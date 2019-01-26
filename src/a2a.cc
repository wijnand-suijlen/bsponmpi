#include "a2a.h"

#include <cstring>
#include <mpi.h>

namespace bsplib {

A2A::A2A( MPI_Comm comm, size_t max_msg_size, size_t small_a2a_size_per_proc)
    : m_pid( )
    , m_nprocs( )
    , m_max_msg_size( max_msg_size )
    , m_small_a2a_size_per_proc( small_a2a_size_per_proc )
    , m_send_sizes()
    , m_send_offsets()
    , m_send_bufs()
    , m_recv_sizes()
    , m_recv_offsets()
    , m_recv_bufs()
    , m_reqs()
    , m_ready()
    , m_small_send_buf()
    , m_small_recv_buf()
    , m_comm( MPI_COMM_NULL )
{ 
    MPI_Comm_dup( comm, &m_comm );
    MPI_Comm_rank(m_comm, &m_pid );
    MPI_Comm_size(m_comm, &m_nprocs );

    m_send_sizes.resize( 2* m_nprocs );
    m_send_offsets.resize( m_nprocs );
    m_send_bufs.resize( m_nprocs );
    m_recv_sizes.resize( 2*m_nprocs );
    m_recv_offsets.resize( m_nprocs );
    m_recv_bufs.resize( m_nprocs );
    m_reqs.resize( 2*m_nprocs, MPI_REQUEST_NULL );
    m_ready.resize( 2*m_nprocs );
    m_small_send_buf.resize( small_a2a_size_per_proc * m_nprocs );
    m_small_recv_buf.resize( small_a2a_size_per_proc * m_nprocs );
}

A2A::~A2A()
{
    MPI_Comm_free( &m_comm );
    m_comm = MPI_COMM_NULL;
}

void * A2A::send( int dst_pid, const void * data, size_t size )
{
    assert( dst_pid >= 0 );
    assert( dst_pid < m_nprocs );

    size_t offset = m_send_sizes[ dst_pid ];
    m_send_sizes[ dst_pid ] += size;
    if (m_send_bufs[ dst_pid ].capacity() < offset + size ) {
        m_send_bufs[ dst_pid ].reserve( 
                std::max( 2ul * m_send_bufs[ dst_pid ].capacity(),
                    offset + size )
                );
    }
    m_send_bufs[ dst_pid ].resize( m_send_sizes[ dst_pid ] );
    void * send_buf = &m_send_bufs[dst_pid][offset];
    std::memcpy( send_buf , data, size );
    return send_buf;
}

void A2A::exchange( )
{
    // exchange data sizes. 
    size_t max_send = *max_element( m_send_sizes.begin(), m_send_sizes.begin() + m_nprocs );
    for (int p = m_nprocs; p > 0; --p ) {
        m_send_sizes[2*p-1] = max_send;
        m_send_sizes[2*p-2] = m_send_sizes[p-1];
    }
    // In normal cases, Bruck's algorithm will be used, 
    // so this will cost about O( log P )
    MPI_Alltoall( m_send_sizes.data(), 2*sizeof(size_t), MPI_BYTE,
            m_recv_sizes.data(), 2*sizeof(size_t), MPI_BYTE,
            m_comm );

    size_t max_recv = *max_element( m_recv_sizes.begin(), m_recv_sizes.end() );
    for (int p = 0; p < m_nprocs; ++p) {
        m_recv_sizes[p] = m_recv_sizes[2*p];
        m_send_sizes[p] = m_send_sizes[2*p];
    }

    // Ensure correct size memory of recv buffers
    for (int p = 0; p < m_nprocs; ++p ) {
        if ( m_recv_sizes[p] > m_recv_bufs[p].capacity() ) {
            m_recv_bufs[p].reserve(
                    std::max( 2ul * m_recv_bufs[p].capacity(),
                              m_recv_sizes[p]
                            )
                    );
        }

        m_recv_bufs[p].resize( m_recv_sizes[p] );
        m_recv_offsets[p] = 0;
        m_send_offsets[p] = 0;
    }

    if ( max_recv <= m_small_a2a_size_per_proc ) {
        for (int p = 0; p < m_nprocs; ++p ) {
            memcpy( m_small_send_buf.data() + p * max_recv,
                    m_send_bufs[p].data(),
                    m_send_sizes[p] );
        }
        // In small exchanges, Bruck's algorithm will be used again
        MPI_Alltoall( m_small_send_buf.data(), max_recv, MPI_BYTE,
                m_small_recv_buf.data(), max_recv, MPI_BYTE,
                m_comm );

        for (int p = 0; p < m_nprocs; ++p ) {
            memcpy( m_recv_bufs[p].data(),
                    m_small_recv_buf.data() + p * max_recv,
                    m_recv_sizes[p] );
        }
    
        for (int p = 0 ; p < m_nprocs; ++p ) {
            m_send_sizes[p] = 0;
            m_send_offsets[p] = 0;
            m_recv_offsets[p] = 0;
            m_send_bufs[p].clear();
        }
}
    else { // we have a large exchange


        // Do a personalized exchange
        int outcount = MPI_UNDEFINED;
        bool first_time = true;
        do {
            for (int p = 0; p < m_nprocs; ++p ) {
                if (m_reqs[p] != MPI_REQUEST_NULL) continue;

                int recv_size = std::min( m_max_msg_size, m_recv_sizes[p] );

                int tag = 0;
                if (recv_size > 0 )
                    MPI_Irecv( m_recv_bufs[p].data() + m_recv_offsets[p],
                            recv_size, MPI_BYTE, p,
                            tag, m_comm, & m_reqs[p] );
                
                m_recv_sizes[p] -= recv_size;
                m_recv_offsets[p] += recv_size;
            }

            if (first_time)
                MPI_Barrier( m_comm ); // Using the barrier the first time
                                       // allows to use ready sends

            for (int p = 0; p < m_nprocs; ++p ) {
                if (m_reqs[m_nprocs + p] != MPI_REQUEST_NULL) continue;

                int send_size = std::min( m_max_msg_size, m_send_sizes[p] );


                int tag = 0;
                if (send_size > 0 ) {
                    if (first_time)
                     MPI_Irsend( m_send_bufs[p].data() + m_send_offsets[p],
                            send_size, MPI_BYTE, p,
                            tag, m_comm, & m_reqs[m_nprocs + p ] );
                    else
                     MPI_Isend( m_send_bufs[p].data() + m_send_offsets[p],
                            send_size, MPI_BYTE, p,
                            tag, m_comm, & m_reqs[m_nprocs + p ] );
                }
                
                m_send_sizes[p] -= send_size;
                m_send_offsets[p] += send_size;
            }
            
            first_time = false;

            MPI_Waitsome( m_reqs.size(), m_reqs.data(), &outcount, 
                    m_ready.data(), MPI_STATUSES_IGNORE );
        } while (outcount != MPI_UNDEFINED );

        for (int p = 0 ; p < m_nprocs; ++p ) {
            assert( m_recv_offsets[p] == m_recv_bufs[p].size() );
            assert( m_recv_sizes[p] == 0 );
            assert( m_send_sizes[p] == 0 );
            assert( m_send_offsets[p] == m_send_bufs[p].size() );
            assert( m_reqs[p] == MPI_REQUEST_NULL );
            assert( m_reqs[m_nprocs+p] == MPI_REQUEST_NULL );
            m_recv_sizes[p] = m_recv_offsets[p];
            m_recv_offsets[p] = 0;
            m_send_offsets[p] = 0;
            m_send_bufs[p].clear();
        }
    } // end of MPI_Isend - Irecv - wait pattern 

}


}
