#include "a2a.h"

#include <cstring>

namespace bsplib {

A2A::A2A( int pid, int nprocs, size_t max_msg_size )
    : m_pid( pid )
    , m_nprocs( nprocs )
    , m_max_msg_size( max_msg_size )
    , m_send_sizes( nprocs )
    , m_send_offsets( nprocs )
    , m_send_bufs( nprocs )
    , m_recv_sizes( nprocs )
    , m_recv_offsets( nprocs )
    , m_recv_bufs( nprocs )
    , m_reqs( 2*nprocs, MPI_REQUEST_NULL )
    , m_ready( 2*nprocs )
{ }

void A2A::send( int dst_pid, const void * data, size_t size )
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
    std::memcpy( &m_send_bufs[dst_pid][offset], data, size );
}

bool A2A::recv( int src_pid, const void *& data, size_t size ) {
    assert( src_pid >= 0 );
    assert( src_pid < m_nprocs );

    size_t o = m_recv_offsets[src_pid];
    if ( m_recv_sizes[src_pid] - o  >= size ) {
        data = static_cast< const void *>( &m_recv_bufs[src_pid][o] );
        m_recv_offsets[src_pid] += size;
        return true;
    }
    else {
        return false;
    }
}

void A2A::exchange( MPI_Comm comm ) 
{
    // exchange data sizes. 
    // In normal cases, Bruck's algorithm will be used, 
    // so this will cost about O( log P )
    MPI_Alltoall( m_send_sizes.data(), sizeof(size_t), MPI_BYTE,
            m_recv_sizes.data(), sizeof(size_t), MPI_BYTE,
            comm );

    // Ensure correct size memory buffers
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


    // Do a personalized exchange
    int outcount = MPI_UNDEFINED;
    do {
        for (int p = 0; p < m_nprocs; ++p ) {
            if (m_reqs[p] != MPI_REQUEST_NULL) continue;

            int recv_size = std::min( m_max_msg_size, m_recv_sizes[p] );

            int tag = 0;
            if (recv_size > 0 )
                MPI_Irecv( m_recv_bufs[p].data() + m_recv_offsets[p],
                        recv_size, MPI_BYTE, p,
                        tag, comm, & m_reqs[p] );
            
            m_recv_sizes[p] -= recv_size;
            m_recv_offsets[p] += recv_size;
        }

        for (int p = 0; p < m_nprocs; ++p ) {
            if (m_reqs[m_nprocs + p] != MPI_REQUEST_NULL) continue;

            int send_size = std::min( m_max_msg_size, m_send_sizes[p] );


            int tag = 0;
            if (send_size > 0 )
                MPI_Isend( m_send_bufs[p].data() + m_send_offsets[p],
                        send_size, MPI_BYTE, p,
                        tag, comm, & m_reqs[m_nprocs + p ] );
            
            m_send_sizes[p] -= send_size;
            m_send_offsets[p] += send_size;
        }

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

}


}
