#include "a2a.h"

#ifdef PROFILE
#include "tictoc.h"
#endif

#include <algorithm>
#include <cstring>
#include <mpi.h>

namespace bsplib {

A2A::A2A( MPI_Comm comm, std::size_t max_msg_size, std::size_t small_a2a_size_per_proc)
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
    , m_new_recv_cap()
    , m_comm( MPI_COMM_NULL )
#ifdef USE_ONESIDED
    , m_renew_win(true)
    , m_recv_win()
#endif
{ 
    MPI_Comm_dup( comm, &m_comm );
    MPI_Comm_rank(m_comm, &m_pid );
    MPI_Comm_size(m_comm, &m_nprocs );

    m_send_sizes.resize( 4* m_nprocs );
    m_send_offsets.resize( m_nprocs );
    m_send_bufs.resize( m_nprocs );
    m_recv_sizes.resize( 4*m_nprocs );
    m_recv_offsets.resize( m_nprocs );
    m_recv_bufs.resize( m_nprocs );
    m_reqs.resize( 2*m_nprocs, MPI_REQUEST_NULL );
    m_ready.resize( 2*m_nprocs );
    m_small_send_buf.resize( small_a2a_size_per_proc * m_nprocs );
    m_small_recv_buf.resize( small_a2a_size_per_proc * m_nprocs );
    m_new_recv_cap.resize( m_nprocs );
#ifdef USE_ONESIDED
    m_recv_win.resize( m_nprocs, MPI_WIN_NULL );
#endif
}

A2A::~A2A()
{
#ifdef USE_ONESIDED
    for (int p = 0; p < m_nprocs; ++p ) {
        if ( m_recv_win[p] != MPI_WIN_NULL )
            MPI_Win_free( &m_recv_win[p] );
        m_recv_win[p] = MPI_WIN_NULL;
    }
#endif
    MPI_Comm_free( &m_comm );
    m_comm = MPI_COMM_NULL;
}

void * A2A::send( int dst_pid, const void * data, std::size_t size )
{
    assert( dst_pid >= 0 );
    assert( dst_pid < m_nprocs );

    std::size_t offset = m_send_sizes[ dst_pid ];
    m_send_sizes[ dst_pid ] += size;
    if (m_send_bufs[ dst_pid ].capacity() < offset + size ) {
#ifdef USE_ONESIDED
        m_renew_win = true;
#endif
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

void A2A :: clear()
{
    for (int p = 0 ; p < m_nprocs; ++p ) {
        m_recv_sizes[p] = 0;
        m_send_sizes[p] = 0;
        m_recv_offsets[p] = 0;
        m_send_offsets[p] = 0;
        m_send_bufs[p].clear();
    }
}


void A2A::exchange( )
{
#ifdef USE_ONESIDED
    bool renew_win = false;
#endif
    std::size_t max_recv = 0;
    {   
#ifdef PROFILE
        TicToc t( TicToc::MPI_META_A2A, 4*sizeof(std::size_t)*m_nprocs );
#endif
        // exchange data sizes. 
        std::size_t max_send = *std::max_element( m_send_sizes.begin(),
                               m_send_sizes.begin() + m_nprocs );
        for (int p = m_nprocs; p > 0; --p ) {
#ifdef USE_ONESIDED
            m_send_sizes[4*p-1] = m_renew_win?1:0;
#endif
            m_send_sizes[4*p-2] = m_send_bufs[p-1].capacity();
            m_send_sizes[4*p-3] = max_send;
            m_send_sizes[4*p-4] = m_send_sizes[p-1];
        }
        // In normal cases, Bruck's algorithm will be used, 
        // so this will cost about O( log P )
        MPI_Alltoall( m_send_sizes.data(), 4*sizeof(std::size_t), MPI_BYTE,
                m_recv_sizes.data(), 4*sizeof(std::size_t), MPI_BYTE,
                m_comm );

        for (int p = 0; p < m_nprocs; ++p) {
#ifdef USE_ONESIDED
            renew_win = renew_win || m_recv_sizes[4*p+3]==1;
#endif
            m_new_recv_cap[p] = m_recv_sizes[4*p+2];
            max_recv = std::max( max_recv, m_recv_sizes[4*p+1] );
            m_recv_sizes[p] = m_recv_sizes[4*p];
            m_send_sizes[p] = m_send_sizes[4*p];
        }
    }

    for (int p = 0; p < m_nprocs; ++p ) {
#ifdef USE_ONESIDED
        if (renew_win && m_recv_win[p] != MPI_WIN_NULL )
            MPI_Win_free( &m_recv_win[p]);
#endif
    // Ensure correct size memory of recv buffers
        if ( m_new_recv_cap[p] > m_recv_bufs[p].capacity() ) {
            m_recv_bufs[p].reserve( m_new_recv_cap[p] );
        }
#ifdef USE_ONESIDED
        if (renew_win) {
            MPI_Win_create( m_recv_bufs[p].data(), 
                    m_recv_bufs[p].capacity(), 1, 
                    MPI_INFO_NULL, m_comm, &m_recv_win[p]);
        }
#endif

        assert( m_recv_bufs[p].capacity() >= m_recv_sizes[p] );
        m_recv_bufs[p].resize( m_recv_sizes[p] );
        m_recv_offsets[p] = 0;
        m_send_offsets[p] = 0;
    }

#ifdef USE_ONESIDED
    m_renew_win = false;
#endif

    if ( max_recv == 0 ) {
        /* no need to do anything */
        clear();
    }
#ifdef USE_ONESIDED
    else {
#ifdef PROFILE
        TicToc t( TicToc::MPI_PUT );
#endif
        for (int p = 0; p < m_nprocs; ++p)
            MPI_Win_fence( 0, m_recv_win[p] );

        for (int p = 0; p < m_nprocs; ++p ) {
            std::size_t size = m_send_sizes[p];
            std::size_t o = 0;
#ifdef PROFILE
            t.addBytes( size );
#endif
            while ( size > 0 ) {
                std::size_t s = std::min( m_max_msg_size, size );
                //size_t s = size;
                MPI_Put( m_send_bufs[p].data() + o, s, MPI_BYTE,
                        p, o, s, MPI_BYTE, m_recv_win[m_pid] );
                size -= s;
                o += s;
            }
        }

        for (int p = 0; p < m_nprocs; ++p)
            MPI_Win_fence( 0, m_recv_win[p] );


        for (int p = 0 ; p < m_nprocs; ++p ) {
            m_send_sizes[p] = 0;
            m_send_offsets[p] = 0;
            m_recv_offsets[p] = 0;
            m_send_bufs[p].clear();
        }
    }
#else
    else if ( max_recv <= m_small_a2a_size_per_proc ) {
#ifdef PROFILE
        TicToc t( TicToc::MPI_SMALL_A2A );
#endif
        for (int p = 0; p < m_nprocs; ++p ) {
            memcpy( m_small_send_buf.data() + p * max_recv,
                    m_send_bufs[p].data(),
                    m_send_sizes[p] );
#ifdef PROFILE
            t.addBytes( m_send_sizes[p] );
#endif
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
#ifdef PROFILE
        TicToc tr( TicToc::MPI_LARGE_RECV );
        TicToc ts( TicToc::MPI_LARGE_SEND );
#endif

        // Do a personalized exchange
        int outcount = MPI_UNDEFINED;
        bool first_time = true;
        do {
            for (int p = 0; p < m_nprocs; ++p ) {
                if (m_reqs[p] != MPI_REQUEST_NULL) continue;

                int recv_size = std::min( m_max_msg_size, m_recv_sizes[p] );
#ifdef PROFILE
                tr.addBytes( recv_size );
#endif

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
#ifdef PROFILE
                ts.addBytes( send_size );
#endif

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
#endif

}


}
