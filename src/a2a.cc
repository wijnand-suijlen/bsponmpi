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
    , m_send_cap(0)
    , m_recv_cap(0)
    , m_send_sizes()
    , m_send_pos()
    , m_send_bufs()
    , m_recv_sizes()
    , m_recv_pos()
    , m_recv_bufs()
    , m_small_send_buf()
    , m_small_recv_buf()
    , m_comm( MPI_COMM_NULL )
#ifdef USE_ONESIDED
    , m_recv_win(MPI_WIN_NULL)
#else
    , m_reqs()
    , m_ready()
#endif
{ 
    MPI_Comm_dup( comm, &m_comm );
    MPI_Comm_rank(m_comm, &m_pid );
    MPI_Comm_size(m_comm, &m_nprocs );

    m_small_a2a_size_per_proc = std::min<size_t>( m_small_a2a_size_per_proc,
            std::numeric_limits<int>::max()/m_nprocs );

    m_send_sizes.resize( 3* m_nprocs );
    m_send_pos.resize( m_nprocs );
    m_recv_sizes.resize( 3*m_nprocs );
    m_recv_pos.resize( m_nprocs );
    m_small_send_buf.resize( small_a2a_size_per_proc * m_nprocs );
    m_small_recv_buf.resize( small_a2a_size_per_proc * m_nprocs );

#if ! defined USE_ONESIDED
    m_reqs.resize( 2*m_nprocs, MPI_REQUEST_NULL );
    m_ready.resize( 2*m_nprocs );
#endif
}

A2A::~A2A()
{
#ifdef USE_ONESIDED
    if ( m_recv_win != MPI_WIN_NULL )
        MPI_Win_free( &m_recv_win );
#endif
    MPI_Comm_free( &m_comm );
    m_comm = MPI_COMM_NULL;
}

void * A2A::send( int dst_pid, const void * data, std::size_t size )
{
    assert( dst_pid >= 0 );
    assert( dst_pid < m_nprocs );

    assert( m_send_cap == m_send_bufs.size() / m_nprocs );

    std::size_t offset = m_send_sizes[ dst_pid ];
    if ( m_send_cap < offset + size ) {
        std::size_t new_cap =
           std::max( 2ul * m_send_cap , offset + size );

        m_send_bufs.resize( m_nprocs * new_cap );
        for ( int p = m_nprocs; p > 0; --p ) {
            std::size_t displ = new_cap - m_send_cap;

            for ( size_t i = p*new_cap; i > (p-1)*new_cap; --i) {
                m_send_bufs[i-1] = m_send_bufs[i-displ*(p-1)-1];
            }
        }  
        m_send_cap = new_cap;
    }
    m_send_sizes[ dst_pid ] += size;
    void * send_buf = &m_send_bufs[ dst_pid * m_send_cap + offset];
    std::memcpy( send_buf , data, size );
    return send_buf;
}

void A2A :: clear()
{
    for (int p = 0 ; p < m_nprocs; ++p ) {
        m_recv_sizes[p] = 0;
        m_send_sizes[p] = 0;
        m_recv_pos[p] = 0;
        m_send_pos[p] = 0;
    }
}


void A2A::exchange( )
{
    std::size_t max_recv = 0;
    std::size_t new_cap = m_send_cap;
    {   
#ifdef PROFILE
        TicToc t( TicToc::MPI_META_A2A, 3*sizeof(std::size_t)*m_nprocs );
#endif
        // exchange data sizes. 
        std::size_t max_send = *std::max_element( m_send_sizes.begin(),
                               m_send_sizes.begin() + m_nprocs );
        for (int p = m_nprocs; p > 0; --p ) {
            m_send_sizes[3*p-1] = m_send_cap;
            m_send_sizes[3*p-2] = max_send;
            m_send_sizes[3*p-3] = m_send_sizes[p-1];
        }
        // In normal cases, Bruck's algorithm will be used, 
        // so this will cost about O( log P )
        MPI_Alltoall( m_send_sizes.data(), 3*sizeof(std::size_t), MPI_BYTE,
                m_recv_sizes.data(), 3*sizeof(std::size_t), MPI_BYTE,
                m_comm );

        for (int p = 0; p < m_nprocs; ++p) {
            new_cap  = std::max( new_cap, m_recv_sizes[3*p+2] );
            max_recv = std::max( max_recv, m_recv_sizes[3*p+1] );
            m_recv_sizes[p] = m_recv_sizes[3*p];
            m_send_sizes[p] = m_send_sizes[3*p];
        }
    }

    // Ensure correct size memory of recv buffers
#ifdef USE_ONESIDED
    if ( new_cap != m_recv_cap && m_recv_win != MPI_WIN_NULL )
        MPI_Win_free( &m_recv_win );
#endif
    m_recv_bufs.resize( new_cap * m_nprocs );

#ifdef USE_ONESIDED
    if ( new_cap != m_recv_cap ) {
        MPI_Win_create( m_recv_bufs.data(), 
                new_cap, 1, 
                MPI_INFO_NULL, m_comm, &m_recv_win );
    }
#endif
    m_recv_cap = new_cap;

    assert( m_recv_cap >= max_recv );

    if ( max_recv == 0 ) {
        /* no need to do anything */
        clear();
    }
    else if ( max_recv <= m_small_a2a_size_per_proc ) {
#ifdef PROFILE
        TicToc t( TicToc::MPI_SMALL_A2A );
#endif
        for (int p = 0; p < m_nprocs; ++p ) {
            memcpy( m_small_send_buf.data() + p * max_recv,
                    m_send_bufs.data() + p * m_send_cap,
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
            memcpy( m_recv_bufs.data() + p * m_recv_cap,
                    m_small_recv_buf.data() + p * max_recv,
                    m_recv_sizes[p] );
        }
    }
#ifdef USE_ONESIDED
    else { // we have a large exchange
#ifdef PROFILE
        TicToc t( TicToc::MPI_PUT );
#endif
        assert( m_recv_cap > 0 );
        assert( m_recv_win != MPI_WIN_NULL );

        MPI_Win_fence( 0, m_recv_win );

        for (int p = 0; p < m_nprocs; ++p ) {
            std::size_t size = m_send_sizes[p];
            std::size_t o1 = m_send_cap * p;
            std::size_t o2 = m_recv_cap * m_pid;
#ifdef PROFILE
            t.addBytes( size );
#endif
            while ( size > 0 ) {
                std::size_t s = std::min( m_max_msg_size, size );
                //size_t s = size;
                MPI_Put( m_send_bufs.data() + o1, s, MPI_BYTE,
                        p, o2, s, MPI_BYTE, m_recv_win );
                size -= s;
                o1 += s;
                o2 += s;
            }
        }

        MPI_Win_fence( 0, m_recv_win );
    }
#else
    else { // we have a large exchange
#ifdef PROFILE
        TicToc tr( TicToc::MPI_LARGE_RECV );
        TicToc ts( TicToc::MPI_LARGE_SEND );
#endif

        for ( int p = 0 ; p < m_nprocs; ++p ) {
            m_send_pos[p] = p * m_send_cap;
            m_recv_pos[p] = p * m_recv_cap;
        }

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
                    MPI_Irecv( m_recv_bufs.data() + m_recv_pos[p],
                            recv_size, MPI_BYTE, p,
                            tag, m_comm, & m_reqs[p] );
                
                m_recv_sizes[p] -= recv_size;
                m_recv_pos[p] += recv_size;
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
                     MPI_Irsend( m_send_bufs.data() + m_send_pos[p],
                            send_size, MPI_BYTE, p,
                            tag, m_comm, & m_reqs[m_nprocs + p ] );
                    else
                     MPI_Isend( m_send_bufs.data() + m_send_pos[p],
                            send_size, MPI_BYTE, p,
                            tag, m_comm, & m_reqs[m_nprocs + p ] );
                }
                
                m_send_sizes[p] -= send_size;
                m_send_pos[p] += send_size;
            }
            
            first_time = false;

            MPI_Waitsome( m_reqs.size(), m_reqs.data(), &outcount, 
                    m_ready.data(), MPI_STATUSES_IGNORE );
        } while (outcount != MPI_UNDEFINED );

        for (int p = 0 ; p < m_nprocs; ++p ) {
            assert( m_recv_sizes[p] == 0 );
            assert( m_send_sizes[p] == 0 );
            assert( m_reqs[p] == MPI_REQUEST_NULL );
            assert( m_reqs[m_nprocs+p] == MPI_REQUEST_NULL );
            m_recv_sizes[p] = m_recv_pos[p] - p * m_recv_cap;
        }
    } // end of MPI_Isend - Irecv - wait pattern 
#endif

    m_send_cap = m_recv_cap;
    m_send_bufs.resize( m_nprocs * m_send_cap );
    for (int p = 0 ; p < m_nprocs; ++p ) {
        m_send_sizes[p] = 0;
        m_send_pos[p] = 0;
        m_recv_pos[p] = 0;
    }
}


}
