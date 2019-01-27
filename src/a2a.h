#ifndef BSPONMPI_A2A_H
#define BSPONMPI_A2A_H

#include <list>
#include <vector>
#include <numeric>
#include <limits>

#include <mpi.h>

#include "uintserialize.h"
#include "dllexport.h"

namespace bsplib {
 
class DLL_LOCAL A2A
{
public:
    A2A( MPI_Comm comm,
            size_t max_msg_size = std::numeric_limits<int>::max(),
            size_t small_a2a_size_per_proc = 1024 );
    ~A2A();

    int pid() const { return m_pid; }
    int nprocs() const { return m_nprocs; }

    void * send( int dst_pid, const void * data, size_t size );
    size_t send_size( int dst_pid ) const
    { return m_send_sizes[ dst_pid ]; }

    const void * recv_top( int src_pid ) const {
        size_t o = m_recv_offsets[src_pid];
        return static_cast< const void *>( &m_recv_bufs[src_pid][o] );
    }

    bool recv_pop( int src_pid, size_t size ){
        size_t o = m_recv_offsets[src_pid];
        if ( m_recv_sizes[src_pid] - o  >= size ) {
            m_recv_offsets[src_pid] += size;
            return true;
        }
        else {
            return false;
        }
    }

    void recv_rewind( int src_pid ) {
        m_recv_offsets[src_pid] = 0;
    }

    size_t recv_pos( int src_pid ) const 
    { return m_recv_offsets[ src_pid ]; }
    
    size_t recv_size( int src_pid ) const
    { return m_recv_sizes[ src_pid ] - m_recv_offsets[ src_pid ]; }
    
    void exchange();

    void clear();

private:
    A2A( const A2A & ); // copying prohibited
    A2A & operator=(const A2A & ); //assignment prohibited

    int m_pid;
    int m_nprocs;
    const size_t m_max_msg_size;
    const size_t m_small_a2a_size_per_proc;

    std::vector< size_t > m_send_sizes, m_send_offsets;
    std::vector< std::vector< char > > m_send_bufs;
    std::vector< size_t > m_recv_sizes, m_recv_offsets;
    std::vector< std::vector< char > > m_recv_bufs;
    std::vector< MPI_Request > m_reqs;
    std::vector< int > m_ready;
    std::vector< char > m_small_send_buf, m_small_recv_buf;
    MPI_Comm m_comm;
};

template <typename UInt>
void serial( A2A & a2a, int pid, UInt x ) {
    typedef UIntSerialize< UInt > S;     
    typename S::Buffer buf;
    { const int n = S::write(x, buf );
      a2a.send( pid, buf, n );
    }
}

template <typename UInt>
void deserial( A2A & a2a, int pid, UInt & x ) {
    typedef UIntSerialize< UInt > S;     
    const unsigned char * m 
        = static_cast<const unsigned char *>(a2a.recv_top(pid));
    { const int n = S::read( m , x );
      a2a.recv_pop( pid, n );
    }
}


}

#endif
