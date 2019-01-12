#ifndef BSPONMPI_A2A_H
#define BSPONMPI_A2A_H

#include <list>
#include <vector>
#include <numeric>
#include <limits>

#include <mpi.h>

namespace bsplib {
 
class A2A
{
public:
    A2A( int pid, int nprocs,
            size_t max_msg_size = std::numeric_limits<int>::max() );
    
    void send( int dst_pid, const void * data, size_t size );
    bool recv( int src_pid, const void ** data, size_t size );
    
    void exchange( MPI_Comm comm ) ;

private:
    const int m_pid;
    const int m_nprocs;
    const size_t m_max_msg_size;

    std::vector< size_t > m_send_sizes, m_send_offsets;
    std::vector< std::vector< char > > m_send_bufs;
    std::vector< size_t > m_recv_sizes, m_recv_offsets;
    std::vector< std::vector< char > > m_recv_bufs;
    std::vector< MPI_Request > m_reqs;
    std::vector< int > m_ready;
};

}

#endif
