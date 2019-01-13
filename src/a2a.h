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
    A2A( MPI_Comm comm,
            size_t max_msg_size = std::numeric_limits<int>::max() );
    ~A2A();

    int pid() const { return m_pid; }
    int nprocs() const { return m_nprocs; }

    void send( int dst_pid, const void * data, size_t size );

    const void * recv_top( int src_pid ) const;
    bool recv_pop( int src_pid, size_t size );
    size_t recv_size( int src_pid ) const
    { return m_recv_sizes[ src_pid ] - m_recv_offsets[ src_pid ]; }
    
    void exchange();

private:
    A2A( const A2A & ); // copying prohibited
    A2A & operator=(const A2A & ); //assignment prohibited

    int m_pid;
    int m_nprocs;
    const size_t m_max_msg_size;

    std::vector< size_t > m_send_sizes, m_send_offsets;
    std::vector< std::vector< char > > m_send_bufs;
    std::vector< size_t > m_recv_sizes, m_recv_offsets;
    std::vector< std::vector< char > > m_recv_bufs;
    std::vector< MPI_Request > m_reqs;
    std::vector< int > m_ready;
    MPI_Comm m_comm;
};

}

#endif
