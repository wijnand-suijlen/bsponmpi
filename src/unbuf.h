#ifndef BSPONMPI_UNBUF_H
#define BSPONMPI_UNBUF_H

#include <vector>
#include <mpi.h>

namespace bsplib {

class Unbuf {
public:
    explicit Unbuf( size_t max_msg_size, MPI_Comm comm );
    ~Unbuf();

    void send( int dst_pid, const void * addr, size_t size );
    void recv( int src_pid, void * addr, size_t size );

    void start();
    void wait();

private:
    Unbuf( const Unbuf & ); //copying prohibited
    Unbuf & operator=( const Unbuf & ); // assignment prohibited

    struct Entry {
        int pid;
        char * addr;
        size_t size;
    };

    MPI_Comm m_comm;
    size_t m_max_msg_size;
    std::vector< Entry > m_sends;
    std::vector< Entry > m_recvs;
    std::vector< MPI_Request > m_reqs;
    std::vector< int > m_ready;
};

}

#endif
