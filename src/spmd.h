#ifndef BSPONMPI_SPMD_H
#define BSPONMPI_SPMD_H

#ifdef HAS_CLOCK_GETTIME
  #include <time.h>
#endif

#include <mpi.h>

#include "dllexport.h"

namespace bsplib {

class DLL_LOCAL Spmd {
public:

    explicit Spmd( int nprocs ); // Collective
    ~Spmd();  // Collective

    bool closed() { return m_closed; }
    bool active() { return m_active; }
    int pid() const { return m_pid; }
    int nprocs() const { return m_nprocs; }

#ifdef HAS_CLOCK_GETTIME
    double time() const { 
      struct timespec now;
      clock_gettime( CLOCK_MONOTONIC, &now );
      return (now.tv_sec - m_time.tv_sec) 
               + 1e-9 * (now.tv_nsec - m_time.tv_nsec );
    }
#else
    double time() const { return MPI_Wtime() - m_time; }
#endif

    int normal_sync(); // return non-zero on error
    int end_sync();   // return non-zero on error
    void close();  // must be called before destructor

    MPI_Comm comm() const { return m_comm; }

private:
    MPI_Comm m_comm;
    bool m_closed;
    bool m_active;
#ifdef HAS_CLOCK_GETTIME
    struct timespec m_time;
#else
    double m_time;
#endif
    int m_pid;
    int m_nprocs;
};

}

#endif
