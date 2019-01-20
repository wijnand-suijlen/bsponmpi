#include "a2a.h"

#include <mpi.h>
#include <cstdio>
#include <cstring>
#include <cstdlib>

using namespace bsplib;

void test_1( MPI_Comm comm, int pid, int nprocs )
{
    A2A a2a( comm, 10 );

    char msg1[] = "Hi, dit is een test";
    char msg2[] = "Boe";
    char msg3[] = "Dit is een hele lange zin om meerdere berichten te testen";

    a2a.send( 0, msg1, sizeof(msg1) );
    a2a.send( 1 % nprocs, msg2, sizeof(msg2) );
    a2a.send( 2 % nprocs, msg3, sizeof(msg3) );

    a2a.exchange( );

    char recv1[ sizeof(msg1) ];
    char recv2[ sizeof(msg2) ];
    char recv3[ sizeof(msg3) ];

    size_t rs[] = { sizeof(msg1), sizeof(msg2), sizeof(msg3) };
    char *recv[] = { &recv1[0], &recv2[0], &recv3[0] };
    std::memcpy( recv[0], msg1, sizeof(msg1) );
    std::memcpy( recv[1], msg2, sizeof(msg2) );
    std::memcpy( recv[2], msg3, sizeof(msg3) );

    if ( pid <= 2 ) {
        for (int p = 0; p < nprocs; ++p ) {
            int i = pid;
            const void * data = NULL;
            int m = 0;
            while ( data = a2a.recv_top(p), a2a.recv_pop( p, rs[i] ) ) {
                if (std::memcmp( data, recv[i], rs[i] ) != 0 ) {
                    std::printf("[%d] wrong data from process %d\n", pid, p );
                    std::printf("[%d] Expected '%s' from %d\n", pid, recv[i], p );
                    std::printf("[%d] Got: '%*s'\n", pid, (int) rs[i], (const char *) data );
                    std::abort();
                }

                i = (i + std::min(3,nprocs)) % 3;
                m += 1;
            }
        }
    }
    else {
        const void * data = NULL;
        for (int p = 0; p < nprocs; ++p ) {
            data = a2a.recv_top(p);
            if ( a2a.recv_pop(p, 1) ) {
                std::printf("[%d] Unexpected message from %d\n", pid, p );
                std::abort();
            }
        }
    }
}

int main( int argc, char ** argv ) 
{
    MPI_Init( &argc, &argv );

    MPI_Comm comm = MPI_COMM_WORLD;

    int pid, nprocs;
    MPI_Comm_rank( comm, & pid );
    MPI_Comm_size( comm, & nprocs );


    test_1( comm, pid, nprocs );


    MPI_Finalize( );
    return 0;
}
