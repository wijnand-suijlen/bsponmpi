#include "unbuf.h"

#include <mpi.h>

using namespace bsplib;

void test_1( MPI_Comm comm, int pid, int nprocs )
{
    Unbuf unbuf( 10, comm );

    char msg1[] = "Hi, dit is een test";
    char msg2[] = "Boe";
    char msg3[] = "Dit is een hele lange zin om meerdere berichten te testen";
    char * msgs[3] = { &msg1[0],&msg2[0], &msg3[0] };

    std::vector< char > recv_buf( sizeof(msg3)*nprocs );
   
    for (int p = 0; p < nprocs; ++p ){
        int tag = pid * nprocs + p;
        switch( p % 3) {
            case 0: unbuf.send( tag, p, msg1, sizeof(msg1) ); break;
            case 1: unbuf.send( tag, p, msg2, sizeof(msg2) ); break;
            case 2: unbuf.send( tag, p, msg3, sizeof(msg3) ); break;
        }
    }

    size_t offset = 0;
    for (int p = 0; p < nprocs; ++p ) {
        int tag = p * nprocs + pid;
        switch( pid % 3) {
            case 0: unbuf.recv( tag, p, recv_buf.data()+offset, sizeof(msg1) );
                    offset += sizeof(msg1);
                    break;
            case 1: unbuf.recv( tag, p, recv_buf.data()+offset, sizeof(msg2) );
                    offset += sizeof(msg2);
                    break;
            case 2: unbuf.recv( tag, p, recv_buf.data()+offset, sizeof(msg3) );
                    offset += sizeof(msg3);
                    break;
        }
    }

    unbuf.start();
    unbuf.wait();

    offset = 0;
    for (int p = 0; p < nprocs; ++p ) {
        int tag = p * nprocs + pid;
        char * ptr = recv_buf.data()+offset;
        size_t size;
        switch( pid % 3) {
            case 0: size = sizeof(msg1); break;
            case 1: size = sizeof(msg2); break;
            case 2: size =sizeof(msg3); break;
        }
        offset += size;
 
        if (memcmp(ptr, msgs[pid%3], size)) {
            printf("[%d] wrong data from process %d with tag %d\n", pid, p, tag );
            printf("[%d] Expected '%s' from %d\n", pid, msgs[p], p );
            printf("[%d] Got: '%*s'\n", pid, (int) size, ptr );
            std::abort();
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
