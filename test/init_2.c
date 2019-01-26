#include <bsp.h>
#include "test.h"

#include <mpi.h>

static int s_some_int;
static bsp_pid_t s_nprocs ;

void spmd(void)
{
    bsp_pid_t s;
    int r;
    bsp_begin( bsp_nprocs() / 2 );

    s = bsp_pid();
    MPI_Comm_rank( MPI_COMM_WORLD, &r);
    EXPECT_IMPLIES( s==0, s_some_int == 5 );

    EXPECT_OP( "%d", r, <, bsp_nprocs() );

    if (s == 0 )
        s_some_int = (int) s - 10;
    else
        s_some_int = 999;

    bsp_end();
}

TEST( init_2, success() )
{
    int pid1, pid2;
    bsp_init( spmd, 0, NULL);

    MPI_Comm_rank( MPI_COMM_WORLD, &pid1 );

    s_some_int = 5;
    s_nprocs = bsp_nprocs();
    EXPECT_OP( "%d", pid1, <, s_nprocs );

    spmd();

    MPI_Comm_rank( MPI_COMM_WORLD, &pid2 );
    EXPECT_EQ( "%d", pid2, 0 );

    EXPECT_EQ("%d", bsp_nprocs(), s_nprocs );

    EXPECT_EQ("%d", s_some_int, -10 );
}
