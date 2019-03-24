#include <bsp.h>
#include <bsc.h>
#include "test.h"

TEST( bsc_alltoall_neighbour, success() )
{
    bsc_group_t group;
    bsp_pid_t neighbours[3];
    int x[3] = { 0, 0, 0 };
    int y[3] = { 0, 0, 0 };
    int left, mid, right;

    bsp_begin( bsp_nprocs() );
    EXPECT_EQ( "%d", 6, bsp_nprocs() ); /* test only makes sense on 6 processes */

    y[0] = bsp_pid() + 1;
    y[1] = bsp_pid() + 7;
    y[2] = bsp_pid() + 13;

    neighbours[0] = ( bsp_nprocs() + bsp_pid() - 2 ) % bsp_nprocs();
    neighbours[1] = bsp_pid();
    neighbours[2] = ( bsp_pid() + 2) % bsp_nprocs();

    group = bsc_group_create_neighbourhood( neighbours, 3 );
    bsp_push_reg( &x[0], sizeof(x) );
    bsp_sync();

    bsc_alltoall( bsc_start, group, y, x, sizeof(x[0]));
    bsc_sync( bsc_flush );


    /* alltoall pattern */
    
     /*
     * Sent
     * [0] ( 1, 7, 13 )
     * [1] ( 2, 8, 14 )
     * [2] ( 3, 9, 15 )
     * [3] ( 4, 10, 16 )
     * [4] ( 5, 11, 17 )
     * [5] ( 6, 12, 18 )
     *
     *
     * Received
     * [0] ( 17, 7, 3 )
     * [1] ( 18, 8, 4 )
     * [2] ( 13, 9, 5 )
     * [3] ( 14, 10, 6 )
     * [4] ( 15, 11, 1 )
     * [5] ( 16, 12, 2 )
     *
     *
     */

    left = (bsp_pid() + 4) % 6+13;
    mid  =  bsp_pid()+ 7;
    right = (bsp_pid() + 2) % 6+1;

    EXPECT_EQ( "%d", x[0], left);
    EXPECT_EQ( "%d", x[1], mid );
    EXPECT_EQ( "%d", x[2], right );
    
    bsp_pop_reg( &x );

    bsp_end();
}
