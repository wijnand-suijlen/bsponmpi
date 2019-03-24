#include <bsp.h>
#include <bsc.h>
#include "test.h"

TEST( bsc_gather_neightbour, success() )
{
    bsc_group_t group;
    bsp_pid_t neighbours[3];
    int x[3] = { 0, 0, 0 };
    int y = 0;

    bsp_begin( bsp_nprocs() );
    EXPECT_EQ( "%d", 6, bsp_nprocs() ); /* test only makes sense on 6 processes */

    y = bsp_pid() + 1;

    neighbours[0] = ( bsp_nprocs() + bsp_pid() - 2 ) % bsp_nprocs();
    neighbours[1] = bsp_pid();
    neighbours[2] = ( bsp_pid() + 2) % bsp_nprocs();

    group = bsc_group_create_neighbourhood( neighbours, 3 );
    bsp_push_reg( &y, sizeof(y) );
    bsp_sync();

    bsc_gather( bsc_start, (bsp_pid()/3)*3, group, &y, x, sizeof(x[0]));
    bsc_sync( bsc_flush );


    /* gather pattern */
    /*   5        2
     *   1    --- 4----
     * --3---/    6   -\-
     *   |  / \   |  /  \
     *   1, 2, 3, 4, 5, 6 */

    if ( bsp_pid() == 0 ) {
        EXPECT_EQ( "%d", x[0], 5 );
        EXPECT_EQ( "%d", x[1], 1 );
        EXPECT_EQ( "%d", x[2], 3 );
    }
    else if (bsp_pid() == 3 ) {
        EXPECT_EQ( "%d", x[0], 2 );
        EXPECT_EQ( "%d", x[1], 4 );
        EXPECT_EQ( "%d", x[2], 6 );
    }
    else {
        EXPECT_EQ( "%d", x[0], 0 );
        EXPECT_EQ( "%d", x[1], 0 );
        EXPECT_EQ( "%d", x[2], 0 );
    }



    bsp_pop_reg( &y );

    bsp_end();
}
