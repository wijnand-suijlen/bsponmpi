#include <bsp.h>
#include <bsc.h>
#include "test.h"

TEST( bsc_allgather_neighbour, success() )
{
    bsc_group_t group;
    bsp_pid_t neighbours[3];
    int x[3] = { 0, 0, 0 };
    int y = 0;
    int i;

    bsp_begin( bsp_nprocs() );
    EXPECT_EQ( "%d", 6, bsp_nprocs() ); /* test only makes sense on 6 processes */

    y = bsp_pid() + 1;

    neighbours[0] = ( bsp_nprocs() + bsp_pid() - 2 ) % bsp_nprocs();
    neighbours[1] = bsp_pid();
    neighbours[2] = ( bsp_pid() + 2) % bsp_nprocs();

    group = bsc_group_create_neighbourhood( neighbours, 3 );
    bsp_push_reg( &y, sizeof(y) );
    bsp_sync();

    bsc_allgather( bsc_start, group, &y, x, sizeof(x[0]));
    bsc_sync( bsc_flush );


    /* gather pattern */
    /*   5        2
     *   1    --- 4----
     * --3---/    6   -\-
     *   |  / \   |  /  \
     *   1, 2, 3, 4, 5, 6 */

    for ( i = 0 ; i < 6; ++i ) {
        int left = (bsp_pid() + 4) % 6+1;
        int mid = bsp_pid()+1;
        int right = (bsp_pid() + 2) % 6+1;

        EXPECT_EQ( "%d", x[0], left);
        EXPECT_EQ( "%d", x[1], mid );
        EXPECT_EQ( "%d", x[2], right );
    }



    bsp_pop_reg( &y );

    bsp_end();
}
