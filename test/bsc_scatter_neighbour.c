#include <bsp.h>
#include <bsc.h>
#include "test.h"

TEST( bsc_scatter_neightbour, success() )
{
    bsc_group_t group;
    bsp_pid_t neighbours[3];
    const int x[3] = { 1, 2, 3 };
    int y = 0;

    bsp_begin( bsp_nprocs() );
    EXPECT_EQ( "%d", 6, bsp_nprocs() ); /* test only makes sense on 6 processes */

    neighbours[0] = ( bsp_nprocs() + bsp_pid() - 1 ) % bsp_nprocs();
    neighbours[1] = bsp_pid();
    neighbours[2] = ( bsp_pid() + 1) % bsp_nprocs();

    group = bsc_group_create_neighbourhood( neighbours, 3 );
    bsp_push_reg( &y, sizeof(y) );
    bsp_sync();

    bsc_scatter( bsc_start, (bsp_pid()/3)*3, group, x, &y, sizeof(x[0]));
    bsc_sync( bsc_flush );


    /* scatter pattern */
    /*   1        1
     *   2        2
     *   3        3 
     * / | \    / | \
     *   2, 3, 1, 2, 3, 1 */

    EXPECT_EQ( "%d", y, x[ (bsp_pid() + 1)%3 ] );

    bsp_pop_reg( &y );

    bsp_end();
}
