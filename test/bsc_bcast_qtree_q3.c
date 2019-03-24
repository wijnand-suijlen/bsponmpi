#include <bsp.h>
#include <bsc.h>
#include "test.h"

bsc_step_t bsc_bcast_qtree( bsc_step_t depends,
                      bsp_pid_t root, bsc_group_t group, 
                      const void * src, void * dst, bsc_size_t size,
                      bsc_pid_t q );

TEST( bsc_bcast_qtree_q3, success() )
{
    int x = 0;
    int y = 0;
    bsp_begin( bsp_nprocs() );

    bsp_push_reg( &y, sizeof(y) );
    bsp_sync();

    if ( bsp_pid() == 2 % bsp_nprocs())
        x = 5;

    bsc_bcast_qtree( bsc_start, 2 % bsp_nprocs(), bsc_all,
            &x, &y, sizeof(x), 3 );

    bsc_sync( bsc_flush );

    EXPECT_EQ( "%d", y, 5 );
    EXPECT_EQ( "%d", x,  (bsp_pid() == 2 % bsp_nprocs()?5:0) );

    bsp_pop_reg( &y );

    bsp_end();
}
