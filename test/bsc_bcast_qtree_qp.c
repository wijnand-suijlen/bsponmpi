#include <bsp.h>
#include <bsc.h>
#include "test.h"

bsc_step_t bsc_bcast_1phase( bsc_step_t depends,
                      bsp_pid_t root, bsc_group_t group, 
                      const void * src, void * dst, bsc_size_t size );

bsc_step_t bsc_bcast_qtree( bsc_step_t depends,
                      bsp_pid_t root, bsc_group_t group, 
                      const void * src, void * dst, bsc_size_t size,
                      bsc_pid_t q );

TEST( bsc_bcast_qtree_qp, success() )
{
    int x1 = 0, x2 = 0;
    int y1 = 0, y2 = 0;
    bsp_begin( bsp_nprocs() );

    bsp_push_reg( &y1, sizeof(y1) );
    bsp_push_reg( &y2, sizeof(y2) );
    bsp_sync();

    if ( bsp_pid() == 2 % bsp_nprocs())
        x1 = 5;

    bsc_bcast_qtree( bsc_start, 2 % bsp_nprocs(), bsc_all,
            &x1, &y1, sizeof(x1), 3 );

    if ( bsp_pid() == 1 % bsp_nprocs())
        x2 = 2;

    bsc_bcast_1phase( bsc_start, 1 % bsp_nprocs(), bsc_all,
            &x2, &y2, sizeof(y2) );
    bsc_sync( bsc_flush );

    EXPECT_EQ( "%d", y1, 5 );
    EXPECT_EQ( "%d", y2, 2 );
    EXPECT_EQ( "%d", x1,  (bsp_pid() == 2 % bsp_nprocs()?5:0) );
    EXPECT_EQ( "%d", x2,  (bsp_pid() == 1 % bsp_nprocs()?2:0) );

    bsp_pop_reg( &y1 );
    bsp_pop_reg( &y2 );

    bsp_end();
}
