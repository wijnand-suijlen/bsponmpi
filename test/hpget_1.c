#include <bsp.h>
#include "test.h"

TEST( hpget_1, success() )
{
    int x, y;
    bsp_begin( bsp_nprocs() );

    x = y =  bsp_pid();

    bsp_push_reg( &x, sizeof(x) );

    bsp_sync();

    bsp_hpget( bsp_nprocs() - 1 - bsp_pid(), &x, 0, &y, sizeof(int) );

    bsp_pop_reg( &x );

    bsp_sync();

    EXPECT_EQ("%d", y, bsp_nprocs()-1-bsp_pid() );

    bsp_end();
}
