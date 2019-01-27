#include <bsp.h>
#include "test.h"

TEST( hpput_1, success() )
{
    int x, y;
    bsp_begin( bsp_nprocs() );

    x = y =  bsp_pid();

    bsp_push_reg( &x, sizeof(x) );

    bsp_sync();

    bsp_hpput( bsp_nprocs() - 1 - bsp_pid(), &y, &x, 0, sizeof(int) );

    bsp_pop_reg( &x );

    bsp_sync();

    EXPECT_EQ("%d", x, bsp_nprocs()-1-bsp_pid() );

    bsp_end();
}
