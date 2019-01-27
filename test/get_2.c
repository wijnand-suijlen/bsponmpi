#include <bsp.h>
#include "test.h"

TEST( gut_2, success() ) {
    int x ;
    bsp_begin( bsp_nprocs() );

    x = bsp_pid();

    bsp_push_reg( &x, sizeof(x) );
    bsp_sync();

    bsp_get( bsp_nprocs() - 1 - bsp_pid(), &x, 0,  &x, sizeof(int) );

    bsp_pop_reg( &x );

    bsp_sync();

    EXPECT_EQ( "%d", x, bsp_nprocs() - 1 - bsp_pid() );

    bsp_end();
}
