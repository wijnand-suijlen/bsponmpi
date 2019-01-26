#include <bsp.h>
#include "test.h"

TEST( put_2, success() ) {
    int x ;
    bsp_begin( bsp_nprocs() );

    x = bsp_pid();

    bsp_push_reg( &x, sizeof(x) );
    bsp_sync();

    bsp_put( bsp_nprocs() - 1 - bsp_pid(), &x, &x, 0, sizeof(int) );

    bsp_pop_reg( &x );

    bsp_sync();

    EXPECT_EQ( "%d", x, bsp_nprocs() - 1 - bsp_pid() );

    bsp_end();
}
