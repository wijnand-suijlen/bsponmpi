#include <bsp.h>
#include "test.h"

TEST( put_6, success() ) {
    int x;
    bsp_begin( bsp_nprocs() );

    x = bsp_pid();
    bsp_push_reg( &x, sizeof(x) );

    bsp_sync();

    bsp_put( bsp_pid(), &x, &x, 0, sizeof(x) );
    x = 8;

    bsp_sync();

    EXPECT_EQ( "%d", x, bsp_pid() );

    bsp_pop_reg( &x );

    bsp_sync();

    bsp_end();
}
