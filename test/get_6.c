#include <bsp.h>
#include "test.h"

TEST( get_6, success() ) {
    int x;
    bsp_begin( bsp_nprocs() );

    x = bsp_pid();
    bsp_push_reg( &x, sizeof(x) );

    bsp_sync();

    bsp_get( bsp_pid(), &x, 0, &x, sizeof(x) );
    x = 8;

    bsp_sync();

    EXPECT_EQ( "%d", x, 8 );

    bsp_pop_reg( &x );

    bsp_sync();

    bsp_end();
}
