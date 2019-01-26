#include <bsp.h>
#include "test.h"

TEST( push_reg_1, success() ) {
    int x;
    bsp_begin( bsp_nprocs() );

    x = bsp_pid();

    bsp_push_reg( &x, sizeof(x) );

    bsp_sync();

    bsp_pop_reg( &x );

    bsp_sync();

    bsp_end();
}
