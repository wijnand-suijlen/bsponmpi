#include <bsp.h>
#include "test.h"

TEST( push_reg_17, success() ) {
    int x, y;
    bsp_begin( bsp_nprocs() );

    bsp_push_reg( &x, sizeof(x) );
    bsp_push_reg( &y, sizeof(y) );

    bsp_pop_reg( &x );
    bsp_sync();

    bsp_pop_reg( &y );

    bsp_sync();

    bsp_end();
}
