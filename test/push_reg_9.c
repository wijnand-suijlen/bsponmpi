#include <bsp.h>
#include <stdlib.h>
#include "test.h"

TEST( push_reg_9, success() ) {
    int x;
    bsp_begin( bsp_nprocs() );

    bsp_push_reg( &x, sizeof(int) );
    bsp_sync();
    bsp_push_reg( &x, sizeof(int) );

    bsp_sync();

    bsp_pop_reg( &x );

    bsp_sync();

    bsp_pop_reg( &x );

    bsp_sync();

    bsp_end();
}
