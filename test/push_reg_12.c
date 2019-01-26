#include <bsp.h>
#include <stdlib.h>
#include "test.h"

TEST( push_reg_12, abort("bsp_pop_reg: memory at address") ) {
    int x;
    bsp_begin( bsp_nprocs() );

    if (bsp_pid()) 
        bsp_push_reg( &x, sizeof(int) );

    bsp_pop_reg( &x );

    bsp_sync();

    bsp_end();
}
