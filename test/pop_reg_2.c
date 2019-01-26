#include <bsp.h>
#include <stdlib.h>
#include "test.h"

TEST( pop_reg_2, abort("bsp_pop_reg: memory at address") ) {
    int x;
    bsp_begin( bsp_nprocs() );

    bsp_push_reg( &x, sizeof(int) );

    bsp_pop_reg( &x );
    bsp_pop_reg( &x );

    bsp_sync();

    bsp_end();
}
