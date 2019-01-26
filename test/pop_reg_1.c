#include <bsp.h>
#include <stdlib.h>
#include "test.h"

TEST( pop_reg_1, abort("bsp_sync/bsp_pop_reg: Not all processes popped the same number") ) {
    int x;
    bsp_begin( bsp_nprocs() );

    bsp_push_reg( &x, sizeof(int) );

    bsp_sync();

    if (bsp_pid())
      bsp_pop_reg( &x );

    bsp_sync();

    bsp_end();
}
