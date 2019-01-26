#include <bsp.h>
#include <stdlib.h>
#include "test.h"

TEST( pop_reg_3, abort("bsp_sync/bsp_pop_reg: Processes 0 and 1 did not pop the same memory registration") ) {
    int x, y;
    bsp_begin( bsp_nprocs() );

    bsp_push_reg( &x, sizeof(int) );
    bsp_push_reg( bsp_pid()?&x:&y, sizeof(int) );
    bsp_pop_reg( &x );

    bsp_sync();

    bsp_end();
}
