#include <bsp.h>
#include <stdlib.h>
#include "test.h"

TEST( push_reg_13, abort("bsp_sync/bsp_push_reg: Not all processes registered the same number of memory blocks") ) {
    int x;
    bsp_begin( bsp_nprocs() );

    if (bsp_pid()==0) 
        bsp_push_reg( &x, sizeof(int) );

    bsp_sync();

    bsp_pop_reg( &x );

    bsp_sync();

    bsp_end();
}
