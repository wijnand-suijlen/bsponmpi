#include <bsp.h>
#include <stdlib.h>
#include "test.h"

TEST( pop_reg_4, abort("bsp_sync/bsp_pop_reg: Tried to deregister NULL on all processes, but could not find a matching regisration") ) {
    int x, y;
    bsp_begin( bsp_nprocs() );

    bsp_push_reg( bsp_pid()?NULL:&x, 0);
    bsp_push_reg( !bsp_pid()?NULL:&x, 0);
    bsp_pop_reg( NULL );

    bsp_sync();

    bsp_end();
}
