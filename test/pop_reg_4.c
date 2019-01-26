#include <bsp.h>
#include <stdlib.h>
#include "test.h"

TEST( pop_reg_4, abort("bsp_sync/bsp_pop_reg: Tried to deregister NULL on all processes, but NULL was never registered") ) {
    bsp_begin( bsp_nprocs() );

    bsp_pop_reg( NULL );

    bsp_sync();

    bsp_end();
}
