#include <bsp.h>
#include <stdlib.h>
#include "test.h"

TEST( push_reg_16 , abort("bsp_sync/bsp_push_reg: Not all processes registered the same number of memory blocks") ) {
    int x, y;
    bsp_begin( bsp_nprocs() );

    if (bsp_pid())
      bsp_push_reg(NULL, 0);

    bsp_pop_reg( NULL );

    bsp_sync();

    bsp_end();
}
