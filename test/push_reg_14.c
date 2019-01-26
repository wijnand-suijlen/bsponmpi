#include <bsp.h>
#include <stdlib.h>
#include "test.h"

TEST( push_reg_14 , success() ) {
    int x, y;
    bsp_begin( bsp_nprocs() );

    bsp_push_reg(NULL, 0);
    bsp_sync();

    bsp_pop_reg( NULL );

    bsp_sync();

    bsp_end();
}
