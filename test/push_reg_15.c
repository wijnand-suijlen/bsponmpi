#include <bsp.h>
#include <stdlib.h>
#include "test.h"

TEST( push_reg_15 , success() ) {
    bsp_begin( bsp_nprocs() );

    bsp_push_reg(NULL, 0);
    bsp_pop_reg( NULL );

    bsp_sync();

    bsp_end();
}
