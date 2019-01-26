#include <bsp.h>
#include <stdlib.h>
#include "test.h"

TEST( pop_reg_6, success() ) {
    int x, y;
    bsp_begin( bsp_nprocs() );

    bsp_push_reg( bsp_pid()?&x:NULL, 0);
    bsp_push_reg( NULL, 0 );

    bsp_sync();
    bsp_pop_reg( bsp_pid()?&x:NULL);

    bsp_sync();

    bsp_end();
}
