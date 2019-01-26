#include <bsp.h>
#include <stdlib.h>
#include "test.h"

TEST( pop_reg_5, success() ) {
    int x;
    bsp_begin( bsp_nprocs() );

    bsp_push_reg( NULL, 0 );
    bsp_push_reg( bsp_pid()?&x:NULL, 0);

    bsp_sync();
    bsp_pop_reg( NULL );

    bsp_sync();

    bsp_end();
}
