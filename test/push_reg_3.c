#include <bsp.h>
#include <stdlib.h>
#include "test.h"

TEST( push_reg_3, success() ) {
    int * x, len, empty;
    bsp_begin( bsp_nprocs() );

    len = bsp_nprocs() - bsp_pid() - 1;
    x = len == 0 ? &empty : calloc( len, sizeof(int) );

    bsp_push_reg( x, len*sizeof(int) );

    bsp_sync();

    bsp_pop_reg( x );

    bsp_sync();
    if (len != 0 ) free(x);

    bsp_end();
}
