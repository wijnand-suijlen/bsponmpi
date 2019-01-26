#include <bsp.h>
#include <stdlib.h>
#include "test.h"

TEST( push_reg_2, success() ) {
    int * x;
    bsp_begin( bsp_nprocs() );

    x = calloc( bsp_pid()+1, sizeof(int) );

    bsp_push_reg( x, (1+bsp_pid())*sizeof(int) );

    bsp_sync();

    bsp_pop_reg( x );

    bsp_sync();
    free(x);

    bsp_end();
}
