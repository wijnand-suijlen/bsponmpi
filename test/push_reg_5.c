#include <bsp.h>
#include <stdlib.h>
#include "test.h"

TEST( push_reg_5, success() ) {
    int x;
    bsp_begin( bsp_nprocs() );

    if ( bsp_pid() == bsp_nprocs() - 1)
        bsp_push_reg( &x, sizeof(int) );
    else
        bsp_push_reg( NULL, 0 );

    bsp_sync();

    if ( bsp_pid() == bsp_nprocs() - 1)
        bsp_pop_reg( &x );
    else
        bsp_pop_reg( NULL );

    bsp_sync();

    bsp_end();
}
