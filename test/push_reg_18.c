#include <bsp.h>
#include "test.h"

TEST( push_reg_18, success() ) {
    int x[100], y[100], i;
    bsp_begin( bsp_nprocs() );

    for ( i = 0; i < 100; ++i ) {
        bsp_push_reg( &x[i], sizeof(x[i]) );
        bsp_push_reg( &y[i], sizeof(y[i]) );
    }

    for ( i = 0; i < 100; ++i ) {
        bsp_pop_reg( &x[i] );
        bsp_sync();
    }

    for ( i = 0; i < 100; ++i )
        bsp_pop_reg( &y[i] );

    bsp_sync();

    bsp_end();
}
