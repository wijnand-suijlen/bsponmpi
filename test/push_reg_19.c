#include <bsp.h>
#include "test.h"

TEST( push_reg_19, success() ) {
    int x[100], y, i, j, k;
    bsp_begin( bsp_nprocs() );

    for ( i = 0; i < 45; ++i ) {
        bsp_push_reg( &y, sizeof(y) );
        bsp_push_reg( &x[i], sizeof(x[i]) );
    }

    k = 0;
    for ( i = 0; i < 10; ++i ) {
        for (j = 0; j < i; ++j ) {
            bsp_pop_reg( &x[k++] );
        }

        bsp_pop_reg( &y );
        bsp_sync();
    }

    for ( i = 0; i < 35; ++i)
        bsp_pop_reg( &y );

    bsp_sync();

    bsp_end();
}
