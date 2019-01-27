#include <bsp.h>
#include "test.h"

/* test gets in a rooftile pattern and tests BSPonMPI's specific behaviour */

TEST( get_4, success() ) {
    int *x, i, P, s, y[2];
    bsp_begin( bsp_nprocs() );
    P = bsp_nprocs();
    s = bsp_pid();

    x = s ? &y[0] : calloc( P + 1, sizeof(int) );
    EXPECT( x );
    for ( i = 0; i < (s?0:P+1); ++i )
        x[i] = i*P + s;

    bsp_push_reg( x, s ? 0 : (P+1)*sizeof(int) );
    bsp_push_reg( &y[0], sizeof(y) );

    bsp_sync();
    bsp_pop_reg( y );
    bsp_pop_reg( x );
    if ( s== 0 ) {
        for (i = 0; i < P; ++i) {
           bsp_get( i, &y[0], 0, &x[i], 2*sizeof(int) );  
        }
    }
    y[0] = s*P + s;
    y[1] = (s+1)*P + s;

    bsp_sync();

    if (s == 0 ) {
        for ( i = 0; i < P; ++i ) {
            EXPECT_EQ( "%d", x[i], i*P+i );
        }
        EXPECT_EQ( "%d", x[P], P*P+P-1 );
        free(x);
    }

    bsp_end();
}
