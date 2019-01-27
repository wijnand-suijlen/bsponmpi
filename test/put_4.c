#include <bsp.h>
#include "test.h"

/* test puts in a rooftile pattern and tests BSPonMPI's specific behaviour */

TEST( put_4, success() ) {
    int *x, i, P, s, y[2];
    bsp_begin( bsp_nprocs() );
    P = bsp_nprocs();
    s = bsp_pid();

    x = s ? &y[0] : calloc( P + 1, sizeof(int) );
    EXPECT( x );
    for ( i = 0; i < (s?0:P+1); ++i )
        x[i] = i*P + s;

    bsp_push_reg( x, s ? 0 : (P+1)*sizeof(int) );

    bsp_sync();

    bsp_pop_reg( x );
    y[0] = s*P + s;
    y[1] = (s+1)*P + s;
    bsp_put( 0, &y[0], &x[0], s*sizeof(int), 2*sizeof(int) );  

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
