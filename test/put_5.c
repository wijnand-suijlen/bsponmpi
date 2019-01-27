#include <bsp.h>
#include "test.h"

/* test puts in a pyramid pattern and tests BSPonMPI's specific behaviour */

TEST( put_5, success() ) {
    int *x, i, P, s;
    bsp_begin( bsp_nprocs() );
    P = bsp_nprocs();
    s = bsp_pid();

    x = calloc( 2*P, sizeof(int) );
    EXPECT( x );
    for ( i = 0; i < 2*P; ++i )
        x[i] = i*P + s;

    bsp_push_reg( x, (2*P)*sizeof(int) );

    bsp_sync();

    bsp_put( P-1, &x[s], &x[0], s*sizeof(int), 2*(P-s)*sizeof(int) );  

    bsp_sync();

    bsp_pop_reg( x );
    if (s == P-1 ) {
        for ( i = 0; i < P; ++i ) {
            EXPECT_EQ( "%d", x[i], i*P+i );
            EXPECT_EQ( "%d", x[2*P-i-1], (2*P-i-1)*P + i );
        }
    }

    free(x);
    bsp_end();
}
