#include <bsp.h>
#include "test.h"

/* test gets in a pyramid pattern and tests BSPonMPI's specific behaviour */

TEST( get_5, success() ) {
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

    if (s == P-1 ) {
        for ( i = 0; i < P; ++i )
           bsp_get( i, &x[0], i*sizeof(int), &x[i], 2*(P-i)*sizeof(int));  
    }

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
