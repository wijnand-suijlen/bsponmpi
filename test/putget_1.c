#include <bsp.h>
#include "test.h"

/* test puts and gets in a rooftile pattern and tests BSPonMPI's specific behaviour */

TEST( putget_1, success() ) {
    int *x, i, P, s, y, z[2];
    bsp_begin( bsp_nprocs() );
    P = bsp_nprocs();
    s = bsp_pid();

    x = s ? &y : calloc( P + 1, sizeof(int) );
    EXPECT( x );
    for ( i = 0; i < (s?0:P+1); ++i )
        x[i] = i*P + s;

    bsp_push_reg( x, s ? 0 : (P+1)*sizeof(int) );
    bsp_push_reg( &z[0], sizeof(z) );

    bsp_sync();

    bsp_pop_reg( &z[0] );
    bsp_pop_reg( x );

    z[0] = s*P + s;
    z[1] = (s+1)*P + s;

    if ( s % 2 == 0 ) {
        bsp_put( 0, &z[0], &x[0], s*sizeof(int), 2*sizeof(int) );  
    } 

    
    if ( s== 0 ) {
        for (i = 1; i < P; i+=2) {
           bsp_get( i, &z[0], 0, &x[i], 2*sizeof(int) );  
        }
    }


    bsp_sync();

    if (s == 0 ) {
        for ( i = 0; i < P/2; i++ ) {
            EXPECT_EQ( "%d", x[2*i]  , (2*i)*P+2*i );
            EXPECT_EQ( "%d", x[2*i+1], (2*i+1)*P+2*i );
        }
        if (P % 2 == 1)
            EXPECT_EQ( "%d", x[P], P*P+P-1 );

        free(x);
    }

    bsp_end();
}
