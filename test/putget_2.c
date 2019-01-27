#include <bsp.h>
#include "test.h"

/* puts and gets in various numbers */

TEST( putget_2, success() ) {
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

    for ( i = 0; i < P; ++i ) {
        if ( s % 2 == 0 && s > 1 )
           bsp_put( i, &x[s], &x[0], s*sizeof(int), 2*(P-s)*sizeof(int) );  
        if ( i % 2 == 1 )
           bsp_get( i, &x[0], i*sizeof(int), &x[i], 2*(P-i)*sizeof(int));
    }

    /*** example write pattern for P == 5 
     *
     * s0: get( 1;    x[1..9]  )
     *     get( 3;    x[3..7]  ) 
     *
     * s1: get( 1;    x[1..9]  )
     *     get( 3;    x[3..7]  ) 
     *
     * s2: put( 0..P; x[2..8] )
     *     get( 1;    x[1..9]  )
     *     get( 3;    x[3..7]  ) 
     *
     * s3: get( 1;    x[1..9]  )
     *     get( 3;    x[3..7]  ) 
     *
     * s4: put( 0..P; x[4..6] )
     *     get( 1;    x[1..9]  )
     *     get( 3;    x[3..7]  ) 
     *
     *     
     */

    bsp_sync();

    bsp_pop_reg( x );
    EXPECT_EQ("%d", x[0], s );
    EXPECT_EQ("%d", x[2*P-1], (2*P-1)*P + s );
    EXPECT_EQ( "%d", x[1], 1*P+1 );
    EXPECT_EQ( "%d", x[2*P-2], (2*P-2)*P+1 );

    for ( i = 2; i < P; i+=2 ) {
        EXPECT_EQ( "%d", x[i], i*P+i);
        EXPECT_EQ( "%d", x[i+1], (i+1)*P+i );

        EXPECT_EQ( "%d", x[2*P-i-1], (2*P-i-1)*P+i);
        EXPECT_EQ( "%d", x[2*P-i-2], (2*P-i-2)*P+i );
    }

    free(x);
    bsp_end();
}
