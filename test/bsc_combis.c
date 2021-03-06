#include <bsp.h>
#include <bsc.h>
#include <string.h>

#include "test.h"

void prefix_sum( void * a, const void * a0, const void * xs, bsc_size_t size )
{
    int * result = a;
    const int * zero = a0;
    const int * nums = xs;
    int i, n = size/sizeof(int);
    int accu = *zero;

    for ( i = 0; i < n ; ++i ) {
        accu += nums[i];
        result[i] = accu;
    }
}

TEST( bsc_combis, success() )
{
    int i, n=6;
    int x[6] = { 888, 888, 888, 888, 888 };
    int y[6] = { 999, 999, 999, 999, 999 };

    int tmp[40] ; 
    int zero = 0;
    bsc_step_t ready = bsc_start, ready2;
    bsp_begin( bsp_nprocs() );
       
     for ( i = 0; (unsigned) i < sizeof(tmp)/sizeof(tmp[0]); ++i)
        tmp[i] = 9999999;
   
    EXPECT_OP( "%d", bsp_nprocs(), <, 20 ); /* because tmp space */

    bsp_push_reg( &tmp[0], (bsp_nprocs()+1)*sizeof(tmp[0]) );
    for (i = 0; i < n; ++i)
        bsp_push_reg( &x[i], sizeof(x[i]) );

    for ( i = 0; i < n; ++i)
        x[i] = n * bsp_pid() + i;

    bsp_sync();

    for (i = 0; i  < n; ++i ) 
        ready2 = bsc_bcast( ready, 0, bsc_all, &x[i], &x[i], sizeof(int) );

    bsc_scan( ready2, bsc_all, x, y, tmp, prefix_sum, &zero, n, sizeof(x[0]));
    
    bsc_sync( bsc_flush );

    for ( i = 0; i < n; ++i ) {
        EXPECT_EQ( "%d", x[i], i ); 
        EXPECT_EQ( "%d", y[i], (i*(i+1)/2) + (n*(n-1)/2*bsp_pid()) );
    }

    bsp_pop_reg( &tmp );

    for (i = 0; i < n; ++i)
        bsp_pop_reg( &x[i] );

    bsp_end();
}
