#include <bsp.h>
#include <bsc.h>
#include <math.h>
#include <stdio.h>
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

bsc_step_t bsc_scan_qtree_single( bsc_step_t depends, bsc_group_t group, 
                    const void * src, void * dst, void * tmp_space,
                    bsc_reduce_t reducer, const void * zero,
                    bsc_size_t nmemb, bsc_size_t size, bsc_pid_t q );

TEST( bsc_scan_3tree, success() )
{
    int i, n1 = 4, n2 = 5, n3 = 6;
    int x[6];
    int y1[4] = { 9999, 9999, 9999, 9999};
    int y2[5] = { 9999, 9999, 9999, 9999, 9999};
    int y3[6] = { 9999, 9999, 9999, 9999, 9999, 9999};

    int tmp[20] ; 
    int zero = 0;
    bsc_step_t ready = bsc_start;
    bsp_begin( bsp_nprocs() );
       
     for ( i = 0; (unsigned) i < sizeof(tmp)/sizeof(tmp[0]); ++i)
        tmp[i] = 9999999;
   
    EXPECT_OP( "%d", bsp_nprocs(), <, 20 ); /* because tmp space */

    for ( i = 0; i < n1; ++i )
        x[i] = n1*bsp_pid()+i;

    bsp_push_reg( &tmp, (bsp_nprocs()+1)*sizeof(tmp) );
    bsp_sync();

    ready = bsc_scan_qtree_single( ready, bsc_all, 
            &x, &y1, &tmp, prefix_sum, &zero, n1, sizeof(x[0]), 
            3 );

    ready = bsc_sync( ready );

    for ( i = 0; i < n1; ++i ) {
        int j = n1*bsp_pid()+i;
        EXPECT_EQ( "%d", y1[i], j*(j+1)/2 );
    }

    for ( i = 0; i < n2; ++i ) {
        x[i] = n2*bsp_pid()+i;
        y2[i] = 9999;
    }
    for ( i = 0; (unsigned) i < sizeof(tmp)/sizeof(tmp[0]); ++i)
        tmp[i] = 9999999;
  
    ready = bsc_scan_qtree_single( ready, bsc_all, 
            &x, &y2, &tmp, prefix_sum, &zero, n2, sizeof(x[0]), 
            2 );

    ready = bsc_sync( ready );

    for ( i = 0; i < n2; ++i ) {
        int j = n2*bsp_pid()+i;
        EXPECT_EQ( "%d", y2[i], j*(j+1)/2 );
    }

    for ( i = 0; i < n3; ++i ) {
        x[i] = n3*bsp_pid()+i;
        y3[i] = 9999;
    }
    for ( i = 0; (unsigned) i < sizeof(tmp)/sizeof(tmp[0]); ++i)
        tmp[i] = 9999999;
 
    ready = bsc_scan_qtree_single( ready, bsc_all, 
            &x, &y3, &tmp, prefix_sum, &zero, n3, sizeof(x[0]), 
            ceil(sqrt(bsp_nprocs())) );

    bsc_sync( ready );

    for ( i = 0; i < n3; ++i ) {
        int j = n3*bsp_pid()+i;
        EXPECT_EQ( "%d", y3[i], j*(j+1)/2 );
    }


    bsp_pop_reg( &tmp );

    bsp_end();
}
