#include <bsp.h>
#include <bsc.h>
#include <math.h>
#include "test.h"

void sum( void * a, const void * a0, const void * xs, bsc_size_t size )
{
    int * result = a;
    const int * zero = a0;
    const int * nums = xs;
    int i, n = size/sizeof(int);
    *result = *zero;
    for ( i = 0; i < n ; ++i )
        *result += nums[i];
}

bsc_step_t bsc_allreduce_qtree_single( bsc_step_t depends, bsc_group_t group, 
                    const void * src, void * dst, void * tmp_space,
                    bsc_reduce_t reducer, const void * zero,
                    bsc_size_t nmemb, bsc_size_t size, bsc_pid_t q );

TEST( bsc_allreduce_qtree_root_p, success() )
{
    int x = 0;
    int y = 0;
    int tmp[20] ;
    int zero = 0;
    bsc_step_t ready;
    bsp_begin( bsp_nprocs() );
    
    EXPECT_OP( "%d", bsp_nprocs(), <, 20 ); /* because tmp space */

    bsp_push_reg( &tmp, (bsp_nprocs()+1)*sizeof(tmp) );
    bsp_sync();

    x = bsp_pid()+1;

    ready = bsc_allreduce_qtree_single( bsc_start, bsc_all, 
            &x, &y, &tmp, sum, &zero, 1, sizeof(x), 
            ceil(sqrt(bsp_nprocs())) );

    bsc_sync( ready );

    EXPECT_EQ( "%d", y, (bsp_nprocs())*(bsp_nprocs()+1)/2 );

    bsp_pop_reg( &tmp );

    bsp_end();
}
