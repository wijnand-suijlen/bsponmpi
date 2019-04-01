#include <bsc.h>
#include <bsp.h>
#include "test.h"

bsc_step_t bsc_bcast_2phase( bsc_step_t depends,
                      bsp_pid_t root, bsc_group_t group, 
                      bsc_coll_params_t * params, bsc_size_t n );

TEST( bsc_bcast_2phase_all, success() )
{
    int i, *xs = NULL, *ys = NULL;
    bsc_coll_params_t params;
    bsc_step_t finished;
    bsp_begin( bsp_nprocs() );

    if (bsp_pid() == 1) {
      xs = calloc( 1000, sizeof(int) );
      EXPECT_OP("%p", (void*) xs, !=, NULL);
      for ( i = 0; i < 1000; ++i) 
          xs[i] = 2 + 3*i;
    }

    ys = calloc( 1000, sizeof(int) );
    EXPECT_OP("%p", (void*) ys, !=, NULL);

    bsp_push_reg( ys, 1000*sizeof(int) );
    bsp_sync();

    
    params.src=xs;
    params.dst=ys;
    params.size=1000*sizeof(int);

    finished = bsc_bcast_2phase( bsc_start, 1, bsc_all, &params, 1);
    EXPECT_EQ("%d", finished-bsc_start, 2 );

    bsc_sync( bsc_flush );
    EXPECT_EQ("%d", 0, bsc_current() );

     for ( i = 0; i < 1000; ++i) 
          EXPECT_EQ( "%d", ys[i],  2 + 3*i );

    bsp_pop_reg( ys );
    free(ys);
    free(xs);

    bsp_end();
} 
