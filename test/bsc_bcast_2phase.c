#include <bsc.h>
#include <bsp.h>
#include "test.h"

bsc_step_t bsc_bcast_2phase( bsc_step_t depends,
                      bsp_pid_t root, bsc_group_t group, 
                      const void * src, void * dst, bsc_size_t size );

TEST( bsc_bcast_2phase_all, success() )
{
    int i, *xs = NULL, *ys = NULL;
    bsp_begin( bsp_nprocs() );

    if (bsp_pid() == 1) {
      xs = calloc( 1000, sizeof(int) );
      for ( i = 0; i < 1000; ++i) 
          xs[i] = 2 + 3*i;
    }

    ys = calloc( 1000, sizeof(int) );

    bsp_push_reg( ys, 1000*sizeof(int) );
    bsp_sync();

    bsc_bcast_2phase( bsc_start, 1, bsc_all,
            xs, ys, 1000*sizeof(int));

    bsc_sync( bsc_flush );

     for ( i = 0; i < 1000; ++i) 
          EXPECT_EQ( "%d", ys[i],  2 + 3*i );

    bsp_pop_reg( ys );
    free(ys);
    free(xs);

    bsp_end();
} 
