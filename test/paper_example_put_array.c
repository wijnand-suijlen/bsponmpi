#include <bsp.h>
#include "test.h"

void put_array( int * xs, int n ) {
    int i, dst_pid, dst_idx;
    int p = bsp_nprocs();
    int n_over_p = n/p;

    if ((n%p) != 0)
        bsp_abort("{put_array} n=%d not divisible by p=%d\n", n, p);
    bsp_push_reg( xs, n_over_p*sizeof(xs[0]));
    bsp_sync();

    for (i = 0; i < n_over_p; ++i) {
        dst_pid = xs[i] / n_over_p;
        dst_idx = xs[i] % n_over_p;
        bsp_put( dst_pid, &xs[i], xs, dst_idx*sizeof(xs[0]), sizeof(xs[0]));
    }
    bsp_sync();
    bsp_pop_reg(xs);
}

TEST( paper_example_put_array, success() )
{
    bsp_pid_t p, s;
    int xs[3];
    bsp_begin( bsp_nprocs() );

    p = bsp_nprocs();
    s = bsp_pid();

    xs[0] = (p - s)*3 - 1;
    xs[1] = (p - s)*3 - 2;
    xs[2] = (p - s)*3 - 3;

    put_array( xs, 3*p);

    EXPECT_EQ( "%d", xs[0], s*3 );
    EXPECT_EQ( "%d", xs[1], s*3 + 1);
    EXPECT_EQ( "%d", xs[2], s*3 + 2);
   
    bsp_end();
}
