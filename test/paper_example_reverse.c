#include <bsp.h>
#include "test.h"

int reverse( int x )
{
    bsp_push_reg( &x, sizeof(x) );
    bsp_sync();

    bsp_put( bsp_nprocs() - bsp_pid() - 1, &x, &x, 0, sizeof(x) );
    bsp_pop_reg(&x);
    bsp_sync();
    return x;
}

TEST( paper_example_reverse, success() )
{
    int x;
    bsp_begin( bsp_nprocs() );
    
    x = reverse( bsp_pid() );
    EXPECT_EQ( "%d", x,  bsp_nprocs() - bsp_pid() - 1 );
   
    bsp_end();
}
