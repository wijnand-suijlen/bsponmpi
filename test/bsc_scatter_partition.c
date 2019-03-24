#include <bsp.h>
#include <bsc.h>
#include "test.h"

TEST( bsc_scatter_partition, success() )
{
    bsc_group_t group;
    const int x[3] = { 1, 2, 3 };
    int y = 0;

    bsp_begin( bsp_nprocs() );

    group = bsc_group_create_partition( bsp_pid() / 3 );
    bsp_push_reg( &y, sizeof(y) );
    bsp_sync();

    bsc_scatter( bsc_start, (bsp_pid()/3)*3, group, x, &y, sizeof(x[0]));
    bsc_sync( bsc_flush );


    EXPECT_EQ( "%d", y, x[bsp_pid() % 3] );

    bsp_pop_reg( &y );

    bsp_end();
}
