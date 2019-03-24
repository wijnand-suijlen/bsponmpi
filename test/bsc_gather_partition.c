#include <bsp.h>
#include <bsc.h>
#include "test.h"

TEST( bsc_gather_partition, success() )
{
    bsc_group_t group;
    int x[3] = { 0, 0, 0 };
    int y;

    bsp_begin( bsp_nprocs() );
    y = bsp_pid()+1;
    EXPECT_EQ( "%d", bsp_nprocs(), 5 ); /* this test must run on 5 processes*/

    group = bsc_group_create_partition( bsp_pid() / 3 );
    bsp_push_reg( &y, sizeof(y) );
    bsp_sync();

    bsc_gather( bsc_start, (bsp_pid()/3)*3+1, group, &y, x, sizeof(x[0]));
    bsc_sync( bsc_flush );

    if ( bsp_pid() == 1 ) {
        EXPECT_EQ( "%d", x[0], 1 );
        EXPECT_EQ( "%d", x[1], 2 );
        EXPECT_EQ( "%d", x[2], 3 );
    }
    else if (bsp_pid() == 4 ) {
        EXPECT_EQ( "%d", x[0], 4 );
        EXPECT_EQ( "%d", x[1], 5 );
        EXPECT_EQ( "%d", x[2], 0 );
    }
    else {
        EXPECT_EQ( "%d", x[0], 0 );
        EXPECT_EQ( "%d", x[1], 0 );
        EXPECT_EQ( "%d", x[2], 0 );
    }

    bsp_pop_reg( &y );

    bsp_end();
}
