#include <bsp.h>
#include "test.h"

TEST( set_tagsize_1, success() )
{
    bsp_size_t x;
    bsp_begin( bsp_nprocs() );

    x = 1;

    bsp_set_tagsize( &x );

    EXPECT_EQ("%d", x, 0 );

    bsp_set_tagsize( &x );

    EXPECT_EQ("%d", x, 1 );

    bsp_set_tagsize( &x );

    EXPECT_EQ("%d", x, 0 );

    bsp_sync();
    x = 5;
    bsp_set_tagsize( &x );

    EXPECT_EQ("%d", x, 1 );

    bsp_set_tagsize( &x );

    EXPECT_EQ("%d", x, 5 );

    bsp_end();
}
