#include <bsp.h>
#include "test.h"

static int s_some_int = 0;

void spmd(void)
{
    bsp_pid_t s;
    bsp_begin( bsp_nprocs() );

    s = bsp_pid();
    EXPECT_IMPLIES( s==0, s_some_int == 5 );

    if (s == 0 )
        s_some_int = (int) s - 10;
    else
        s_some_int = 999;

    bsp_end();
}

TEST( init_1, success() )
{
    bsp_init( spmd, 0, NULL );

    s_some_int = 5;

    spmd();

    EXPECT_EQ("%d", s_some_int, -10 );
}
