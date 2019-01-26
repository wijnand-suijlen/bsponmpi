#include <bsp.h>
#include "test.h"

TEST( nprocs_1, success() )
{
    bsp_pid_t p;

    p = bsp_nprocs();

    EXPECT_OP( "%d", p, >, 0 );
}
