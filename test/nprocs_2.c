#include <bsp.h>
#include "test.h"

TEST( nprocs_2, success() )
{
    bsp_pid_t p, s;
    bsp_begin( bsp_nprocs() );
    
    p = bsp_nprocs();
    s = bsp_pid();

    EXPECT_OP( "%d", p, >, 0 );
    EXPECT_OP( "%d", s, <, p );
    EXPECT_OP( "%d", s, >=, 0);

    bsp_end();
}
