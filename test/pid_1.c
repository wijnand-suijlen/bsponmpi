#include <bsp.h>
#include "test.h"

TEST( pid_1, success() )
{
    bsp_begin( bsp_nprocs() );

    EXPECT_OP("%d", 0, <=, bsp_pid() );
    EXPECT_OP("%d", bsp_pid(), <, bsp_nprocs() ); 

    bsp_end();
}
