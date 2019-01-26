#include <bsp.h>
#include "test.h"

TEST( sync_2, abort("bsp_sync/bsp_end: Some processes have called bsp_sync, "
            "while others have called bsp_end instead") )
{
    bsp_begin( bsp_nprocs() );
    
    if (bsp_pid() != 0 ) bsp_sync();

    bsp_end();
}
