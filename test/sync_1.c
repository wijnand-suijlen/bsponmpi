#include <bsp.h>
#include "test.h"

TEST( sync_1, success() )
{
    bsp_begin( bsp_nprocs() );
    
    bsp_sync();

    bsp_end();
}
