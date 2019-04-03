#include <bsp.h>
#include <bsc.h>
#include "test.h"

TEST( bsc_sync_stationary, abort("bsp_put: Remote address was just registered") )
{
    int x;
    bsp_begin( bsp_nprocs() );
    
    bsc_sync( 0 );

    bsp_push_reg( &x, sizeof(x) );
    bsc_sync( 0 );
    /* because no sync have been executed, the following will cause an error */
    bsp_put( (bsp_pid() + 1)%bsp_nprocs(), &x, &x, 0, sizeof(int) );

    bsp_end();
}
