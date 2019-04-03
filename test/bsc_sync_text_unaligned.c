#include <bsp.h>
#include <bsc.h>
#include "test.h"

TEST( bsc_sync_text_unaligned, success() )
{
    bsc_step_t d ;
    bsp_begin( bsp_nprocs() );

    d = bsc_current();

    if ( bsp_pid() ) {
        bsc_sync( d );
    }

    EXPECT_EQ("%d", bsc_current(), bsp_pid()?1:0 );

    bsc_sync( bsc_flush );

    
    bsp_end();
}
