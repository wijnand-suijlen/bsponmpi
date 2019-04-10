#include <bsp.h>
#include <bsc.h>
#include "test.h"

TEST( bsc_sync_ring_comm, success())
{
    int y;
    bsp_begin( bsp_nprocs() );
    
    bsp_push_reg( &y, sizeof(y) );
    bsp_sync();
    y = bsp_pid()?0:0xABCDEF;

    bsc_put( bsp_pid(), (bsp_pid() + 1)%bsp_nprocs(), &y, &y, 0, sizeof(int) );
    bsc_sync( bsp_nprocs() );

    EXPECT_EQ("%d", y, 0xABCDEF );
    EXPECT_EQ("%d", (int) bsc_current(), bsp_nprocs()  );

    bsc_sync( bsc_flush );

    EXPECT_EQ("%d", bsc_current(), 0 );

    bsp_end();
}
