#include <bsp.h>
#include <bsc.h>
#include "test.h"

TEST( bsc_sync_ring_comm_with_gaps, success())
{
    int y, gap;
    bsp_begin( bsp_nprocs() );
    
    bsp_push_reg( &y, sizeof(y) );
    bsp_sync();

    y = bsp_pid()?0:0xABCDEF;
    gap = 10;

    bsc_put( gap*bsp_pid(), (bsp_pid() + 1)%bsp_nprocs(), &y, &y, 0, sizeof(int) );
    bsc_sync( gap*(bsp_nprocs() - 1)+1);

    EXPECT_EQ("%d", (int) bsc_current(), gap*(bsp_nprocs()-1)+1  );
    EXPECT_EQ("%d", y, 0xABCDEF );

    bsc_sync( bsc_flush );

    EXPECT_EQ("%d", bsc_current(), 0 );

    bsp_end();
}
