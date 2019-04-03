#include <bsp.h>
#include <bsc.h>
#include "test.h"

TEST( bsc_sync_ring_double_comm_with_gaps, success())
{
    int y, z, gap;
    bsp_begin( bsp_nprocs() );
    
    bsp_push_reg( &y, sizeof(y) );
    bsp_push_reg( &z, sizeof(z) );
    bsp_sync();
    y = bsp_pid()?0:0xABCDEF;
    z = 0;
    gap = 10;

    /* ring forward */
    bsc_put( gap*bsp_pid(), (bsp_pid() + 1)%bsp_nprocs(), &y, &y, 0, sizeof(int) );
    /* and a ring backward */
    bsc_put( gap*(2*bsp_nprocs()-bsp_pid()-1),
            (bsp_pid() + bsp_nprocs() - 1)%bsp_nprocs(), &y, &z, 0, sizeof(int) );
    bsc_sync( gap*(2*bsp_nprocs() - 1));

    EXPECT_EQ("%d", y, 0xABCDEF );
    EXPECT_EQ("%d", z, 0xABCDEF );
    EXPECT_EQ("%d", (int) bsc_current(), gap*(2*bsp_nprocs()-1)+1  );

    bsc_sync( bsc_flush );

    EXPECT_EQ("%d", bsc_current(), 0 );

    bsp_end();
}
