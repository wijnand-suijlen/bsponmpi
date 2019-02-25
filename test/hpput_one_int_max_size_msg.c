#include <bsp.h>
#include <limits.h>
#include <string.h>
#include "test.h"

TEST( hpput_1, success() )
{
    int i;
    char * buf = NULL;
    bsp_begin( bsp_nprocs() );

    buf = malloc( INT_MAX );
    
    for ( i = 0; i < INT_MAX; ++i ) {
        buf[i] = (char) (bsp_pid() + bsp_nprocs() * i );
    }

    bsp_push_reg( buf, INT_MAX );

    bsp_sync();

    if (bsp_pid() == 0)
      bsp_hpput( ( 1 + bsp_pid() ) % bsp_nprocs() , buf, buf, 0, INT_MAX );

    bsp_pop_reg( buf );

    bsp_sync();

    if ( 1 % bsp_nprocs() == bsp_pid() ) {
        for ( i = 0; i < INT_MAX; ++i ) {
                EXPECT_EQ( "%d", (int) buf[i], (int) (char) (bsp_nprocs() * i ) );
        }
    } else {
        for ( i = 0; i < INT_MAX; ++i ) {
                EXPECT_EQ( "%d", (int) buf[i], (int) (char) (bsp_pid() + bsp_nprocs() * i ) );
        }
    }

    bsp_end();
}
