#include <bsp.h>
#include "test.h"

TEST( put_1, abort("bsp_put: Writes 4 bytes beyond registered range")) {
    int x;
    bsp_begin( bsp_nprocs() );

    x = bsp_pid();

    bsp_push_reg( bsp_pid()?&x:NULL, bsp_pid()?sizeof(x):0 );

    bsp_sync();

    if (bsp_pid())
      bsp_put( bsp_nprocs() - 1 - bsp_pid(), &x, &x, 0, sizeof(int) );

    bsp_sync();

    bsp_pop_reg( &x );

    bsp_sync();

    bsp_end();
}
