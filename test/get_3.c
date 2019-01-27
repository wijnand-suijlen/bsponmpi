#include <bsp.h>
#include "test.h"

TEST( get_3, abort("bsp_get: Reads 4 bytes beyond registered range")) {
    int x;
    bsp_begin( bsp_nprocs() );

    x = bsp_pid();

    bsp_push_reg( bsp_pid()?&x:NULL, bsp_pid()?sizeof(x):0 );

    bsp_sync();

    if (bsp_pid())
      bsp_get( bsp_nprocs() - 1 - bsp_pid(), &x, 0, &x, sizeof(int) );

    bsp_sync();

    bsp_pop_reg( &x );

    bsp_sync();

    bsp_end();
}
