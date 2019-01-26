#include <bsp.h>
#include "test.h"

TEST( put_1, abort("bsp_put: Remote address may not be used because it was registered with NULL")) {
    int x;
    bsp_begin( bsp_nprocs() );

    x = bsp_pid();

    bsp_push_reg( bsp_pid()?&x:NULL, sizeof(x) );

    bsp_sync();

    if (bsp_pid())
      bsp_put( bsp_nprocs() - 1 - bsp_pid(), &x, &x, 0, sizeof(int) );

    bsp_sync();

    bsp_pop_reg( &x );

    bsp_sync();

    bsp_end();
}
