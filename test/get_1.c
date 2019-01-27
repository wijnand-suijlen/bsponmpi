#include <bsp.h>
#include "test.h"

TEST( get_1, abort("bsp_get: Remote address was just registered")) {
    int x;
    bsp_begin( bsp_nprocs() );

    x = bsp_pid();

    bsp_push_reg( &x, sizeof(x) );

    bsp_get( bsp_nprocs() - 1 - bsp_pid(), &x, 0, &x, sizeof(int) );

    bsp_sync();

    bsp_pop_reg( &x );

    bsp_sync();

    bsp_end();
}
