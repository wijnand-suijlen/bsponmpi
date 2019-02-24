#include <bsp.h>
#include "test.h"

TEST( move_nbytes_negative,
      abort("bsp_move: size argument may not be negative") )
{
    bsp_begin(bsp_nprocs());

    int x = bsp_pid();

    bsp_send( (bsp_pid() + 1)%bsp_nprocs(), NULL, &x, sizeof(int));

    bsp_sync();
    bsp_size_t s;
    bsp_move( &x, (bsp_size_t) -sizeof(x) );

    bsp_end();
}
