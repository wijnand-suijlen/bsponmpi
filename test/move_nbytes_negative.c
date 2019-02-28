#include <bsp.h>
#include "test.h"

TEST( move_nbytes_negative,
      abort("bsp_move: size argument may not be negative") )
{
    int x;
    int int_size = (int) sizeof(x);
    bsp_begin(bsp_nprocs());

    x = bsp_pid();

    bsp_send( (bsp_pid() + 1)%bsp_nprocs(), NULL, &x, sizeof(int));

    bsp_sync();
    bsp_move( &x, 0 - int_size );

    bsp_end();
}
