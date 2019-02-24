#include <bsp.h>
#include "test.h"

TEST( move_payload_null,
      abort("bsp_move: payload may not be NULL if reception_nbytes is non-zero") )
{
    bsp_begin(bsp_nprocs());

    int x = bsp_pid();

    bsp_send( (bsp_pid() + 1)%bsp_nprocs(), NULL, &x, sizeof(int));

    bsp_sync();

    bsp_move( NULL, sizeof(x) );

    bsp_end();
}
