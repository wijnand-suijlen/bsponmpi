#include <bsp.h>
#include "test.h"

TEST( move_1, abort("bsp_move: Message queue was empty") )
{
    bsp_begin(bsp_nprocs());

    bsp_send( (bsp_pid() + 1)%bsp_nprocs(), NULL, NULL, 0);

    bsp_sync();

    bsp_sync();

    bsp_move( NULL, 0 );

    bsp_end();
}
