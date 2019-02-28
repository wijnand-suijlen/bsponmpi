#include <bsp.h>
#include "test.h"

TEST( move_queue_empty, abort("bsp_move: Message queue was empty") )
{
    int x;
    bsp_begin(bsp_nprocs());

    bsp_move( &x, sizeof(int) );

    bsp_end();
}
