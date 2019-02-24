#include <bsp.h>
#include "test.h"

TEST( move_queue_empty, abort("bsp_move: Message queue was empty") )
{
    bsp_begin(bsp_nprocs());

    int x;

    bsp_move( &x, sizeof(int) );

    bsp_end();
}
