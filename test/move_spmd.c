#include <bsp.h>
#include "test.h"

TEST( move_spmd, abort("bsp_move: can only be called within SPMD section") )
{
    bsp_move( NULL, 0 );
}
