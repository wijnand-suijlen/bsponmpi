#include <bsp.h>
#include "test.h"

TEST( hpmove_spmd, abort("bsp_hpmove: can only be called within SPMD section") )
{
    void * tag;
    void * payload;

    bsp_hpmove( &tag, &payload );
}
