#include <bsp.h>
#include "test.h"

TEST( hpmove_spmd, abort("bsp_hpmove: can only be called within SPMD section") )
{
    char * tag;
    char * payload;

    bsp_hpmove( &tag, &payload );
}
