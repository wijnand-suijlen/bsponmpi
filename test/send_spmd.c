#include <bsp.h>
#include "test.h"

#include <stdint.h>

TEST( send_spmd,
      abort("bsp_send: can only be called within SPMD section") )
{
    bsp_send( 0, NULL, NULL, 0 );
}
