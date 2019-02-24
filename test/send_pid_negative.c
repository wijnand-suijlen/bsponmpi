#include <bsp.h>
#include "test.h"

#include <stdint.h>

TEST( send_pid_negative,
      abort("bsp_send: The source process ID does not exist") )
{
    bsp_begin(bsp_nprocs());

    bsp_send( -1, NULL, NULL, 0 );
    bsp_sync();

    bsp_end();
}
