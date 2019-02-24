#include <bsp.h>
#include "test.h"

#include <stdint.h>

TEST( send_size_negative,
      abort("bsp_send: Payload size may not be negative") )
{
    bsp_begin(bsp_nprocs());

    int x = 42;

    bsp_send( (bsp_pid() + 1) % bsp_nprocs(), NULL, NULL, -1 );
    bsp_send( (bsp_pid() + 1) % bsp_nprocs(), NULL, &x, -1 );
    bsp_sync();

    bsp_end();
}
