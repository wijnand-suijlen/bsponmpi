#include <bsp.h>
#include "test.h"

#include <stdint.h>
#include <stdlib.h>

TEST( send_size_endtag,
      abort("bsp_send: Invalid payload size, becaues it is also used as end marker") )
{
    bsp_begin(bsp_nprocs());

    bsp_size_t size = -1;

    char * x = malloc(size);

    // Test does not matter if that size cannot even be allocated ?
    if (x != NULL)
    {
        bsp_send( (bsp_pid() + 1) % bsp_nprocs(), NULL, x, size );
    }

    bsp_sync();

    if (x != NULL)
    {
        free(x);
    }

    bsp_end();
}
