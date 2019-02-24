#include <bsp.h>
#include "test.h"

TEST( hpmove_payload_null,
      abort("bsp_hpmove: pointer arugments may not be NULL") )
{
    bsp_begin(bsp_nprocs());

    char * tag;

    bsp_send((bsp_pid() + 1) % bsp_nprocs(), NULL, NULL, 0);

    bsp_sync();

    bsp_hpmove( (void **) &tag, NULL );

    bsp_end();
}
