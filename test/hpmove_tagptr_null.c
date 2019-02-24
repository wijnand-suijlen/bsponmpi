#include <bsp.h>
#include "test.h"

TEST( hpmove_tagptr_null,
      abort("bsp_hpmove: pointer arugments may not be NULL") )
{
    bsp_begin(bsp_nprocs());

    char * payload;

    bsp_send((bsp_pid() + 1) % bsp_nprocs(), NULL, NULL, 0);

    bsp_sync();

    bsp_hpmove( NULL, (void **) &payload );

    bsp_end();
}
