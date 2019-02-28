#include <bsp.h>
#include "test.h"

TEST( hpmove_tagptr_null,
      abort("bsp_hpmove: pointer arugments may not be NULL") )
{
    void * payload;

    bsp_begin(bsp_nprocs());

    bsp_send((bsp_pid() + 1) % bsp_nprocs(), NULL, NULL, 0);

    bsp_sync();

    bsp_hpmove( NULL, &payload );

    bsp_end();
}
