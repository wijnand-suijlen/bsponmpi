#include <bsp.h>
#include "test.h"

TEST( set_tagsize_4, abort("bsp_sync/bsp_set_tagsize: Not all processes "
           "announced the same tag size") )
{
    bsp_size_t x;
    bsp_begin( bsp_nprocs() );

    x = bsp_pid();
    bsp_set_tagsize( &x );

    bsp_sync();
   
    bsp_end();
}
