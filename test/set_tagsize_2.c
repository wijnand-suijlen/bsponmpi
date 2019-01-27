#include <bsp.h>
#include "test.h"

TEST( set_tagsize_2, abort("bsp_sync/bsp_set_tagsize: Not all processes "
           "set the tag size same number of times") )
{
    bsp_size_t x;
    bsp_begin( bsp_nprocs() );

    x = 1;

    if (bsp_pid())
      bsp_set_tagsize( &x );

    bsp_sync();
   
    bsp_end();
}
