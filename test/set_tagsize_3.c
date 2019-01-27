#include <bsp.h>
#include "test.h"

TEST( set_tagsize_3, abort("bsp_sync/bsp_set_tagsize: Not all processes" 
            " set the tag size same number of times") )
{
    bsp_size_t x = 0;
    bsp_begin( bsp_nprocs() );


    if (bsp_pid())
      bsp_set_tagsize( &x );

    bsp_sync();
   
    bsp_end();
}
