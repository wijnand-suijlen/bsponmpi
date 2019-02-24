#include <bsp.h>
#include "test.h"

TEST( set_tagsize_negative_size,
      abort("bsp_set_tagsize: Tag size may not be negative") )
{
    bsp_size_t x;
    bsp_begin( bsp_nprocs() );

    x = -1;
    bsp_set_tagsize( &x );

    bsp_sync();

    bsp_end();
}
