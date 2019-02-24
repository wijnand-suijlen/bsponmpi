#include <bsp.h>
#include "test.h"

TEST( set_tagsize_null,
      abort("bsp_set_tagsize: NULL ptr as argument is now allowed") )
{
    bsp_begin( bsp_nprocs() );

    bsp_set_tagsize( NULL );

    bsp_sync();

    bsp_end();
}
