#include <bsp.h>
#include "test.h"

TEST( set_tagsize_spmd,
      abort("bsp_set_tagsize: can only be called within SPMD section") )
{
    bsp_size_t x;

    x = 1;
    bsp_set_tagsize( &x );
}
