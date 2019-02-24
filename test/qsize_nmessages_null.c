#include <bsp.h>
#include "test.h"

#include <stdint.h>

TEST( qsize_nmessages_null,
      abort("bsp_qsize: both arguments may not be NULL") )
{
    bsp_begin(bsp_nprocs());

    bsp_size_t nbytes;
    bsp_qsize(NULL, &nbytes);

    bsp_end();
}
