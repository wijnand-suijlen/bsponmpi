#include <bsp.h>
#include "test.h"

#include <stdint.h>

TEST( qsize_accum_nbytes_null,
      abort("bsp_qsize: both arguments may not be NULL") )
{
    bsp_begin(bsp_nprocs());

    bsp_nprocs_t nmessages;
    bsp_qsize(&nmessages, NULL);

    bsp_end();
}
