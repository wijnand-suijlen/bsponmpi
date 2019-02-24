#include <bsp.h>
#include "test.h"

#include <stdint.h>

TEST( qsize_spmd,
      abort("bsp_qsize: can only be called within SPMD section") )
{
    bsp_nprocs_t nmessages;
    bsp_size_t nbytes;
    bsp_qsize(&nmessages, &nbytes);
}
