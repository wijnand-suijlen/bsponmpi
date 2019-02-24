#include <bsp.h>
#include "test.h"

TEST( get_tag_spmd, abort("bsp_get_tag: can only be called within SPMD section") )
{
    bsp_size_t status;

    bsp_get_tag(&status, NULL);
}
