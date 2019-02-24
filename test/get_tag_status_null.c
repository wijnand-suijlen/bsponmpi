#include <bsp.h>
#include "test.h"

TEST( get_tag_status_null,
      abort("bsp_get_tag: first arguments may not be NULL") )
{
    bsp_begin(bsp_nprocs());

    bsp_get_tag(NULL, NULL);

    bsp_end();
}
