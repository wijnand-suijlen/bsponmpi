#include <bsp.h>
#include "test.h"

TEST( get_tag_tag_null,
      abort("bsp_get_tag: tag must point to a block of available "
            "memory of at least 1 bytes long\n") )
{
    bsp_size_t status, tagsize;
    char tag;

    bsp_begin(bsp_nprocs());

    tagsize = sizeof(char);
    tag = 42;

    bsp_set_tagsize(&tagsize);

    bsp_sync();

    bsp_send((bsp_pid() + 1) % bsp_nprocs(), &tag, NULL, 0);

    bsp_sync();

    bsp_get_tag(&status, NULL);

    bsp_end();
}
