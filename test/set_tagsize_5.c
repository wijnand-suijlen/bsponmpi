#include <bsp.h>
#include "test.h"

TEST( set_tagsize_5, success() )
{
    const int tag = 42, tagr_init = 24;
    int tagr = tagr_init;
    bsp_size_t status, tagsize = 1;

    bsp_begin( bsp_nprocs() );

    bsp_set_tagsize(&tagsize);
    bsp_send((bsp_pid() + 1) % bsp_nprocs(), &tag, NULL, 0);

    bsp_sync();
    bsp_set_tagsize(&tagsize);

    bsp_get_tag(&status, &tagr);

    /* Default tag size is zero, initial value of tagr should not change.
       Even if tag size is changed at the same time or even after the sync */
    EXPECT_EQ("%d", tagr, tagr_init);

    bsp_send((bsp_pid() + 1) % bsp_nprocs(), &tag, NULL, 0);
    bsp_sync();

    bsp_get_tag(&status, &tagr);

    EXPECT_EQ("%d", tagr, tag);

    bsp_end();
}
