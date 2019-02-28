#include <bsp.h>
#include "test.h"

TEST( default_tagsize_1, success() )
{
    const char tag = 42;
    const char tagr_init = 24;
    char tagr;
    bsp_size_t status;

    bsp_begin( bsp_nprocs() );

    tagr = tagr_init;

    bsp_send((bsp_pid() + 1) % bsp_nprocs(), &tag, NULL, 0);

    bsp_sync();

    bsp_get_tag(&status, &tagr);

    /* Default tag size is zero, initial value of tagr should not change */
    EXPECT_EQ("%d", tagr, tagr_init);

    bsp_end();
}
