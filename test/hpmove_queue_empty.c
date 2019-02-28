#include <bsp.h>
#include "test.h"

TEST( hpmove_queue_empty,
      success() )
{
    void * tag;
    void * payload;
    bsp_size_t ret;

    bsp_begin(bsp_nprocs());

    ret = bsp_hpmove( &tag, &payload );

    EXPECT_EQ("%d", -1, ret);

    bsp_end();
}
