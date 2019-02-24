#include <bsp.h>
#include "test.h"

TEST( hpmove_queue_empty,
      success() )
{
    bsp_begin(bsp_nprocs());

    char * tag;
    char * payload;
    bsp_size_t ret;

    ret = bsp_hpmove( (void **) &tag, (void **) &payload );

    EXPECT_EQ("%d", -1, ret);

    bsp_end();
}
