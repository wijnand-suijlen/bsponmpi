#include <bsp.h>
#include "test.h"

TEST( begin_1, abort("bsp_begin: May be called only once") )
{
    bsp_begin( bsp_nprocs() );
    
    bsp_begin( 2 );

    bsp_end();
}
