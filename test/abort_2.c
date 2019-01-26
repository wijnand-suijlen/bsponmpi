#include <bsp.h>
#include "test.h"

TEST( abort_2, abort("Dit is een test abort. Twee punt vijfenveertig = 2.45"))
{
    bsp_begin( bsp_nprocs() );
    bsp_abort("Dit is een test abort. Twee punt vijfenveertig = %g\n", 2.45);
    bsp_end();
}
