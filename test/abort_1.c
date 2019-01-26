#include <bsp.h>
#include "test.h"


TEST( abort_1, abort("Dit is een test abort. Twee = 2") )
{
    bsp_abort("Dit is een test abort. Twee = %d\n", 2);
}
