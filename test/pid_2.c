#include <bsp.h>
#include "test.h"

TEST( pid_2, abort("bsp_pid: can only be called within SPMD section") )
{
    bsp_pid();
}
