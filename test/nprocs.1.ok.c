#include <bsp.h>
#include <assert.h>

int main( int argc, char ** argv )
{
    bsp_pid_t p;
    (void) argc; (void) argv;

    p = bsp_nprocs();

    assert( p > 0 );

    return 0;
}
