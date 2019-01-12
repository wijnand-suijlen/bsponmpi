#include <bsp.h>
#include <assert.h>

int main( int argc, char ** argv )
{
    (void) argc; (void) argv;

    bsp_pid_t p = bsp_nprocs();

    assert( p > 0 );

    return 0;
}
