#include <bsp.h>
#include <assert.h>

int main( int argc, char ** argv )
{
    bsp_pid_t p, s;
    (void) argc; (void) argv;
    bsp_begin( bsp_nprocs() );
    
    p = bsp_nprocs();
    s = bsp_pid();

    assert( p > 0 );
    assert( s < p );
    assert( s >= 0);

    bsp_end();
    return 0;
}
