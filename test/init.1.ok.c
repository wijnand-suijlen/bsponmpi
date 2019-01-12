#include <bsp.h>
#include <assert.h>

static int s_some_int = 0;

void spmd(void)
{
    bsp_begin( bsp_nprocs() );

    bsp_pid_t s = bsp_pid();
    assert( !(s==0) || s_some_int == 5 );

    if (s == 0 )
        s_some_int = (int) s - 10;
    else
        s_some_int = 999;

    bsp_end();
}

int main( int argc, char ** argv )
{
    (void) argc; (void) argv;
    bsp_init( spmd, argc, argv );

    s_some_int = 5;

    spmd();

    assert( s_some_int == -10 );

    return 0;
}
