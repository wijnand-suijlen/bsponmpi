#include <bsp.h>
#include <stdio.h>

int main( int argc, char ** argv )
{
    (void) argc; (void) argv;
    bsp_begin( bsp_nprocs( ) );

    printf("Hi from %d/%d\n", bsp_pid(), bsp_nprocs() );

    bsp_end();
    return 0;
}

