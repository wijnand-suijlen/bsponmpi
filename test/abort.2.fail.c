#include <bsp.h>
#include <assert.h>

int main( int argc, char ** argv )
{
    (void) argc; (void) argv;
    bsp_begin( bsp_nprocs() );
    
    bsp_abort("Dit is een test abort. Twee punt vijfenveertig = %g\n", 2.45);

    bsp_end();
    return 0;
}
