#include <bsp.h>
#include <assert.h>

int reverse( int x )
{
    bsp_push_reg( &x, sizeof(x) );
    bsp_sync();

    bsp_put( bsp_nprocs() - bsp_pid() - 1, &x, &x, 0, sizeof(x) );
    bsp_pop_reg(&x);
    bsp_sync();
    return x;
}

int main( int argc, char ** argv )
{
    int x;
    (void) argc; (void) argv;
    bsp_begin( bsp_nprocs() );
    
    x = reverse( bsp_pid() );
    assert( bsp_nprocs() - bsp_pid() - 1 == x );
   
    bsp_end();
    return 0;
}
