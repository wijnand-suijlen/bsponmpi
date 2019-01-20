#include <stdlib.h>
#include <assert.h>
#include <stdio.h>
#include <bsp.h>

int bsp_sum( int *xs, int nelem ) {
    int * local_sums, i, j, result, p;
    result = 0;
    p = bsp_nprocs();

    for (j = 0; j < nelem; ++j)
        result += xs[j];

    bsp_push_reg( &result, sizeof(result) );
    bsp_sync();

    local_sums = calloc(p, sizeof(local_sums[0]));
    if (local_sums == NULL)
        bsp_abort("bsp_sum: no memory for %d int", p);

    for (i = 0 ; i < p; ++i) 
        bsp_hpget( i, &result, 0, &local_sums[i], sizeof(int));
    
printf("[%d] e\n", bsp_pid() );
    bsp_sync();
printf("[%d] f\n", bsp_pid() );


    result = 0;
    for (i = 0 ; i < p; ++i)
        result += local_sums[i];

    bsp_pop_reg(&result);
    free(local_sums);
    return result;
}

int main( int argc, char * argv[]) 
{
    int p, s, y, i;
    int xs[6];
    (void) argc; (void) argv;
    bsp_begin( bsp_nprocs() );
    p = bsp_nprocs();
    s = bsp_pid();
    for ( i = 0; i < 6 ; ++i )
        xs[i] = i * s;

    y = bsp_sum( xs, 6 );

    assert( (0 + 1 + 2 + 3 + 4+ 5)*p*(p-1)/2 == y );

    bsp_end();
    return 0;
}
