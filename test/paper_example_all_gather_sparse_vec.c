#include <bsp.h>
#include <stdlib.h>
#include "test.h"
#include <stdio.h>

int all_gather_sparse_vec( float * dense, int n_over_p,
                            float ** sparse_out,
                            int ** sparse_ivec_out) {
   int global_idx, i, j, tag_size, p, 
       nonzeros, nonzeros_size, status;
   int *sparse_ivec = NULL;
   float *sparse = NULL;

   p = bsp_nprocs();
   tag_size = sizeof(int);
   bsp_set_tagsize(&tag_size);
   bsp_sync();

   for ( i = 0; i < n_over_p; ++i ) {
     if ( dense[i] != 0.0 ) {
       global_idx = n_over_p * bsp_pid() + i;
       for (j = 0; j < p; ++j)
         bsp_send(j, &global_idx, &dense[i], sizeof(float) );
     }
   }
   bsp_sync();

   bsp_qsize( &nonzeros, &nonzeros_size );
   if (nonzeros > 0 ) {
     sparse      = calloc( nonzeros, sizeof(float) );
     sparse_ivec = calloc( nonzeros, sizeof(int) );
     if (sparse == NULL || sparse_ivec == NULL)
       bsp_abort("Unable to allocate memory\n");
 
     for ( i = 0 ; i < nonzeros; ++i) {
       bsp_get_tag( &status, &sparse_ivec[i] );
       EXPECT_EQ( "%d", status, (int) sizeof(float));
       bsp_move( &sparse[i], sizeof(float) );
     }
   }
   bsp_set_tagsize(&tag_size);
   *sparse_out = sparse;
   *sparse_ivec_out = sparse_ivec;
   return nonzeros;
}

TEST( paper_all_gather_sparce_vec, success() )
{
    int p, s, N, i, j, k, nz;
    float x;
    float * dense, *dense_all;
    float * sparse = NULL;
    int * sparse_ivec = NULL;
    bsp_begin(bsp_nprocs());

    p = bsp_nprocs();
    s = bsp_pid();
    N = 100*p;
    nz = 10*p;

    
    dense = calloc( N/p, sizeof(dense[0]) );
    dense_all = calloc( N, sizeof(dense_all[0])); 

    EXPECT( dense );
    EXPECT( dense_all );
    bsp_push_reg( dense, N/p*sizeof(dense[0]));
    bsp_sync();
    
    for ( i = 0; i < 2*s; ++i)
        rand();

    for ( i = 0; i < nz/p; ++i ) {
        j = (int) (rand() / (1.0 + RAND_MAX) * N/p);
        while (dense[j] != 0.0) j = (j+1)%(N/p);

        x = (float) (rand() / (1.0 + RAND_MAX) );

        EXPECT_OP( "%d", j, >=, 0);
        EXPECT_OP( "%d", j, <, N/p );
        dense[j] = x;

        for ( k = 0; k < 2*p; ++k)
           rand();
    }
    
    for ( i = 0 ; i < p; ++i) 
        bsp_get( i, dense, 0, dense_all + i*N/p, sizeof(float)*N/p);

    nz = all_gather_sparse_vec( dense, N/p, &sparse, &sparse_ivec );

    for ( i = 0; i < nz; ++i) {
        EXPECT_EQ("%g", sparse[i], dense_all[ sparse_ivec[i] ] );
    }

    bsp_pop_reg( dense );

    bsp_end();
}

