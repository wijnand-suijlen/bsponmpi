#ifdef _MSC_VER
#define _CRT_SECURE_NO_WARNINGS
#endif
#include <bsp.h>
#include <bsc.h>
#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <signal.h>
#include <limits.h>
#include "util.h"

typedef int (*conf_interval_method_t)(
        size_t * rng, void * samples, int sample_bytes, 
        int (*get_statistic)(const void * samples, int offset, int n, 
                double * result),
        const int * category_index, int n_categories,
        double conf_level,
        int n_avgs, 
        double * T_avg, double * T_min, double * T_max );

double measure_hrel_memcpy( MPI_Comm comm, char * sendbuf, char * recvbuf,
                            int h, int size, int repeat ) 
{
    double t0, t1;
    int i, j;
    MPI_Barrier( comm );
    t0 = bsp_time();
    for ( j = 0; j < repeat; ++j ) {
        for ( i = 0; i < h; ++i ) {
            memcpy( recvbuf + i*size, sendbuf + i*size, size );
        }
    }
    t1 = bsp_time() - t0;
    MPI_Allreduce( &t1, &t0, 1, MPI_DOUBLE, MPI_MAX, comm );
    return t0 / repeat;
}


double measure_hrel_rma( MPI_Comm comm, MPI_Win win, 
        char * sendbuf, int h, int size, int repeat ) 
{
    double t0, t1;
    int i, j, pid, nprocs;
    MPI_Comm_size( comm, &nprocs);
    MPI_Comm_rank( comm, &pid );
    MPI_Barrier( comm );
    t0 = bsp_time();
    for ( j = 0; j < repeat; ++j ) {
        MPI_Win_fence( 0, win );
        for ( i = 0; i < h; ++i ) {
            MPI_Put( sendbuf + i*size, size, MPI_BYTE,
                    (i+pid) % nprocs, i*size, size, MPI_BYTE,
                    win );
        }
        MPI_Win_fence( 0, win );
    }
    t1 = bsp_time() - t0;
    MPI_Allreduce( &t1, &t0, 1, MPI_DOUBLE, MPI_MAX, comm );
    return t0 / repeat;
}

double measure_hrel_msg( MPI_Comm comm, 
        char * sendbuf, char * recvbuf, MPI_Request * reqs, 
        int h, int size, int repeat ) 
{
    double t0, t1;
    int i, j, pid, nprocs, tag = 0;
    MPI_Comm_size( comm, &nprocs );
    MPI_Comm_rank( comm, &pid );

    MPI_Barrier( comm );
    t0 = bsp_time();
    for ( j = 0; j < repeat; ++j ) {

        for ( i = 0; i < h; ++i ) {
            MPI_Irecv( recvbuf + i*size, size, MPI_BYTE,
                    (pid + nprocs - (i%nprocs)) % nprocs, 
                    tag, comm, &reqs[i] );
        }
        MPI_Barrier( comm );

        for ( i = 0; i < h; ++i ) {
            MPI_Irsend( sendbuf + i*size, size, MPI_BYTE,
                    (i+pid) % nprocs, tag, comm, &reqs[h+i] );
        }
        MPI_Waitall( 2*h, reqs, MPI_STATUSES_IGNORE );
        MPI_Barrier( comm );
    }

    t1 = bsp_time() - t0;
    MPI_Allreduce( &t1, &t0, 1, MPI_DOUBLE, MPI_MAX, comm );
    return t0 / repeat;
}


typedef struct lincost_sample {
    double time;   /* time per msg */
    int    h;      /* # msgs */
    int    msg_size; /* size per msg */
    int    comm;   /* topology */
    int    method; /* 0 == rma, 1 == msg */

    /* info on the order of samples taken */
    int    serial; /* serial number of sample */
    double timestamp; /* moment of sample */
} lincost_sample_t;

void permute( size_t * rng, void * data, 
        int nmemb, int size )
{
    void * tmp = malloc( size );
    char * xs = data;
    while (nmemb > 1 ){
        int i = (int) (nmemb * rand_next(rng));
        if (i != nmemb - 1 ) {
            memcpy( tmp, xs + i * size, size);
            memcpy( xs + i * size, xs + (nmemb-1)*size, size );
            memcpy( xs + (nmemb-1)*size, tmp, size );
        }
        nmemb -= 1;
    }
    free( tmp );
}

int cmp_lincost_sample( const void * a, const void * b )
{
    const lincost_sample_t * x = a, * y = b;

    int i = x->method - y->method;
    int j = x->msg_size - y->msg_size;

    return i?i:j;
}

int cmp_double( const void * a, const void * b )
{
    const double *  x = a, *y = b;
    return *x < *y ? -1 : *x > *y ? 1 : 0;
}

/* Compute alpha and beta s.t. 
 * \sum_{i=0}^n ( y[i] - alpha + beta x[i] )^2 is minimal 
 */
void linear_least_squares( const double *x, const double * y, int n,
        double * alpha, double * beta )
{
    double avg_x = 0.0, avg_y = 0.0;
    double cov_xy = 0.0, var_x = 0.0;
    int i;

    for ( i = 0; i < n; ++i) {
        avg_x += x[i];
        avg_y += y[i];
    }
    avg_x /= n;
    avg_y /= n;

    for ( i = 0; i < n; ++i)  {
        var_x += (x[i] - avg_x)*(x[i] - avg_x);
        cov_xy += (x[i] - avg_x)*(y[i] - avg_y);
    }

    *beta = cov_xy / var_x;
    *alpha = avg_y - (*beta) * avg_x;
}

/* measure_hrel_* measure the time it requires to complete
   a barrier and n message of size m: 
     T(n, m) = gamma + alpha n + beta n m
   Because we are not interested in the time of a barrier
   we use least squares to compute: alpha + beta m
 */
int lincost_samples_get_hrel_time( const void * s,
        int offset, int n, double * result )
{
    int i;
    double *xs = NULL, *ys = NULL;
    double ignore;
    const lincost_sample_t * samples = s;
    int error = 1;

    if (n < 2 ) {
        fprintf(stderr, "ERROR: Insufficient number of samples\n");
        goto exit;
    }

    xs = calloc( n, sizeof(xs[0]));
    ys = calloc( n, sizeof(ys[0]));

    if (!xs || !ys) {
        fprintf(stderr, "ERROR: Insufficient memory\n");
        goto exit;
    }

    for ( i = 0; i < n; ++i ) {
        xs[i] = (double) samples[offset+i].h;
        ys[i] = samples[offset+i].time;
    }

    linear_least_squares( xs, ys, n, &ignore, result );
    error = 0;

exit:
    free(ys);
    free(xs);
    return error;
}

int subsample( size_t * rng, void * samples, int sample_bytes, 
        int (*get_statistic)(const void * samples, int offset, int n, 
                double * result),
        const int * category_index, int n_categories,
        double conf_level, int n_avgs, 
        double * T_avg, double * T_min, double * T_max )
{
    int i, k;
    int error = 1;
    double * avgs = calloc( n_avgs, sizeof(double) );
    if (!avgs) { 
        fprintf(stderr, "ERROR: Insufficient memory for subsampling\n");
        goto exit;
    }

    for ( k = 0; k < n_categories; ++k ) {
        int n = category_index[k+1] - category_index[k]; 
        int l = 0;
        double total_sum = 0.0;
        /* now do to the subsampling: take random samples from the 
         * current samples withour replacement. 
         */
        permute( rng, 
                (char *) samples + sample_bytes * category_index[k], 
                 n, sample_bytes ); 
        for ( i = 0; i < n_avgs; i++) {
            int m = (n + n_avgs - i - 1)/n_avgs;

            if ( (*get_statistic)( samples, category_index[k] + l,
                                   m, &avgs[i] )) {
                goto exit;
            }

            l += m;
            total_sum += avgs[i];
        }
        T_avg[k] = total_sum / n_avgs; /* The average time for the h-rel */

        /* generate CFD of averages */
        qsort( avgs, n_avgs, sizeof(double), cmp_double );

        /* compute the confidence interval */
        T_min[k] = avgs[ (int) (n_avgs * ( 1.0 - conf_level) * 0.5) ];
        T_max[k] = avgs[ (int) (n_avgs * ( 1.0 - (1.0 - conf_level) * 0.5)) ];
    }
    error = 0;

exit:
    free(avgs);
    return error;
}

int bootstrap( size_t * rng, void * samples, int sample_bytes, 
        int (*get_statistic)(const void * samples, int offset, int n, 
                double * result),
        const int * category_index, int n_categories,
        double conf_level,
        int n_avgs, 
        double * T_avg, double * T_min, double * T_max )
{
    int i, j, k;
    int error = 1;
    double * avgs = calloc( n_avgs, sizeof(double) );
    char * bootstrap_sample = NULL;
    int max_n = 0;
    for ( k = 0; k < n_categories; ++k ) { 
        int n = category_index[k+1] - category_index[k]; 
        if (n > max_n) max_n = n;
    }
     
    bootstrap_sample = calloc( max_n, sample_bytes );

    if (!avgs || !bootstrap_sample) {
        fprintf(stderr, "ERROR: Insufficient memory for computing bootstrap\n");
        goto exit;
    }

    for ( k = 0; k < n_categories; ++k ) {
        int n = category_index[k+1] - category_index[k]; 
        double total_sum = 0.0;
        const char * offset_samples 
            = (char *) samples + sample_bytes * category_index[k];

        /* now do to the bootstrapping: take random samples from the 
         * current samples with replacement. 
         */

        for ( i = 0; i < n_avgs; i++) {
            for ( j = 0; j < n; ++j ) {
                int r = (int) (rand_next(rng) * n);
                memcpy( bootstrap_sample + sample_bytes * j,
                        offset_samples + sample_bytes * r,
                        sample_bytes );
            }

            if ( (*get_statistic)( bootstrap_sample, 0, n, &avgs[i] )) {
                goto exit;
            }

            total_sum += avgs[i];
        }
        T_avg[k] = total_sum / n_avgs; /* The average time for the h-rel */

        /* generate CFD of averages */
        qsort( avgs, n_avgs, sizeof(double), cmp_double );

        /* compute the confidence interval */
        T_min[k] = avgs[ (int) (n_avgs * ( 1.0 - conf_level) * 0.5) ];
        T_max[k] = avgs[ (int) (n_avgs * ( 1.0 - (1.0 - conf_level) * 0.5)) ];
    }
    error = 0;
exit:
    free( bootstrap_sample );
    free( avgs );
    return error;
}



static int s_interrupted = 0;

int measure_lincost( int pid, int nprocs, int repeat,
        int msg_size, int total_size, int min_niters, int niters, int ncomms,
        conf_interval_method_t conf_interval_method,
        int * n_samples, 
        const char * sample_file_name,
        double * alpha_msg, double * alpha_msg_min, double * alpha_msg_max, 
        double * alpha_rma, double * alpha_rma_min, double * alpha_rma_max, 
        double * beta_msg, double * beta_msg_min, double * beta_msg_max,
        double * beta_rma, double * beta_rma_min, double * beta_rma_max,
        double * beta_memcpy, double *beta_memcpy_min, double *beta_memcpy_max,
        double conf_level )
{
    char * sendbuf, * recvbuf;
    lincost_sample_t * samples;
    size_t rng = 0;
    int i, k, my_error=0, glob_error=0;
    MPI_Request * reqs = 0; 
    int * pid_perms;
    MPI_Win * wins;  MPI_Comm * comms;
    const int n_avgs = (int) (10.0 / (1.0 - conf_level));
    const int max_h = total_size / msg_size;
    int o[7];  /* two methods x two msg sizes + 1 */
    double T[6], T_min[6], T_max[6]; /* two methods x two msg sizes */
    FILE * sample_file;
    double t0, t1;

    * n_samples = 0;

    if ( pid == 0 && sample_file_name ) {
        sample_file = fopen( sample_file_name, "w" );
        if (sample_file) {
            fprintf( sample_file, "#ID\tSTAMP\tMETHOD\tTOPOLOGY\t"
                                  "MSG_SIZE\tH_REL\tTIME\n" );
            my_error = 0;
        }
        else {
            my_error = 1;
        }
    }
    else {
        sample_file = NULL;
        my_error = 0;
    }

    MPI_Allreduce( &my_error, &glob_error, 1, MPI_INT, 
            MPI_MAX, MPI_COMM_WORLD );
    if ( glob_error ) {
        if (!pid) fprintf( stderr,
                "Could not open file '%s' for writing sample data\n",
                sample_file_name );
        return 1;
    }
 

    /* allocate memory */
    sendbuf = calloc( 2* total_size, 1 );
    recvbuf = calloc( 2* total_size, 1 );
    samples =  calloc( niters, sizeof(samples[0]) );
    reqs = calloc( 2*max_h, sizeof(MPI_Request) );
    pid_perms = calloc( nprocs, sizeof(pid_perms[0]) );
    wins = calloc( ncomms, sizeof(MPI_Win) );
    comms = calloc( ncomms, sizeof(MPI_Comm) ); 

    if (comms) {
        for ( i = 0; i < ncomms; ++i )
            comms[i] = MPI_COMM_NULL;
    }
    if (wins) { 
        for ( i = 0; i < ncomms; ++i) 
            wins[i] = MPI_WIN_NULL;
    }

    my_error = !sendbuf || !recvbuf || !samples || !reqs
            || !pid_perms || !wins || !comms;
    MPI_Allreduce( &my_error, &glob_error, 1, MPI_INT, 
            MPI_MAX, MPI_COMM_WORLD );
    if ( glob_error ) {
        if (!pid) fprintf( stderr,
                "Insufficient memory. For a benchmark on "
                "%d processes I require 2x %d bytes of memory "
                "for send and receive buffers, %ld bytes for "
                "storing the measurements, %ld bytes for "
                "the request queue, %ld bytes for "
                "communication infrastructure\n",
                nprocs,  2*total_size,
                (long) niters * (long) sizeof(samples[0]),
                2*max_h * (long) sizeof(MPI_Request),
               ncomms * (long) sizeof(MPI_Win) +
               ncomms * (long) sizeof(MPI_Comm) +
               nprocs * (long) sizeof(pid_perms[0]) );
        goto exit;
    }

    memset( sendbuf, 1, 2*total_size );
    memset( recvbuf, 2, 2*total_size );

    /* Create random topologies */
    for ( i = 0 ; i < nprocs; ++i )
        pid_perms[i] = i;

    MPI_Comm_dup( MPI_COMM_WORLD, &comms[0] );
    MPI_Win_create( recvbuf, 2*total_size, 1,
            MPI_INFO_NULL, comms[0], & wins[0] );

    for ( i = 1; i < ncomms; ++i ) {
        permute( &rng, pid_perms, nprocs, sizeof(pid_perms[0]) );

        MPI_Comm_split( comms[0], 0, pid_perms[pid], &comms[i] );
        MPI_Win_create( recvbuf, 2*total_size, 1,
            MPI_INFO_NULL, comms[i], & wins[i] );
    }

    /* Create measurement points */
    for ( i = 0 ; i < niters; ++i )
        samples[i].msg_size = (i%2 + 1) * msg_size;

    permute( &rng, samples, niters, sizeof(samples[0]) );

    for ( i = 0; i < niters; ++i ) 
        samples[i].h = i % ( max_h + 1);

    permute( &rng, samples, niters, sizeof(samples[0]) );

    for ( i = 0; i < niters; ++i )
       samples[i].comm = i % ncomms; 

    permute( &rng, samples, niters, sizeof(samples[0]) );

    for ( i = 0 ; i < niters; ++i )
        samples[i].method = i % 3;

    permute( &rng, samples, niters, sizeof(samples[0]) );


    /* Warm up */
    for ( i = 0; i < ncomms; ++i ) {
        measure_hrel_msg( comms[i], sendbuf, recvbuf, reqs, 
                max_h, 2*msg_size, repeat );
        measure_hrel_rma( comms[i], wins[i], sendbuf, 
                max_h, 2*msg_size, repeat );
        measure_hrel_memcpy( comms[i], sendbuf, recvbuf, 
                max_h, 2*msg_size, repeat );
    }

    /* And now the measurements */
    t0 = t1 = bsp_time();
    for ( i = 0; i < niters; ++i ) {
        double t = HUGE_VAL;
        int h    = samples[i].h;
        int size = samples[i].msg_size;
        MPI_Comm comm = comms[ samples[i].comm ];
        MPI_Win win = wins[ samples[i].comm ];

        switch( samples[i].method ) {
            case 0: t = measure_hrel_rma( comm, win, sendbuf, 
                                          h, size, repeat );
                    break;
            case 1: t = measure_hrel_msg( comm, sendbuf, recvbuf, reqs, 
                                          h, size, repeat );
                    break;
            case 2: t = measure_hrel_memcpy( comm, sendbuf, recvbuf,
                                          h, size, repeat );
                    break;
        }

        samples[i].time = t;
        samples[i].timestamp = bsp_time() - t0;
        samples[i].serial = i;

        if ( !pid && bsp_time() - t1 > 10.0 ) {
            printf("# Sample time remaining: %.0f seconds\n"
                   "#    (Press CTRL+C to interrupt)\n",
                    (bsp_time() - t0) / i *(niters - i - 1) );
            t1 = bsp_time();
        }

        MPI_Allreduce( &s_interrupted, &glob_error, 1, MPI_INT,
                MPI_MAX, MPI_COMM_WORLD );
        if ( glob_error ) { 
            if (i > min_niters) {
                niters = i;
                break;
            }
            else {
                goto exit;
            }
        }
    } 

    /* Saving data */
    if (sample_file) {
        for ( i = 0; i < niters; ++i ) {
            fprintf( sample_file, 
                      "%d\t%.15e\t%d\t%d\t%d\t%d\t%.15e\n",
                      samples[i].serial, samples[i].timestamp,
                      samples[i].method, samples[i].comm,
                      samples[i].msg_size, samples[i].h, samples[i].time );
        }
        fclose( sample_file );
        sample_file = NULL;
    }

    /* Analysis */

    /* sort on (msg_size, method) */
    qsort( samples, niters, sizeof(samples[0]), cmp_lincost_sample );
    
    /* create an index to each (msg_size, method) in array 'o' */
    o[0] = 0;
    i = 0;
    for ( k = 0; k < 6; ++k ) {
        for ( ; i < niters; ++i )
            if (samples[i].msg_size != (k%2+1)*msg_size
                || samples[i].method != (k/2) )
                break;

        o[k+1] = i;
    }
 
    my_error = (*conf_interval_method)( &rng, samples, sizeof(samples[0]), 
           lincost_samples_get_hrel_time, 
           o, 6, conf_level, n_avgs,
           T, T_min, T_max );
    MPI_Allreduce( &my_error, &glob_error, 1, MPI_INT, 
            MPI_MAX, MPI_COMM_WORLD );
    if (glob_error) {
        goto exit;
    }
   
    /* saving results */
    *alpha_rma     = 2*T[0] - T[1];
    *alpha_rma_min = 2*T_min[0] - T_max[1];
    *alpha_rma_max = 2*T_max[0] - T_min[1];

    *alpha_msg     = 2*T[2] - T[3];
    *alpha_msg_min = 2*T_min[2] - T_max[3];
    *alpha_msg_max = 2*T_max[2] - T_min[3];

    *beta_rma = (T[1] - T[0])/msg_size ;
    *beta_rma_min = (T_min[1] - T_max[0])/msg_size;
    *beta_rma_max = (T_max[1] - T_min[0])/msg_size;

    *beta_msg = (T[3] - T[2])/msg_size ;
    *beta_msg_min = (T_min[3] - T_max[2])/msg_size;
    *beta_msg_max = (T_max[3] - T_min[2])/msg_size;

    *beta_memcpy = (T[5] - T[4])/msg_size ;
    *beta_memcpy_min = (T_min[5] - T_max[4])/msg_size;
    *beta_memcpy_max = (T_max[5] - T_min[4])/msg_size;

    *n_samples = niters;

exit:
    /* free resources */
    for ( i = 0; i < ncomms; ++i ) {
        if (wins[i] != MPI_WIN_NULL)
            MPI_Win_free( &wins[i] );

        if (comms[i] != MPI_COMM_NULL )
            MPI_Comm_free( &comms[i] );
    }

    free( sendbuf );
    free( recvbuf );
    free( samples );
    free( reqs );
    free( comms );
    free( pid_perms );
    free( wins );
    
    return glob_error;
}


typedef enum bsp_method { PUT, HPPUT, GET, HPGET, SEND, 
    N_BSP_METHODS } bsp_method_t;
static const char * bsp_method_labels[] = { "PUT", "HPPUT",
    "GET", "HPGET", "SEND", "UNKNOWN" };

typedef struct bsp_sample {
    double time;   /* time per h-relation */
    int    h;  /* size of h-relation */
    int    word_size;  /* word size  */
    int    comm;   /* toplogy number */
    bsp_method_t method; /* method */

    /* info on the order of samples taken */
    int    serial; /* serial number of sample */
    double timestamp; /* moment of sample */
} bsp_sample_t;

int cmp_bsp_sample( const void * a, const void * b )
{
    const bsp_sample_t * x = a, * y = b;

    int i = x->method - y->method;
    int j = x->word_size - y->word_size;
    int k = x->h - y->h;
    return i?i:j?j:k;
}

int cmp_int( const void * a, const void * b )
{
    const int * x = a, * y = b;
    return *x - *y;
}

static void bsp_bcast( bsc_pid_t root, void * data, bsc_size_t size ) 
{
    bsp_push_reg( data, size );
    bsp_sync();
    bsc_bcast( bsc_start, root, bsc_all, data, data, size);
    bsc_sync( bsc_flush );
    bsp_pop_reg( data );
}

static void int_or( void * a, const void *a0, const void* xs, int size)
{
    int * result = a;
    const int * zero = a0;
    const int * item = xs;
    const int n = size/sizeof(int);
    int i;

    *result = *zero;
    for ( i = 0 ; i < n; ++i ) {
        *result = *result || item[i];
    }
}


static int bsp_or( int x )
{
    int result = 0;
    const int zero = 0;
    void * tmp = calloc( bsp_nprocs() + 1, sizeof(int) );
    bsp_push_reg(tmp, (bsp_nprocs() + 1)*sizeof(int));
    bsp_sync();
    bsc_allreduce( bsc_start, bsc_all, &x, &result, tmp,
            int_or, &zero, 1, sizeof(int) );
    bsc_sync( bsc_flush );

    bsp_pop_reg(tmp);
    free(tmp);
    return result;
}

static void dbl_max( void * a, const void *a0, const void * xs, int size)
{
    double * result = a;
    const double * neutral = a0;
    const double * item = xs;
    const int n = size/sizeof(double);
    int i;

    *result = *neutral;
    for ( i = 0 ; i < n; ++i ) {
        *result = MAX( *result, item[i] );
    }
}


static double bsp_max( double x )
{
    const double neutral = -HUGE_VAL;
    double result = neutral;
    void * tmp = calloc( bsp_nprocs() + 1, sizeof(double) );
    bsp_push_reg(tmp, (bsp_nprocs() + 1)*sizeof(double));
    bsp_sync();
    bsc_allreduce( bsc_start, bsc_all, &x, &result, tmp,
            dbl_max, &neutral, 1, sizeof(double) );
    bsc_sync( bsc_flush );

    bsp_pop_reg(tmp);
    free(tmp);
    return result;
}


double measure_bsp_hrel( const int * pid_perm, const int * pid_perm_inv,
        char * sendbuf,  char * recvbuf, bsp_method_t method, int word_size, 
        int h, int repeat ) 
{
    double t0, t1;
    int i, k;
    const int pid = pid_perm[bsp_pid()];
    const int nprocs = bsp_nprocs();

    bsp_sync();
    t0 = bsp_time();
    for ( k = 0; k < repeat; ++k ) {
        for ( i = 0; i < h; ++i ) {
            int dst_pid = pid_perm_inv[ (i+pid) % nprocs ];
            int offset = i * word_size;
            switch( method ) {
                case PUT: bsp_put( dst_pid, sendbuf + i * word_size, 
                                   recvbuf, offset, word_size );
                          break;
                case HPPUT: bsp_hpput( dst_pid, sendbuf + i * word_size, 
                                   recvbuf, offset, word_size );
                          break;
                case GET: bsp_get( dst_pid, sendbuf, offset,
                                   recvbuf+ i * word_size, 
                                   word_size );
                          break;
                case HPGET: bsp_hpget( dst_pid, sendbuf, offset,
                                   recvbuf + i * word_size, 
                                   word_size );
                          break;

                case SEND: bsp_send( dst_pid, NULL, 
                                   sendbuf + i * word_size, word_size );
                          break;
                case N_BSP_METHODS:
                          abort();
                          break;
            }
        } /* for i = 0 to h */
        bsp_sync();
    } 
    t1 = bsp_time() ;
    
    return bsp_max(t1-t0) / repeat;
}

int bsp_sample_get_avg_time( const void * s,
        int offset, int n, double * result )
{
    int i;
    double sum = 0.0;
    const bsp_sample_t * samples = s;

    if (n < 1) {
        fprintf(stderr, "ERROR: Insufficient number of samples\n" );
        return EXIT_FAILURE;
    }

    for ( i = 0; i < n ; ++i ) {
        sum += samples[ offset + i ].time;
    }

    *result = sum / n;
    return EXIT_SUCCESS;
}


int measure_bsp_params( int pid, int nprocs, int repeat,
        int word_size, int total_size, int min_niters, int niters, 
        int ncomms, conf_interval_method_t conf_interval_method,
        int * n_samples,
        const char * sample_file_name, 
        double * L, double * L_min, double * L_max,
        double * o, double * o_min, double * o_max,
        double * g, double * g_min, double * g_max, 
        double conf_level )
{
    char * sendbuf, * recvbuf;
    bsp_sample_t * samples;
    int * pid_perm, * pid_perm_inv;
    size_t rng = 0;
    FILE * sample_file = NULL;
    int i, j, k, my_error=0, glob_error=0;
    const int n_avgs = (int) (10.0 / (1.0 - conf_level));
    int *offsets;  /* index into samples */
    double *T, *T_min, *T_max; /* averages + confidence interval */
    double t0, t1 ;

    *n_samples = 0;

    if ( pid == 0 && sample_file_name ) {
        sample_file = fopen( sample_file_name, "w" );
        if (sample_file) {
            fprintf( sample_file, "#ID\tSTAMP\tMETHOD\tTOPOLOGY\t"
                                  "WORD_SIZE\tH_REL\tTIME\n" );
            my_error = 0;
        }
        else {
            my_error = 1;
        }
    }
    else {
        sample_file = NULL;
        my_error = 0;
    }

    glob_error = bsp_or( my_error );
    if ( glob_error ) {
        if (!pid) fprintf( stderr,
                "Could not open file '%s' for writing sample data\n",
                sample_file_name );
        return 1;
    }

    /* allocate memory */
    sendbuf = calloc( 4*total_size, 1 );
    recvbuf = calloc( 4*total_size, 1 );
    samples =  calloc( niters, sizeof(samples[0]) );
    pid_perm = calloc( ncomms * nprocs, sizeof(pid_perm[0]) );
    pid_perm_inv = calloc( ncomms * nprocs, sizeof(pid_perm_inv[0]) );
    offsets = calloc( 1 + 6*N_BSP_METHODS, sizeof(offsets[0]) );
    T = calloc( 6*N_BSP_METHODS, sizeof(T[0]));
    T_min = calloc( 6*N_BSP_METHODS, sizeof(T_min[0]) );
    T_max = calloc( 6*N_BSP_METHODS, sizeof(T_max[0]) );

    my_error = !sendbuf || !recvbuf || !samples 
                || !pid_perm || !pid_perm_inv 
                || !offsets || !T || !T_min || !T_max;

    glob_error = bsp_or( my_error );
    if ( glob_error ) {
        if (!pid) fprintf( stderr,
                "Insufficient memory. For a benchmark on "
                "%d processes I require 2x %ld bytes of memory "
                "for send and receive buffers, %ld bytes for "
                "storing the measurements, %ld bytes for "
                "preparing sample points, and %ld bytes for "
                "analysis.\n",
                nprocs, (long) 4*total_size,
                (long) niters * (long) sizeof(samples[0]),
                (long) 2*ncomms * nprocs * (long) sizeof(pid_perm[0]),
               (long) (1+ 6 * N_BSP_METHODS) * (long) sizeof(offsets[0]) +
               (long) 18 * N_BSP_METHODS * (long) sizeof(T[0])  );
        goto exit;
    }

    memset( sendbuf, 1, 4*total_size );
    memset( recvbuf, 2, 4*total_size );

    bsp_push_reg( sendbuf, 4*total_size );
    bsp_push_reg( recvbuf, 4*total_size );
   
    /* Create random topologies */
    for ( j = 0; j < ncomms; ++j ) {
        for ( i = 0; i < nprocs; ++i ) {
            pid_perm[j*nprocs + i ] = i;
        }
        permute( &rng, pid_perm + j*nprocs, nprocs, sizeof(pid_perm[0]) );

        for ( i = 0; i < nprocs; ++i ) {
            pid_perm_inv[j*nprocs + pid_perm[j*nprocs+ i] ] = i;
        }
    }

    /* Assign word size randomly to each sample */
    for ( i = 0; i < niters; ++i ) 
        samples[i].word_size = (i%2 + 1) * word_size;

    permute( &rng, samples, niters, sizeof(samples[0]) );

    /* Assign h randomly to each sample */
    for ( i = 0; i < niters; ++i ) 
        samples[i].h = total_size / word_size * (i%3);

    permute( &rng, samples, niters, sizeof(samples[0]) );

    /* Assign a topology randomly to each sample */
    for ( i = 0; i < niters; ++i ) 
        samples[i].comm = i % ncomms;

    permute( &rng, samples, niters, sizeof(samples[0]) );

    /* Assign a random bsp method (PUT, GET, etc...) */
    for ( i = 0; i < niters ; ++i ) 
        samples[i].method = i % N_BSP_METHODS;

    permute( &rng, samples, niters, sizeof(samples[0]) );

    /* Warm up */
    if (!pid) {printf("# Warm up phase: "); fflush(stdout); }
    for ( i = 0; i < ncomms; ++i ) {
      if (!pid) printf(" %d/%d; ", i+1, ncomms);
      for ( j = 0; j < N_BSP_METHODS; ++j ) {
        measure_bsp_hrel( pid_perm + i*nprocs, pid_perm_inv + i*nprocs,
                sendbuf, recvbuf, (bsp_method_t) j, word_size, 
                total_size / word_size * 2, repeat );
        measure_bsp_hrel( pid_perm + i*nprocs, pid_perm_inv + i*nprocs,
                sendbuf, recvbuf, (bsp_method_t) j, 2*word_size, 
                total_size / word_size * 2, repeat );
      }
    }
    if (!pid) printf("\n");


    /* And now the measurements */
    t1 = t0 = bsp_time();
    for ( i = 0; i < niters; ++i ) {
        int comm = samples[i].comm;
        int h    = samples[i].h;
        int ws   = samples[i].word_size;
        bsp_method_t m = samples[i].method;

        samples[i].time =
            measure_bsp_hrel( pid_perm + comm*nprocs, 
                    pid_perm_inv + comm*nprocs,
                    sendbuf, recvbuf, m, ws, h, repeat );

        samples[i].timestamp = bsp_time() - t0;
        samples[i].serial = i;
        
        if ( !pid && bsp_time() - t1 > 10.0 ) {
            printf("# Sample time remaining: %.0f seconds\n"
                   "#    (Press CTRL+C to interrupt)\n",
                    (bsp_time() - t0) / i *(niters - i - 1) );
            t1 = bsp_time();
        }

        glob_error = bsp_or( s_interrupted );
        if ( glob_error ) { 
            if (i > min_niters) {
                niters = i;
                break;
            }
            else {
                goto exit;
            }
        }
    } 

    /* Saving data */
    if (sample_file) {
        for ( i = 0; i < niters; ++i ) {
            fprintf( sample_file, "%d\t%.15e\t%s\t%d\t%d\t%d\t%.15e\n",
                      samples[i].serial, samples[i].timestamp,
                      bsp_method_labels[ samples[i].method ],
                      samples[i].comm, samples[i].word_size, 
                      samples[i].h, samples[i].time );
        }
        fclose( sample_file );
        sample_file = NULL;
    }

    /* Analysis */

    /* sort on (msg_size, method) */
    qsort( samples, niters, sizeof(samples[0]), cmp_bsp_sample );
    
    /* create an index 'offsets' to each subarray with the same .bytes value
     * in samples */
    offsets[0] = 0;
    i = 0;
    for ( k = 0; k < 6*N_BSP_METHODS; ++k ) {
        for ( ; i < niters; ++i ) {
            if ( samples[i].h != total_size / word_size * (k%3)
                || samples[i].word_size != word_size * (1+(k/3)%2)
                || samples[i].method != (bsp_method_t) k / 6
               )  break;
        }

        offsets[k+1] = i;
    }
    
    my_error = (*conf_interval_method)( &rng, samples, sizeof(samples[0]), 
           bsp_sample_get_avg_time, 
           offsets, 6*N_BSP_METHODS, conf_level, n_avgs,
           T, T_min, T_max );
    if ( bsp_or( my_error) ) {
        glob_error = 1;
        goto exit;
    }

    for ( k = 0; k < N_BSP_METHODS; ++k ) {
        /* We assume that sending 'n' words of size 'w' requires time
         *
         * T(n, w) = L + n o + n w g
         *
         * where L is the latency, o the overhead per message, and g the
         * reciprocal throughput.
         */

        /* reminder: We now has values of T(n, w) at the following indices
         *
         * T[6 k + 0] = T(0, w)
         * T[6 k + 1] = T(n, w)
         * T[6 k + 2] = T(2n, w)
         * T[6 k + 3] = T(0, 2 w)
         * T[6 k + 4] = T(n, 2 w)
         * T[6 k + 5] = T(2n, 2 w)
         */

        double L_cand;
        int n = total_size / word_size; 
        int nw = n * word_size ;
        int _0_w = 6*k;
        int _n_w = 6*k+1;
        int _2n_w = 6*k+2;
        int _0_2w = 6*k+3;
        int _n_2w = 6*k+4;
        /*  _2n_2w = 6*k+5; unused */
        
        /* get minimum estimate for L >= max{ T(0, w), T(0, 2w) } */
        if ( T[_0_w] < T[_0_2w]) { 
            L[k] = T[_0_2w];
            L_min[k] = T_min[_0_2w];
            L_max[k] = T_max[_0_2w];
        }
        else {
            L[k] = T[_0_w];
            L_min[k] = T_min[_0_w];
            L_max[k] = T_max[_0_w];
        }

        /* estimate assymptotic g = (T(n, 2 w ) - T(n , w)) / (n w)*/
        g[k] = ( T[_n_2w] - T[_n_w] ) / nw;
        g_min[k] = (T_min[_n_2w] - T_max[_n_w]) / nw;
        g_max[k] = (T_max[_n_2w] - T_min[_n_w]) / nw;

        /* correct L if its previous estimate was too low, because
         * 
         *     L >= 2 T(n, w) - T(2 n , w)
         * 
         **/

        L_cand = 2*T[_n_w] - T[_2n_w];
        if (L_cand > L[k]) {
            L[k] = L_cand;
            L_min[k] = 2*T_min[_n_w] - T_max[_2n_w];
            L_max[k] = 2*T_max[_n_w] - T_min[_2n_w];
        }

        /* Estimate the message overhead
         *
         *    o = ( T( 2n, w ) - T(n, 2w ) ) / n
         *
         */
        o[k] = ( T[_2n_w] - T[_n_2w] ) / n;
        o_min[k] = ( T_min[_2n_w] - T_max[_n_2w]) / n;
        o_max[k] = ( T_max[_2n_w] - T_min[_n_2w]) / n;
    }
    *n_samples = niters;

exit:
    /* free resources */
    free( sendbuf ); 
    free( recvbuf ); 
    free( samples );
    free( pid_perm ); 
    free( pid_perm_inv ); 
    free( offsets ); 
    free( T ); 
    free( T_min );
    free( T_max );

    return glob_error;
}


const char * get_param( const char * arg, const char * find ) {
    while ( *arg && *find ) {
        if (*arg != *find ) 
            return NULL;
        arg++;
        find++;
    }

    return arg;
}

typedef enum ci_method { SUBSAMPLE, BOOTSTRAP } ci_method_t;

int parse_params( int argc, char ** argv,
        int * ncomms, int * msg_size,
        double * conf_level, double * ci_rel, 
        ci_method_t * ci_method ,
        int * min_niters, int * max_niters,
        int * repeat, 
        int * avg_total_bytes, 
        int * mes_lincost, 
        char * sample_file_name, int max_file_name_length )
{
    int i;
    const char * val = NULL;

    for ( i = 1; i < argc; ++i ) {
        if ( get_param( argv[i], "--help") ) {
            fprintf(stderr, "Usage %s [--ntopos=<number>] "
                    "[--granularity=<number>]"
                    "[--msg-size=<size>] [--conf-level=<percentage>] "
                    "[--conf-interval=<percentage>] "
                    "[--conf-interval-method=<bootstrap|sample>] "
                    "[--min-samples=<number>] [--max-samples=<number>] "
                    "[--word-size=<bytes>] "
                    "[--avg-total-bytes=<bytes>] "
                    "[--save-sample-data=<file name>] "
                    "[--bsp|--lincost]\n",
                    argv[0]);
            return 1;
        }
        else if ( (val = get_param( argv[i], "--save-sample-data=") ) ) {
            strncpy( sample_file_name, val, max_file_name_length );
            sample_file_name[ max_file_name_length - 1 ] = '\0';
        }
        else if ( (val = get_param( argv[i], "--ntopos=" ) ) ) {
            *ncomms = atoi( val );
        }
        else if ( (val = get_param(argv[i], "--msg-size=") ) ) {
            *msg_size = atoi( val );
        }
        else if ( (val = get_param(argv[i], "--conf-level=") ) ) {
            *conf_level = atof( val ) * 0.01;
        }
        else if ( (val = get_param(argv[i], "--conf-interval=") ) ) {
            *ci_rel = atof( val ) * 0.01;
        }
        else if ( (val = get_param(argv[i], "--conf-interval-method=") ) ) {
            if ( strcmp(val, "subsample") == 0 )
                *ci_method = SUBSAMPLE;
            else if (strcmp(val, "bootstrap") == 0 )
                *ci_method = BOOTSTRAP;
            else
                fprintf(stderr, "Unrecognized confidence level method estimation "
                       "method '%s'\n", val );
        }
        else if ( (val = get_param(argv[i], "--min-samples=") ) ) {
            *min_niters = atoi( val );
        }
        else if ( (val = get_param(argv[i], "--max-samples=") ) ) {
            *max_niters = atoi( val );
        }
        else if ( (val = get_param(argv[i], "--granularity=") ) ) {
            *repeat = atoi( val );
        }
        else if ( (val = get_param(argv[i], "--word-size=") ) ) {
            *msg_size = atoi( val );
        }
        else if ( (val = get_param(argv[i], "--avg-total-bytes=") ) ) {
            *avg_total_bytes = atoi( val );
        }
        else if ( (val = get_param(argv[i], "--bsp") ) ) {
            *mes_lincost = 0;
        }
        else if ( (val = get_param(argv[i], "--lincost") ) ) {
            *mes_lincost = 1;
        }
        else
            fprintf(stderr, "Ignoring unrecognized parameter '%s'\n",
                    argv[i]);
    }
    return 0;
}

static void interrupt( int signal )
{
    (void) signal;
    s_interrupted = 1;
}


int main( int argc, char ** argv )
{
    int nprocs, pid, i;
    double alpha_rma=1e+4, alpha_msg=1e+4;
    double alpha_rma_min=1e+4, alpha_rma_max=1e+4;
    double alpha_msg_min=1e+4, alpha_msg_max=1e+4;
    double beta_rma=1.0, beta_msg=1.0;
    double beta_rma_min=1.0, beta_rma_max=1.0;
    double beta_msg_min=1.0, beta_msg_max=1.0;
    double beta_memcpy = 0.1, beta_memcpy_min = 0.1, beta_memcpy_max = 0.1;
    double L[N_BSP_METHODS], L_min[N_BSP_METHODS], L_max[N_BSP_METHODS];
    double g[N_BSP_METHODS], g_min[N_BSP_METHODS], g_max[N_BSP_METHODS];
    double o[N_BSP_METHODS], o_min[N_BSP_METHODS], o_max[N_BSP_METHODS];
    double conf_level = 0.95;
    double ci_rel = 0.5;
    ci_method_t ci_method = BOOTSTRAP;
    double grow_factor = 2.0;
    double speed_estimate = HUGE_VAL;
    char sample_file_name[256]="";
    int lincost_mode = 1;
    int min_niters = 0, max_niters = INT_MAX, repeat=10;
    int niters = 1000;
    int ncomms = 10;
    int msg_size = 10000;
    int avg_total_bytes = 2500000;
    int error = 0;
    bsp_begin( bsp_nprocs() );
    signal( SIGINT, interrupt );

    pid = bsp_pid();
    nprocs = bsp_nprocs();

    for ( i = 0; i < N_BSP_METHODS; ++i ) {
        L[i]=bsp_nprocs()*1e+4; L_min[i] = 0.0; L_max[i] = HUGE_VAL;
        g[i]=1.0; g_min[i]=0.0; g_max[i]=HUGE_VAL; 
        o[i]=0.0, o_min[i]=0.0; o_max[i]=HUGE_VAL;
    }

    if (!pid) error = parse_params( argc, argv, 
            &ncomms, &msg_size, &conf_level, &ci_rel, &ci_method,
            &niters, &max_niters, &repeat,
            & avg_total_bytes,
            & lincost_mode,
            & sample_file_name[0], sizeof(sample_file_name)  );

    bsp_bcast( 0, &error, sizeof(error) );

    if (error) { 
        goto exit;
    }

    bsp_bcast( 0, &ncomms, sizeof(ncomms) );
    bsp_bcast( 0, &msg_size, sizeof(msg_size));
    bsp_bcast( 0, &conf_level, sizeof(conf_level));
    bsp_bcast( 0, &ci_rel, sizeof(ci_rel));
    bsp_bcast( 0, &ci_method, sizeof(ci_method));
    bsp_bcast( 0, &niters, sizeof(niters));
    bsp_bcast( 0, &max_niters, sizeof(max_niters));
    bsp_bcast( 0, &repeat, sizeof(repeat));
    bsp_bcast( 0, &avg_total_bytes, sizeof(avg_total_bytes) );
    bsp_bcast( 0, &lincost_mode, sizeof(lincost_mode) );
    /* do not broadcast sample_file_name, because it is not necessary */

    if (!pid) printf(
           "# Number of processes       = %10d (bsprun -n)\n"
           "# Confidence level          = %10.2f%% (--conf-level)\n"
           "# Max confidence interval   = %10.2f%% (--conf-interval)\n"
           "# Confidence interval method= %10s (--conf-interval-method)\n"
           "# Random process layouts    = %10d (--ntopos)\n"
           "# Min number of samples     = %10d (--min-samples)\n"
           "# Max number of samples     = %10d (--max-samples)\n"
           "# Granularity               = %10d (--granularity)\n"
           "# Save samples to           = %s (--save-sample-data)\n"
           "# Average communication vol.= %10d bytes (--avg-total-bytes)\n"
           "# Typical message size      = %10d bytes (--msg-size or --word-size)\n",
           nprocs, 100.0*conf_level, 100.0*ci_rel,
           ci_method == SUBSAMPLE? "subsample" :
           ci_method == BOOTSTRAP ? "bootstrap" : "UNKNOWN",
           ncomms,
           niters, max_niters, repeat, 
           sample_file_name[0]=='\0'?"NA":sample_file_name,
           avg_total_bytes, msg_size );

    if (!pid && lincost_mode )
        printf( "\n"
                "# Measuring parameters of linear cost model\n" );

    if (!pid && !lincost_mode)
        printf("\n"
               "# Measuring BSP machine parameters\n" );

                

    do { 
        double t0, t1, dbl_niters;
        int error;
        int performed_samples = 0;

        if (!pid)
            printf("#\n"
                   "#     Taking %d samples ...\n"
                   "#     Expected sampling time: %g seconds\n" 
                   "#     (Press CTRL-C to interrupt)\n",
                   niters, speed_estimate * niters );

        t0 = bsp_time();
        if (lincost_mode) {
            error = measure_lincost( pid, nprocs, repeat, msg_size,
                avg_total_bytes, min_niters, niters, ncomms,
                ci_method==SUBSAMPLE?subsample:bootstrap,
                & performed_samples,
                sample_file_name[0] == '\0' ? NULL : sample_file_name,
                &alpha_msg, &alpha_msg_min, &alpha_msg_max,
                &alpha_rma, &alpha_rma_min, &alpha_rma_max,
                &beta_msg, & beta_msg_min, &beta_msg_max,
                &beta_rma, & beta_rma_min, &beta_rma_max,
                &beta_memcpy, &beta_memcpy_min, &beta_memcpy_max,
                conf_level );
        }
        else {
            if (!pid) {
                const char * min_n_hp_str = getenv( "BSPONMPI_P2P_N_HP");
                if (min_n_hp_str) {
                    int n = atoi( min_n_hp_str );
                    if (n > msg_size ) {
                        fprintf(stderr, "WARNING: Word size should be taken "
                                "greater than %d\n", n );
                    }
                }
            }

            error = measure_bsp_params( pid, nprocs, repeat, msg_size,
                avg_total_bytes, min_niters, niters, ncomms, 
                ci_method==SUBSAMPLE?subsample:bootstrap,
                & performed_samples,
                sample_file_name[0] == '\0' ? NULL : sample_file_name,
                L, L_min, L_max, o, o_min, o_max, g, g_min, g_max,
                conf_level );
        }

        t1 = bsp_time();

        if (!pid && lincost_mode && performed_samples > 0) {
            printf("#     Sampling time was %g seconds; Number of samples %d\n"
                    "# MSG:  alpha = %12.4g [ %12.4g - %12.4g ] ;"
                    "beta =  %12.4g [ %12.4g - %12.4g ] \n"
                    "# RMA:  alpha = %12.4g [ %12.4g - %12.4g ] ;"
                    "beta =  %12.4g [ %12.4g - %12.4g ] \n"
                    "# MEMCPY: beta= %12.4g [ %12.4g - %12.4g ]\n",
                    t1-t0, performed_samples,
                    alpha_msg, alpha_msg_min, alpha_msg_max,
                    beta_msg, beta_msg_min, beta_msg_max,
                    alpha_rma, alpha_rma_min, alpha_rma_max,
                    beta_rma, beta_rma_min, beta_rma_max,
                    beta_memcpy, beta_memcpy_min, beta_memcpy_max );
        }

        if (!pid && !lincost_mode && performed_samples > 0) {
            printf("#     Sampling time was %g seconds; "
                   "Number of samples %d\n", t1-t0, performed_samples );
            for ( i = 0; i < N_BSP_METHODS; ++i ) {
                printf( "# %5s :: L = %12.4g [ %12.4g - %12.4g ] ;"
                        " o = %12.4g [ %12.4g - %12.4g ] ;"
                        " g = %12.4g [ %12.4g - %12.4g ]\n",
                        bsp_method_labels[ i ],
                        L[i], L_min[i], L_max[i], 
                        o[i], o_min[i], o_max[i],
                        g[i], g_min[i], g_max[i]);
            }
        }
        if (error) {
            if (!pid)
                printf("# Interrupted: target precision was not achieved\n");
             break;
        }
            
        speed_estimate = (t1-t0)/niters;

        if ( lincost_mode ) {
            grow_factor = 
                MAX( MAX( MAX( 
                    fabs((alpha_rma_max - alpha_rma_min) / alpha_rma),
                    fabs((alpha_msg_max - alpha_msg_min) / alpha_msg) ),  
                    fabs((beta_msg_max - beta_msg_min) / beta_msg) ),
                    fabs((beta_rma_max - beta_rma_min) / beta_rma ) ) 
                    / ci_rel;
        } else {
            for ( i = 0; i < N_BSP_METHODS; ++i ) {
                grow_factor = MAX( MAX( 
                        fabs(( L_max[i] - L_min[i]) / L[i]), 
                        fabs(( o_max[i] - o_min[i]) / o[i])),
                        fabs(( g_max[i] - g_min[i]) / g[i]))
                     / ci_rel;
          }
        }

        if ( niters == max_niters ) {
            if (!pid)
                printf("# Warning: target precision was not achieved\n");
            break;
        }

        min_niters = niters;
        dbl_niters = MIN(100.0,MAX(2.0, grow_factor)) * niters ;
        if ( dbl_niters > (double) max_niters )
            niters = max_niters;
        else
            niters = (int) dbl_niters;

    } while ( grow_factor > 1.0 );

    if (!pid && lincost_mode ) {
        /* n_hp is the minimum message size at which bsp_hpget or bsp_hpput
         * is faster than buffered bsp_put and bsp_get 
         *
         * n_hp = (alpha - bsp_msg_overhead) / ( 2*beta_memcpy )
         *
         * in reality bsp_msg_overhead is negligible compared to alpha,
         * so we can assume it to be zero.
         * */

        long n_hp_rma = (long) (alpha_rma / (2 * beta_memcpy ));
        long n_hp_msg = (long) (alpha_msg / (2 * beta_memcpy ));

        if ( alpha_rma < 0 || alpha_msg < 0 ) {
            printf("# WARNING Point-to-point latency was measured to be negative.\n"
                   "#         This may be caused by super-linear cost of sending a\n"
                   "#         message. Reduce communication volume with --avg-total-bytes\n");
        }
        if ( beta_rma < 0 || beta_msg < 0 || beta_memcpy < 0 ) {
            printf("# WARNING Reciprocal throughput was measured to be negative.\n"
                   "#         Increase the message size with --msg-size\n");
        }
        if ( alpha_rma < 0 || alpha_msg < 0 || beta_rma < 0 || beta_msg < 0 ) {
            printf("# WARNING One of the machine cost parameters is negative, which is\n"
                    "         nonsense. Please inspect the sample data manually\n"
                    "         (use --save-sample-data)\n");
        }

        printf(
           "if [ x${BSPONMPI_A2A_METHOD} = xrma ]; then\n"
           "   export BSPONMPI_P2P_LATENCY=\"%.2g\"\n"
           "   export BSPONMPI_P2P_MSGGAP=\"%.2g\"\n"
           "   export BSPONMPI_P2P_N_HP=%ld\n"
           "elif [ x${BSPONMPI_A2A_METHOD} = xmsg ]; then\n"
           "   export BSPONMPI_P2P_LATENCY=\"%.2g\"\n"
           "   export BSPONMPI_P2P_MSGGAP=\"%.2g\"\n"
           "   export BSPONMPI_P2P_N_HP=%ld\n"
           "fi\n", 
           MAX(0.0, alpha_rma),
           MAX(0.0, beta_rma),
           MAX(0l, n_hp_rma),
           MAX(0.0, alpha_msg), 
           MAX(0.0, beta_msg),
           MAX(0l, n_hp_msg) );

        printf("\n# Save this output to a file, and use it with bsprun's\n"
                 "# --use-benchmark-file parameter\n");
    }

    if (!pid && !lincost_mode ) {
        const char * mode = getenv("BSPONMPI_A2A_METHOD");
        if (!mode) mode = "rma";

        for ( i = 0 ; i < N_BSP_METHODS; ++i ) {
            if ( g[i] < 0 ) {
                printf("# WARNING Reciprocal throughput was measured to be negative.\n"
                       "#         Increase the word size with --word-size\n");
            }
            if ( o[i] < 0 ) {
                printf("# WARNING Message overhead was measured to be negative.\n"
                       "#         Decrease the word size with --word-size and/or\n"
                       "#         increase the total volume with --avg-total-bytes.\n");
            }
        }

        printf("if [ x${BSPONMPI_A2A_METHOD} = x%s ]; then\n", mode);
        for ( i = 0; i < N_BSP_METHODS; ++i ) {
            printf(
               "  export BSC_%s_L=\"%.2g\"\n"
               "  export BSC_%s_G=\"%.2g\"\n"
               "  export BSC_%s_O=\"%.2g\"\n",
                bsp_method_labels[i], MAX(0.0, L[i]), 
                bsp_method_labels[i], MAX(0.0, g[i]),
                bsp_method_labels[i], MAX(0.0, o[i]) );
        }
        printf("fi\n");
        printf("\n# Save this output to a file, and use it with bsprun's\n"
                 "# --use-benchmark-file parameter\n");
    }

exit:
    bsp_end();
    return EXIT_SUCCESS;
}
