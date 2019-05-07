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
                    (i+pid) % nprocs, pid*size, size, MPI_BYTE,
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
                    (pid + nprocs - i) % nprocs, tag, comm, &reqs[i] );
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
        int i = nmemb * rand_next(rng);
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
    return x < y ? -1 : x > y ? 1 : 0;
}

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


static int s_interrupted = 0;

int measure_lincost( int pid, int nprocs, int repeat,
        int msg_size, int min_niters, int niters, int ncomms,
        const char * sample_file_name,
        double * alpha_msg, double * alpha_msg_ci, 
        double * alpha_rma, double * alpha_rma_ci, 
        double * beta_msg, double * beta_msg_ci,
        double * beta_rma, double * beta_rma_ci,
        double * beta_memcpy, double * beta_memcpy_ci,
        double conf_level )
{
    char * sendbuf, * recvbuf;
    lincost_sample_t * samples;
    size_t rng = 0;
    int i, j, k, my_error=0, glob_error=0;
    MPI_Request * reqs = 0; 
    int * pid_perms;
    MPI_Win * wins;  MPI_Comm * comms;
    const int n_avgs = 10.0 / (1.0 - conf_level);
    double * avgs, *xs, *ys;
    int o[7];  /* two methods x two msg sizes + 1 */
    double T[6], T_ci[6]; /* two methods x two msg sizes */
    FILE * sample_file;
    double t0, t1;

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
    sendbuf = calloc( (long) 2 * msg_size * nprocs, 1 );
    recvbuf = calloc( (long) 2 * msg_size * nprocs, 1 );
    samples =  calloc( niters, sizeof(samples[0]) );
    reqs = calloc( 2*nprocs, sizeof(MPI_Request) );
    pid_perms = calloc( nprocs, sizeof(pid_perms[0]) );
    wins = calloc( ncomms, sizeof(MPI_Win) );
    comms = calloc( ncomms, sizeof(MPI_Comm) ); 
    avgs = calloc( n_avgs, sizeof(avgs[0]) );
    xs = calloc( niters, sizeof(xs[0]));
    ys = calloc( niters, sizeof(ys[0]));

    if (comms) {
        for ( i = 0; i < ncomms; ++i )
            comms[i] = MPI_COMM_NULL;
    }
    if (wins) { 
        for ( i = 0; i < ncomms; ++i) 
            wins[i] = MPI_WIN_NULL;
    }

    my_error = !sendbuf || !recvbuf || !samples || !reqs
            || !pid_perms || !wins || !comms || !avgs || !xs || !ys ;
    MPI_Allreduce( &my_error, &glob_error, 1, MPI_INT, 
            MPI_MAX, MPI_COMM_WORLD );
    if ( glob_error ) {
        if (!pid) fprintf( stderr,
                "Insufficient memory. For a benchmark on "
                "%d processes I require 2x %ld bytes of memory "
                "for send and receive buffers, %ld bytes for "
                "storing the measurements, %ld bytes for "
                "the request queue, %ld bytes for "
                "communication infrastructure, and %ld bytes for "
                "analysis.\n",
                nprocs, (long) 2*msg_size * nprocs,
                (long) niters * (long) sizeof(samples[0]),
                (long) 2*nprocs*sizeof(MPI_Request),
               ncomms * (long) sizeof(MPI_Win) +
               ncomms * (long) sizeof(MPI_Comm) +
               nprocs * (long) sizeof(pid_perms[0]), 
               (long) n_avgs * (long) sizeof(avgs[0]) + 
               (long) 2*niters*sizeof(double) );
        goto exit;
    }

    memset( sendbuf, 1, (long) 2*msg_size*nprocs );
    memset( recvbuf, 2, (long) 2*msg_size*nprocs );

    /* Create random topologies */
    for ( i = 0 ; i < nprocs; ++i )
        pid_perms[i] = i;

    MPI_Comm_dup( MPI_COMM_WORLD, &comms[0] );
    MPI_Win_create( recvbuf, (long) 2*msg_size * nprocs, 1,
            MPI_INFO_NULL, comms[0], & wins[0] );

    for ( i = 1; i < ncomms; ++i ) {
        permute( &rng, pid_perms, nprocs, sizeof(pid_perms[0]) );

        MPI_Comm_split( comms[0], 0, pid_perms[pid], &comms[i] );
        MPI_Win_create( recvbuf, (long) 2*msg_size * nprocs, 1,
            MPI_INFO_NULL, comms[i], & wins[i] );
    }

    /* Create measurement points */
    for ( i = 0; i < niters; ++i ) 
        samples[i].h = i % (nprocs+1);

    permute( &rng, samples, niters, sizeof(samples[0]) );

    for ( i = 0; i < niters; ++i )
       samples[i].comm = i % ncomms; 

    permute( &rng, samples, niters, sizeof(samples[0]) );

    for ( i = 0 ; i < niters; ++i )
        samples[i].msg_size = (i%2 + 1) * msg_size;

    permute( &rng, samples, niters, sizeof(samples[0]) );

    for ( i = 0 ; i < niters; ++i )
        samples[i].method = i % 3;

    permute( &rng, samples, niters, sizeof(samples[0]) );


    /* Warm up */
    for ( i = 0; i < ncomms; ++i ) {
        measure_hrel_msg( comms[i], 
                sendbuf, recvbuf, reqs, nprocs, 2*msg_size, repeat );
        measure_hrel_rma( comms[i], wins[i], 
                sendbuf, nprocs, 2*msg_size, repeat );
        measure_hrel_memcpy( comms[i], sendbuf, recvbuf, nprocs,
                2*msg_size, repeat );
    }

    /* And now the measurements */
    t0 = t1 = bsp_time();
    for ( i = 0; i < niters; ++i ) {
        double t;
        int h    = samples[i].h;
        int size = samples[i].msg_size;
        MPI_Comm comm = comms[ samples[i].comm ];
        MPI_Win win = wins[ samples[i].comm ];

        switch( samples[i].method ) {
            case 0: t = measure_hrel_rma( comm, win, sendbuf, h, size, 
                                        repeat );
                    break;
            case 1: t = measure_hrel_msg( comm, sendbuf, recvbuf, reqs, h, 
                                        size, repeat );
                    break;
            case 2: t = measure_hrel_memcpy( comm, sendbuf, recvbuf, h, 
                                        size, repeat );
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
    
    for ( k = 0; k < 6; ++k ) {
        int n = o[k+1] - o[k]; 
        int l = 0;
        double total_sum = 0.0;
        double left, right;
        /* now do to the subsampling: take random samples from the 
         * current samples. 
         */
        permute( &rng, samples + o[k], n, sizeof(samples[0])); 
        for ( i = 0; i < n_avgs; i++) {
            double T_barrier, T_hrel;
            int m = (n + n_avgs - i - 1)/n_avgs;

            if ( m < 2 ) {
                glob_error = 1;
                goto exit;
            }

            /* Compute a linear least squares fit.
             * Denote f(h) the recorded time required by a h-relation,
             * then f(h) = T_barrier + h ( alpha + n beta ).
             * We want to know the: alpha + n beta */
            for ( j = 0; j < m; ++j ) {
                xs[j] = (double) samples[o[k]+l].h;
                ys[j] = samples[o[k]+l].time;
                l++;
            }

            linear_least_squares( xs, ys, m, &T_barrier, &T_hrel );

            avgs[i] = T_hrel;
            total_sum += T_hrel; /* = alpha + n beta */
        }

        T[k] = total_sum / n_avgs; /* The average time for alpha + n beta */

        /* generate CFD of averages */
        qsort( avgs, n_avgs, sizeof(double), cmp_double );

        left = avgs[ (int) (n_avgs * conf_level * 0.5) ];
        right = avgs[ (int) (n_avgs * ( 1.0 - conf_level * 0.5)) ];

        /* compute the confidence interval */
        T_ci[k] = MAX( fabs(T[k] - left),  fabs(right - T[k]) );
    }

    /* saving results */
    *alpha_rma_ci = fabs(2*T_ci[0] + T_ci[1] );
    *alpha_rma = MAX( *alpha_rma_ci, 2*T[0] - T[1]);

    *alpha_msg_ci = fabs(2*T_ci[2] + T_ci[3]);
    *alpha_msg = MAX( *alpha_msg_ci, 2*T[2] - T[3]);

    *beta_rma_ci = (T_ci[1] + T_ci[0]) / msg_size;
    *beta_rma = MAX( *beta_rma_ci, (T[1] - T[0])/msg_size );

    *beta_msg_ci = (T_ci[3] + T_ci[2]) / msg_size;
    *beta_msg = MAX( *beta_msg_ci, (T[3] - T[2])/msg_size );

    *beta_memcpy_ci = (T_ci[5] + T_ci[4]) / msg_size;
    *beta_memcpy = MAX( *beta_msg_ci, (T[5] - T[4])/msg_size );

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
    free( avgs );
    free( xs );
    free( ys );

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
        int h, int max_total_bytes, int repeat ) 
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
            int offset = ((i/nprocs)*nprocs + pid) * word_size;
            if ( offset + word_size <= max_total_bytes ) {
                switch( method ) {
                    case PUT: bsp_put( dst_pid, sendbuf + i * word_size, 
                                       recvbuf, offset, word_size );
                              break;
                    case HPPUT: bsp_hpput( dst_pid, sendbuf + i * word_size, 
                                       recvbuf, offset, word_size );
                              break;
                    case GET: bsp_get( dst_pid, recvbuf, offset,
                                       sendbuf + i * word_size, 
                                       word_size );
                              break;
                    case HPGET: bsp_hpget( dst_pid, recvbuf, offset,
                                       sendbuf + i * word_size, 
                                       word_size );
                              break;

                    case SEND: bsp_send( dst_pid, NULL, 
                                       sendbuf + i * word_size, word_size );
                              break;
                    case N_BSP_METHODS:
                              abort();
                              break;
                }
            }
        } /* for i = 0 to h */
        bsp_sync();
    } 
    t1 = bsp_time() ;
    
    return bsp_max(t1-t0) / repeat;
}



int measure_bsp_params( int pid, int nprocs, int repeat,
        int pessimistic_word_size, int optimistic_word_size,
        int max_total_bytes, int nhrels, int min_niters, int niters, 
        int ncomms, const char * sample_file_name,
        double * L, double * L_ci,
        double * g_pes, double * g_pes_ci, 
        double * g_opt, double * g_opt_ci,
        double conf_level )
{
    char * sendbuf, * recvbuf;
    bsp_sample_t * samples;
    int * pid_perm, * pid_perm_inv;
    size_t rng = 0;
    FILE * sample_file = NULL;
    int i, j, k, my_error=0, glob_error=0;
    int * hrels;
    const int n_avgs = 10.0 / (1.0 - conf_level);
    double * avgs, *xs, *ys;
    int *o;  /* index into samples */
    double *T, *T_ci; /* averages + confidence interval */
    double t0, t1 ;

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
    sendbuf = calloc( max_total_bytes, 1 );
    recvbuf = calloc( max_total_bytes, 1 );
    samples =  calloc( niters, sizeof(samples[0]) );
    pid_perm = calloc( ncomms * nprocs, sizeof(pid_perm[0]) );
    pid_perm_inv = calloc( ncomms * nprocs, sizeof(pid_perm_inv[0]) );
    hrels = calloc( nhrels, sizeof(hrels[0]) );
    avgs = calloc( n_avgs, sizeof(avgs[0]) );
    xs = calloc( niters, sizeof(xs[0]));
    ys = calloc( niters, sizeof(ys[0]));
    o = calloc( (2*nhrels+2)*N_BSP_METHODS, sizeof(o[0]) );
    T = calloc( 2*nhrels*N_BSP_METHODS, sizeof(T[0]));
    T_ci = calloc( 2*nhrels*N_BSP_METHODS, sizeof(T_ci[0]) );

    my_error = !sendbuf || !recvbuf || !samples 
                || !pid_perm || !pid_perm_inv || !hrels
                || !avgs || !xs || !ys || !o || !T || !T_ci;

    glob_error = bsp_or( my_error );
    if ( glob_error ) {
        if (!pid) fprintf( stderr,
                "Insufficient memory. For a benchmark on "
                "%d processes I require 2x %ld bytes of memory "
                "for send and receive buffers, %ld bytes for "
                "storing the measurements, %ld bytes for "
                "preparing sample points, and %ld bytes for "
                "analysis.\n",
                nprocs, (long) max_total_bytes,
                (long) niters * (long) sizeof(samples[0]),
                (long) 2*ncomms * nprocs * (long) sizeof(pid_perm[0])+
                (long) nhrels * (long) sizeof(hrels[0]),
               (long) n_avgs * (long) sizeof(avgs[0]) + 
               (long) 2*niters* (long) sizeof(double) +
               (long) (2*nhrels+2) * N_BSP_METHODS * (long) sizeof(o[0]) +
               (long) 4*nhrels * N_BSP_METHODS * (long) sizeof(T[0])  );
        goto exit;
    }

    memset( sendbuf, 1, max_total_bytes );
    memset( recvbuf, 2, max_total_bytes );

    bsp_push_reg( sendbuf, max_total_bytes );
    bsp_push_reg( recvbuf, max_total_bytes );
   
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

    /* Define nhrels measurements at equidistant intervals */
    hrels[0] = 0;
    for ( i = 1; i < nhrels ; ++i ) {
        hrels[i] = (int) ( (double) i / nhrels * max_total_bytes );
    }

    /* Assign word size randomly to each sample */
    for ( i = 0; i < niters; ++i ) 
        samples[i].word_size = 
            (i%2) ? pessimistic_word_size : optimistic_word_size;

    permute( &rng, samples, niters, sizeof(samples[0]) );

    /* Assign h randomly to each sample */
    for ( i = 0; i < niters; ++i ) 
        samples[i].h = hrels[i % nhrels] / samples[i].word_size;

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
                sendbuf, recvbuf, (bsp_method_t) j, pessimistic_word_size, 
                max_total_bytes / pessimistic_word_size, 
                max_total_bytes, repeat );
        measure_bsp_hrel( pid_perm + i*nprocs, pid_perm_inv + i*nprocs,
                sendbuf, recvbuf, (bsp_method_t) j, optimistic_word_size, 
                max_total_bytes / optimistic_word_size,
                max_total_bytes, repeat );
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
                    sendbuf, recvbuf, m, ws, h,
                    max_total_bytes, repeat );

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
    
    /* create an index 'o' to each subarray with the same .bytes value
     * in samples */
    o[0] = 0;
    i = 0;
    for ( k = 0; k < 2*nhrels*N_BSP_METHODS; ++k ) {
        for ( ; i < niters; ++i ) {
            if ( samples[i].h != hrels[k%nhrels] / samples[i].word_size 
                || samples[i].word_size != 
                 ((k/nhrels%2)? optimistic_word_size : pessimistic_word_size)
                || samples[i].method != (bsp_method_t) k / (2 * nhrels)
               )  break;
        }

        o[k+1] = i;
    }
    
    for ( k = 0; k < 2*nhrels*N_BSP_METHODS; ++k ) {
        int n = o[k+1] - o[k]; 
        int l = 0;
        double total_sum = 0.0;
        double left, right;
        /* now do to the subsampling: take random samples from the 
         * current samples. 
         */
        permute( &rng, samples + o[k], n, sizeof(samples[0])); 
        for ( i = 0; i < n_avgs; i++) {
            double sum = 0.0;
            int m = (n + n_avgs - i - 1)/n_avgs;

            if ( m < 1 ) {
                glob_error = 1;
                goto exit;
            }

            for ( j = 0; j < m ; ++ j ) {
                sum += samples[ o[k] + l ].time;
                l++;
            }

            avgs[i] = sum / m;
            total_sum += sum;
        }
        T[k] = total_sum / n; /* The average time for the h-rel */

        /* generate CFD of averages */
        qsort( avgs, n_avgs, sizeof(double), cmp_double );

        left = avgs[ (int) (n_avgs * conf_level * 0.5) ];
        right = avgs[ (int) (n_avgs * ( 1.0 - conf_level * 0.5)) ];

        /* compute the confidence interval */
        T_ci[k] = MAX( fabs(T[k] - left),  fabs(right - T[k]) );
    }

    for ( k = 0; k < N_BSP_METHODS; ++k ) {
        /* now infer L and g results */
        if ( T[2*k*nhrels] < T[(2*k+1)*nhrels] ) { 
            L[k] = T[2*k*nhrels];
            L_ci[k] = T_ci[2*k*nhrels];
        }
        else {
            L[k] = T[(2*k+1)*nhrels];
            L_ci[k] = T_ci[(2*k+1)*nhrels];
        }

        g_opt[k] = 0.0;
        g_pes[k] = 0.0;
        g_opt_ci[k] = HUGE_VAL;
        g_pes_ci[k] = HUGE_VAL;
        for ( i = 1 ; i < nhrels-1; ++i) {
            const int ws_pes = pessimistic_word_size;
            const int ws_opt = optimistic_word_size;
            const int h_pes_a = hrels[i] / ws_pes, h_pes_b = hrels[i+1]/ws_pes;
            const int h_opt_a = hrels[i] / ws_opt, h_opt_b = hrels[i+1]/ws_opt;
            const int bytes_pes_a = ws_pes * h_pes_a;
            const int bytes_pes_b = ws_pes * h_pes_b;
            const int bytes_opt_a = ws_opt * h_opt_a;
            const int bytes_opt_b = ws_opt * h_opt_b;

            double g_opt_cand = (T[(2*k+1)*nhrels+i+1]-T[(2*k+1)*nhrels+i]) 
                                / ( bytes_opt_b - bytes_opt_a);
            double g_pes_cand = (T[2*k*nhrels + i+1]-T[2*k*nhrels+i]) 
                                / ( bytes_pes_b - bytes_pes_a);

            if ( g_opt[k] < g_opt_cand ) {
                double L_low;
                g_opt[k] = g_opt_cand;
                g_opt_ci[k] = (T_ci[(2*k+1)*nhrels+i-1]+T_ci[(2*k+1)*nhrels+i]) 
                                / (bytes_opt_b - bytes_opt_a);

                L_low = T[(2*k+1)*nhrels+i] - T_ci[(2*k+1)*nhrels+i] 
                                - ( g_opt[k] + g_opt_ci[k] ) * bytes_opt_a ;
                if (L_low > L[k]) {
                    L[k] = T[(2*k+1)*nhrels+i] - g_opt[k] * bytes_opt_a;
                    L_ci[k] = L[k] - L_low;
                }
            }

            if ( g_pes[k] < g_pes_cand ) {
                double L_low;
                g_pes[k] = g_pes_cand;
                g_pes_ci[k] = (T_ci[2*k*nhrels+i-1]+T_ci[2*k*nhrels+i]) 
                                / (bytes_pes_b - bytes_pes_a);

                L_low = T[2*k*nhrels+i] - T_ci[2*k*nhrels+i] 
                        - ( g_pes[k] + g_pes_ci[k] ) * bytes_pes_a ;
                if (L_low > L[k]) {
                    L[k] = T[2*k*nhrels+i] - g_pes[k] * bytes_pes_a;
                    L_ci[k] = L[k] - L_low;
                }
            }
        }
    }


exit:
    /* free resources */
    free( sendbuf ); 
    free( recvbuf ); 
    free( samples );
    free( pid_perm ); 
    free( pid_perm_inv ); 
    free( hrels ); 
    free( avgs ); 
    free( xs ); 
    free( ys );
    free( o ); 
    free( T ); 
    free( T_ci );

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

int parse_params( int argc, char ** argv,
        int * ncomms, int * msg_size,
        double * conf_level, int * digits, 
        int * min_niters, int * max_niters,
        int * repeat, int * pessimistic_word_size, 
        int * optimistic_word_size,
        int * max_total_bytes, int * nhrels,
        int * mes_lincost, 
        char * sample_file_name, int max_file_name_length )
{
    int i;
    const char * val = NULL;

    for ( i = 1; i < argc; ++i ) {
        if ( get_param( argv[i], "--help") ) {
            fprintf(stderr, "Usage %s [--ncomms=<number>] "
                    "[--granularity=<number>]"
                    "[--msg-size=<size>] [--conf-level=<percentage>] "
                    "[--digits=<number>] "
                    "[--min-samples=<number>] [--max-samples=<number>] "
                    "[--pessimistic-word-size=<bytes>] [--optimistic-word-size=<bytes>] "
                    "[--max-total-bytes=<bytes>] [--nhrels=<number>] "
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
        else if ( (val = get_param(argv[i], "--digits=") ) ) {
            *digits = atoi( val );
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
        else if ( (val = get_param(argv[i], "--pessimistic-word-size=") ) ) {
            *pessimistic_word_size = atoi( val );
        }
        else if ( (val = get_param(argv[i], "--optimistic-word-size=") ) ) {
            *optimistic_word_size = atoi( val );
        }
        else if ( (val = get_param(argv[i], "--max-total-bytes=") ) ) {
            *max_total_bytes = atoi( val );
        }
        else if ( (val = get_param(argv[i], "--nhrels=") ) ) {
            *nhrels = atoi( val );
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
    double alpha_rma_ci=1e+4, alpha_msg_ci=1e+4;
    double beta_rma=1.0, beta_msg=1.0;
    double beta_rma_ci=1.0, beta_msg_ci=1.0;
    double beta_memcpy = 0.1, beta_memcpy_ci = 0.1;
    double L[N_BSP_METHODS], L_ci[N_BSP_METHODS];
    double g_pes[N_BSP_METHODS], g_pes_ci[N_BSP_METHODS];
    double g_opt[N_BSP_METHODS], g_opt_ci[N_BSP_METHODS];
    double conf_level = 0.95;
    double grow_factor = 2.0;
    double speed_estimate = HUGE_VAL;
    char sample_file_name[256]="";
    int lincost_mode = 1;
    int min_niters = 0, max_niters = INT_MAX, repeat=10;
    int niters = 20000;
    int ncomms = 10;
    int msg_size = 2500000;
    int pessimistic_word_size = 8;
    int optimistic_word_size = 1000000;
    int max_total_bytes = 4000000;
    int nhrels = 3; 
    int digits = 1, pow_digits=10;
    int error = 0;
    bsp_begin( bsp_nprocs() );
    signal( SIGINT, interrupt );

    pid = bsp_pid();
    nprocs = bsp_nprocs();
    max_total_bytes = nprocs * optimistic_word_size;

    for ( i = 0; i < N_BSP_METHODS; ++i ) {
        L[i]=bsp_nprocs()*1e+4;
        L_ci[i] = bsp_nprocs()*1e+4;
        g_pes[i]=1.0; g_pes_ci[i]=1.0; 
        g_opt[i]=1.0, g_opt_ci[i]=1.0;
    }

    if (!pid) error = parse_params( argc, argv, 
            &ncomms, &msg_size, &conf_level, &digits,
            &niters, &max_niters, &repeat,
            & pessimistic_word_size, & optimistic_word_size, 
            & max_total_bytes, &nhrels,
            & lincost_mode,
            & sample_file_name[0], sizeof(sample_file_name)  );

    bsp_bcast( 0, &error, sizeof(error) );

    if (error) { 
        goto exit;
    }

    bsp_bcast( 0, &ncomms, sizeof(ncomms) );
    bsp_bcast( 0, &msg_size, sizeof(msg_size));
    bsp_bcast( 0, &conf_level, sizeof(conf_level));
    bsp_bcast( 0, &digits, sizeof(digits));
    bsp_bcast( 0, &niters, sizeof(niters));
    bsp_bcast( 0, &max_niters, sizeof(max_niters));
    bsp_bcast( 0, &repeat, sizeof(repeat));
    bsp_bcast( 0, &optimistic_word_size, sizeof(optimistic_word_size) );
    bsp_bcast( 0, &pessimistic_word_size, sizeof(pessimistic_word_size) );
    bsp_bcast( 0, &max_total_bytes, sizeof(max_total_bytes) );
    bsp_bcast( 0, &nhrels, sizeof(nhrels) );
    bsp_bcast( 0, &lincost_mode, sizeof(lincost_mode) );
    /* do not broadcast sample_file_name, because it is not necessary */

    if (!pid) printf(
           "# Number of processes    = %d\n"
           "# Confidence level       = %.2f%%\n"
           "# Significant digits     = %d\n"
           "# Random process layouts = %d\n"
           "# Max number of samples  = %d\n"
           "# Granularity            = %d\n"
           "# Save samples to        = %s\n",
           nprocs, 100.0*conf_level, digits, ncomms,
           max_niters, repeat, 
           sample_file_name[0]=='\0'?"NA":sample_file_name );

    if (!pid && lincost_mode )
        printf( "\n"
                "# Measuring parameters of linear cost model\n"
                "# Typical message size   = %d bytes\n",
                msg_size );

    if (!pid && !lincost_mode)
        printf("\n"
               "# Measuring BSP machine parameters\n"
               "# Optimistic word size    = %d bytes\n"
               "# Pessimistic word size   = %d bytes\n"
               "# Size of max h-relation  = %d bytes\n"
               "# Number of h-relations   = %d\n",
               optimistic_word_size, pessimistic_word_size,
               max_total_bytes, nhrels );

                

    do { 
        double t0, t1, dbl_niters;
        int error;

        if (!pid)
            printf("#\n"
                   "#     Taking %d samples ...\n"
                   "#     Expected sampling time: %g seconds\n" 
                   "#     (Press CTRL-C to interrupt)\n",
                   niters, speed_estimate * niters );

        t0 = bsp_time();
        if (lincost_mode) {
            error = measure_lincost( pid, nprocs, repeat,
                msg_size,  min_niters, niters, ncomms,
                sample_file_name[0] == '\0' ? NULL : sample_file_name,
                &alpha_msg, &alpha_msg_ci, 
                &alpha_rma, &alpha_rma_ci, 
                &beta_msg, & beta_msg_ci,
                &beta_rma, & beta_rma_ci,
                &beta_memcpy, &beta_memcpy_ci,
                conf_level );
        }
        else {
            error = measure_bsp_params( pid, nprocs, repeat,
                pessimistic_word_size, optimistic_word_size,
                max_total_bytes, nhrels, min_niters, niters, ncomms, 
                sample_file_name[0] == '\0' ? NULL : sample_file_name,
                L, L_ci, g_pes, g_pes_ci, g_opt, g_opt_ci,
                conf_level );
        }

        t1 = bsp_time();

        if (error) {
            if (!pid)
                printf("# Interrupted: target precision was not achieved\n");
             break;
        }
            
        if (!pid && lincost_mode) {
            printf("#     Sampling time was %g seconds\n"
                    "# MSG:  alpha = %12.4g +/- %12.4g ;"
                    "beta =  %12.4g +/- %12.4g \n"
                    "# RMA:  alpha = %12.4g +/- %12.4g ;"
                    "beta =  %12.4g +/- %12.4g\n"
                    "# MEMCPY: beta= %12.4g +/- %12.4g\n",
                    t1-t0,
                    alpha_msg, alpha_msg_ci,
                    beta_msg, beta_msg_ci,
                    alpha_rma, alpha_rma_ci,
                    beta_rma, beta_rma_ci,
                    beta_memcpy, beta_memcpy_ci );
        }

        if (!pid && !lincost_mode) {
            printf("#     Sampling time was %g seconds\n", t1-t0);
            for ( i = 0; i < N_BSP_METHODS; ++i ) {
                printf( "# %5s :: L = %12.4g +/- %12.4g ;"
                        " g_opt =  %12.4g +/- %12.4g ;"
                        " g_pes = %12.4g +/- %12.4g\n",
                        bsp_method_labels[ i ],
                        L[i], L_ci[i], g_opt[i], g_opt_ci[i], g_pes[i], g_pes_ci[i]);
            }
        }
        speed_estimate = (t1-t0)/niters;
        pow_digits = int_pow(10,digits);

        if ( lincost_mode ) {
            grow_factor = 
                MAX( MAX( MAX( alpha_rma_ci * pow_digits  / alpha_rma
                         ,  alpha_msg_ci * pow_digits / alpha_msg) 
                    , beta_msg_ci * pow_digits / beta_msg )
                ,  beta_rma_ci * pow_digits / beta_rma );
        } else {
            grow_factor = 1.0;
            for ( i = 0; i < N_BSP_METHODS; ++i ) {
                grow_factor = MAX( MAX( MAX( grow_factor,
                                    L_ci[i] * pow_digits / L[i]) 
                                 ,  g_opt_ci[i] * pow_digits / g_opt[i])
                                 ,  g_pes_ci[i] * pow_digits / g_pes[i]);
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
           alpha_rma, beta_rma, n_hp_rma, alpha_msg, beta_msg, n_hp_msg );

        printf("\n# Save this output to a file, and use it with bsprun's\n"
                 "# --use-p2p-benchmark-file parameter\n");
    }

    if (!pid && !lincost_mode ) {
        for ( i = 0; i < N_BSP_METHODS; ++i ) {
            double w_opt = optimistic_word_size;
            double w_pes = pessimistic_word_size;
            double g_assymp 
                = (g_opt[i]*w_opt - g_pes[i]*w_pes) / (w_opt - w_pes );
            double req_overhead = g_pes[i] * w_pes - g_assymp * w_pes;

            g_assymp = MAX( 0.0, g_assymp );
            req_overhead = MAX( 0.0, req_overhead );
            printf(
               "export BSC_%s_L=\"%.2g\"\n"
               "export BSC_%s_G=\"%.2g\"\n"
               "export BSC_%s_O=\"%.2g\"\n",
                bsp_method_labels[i], L[i], 
                bsp_method_labels[i], g_assymp,
                bsp_method_labels[i], req_overhead );
        }
           

        printf("\n# Save this output to a file, and use it with bsprun's\n"
                 "# --use-p2p-benchmark-file parameter\n");
    }

exit:
    bsp_end();
    return EXIT_SUCCESS;
}
