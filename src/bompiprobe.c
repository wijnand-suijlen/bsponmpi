#ifdef HAS_CLOCK_GETTIME
  #define _POSIX_C_SOURCE 199309L
  #include <time.h>
#endif

#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <signal.h>
#include "util.h"

#ifdef HAS_CLOCK_GETTIME
double get_time() { 
  struct timespec now;
  clock_gettime( CLOCK_MONOTONIC, &now );
  return now.tv_sec + 1e-9 * now.tv_nsec;
}
#else
double get_time() { return MPI_Wtime() }
#endif


double measure_hrel_rma( MPI_Comm comm, MPI_Win win, 
        char * sendbuf, int h, int size ) 
{
    double t0, t1;
    int i, pid, nprocs;
    MPI_Comm_size( comm, &nprocs);
    MPI_Comm_rank( comm, &pid );
    MPI_Win_fence( 0, win );
    MPI_Barrier( comm );
    t0 = get_time();
    for ( i = 0; i < h; ++i ) {
        MPI_Put( sendbuf + i*size, size, MPI_BYTE,
                (i+pid) % nprocs, pid*size, size, MPI_BYTE,
                win );
    }
    MPI_Win_fence( 0, win );
    t1 = get_time() - t0;
    MPI_Allreduce( &t1, &t0, 1, MPI_DOUBLE, MPI_MAX, comm );
    return t0;
}

double measure_hrel_msg( MPI_Comm comm, 
        char * sendbuf, char * recvbuf, MPI_Request * reqs, 
        int h, int size ) 
{
    double t0, t1;
    int i, pid, nprocs, tag = 0;
    MPI_Comm_size( comm, &nprocs );
    MPI_Comm_rank( comm, &pid );

    MPI_Barrier( comm );
    t0 = get_time();

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

    t1 = get_time() - t0;
    MPI_Allreduce( &t1, &t0, 1, MPI_DOUBLE, MPI_MAX, comm );
    return t0;
}


typedef struct sample {
    double time;   /* time per msg */
    int    h;      /* # msgs */
    int    msg_size; /* size per msg */
    int    comm;   /* topology */
    int    method; /* 0 == rma, 1 == msg */
} sample_t;

void permute( size_t * rng, sample_t * samples, int n_samples )
{
    while (n_samples > 1 ){
        int i = n_samples * rand_next(rng);
        sample_t tmp = samples[i];
        samples[i] = samples[n_samples-1];
        samples[n_samples-1] = tmp;
        n_samples -= 1;
    }
}

int cmp_sample( const void * a, const void * b )
{
    const sample_t * x = a, * y = b;

    int i = x->method - y->method;
    int j = x->msg_size - y->msg_size;

    return i?i:j;
}

int cmp_double( const void * a, const void * b )
{
    const double *  x = a, *y = b;
    return x < y ? -1 : x > y ? 1 : 0;
}

void measure_lincost( int pid, int nprocs, 
        int msg_size, int niters, int ncomms,
        double * alpha_msg, double * alpha_msg_ci, 
        double * alpha_rma, double * alpha_rma_ci, 
        double * beta_msg, double * beta_msg_ci,
        double * beta_rma, double * beta_rma_ci,
        double conf_level )
{
    char * sendbuf, * recvbuf;
    sample_t * samples;
    size_t rng = 0;
    int i, j, k;
    MPI_Request * reqs = 0; 
    MPI_Win * wins;  MPI_Comm * comms;
    const int n_avgs = 10.0 / conf_level;
    double * avgs;
    int o[5];  /* two methods x two msg sizes + 1 */
    double T[4], T_ci[4]; /* two methods x two msg sizes */

    /* allocate memory */
    sendbuf = calloc( (long) 2 * msg_size * nprocs, 1 );
    recvbuf = calloc( (long) 2 * msg_size * nprocs, 1 );
    samples =  calloc( niters, sizeof(samples[0]) );
    reqs = calloc( 2*nprocs, sizeof(MPI_Request) );
    wins = calloc( ncomms, sizeof(MPI_Win) );
    comms = calloc( ncomms, sizeof(MPI_Comm) ); 
    avgs = calloc( n_avgs, sizeof(avgs[0]) );

    if (!sendbuf || !recvbuf || !samples || !reqs
            || !wins || !comms || !avgs  ) {
        fprintf(stderr, "Insufficient memory. For a benchmark on "
                "%d processes I require 2x %ld bytes of memory "
                "for send and receive buffers, %ld bytes for "
                "storing the measurements, %ld bytes for "
                "the request queue, %ld + %ld bytes for "
                "communication infrastructure, and %ld bytes for "
                "analysis.\n",
                nprocs, (long) 2*msg_size * nprocs,
                (long) niters * (long) sizeof(samples[0]),
                (long) 2*nprocs*sizeof(MPI_Request),
               ncomms * (long) sizeof(MPI_Win), 
               ncomms * (long) sizeof(MPI_Comm), 
               (long) n_avgs * (long) sizeof(avgs[0]) );
        MPI_Abort( MPI_COMM_WORLD, EXIT_FAILURE );
    }

    memset( sendbuf, 1, (long) 2*msg_size*nprocs );
    memset( recvbuf, 2, (long) 2*msg_size*nprocs );

    /* Create random topologies */
    MPI_Comm_dup( MPI_COMM_WORLD, &comms[0] );
    MPI_Win_create( recvbuf, (long) 2*msg_size * nprocs, 1,
            MPI_INFO_NULL, comms[0], & wins[0] );

    for ( i = 1; i < ncomms; ++i ) {
        int key = 0;
        /* Note: the random key can be chosen much quicker */
        for ( j = 0 ; j < nprocs; ++j ) {
            double r = rand_next( &rng ) ;
            if (j == pid ) key = 10*nprocs * r;
        }
    
        MPI_Comm_split( comms[0], 0, key, &comms[i] );
        MPI_Win_create( recvbuf, (long) 2*msg_size * nprocs, 1,
            MPI_INFO_NULL, comms[i], & wins[i] );
    }

    /* Create measurement points */
    for ( i = 0; i < niters; ++i ) 
        samples[i].h = (1 + i % nprocs);

    permute( &rng, samples, niters );

    for ( i = 0; i < niters; ++i )
       samples[i].comm = i % ncomms; 

    permute( &rng, samples, niters );

    for ( i = 0 ; i < niters; ++i )
        samples[i].msg_size = i < niters/2 ? msg_size : 2*msg_size;

    permute( &rng, samples, niters );

    for ( i = 0 ; i < niters; ++i )
        samples[i].method = i < niters/2 ;

    permute( &rng, samples, niters );


    /* Warm up */
    for ( i = 0; i < ncomms; ++i ) {
        measure_hrel_msg( comms[i], 
                sendbuf, recvbuf, reqs, nprocs, 2*msg_size );
        measure_hrel_rma( comms[i], wins[i], 
                sendbuf, nprocs, 2*msg_size );
    }

    /* And now the measurements */
    for ( i = 0; i < niters; ++i ) {
        double t;
        int h    = samples[i].h;
        int size = samples[i].msg_size;
        MPI_Comm comm = comms[ samples[i].comm ];
        MPI_Win win = wins[ samples[i].comm ];

        if ( samples[i].method )
            t = measure_hrel_msg( comm, sendbuf, recvbuf, reqs, h, size );
        else
            t = measure_hrel_rma( comm, win, sendbuf, h, size );

        samples[i].time = t / h;
    } 

    /* Analysis */

    /* sort on (msg_size, method) */
    qsort( samples, niters, sizeof(sample_t), cmp_sample );
    
    o[0] = 0;
    i = 0;
    for ( k = 0; k < 4; ++k ) {
        for ( ; i < niters; ++i )
            if (samples[i].msg_size != (k%2 + 1)*msg_size
                || samples[i].method != (k/2) )
                break;

        o[k+1] = i;
    }
    
    for ( k = 0; k < 4; ++k ) {
        int n = o[k+1] - o[k]; 
        int l = 0;
        double total_sum = 0.0;
        double left, right;
        /* permute to randomly sample over time */
        permute( &rng, samples + o[k], n); 
        for ( i = 0; i < n_avgs; i++) {
            int m = (n + n_avgs - i - 1)/n_avgs;
            double sum = 0.0;
            for ( j = 0; j < m; ++j ) {
                sum += samples[o[k]+l].time;
                l++;
            }
            avgs[i] = sum / m;
            total_sum += sum;
        }

        T[k] = total_sum / n;

        /* generate CFD of averages */
        qsort( avgs, n_avgs, sizeof(double), cmp_double );

        left = avgs[ (int) (n_avgs * conf_level * 0.5) ];
        right = avgs[ (int) (n_avgs * ( 1.0 - conf_level * 0.5)) ];

        /* compute the confidence interval */
        T_ci[k] = MAX( fabs(T[k] - left),  fabs(right - T[k]) );
    }

    /* saving results */
    *alpha_rma = 2*T[0] - T[1];
    *alpha_rma_ci = fabs(2*T_ci[0] + T_ci[1] );

    *alpha_msg = 2*T[2] - T[3];
    *alpha_msg_ci = fabs(2*T_ci[2] + T_ci[3] );

    *beta_rma = (T[1] - T[0])/msg_size;
    *beta_rma_ci = (T_ci[1] + T_ci[0]) / msg_size;

    *beta_msg = (T[3] - T[2])/msg_size;
    *beta_msg_ci = (T_ci[3] + T_ci[2]) / msg_size;

    /* free resources */
    for ( i = 0; i < ncomms; ++i ) {
        MPI_Win_free( &wins[i] );
        MPI_Comm_free( &comms[i] );
    }

    free( sendbuf );
    free( recvbuf );
    free( samples );
    free( reqs );
    free( comms );
    free( wins );
    free( avgs );
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
        double * conf_level, int * digits )
{
    int i;
    const char * val = NULL;

    for ( i = 1; i < argc; ++i ) {
        if ( get_param( argv[i], "--help") ) {
            fprintf(stderr, "Usage %s [--ncomms=<number>] "
                    "[--msg-size=<size>] [--conf-level=<percentage>] "
                    "[--digits=<number>]\n", argv[0]);
            return 1;
        }
        else if ( (val = get_param( argv[i], "--ncomms=" ) ) ) {
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
        else
            fprintf(stderr, "Ignoring unrecognized parameter '%s'\n",
                    argv[i]);
    }
    return 0;
}

int main( int argc, char ** argv )
{
    int nprocs, pid;
    double alpha_rma=1e+4, alpha_msg=1e+4;
    double alpha_rma_ci=1e+4, alpha_msg_ci=1e+4;
    double beta_rma=1.0, beta_msg=1.0;
    double beta_rma_ci=1.0, beta_msg_ci=1.0;
    double conf_level = 0.95;
    int niters = 1000;
    int ncomms = 10;
    int msg_size = 10000;
    int digits = 1;
    int error = 0;
    MPI_Init( &argc, & argv );

    MPI_Comm_rank( MPI_COMM_WORLD, &pid );
    MPI_Comm_size( MPI_COMM_WORLD, &nprocs );

    if (!pid) error = parse_params( argc, argv, 
            &ncomms, &msg_size, &conf_level, &digits );

    MPI_Bcast( &error, 1, MPI_INT, 0, MPI_COMM_WORLD );

    if (error) { 
        MPI_Finalize();
        return EXIT_FAILURE;
    }

    MPI_Bcast( &ncomms, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast( &msg_size, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast( &conf_level, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast( &digits, 1, MPI_INT, 0, MPI_COMM_WORLD);

    if (!pid) printf(
           "# Number of MPI processes= %d\n"
           "# Typical message size   = %d\n"
           "# Confidence level       = %.2f%%\n"
           "# Significant digits     = %d\n"
           "# Random process layouts = %d\n",
           nprocs, msg_size, 100.0*conf_level, digits, ncomms );


    do { 
        double t0, t1;
        niters *= 10;

        if (!pid)
            printf("#     Samples = %d\n", niters );

        t0 = get_time();
        measure_lincost( pid, nprocs, 
            msg_size,  niters, ncomms,
            &alpha_msg, &alpha_msg_ci, 
            &alpha_rma, &alpha_rma_ci, 
            &beta_msg, & beta_msg_ci,
            &beta_rma, & beta_rma_ci,
            conf_level );
        t1 = get_time();

        if (!pid) {
            printf("#     Sampling time =  %g seconds\n"
                    "# MSG: alpha = %12.4g +/- %12.4g ;"
                    "beta =  %12.4g +/- %12.4g \n"
                    "# RMA: alpha = %12.4g +/- %12.4g ;"
                    "beta =  %12.4g +/- %12.4g\n",
                    t1-t0,
                    alpha_msg, alpha_msg_ci,
                    beta_msg, beta_msg_ci,
                    alpha_rma, alpha_rma_ci,
                    beta_rma, beta_rma_ci );

        }
    } while ( alpha_rma_ci * int_pow(10,digits) > alpha_rma
         || alpha_msg_ci * int_pow(10, digits) > alpha_msg
         || beta_msg_ci * int_pow(10, digits) > beta_msg
         || beta_rma_ci * int_pow(10, digits) > beta_rma );

    if (!pid) {
        printf(
           "if [ x${BSPONMPI_A2A_METHOD} = xrma ]; then\n"
           "   export BSPONMPI_P2P_LATENCY=\"%.2g\"\n"
           "   export BSPONMPI_P2P_MSGGAP=\"%.2g\"\n"
           "elif [ x${BSPONMPI_A2A_METHOD} = xmsg ]; then\n"
           "   export BSPONMPI_P2P_LATENCY=\"%.2g\"\n"
           "   export BSPONMPI_P2P_MSGGAP=\"%.2g\"\n"
           "fi\n", 
           alpha_rma, beta_rma, alpha_msg, beta_msg);

        printf("\n# Save this output to a file, and use it with bsprun's\n"
                 "# --use-p2p-benchmark-file parameter\n");
    }

    MPI_Finalize();
    return EXIT_SUCCESS;
}
