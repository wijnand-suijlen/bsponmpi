#include "bsc.h"
#include <limits.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#define MAX( a, b ) ( (a) < (b) ? (b) : (a) )
#define MIN( a, b ) ( (a) > (b) ? (b) : (a) )

const bsc_step_t bsc_flush = UINT_MAX;
const bsc_step_t bsc_start = 0;

typedef struct {
    const void * src_addr;
    void *dst_addr;
    bsp_pid_t src_pid, dst_pid;
    bsp_size_t size, offset;
} comm_request_t ;

typedef struct {
    bsc_reduce_t func;
    void * a;
    const void * a0;
    const void * xs;
    bsc_size_t size;
} exec_request_t;

typedef struct request {
    enum { COMM, EXEC } kind;
    union { 
        comm_request_t comm; 
        exec_request_t exec;
    } payload;
} request_t;

typedef struct state { 
    bsc_step_t current;
    unsigned   horizon;
    unsigned   max_requests_per_step;
    unsigned   queue_start;
    unsigned   * n_requests;
    request_t  * queue;
    unsigned   total_n_outstanding_requests;
    int        * flushed;
} state_t;

static state_t s_bsc;

static request_t * new_request( bsc_step_t step ) 
{
    unsigned i, j, delay;
    request_t * entry;

    if ( step < s_bsc.current )
        bsp_abort("bsc: a flush sync is required before posting"
                  " another request\n");

    delay = step - s_bsc.current;
    if ( delay >= s_bsc.horizon || 
          s_bsc.n_requests[ delay ] >= s_bsc.max_requests_per_step ) {

        /* new memory allocation is necessary */
        bsc_step_t new_horizon;
        unsigned * new_n_requests;
        unsigned new_max_requests;
        request_t * new_queue;

        new_horizon = s_bsc.horizon;
        if (delay >= s_bsc.horizon) {
            new_horizon = MAX( 2*s_bsc.horizon, delay + 1);
        }

        if ( new_horizon <= delay )
            bsp_abort("bsc: cannot plan that far in the future"
                      "due to integer overflow\n");
        
        new_n_requests = calloc( new_horizon, sizeof(unsigned));
        if ( ! new_n_requests ) 
            bsp_abort("bsc: insufficient memory; "
                      "cannot allocate %u x %u bytes\n",
                new_horizon, (unsigned) sizeof(unsigned) );

        for ( i = 0; i < s_bsc.horizon; ++i ) {
            j = (s_bsc.queue_start + i) % s_bsc.horizon;
            new_n_requests[ i ] = s_bsc.n_requests[j];
        }

        new_max_requests = s_bsc.max_requests_per_step; 
        if (new_max_requests <= new_n_requests[delay]) 
            new_max_requests = MAX( 2 * s_bsc.max_requests_per_step, 
                    new_n_requests[delay] + 1  );

        if ( new_max_requests <= new_n_requests[delay] )
            bsp_abort("bsc: cannot plan that far in the future"
                      "due to integer overflow\n");
     
        new_queue = calloc( (size_t) new_horizon * new_max_requests, 
                            sizeof(request_t) );

        for ( i = 0; i < s_bsc.horizon; ++i ) {
            for ( j = 0; j < s_bsc.max_requests_per_step; ++j ){
                unsigned a = i * new_max_requests + j;
                unsigned b = (i + s_bsc.queue_start) % s_bsc.horizon 
                              * s_bsc.max_requests_per_step + j;
                new_queue[ a ] = s_bsc.queue[ b ];
            }
        }

        free( s_bsc.n_requests );
        free( s_bsc.queue );
        s_bsc.n_requests = new_n_requests;
        s_bsc.queue = new_queue;
        s_bsc.horizon = new_horizon;
        s_bsc.max_requests_per_step = new_max_requests;
        s_bsc.queue_start = 0;
    }

    entry = s_bsc.queue + 
        (delay + s_bsc.queue_start) % s_bsc.horizon * 
        s_bsc.max_requests_per_step + s_bsc.n_requests[delay] ;
    
    s_bsc.n_requests[delay] += 1;
    s_bsc.total_n_outstanding_requests += 1;

    return entry;
}


bsc_step_t bsc_put( bsc_step_t depends, bsc_pid_t dst_pid, 
        const void * src, void * dst, bsc_size_t offset, bsc_size_t size)
{
    request_t * req;
    comm_request_t r;
    r.src_addr = src;
    r.dst_addr = dst;
    r.src_pid = bsp_pid();
    r.dst_pid = dst_pid;
    r.size = size;
    r.offset = offset;
        
    req = new_request( depends );
    req->kind = COMM;
    req->payload.comm = r;
    return depends + 1;
}

bsc_step_t bsc_get( bsc_step_t depends, bsc_pid_t src_pid, 
        const void * src, bsc_size_t offset, void * dst, bsc_size_t size )
{
    request_t * req;
    comm_request_t r;
    r.src_addr = src;
    r.dst_addr = dst;
    r.dst_pid = bsp_pid();
    r.src_pid = src_pid;
    r.size = size;
    r.offset = offset;
        
    req = new_request( depends );
    req->kind = COMM;
    req->payload.comm = r;
    return depends + 1;
}

bsc_step_t bsc_exec_reduce( bsc_step_t depends, 
        bsc_reduce_t reducer, void * a, const void * a0,
        const void * xs, bsc_size_t size  )
{
    request_t * req;
    exec_request_t r;
    r.func = reducer;
    r.a= a;
    r.a0 = a0;
    r.xs = xs;
    r.size = size;
    
    req = new_request( depends );
    req->kind = EXEC;
    req->payload.exec = r;
    return depends + 1;
}

bsc_step_t bsc_current(void)
{
    return s_bsc.current;
}

bsc_step_t bsc_sync( bsc_step_t until )
{
    unsigned i;
    if ( !s_bsc.flushed ) {
        /* do some initialization */
        s_bsc.flushed = calloc( bsp_nprocs(), sizeof(s_bsc.flushed[0]));
        bsp_push_reg( s_bsc.flushed, 
                bsp_nprocs() * sizeof(s_bsc.flushed[0]) );
        bsp_sync();
    }

    while ( s_bsc.current <= until || until == bsc_flush ) {
        if (s_bsc.horizon > 0) {
            for ( i = 0; i < s_bsc.n_requests[s_bsc.queue_start]; ++i ) {
                request_t * r = s_bsc.queue + 
                    s_bsc.queue_start * s_bsc.max_requests_per_step + i ;
                 if ( r->kind == COMM ) {
                    if ( r->payload.comm.src_pid == bsp_pid() ) {
                        bsp_put( r->payload.comm.dst_pid, 
                                r->payload.comm.src_addr, r->payload.comm.dst_addr, 
                                r->payload.comm.offset, r->payload.comm.size );
                    }
                    else
                    {
                        bsp_get( r->payload.comm.src_pid, r->payload.comm.src_addr,
                                r->payload.comm.offset, r->payload.comm.dst_addr, 
                                r->payload.comm.size );
                    }
                }
            }
            s_bsc.total_n_outstanding_requests 
                -= s_bsc.n_requests[ s_bsc.queue_start ];
        }
        if (until == bsc_flush && s_bsc.total_n_outstanding_requests == 0)
        {
            unsigned P = bsp_nprocs(), sz = sizeof(s_bsc.flushed[0]);
            for ( i = 0; i < P; ++i ) {
                int flushed = 1;
                bsp_put( i, &flushed, s_bsc.flushed, bsp_pid()*sz, sz );
            }
        }
        bsp_sync();
        if (s_bsc.horizon > 0) {
            for ( i = 0; i < s_bsc.n_requests[s_bsc.queue_start]; ++i ) {
                request_t * r = s_bsc.queue + 
                    s_bsc.queue_start * s_bsc.max_requests_per_step + i ;
                if ( r->kind == EXEC ) {
                    bsc_reduce_t f = r->payload.exec.func;
                    void * a = r->payload.exec.a;
                    const void * a0 = r->payload.exec.a0 ;
                    const void * xs = r->payload.exec.xs;
                    bsc_size_t size = r->payload.exec.size;
                    (*f)( a, a0, xs, size );
                }
            }
            s_bsc.n_requests[ s_bsc.queue_start ] = 0;
            s_bsc.queue_start = (s_bsc.queue_start + 1) % s_bsc.horizon;
        }
        s_bsc.current += 1;

        if ( until == bsc_flush ) {
            int flushed = 1;
            unsigned P = bsp_nprocs();
            for ( i = 0 ; i < P; ++i ) {
                flushed = flushed && s_bsc.flushed[i];
            }
            if (flushed) { 
                memset( s_bsc.flushed, 0, 
                        sizeof(s_bsc.flushed[0]) * bsp_nprocs() );
                break;
            }
        }
    }
    if (s_bsc.total_n_outstanding_requests == 0 )
        s_bsc.current = 0;

    return s_bsc.current;
}

double bsc_g()
{
    return 1.0;
}

double bsc_L()
{
    return 1.0e+4;
}


typedef struct group {
    bsp_pid_t size;
    bsp_pid_t local_pid;
    bsp_pid_t * members;
} group_t;

static void group_create( group_t * g , bsp_pid_t local_pid, bsp_pid_t members )
{
    g->size = members;
    g->local_pid = local_pid;
    g->members = calloc( members, sizeof(g->members[0]));
    if (!g->members)
        bsp_abort("bsc: insufficient memory for creating process group\n");
}

static void group_destroy( group_t * g )
{
    g->size = 0;
    g->local_pid = 0;
    free(g->members);
    g->members = NULL;
}

const bsc_group_t bsc_all = NULL;

bsc_group_t bsc_group_create( unsigned color )
{
    bsp_pid_t i, j, P = bsp_nprocs(), s = bsp_pid();
    bsp_pid_t n_members = 0;
    bsp_pid_t local_pid = 0;
    unsigned * all_colors = calloc( P, sizeof(all_colors[0]) );
    group_t * g = calloc( 1, sizeof(*g) );
    if (!all_colors || !g)
        bsp_abort("bsc_group_create: Insufficient memory");
    bsp_push_reg( all_colors, P * sizeof(all_colors[0]) );
    bsp_sync();

    /* all-gather the colors */
    for ( i = 0 ; i < P; ++i ) {
        bsp_put( i, &color, all_colors, sizeof(all_colors[0])*s,
                sizeof(all_colors[0]) );
    }
    bsp_sync();

    /* count number of members */
    for ( i = 0; i < P; ++i )
        if ( all_colors[i] == color ) {
            if (i == s ) local_pid = n_members;
            n_members += 1;
        }

    group_create( g, local_pid, n_members );
    j = 0;
    for ( i = 0; i < P; ++i )
        if ( all_colors[i] == color )
            g->members[ j++ ] = i;

    bsp_pop_reg(all_colors);
    bsp_sync();
    free(all_colors);

    return g;
}

void bsc_group_destroy( bsc_group_t group )
{
    group_t * g = group;
    group_destroy( g );
    free( g );
}

bsc_step_t bsc_scatter( bsc_step_t depends,
                        bsc_pid_t root, bsc_group_t group,
                       const void * src, void *dst, bsc_size_t size )
{
    const group_t * g = group;
    bsp_pid_t i, n = g?g->size:bsp_nprocs() ;
    if ( bsp_pid() == root ) {
        const char * src_bytes = src;
        for ( i = 0; i < n; ++i ) {
            bsc_put( depends, g?g->members[i]:i, 
                    src_bytes + size * i, dst, 0, size );
        }
    }
    return depends+1;
}

bsc_step_t bsc_gather( bsc_step_t depends, 
                       bsp_pid_t root, bsc_group_t group,
                       const void * src, void *dst, bsc_size_t size )
{
    const group_t * g = group;
    bsp_pid_t s = g?g->local_pid:bsp_pid();
    return bsc_put( depends, root, src, dst, size*s, size);
}

bsc_step_t bsc_allgather( bsc_step_t depends, bsc_group_t group,
                        const void * src, void * dst, bsc_size_t size )
{
    const group_t * g = group;
    bsp_pid_t i, n = g?g->size:bsp_nprocs() ;
    bsp_pid_t s = g?g->local_pid:bsp_pid();
    for ( i = 0; i < n; ++i ) {
        bsc_put( depends, g?g->members[i]:i, src, dst, s * size, size );
    }
    return depends + 1;
}

bsc_step_t bsc_alltoall( bsc_step_t depends, bsc_group_t group,
                      const void * src, void * dst, bsc_size_t size )
{
    const group_t * g = group;
    bsp_pid_t i, n = g?g->size:bsp_nprocs() ;
    bsp_pid_t s = g?g->local_pid:bsp_pid();
    const char * src_bytes = (const char *) src;
    for ( i = 0; i < n; ++i ) {
        bsc_put( depends, g?g->members[i]:i, src_bytes + i * size, 
                dst, s * size, size );
    }
    return depends + 1;
}


bsc_step_t bsc_bcast( bsc_step_t depends,
                      bsp_pid_t root, bsc_group_t group, 
                      const void * src, void * dst, bsc_size_t size )
{
    const group_t * g = group;
    bsp_pid_t i, n = g?g->size:bsp_nprocs() ;
    double one_phase = n * size * bsc_g() + bsc_L();
    double two_phase = (2*size + n * (size % n)) * bsc_g() + 2 * bsc_L();

    if ( one_phase <= two_phase ) {
        if ( bsp_pid() == root ) {
            for ( i = 0; i <  n; ++i ) 
                bsc_put( depends, g?g->members[i]:i, src, dst, 0, size );
        }
        return depends + 1;
    }
    else {
        const char * src_bytes = src;
        bsc_size_t rest = size % n;
        bsc_size_t rest_offset = size - rest;
        if ( bsp_pid() == root ) {
            for ( i = 0; i < n; ++i ) 
                bsc_put( depends, g?g->members[i]:i, 
                        src_bytes + rest_offset, dst, rest_offset, rest );
        }
        depends = bsc_scatter( depends, root, group, src, dst, size / n );
        depends = bsc_allgather( depends, group, dst, dst, size / n );
        return depends;
    }
}


bsc_step_t bsc_reduce( bsc_step_t depends, 
        bsc_pid_t root, bsc_group_t group,
        const void * src, void * dst, void * tmp_space,
        bsc_reduce_t reducer, const void * zero,
        bsc_size_t nmemb, bsc_size_t size )
{
    const group_t * g = group;
    (*reducer)( dst, zero, src, size * nmemb );
    depends = bsc_gather( depends, root, group, dst, tmp_space, size );
    if ( root == bsp_pid() )
        bsc_exec_reduce( depends-1, reducer, dst, zero, tmp_space, 
                size * (g?g->size:bsp_nprocs())  );

    return depends;
}

bsc_step_t bsc_allreduce( bsc_step_t depends, bsc_group_t group,
        const void * src, void * dst, void * tmp_space,
        bsc_reduce_t reducer, const void * zero,
        bsc_size_t nmemb, bsc_size_t size )
{
    const group_t * g = group;
    (*reducer)( dst, zero, src, size * nmemb );
    depends = bsc_allgather( depends, group, dst, tmp_space, size );
    return bsc_exec_reduce( depends-1, reducer, dst, zero, tmp_space, 
                size * (g?g->size:bsp_nprocs())  );
}

bsc_step_t bsc_scan( bsc_step_t depends, bsc_group_t group,
        const void * src, void * dst, 
        void * tmp_space1, void * tmp_space2,
        bsc_reduce_t reducer, const void * zero,
        bsc_size_t nmemb, bsc_size_t size )
{
    const group_t * g = group;
    bsp_pid_t n = g?g->size:bsp_nprocs() ;
    bsp_pid_t s = g?g->local_pid:bsp_pid();
    void * local_last = (char *) dst + (nmemb-1) * size;
    const void * start = s == 0? zero : (const char * ) tmp_space2 + (s-1) * size;

    /* compute local prefix sum */
    (*reducer)( dst, zero, src, size * nmemb );

    /* get all last elements of on all processes */
    depends = bsc_allgather( depends, group, local_last, tmp_space1, size );

    /* compute the prefix sum of all sums in tmp_space2 */
    bsc_exec_reduce( depends-1, reducer, tmp_space2, zero, tmp_space1, size * n );

    /* Compute the prefix sum of src and use tmp_space2[s-1] as start */
    bsc_exec_reduce( depends-1, reducer, dst, start, src, size * nmemb );

    return depends;
}




