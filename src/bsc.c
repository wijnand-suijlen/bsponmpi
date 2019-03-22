#include "bsc.h"
#include <limits.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
#include <assert.h>

#define MAX( a, b ) ( (a) < (b) ? (b) : (a) )
#define MIN( a, b ) ( (a) > (b) ? (b) : (a) )

static uint32_t log_2( uint32_t x )
{
   uint32_t a = (x & 0xffff0000)?16:0;
   uint32_t b = (x & 0xff00ff00)?8:0;
   uint32_t c = (x & 0xf0f0f0f0)?4:0;
   uint32_t d = (x & 0xCCCCCCCC)?2:0;
   uint32_t e = (x & 0xAAAAAAAA)?1:0;
   return a + b + c + d + e;
}

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
    enum { PUT, GET, EXEC } kind;
    union { 
        comm_request_t put; 
        comm_request_t get;
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
    req->kind = PUT;
    req->payload.put = r;
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
    req->kind = GET;
    req->payload.get = r;
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
                 if ( r->kind == PUT ) {
                    assert( r->payload.put.src_pid == bsp_pid() ) ;
                    bsp_put( r->payload.put.dst_pid, 
                            r->payload.put.src_addr, r->payload.put.dst_addr, 
                            r->payload.put.offset, r->payload.put.size );
                 } else if ( r->kind == GET ) {
                    assert( r->payload.get.dst_pid == bsp_pid() ) ;
                    bsp_get( r->payload.get.src_pid, r->payload.get.src_addr,
                            r->payload.get.offset, r->payload.get.dst_addr, 
                            r->payload.get.size );
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
    /** the number of members in this group */
    bsp_pid_t size;

    /** Array of size nprocs with the local pid for each global pid */
    bsp_pid_t * lid;

    /** the global PIDs of the group members */
    bsp_pid_t * gid;

    /** the memory offset this process in a remote arary */
    bsp_size_t * rof;

} group_t;

static void group_create( group_t * g , bsp_pid_t members )
{
    g->size = members;
    g->lid = calloc( bsp_nprocs(), sizeof(g->lid[0]));
    g->gid = calloc( members, sizeof(g->gid[0]));
    g->rof = calloc( members, sizeof(g->rof[0]));
    if (!g->lid || !g->rof || !g->gid )
        bsp_abort("bsc: insufficient memory for creating process group\n");
}

static void group_destroy( group_t * g )
{
    g->size = 0;
    free(g->lid);
    free(g->rof);
    free(g->gid);
    g->lid = NULL;
    g->rof = NULL;
    g->gid = NULL;
}

const bsc_group_t bsc_all = NULL;

bsc_group_t bsc_group_create_partition( unsigned color )
{
    bsp_pid_t i, j, P = bsp_nprocs(), s = bsp_pid();
    bsp_pid_t n_members = 0;
    bsp_pid_t rof = 0;
    unsigned * all_colors = calloc( P, sizeof(all_colors[0]) );
    group_t * g = calloc( 1, sizeof(*g) );
    if (!all_colors || !g)
        bsp_abort("bsc_group_create_partition: Insufficient memory");
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
            if (i == s ) rof = n_members;
            n_members += 1;
        }

    group_create( g, n_members );
    j = 0;
    for ( i = 0; i < P; ++i )
        if ( all_colors[i] == color ) {
            g->lid[ i ] = j;
            g->rof[ j ] = rof;
            g->gid[ j ] = i;
            j++;
        }
        else {
            g->lid[i] = -1;
        }

    bsp_pop_reg(all_colors);
    bsp_sync();
    free(all_colors);

    return g;
}

bsc_group_t bsc_group_create_neighbourhood( bsc_pid_t * neighbours, 
        bsc_size_t size )
{
    bsc_pid_t i, P = bsp_nprocs(), s = bsp_pid();
    bsc_pid_t n_members = 0;
    bsc_size_t * incoming = calloc( P, sizeof(incoming[0]));
    bsc_size_t * outgoing = calloc( P, sizeof(outgoing[0]));
    group_t * g = calloc( 1, sizeof(*g) );
    if (!incoming || !outgoing || !g)
        bsp_abort("bsc_group_create_neighbourhood: Insufficient memory");
    bsp_push_reg( incoming, P );
    bsp_sync();

    /* fill the communication matrix with a 1 for each connection */
    for ( i = 0 ; i < size ; ++i ) {
        bsc_size_t offset = i+1;
        bsp_put( neighbours[i], &offset, incoming, s, sizeof(offset) );
        outgoing[ neighbours[i] ] = i+1 ;
    }
    bsp_sync();

    if ( !outgoing[s] )
        bsp_abort("bsc_group_create_neighbourhood: Group must include calling process");

    /* take the intersection of incoming and outgoing */
    for ( i = 0; i < P; ++i ) {
        if (incoming[i] && !outgoing[i])  {
            bsp_abort("bsc_group_create_neighourhood: Processes %d thinks %d"
                    " is a neighbour, but that feeling is not mutual\n",
                    i, s );
        }
        else if (!incoming[i] && outgoing[i]) {
            bsp_abort("bsc_group_create_neighourhood: Processes %d thinks %d"
                    " is a neighbour, but that feeling is not mutual\n",
                    s, i );
        }
        else if (incoming[i] && outgoing[i])
            n_members += 1;
    }

    group_create( g, n_members );
    for ( i = 0; i < P; ++i )
        g->lid[i] = -1;

    for ( i = 0; i < size; ++i ) {
        g->lid[ neighbours[i] ] = i;
        g->rof[ i ] = incoming[ neighbours[i] ] - 1;
        g->gid[ i ] = neighbours[i];
    }

    bsp_pop_reg(incoming);
    bsp_sync();
    free(incoming);
    free(outgoing);

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
    bsc_pid_t i, n = g?g->size:bsp_nprocs() ;
    if ( bsp_pid() == root ) {
        for ( i = 0; i < n; ++i ) {
            bsc_pid_t m = g?g->gid[i]:i;
            const char * src_bytes = src;
            bsc_put( depends, m, src_bytes + i * size, dst, 0, size );
        }
    }
    return depends + 1;
}

bsc_step_t bsc_gather( bsc_step_t depends, 
                       bsc_pid_t root, bsc_group_t group,
                       const void * src, void *dst, bsc_size_t size )
{
    const group_t * g = group;
    bsc_pid_t i, n = g?g->size:bsp_nprocs() ;
    if ( bsp_pid() == root ) {
        for ( i = 0; i < n; ++i ) {
            bsc_pid_t m = g?g->gid[i]:i;
            char * dst_bytes = dst;
            bsc_get( depends, m, src, 0, dst_bytes + i * size, size );
        }
    }

    return depends + 1;
}

bsc_step_t bsc_allgather( bsc_step_t depends, bsc_group_t group,
                        const void * src, void * dst, bsc_size_t size )
{
    const group_t * g = group;
    bsc_pid_t i, n = g?g->size:bsp_nprocs() ;
    for ( i = 0; i < n; ++i ) {
        char * dst_bytes = dst;
        bsc_get( depends, g?g->gid[i]:i, src, 0, dst_bytes + i * size, size );
    }
    return depends + 1;
}

bsc_step_t bsc_alltoall( bsc_step_t depends, bsc_group_t group,
                      const void * src, void * dst, bsc_size_t size )
{
    const group_t * g = group;
    bsp_pid_t i, n = g?g->size:bsp_nprocs() ;
    const char * src_bytes = (const char *) src;
    for ( i = 0; i < n; ++i ) {
        bsc_pid_t m = g?g->gid[i]:i;
        bsc_size_t doff = g?g->rof[i]:bsp_pid();
        bsc_put( depends, m, src_bytes + i * size, dst, doff * size, size );
    }
    return depends + 1;
}

static bsc_step_t bcast_1phase( bsc_step_t depends,
                      bsp_pid_t root, bsc_group_t group, 
                      const void * src, void * dst, bsc_size_t size )
{
    const group_t * g = group;
    bsc_pid_t i, n = g?g->size:bsp_nprocs() ;
    if ( bsp_pid() == root ) {
        for ( i = 0 ; i < n; ++i) {
            bsc_put( depends, g?g->gid[i]:i, src, dst, 0, size );
        }
    }
    return depends + 1;
}


static bsc_step_t bcast_qtree( bsc_step_t depends,
                      bsp_pid_t root, bsc_group_t group, 
                      const void * src, void * dst, bsc_size_t size,
                      bsc_pid_t q )
{
    const group_t * g = group;
    bsc_pid_t i, j, n = g?g->size:bsp_nprocs() ;
    bsc_pid_t range, r = g?g->lid[root]:root;
    bsc_pid_t me = ((g?g->lid[ bsp_pid() ] : bsp_pid()) + n - r) %n;

    if ( me == 0 ) 
        memcpy( dst, src, size ); 

    for (range = 1; range < n; range *= q ) {
        if ( me < range ) {
            for ( i = 1 ; i < q; ++i) {
                j = i * range + me; 
                if (j < n) {
                    j = ( j + r ) % n;
                    bsc_put( depends, g?g->gid[j]:j, dst, dst, 0, size );
                }
            }
        }
        depends += 1;
    }
    return depends;
}

static bsc_step_t bcast_2phase( bsc_step_t depends,
                      bsp_pid_t root, bsc_group_t group, 
                      const void * src, void * dst, bsc_size_t size )
{
    const group_t * g = group;
    bsc_pid_t n = g?g->size:bsp_nprocs() ;
    depends = bsc_scatter( depends, root, group, src, dst, size / n );
    depends = bsc_allgather( depends, group, dst, dst, size / n );
    return depends;
}


bsc_step_t bsc_bcast( bsc_step_t depends,
                      bsp_pid_t root, bsc_group_t group, 
                      const void * src, void * dst, bsc_size_t size )
{
    const group_t * g = group;
    int i;
    bsc_pid_t n = g?g->size:bsp_nprocs() ;
    bsc_size_t rest = size % n;
    bsc_pid_t root_n = ceil( sqrt( (double) n ) );
    enum ALG { ONE_PHASE, TWO_PHASE, TWO_TREE, ROOT_TREE, N_ALGS } ;
    double costs[N_ALGS];
    enum ALG min_alg = ONE_PHASE ;
    double min_cost = HUGE_VAL;
    costs[ ONE_PHASE] = n * size * bsc_g() + bsc_L();
    costs[ TWO_PHASE] = (2*size + n * rest ) * bsc_g() + 2 * bsc_L();
    costs[ TWO_TREE]  = ((log_2(n-1)+1) * ( size * bsc_g() +  bsc_L() ) );
    costs[ ROOT_TREE] = 2*(root_n - 1) * size * bsc_g() + 2 * bsc_L();

    for ( i = 0; i < N_ALGS; ++i )
        if ( costs[i] < min_cost ) { 
            min_cost = costs[i];
            min_alg = (enum ALG) i;
        }

    switch (min_alg) {
        case ONE_PHASE:
            return bcast_1phase( depends, root, group, src, dst, size );

        case TWO_PHASE: {
            const char * src_bytes = src;
            char * dst_bytes = dst;
            bsc_size_t rest_offset = size - rest;
            bcast_1phase( depends, root, group, 
                    src_bytes + rest_offset, dst_bytes + rest_offset, rest );
            bcast_2phase( depends, root, group, src, dst, rest_offset );
            return depends;
        }

        case TWO_TREE: 
            return bcast_qtree( depends, root, group, src, dst, size, 2 );

        case ROOT_TREE:
            return bcast_qtree( depends, root, group, src, dst, size, root_n);

        case N_ALGS:
            bsp_abort("bsc_bcast: missing algorithm\n");
            return depends;
    }
    return depends;
}

static bsc_step_t reduce_1phase( bsc_step_t depends, 
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

static bsc_step_t reduce_qtree( bsc_step_t depends, 
        bsc_pid_t root, bsc_group_t group,
        const void * src, void * dst, void * tmp_space,
        bsc_reduce_t reducer, const void * zero,
        bsc_size_t nmemb, bsc_size_t size, bsc_pid_t q )
{
    const group_t * g = group;
    bsc_pid_t i, j, n = g?g->size:bsp_nprocs() ;
    bsc_pid_t range, r = g?g->lid[root]:root;
    bsc_pid_t me = ((g?g->lid[ bsp_pid() ] : bsp_pid()) + n - r) %n;
    char * tmp_bytes = tmp_space;

    (*reducer)( dst, zero, src, size * nmemb );
    for (range = 1; range < n; range *= q) {
        /* search smallest power of q >= n */
    }

    for (range = range/q; range >= 1; range /= q ) {
        if ( me < range ) {
            for ( i = 1 ; i < q; ++i) {
                j = i * range + me; 
                if (j < n) {
                    j = (j + r )% n;
                    bsc_get( depends, g?g->gid[j]:j, dst, 0, tmp_bytes + i*size, size );
                }
            }
            j = (me + r)%n;
            bsc_get( depends, g?g->gid[j]:j, dst, 0, tmp_space, size );
            bsc_exec_reduce( depends, reducer, dst, zero, tmp_space, size * q );
        }
        depends += 1;
    }
    return depends;

}



bsc_step_t bsc_reduce( bsc_step_t depends, 
        bsc_pid_t root, bsc_group_t group,
        const void * src, void * dst, void * tmp_space,
        bsc_reduce_t reducer, const void * zero,
        bsc_size_t nmemb, bsc_size_t size )
{
    const group_t * g = group;
    int i;
    bsc_pid_t n = g?g->size:bsp_nprocs() ;
    bsc_pid_t root_n = ceil( sqrt( (double) n ) );
    enum ALG { ONE_PHASE, TWO_TREE, ROOT_TREE, N_ALGS } ;
    double costs[N_ALGS];
    enum ALG min_alg = ONE_PHASE ;
    double min_cost = HUGE_VAL;
    costs[ ONE_PHASE] = n * size * bsc_g() + bsc_L();
    costs[ TWO_TREE]  = ((log_2(n-1)+1) * ( size * bsc_g() +  bsc_L() ) );
    costs[ ROOT_TREE] = 2*(root_n - 1) * size * bsc_g() + 2 * bsc_L();

    for ( i = 0; i < N_ALGS; ++i )
        if ( costs[i] < min_cost ) { 
            min_cost = costs[i];
            min_alg = (enum ALG) i;
        }

    switch (min_alg) {
        case ONE_PHASE:
            return reduce_1phase( depends, root, group, src, dst, tmp_space,
                      reducer, zero, nmemb, size );

        case TWO_TREE: 
            return reduce_qtree( depends, root, group, src, dst, tmp_space,
                      reducer, zero, nmemb, size, 2 );

        case ROOT_TREE: 
            return reduce_qtree( depends, root, group, src, dst, tmp_space,
                      reducer, zero, nmemb, size, root_n );
        case N_ALGS:
            bsp_abort("bsc_reduce: missing algorithm\n");
            return depends;
    }
    return depends;
}

static bsc_step_t allreduce_1phase( bsc_step_t depends, bsc_group_t group,
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


bsc_step_t bsc_allreduce( bsc_step_t depends, bsc_group_t group,
        const void * src, void * dst, void * tmp_space,
        bsc_reduce_t reducer, const void * zero,
        bsc_size_t nmemb, bsc_size_t size )
{
    const group_t * g = group;
    int i;
    bsc_pid_t n = g?g->size:bsp_nprocs(), root=0 ;
    bsc_pid_t root_n = ceil( sqrt( (double) n ) );
    enum ALG { ONE_PHASE, TWO_TREE, ROOT_TREE, N_ALGS } ;
    double costs[N_ALGS];
    enum ALG min_alg = ONE_PHASE ;
    double min_cost = HUGE_VAL;
    costs[ ONE_PHASE] = n * size * bsc_g() + bsc_L();
    costs[ TWO_TREE]  = 2*((log_2(n-1)+1) * ( size * bsc_g() +  bsc_L() ) );
    costs[ ROOT_TREE] = 4*(root_n - 1) * size * bsc_g() + 4 * bsc_L();

    for ( i = 0; i < N_ALGS; ++i )
        if ( costs[i] < min_cost ) { 
            min_cost = costs[i];
            min_alg = (enum ALG) i;
        }

    switch (min_alg) {
        case ONE_PHASE:
            return allreduce_1phase( depends, group, src, dst, tmp_space,
                      reducer, zero, nmemb, size );

        case TWO_TREE: 
            depends = reduce_qtree( depends, root, group, src, dst, tmp_space,
                      reducer, zero, nmemb, size, 2 );
            return bcast_qtree( depends, root, group, dst, dst, size, 2 );

        case ROOT_TREE: 
            depends = reduce_qtree( depends, root, group, src, dst, tmp_space,
                      reducer, zero, nmemb, size, root_n );
            return bcast_qtree( depends, root, group, dst, dst, size, root_n );
        
        case N_ALGS:
            bsp_abort("bsc_allreduce: missing algorithm\n");
            return depends;
    }
    return depends;
}

bsc_step_t bsc_scan( bsc_step_t depends, bsc_group_t group,
        const void * src, void * dst, 
        void * tmp_space1, void * tmp_space2,
        bsc_reduce_t reducer, const void * zero,
        bsc_size_t nmemb, bsc_size_t size )
{
    const group_t * g = group;
    bsp_pid_t n = g?g->size:bsp_nprocs() ;
    bsc_size_t o = g->lid[bsp_pid()];
    void * local_last = (char *) dst + (nmemb-1) * size;
    const void * start = o == 0? zero : (const char * ) tmp_space2 + (o-1) * size;

    /* compute local prefix sum */
    (*reducer)( dst, zero, src, size * nmemb );

    /* get all last elements of on all processes */
    depends = bsc_allgather( depends, group, local_last, tmp_space1, size );

    /* compute the prefix sum of all sums in tmp_space2 */
    bsc_exec_reduce( depends-1, reducer, tmp_space2, zero, tmp_space1, size * n );

    /* Compute the prefix sum of src and use tmp_space2[o-1] as start */
    bsc_exec_reduce( depends-1, reducer, dst, start, src, size * nmemb );

    return depends;
}




