#include "bsc.h"
#include "util.h"
#include <limits.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
#include <assert.h>

typedef struct group group_t;

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

typedef struct coll_request {
    const bsc_collective_t * coll;
    bsc_pid_t          root;
    group_t *          group;
    bsc_coll_params_t  params;
    struct coll_request * next;
} coll_request_t;


typedef struct request {
    enum { PUT, GET, EXEC, COLL } kind;
    union { 
        comm_request_t put; 
        comm_request_t get;
        exec_request_t exec;
        coll_request_t coll;
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

    /* Collective coalescing table. The items are the
     * coll_request_t structs. Keys are made up of
     * (coll, root, group )
     */
    hash_table_t collcoa;
    coll_request_t **collcoa_set;
    bsc_size_t collcoa_set_reserved, collcoa_set_size;
    bsc_coll_params_t * collcoa_param_list;
    bsc_size_t collcoa_param_list_size;
} state_t;

static state_t s_bsc;

int collcoa_is_equal( const void * a, const void * b )
{
    const coll_request_t * x = a, * y = b;
    return x->coll == y->coll && x->root == y->root
        && x->group == y->group ;
}

size_t collcoa_hash( const void * a )
{
    const coll_request_t * x = a;
    return ((size_t) x->coll)*2 +
           ((size_t) x->group)*3 +
            5*x->root ; 
}


static void expand_coll( coll_request_t * c )
{
    coll_request_t * item; 
    if ( s_bsc.collcoa.n_buckets == 0 ) {
        hash_table_create( & s_bsc.collcoa, 1000,
                collcoa_is_equal, collcoa_hash );
    }
    
    if ( s_bsc.collcoa_set_size == s_bsc.collcoa_set_reserved )
    {
        bsc_size_t new_collcoa_set_reserved =
            MAX(1000, 2*s_bsc.collcoa_set_reserved) ;
        coll_request_t ** new_collcoa_set = calloc( 
                new_collcoa_set_reserved, sizeof(coll_request_t*) );
        if ( !new_collcoa_set )
            bsp_abort("bsc_sync: insufficient memory\n");
        memcpy( new_collcoa_set, s_bsc.collcoa_set, 
                sizeof(coll_request_t*)*s_bsc.collcoa_set_size );
        free(s_bsc.collcoa_set);
        s_bsc.collcoa_set = new_collcoa_set;
        s_bsc.collcoa_set_reserved = new_collcoa_set_reserved;
    }
               
    item = hash_table_new_item( &s_bsc.collcoa, c );
    if (item != c) {
        c->next = item->next;
        item->next = c;
    }
    else {
        assert( s_bsc.collcoa_set_size < s_bsc.collcoa_set_reserved );
        s_bsc.collcoa_set[ s_bsc.collcoa_set_size++ ] = item;
    }
}

static void schedule_colls(void)
{
    bsc_size_t i, j;
    for ( i = 0; i < s_bsc.collcoa_set_size; ++i ) {
        coll_request_t * ptr = s_bsc.collcoa_set[i];
        const bsc_collective_t * coll = ptr->coll;
        group_t * group = ptr->group;
        bsc_pid_t root = ptr->root;
        bsc_size_t n = 0;
        double min_cost = HUGE_VAL;
        bsc_size_t min_alg = 0;

        /* put all parameters in an array */
        while (ptr) {
            if (s_bsc.collcoa_param_list_size <= n ) {
                free( s_bsc.collcoa_param_list );
                s_bsc.collcoa_param_list_size = MAX(100,
                        s_bsc.collcoa_param_list_size*2 );
                s_bsc.collcoa_param_list = calloc( 
                        s_bsc.collcoa_param_list_size,
                        sizeof(bsc_coll_params_t) );
                if (!s_bsc.collcoa_param_list)
                    bsp_abort("bsc_sync: insufficient memory\n");
                
                /* start copying from scratch again */
                ptr = s_bsc.collcoa_set[i];
                n = 0;
            } else {
                s_bsc.collcoa_param_list[n] = ptr->params;
                ptr = ptr->next;
                n += 1;
            }
        }

        /* compute cheapest algorithm */
        for ( j = 0; j < coll->n_algs; ++j ) {
            double c = (*coll->costfuncs[j])(group, s_bsc.collcoa_param_list, n);
            if (c < min_cost ) {
                min_cost = c;
                min_alg = j;
            }
        }

        /* execute cheapest algorithm */
        (*coll->algorithms[min_alg])( s_bsc.current, root, group,
                s_bsc.collcoa_param_list, n );
    }
}

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

bsc_step_t bsc_collective( bsc_step_t depends,
       const bsc_collective_t * collective,
       bsc_pid_t root, bsc_group_t group, 
       bsc_coll_params_t params ) 
{
    request_t * req;
    coll_request_t r;
    bsc_step_t maxsteps = 0, steps;
    bsc_size_t i ;
    for ( i = 0; i < collective->n_algs; ++i ) 
    {
        steps = (*collective->maxsteps[i])( group );
        if ( maxsteps < steps )
            maxsteps = steps;
    }

    r.coll = collective;
    r.root = root;
    r.group = (group_t *) group;
    r.params = params;
    r.next = NULL;

    req = new_request( depends );
    req->kind = COLL;
    req->payload.coll = r;
    return depends + maxsteps;
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
            /* expand all collectives */
            for ( i = 0; i < s_bsc.n_requests[s_bsc.queue_start]; ++i ) {
                 request_t * r = s_bsc.queue + 
                    s_bsc.queue_start * s_bsc.max_requests_per_step + i ;
                 if ( r->kind == COLL ) {
                     expand_coll( & r->payload.coll );
                 }
            }
            /* and insert them in small pieces into the queue */
            schedule_colls();

            /* execute the delayed puts and gets */
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
            /* execute the delayed local reduce operations */
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


struct group {
    /** the number of members in this group */
    bsp_pid_t size;

    /** Array of size nprocs with the local pid for each global pid */
    bsp_pid_t * lid;

    /** the global PIDs of the group members */
    bsp_pid_t * gid;

    /** the memory offset this process in a remote arary */
    bsp_size_t * rof;
} ; 

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
    bsp_push_reg( incoming, P*sizeof(incoming[0]) );
    bsp_sync();

    /* fill the communication matrix with a 1 for each connection */
    for ( i = 0 ; i < size ; ++i ) {
        bsc_size_t offset = i+1;
        bsp_put( neighbours[i], &offset, incoming, s*sizeof(offset), sizeof(offset) );
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

bsc_step_t bsc_steps_1phase( bsc_group_t group )
{
    (void) group;
    return 1;
}

bsc_step_t bsc_steps_2phase( bsc_group_t group )
{
    (void) group;
    return 2;
}

bsc_step_t bsc_steps_2tree( bsc_group_t group )
{
    group_t * g = group;
    bsc_pid_t P = g?g->size:bsp_nprocs();
    return int_log2(P-1)+1;
}


bsc_step_t bsc_scatter_single( bsc_step_t depends,
                        bsc_pid_t root, bsc_group_t group,
                       const void * src, void *dst, bsc_size_t size )
{
    const group_t * g = group;
    bsc_pid_t i, P = g?g->size:bsp_nprocs() ;
    if ( bsp_pid() == root ) {
        for ( i = 0; i < P; ++i ) {
            bsc_pid_t m = g?g->gid[i]:i;
            const char * src_bytes = src;
            bsc_put( depends, m, src_bytes + i * size, dst, 0, size );
        }
    }
    return depends + 1;
}

bsc_step_t bsc_scatter_multiple( bsc_step_t depends,
                        bsc_pid_t root, bsc_group_t group,
                        bsc_coll_params_t * set, bsc_size_t n )
{
    bsc_size_t i;
    for ( i = 0; i < n ; ++i )
        bsc_scatter_single( depends, root, group, 
                set[i].src, set[i].dst, set[i].size );
    return depends + 1;
}

double bsc_scatter_costs( bsc_group_t group,
                        bsc_coll_params_t * set, bsc_size_t n )
{
    double costs = 0.0;
    bsc_size_t i;
    group_t * g = group;
    bsc_pid_t P = g?g->size:bsp_nprocs() ;
    for ( i = 0; i < n ; ++i )
        costs += set[i].size * P * bsc_g();
    costs += bsc_L();
    return costs;
}

const bsc_collective_t bsc_scatter_algs = {
    { bsc_scatter_multiple },
    { bsc_scatter_costs },
    { bsc_steps_1phase },
    1 
};
const bsc_collective_t * const bsc_coll_scatter = & bsc_scatter_algs;

bsc_step_t bsc_scatter( bsc_step_t depends,
                        bsc_pid_t root, bsc_group_t group,
                       const void * src, void *dst, bsc_size_t size )
{
    bsc_coll_params_t p; memset(&p, 0, sizeof(p));
    p.src = src;
    p.dst = dst;
    p.size = size;
    return bsc_collective( depends, bsc_coll_scatter, root, group, p );
}



bsc_step_t bsc_gather_single( bsc_step_t depends, 
                       bsc_pid_t root, bsc_group_t group,
                       const void * src, void *dst, bsc_size_t size )
{
    const group_t * g = group;
    bsc_pid_t i, P = g?g->size:bsp_nprocs() ;
    if ( bsp_pid() == root ) {
        for ( i = 0; i < P; ++i ) {
            bsc_pid_t m = g?g->gid[i]:i;
            char * dst_bytes = dst;
            bsc_get( depends, m, src, 0, dst_bytes + i * size, size );
        }
    }

    return depends + 1;
}

bsc_step_t bsc_gather_multiple( bsc_step_t depends,
                        bsc_pid_t root, bsc_group_t group,
                        bsc_coll_params_t * set, bsc_size_t n )
{
    bsc_size_t i;
    for ( i = 0; i < n ; ++i )
        bsc_gather_single( depends, root, group, 
                set[i].src, set[i].dst, set[i].size );
    return depends + 1;
}

double bsc_gather_costs( bsc_group_t group,
                        bsc_coll_params_t * set, bsc_size_t n )
{
    double costs = 0.0;
    bsc_size_t i;
    group_t * g = group;
    bsc_pid_t P = g?g->size:bsp_nprocs() ;
    for ( i = 0; i < n ; ++i )
        costs += set[i].size * P * bsc_g();
    costs += bsc_L();
    return costs;
}

const bsc_collective_t bsc_gather_algs = {
    { bsc_gather_multiple },
    { bsc_gather_costs },
    { bsc_steps_1phase },
    1 
};
const bsc_collective_t * const bsc_coll_gather = & bsc_gather_algs;

bsc_step_t bsc_gather( bsc_step_t depends,
                        bsc_pid_t root, bsc_group_t group,
                       const void * src, void *dst, bsc_size_t size )
{
    bsc_coll_params_t p; memset(&p, 0, sizeof(p));
    p.src = src;
    p.dst = dst;
    p.size = size;
    return bsc_collective( depends, bsc_coll_gather, root, group, p );
}



bsc_step_t bsc_allgather_single( bsc_step_t depends, bsc_group_t group,
                        const void * src, void * dst, bsc_size_t size )
{
    const group_t * g = group;
    bsc_pid_t i, P = g?g->size:bsp_nprocs() ;
    for ( i = 0; i < P; ++i ) {
        char * dst_bytes = dst;
        bsc_get( depends, g?g->gid[i]:i, src, 0, dst_bytes + i * size, size );
    }
    return depends + 1;
}

bsc_step_t bsc_allgather_multiple( bsc_step_t depends, 
                        bsc_pid_t root, bsc_group_t group,
                        bsc_coll_params_t * set, bsc_size_t n )
{
    bsc_size_t i;
    (void) root;
    for ( i = 0; i < n ; ++i )
        bsc_allgather_single( depends, group, 
                set[i].src, set[i].dst, set[i].size );
    return depends + 1;
}

double bsc_allgather_costs( bsc_group_t group,
                        bsc_coll_params_t * set, bsc_size_t n )
{
    double costs = 0.0;
    bsc_size_t i;
    group_t * g = group;
    bsc_pid_t P = g?g->size:bsp_nprocs() ;
    for ( i = 0; i < n ; ++i )
        costs += set[i].size * P * bsc_g();
    costs += bsc_L();
    return costs;
}

const bsc_collective_t bsc_allgather_algs = {
    { bsc_allgather_multiple },
    { bsc_allgather_costs },
    { bsc_steps_1phase },
    1 
};
const bsc_collective_t * const bsc_coll_allgather = & bsc_allgather_algs;

bsc_step_t bsc_allgather( bsc_step_t depends, bsc_group_t group,
                       const void * src, void *dst, bsc_size_t size )
{
    bsc_coll_params_t p; memset(&p, 0, sizeof(p));
    p.src = src;
    p.dst = dst;
    p.size = size;
    return bsc_collective( depends, bsc_coll_allgather, 0, group, p );
}



bsc_step_t bsc_alltoall_single( bsc_step_t depends, bsc_group_t group,
                      const void * src, void * dst, bsc_size_t size )
{
    const group_t * g = group;
    bsp_pid_t i, P = g?g->size:bsp_nprocs() ;
    const char * src_bytes = (const char *) src;
    for ( i = 0; i < P; ++i ) {
        bsc_pid_t m = g?g->gid[i]:i;
        bsc_size_t doff = g?g->rof[i]:bsp_pid();
        bsc_put( depends, m, src_bytes + i * size, dst, doff * size, size );
    }
    return depends + 1;
}

bsc_step_t bsc_alltoall_multiple( bsc_step_t depends,
                        bsc_pid_t root, bsc_group_t group,
                        bsc_coll_params_t * set, bsc_size_t n )
{
    bsc_size_t i;
    (void) root;
    for ( i = 0; i < n ; ++i )
        bsc_alltoall_single( depends, group,
                set[i].src, set[i].dst, set[i].size );
    return depends + 1;
}

double bsc_alltoall_costs( bsc_group_t group,
                        bsc_coll_params_t * set, bsc_size_t n )
{
    double costs = 0.0;
    bsc_size_t i;
    group_t * g = group;
    bsc_pid_t P = g?g->size:bsp_nprocs() ;
    for ( i = 0; i < n ; ++i )
        costs += set[i].size * P * bsc_g();
    costs += bsc_L();
    return costs;
}

const bsc_collective_t bsc_alltoall_algs = {
    { bsc_alltoall_multiple },
    { bsc_alltoall_costs },
    { bsc_steps_1phase },
    1 
};
const bsc_collective_t * const bsc_coll_alltoall = & bsc_alltoall_algs;

bsc_step_t bsc_alltoall( bsc_step_t depends, bsc_group_t group,
                       const void * src, void *dst, bsc_size_t size )
{
    bsc_coll_params_t p; memset(&p, 0, sizeof(p));
    p.src = src;
    p.dst = dst;
    p.size = size;
    return bsc_collective( depends, bsc_coll_alltoall, 0, group, p );
}


bsc_step_t bsc_bcast_qtree_single( bsc_step_t depends,
                      bsp_pid_t root, bsc_group_t group, 
                      const void * src, void * dst, bsc_size_t size,
                      bsc_pid_t q )
{
    const group_t * g = group;
    bsc_pid_t i, j, P = g?g->size:bsp_nprocs() ;
    bsc_pid_t range, r = g?g->lid[root]:root;
    bsc_pid_t me = ((g?g->lid[ bsp_pid() ] : bsp_pid()) + P - r) %P;

    if ( me == 0 ) 
        memcpy( dst, src, size ); 

    for (range = 1; range < P; range *= q ) {
        if ( me < range ) {
            for ( i = 1 ; i < q; ++i) {
                j = i * range + me; 
                if (j < P) {
                    j = ( j + r ) % P;
                    bsc_put( depends, g?g->gid[j]:j, dst, dst, 0, size );
                }
            }
        }
        depends += 1;
    }
    return depends;
}

bsc_step_t bsc_bcast_1phase_multiple( bsc_step_t depends,
                      bsp_pid_t root, bsc_group_t group, 
                      bsc_coll_params_t * params, bsc_size_t n)
{
    group_t * g = group;
    bsc_pid_t P = g?g->size:bsp_nprocs() ;
    bsc_size_t i;
    bsc_step_t finished = depends;
    for ( i = 0; i < n ; ++i ) {
        finished = bsc_bcast_qtree_single( depends, root, group,
                params[i].src, params[i].dst, params[i].size, P );
    }

    return finished;
}

bsc_step_t bsc_bcast_2tree_multiple( bsc_step_t depends,
                      bsp_pid_t root, bsc_group_t group, 
                      bsc_coll_params_t * params, bsc_size_t n)
{
    bsc_size_t i;
    bsc_step_t finished = depends;
    for ( i = 0; i < n ; ++i ) {
        finished = bsc_bcast_qtree_single( depends, root, group,
                params[i].src, params[i].dst, params[i].size, 2 );
    }

    return finished;
}

bsc_step_t bsc_bcast_root_p_tree_multiple( bsc_step_t depends,
                      bsp_pid_t root, bsc_group_t group, 
                      bsc_coll_params_t * params, bsc_size_t n)
{
    group_t * g = group;
    bsc_pid_t P = g?g->size:bsp_nprocs(), root_n = ceil(sqrt(P)) ;
    bsc_size_t i;
    bsc_step_t finished = depends;
    for ( i = 0; i < n ; ++i ) {
        finished = bsc_bcast_qtree_single( depends, root, group,
                params[i].src, params[i].dst, params[i].size, root_n );
    }

    return finished;
}



bsc_step_t bsc_bcast_2phase( bsc_step_t depends,
                      bsp_pid_t root, bsc_group_t group, 
                      bsc_coll_params_t * params, bsc_size_t n )
{
    const group_t * g = group;
    bsc_pid_t P = g?g->size:bsp_nprocs() ;
    bsc_pid_t s = g?g->lid[bsp_pid()]:bsp_pid();
    bsc_size_t i, total_size = 0;
    bsc_size_t block;
    bsc_size_t b, my_start_b, my_block_size;
    bsc_size_t j, my_start_j;
    const char * src_bytes;
    char * dst_bytes;
 
    /* compute preferred block size */
    for ( i = 0 ; i < n; ++i )
        total_size += params[i].size;
    block = (total_size + P - 1)/P;
   
    /* the root scatters, all compute the block for the
     * following all-gather */
    b = 0; j = 0;
    for ( i = 0 ; i < P; ++i) {
        bsc_size_t block_begin = MIN(     i*block, total_size);
        bsc_size_t block_end   = MIN( (i+1)*block, total_size);
        bsc_size_t block_remaining = block_end - block_begin;

        if (i == s ) {
            my_start_j = j;
            my_start_b = b;
            my_block_size = block_remaining;
        }
        
        do {
            bsc_size_t e = b + block_remaining;
            bsc_size_t next_b = e;
            bsc_size_t next_j = j;
            if (e > params[j].size) {
                e = params[j].size;
                next_j = j+1;
                next_b = 0;
            }
            src_bytes = params[j].src;
            dst_bytes = params[j].dst;
            if (bsp_pid() == root ) 
                bsc_put( depends, g?g->gid[i]:i, 
                        src_bytes + b, dst_bytes, b, e - b);
            block_remaining -= (e-b);
            j = next_j;
            b = next_b;
        } while ( block_remaining > 0 ) ;

    }
    depends += 1;
    /* the all-gather */
    for ( i = 0; i < P; ++i ) {
        bsc_size_t block_remaining = my_block_size;
        b = my_start_b;
        j = my_start_j;

        do {
            bsc_size_t e = b + block_remaining;
            bsc_size_t next_b = e;
            bsc_size_t next_j = j;
            if (e > params[j].size) {
                e = params[j].size;
                next_j = j+1;
                next_b = 0;
            }
            src_bytes = params[j].src;
            dst_bytes = params[j].dst;
            bsc_put( depends, g?g->gid[i]:i, 
                    dst_bytes + b, dst_bytes, b, e - b);
            block_remaining -= (e-b);
            j = next_j;
            b = next_b;
        } while ( block_remaining > 0 ) ;
    }
    return depends+1;
}

double bsc_bcast_1phase_costs( bsc_group_t group,
         bsc_coll_params_t * set, bsc_size_t n )
{
    bsc_size_t i, total_size = 0;
    group_t * g = group;
    bsc_pid_t P = g?g->size:bsp_nprocs() ;
    for ( i = 0; i < n ; ++i )
        total_size += set[i].size;

    return bsc_g() * P * total_size + bsc_L();
}


double bsc_bcast_2tree_costs(  bsc_group_t group,
         bsc_coll_params_t * set, bsc_size_t n )
{
    bsc_size_t i, total_size = 0;
    group_t * g = group;
    bsc_pid_t P = g?g->size:bsp_nprocs() ;
    for ( i = 0; i < n ; ++i )
        total_size += set[i].size;

    return (int_log2(P-1)+1) * (total_size * bsc_g() +  bsc_L() ) ;
}

double bsc_bcast_root_p_tree_costs(  bsc_group_t group,
         bsc_coll_params_t * set, bsc_size_t n )
{
    bsc_size_t i, total_size = 0;
    group_t * g = group;
    bsc_pid_t P = g?g->size:bsp_nprocs() ;
    bsc_pid_t root_P = ceil(sqrt(P));
    for ( i = 0; i < n ; ++i )
        total_size += set[i].size;

    return 2*(root_P - 1) * total_size * bsc_g() + 2 * bsc_L();
}

double bsc_bcast_2phase_costs(  bsc_group_t group,
         bsc_coll_params_t * set, bsc_size_t n )
{
    bsc_size_t i, total_size = 0;
    group_t * g = group;
    bsc_pid_t P = g?g->size:bsp_nprocs() ;
    for ( i = 0; i < n ; ++i )
        total_size += set[i].size;

    return (total_size + ((total_size+P-1)/P)*(P-2)) * bsc_g() + 2 * bsc_L();
}


const bsc_collective_t bsc_bcast_algs = {
    { bsc_bcast_1phase_multiple, 
      bsc_bcast_2tree_multiple,
      bsc_bcast_root_p_tree_multiple,
      bsc_bcast_2phase },
    { bsc_bcast_1phase_costs,
      bsc_bcast_2tree_costs,
      bsc_bcast_root_p_tree_costs,
      bsc_bcast_2phase_costs },
    { bsc_steps_1phase,
      bsc_steps_2tree ,
      bsc_steps_2phase,
      bsc_steps_2phase },
    4 
};
const bsc_collective_t * const bsc_coll_bcast = & bsc_bcast_algs;

bsc_step_t bsc_bcast( bsc_step_t depends, bsc_pid_t root, bsc_group_t group, 
                       const void * src, void *dst, bsc_size_t size )
{
    bsc_coll_params_t p; memset(&p, 0, sizeof(p));
    p.src = src;
    p.dst = dst;
    p.size = size;
    return bsc_collective( depends, bsc_coll_bcast, root, group, p );
}



bsc_step_t bsc_reduce_qtree_single( bsc_step_t depends, 
        bsc_pid_t root, bsc_group_t group,
        const void * src, void * dst, void * tmp_space,
        bsc_reduce_t reducer, const void * zero,
        bsc_size_t nmemb, bsc_size_t size, bsc_pid_t q )
{
    const group_t * g = group;
    bsc_pid_t i, j, k, P = g?g->size:bsp_nprocs() ;
    bsc_pid_t range, r = g?g->lid[root]:root;
    bsc_pid_t me = ((g?g->lid[ bsp_pid() ] : bsp_pid()) + P - r) %P;
    char * tmp_bytes = tmp_space;

    (*reducer)( dst, zero, src, size * nmemb );
    if (P < 2)
        return depends;
    
    for (range = int_pow( q, int_log( q, P-1));
            range >= 1; range /= q ) {
        if ( me < range ) {
            k = 1;
            for ( i = 1 ; i < q; ++i) {
                j = i * range + me; 
                if (j < P) {
                    j = (j + r )% P;
                    bsc_get( depends, g?g->gid[j]:j, dst, 0, 
                            tmp_bytes + k*size, size );
                    k += 1;
                }
            }
            bsc_get( depends, bsp_pid(), dst, 0, tmp_space, size );
            bsc_exec_reduce( depends, reducer, dst, zero, 
                    tmp_space, size * k );
        }
        depends += 1;
    }
    return depends;

}

bsc_step_t bsc_reduce_1phase_multiple( bsc_step_t depends,
                      bsp_pid_t root, bsc_group_t group, 
                      bsc_coll_params_t * params, bsc_size_t n)
{
    group_t * g = group;
    bsc_pid_t P = g?g->size:bsp_nprocs() ;
    bsc_size_t i;
    bsc_step_t finished = depends;
    for ( i = 0; i < n ; ++i ) {
        finished = bsc_reduce_qtree_single( depends, root, group,
                   params[i].src, params[i].dst, params[i].tmp,
                   params[i].reducer, params[i].zero, 
                   params[i].nmemb, params[i].size, P );
    }

    return finished;
}

bsc_step_t bsc_reduce_2tree_multiple( bsc_step_t depends,
                      bsp_pid_t root, bsc_group_t group, 
                      bsc_coll_params_t * params, bsc_size_t n)
{
    bsc_size_t i;
    bsc_step_t finished = depends;
    for ( i = 0; i < n ; ++i ) {
        finished = bsc_reduce_qtree_single( depends, root, group,
                   params[i].src, params[i].dst, params[i].tmp,
                   params[i].reducer, params[i].zero, 
                   params[i].nmemb, params[i].size, 2 );
    }

    return finished;
}

bsc_step_t bsc_reduce_root_p_tree_multiple( bsc_step_t depends,
                      bsp_pid_t root, bsc_group_t group, 
                      bsc_coll_params_t * params, bsc_size_t n)
{
    group_t * g = group;
    bsc_pid_t P = g?g->size:bsp_nprocs(), root_P = ceil(sqrt(P)) ;
    bsc_size_t i;
    bsc_step_t finished = depends;
    for ( i = 0; i < n ; ++i ) {
        finished = bsc_reduce_qtree_single( depends, root, group,
                   params[i].src, params[i].dst, params[i].tmp,
                   params[i].reducer, params[i].zero, 
                   params[i].nmemb, params[i].size, root_P );
    }

    return finished;
}

double bsc_reduce_1phase_costs( bsc_group_t group,
         bsc_coll_params_t * set, bsc_size_t n )
{
    bsc_size_t i, total_size = 0;
    group_t * g = group;
    bsc_pid_t P = g?g->size:bsp_nprocs() ;
    for ( i = 0; i < n ; ++i )
        total_size += set[i].size;

    return bsc_g() * P * total_size + bsc_L();
}

double bsc_reduce_2tree_costs(  bsc_group_t group,
         bsc_coll_params_t * set, bsc_size_t n )
{
    bsc_size_t i, total_size = 0;
    group_t * g = group;
    bsc_pid_t P = g?g->size:bsp_nprocs() ;
    for ( i = 0; i < n ; ++i )
        total_size += set[i].size;

    return (int_log2(P-1)+1) * (total_size * bsc_g() +  bsc_L() ) ;
}

double bsc_reduce_root_p_tree_costs(  bsc_group_t group,
         bsc_coll_params_t * set, bsc_size_t n )
{
    bsc_size_t i, total_size = 0;
    group_t * g = group;
    bsc_pid_t P = g?g->size:bsp_nprocs() ;
    bsc_pid_t root_P = ceil(sqrt(P));
    for ( i = 0; i < n ; ++i )
        total_size += set[i].size;

    return 2*(root_P - 1) * total_size * bsc_g() + 2 * bsc_L();
}


const bsc_collective_t bsc_reduce_algs = {
    { bsc_reduce_1phase_multiple, 
      bsc_reduce_2tree_multiple,
      bsc_reduce_root_p_tree_multiple,
    },      
    { bsc_reduce_1phase_costs,
      bsc_reduce_2tree_costs,
      bsc_reduce_root_p_tree_costs,
    },
    { bsc_steps_1phase,
      bsc_steps_2tree ,
      bsc_steps_2phase,
    },
    3 
};
const bsc_collective_t * const bsc_coll_reduce = & bsc_reduce_algs;

bsc_step_t bsc_reduce( bsc_step_t depends, 
        bsc_pid_t root, bsc_group_t group,
        const void * src, void * dst, void * tmp_space,
        bsc_reduce_t reducer, const void * zero,
        bsc_size_t nmemb, bsc_size_t size )
{
    bsc_coll_params_t p; memset(&p, 0, sizeof(p));
    p.src = src;
    p.dst = dst;
    p.tmp = tmp_space;
    p.reducer = reducer;
    p.zero = zero;
    p.nmemb = nmemb;
    p.size = size;
    return bsc_collective( depends, bsc_coll_reduce, root, group, p );
}


bsc_step_t bsc_allreduce_1phase( bsc_step_t depends, bsc_group_t group,
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

bsc_step_t bsc_allreduce_qtree( bsc_step_t depends, bsc_group_t group,
        const void * src, void * dst, void * tmp_space,
        bsc_reduce_t reducer, const void * zero,
        bsc_size_t nmemb, bsc_size_t size, 
        bsc_pid_t q )
{
    const group_t * g = group;
    bsc_pid_t i, j, k, P = g?g->size:bsp_nprocs() ;
    bsc_pid_t range, s = g?g->lid[bsp_pid()]:bsp_pid();
    char * tmp_bytes = tmp_space;

    (*reducer)( dst, zero, src, size * nmemb );

    for (range = 1; range < P; range *= q ) {
        k = 1;
        for ( i = 1 ; i < q; ++i) {
            if (i * range < P) {
                j = i * range + s; 
                bsc_get( depends, g?g->gid[j]:j, dst, 0,
                        tmp_bytes + k * size, size);
                k+=1;
            }
        }
        bsc_get( depends, bsp_pid(), dst, 0, tmp_bytes, size );
        bsc_exec_reduce( depends, reducer, dst, zero, tmp_space, size*k );
        depends += 1;
    }
    return depends;
}



bsc_step_t bsc_allreduce( bsc_step_t depends, bsc_group_t group,
        const void * src, void * dst, void * tmp_space,
        bsc_reduce_t reducer, const void * zero,
        bsc_size_t nmemb, bsc_size_t size )
{
    const group_t * g = group;
    int i;
    bsc_pid_t P = g?g->size:bsp_nprocs();
    bsc_pid_t root_n = ceil( sqrt( (double) P ) );
    enum ALG { ONE_PHASE, TWO_TREE, ROOT_TREE, N_ALGS } ;
    double costs[N_ALGS];
    enum ALG min_alg = ONE_PHASE ;
    double min_cost = HUGE_VAL;
    costs[ ONE_PHASE] = P * size * bsc_g() + bsc_L();
    costs[ TWO_TREE]  = (int_log2(P-1)+1) * (size * bsc_g() +  bsc_L());
    costs[ ROOT_TREE] = 2*(root_n - 1) * size * bsc_g() + 2 * bsc_L();

    for ( i = 0; i < N_ALGS; ++i )
        if ( costs[i] < min_cost ) { 
            min_cost = costs[i];
            min_alg = (enum ALG) i;
        }

    switch (min_alg) {
        case ONE_PHASE:
            return bsc_allreduce_1phase( depends, group, src, dst, tmp_space,
                      reducer, zero, nmemb, size );

        case TWO_TREE: 
            return bsc_allreduce_qtree( depends, group, src, dst, tmp_space,
                      reducer, zero, nmemb, size, 2 );

        case ROOT_TREE: 
            return bsc_allreduce_qtree( depends, group, src, dst, tmp_space,
                      reducer, zero, nmemb, size, root_n );
 
        case N_ALGS:
            bsp_abort("bsc_allreduce: missing algorithm\n");
            return depends;
    }
    return depends;
}

bsc_step_t bsc_scan_qtree( bsc_step_t depends, bsc_group_t group,
        const void * src, void * dst, void * tmp_space, 
        bsc_reduce_t reducer, const void * zero,
        bsc_size_t nmemb, bsc_size_t size,
        int q )
{
    const group_t * g = group;
    bsc_pid_t i, j, k, P = g?g->size:bsp_nprocs() ;
    bsc_pid_t range, s = g?g->lid[bsp_pid()]:bsp_pid();
    char * tmp_bytes = tmp_space;
    char * dst_bytes = dst;
    const void * prefix;

    (*reducer)( dst, zero, src, size * nmemb );

    for (range = 1; range < P; range *= q ) {
        const bsc_size_t shift = (range==1);
        const bsc_size_t last = range == 1 ? (nmemb-1)*size : q*size;
        for ( i = 1; i < q; ++i ) {
            j = s + range * i + shift;
            if ( j < P ) {
                bsc_put( depends, g?g->gid[j]:j, dst_bytes + last,
                        tmp_bytes , (q-i-1) * size, size);
            }
        }

        k = (s >= shift);
        for ( i = 1; i < q; ++i )
            k += (s >= range*i + shift);

        if (s < P - shift) {
            j = g?g->gid[s+shift]:s+shift;
            bsc_put(depends, j, dst_bytes + last, tmp_bytes, last, size);
        }
        bsc_exec_reduce( depends, reducer, dst, zero, tmp_space, size*k );
        depends += 1;
    }

    prefix = k > 0 ? tmp_bytes + (k-1)*size : zero ;
    bsc_exec_reduce( depends, reducer, dst, prefix, src, size*nmemb );
    return depends;
}



bsc_step_t bsc_scan( bsc_step_t depends, bsc_group_t group,
        const void * src, void * dst, void * tmp_space,
        bsc_reduce_t reducer, const void * zero,
        bsc_size_t nmemb, bsc_size_t size )
{
    const group_t * g = group;
    int i;
    bsc_pid_t P = g?g->size:bsp_nprocs();
    bsc_pid_t root_P = ceil( sqrt( (double) P ) );
    enum ALG { ONE_PHASE, TWO_TREE, ROOT_TREE, N_ALGS } ;
    double costs[N_ALGS];
    enum ALG min_alg = ONE_PHASE ;
    double min_cost = HUGE_VAL;

    costs[ ONE_PHASE] = (P-1) * size * bsc_g() + bsc_L();
    costs[ TWO_TREE]  = int_log2(P-1) * (size * bsc_g() +  bsc_L())
                         + size * bsc_g();
    costs[ ROOT_TREE] = (2*root_P - 1) * size * bsc_g() + 2 * bsc_L();

    for ( i = 0; i < N_ALGS; ++i )
        if ( costs[i] < min_cost ) { 
            min_cost = costs[i];
            min_alg = (enum ALG) i;
        }

    switch (min_alg) {
        case ONE_PHASE:
            return bsc_scan_qtree( depends, group, src, dst, tmp_space,
                      reducer, zero, nmemb, size, P );

        case TWO_TREE: 
            return bsc_scan_qtree( depends, group, src, dst, tmp_space,
                      reducer, zero, nmemb, size, 2 );

        case ROOT_TREE: 
            return bsc_scan_qtree( depends, group, src, dst, tmp_space,
                      reducer, zero, nmemb, size, root_P );

        case N_ALGS:
            bsp_abort("bsc_scan: missing algorithm\n");
            return depends;
    }
    return depends;
}




