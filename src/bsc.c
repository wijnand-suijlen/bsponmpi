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

typedef struct coll_header {
    const bsc_collective_t * coll ;
    group_t * group ;
    bsc_pid_t root ;
    bsc_size_t n_params;
} coll_header_t;


typedef struct state { 
    bsc_step_t current;
    unsigned   horizon;
    unsigned   max_requests_per_step;
    unsigned   queue_start;
    unsigned   * n_requests;
    request_t  * queue;
    unsigned   total_n_outstanding_requests;
    bsc_step_t * next_step;

    /* Collective coalescing table. The items are the
     * coll_request_t structs. Keys are made up of
     * (coll, root, group )
     */
    hash_table_t collcoa;
    coll_request_t **collcoa_set;
    bsc_size_t collcoa_set_reserved, collcoa_set_size;
    bsc_coll_params_t * collcoa_param_list;
    coll_header_t *collcoa_headers;
    bsc_size_t collcoa_param_list_size;
} state_t;


static double load_param(const char * param, double default_val, 
        const char * func ) 
{
    char * numend = NULL;
    char * str = getenv(param);
    double result = str?strtod( str, &numend ):default_val;
    if ( (str && numend == str) || result == HUGE_VAL 
            || result == -HUGE_VAL ) {
        result = default_val;
        fprintf(stderr, "%s: %s was not a finite positve "
                        "floating point number."
                        " Using default: %g\n", 
                        func, param, result );
    }
    return result;
}

static double s_g[BSC_N_REQUEST_TYPES];
static double s_o[BSC_N_REQUEST_TYPES];
static double s_L[BSC_N_REQUEST_TYPES];

void load_param_g_once(void)
{
    const char * env_names_g[] = { 
        "BSC_PUT_G", "BSP_HPPUT_G", "BSP_GET_G",
        "BSP_HPGET_G", "BSP_SEND_G" } ;
    const char * env_names_o[] = { 
        "BSC_PUT_O", "BSP_HPPUT_O", "BSP_GET_O",
        "BSP_HPGET_O", "BSP_SEND_O" } ;
    int i;
    for ( i = 0; i < BSC_N_REQUEST_TYPES; ++i )
        if (env_names_g[i]) {
                s_g[i] = load_param(env_names_g[i], 1.0, "bsc_g");
                s_o[i] = load_param(env_names_o[i], 0.0, "bsc_g");
        }
}

void load_param_L_once(void) {
    const char * env_names[] = { 
        "BSC_PUT_L", "BSP_HPPUT_L", "BSP_GET_L",
        "BSP_HPGET_L", "BSP_SEND_L" } ;
    int i;
    for ( i = 0; i < BSC_N_REQUEST_TYPES; ++i )
        if (env_names[i])
                s_L[i] = load_param(env_names[i], 1.0e+4, "bsc_L");
}

#ifdef PTHREADSAFE
#include <pthread.h>
static pthread_key_t s_tls;
static pthread_once_t s_tls_init = PTHREAD_ONCE_INIT;

static state_t * state_create( void ) {
    return calloc( 1, sizeof(state_t) );
}

static void state_destroy( void * s ) {
    state_t * state = s;
    free(state->n_requests );
    free(state->queue );
    free(state->next_step );

    hash_table_destroy( &state->collcoa );
    free(state->collcoa_set);
    free(state->collcoa_headers);
    free(state->collcoa_param_list );
    memset( state, 0, sizeof(state_t) );
    free(state);
}

static void init_tls(void) {
    if ( pthread_key_create( &s_tls, state_destroy ) ) {
        fprintf(stderr, "FATAL ERROR: Cannot initalize thread local storage\n");
        exit(EXIT_FAILURE);
    }
}

static state_t * bsc()  {
    void * data = NULL;
    (void) pthread_once( &s_tls_init, init_tls );
    data = pthread_getspecific( s_tls );
    if ( NULL == data ) {
        data = state_create();
        if (!data || pthread_setspecific( s_tls, data )) {
            fprintf(stderr, "FATAL ERROR: Cannot initialize thread local storage\n");
            exit(EXIT_FAILURE);
        }
    }
    return data;
}


double bsc_g(bsc_request_t type, int word_size)
{
    static pthread_once_t once = PTHREAD_ONCE_INIT;
    (void) pthread_once( &once, load_param_g_once );
    return s_g[type] + s_o[type] / word_size;
}

double bsc_L(bsc_request_t type)
{
    static pthread_once_t once = PTHREAD_ONCE_INIT;
    (void) pthread_once( &once, load_param_L_once );
    return s_L[type];
}

#else

static state_t * bsc(void) {
    static state_t s_bsc;
    return & s_bsc; 
}

double bsc_g(bsc_request_t type, int word_size)
{
    if ( s_g[0] == 0.0 ) {
        load_param_g_once();
    }
    return s_g[type] + s_o[type] / word_size;
}

double bsc_L(bsc_request_t type)
{
    if ( s_L[0] == 0.0 ) {
        load_param_L_once();
    }
    return s_L[type];
}
#endif

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
    state_t * const sbsc = bsc();
    if ( sbsc->collcoa.n_buckets == 0 ) {
        hash_table_create( & sbsc->collcoa, 1000,
                collcoa_is_equal, collcoa_hash );
    }
    
    if ( sbsc->collcoa_set_size == sbsc->collcoa_set_reserved )
    {
        bsc_size_t new_collcoa_set_reserved =
            MAX(1000, 2*sbsc->collcoa_set_reserved) ;
        coll_request_t ** new_collcoa_set = calloc( 
                new_collcoa_set_reserved, sizeof(coll_request_t*) );
        coll_header_t * new_headers = calloc( 
                new_collcoa_set_reserved, sizeof(coll_header_t) );
        if ( !new_collcoa_set || ! new_headers )
            bsp_abort("bsc_sync: insufficient memory\n");
        memcpy( new_collcoa_set, sbsc->collcoa_set, 
                sizeof(coll_request_t*)*sbsc->collcoa_set_size );
        memcpy( new_headers, sbsc->collcoa_headers,
                sizeof(coll_header_t)*sbsc->collcoa_set_size );
        free(sbsc->collcoa_set);
        free(sbsc->collcoa_headers);
        sbsc->collcoa_set = new_collcoa_set;
        sbsc->collcoa_headers = new_headers;
        sbsc->collcoa_set_reserved = new_collcoa_set_reserved;
    }
               
    item = hash_table_new_item( &sbsc->collcoa, c );
    if (item != c) {
        c->next = item->next;
        item->next = c;
    }
    else {
        assert( sbsc->collcoa_set_size < sbsc->collcoa_set_reserved );
        sbsc->collcoa_set[ sbsc->collcoa_set_size++ ] = item;
    }
}

static void schedule_colls(void)
{
    bsc_size_t i, j, k;
    state_t * const sbsc = bsc();

    /* first copy all parameter lists from the request queue 
     * because replacing collectives with put, get requests
     * will probably cause reallocation so that our 
     * hash table structure gets invalidated */
    for ( i = 0, k = 0; i < sbsc->collcoa_set_size; ++i ) {
        coll_request_t * ptr = sbsc->collcoa_set[i];
        bsc_size_t n = 0;
        
        sbsc->collcoa_headers[i].coll = ptr->coll;
        sbsc->collcoa_headers[i].group = ptr->group;
        sbsc->collcoa_headers[i].root = ptr->root;

        /* put all parameters in an array */
        while (ptr) {
            if (sbsc->collcoa_param_list_size <= k+n ) {
                bsc_size_t new_size = MAX(100,
                        sbsc->collcoa_param_list_size*2 );
                bsc_coll_params_t * new_params = calloc( 
                        new_size, sizeof(bsc_coll_params_t) );

                if (!new_params)
                    bsp_abort("bsc_sync: insufficient memory\n");
                
                memcpy( new_params, sbsc->collcoa_param_list, 
                  sbsc->collcoa_param_list_size * sizeof(bsc_coll_params_t) ); 

                free( sbsc->collcoa_param_list );
                sbsc->collcoa_param_list_size = new_size;
                sbsc->collcoa_param_list = new_params;
            } else {
                sbsc->collcoa_param_list[k+n] = ptr->params;
                ptr = ptr->next;
                n += 1;
            }
        }
        sbsc->collcoa_headers[i].n_params = n;
        k += n;
    }

    for ( i = 0, k = 0 ; i < sbsc->collcoa_set_size; ++i ) {
        coll_header_t  hdr = sbsc->collcoa_headers[i];
        const bsc_collective_t * coll = hdr.coll;
        group_t * group = hdr.group ;
        bsc_pid_t root = hdr.root  ;
        bsc_size_t n = hdr.n_params;

        double min_cost = HUGE_VAL;
        bsc_size_t min_alg = 0;
        /* compute cheapest algorithm */
        for ( j = 0; j < coll->n_algs; ++j ) {
            double c = (*coll->costfuncs[j])(group, sbsc->collcoa_param_list + k, n);
            if (c < min_cost ) {
                min_cost = c;
                min_alg = j;
            }
        }

        /* execute cheapest algorithm */
        (*coll->algorithms[min_alg])( sbsc->current, root, group,
                sbsc->collcoa_param_list + k, n );

        k += n;
    }

    /* clear hash table for next use */
    hash_table_clear( &sbsc->collcoa );
    sbsc->collcoa_set_size = 0;
}

static request_t * new_request( bsc_step_t step ) 
{
    unsigned i, j, delay;
    request_t * entry;
    state_t * const sbsc = bsc();

    if ( step < sbsc->current )
        bsp_abort("bsc: a flush sync is required before posting"
                  " another request. Requested superstep %d comes before current step %d\n",
                  step, sbsc->current );

    delay = step - sbsc->current;
    if ( delay >= sbsc->horizon || 
          sbsc->n_requests[ (sbsc->queue_start + delay) % sbsc->horizon ] 
            >= sbsc->max_requests_per_step ) {

        /* new memory allocation is necessary */
        bsc_step_t new_horizon;
        unsigned * new_n_requests;
        unsigned new_max_requests;
        request_t * new_queue;

        new_horizon = sbsc->horizon;
        if (delay >= sbsc->horizon) {
            new_horizon = MAX( 2*sbsc->horizon, delay + 1);
        }

        if ( new_horizon <= delay )
            bsp_abort("bsc: cannot plan that far in the future"
                      "due to integer overflow\n");
        
        new_n_requests = calloc( new_horizon, sizeof(unsigned));
        if ( ! new_n_requests ) 
            bsp_abort("bsc: insufficient memory; "
                      "cannot allocate %u x %u bytes\n",
                new_horizon, (unsigned) sizeof(unsigned) );

        for ( i = 0; i < sbsc->horizon; ++i ) {
            j = (sbsc->queue_start + i) % sbsc->horizon;
            new_n_requests[ i ] = sbsc->n_requests[j];
        }

        new_max_requests = sbsc->max_requests_per_step; 
        if (new_max_requests <= new_n_requests[delay]) 
            new_max_requests = MAX( 2 * sbsc->max_requests_per_step, 
                    new_n_requests[delay] + 1  );

        if ( new_max_requests <= new_n_requests[delay] )
            bsp_abort("bsc: cannot plan that far in the future"
                      "due to integer overflow\n");
     
        new_queue = calloc( (size_t) new_horizon * new_max_requests, 
                            sizeof(request_t) );

        for ( i = 0; i < sbsc->horizon; ++i ) {
            for ( j = 0; j < sbsc->max_requests_per_step; ++j ){
                unsigned a = i * new_max_requests + j;
                unsigned b = (i + sbsc->queue_start) % sbsc->horizon 
                              * sbsc->max_requests_per_step + j;
                new_queue[ a ] = sbsc->queue[ b ];
            }
        }

        free( sbsc->n_requests );
        free( sbsc->queue );
        sbsc->n_requests = new_n_requests;
        sbsc->queue = new_queue;
        sbsc->horizon = new_horizon;
        sbsc->max_requests_per_step = new_max_requests;
        sbsc->queue_start = 0;
    }

    i = (delay + sbsc->queue_start) % sbsc->horizon;
    entry = sbsc->queue + i * 
        sbsc->max_requests_per_step + sbsc->n_requests[i] ;
    
    sbsc->n_requests[i] += 1;
    sbsc->total_n_outstanding_requests += 1;

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
    return bsc()->current;
}

bsc_step_t bsc_sync( bsc_step_t until )
{
    state_t * const sbsc = bsc();
    unsigned i, P = bsp_nprocs(), sz=sizeof(sbsc->next_step[0]);
    if ( !sbsc->next_step) {
        /* do some initialization */
        sbsc->next_step = calloc( bsp_nprocs(), sizeof(sbsc->next_step[0]));
        bsp_push_reg( sbsc->next_step, 
                bsp_nprocs() * sizeof(sbsc->next_step[0]) );
        bsp_sync();
    }

    while ( sbsc->current < until || until == bsc_flush ) {
        bsc_step_t common_step = until;
        for ( i = 0; i < P; ++i ) 
            sbsc->next_step[i] = bsc_flush;

        if (sbsc->horizon > 0) {
            /* expand all collectives */
            for ( i = 0; i < sbsc->n_requests[sbsc->queue_start]; ++i ) {
                 request_t * r = sbsc->queue + 
                    sbsc->queue_start * sbsc->max_requests_per_step + i ;
                 if ( r->kind == COLL ) {
                     expand_coll( & r->payload.coll );
                 }
            }
            /* and insert them in small pieces into the queue */
            schedule_colls();

            /* execute the delayed puts and gets */
            for ( i = 0; i < sbsc->n_requests[sbsc->queue_start]; ++i ) {
                 request_t * r = sbsc->queue + 
                    sbsc->queue_start * sbsc->max_requests_per_step + i ;
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
            sbsc->total_n_outstanding_requests 
                -= sbsc->n_requests[ sbsc->queue_start ];

            /* search for the step this process can skip to */
            if ( 0 == sbsc->total_n_outstanding_requests ) {
                common_step = until;
            } else {
                for (common_step = 1; common_step < sbsc->horizon ; common_step++){
                    unsigned j = (sbsc->queue_start + common_step) % sbsc->horizon;
                    if ( sbsc->n_requests[j]  != 0 ) {
                        break;
                    }
                } 
                common_step += sbsc->current;
            }
        }

        if ( common_step != bsc_flush ) {
            for ( i = 0; i < P; ++i ) 
                bsp_put( i, &common_step, sbsc->next_step, bsp_pid()*sz, sz );
        }
        bsp_sync();
        common_step = until;
        for ( i = 0; i < P; ++i ) {
            if ( sbsc->next_step[i] < common_step )
                common_step = sbsc->next_step[i];
        }

        if (sbsc->horizon > 0) {
            /* execute the delayed local reduce operations */
            for ( i = 0; i < sbsc->n_requests[sbsc->queue_start]; ++i ) {
                request_t * r = sbsc->queue + 
                    sbsc->queue_start * sbsc->max_requests_per_step + i ;
                if ( r->kind == EXEC ) {
                    bsc_reduce_t f = r->payload.exec.func;
                    void * a = r->payload.exec.a;
                    const void * a0 = r->payload.exec.a0 ;
                    const void * xs = r->payload.exec.xs;
                    bsc_size_t size = r->payload.exec.size;
                    (*f)( a, a0, xs, size );
                }
            }
            sbsc->n_requests[ sbsc->queue_start ] = 0;
            sbsc->queue_start = (sbsc->queue_start+common_step-sbsc->current) % sbsc->horizon;
        }
        if ( common_step == bsc_flush ) {
            assert( sbsc->total_n_outstanding_requests == 0 );
            sbsc->current = 0;
            sbsc->queue_start = 0;
            break;
        }
        sbsc->current = common_step ;
    }

    return sbsc->current;
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

double bsc_costs_std_1phase_put( bsc_group_t group,
         bsc_coll_params_t * set, bsc_size_t n )
{
    bsc_size_t i;
    group_t * g = group;
    bsc_pid_t P = g?g->size:bsp_nprocs() ;
    double costs = bsc_L(BSC_PUT);
    for ( i = 0; i < n ; ++i )
        costs += bsc_g(BSC_PUT, set[i].size) * P * set[i].size;

    return costs;
}

double bsc_costs_std_1phase_get( bsc_group_t group,
         bsc_coll_params_t * set, bsc_size_t n )
{
    bsc_size_t i;
    group_t * g = group;
    bsc_pid_t P = g?g->size:bsp_nprocs() ;
    double costs = bsc_L(BSC_GET);
    for ( i = 0; i < n ; ++i )
        costs += bsc_g(BSC_GET, set[i].size) * P * set[i].size;

    return costs;
}


double bsc_costs_std_2tree_put(  bsc_group_t group,
         bsc_coll_params_t * set, bsc_size_t n )
{
    bsc_size_t i;
    group_t * g = group;
    bsc_pid_t P = g?g->size:bsp_nprocs() ;
    double one_phase = bsc_L(BSC_PUT);
    for ( i = 0; i < n ; ++i )
        one_phase += bsc_g(BSC_PUT, set[i].size) * set[i].size;

    return P==1?0:(int_log2(P-1)+1) * one_phase ;
}

double bsc_costs_std_2tree_get(  bsc_group_t group,
         bsc_coll_params_t * set, bsc_size_t n )
{
    bsc_size_t i;
    group_t * g = group;
    bsc_pid_t P = g?g->size:bsp_nprocs() ;
    double one_phase = bsc_L(BSC_GET);
    for ( i = 0; i < n ; ++i )
        one_phase += bsc_g(BSC_GET, set[i].size) * set[i].size;

    return P==1?0:(int_log2(P-1)+1) * one_phase ;
}


double bsc_costs_std_root_p_tree_put(  bsc_group_t group,
         bsc_coll_params_t * set, bsc_size_t n )
{
    bsc_size_t i;
    group_t * g = group;
    bsc_pid_t P = g?g->size:bsp_nprocs() ;
    bsc_pid_t root_P = ceil(sqrt(P));
    double one_phase = bsc_L(BSC_PUT); 
    for ( i = 0; i < n ; ++i )
        one_phase += (root_P - 1) * bsc_g(BSC_PUT, set[i].size) * set[i].size;

    return 2*one_phase;
}

double bsc_costs_std_root_p_tree_get(  bsc_group_t group,
         bsc_coll_params_t * set, bsc_size_t n )
{
    bsc_size_t i;
    group_t * g = group;
    bsc_pid_t P = g?g->size:bsp_nprocs() ;
    bsc_pid_t root_P = ceil(sqrt(P));
    double one_phase = bsc_L(BSC_GET); 
    for ( i = 0; i < n ; ++i )
        one_phase += (root_P - 1) * bsc_g(BSC_GET, set[i].size) * set[i].size;

    return 2*one_phase;
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
    return P==1?1:int_log2(P-1)+1;
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

const bsc_collective_t bsc_scatter_algs = {
    { bsc_scatter_multiple },
    { bsc_costs_std_1phase_put},
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

const bsc_collective_t bsc_gather_algs = {
    { bsc_gather_multiple },
    { bsc_costs_std_1phase_get },
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

const bsc_collective_t bsc_allgather_algs = {
    { bsc_allgather_multiple },
    { bsc_costs_std_1phase_get},
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

const bsc_collective_t bsc_alltoall_algs = {
    { bsc_alltoall_multiple },
    { bsc_costs_std_1phase_put},
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

    if ( me == 0 && dst != src ) 
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
    bsc_size_t b, my_start_b=-1, my_block_size=-1;
    bsc_size_t j, my_start_j=-1;
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

double bsc_bcast_2phase_costs(  bsc_group_t group,
         bsc_coll_params_t * set, bsc_size_t n )
{
    bsc_size_t i;
    group_t * g = group;
    bsc_pid_t P = g?g->size:bsp_nprocs() ;
    double one_phase = bsc_L(BSC_PUT);

    /* Note: the following is only an estimate */

    if ( n < P ) {
        for ( i = 0; i < n ; ++i ) 
            one_phase += bsc_g(BSC_PUT, (set[i].size+P-1)/P) * set[i].size ;
    } else {
        for ( i = 0 ; i < n; ++i )
            one_phase += bsc_g(BSC_PUT, set[i].size) * set[i].size;
    }

    return 2.0 * one_phase;
}


const bsc_collective_t bsc_bcast_algs = {
    { bsc_bcast_1phase_multiple, 
      bsc_bcast_2tree_multiple,
      bsc_bcast_root_p_tree_multiple,
      bsc_bcast_2phase },
    { bsc_costs_std_1phase_put,
      bsc_costs_std_2tree_put,
      bsc_costs_std_root_p_tree_put,
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
    void * accum = tmp_bytes;
    char * gather = tmp_bytes + size;

    (*reducer)( P==1?dst:accum, zero, src, size * nmemb );
    if (P == 1) {
        return depends;
    }
    
    for (range = int_pow( q, int_log( q, P-1));
            range >= 1; range /= q ) {
        if ( me < range ) {
            k = 1;
            for ( i = 1 ; i < q; ++i) {
                j = i * range + me; 
                if (j < P) {
                    j = (j + r )% P;
                    bsc_get( depends, g?g->gid[j]:j, accum, 0, 
                            &gather[k*size], size );
                    k += 1;
                }
            }
            bsc_get( depends, bsp_pid(), accum, 0, &gather[0], size );
            bsc_exec_reduce( depends, reducer, 
                    range==1?dst:accum, zero, gather, size * k );
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

const bsc_collective_t bsc_reduce_algs = {
    { bsc_reduce_1phase_multiple, 
      bsc_reduce_2tree_multiple,
      bsc_reduce_root_p_tree_multiple,
    },      
    { bsc_costs_std_1phase_get,
      bsc_costs_std_2tree_get,
      bsc_costs_std_root_p_tree_get,
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


bsc_step_t bsc_allreduce_qtree_single( 
        bsc_step_t depends, bsc_group_t group,
        const void * src, void * dst, void * tmp_space,
        bsc_reduce_t reducer, const void * zero,
        bsc_size_t nmemb, bsc_size_t size, 
        bsc_pid_t q )
{
    const group_t * g = group;
    bsc_pid_t i, j, k, P = g?g->size:bsp_nprocs() ;
    bsc_pid_t range, s = g?g->lid[bsp_pid()]:bsp_pid();
    char * tmp_bytes = tmp_space;
    void * accum = tmp_bytes;
    char * gather = tmp_bytes + size;

    (*reducer)( P==1?dst:accum, zero, src, size * nmemb );

    for (range = 1; range < P; range *= q ) {
        k = 1;
        for ( i = 1 ; i < q; ++i) {
            if (i * range < P) {
                j = (i * range + s)%P; 
                bsc_get( depends, g?g->gid[j]:j, accum, 0,
                        &gather[k * size], size);
                k+=1;
            }
        }
        bsc_get( depends, bsp_pid(), accum, 0, gather, size );
        bsc_exec_reduce( depends, reducer, range*q>=P?dst:accum, zero, gather, size*k );
        depends += 1;
    }
    return depends;
}

bsc_step_t bsc_allreduce_1phase_multiple( bsc_step_t depends,
                      bsc_pid_t root, bsc_group_t group, 
                      bsc_coll_params_t * params, bsc_size_t n)
{
    group_t * g = group;
    bsc_pid_t P = g?g->size:bsp_nprocs() ;
    bsc_size_t i;
    bsc_step_t finished = depends;
    (void) root;
    for ( i = 0; i < n ; ++i ) {
        finished = bsc_allreduce_qtree_single( depends, group,
                   params[i].src, params[i].dst, params[i].tmp,
                   params[i].reducer, params[i].zero, 
                   params[i].nmemb, params[i].size, P );
    }

    return finished;
}

bsc_step_t bsc_allreduce_2tree_multiple( bsc_step_t depends,
                      bsc_pid_t root, bsc_group_t group, 
                      bsc_coll_params_t * params, bsc_size_t n)
{
    bsc_size_t i;
    bsc_step_t finished = depends;
    (void) root;
    for ( i = 0; i < n ; ++i ) {
        finished = bsc_allreduce_qtree_single( depends, group,
                   params[i].src, params[i].dst, params[i].tmp,
                   params[i].reducer, params[i].zero, 
                   params[i].nmemb, params[i].size, 2 );
    }

    return finished;
}

bsc_step_t bsc_allreduce_root_p_tree_multiple( bsc_step_t depends, 
                      bsc_pid_t root, bsc_group_t group, 
                      bsc_coll_params_t * params, bsc_size_t n)
{
    group_t * g = group;
    bsc_pid_t P = g?g->size:bsp_nprocs(), root_P = ceil(sqrt(P)) ;
    bsc_size_t i;
    bsc_step_t finished = depends;
    (void) root;
    for ( i = 0; i < n ; ++i ) {
        finished = bsc_allreduce_qtree_single( depends, group,
                   params[i].src, params[i].dst, params[i].tmp,
                   params[i].reducer, params[i].zero, 
                   params[i].nmemb, params[i].size, root_P );
    }

    return finished;
}

const bsc_collective_t bsc_allreduce_algs = {
    { bsc_allreduce_1phase_multiple, 
      bsc_allreduce_2tree_multiple,
      bsc_allreduce_root_p_tree_multiple,
    },      
    { bsc_costs_std_1phase_get,
      bsc_costs_std_2tree_get,
      bsc_costs_std_root_p_tree_get,
    },
    { bsc_steps_1phase,
      bsc_steps_2tree ,
      bsc_steps_2phase,
    },
    3 
};
const bsc_collective_t * const bsc_coll_allreduce = & bsc_allreduce_algs;



bsc_step_t bsc_allreduce( bsc_step_t depends, bsc_group_t group,
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
    return bsc_collective( depends, bsc_coll_allreduce, 0, group, p );
}

bsc_step_t bsc_scan_qtree_single( bsc_step_t depends, bsc_group_t group,
        const void * src, void * dst, void * tmp_space, 
        bsc_reduce_t reducer, const void * zero,
        bsc_size_t nmemb, bsc_size_t size,
        int q )
{
    const group_t * g = group;
    bsc_pid_t i, j, k=0, P = g?g->size:bsp_nprocs() ;
    bsc_pid_t range, s = g?g->lid[bsp_pid()]:bsp_pid();
    char * dst_bytes = dst;
    char * tmp1 = tmp_space;
    char * tmp2 = tmp1 + P * size;
    const void * prefix;

    (*reducer)( dst, zero, src, size * nmemb );

    for (range = 1; range < P; range *= q ) {
        const bsc_size_t shift = (range==1);
        const void * last = NULL;
        if (range == 1)
            last = (nmemb==0 ? zero : &dst_bytes[(nmemb-1)*size]);
        else
            last = (k==0 ? zero : &tmp2[(k-1)*size]);

        for ( i = 1; i < q; ++i ) {
            j = s + range * i + shift;
            if ( j < P ) {
                bsc_put( depends, g?g->gid[j]:j, last, tmp1, i * size, size);
            }
        }
    
        /* k is number of received elements */
        k = MIN( q, (range+s-shift)/range);

        if (s < P - shift) {
            j = g?g->gid[s+shift]:s+shift;
            bsc_put(depends, j, last, tmp1, 0, size);
        }
        bsc_exec_reduce( depends, reducer, tmp2, zero, tmp1, size*k );
        depends += 1;
    }

    prefix = k > 0 ? &tmp2[(k-1)*size] : zero ;
    bsc_put(depends, bsp_pid(), prefix, tmp1, 0, size);

    bsc_exec_reduce( depends, reducer, dst, tmp1, src, size*nmemb );
    return depends+1;
}

bsc_step_t bsc_scan_1phase_multiple( bsc_step_t depends,
                      bsc_pid_t root, bsc_group_t group, 
                      bsc_coll_params_t * params, bsc_size_t n)
{
    group_t * g = group;
    bsc_pid_t P = g?g->size:bsp_nprocs() ;
    bsc_size_t i;
    bsc_step_t finished = depends;
    (void) root;
    for ( i = 0; i < n ; ++i ) {
        finished = bsc_scan_qtree_single( depends, group,
                   params[i].src, params[i].dst, params[i].tmp,
                   params[i].reducer, params[i].zero, 
                   params[i].nmemb, params[i].size, P );
    }

    return finished;
}

bsc_step_t bsc_scan_2tree_multiple( bsc_step_t depends,
                      bsc_pid_t root, bsc_group_t group, 
                      bsc_coll_params_t * params, bsc_size_t n)
{
    bsc_size_t i;
    bsc_step_t finished = depends;
    (void) root;
    for ( i = 0; i < n ; ++i ) {
        finished = bsc_scan_qtree_single( depends, group,
                   params[i].src, params[i].dst, params[i].tmp,
                   params[i].reducer, params[i].zero, 
                   params[i].nmemb, params[i].size, 2 );
    }

    return finished;
}

bsc_step_t bsc_scan_root_p_tree_multiple( bsc_step_t depends,
                      bsc_pid_t root, bsc_group_t group, 
                      bsc_coll_params_t * params, bsc_size_t n)
{
    group_t * g = group;
    bsc_pid_t P = g?g->size:bsp_nprocs(), root_P = ceil(sqrt(P)) ;
    bsc_size_t i;
    bsc_step_t finished = depends;
    (void) root;
    for ( i = 0; i < n ; ++i ) {
        finished = bsc_scan_qtree_single( depends, group,
                   params[i].src, params[i].dst, params[i].tmp,
                   params[i].reducer, params[i].zero, 
                   params[i].nmemb, params[i].size, root_P );
    }

    return finished;
}

const bsc_collective_t bsc_scan_algs = {
    { bsc_scan_1phase_multiple, 
      bsc_scan_2tree_multiple,
      bsc_scan_root_p_tree_multiple,
    },      
    { bsc_costs_std_1phase_put,
      bsc_costs_std_2tree_put,
      bsc_costs_std_root_p_tree_put,
    },
    { bsc_steps_1phase,
      bsc_steps_2tree ,
      bsc_steps_2phase,
    },
    3 
};
const bsc_collective_t * const bsc_coll_scan = & bsc_scan_algs;



bsc_step_t bsc_scan( bsc_step_t depends, bsc_group_t group,
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
    return bsc_collective( depends, bsc_coll_scan, 0, group, p );
}
