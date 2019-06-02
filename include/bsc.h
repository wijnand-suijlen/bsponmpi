#ifndef BSPONMPI_BSC_H
#define BSPONMPI_BSC_H

#include <bsp.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

/* copied from https://gcc.gnu.org/wiki/Visibility */
#if defined DOXYGEN
  #define DLL_PUBLIC
  #define DLL_LOCAL
#elif defined _WIN32 || defined __CYGWIN__
  #ifdef BUILDING_DLL
    #ifdef __GNUC__
      #define DLL_PUBLIC __attribute__ ((dllexport))
    #elif defined __clang__
      #define DLL_PUBLIC __attribute__ ((dllexport))
    #else
      #define DLL_PUBLIC __declspec(dllexport)
    #endif
  #else
    #ifdef __GNUC__
      #define DLL_PUBLIC __attribute__ ((dllimport))
    #elif defined __clang__
      #define DLL_PUBLIC __attribute__ ((dllimport))
    #else
      #define DLL_PUBLIC __declspec(dllimport) 
    #endif
  #endif
  #define DLL_LOCAL
#else
  #if __GNUC__ >= 4 || defined __clang__
    #define DLL_PUBLIC __attribute__ ((visibility ("default")))
    #define DLL_LOCAL  __attribute__ ((visibility ("hidden")))
  #else
    #define DLL_PUBLIC
    #define DLL_LOCAL
  #endif
#endif


/** \defgroup BSC Bulk Synchronous Collectives 
 *
 * This library defines a collectives library for BSPlib. It differs from the
 * level-1 library from the Oxford Toolset in that it allows collective
 * operations on process subgroups. This requires some sophistication, however,
 * because BSPlib's global barrier synchronisations require active
 * participation from all processes, including the processes that do not
 * participate in the collective. As a solution, this library choses to augment
 * BSPlib with <i>delayed</i> requests. Instead of requiring execution of the
 * primitives bsp_put() and bsp_get() during this superstep, this library
 * defines a shell that allows scheduling of their execution in later
 * supersteps. This enables each process in each subgroup to schedule group
 * private communications in a local queue completely asynchronously. A special
 * bsc_sync() will then play the queued requests and call bsp_sync()
 * accordingly. Depending on the user's wishes the bsc_sync() can either
 * synchronize globally or synchronize only the subgroup engaged in a
 * particular collective. 
 *
 * Collectives are built by building a dependency DAG at runtime, shortly
 * before the collective require execution. Each node in the graph is superstep
 * number. Edges are formed by calls to primitives bsc_put() and bsc_get(). As
 * the first parameter they take the start superstep while their return value
 * is the end superstep. In case of the primitives the end superstep is always
 * equal to the start superstep incremented by one. Collectives follow the same
 * calling convention, but their end superstep may be farther away that one. This
 * convention allows simultaneous execution of collectives by using the 
 * same start superstep <i>or</i> forces serial execution by using the return value
 * of the preceding collective as the first parameter in the next collective.
 * 
 * The implementation of the two-phase broadcast illustrates how complex
 * collectives can be built from this.
 *
 * \code
    bsc_step_t bsc_bcast( bsc_step_t depends,
                          bsc_pid_t root, bsc_group_t group, 
                          const void * src, void * dst, bsc_size_t size )
    {
        const group_t * g = group;
        bsp_pid_t i, n = g?g->size:bsp_nprocs() ;
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
 * \endcode
 *
 * The newly defined bsc_sync() will execute the request queue until the step
 * that is given as parameter. When the given time step has the special value
 * #bsc_flush, it will continue processing until all processes have finished
 * all enqueued requests and have also reached a bsc_sync() with a #bsc_flush.
 * It behaves therefore as a global barrier synchronisation, requiring textual
 * alignment. The following example comes from BSPedupack's bsplu. 
 *
 * \code
    bsc_group_t col_procs = bsc_group_create_partition( t );
    bsc_group_t row_procs = bsc_group_create_partition( s );
    if (k%N==t){  
        ...
        bsc_step_t pivot_known = bsc_allreduce( bsc_start, col_procs, ...  );
        bsc_sync( pivot_known );
        ...
    }
    bsc_sync( bsc_flush ); 
 * \endcode
 *
 *
 * @{
 *
 */

/** \defgroup BSC_TYPEDEFS Elementary types 
 *
 * @{
 * 
 */

/** An unsigned integer that can hold the delay in supersteps */
typedef unsigned bsc_step_t;

/** An integer to hold the size of data communication requests */
typedef bsp_size_t bsc_size_t;

/** An integer to hold process identifiers */
typedef bsp_pid_t bsc_pid_t;

/**
 * @}
 */

/** \defgroup BSC_GROUP Process group management
 *
 * Process groups restrict a collective to a subset of processes
 *
 * @{
 *
 */

/** An object representing a group of processes */
typedef void * bsc_group_t;

/** The group that holds all processes. */
DLL_PUBLIC extern const bsc_group_t bsc_all;

/** A collective function that partitions all processes in disjoint groups.
 * All processes must collectively engage in the call to this function, while
 * chosing a \a symbol of their liking. Processes that chose the same symbol
 * are put together in the same group. 
 *
 * \note This function synchronizes
 *
 * \param symbol The symbol identifying the subgroup to be partitioned in
 * 
 * \returns An object representing the group
 */
DLL_PUBLIC bsc_group_t bsc_group_create_partition( unsigned symbol );

/** A collective function that groups all processes in neighbourhoods. 
 * All processes must collectively engage in the call to this function, while
 * chosing their neighbours. Only neighbours that acknowledge each other
 * are put in the same group. The order in which the neighbours appear is significant,
 * as it will be the same order used for bsc_scatter(), bsc_gather(), bsc_allgather(),
 * an bsc_scan(). The neighbour list must include the calling process too.
 * 
 * \note This function synchronizes
 *
 * \param neighbours An array with the process IDs of the neighbours of this process
 *                   and itself.
 * \param size The number of neighbours.
 * 
 * \returns An object representing the group
 */
DLL_PUBLIC bsc_group_t bsc_group_create_neighbourhood( bsc_pid_t * neighbours, 
        bsc_size_t size );

/** Frees resources held to maintain bookkeeping for this group. 
 *
 * \param group The group to be destroyed.
 * 
 */
DLL_PUBLIC void bsc_group_destroy( bsc_group_t group );

/**
 * @}
 */

/** \defgroup BSC_INTROSP  BSP machine performance parameters
 *
 * These functions can be used by collectives implementations to decide which
 * BSP algorithm to use
 *
 * @{
 *
 */

typedef enum bsc_request { 
    BSC_PUT, BSC_HPPUT, BSC_GET, BSC_HPGET, BSC_SEND,
    BSC_N_REQUEST_TYPES
} bsc_request_t;

/** BSP cost parameter \f$g\f$, denoting the message gap or reciprocal
 * throughput. This can be used to choose the best algorithm, depending
 * on the input data size. The value is relative to bsp_L().
 *
 * \param type      The type of requests 
 * \param word_size The value of the \a size argument in the request.
 * \returns The BSP cost parameter \f$g\f$ for the specified type of requests
 *          and specified size.
 */
DLL_PUBLIC double bsc_g(bsc_request_t type, int word_size);

/** BSP cost parameter \f$\ell\f$, denoting the synchronisation latency.
 * This can be used to choose the best algorithm, depending
 * on the input data size. The value is relative to bsp_g().
 * \param type      The type of requests 
 * \returns The BSP cost parameter \f$\ell\f$ for the specified type of
 *          requests.
 */
DLL_PUBLIC double bsc_L(bsc_request_t type);

/**
 * @}
 */


/** \defgroup BSC_PRIM Delayed one-sided communication primitives
 *
 * All collectives are built from only these primitives
 *
 * @{
 *
 */


/** The special delay value to use with bsc_sync() to proceed execution of all
 * outstanding requests and to perform a global barrier synchronisation. This
 * ends the current epoch and starts a new one. */
DLL_PUBLIC extern const bsc_step_t bsc_flush;

/** The special delay value that denotes the first superstep of an epoch.
 * This number is equal to zero.
 *
 * \see #bsc_flush */
DLL_PUBLIC extern const bsc_step_t bsc_start;

/** The current superstep number. This is equal to zero at the start of an epoch.
 */
DLL_PUBLIC bsc_step_t bsc_current(void);

/** Completes requests until superstep \a until is started. If #bsc_flush is given,
 * it ends an epoch and starts a new one by completing all outstanding requests
 * and waiting for the other processes to also call bsc_sync() with #bsc_flush,
 * hence acting as a global barrier synchronisation. Using bsc_current() as
 * value, results in a no-op.  Any other value will proceed to that specific
 * superstep. In that case also a P-relation is incurred as extra BSP cost.
 *
 * \param until The superstep number to proceed to. If #bsc_flush is given, it
 *              will proceed until no more requests remain in the queue and
 *              other processes also have called bsc_sync() with #bsc_flush.
 *
 * \returns The superstep number after the synchronisation. This will be equal
 *          to \a until unless #bsc_flush was used, in which case it is zero.
 */
DLL_PUBLIC bsc_step_t bsc_sync( bsc_step_t until );


/** Requests a bsp_put() to be executed in \a depends supersteps from now.
 *
 * \param depends The superstep number in which to request the bsp_put().
 * \param dst_pid ID of the destination process
 * \param src     Pointer to local memory where the source data will be in 
 *                \a depends supersteps from now
 * \param dst     Pointer of the registered memory area where the data must be 
 *                sent to.
 * \param offset  Offset in bytes from \a dst where the data must be put.
 * \param size    Size in bytes of the data.
 *
 * \returns  The superstep number when this bsp_put() has completed
 */
DLL_PUBLIC bsc_step_t bsc_put( bsc_step_t depends, bsc_pid_t dst_pid, 
        const void * src, void * dst, bsc_size_t offset, bsc_size_t size);

/** Requests a bsp_get() to be executed in \a depends supersteps from now.
 *
 * \param depends The superstep number in which to request the bsp_get().
 * \param src_pid ID of the source process
 * \param src     Pointer to the registered memory area where the source data
 *                will be in \a depends supersteps from now
 * \param offset  Offset in bytes from \a src where the data must be read from.
 * \param dst     Pointer of local memory area where the data must be sent to.
 * \param size    Size in bytes of the data.
 
 * \returns  The superstep number when this bsp_get() has completed
*/
DLL_PUBLIC bsc_step_t bsc_get( bsc_step_t depends, bsc_pid_t src_pid, 
        const void * src, bsc_size_t offset, void * dst, bsc_size_t size);

/** A function pointer that points to an implementation of an associative
 * operator \f$\oplus\f$. It must perform the reduction on the array \a xs of
 * \a size bytes in total and store the result in \a a, either as a single value
 * or a prefix calculation (scan).  The value pointer to by \a a0 is used as
 * start value for \a a. The pointers \a a, \a a0, and \a xs may not refer to
 * overlapping memory blocks.
 *
 * \param a     Pointer to memory where to store the result
 * \param a0    Start value for \a a
 * \param xs    Pointer to the start of the array
 * \param size  Number of bytes that the array holds
 */
typedef void (*bsc_reduce_t)(void * a, const void * a0,
        const void * xs, bsc_size_t size );

/** Requests execution of a \a bsp_reduce_t function after all communication 
 * of the superstep of \a depends steps from now has been completed.
 *
 * \param depends The superstep number just before the step that executes this
 * \param reducer The reduction function to execute on the local array
 * \param a       Pointer to local memory where the result must be stored
 * \param a0      Pointer to local memory where start value for \a a is stored
 * \param xs      Pointer to local memory where the array is stored
 * \param size    Size in bytes of the local array.
 *
 * \returns  The superstep number when this reduction has completed
 */
DLL_PUBLIC bsc_step_t bsc_exec_reduce( bsc_step_t depends, 
        bsc_reduce_t reducer, void * a, const void * a0,
        const void * xs, bsc_size_t size );

/** Stores the data and operation to be executed by the collective.  This is
 * used to define collectives whose operation can be coalesced.
 */
typedef struct bsc_coll_params {
    const void * src; /**< pointer to source data */
    void * dst;       /**< pointer to destination data */
    void * tmp;       /**< pointer to temporary memory */
    bsc_reduce_t reducer; /**< function computing reduction or prefix */
    const void * zero; /**< the neutral element for a reduction */
    bsc_size_t nmemb; /**< number of elements */
    bsc_size_t size; /**< size of each element */
} bsc_coll_params_t;

/** A parallel collective algorithm for coalesced collective operations.
 *
 * \param depends Collective is executed at the end of this superstep
 * \param root    A root process, if relevant
 * \param group   The process group
 * \param set     A pointer to an array with the coalesced set of calls
 * \param n       The number of elements in \a set.
 *
 * \returns The superstep of completion
 *
 */
typedef bsc_step_t (*bsc_coll_alg_t)( bsc_step_t depends, 
        bsc_pid_t root, bsc_group_t group, 
        bsc_coll_params_t * set, bsc_size_t n );

/** The cost function for a parallel, coalesced collective algorithm.
 *
 * \param group The process group
 * \param set   A pointer to an array with the coalesced set of calls
 * \param n     The number of elements in \a set.
 *
 * \returns The cost of this algorithm
 */
typedef double (*bsc_coll_cost_t)( bsc_group_t group,
        bsc_coll_params_t * set, bsc_size_t n );

/** The number of supersteps required by a parallel, coalesced collective
 * algorithm
 *
 * \param group The process group
 *
 * \return The number of supersteps required
 */
typedef bsc_step_t (*bsc_coll_steps_t)(bsc_group_t group);

/** The maximum number of algorithms allowed for each collective */
#define BSC_MAX_COLL_ALGS 10

/** Defines a collective */
typedef struct bsc_collective { 
    
    /** Each collective algorithm */
    bsc_coll_alg_t   algorithms[BSC_MAX_COLL_ALGS];

    /** The cost function of each algorithm */
    bsc_coll_cost_t  costfuncs[BSC_MAX_COLL_ALGS];

    /** The number of supersteps required by each algorithm */
    bsc_coll_steps_t maxsteps[BSC_MAX_COLL_ALGS];

    /** The number of algorithms */
    bsc_size_t       n_algs;
} bsc_collective_t;

/** Schedules \a collective operation for superstep \a depends. This call
 * can be coalesced with other calls that have the same \a root, \a group, 
 * and \a depends.
 *
 * \param  depends The collective operation is schedules at the end of
 *                 superstep \a depends.
 * \param  collective The collective operation that is scheduled
 * \param  root    The root process of the collective operation, if applicable
 * \param  group   The process group 
 * \param  params  The collective operation parameters
 */
DLL_PUBLIC bsc_step_t bsc_collective( bsc_step_t depends, 
       const bsc_collective_t * collective,
       bsc_pid_t root, bsc_group_t group, 
       bsc_coll_params_t params ) ;


/**
 * @}
 *
 */



/** \defgroup BSC_COLLS Collectives 
 *
 * The library of communication collectives
 *
 * @{
 */


/** Schedules a collective scatter operation. This is a collective call
 * on each process in the given \a group.
 *
 * \f[
 * \forall_{i \in \texttt{group}} \quad
 * \texttt{dst}_i[ 0 \ldots \texttt{size} ] 
 *      \leftarrow \texttt{src}_\texttt{root}[ i \texttt{size} \ldots (i+1) \texttt{size} ]
 * \f]
 *
 * \param       depends <i>collective</i> The superstep number to schedule the
 *                      first communications
 * \param       group   <i>collective</i> The group to perform this collective
 * \param       root    <i>collective</i> The ID of the process that holds the
 *                      data to be scattered.
 *
 * \param       src     Pointer to the source memory, local to the \a root.
 * \param       dst     Pointer to <i>registered</i> destination memory.
 * \param       size    <i>collective</i> Size of each chunk to be scattered. 
 *
 * \returns The superstep number that this collective will have been completed.
 */
DLL_PUBLIC bsc_step_t bsc_scatter( bsc_step_t depends,
                        bsc_pid_t root, bsc_group_t group,
                       const void * src, void *dst, bsc_size_t size );

/** Schedules a collective gather operation. This is a collective call
 * on each process in the given \a group.
 *
 * \f[
 * \forall_{i \in \texttt{group}} \quad
 * \texttt{dst}_\texttt{root}[ i \texttt{size} \ldots (i+1) \texttt{size} ]
 *      \leftarrow \texttt{src}_i[ 0 \ldots \texttt{size} ] 
 * \f]
 *
 * \param depends <i>collective</i> The superstep number to schedule the first
 *                 communications
 * \param group   <i>collective</i> The group to perform this collective
 * \param root    <i>collective</i> The ID of the process on which all data is
 *                gathered.
 * \param src     Pointer to <i>registered</i> source memory.
 * \param dst     Pointer to destinatinon memory, local to the \a root.
 * \param size    <i>collective</i> Size of each chunk to be gathered. 
 *
 * \returns The superstep number that this collective will have been completed.
 */
DLL_PUBLIC bsc_step_t bsc_gather( bsc_step_t depends, 
                       bsp_pid_t root, bsc_group_t group,
                       const void * src, void *dst, bsc_size_t size );

/** Schedules a collective all-gather operation. This is a collective call
 * on each process in the given \a group.
 *
 * \f[
 * \forall_{i,j \in \texttt{group}} \quad
 * \texttt{dst}_j[ i \texttt{size} \ldots (i+1) \texttt{size} ]
 *      \leftarrow \texttt{src}_i[ 0 \ldots \texttt{size} ] 
 * \f]
 *
 * \param depends <i>collective</i> The superstep number to schedule the first
 *                communications
 * \param group   <i>collective</i> The group to perform this collective
 * \param src     Pointer to <i>registered</i> source memory.
 * \param dst     Pointer to the destination memory.
 * \param size    <i>collective</i> Size of each chunk to be gathered. 
 *
 * \returns The superstep number that this collective will have been completed.
 */
DLL_PUBLIC bsc_step_t bsc_allgather( bsc_step_t depends, bsc_group_t group,
                        const void * src, void * dst, bsc_size_t size );

/** Schedules a collective total-exchange operation. This is a collective call
 * on each process in the given \a group.
 *
 * \f[
 * \forall_{i,j \in \texttt{group}} \quad
 * \texttt{dst}_j[ i \texttt{size} \ldots (i+1) \texttt{size} ]
 *      \leftarrow \texttt{src}_i[ j \texttt{size} \ldots (j+1) \texttt{size}] 
 * \f]
 *
 * \param depends <i>collective</i> The superstep number to schedule the first
 *                communications
 * \param group   <i>collective</i> The group to perform this collective
 * \param src     Pointer to source memory, local to each process.
 * \param dst     Pointer to the <i>registered</i> destination memory.
 * \param size    <i>collective</i> Size of each chunk to be gathered. 
 *
 * \returns The superstep number that this collective will have been completed.
 */
DLL_PUBLIC bsc_step_t bsc_alltoall( bsc_step_t depends, bsc_group_t group,
                        const void * src, void * dst, bsc_size_t size );

/** Schedules a collective broadcast operation. This is a collective call
 * on each process in the given \a group.
 *
 * \f[
 * \forall_{i \in \texttt{group}} \quad
 * \texttt{dst}_i[ 0 \ldots  \texttt{size} ]
 *      \leftarrow \texttt{src}_\texttt{root}[ 0 \ldots \texttt{size} ] 
 * \f]
 *
 * \param depends <i>collective</i> The superstep number to schedule the first
 *                communications
 * \param root    <i>collective</i> The process that broadcasts 
 * \param group   <i>collective</i> The group to perform this collective
 * \param src     Pointer to source memory, local to the \a root process.
 * \param dst     Pointer to the <i>registered</i> destination memory.
 * \param size    <i>collective</i> Size of the data.
 *
 * \returns The superstep number that this collective will have been completed.
 */
DLL_PUBLIC bsc_step_t bsc_bcast( bsc_step_t depends,
                      bsp_pid_t root, bsc_group_t group, 
                      const void * src, void * dst, bsc_size_t size );


/** Schedules a collective reduction operation and stores the result on the 
 * \a root process. This is a collective call on each process in the given \a
 * group.
 *
 * \f[
 * \texttt{dst}_\texttt{root}
 *      \leftarrow \bigoplus_{s \in \texttt{group}} 
 *      \bigoplus_{k=1}^{\texttt{nmemb}_s}
 *      \texttt{src}_{s,k-1} 
 * \f]
 *
 * \param depends <i>collective</i> The superstep number to schedule the first
 *                      communications
 * \param root    <i>collective</i> The process that will get the result in \a dst.
 * \param group   <i>collective</i> The group to perform this collective
 * \param src     Pointer to locally stored array with \a nmemb objects 
 *                      of  \a size bytes size.
 * \param dst     Pointer to memory of \a size bytes, local to root process.
 * \param tmp_space Pointer to <i>registered</i> scratch space 
 *                of \f$(p+1) \cdot \texttt{size}\f$ bytes size.
 * \param reducer <i>collective</i> The reduction function.
 * \param zero    <i>collective</i> The zero for the reduction function.
 *                It must hold that \f$ \texttt{zero} \oplus x = x, \ \forall x \f$.
 * \param nmemb   The number of elements in the local input array \a src
 *                      on this process
 * \param size    <i>collective</i> The size of each element in bytes.
 *
 * \returns The superstep number that this collective will have been completed.
 */
DLL_PUBLIC bsc_step_t bsc_reduce( bsc_step_t depends, 
        bsc_pid_t root, bsc_group_t group,
        const void * src, void * dst, void * tmp_space,
        bsc_reduce_t reducer, const void * zero, 
        bsc_size_t nmemb, bsc_size_t size );

/** Schedules a collective reduction operation and stores the result on each
 * process. This is a collective call on each process in the given \a group.
 *
 * \f[
 * \forall_{j \in \texttt{group}} \quad
 * \texttt{dst}_j
 *      \leftarrow \bigoplus_{s \in \texttt{group}} 
 *      \bigoplus_{k=1}^{\texttt{nmemb}_s}
 *      \texttt{src}_{s,k-1} 
 * \f]
 *
 * \param depends <i>collective</i>The superstep number to schedule the first
 *                      communications
 * \param group   <i>collective</i> The group to perform this collective
 * \param src     Pointer to locally stored array with \a nmemb objects 
 *                      of  \a size bytes size.
 * \param dst     Pointer to memory of \a size bytes
 *                where the result should be stored.
 * \param tmp_space Pointer to <i>registered</i> scratch space of
 *                      \f$(p + 1) \cdot \texttt{size}\f$ bytes size.  
 * \param reducer <i>collective</i> The reduction function.
 * \param zero    <i>collective</i> The zero for the reduction function.
 *                It must hold that \f$ \texttt{zero} \oplus x = x, \ \forall x \f$.
 * \param nmemb   The number of elements in the local input array \a src
 *                      on this process
 * \param size    <i>collective</i> The size of each element in bytes.
 *
 * \returns The superstep number that this collective will have been completed.
 */
DLL_PUBLIC bsc_step_t bsc_allreduce( bsc_step_t depends, bsc_group_t group,
        const void * src, void * dst, void * tmp_space,
        bsc_reduce_t reducer, const void * zero,
        bsc_size_t nmemb, bsc_size_t size );

/** Schedules a collective scan operation. This is a collective call on each
 * process in the given \a group.
 *
 * \f[
 * \forall_{j \in \texttt{group}} \quad \forall_i: 0 \leq i < \texttt{nmemb}_j \quad
 * \texttt{dst}_{j, i}
 *      \leftarrow \left[ \bigoplus_{s=0}^{j-1} 
 *      \bigoplus_{k=1}^{\texttt{nmemb}_s} 
 *      \texttt{src}_{s,k-1}  \right]
 *      \
 *      \oplus
 *      \
 *      \left[
 *      \bigoplus_{k=0}^{i} 
 *      \texttt{src}_{j,k} 
 *      \right]
 * \f]
 *
 * \param depends <i>collective</i>The superstep number to schedule the first
 *                      communications
 * \param group   <i>collective</i> The group to perform this collective
 * \param src     Pointer to locally stored array with \a nmemb objects 
 *                      of  \a size bytes size.
 * \param dst     Pointer to memory of array \a nmemb elements 
 *                of \a size bytes, where the result should be stored.
 * \param tmp_space Pointer to <i>registered</i> scratch space of
 *                      \f$2*p \cdot \texttt{size}\f$ bytes size.  
 * \param reducer <i>collective</i> The reduction function that computes a scan.
 * \param zero    <i>collective</i> The zero for the reduction function.
 *                It must hold that \f$ \texttt{zero} \oplus x = x, \ \forall x \f$.
 * \param nmemb   The number of elements in the input array \a src and output
 *                array \a dst  on this process
 * \param size    <i>collective</i> The size of each element in bytes.
 *
 * \returns The superstep number that this collective will have been completed.
 */
DLL_PUBLIC bsc_step_t bsc_scan( bsc_step_t depends, bsc_group_t group,
        const void * src, void * dst, void * tmp_space, 
        bsc_reduce_t reducer, const void * zero,
        bsc_size_t nmemb, bsc_size_t size );


/**
 * @}
 */


/** \defgroup BSC_COLLS_BLOCKS Parallel collective operation algorithms 
 *
 * The algorithms used to implement the collective operations
 *
 * @{
 */



/** Definition for bsc_scatter() */
DLL_PUBLIC extern const bsc_collective_t * const bsc_coll_scatter;

/** Definition for bsc_gather() */
DLL_PUBLIC extern const bsc_collective_t * const bsc_coll_gather;

/** Definition for bsc_allgather() */
DLL_PUBLIC extern const bsc_collective_t * const bsc_coll_allgather;

/** Definition for bsc_alltoall() */
DLL_PUBLIC extern const bsc_collective_t * const bsc_coll_alltoall;

/** Definition for bsc_bcast */
DLL_PUBLIC extern const bsc_collective_t * const bsc_coll_bcast;

/** Definition for bsc_reduce */
DLL_PUBLIC extern const bsc_collective_t * const bsc_coll_reduce;

/** Definition for bsc_allreduce */
DLL_PUBLIC extern const bsc_collective_t * const bsc_coll_allreduce;

/** Definition for bsc_scan */
DLL_PUBLIC extern const bsc_collective_t * const bsc_coll_scan;


/**
 * @}
 */


/**
 * @}
 */


#ifdef __cplusplus
}
#endif


#endif
