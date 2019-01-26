#ifndef BSPONMPI_BSP_H
#define BSPONMPI_BSP_H

#include <stdarg.h>
#include <stddef.h>

#ifdef __cplusplus 
extern "C" {
#endif

/** \mainpage BSPonMPI - A \ref BSPlib implementation on top of MPI
 *
 * \author Wijnand J. Suijlen
 * \copyright MIT License
 *
 * \section Introduction Introduction
 * BSPonMPI offers on any parallel machine supporting MPI a parallel
 * programming API that adheres to the BSP cost-model. This API is better known
 * as \ref BSPlib and supports Single Program Multiple Data (\ref BSP_SPMD "SPMD")
 * style programs that communicate using Direct Remote Memory Access (\ref
 * BSP_DRMA "DRMA") and Bulk Synchronous Message Passing (\ref BSP_BSMP
 * "BSMP"). This library enables <i>immortal algorithms</i> that once written
 * are scalable everywhere, since the BSP-cost model prescribes that a
 * program's run-time is bounded by 
 * \f[ W + g H + \ell S, \f]
 * where \f$W\f$, the local computation, \f$H\f$, the total communiation
 * volume, and \f$S\f$, the number of supersteps, are determined by the program
 * and \f$g\f$, the <i>message gap</i> (i.e. reciprocal throughput), and
 * \f$\ell\f$, the latency, are determined by the machine. 
 *
 * The API consists of three parts:
 *   - \ref BSP_SPMD "The SPMD Framework"
 *   - \ref BSP_DRMA "Direct Remote Memory Access communications"
 *   - \ref BSP_BSMP "Bulk Sychronous Message Passing communications"
 *
 * In each of these sections, examples are provided how they can be used.
 *
 * Additionally a \ref MCBSP "MulticoreBSP for C compatibility" layer is
 * provided. Use the --mcbsp parameter using the \c bspcc frontend or
 * define the macro \c BSPONMPI_MCBSP_COMPAT to enable this.
 *
 * \section Example Hello World
 * This example is shamelessly copied from [1].
 * \code
   #include <stdio.h>
   #include <bsp.h>
  
   void main(void) { 
     bsp_begin( bsp_nprocs() );
       printf("Hello BSP Worldwide from process %d of %d\n",
            bsp_pid(), bsp_nprocs() );
     bsp_end();
   }
   \endcode
 *
 * \section References References
 * \anchor BSPlib [1] "BSPlib: The BSP programming library," by J. M. D. Hill,
 * W. F. McColl, D. C. Stefanescu, M. W. Goudreau, K. Lang, S. B. Rao, T. Suel,
 * Th. Tsantilas, R. H. Bisseling, Elsevier, Parallel Computing, Volume 24, 
 * Issue 14, December 1998, pages 1947--1980. 
 */

/** \defgroup BSP_TYPES BSPlib convenience typedefs
 *
 * These are used as convenience for the programmer to maker her/his
 * program BSPlib dialect oblivious.
 *
 * By default this implementation is compatible with \ref BSPlib 1998.
 * See the section about the \ref MCBSP "MulticoreBSP for C compatibilty"  
 * layer how to switch modes. 
 *
 * @{ 
 */

/** Data type large enough to store the number of processes */
typedef int bsp_pid_t;

/** Data type large enough to express the volume in any communication request */
typedef int bsp_size_t;

/**
  * @}
  *
  */


#ifdef _GNUC_
  #define BSPONMPI_PRINTF_FORMAT_ATTRIBUTE(format_idx, first_va)\
                __attribute__(( format(printf, format_idx, first_va) ))
#else

  #define BSPONMPI_PRINTF_FORMAT_ATTRIBUTE(format_idx, first_va) /* empty */

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

/** \defgroup BSP_SPMD BSPlib SPMD Framework
 *  The SPMD framework enables writing programs in Single Program Multiple Data
 * style. This means that the runtime system will start the same program
 * multiple times and will connect the resulting processes through a communications
 * system (see \ref BSP_DRMA and \ref BSP_BSMP). A process can identify itself
 * only through the the bsp_pid() function, which returns a unique value \f$ 0
 * \leq \texttt{bsp\_pid()} < \texttt{bsp\_nprocs()}\f$ for each process.
 *
 * Because the whole program need not to be SPMD, the programmer should enclose
 * the SPMD part of her program in a function which has as first statement
 * bsp_begin() and as last statement a call to bsp_end(). This function may be
 * main(), but if it is not, the main() function must call bsp_init() as first
 * statement with a reference to the SPMD function. 
 *
 * Example 1
   \code
   #include <stdio.h>
   #include <bsp.h>
  
   void main(void) { 
     bsp_begin( bsp_nprocs() );
       printf("Hello BSP Worldwide from process %d of %d\n",
            bsp_pid(), bsp_nprocs() );
     bsp_end();
   }
   \endcode
  
  
   Example 2
   \code
   #include <stdio.h>
   #include <bsp.h>
  
   static int s_nprocs;
   
   void spmd_part( void ) {
     bsp_begin( s_nprocs );
       printf("Hello BSP Worldwide from process %d of %d\n",
              bsp_pid(), bsp_nprocs() );
     bsp_end();
   }
  
   void main(int argc, char *argv[]) {
      bsp_init( spmd_part, argc, argv);
      printf("Please enter number of processes: ");
      scanf("%d", &s_nprocs );
      spmd_part();
   }
   \endcode
 *
 *  @{
 */


/** Starts an SPMD section with at most \a maxprocs parallel processes. This
 * must be called as the first statement (except for a call to bsp_nprocs()) of
 * a function, which may be main(). In that same function there must also be
 * call to bsp_end() as the last statement. There may be only one instance of a
 * bsp_begin()/bsp_end() pair within a program. If the enclosing function is
 * not main(), then its reference must be passed as parameter to a call to
 * bsp_init(), which in its turn must be the first statement in main().
 *
 * In the SPMD section normal I/O operations are only guaranteed to work on the
 * process with ID 0. Only output to standard output and standard error are
 * guaranteed to work from all other processes. These outputs will be
 * multiplexed in arbitrary order.
 *
 * Additionally, only the process with ID 0 will be the continuation of the
 * calling process, which means that all variables on the other processes will
 * be undefined  after the call to bsp_begin(). The same holds for the (Unix)
 * environment. 
 *
 * \param maxprocs The requested number of processes by the user.
 * 
 * \throws bsp_abort When a 0 or less processes are requested.
 * \throws bsp_abort When bsp_begin() is called for the second time
 * 
 * \warning Undefined behaviour results when bsp_begin() is not the first
 *          statement. Note that variable or function declarations are not
 *          statements.
*/
DLL_PUBLIC void bsp_begin(bsp_pid_t maxprocs);

/** Ends an SPMD section. This must be called as the last statement of the
 * function that calls bsp_begin(). There can only be one SPMD section in a
 * program and, hence, only one instance of a bsp_begin()/bsp_end() pair is
 * allowed.
 *
 * \throws bsp_abort When bsp_end() is called without a preceding bsp_begin()
 * \throws bsp_abort When bsp_end() is not called before the end of the program.
 * 
 * \warning Undefined behaviour results when bsp_begin() is not the last 
 *          statement. Note that variable or function declarations are not
 *          statements.
 * \warning Undefined behaviour results when bsp_end() is called in another
 *          function than where bsp_begin() was called.
*/
DLL_PUBLIC void bsp_end(void);

/** Declares that \a spmd_part will be the function containing the SPMD section. 
    Whenever a sequential section has to precede the SPMD section, this function
    must be called as the first statement to main(). 

    \param spmd_part The pointer to the function containing the SPMD function.
    \param argc      The argc parameter as was passed to main()
    \param argv      The argv parameter as was passed to main()

    \throws bsp_abort When called for the second time
    \throws bsp_abort When called after or from within an SPMD section

    \warning Undefined behaviour results when bsp_init() or bsp_begin() is not
    the first statement of main()
  */
DLL_PUBLIC void bsp_init( void (*spmd_part)(void), int argc, char *argv[]) ;

/** Terminates the program abnormally from any place in the program by 
  any process, while mentioning the message as produced by the \c printf style
  format. This does not require a call to bsp_sync(). If more than one process
  calls bsp_abort() then either one, all, or a subset of processes may print
  their error message.

  \param format A printf style format string. 
  \param ...    The values as referenced by the format string

  \returns never
*/
DLL_PUBLIC void bsp_abort( const char * format, ... )
#ifndef DOXYGEN
   BSPONMPI_PRINTF_FORMAT_ATTRIBUTE(1,2)
#endif
   ; 


/** Terminates the program abnormally from any place in the program by 
  any process, while mentioning the message as produced by the \c printf style
  format. This does not require a call to bsp_sync(). If more than one process
  calls bsp_abort() then either one, all, or a subset of processes may print
  their error message.

  \param format A printf style format string. 
  \param ap The \c va_list list of arguments to the format string.

  \returns never
*/
DLL_PUBLIC void bsp_abort_va( const char * format, va_list ap );

/** Returns the number of processes available to the SPMD section. Only
  * when bsp_nprocs() is called from within the SPMD section, the returned
  * value is the actual number of parallel processes.
  *
  * \returns The number of processes. This value will be strictly greater than 0.
  */
DLL_PUBLIC bsp_pid_t bsp_nprocs(void);

/** Returns the process ID \f$ s \f$ of a particular process in the SPMD
 * section. If \f$ p \f$ is the number of processes, then \f$ 0 \leq s < p\f$.
 * This value is unique for this process in that SPMD section.
 *
 * \returns The number of processes. This value will be strictly greater than 0.
 *
 * \throws bsp_abort When called outside SPMD section.
 */
DLL_PUBLIC bsp_pid_t bsp_pid(void);

/** Returns the number of seconds elapsed since bsp_begin(). The clocks will
  * be monotonic for each process individually, but are not guaranteed to be
  * monotonic collectively, because the clocks are not required to be
  * synchronous. 
  * 
  * \warning Because the clocks are asynchronous, the value will depend on the
  * process ID. Such values should be handled carefully in branch statement
  * conditions.

  * \returns The number of elapsed seconds since bsp_begin().
  *
  * \throws bsp_abort When called outside SPMD section
 */
DLL_PUBLIC double bsp_time(void);


/** Ends the computation phase and completes the superstep by performing the
 * communication phase as prescribed by preceding calls to bsp_put(),
 * bsp_hpput(), bsp_get(), bsp_hpget(), and bsp_send().  A side-effect is that
 * all processes synchronise. After the call, a new superstep begins.
 *
 * \throws bsp_abort 
 */
DLL_PUBLIC void bsp_sync(void);


/** 
  * @}
  */


/** \defgroup BSP_DRMA BSPlib Direct Remote Memory Access methods
 *
 * Direct Remote Memory Access (DRMA) or RDMA allows processes to write
 * directly into each other's memory with a minimum of cooperation from 
 * the remote party.
 *
 * Example (from section 3.2.2 in \ref BSPlib "[1]" )
 *
 * \code
   #include <stdio.h>
   #include <bsp.h>
  
   int reverse( int x )
   {
       bsp_push_reg( &x, sizeof(x) );
       bsp_sync();
  
       bsp_put( bsp_nprocs() - bsp_pid() - 1, &x, &x, 0, sizeof(x) );
       bsp_pop_reg(&x);
       bsp_sync();
  
       return x;
   }
  
   int main( int argc, char ** argv )
   {
       int x;
       (void) argc; (void) argv;
       bsp_begin( bsp_nprocs() );
       
       x = reverse( bsp_pid() );
       printf("Process %d/%d has %d after calling reverse\n",
            bsp_pid(), bsp_nprocs(), x );
       
       bsp_end();
       return 0;
   }
   \endcode
 *
 * @{
 *
 */

/** By calling this function collectively from all processes, it creates an 
 * association between the referred memory blocks as if it were 
 * a disributed array, which is distributed over all processes and locally
 * identified by the pointer to the locally held block. 
 *
 * All processes must call this function collectively and in the same order
 * w.r.t other calls to bsp_push_reg() and bsp_pop_reg(), while supplying the
 * base address and size of the locally held part. The registration takes
 * effect after the next call to bsp_sync().
 * 
 * If this process does not have part in this association, i.e it does not have
 * data to share nor does not need data to read, it may register \c NULL as a
 * zero size memory block. If, on the contrary, this process does require to
 * access remote pieces of this array, it must offer at least a unique address
 * even though the local memory block size is zero. Such unique addresses may
 * be obtained, for example, by declaring a automatic variable on the stack.
 *
 * \param addr The base address of the local memory block
 * \param size The size of the local memory block in bytes
 *
 * \throws bsp_abort When called outside SPMD section
 * \throws bsp_abort When NULL is registered with non-zero size
 * \throws bsp_abort When a negative size is given.
 *
 * \warning Registering memory that is not owned by this process will in
 * general not lead to memory faults (usually called segmentation fault or
 * general exception errror, depending on the OS) right away. These faults will
 * usually only surface during a bsp_sync.
 *
 */
DLL_PUBLIC void bsp_push_reg( const void * addr, bsp_size_t size );

/** Removes the most recent memory association that is identified by
  * local memory block based at \a addr.  All processes must call this function
  * collectively and in the same order w.r.t to other calls to bsp_pop_reg()
  * and bsp_push_reg, with the obvious requirement that the referenced
  * memory block is the same one on all processes. 
  * 
  * \param addr The base of the local memory block
  *
  * \throws bsp_abort When call outside SPMD section
  * \throws bsp_Abort When no association exists of a memory block based at \a addr
  *
  * \warning Non-respect of collectivity and call order is not always
  * detectable by the run-time, because logically different memory areas may
  * have have the same memory address, due to memory reuse by malloc() or by
  * the C runtime.  When an error is detected, the coinciding error can only
  * be emitted at the next bsp_sync().
  * 
  */
DLL_PUBLIC void bsp_pop_reg( const void * addr );

/** Copies locally held data to the memory on another process. The destination
 * memory must have been registered through bsp_push_reg() before.  The memory
 * block with the source data can be reused right after the call has finished,
 * because the data is buffered. The data will be written at the remote process
 * after all data has been read and written by bsp_get() calls, at the end of
 * the next bsp_sync(). If the destinations of multiple bsp_put() requests,
 * which can be issued from multiple processes, overlap, then it will appear as
 * if the data has been written in some arbitrary, sequential order. 
 * 
 * \param pid The destination process
 * \param src The pointer to the locally held data
 * \param dst The local base pointer associated to the remote memory. This must
 * be the same pointer as was used to in the bsp_push_reg() call.
 * \param offset The offset in number of bytes w.r.t to the base pointer at 
 *               the remotely registered memory block on process \a pid where
 *               the data should be written
 * \param nbytes Number of bytes of data that should be copied. The call
 *               has no effect when this is zero.
 *
 * \throws bsp_abort When called outside SPMD section, unless \a nbytes is zero.
 * \throws bsp_abort When \a pid is not valid, unless \a nbytes is zero
 * \throws bsp_abort When \a offset is negative, unless \a nbytes is zero
 * \throws bsp_abort When \a dst is \c NULL, unless \a nbytes is zero
 * \throws bsp_abort When \a dst has not been registered before, unless 
 *                   \a nbytes is zero.
 * \throws bsp_abort When \a dst refers to memory that has been registered with
 *                   \c NULL on the remote process \a pid, unless \a nbytes is zero.
 * \throws bsp_abort When the remote write will be beyond the registered bounds,
 *                   unless \a nbytes is zero.
 */
DLL_PUBLIC void bsp_put( bsp_pid_t pid, const void * src, void * dst,
        bsp_size_t offset, bsp_size_t nbytes );

/** Copies locally held data to the memory on another process. The destination
 * memory must have been registered through bsp_push_reg() before. The memory
 * will not be buffered at the local process nor at the remote process. This
 * also means that neither the local memory nor the remote memory may be accessed
 * by any local process nor any other BSPlib primitive call, local or remote, 
 * until the next call to bsp_sync() . Ignoring that advice, may lead to memory
 * corruption inside the referenced memory areas.
 *
 * \note The \c hp prefix to put may lead to the expectation that there is also
 * a speed improvement. That is not necessarily true. The BSP cost model does
 * not assign different costs to buffered and unbuffered communication. Also, from
 * the perspective of using MPI as communications library, it is complicated to
 * make unbuffered communication fast. Therefore, the \c hp primitives should only
 * be used if buffering would require too much memory.
 * 
 * \param pid The destination process
 * \param src The pointer to the locally held data
 * \param dst The local base pointer associated to the remote memory. This must
 * be the same pointer as was used to in the bsp_push_reg() call.
 * \param offset The offset in number of bytes w.r.t to the base pointer at 
 *               the remotely registered memory block on process \a pid where
 *               the data should be written
 * \param nbytes Number of bytes of data that should be copied. The call
 *               has no effect when this is zero.
 *
 * \throws bsp_abort When called outside SPMD section, unless \a nbytes is zero.
 * \throws bsp_abort When \a pid is not valid, unless \a nbytes is zero
 * \throws bsp_abort When \a offset is negative, unless \a nbytes is zero
 * \throws bsp_abort When \a dst is \c NULL, unless \a nbytes is zero
 * \throws bsp_abort When \a dst has not been registered before, unless 
 *                   \a nbytes is zero.
 * \throws bsp_abort When \a dst refers to memory that has been registered with
 *                   \c NULL on the remote process \a pid, unless \a nbytes is zero.
 * \throws bsp_abort When the remote write will be beyond the registered bounds,
 *                   unless \a nbytes is zero.
 */
DLL_PUBLIC void bsp_hpput( bsp_pid_t pid, const void * src, void * dst,
        bsp_size_t offset, bsp_size_t nbytes );

/** Copies data held on a remote process to the local process. The source 
 * memory must have been registered through bsp_push_reg() before. The data
 * will be read at the start of the next bsp_sync() and written before any
 * data is written by bsp_put(). If the destination memory blocks overlap, 
 * it will appear that they have been written in arbitrary sequential order.
 * 
 * \param pid The source process
 * \param src The local base pointer associated to the remote memory area.
 *            This must the same pointer as was used to in the bsp_push_reg() call.
 * \param offset The offset in number of bytes w.r.t to the base pointer at 
 *               the remotely registered memory block on process \a pid where
 *               the data should be read from 
 * \param dst The local address where the retrieved data should be written.
 * \param nbytes Number of bytes of data that should be copied. The call
 *               has no effect when this is zero.
 *
 * \throws bsp_abort When called outside SPMD section, unless \a nbytes is zero.
 * \throws bsp_abort When \a pid is not valid, unless \a nbytes is zero
 * \throws bsp_abort When \a offset is negative, unless \a nbytes is zero
 * \throws bsp_abort When \a src is \c NULL, unless \a nbytes is zero
 * \throws bsp_abort When \a src has not been registered before, unless 
 *                   \a nbytes is zero.
 * \throws bsp_abort When \a src refers to memory that has been registered with
 *                   \c NULL on the remote process \a pid, unless \a nbytes is zero.
 * \throws bsp_abort When the remote read will be beyond the registered bounds,
 *                   unless \a nbytes is zero.
 */
DLL_PUBLIC void bsp_get( bsp_pid_t pid, const void * src, bsp_size_t offset,
        void * dst, bsp_size_t nbytes );

/** Copies data held on a remote process to the local process. The source 
 * memory must have been registered through bsp_push_reg() before. The memory
 * will not be buffered at the local process nor at the remote process. This
 * also means that neither the local memory nor the remote memory may be accessed
 * by the local process nor any other BSPlib primitive call, local or remote, 
 * until the next call to bsp_sync(). Ignoring that advice, may lead to memory
 * corruption inside the referenced memory areas.
 *
 * \note The \c hp prefix to put may lead to the expectation that there is also
 * a speed improvement. That is not necessarily true. The BSP cost model does
 * not assign different costs to buffered and unbuffered communication. Also, from
 * the perspective of using MPI as communications library, it is complicated to
 * make unbuffered communication fast. Therefore, the \c hp primitives should only
 * be used if buffering would require too much memory.
 * 
 * \param pid The source process
 * \param src The local base pointer associated to the remote memory area.
 *            This must the same pointer as was used to in the bsp_push_reg() call.
 * \param offset The offset in number of bytes w.r.t to the base pointer at 
 *               the remotely registered memory block on process \a pid where
 *               the data should be read from 
 * \param dst The local address where the retrieved data should be written.
 * \param nbytes Number of bytes of data that should be copied. The call
 *               has no effect when this is zero.
 *
 * \throws bsp_abort When called outside SPMD section, unless \a nbytes is zero.
 * \throws bsp_abort When \a pid is not valid, unless \a nbytes is zero
 * \throws bsp_abort When \a offset is negative, unless \a nbytes is zero
 * \throws bsp_abort When \a src is \c NULL, unless \a nbytes is zero
 * \throws bsp_abort When \a src has not been registered before, unless 
 *                   \a nbytes is zero.
 * \throws bsp_abort When \a src refers to memory that has been registered with
 *                   \c NULL on the remote process \a pid, unless \a nbytes is zero.
 * \throws bsp_abort When the remote read will be beyond the registered bounds,
 *                   unless \a nbytes is zero.
 */
DLL_PUBLIC void bsp_hpget( bsp_pid_t pid, const void * src, bsp_size_t offset,
        void * dst, bsp_size_t nbytes );


/** 
 * @}
 */

/** \defgroup BSP_BSMP BSPlib Bulk Synchronous Message Passing methods
  * Bulk Synchronous Message Passing allows processes to send each other
  * messages. It differs from message-passing like in e.g. MPI in that there is
  * still only bsp_sync() to synchronize. Point-to-point synchronisations are
  * not possible, since use of those could possibly lead to unacceptable
  * slow-down.
  *
  * Example (from section 4.5.2 in \ref BSPlib "[1]")
  * \code
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
          assert(status == sizeof(float));
          bsp_move( &sparse[i], sizeof(float) );
        }
      }
      bsp_set_tagsize(&tag_size);
      *sparse_out = sparse;
      *sparse_ivec_out = sparse_ivec;
      return nonzeros;
    }
    \endcode
  *
  * @{
  */

/** Set the tag size of messages for the next superstep and queries the previously
 * set value. The message tag is the fixed size part of a message. By default
 * this value is zero. The function must be collectively called by all processes
 * with the same value. The value returned through \a tag_nbytes is the 
 * value that was given in the previous call to bsp_set_tagsize() or \c 0 if 
 * there were no earlier calls.
 *
 * \note The definitions in the original BSPlib paper and the man page in the
 * Oxford BSP Toolset are rather vague. They both declare that the tag size
 * value becomes valid after the next bsp_sync(). This implies that it could
 * have been invalid before the bsp_sync() and that during that time it was not
 * allowed to send and receive messages. The Oxford BSP Toolset v1.4 implementation
 * implements slightly inconsistent behaviour, because it allows bsp_send() calls
 * immediately following the bsp_set_tagsize() while it does not define whether
 * or how to receive messages after the tag size has changed.
 *
 * \note A much preciser description would be to model this behaviour with a
 * 3-tuple \f$(A, B, C)\f$, where \f$A\f$ is the tag size of the messages
 * in the current receive queue, \f$B\f$ is the tag size of the messages that
 * are now being sent, \f$C\f$ the tag size that was given in the last call
 * to bsp_set_tagsize(). Directly after the bsp_begin() the state is
 * \f$(0, 0, 0)\f$. A call to \c bsp_set_tagsize(X) changes the state 
 * \f$(A, B, C) \rightarrow (A, B, X)\f$ and returns C. A bsp_sync() changes
 * the state \f$(A, B, C) \rightarrow (B, C, C)\f$.
 *
 * \param tag_nbytes On entry the value is used to set the tag size for after
 * the next bsp_sync(). On exit the value is changed to the value that was
 * given in the previous call to bsp_set_tagsize().
 *
 * \throws bsp_abort When called outside SPMD section
 * \throws bsp_abort When NULL is given
 * \throws bsp_abort When a negative size is given.
 */
DLL_PUBLIC void bsp_set_tagsize( bsp_size_t *tag_nbytes );

/** Send a message to another process. The tag and payload are copied right
 * away, so that the memory pointed to by \a tag and \a payload can be reused
 * right after the call. The message order is not guaranteed to be maintained
 * at reception, not even for messages between the same process pair.
 *
 * \note A message with a zero size tag and no payload will be sent.
 *  
 * \param pid The destination process ID
 * \param tag A pointer to the fixed size part of the message. This may be \c NULL
 *            if active tag size zero as well.
 * \param payload A pointer to the variable size part of the message. This may be 
              \c NULL if the \a payload_nbytes is zero.
 * \param payload_nbytes The size of the variable size part.
 *
 * \throws bsp_abort When called outside SPMD section
 * \throws bsp_abort When \a pid is not valid.
 * \throws bsp_abort When a negative size is given.
 * \throws bsp_abort When \a payload_nbytes equals \c bsp_size_t(-1), which may
 *                   happen if #bsp_size_t is an unsigned data type. 
 *
 */
DLL_PUBLIC void bsp_send( bsp_pid_t pid, const void * tag, const void * payload, 
        bsp_size_t payload_nbytes);

/** Returns the number of remaining messages in the message reception queue and
 * their total size of their payloads. 
 *
 * \param nmessages On exit will be set to the number of remaining messages in
 *                  recepetion queue. 
 * \param accum_nbytes On exit will be set to the total payload size of all
 *                  remaining messages in the reception queue. 
 *
 * \throws bsp_abort When called outside SPMD section
 * \throws bsp_abort When \a nmessages or \a accum_nbytes is set to \c NULL.
 * \throws bsp_abort When the number of messages exceeds capacity of #bsp_size_t
 * \throws bsp_abort When the total payload size exceeds capacity of #bsp_size_t
 */
DLL_PUBLIC void bsp_qsize( bsp_size_t * nmessages, bsp_size_t * accum_nbytes );

/** Copies the tag of the first message in the reception queue and sets
  * \a status to the length of the payload. If no messages is available, 
  * the tag will not be copied and \a status will be set to -1.
  *
  * \param status On exit will be set to the size of the payload of the first
  * message in the reception  queue. If the queue is empty, will be set to -1.
  * \param tag On exit the memory this pointer refers to will be overwritten
  *            by the message tag. 
  *
  * \throws bsp_abort When called outside SPMD section
  * \throws bsp_abort When \a tag is \c NULL and the tag size was not zero.
  * \throws bsp_abort When \a status is \c NULL.
 */
DLL_PUBLIC void bsp_get_tag( bsp_size_t * status, void * tag );

/** Copies the first \a reception_nbytes bytes from the payload of the first
 * message in the queue to \a payload and removes it from the queue.
 *
 * \param payload the location to copy the payload to.
 * \param reception_nbytes the number of bytes to copy from the payload
 * 
 * \throws bsp_abort When called outside SPMD section
 * \throws bsp_abort When \a payload is \c NULL and \a reception_nbytes was not zero.
 * \throws bsp_abort When \a reception_nbytes was negative
 * \throws bsp_abort When the receive queue was empty.
 */
DLL_PUBLIC void bsp_move( void * payload, bsp_size_t reception_nbytes );

/** Returns the length of the payload of the first message in the queue,
 * sets a pointer to the tag and to the payload as held in the buffer, and
 * removes the message. The resulting pointers will be sufficiently aligned to
 * hold any data type. If the queue was empty, -1 is returned.
 *
 * \param tag_ptr On exit this is set to the location where the tag is stored
 *                in the receive buffer. The location is sufficiently aligned
 *                to reference any data type.
 * \param payload_ptr On exit this is set to the location where the payload is
 *                stored in the receive buffer. The location is sufficently aligned
 *                to reference any data type.
 * \returns The size of the payload
 *
 * \throws bsp_abort When called outside SPMD section
 * \throws bsp_abort When \a tag_ptr or \a payload_ptr are NULL. 
*/
DLL_PUBLIC bsp_size_t bsp_hpmove( void ** tag_ptr, void ** payload_ptr );

/**
 * @}
 */


/** \defgroup MCBSP MulticoreBSP for C compatibility
 * This module defines symbols for MulticoreBSP for C compatibility.
 * Use the --mcbsp argument to the bspcc frontend to use these
 * signatures instead or define BSPONMPI_MCBSP_COMPAT at compile time.
 * This enables macros to change the bsp_* into mcbsp_* calls. 
 *
 * \note Only the functions with different signatures are listed here.
 * 
 * @{
 */

/** Type to store a process ID in MulticoreBSP for C.
 *  @see bsp_pid_t
 */
typedef unsigned int mcbsp_pid_t;

/** Type to store the number of processes in MulticoreBSP for C.
 *  @see bsp_pid_t
 */
typedef unsigned int mcbsp_nprocs_t;


/** Type to store the communication volume in MulticoreBSP for C.
 *  @see bsp_size_t
 */
typedef size_t       mcbsp_size_t;

/** @see bsp_begin */
DLL_PUBLIC void mcbsp_begin( mcbsp_pid_t P );

/** @see bsp_push_reg */
DLL_PUBLIC void mcbsp_push_reg( void * address, mcbsp_size_t size );

/** @see bsp_pop_reg */
DLL_PUBLIC void mcbsp_pop_reg( void * address );

/** @see bsp_put */
DLL_PUBLIC void mcbsp_put( mcbsp_pid_t pid, const void * src,
             const void * dst, mcbsp_size_t offset, mcbsp_size_t size );

/** @see bsp_get */
DLL_PUBLIC void mcbsp_get( mcbsp_pid_t pid, const void * src,
    mcbsp_size_t offset, const void * dst, mcbsp_size_t size );

/** @see bsp_hpput */
DLL_PUBLIC void mcbsp_hpput( mcbsp_pid_t pid, const void * src,
             const void * dst, mcbsp_size_t offset, mcbsp_size_t size );

/** @see bsp_hpget */
DLL_PUBLIC void mcbsp_hpget( mcbsp_pid_t pid, const void * src, 
    mcbsp_size_t offset, const void * dst, mcbsp_size_t size );

/** @see bsp_set_tagsize */
DLL_PUBLIC void mcbsp_set_tagsize( mcbsp_size_t * size );

/** @see bsp_send */
DLL_PUBLIC void mcbsp_send( mcbsp_pid_t pid, const void * tag, 
             const void * payload, mcbsp_size_t size );

/** @see bsp_qsize */
DLL_PUBLIC void mcbsp_qsize( mcbsp_nprocs_t * packets, 
                             mcbsp_size_t * total_size );

/** @see bsp_get_tag */
DLL_PUBLIC void mcbsp_get_tag( mcbsp_size_t * status, void * tag );

/** @see bsp_move */
DLL_PUBLIC void mcbsp_move( void * payload, mcbsp_size_t size );
 

/**
 *  @}
 *  */

#ifdef BSPONMPI_MCBSP_COMPAT

#define bsp_pid_t    mcbsp_pid_t
#define bsp_nprocs_t mcbsp_nprocs_t
#define bsp_size_t   mcbsp_size_t

#define bsp_direct_get  \
   #error "MulticoreBSP for C function bsp_direct_get is not supported"

#define bsp_begin         mcbsp_begin
#define bsp_push_reg      mcbsp_push_reg
#define bsp_pop_reg       mcbsp_pop_reg
#define bsp_put           mcbsp_put
#define bsp_get           mcbsp_get
#define bsp_hpput         mcbsp_hpput
#define bsp_hpget         mcbsp_hpget
#define bsp_set_tagsize   mcbsp_set_tagsize
#define bsp_send          mcbsp_send
#define bsp_hpsend        mcbsp_send
#define bsp_qsize         mcbsp_qsize
#define bsp_get_tag       mcbsp_get_tag
#define bsp_move          mcbsp_move

#endif

 
#ifdef __cplusplus
}
#endif

#endif
