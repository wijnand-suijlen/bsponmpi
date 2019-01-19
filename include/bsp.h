#ifndef BSPONMPI_BSP_H
#define BSPONMPI_BSP_H

#include <stdarg.h>

#ifdef __cplusplus 
extern "C" {
#endif

typedef int bsp_pid_t;
typedef int bsp_size_t;

#ifdef _GNUC_
  #define BSPONMPI_PRINTF_FORMAT_ATTRIBUTE(format_idx, first_va)\
                __attribute__(( format(printf, format_idx, first_va) ))
#else

  #define BSPONMPI_PRINTF_FORMAT_ATTRIBUTE(format_idx, first_va) /* empty */

#endif


/** \defgroup SPMD SPMD Framework
 *
  *  @{
 */


/** Starts an SPMD section with at most \a maxprocs parallel processes. This
 * must be called as the first statement of a function, which may be main().
 * In that same function there must also be call to bsp_end() as the last
 * statement.  There may be only one instance of a bsp_begin()/bsp_end() pair
 * within a program.  If the enclosing function is not  main(), then its
 * reference must be passed as parameter to a call to bsp_init(), which in its
 * turn must be the first statement in main().
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
void bsp_begin(bsp_pid_t maxprocs);

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
void bsp_end(void);

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
void bsp_init( void (*spmd_part)(void), int argc, char *argv[]) ;

/** Terminates the program abnormally from any place in the program by 
  any process, while mentioning the message as produced by the \c printf style
  format. This does not require a call to bsp_sync(). If more than one process
  calls bsp_abort() then either one, all, or a subset of processes may print
  their error message.

  \param format A printf style format string. 
  \param ...    The values as referenced by the format string

  */
void bsp_abort( const char * format, ... )
   BSPONMPI_PRINTF_FORMAT_ATTRIBUTE(1,2); 

void bsp_abort_va( const char * format, va_list ap );

/** Returns the number of processes available to the SPMD section. Only
  * when bsp_nprocs() is called from within the SPMD section, the returned
  * value is the actual number of parallel processes.
  *
  * \returns The number of processes. This value will be strictly greater than 0.
  */
bsp_pid_t bsp_nprocs(void);

/** Returns the process ID \f$ s \f$ of a particular process in the SPMD
 * section. If \f$ p \f$ is the number of processes, then \f$ 0 \leq s < p$.
 * This value is unique for this process in that SPMD section.
 *
 * \returns The number of processes. This value will be strictly greater than 0.
 *
 * \throws bsp_abort When called outside SPMD section.
 */
bsp_pid_t bsp_pid(void);

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
double bsp_time(void);


/* Ends the computation phase and completes the superstep by performing the
 * communication phase as prescribed by preceding calls to bsp_put(),
 * bsp_hpput(), bsp_get(), bsp_hpget(), and bsp_send().  A side-effect is that
 * all processes synchronise. After the call, a new superstep begins.
 */
void bsp_sync(void);


/** 
  * @}
  */


/** \defgroup DRMA Direct Remote Memory Access
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
 * If this process does not have part in this associaten , i.e it does not have
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
void bsp_push_reg( const void * addr, bsp_size_t size );

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
void bsp_pop_reg( const void * addr );

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
void bsp_put( bsp_pid_t pid, const void * src, void * dst,
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
 */void bsp_hpput( bsp_pid_t pid, const void * src, void * dst,
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
void bsp_get( bsp_pid_t pid, const void * src, bsp_size_t offset,
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
void bsp_hpget( bsp_pid_t pid, const void * src, bsp_size_t offset,
        void * dst, bsp_size_t nbytes );


/** 
 * @}
 */

/** \defgroup BSMP Bulk Synchronous Message Passing
  * @{
  */

/** Set the tag size of messages for the next superstep and queries the previously
  set value. */
void bsp_set_tagsize( int *tag_nbytes );

/* @} */

#ifdef __cplusplus
}
#endif

#endif
