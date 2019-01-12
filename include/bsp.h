#ifndef BSPONMPI_BSP_H
#define BSPONMPI_BSP_H

#include <stdarg.h>

#ifdef __cplusplus 
extern "C" {
#endif

typedef int bsp_pid_t;

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


#ifdef __cplusplus
}
#endif

#endif
