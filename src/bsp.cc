#include "bsp.h"
#include "spmd.h"
#include "rdma.h"
#include "bsmp.h"
#include "exception.h"
#include "config.h"

#include <cstdarg>
#include <cstdio>
#include <cstdlib>
#include <cassert>
#include <limits>
#include <sstream>

#include <mpi.h>

#ifdef PROFILE
#include "tictoc.h"
using bsplib::TicToc;
#endif

static bsplib::Spmd * s_spmd = NULL;
static bsplib::Rdma * s_rdma = NULL;
static bsplib::Bsmp * s_bsmp = NULL;

static const char * s_expect_abort_msg = NULL;
static bool s_aborting = false ;

extern "C" {
DLL_PUBLIC void bsp_intern_expect_success(void)
{
    /* empty */
}

DLL_PUBLIC void bsp_intern_expect_abort( const char * message )
{
    s_expect_abort_msg = message;
}
}

static void bsp_check_exit(void) 
{
    if ( s_expect_abort_msg ) {
        std::fprintf(stderr, "An error should have been detected '%s'\n",
                s_expect_abort_msg );
        std::abort();
    }

    int init = 0, fini = 0;
    MPI_Initialized(&init);
    MPI_Finalized(&fini);
    if (init && !fini && !s_aborting ) {
        // calling mpi_finalize() in atexit handler, might
        // deadlock some MPI implementations (Microsoft, IBM)
        if (s_spmd) {
            std::fprintf( stderr, "bsp_end: was not called before end of program\n");
            MPI_Abort( MPI_COMM_WORLD, EXIT_FAILURE );
        }
        MPI_Abort( MPI_COMM_WORLD, EXIT_SUCCESS );
    }
}

void bsp_abort( const char * format, ... )
{
    va_list ap;
    va_start( ap, format );
    bsp_abort_va( format, ap );
    va_end( ap );
}

void bsp_abort_va( const char * format, va_list ap )
{
    int mpi_init = 0;
    MPI_Initialized( & mpi_init );

    s_aborting = true;

    if (s_expect_abort_msg) {
        char buffer[200];
        std::vsnprintf(buffer, sizeof(buffer), format, ap );
        int len = strlen( s_expect_abort_msg );
        if (len > int(sizeof(buffer))) len = sizeof(buffer);
        // test equality of prefix
        if ( strncmp( buffer, s_expect_abort_msg, len ) == 0 ) {
            s_expect_abort_msg = NULL;
            if (mpi_init)
                MPI_Abort(MPI_COMM_WORLD, EXIT_SUCCESS );
            else
                std::exit(EXIT_SUCCESS);
        }
        std::fprintf( stderr, "Test did indeed abort, but produced wrong message.\n"
                "  Expected: '%s...'\n"
                "    Actual: '%s'\n", s_expect_abort_msg, buffer  );
    }
    std::vfprintf(stderr, format, ap );

    if (mpi_init)
        MPI_Abort(MPI_COMM_WORLD, 6 );
    else
        std::abort();
}

void bsp_init( void (*spmd_part)(void), int argc, char *argv[])
{
    int pid = -1;
    int mpi_init = 0;
    MPI_Initialized(&mpi_init);
    if (!mpi_init) {
        MPI_Init(&argc, &argv);
        atexit( bsp_check_exit );
    }

    if (s_spmd)
        bsp_abort( "bsp_init: May be called only once\n" );

    if (!spmd_part)
        bsp_abort( "bsp_init: Pointer to function must not point to NULL\n");

    MPI_Comm_rank( MPI_COMM_WORLD, &pid );
    if (pid > 0) {
        (*spmd_part)();

        if (!s_spmd)
            bsp_abort("bsp_init: spmd function did not contain SPMD section\n");

        if (!s_spmd->closed())
            bsp_abort("bsp_init: spmd function did not end with call to bsp_end\n");

        delete s_spmd;
        s_spmd = NULL;
        std::exit(EXIT_SUCCESS);
    }
}

void bsp_begin( bsp_pid_t maxprocs )
{
    int mpi_init = 0;
    MPI_Initialized(&mpi_init);
    if (!mpi_init) {
        MPI_Init(NULL, NULL);
        atexit( bsp_check_exit );
    }

    int mpi_pid = 0;
    MPI_Comm_rank( MPI_COMM_WORLD, &mpi_pid );

    if (s_spmd)
        bsp_abort( "bsp_begin: May be called only once\n" );
 
    if (mpi_pid == 0 && maxprocs <= 0)
        bsp_abort("bsp_begin: The requested number of processes must be "
                "strictly greater than 0. Instead, it was %d\n", maxprocs );

    s_spmd = new bsplib::Spmd( maxprocs );
    if ( ! s_spmd->active() ) {
        if (s_spmd->end_sync())
            bsp_abort( "bsp_end: Unexpected internal error\n");
        s_spmd->close();
        delete s_spmd;
        s_spmd = NULL;
        std::exit(EXIT_SUCCESS);
    }

    size_t max_msg_size = INT_MAX; 
    if ( bsplib::read_env( "BSPONMPI_MAX_MSG_SIZE", max_msg_size ) 
        || max_msg_size == 0 ) {
        fprintf(stderr, "bsp_begin: WARNING! BSPONMPI_MAX_MSG_SIZE "
                "was not set to a positive, non-zero integer\n");
    }

    int max_small_exchange_size = 
        std::numeric_limits<int>::max() / s_spmd->nprocs() ;
    int small_exch_size = std::min( 1024, max_small_exchange_size ) ;
    if (bsplib::read_env( "BSPONMPI_SMALL_EXCHANGE_SIZE", small_exch_size ) 
        || small_exch_size == 0 
        || small_exch_size > max_small_exchange_size ) {
        fprintf(stderr, "bsp_begin: WARNING! BSPONMPI_SMALL_EXCHANGE_SIZE"
                        " was not an integer between 1 and %d \n",
                        max_small_exchange_size );
    }
    s_rdma = new bsplib::Rdma( s_spmd->comm(), max_msg_size, 
                               small_exch_size);
    s_bsmp = new bsplib::Bsmp( s_spmd->comm(), max_msg_size,
                               small_exch_size);
}

void bsp_end() 
{
    if (!s_spmd)
        bsp_abort("bsp_end: called without calling bsp_begin() before\n");

    if (s_spmd->closed())
        bsp_abort("bsp_end: may not be called multiple times\n");

    if (s_spmd->end_sync())
        bsp_abort("bsp_sync/bsp_end: Some processes have called bsp_sync, "
               "while others have called bsp_end instead\n");

    
    delete s_rdma;
    s_rdma = NULL;
    delete s_bsmp;
    s_bsmp = NULL;

#ifdef PROFILE
    std::ostringstream out;
    bsplib::TicToc::print_stats( s_spmd->comm(), out );
    if ( s_spmd->pid() == 0 )
        printf("--- PROFILE between bsp_begin() - bsp_end() ---\n%s\n",
                out.str().c_str() );
#endif

    s_spmd->close();
    assert( s_spmd->closed() );
}

bsp_pid_t bsp_nprocs()
{
    if (!s_spmd) {
        int init, fini;
        MPI_Initialized( & init );
        MPI_Finalized( &fini );
        if ( fini )
            bsp_abort("bsp_nprocs: internal error\n");
        
        if (!init ) {
            MPI_Init(NULL, NULL);
            atexit( bsp_check_exit );
        }

        int nprocs = -1;
        MPI_Comm_size( MPI_COMM_WORLD, &nprocs);

        return nprocs;
    }
    else {
        return s_spmd->nprocs();
    }
}

bsp_pid_t bsp_pid() 
{
    if (!s_spmd && !s_spmd->closed())
        bsp_abort("bsp_pid: can only be called within SPMD section\n");

    return s_spmd->pid();
}

double bsp_time()
{
    if (!s_spmd && !s_spmd->closed())
        bsp_abort("bsp_time: can only be called within SPMD section\n");

    return s_spmd->time();
}

  
void bsp_sync()
{
    if (!s_spmd && !s_spmd->closed())
        bsp_abort("bsp_sync: can only be called within SPMD section\n");

    if (s_spmd->normal_sync())
        bsp_abort("bsp_sync/bsp_end: Some processes have called bsp_sync, "
               "while others have called bsp_end instead\n");

#ifdef PROFILE
    { TicToc imb( TicToc::IMBALANCE );
      MPI_Barrier( s_spmd->comm() );
    }
#endif

#ifdef PROFILE
    TicToc t( TicToc::SYNC );
#endif
    bool dummy_bsmp = s_bsmp->is_dummy();

    try {
        dummy_bsmp = s_rdma->sync( dummy_bsmp );
    }
    catch( bsplib :: exception & e ) {
        bsp_abort("bsp_sync/%s\n", e.str().c_str() );
    }

    try {
        s_bsmp->sync( dummy_bsmp );
    }
    catch( bsplib  :: exception & e ) {
        bsp_abort("bsp_sync/%s\n", e.str().c_str() );
    }
}

void bsp_push_reg( const void * addr, bsp_size_t size )
{
    if (!s_spmd && !s_spmd->closed())
        bsp_abort("bsp_push_reg: can only be called within SPMD section\n");

    if (size < 0)
        bsp_abort("bsp_push_reg: memory size must be positive\n");

    s_rdma->push_reg( const_cast<void*>(addr), size );
}

void bsp_pop_reg( const void * addr )
{
    if (!s_spmd && !s_spmd->closed())
        bsp_abort("bsp_pop_reg: can only be called within SPMD section\n");

    bsplib::Rdma::Memslot slot = s_rdma->lookup_reg( addr, true, false );
    if ( s_rdma->no_slot() == slot )
        bsp_abort("bsp_pop_reg: memory at address %p was not registered\n",
                addr );

    s_rdma->pop_reg( slot );
}

static bsplib::Rdma::Memslot lookup_usable_reg( const void * addr, const char * func )
{
    assert( addr != NULL );

    bsplib::Rdma::Memslot id = 
        s_rdma->lookup_reg( const_cast<void*>(addr), false, true);
    if ( id == s_rdma->no_slot() ) {
        id = s_rdma->lookup_reg( addr, true, false );
        unsigned PUSHED = bsplib::Rdma::Memblock::PUSHED;
        if ( id != s_rdma->no_slot() && 
             ( s_rdma->slot(s_spmd->pid(), id).status & PUSHED) )
          bsp_abort("%s: Remote address was just registered "
                 " with a bsp_push_reg(%p), but it hasn't become effective yet, "
                " because bsp_sync() hasn't been performed yet\n", func, addr);
        else
          bsp_abort("%s: Remote address was not registered\n", func );
    }

    assert( id != s_rdma->null_slot() );

    return id;
}


void bsp_put( bsp_pid_t pid, const void * src, void * dst,
        bsp_size_t offset, bsp_size_t nbytes )
{
#ifdef PROFILE
    TicToc t( TicToc::PUT, nbytes );
#endif

    if (nbytes == 0) // ignore any empty writes
        return;

    if (!s_spmd && !s_spmd->closed())
        bsp_abort("bsp_put: can only be called within SPMD section\n");

    if (pid < 0 || pid > s_spmd->nprocs())
        bsp_abort("bsp_put: The destination process ID does not exist\n");

    if (offset < 0)
        bsp_abort("bsp_put: The offset was negative, which is illegal\n");
    
    if (nbytes< 0)
        bsp_abort("bsp_put: The size was negative, which is illegal\n");

    if ( dst == NULL )
        bsp_abort("bsp_put: Destination address cannot be identified by NULL\n");

    bsplib::Rdma::Memslot dst_slot_id = lookup_usable_reg( dst, "bsp_put"); 

    bsplib::Rdma::Memblock dst_slot = s_rdma->slot( pid, dst_slot_id );

    if ( size_t(offset + nbytes) > dst_slot.size )
        bsp_abort("bsp_put: Writes %zu bytes beyond registered "
                  "range [%p,%p+%zu) at process %d\n",
                  offset + nbytes - dst_slot.size,
                  dst_slot.addr, dst_slot.addr, dst_slot.size, pid );

    s_rdma->put( src, pid, dst_slot_id, offset, nbytes );
}

void bsp_hpput( bsp_pid_t pid, const void * src, void * dst,
        bsp_size_t offset, bsp_size_t nbytes )
{
#ifdef PROFILE
    TicToc t( TicToc::HPPUT, nbytes );
#endif

    if (nbytes == 0) // ignore any empty writes
        return;

    if (!s_spmd && !s_spmd->closed())
        bsp_abort("bsp_put: can only be called within SPMD section\n");

    if (pid < 0 || pid > s_spmd->nprocs())
        bsp_abort("bsp_hpput: The destination process ID does not exist\n");

    if (offset < 0)
        bsp_abort("bsp_hpput: The offset was negative, which is illegal\n");
    
    if (nbytes< 0)
        bsp_abort("bsp_hpput: The size was negative, which is illegal\n");

    if ( dst == NULL )
        bsp_abort("bsp_hpput: Destination address cannot be identified by NULL\n");

    bsplib::Rdma::Memslot dst_slot_id = lookup_usable_reg( dst, "bsp_hpput" );
    bsplib::Rdma::Memblock dst_slot = s_rdma->slot( pid, dst_slot_id );

    if ( size_t(offset + nbytes) > dst_slot.size )
        bsp_abort("bsp_hpput: Writes %zu bytes beyond registered "
                  "range [%p,%p+%zu) at process %d\n",
                  offset + nbytes - dst_slot.size,
                  dst_slot.addr, dst_slot.addr, dst_slot.size, pid );

    if (pid == s_spmd->pid() )
        memcpy( dst, src, nbytes);
    else
        s_rdma->hpput( src, pid, dst_slot_id, offset, nbytes );

}

void bsp_get( bsp_pid_t pid, const void * src, bsp_size_t offset,
        void * dst, bsp_size_t nbytes )
{
#ifdef PROFILE
    TicToc t( TicToc::GET, nbytes );
#endif

    if (nbytes == 0) // ignore any empty reads
        return;

    if (!s_spmd && !s_spmd->closed())
        bsp_abort("bsp_get: can only be called within SPMD section\n");

    if (pid < 0 || pid > s_spmd->nprocs())
        bsp_abort("bsp_get: The source process ID does not exist\n");

    if (offset < 0)
        bsp_abort("bsp_get: The offset was negative, which is illegal\n");
    
    if (nbytes< 0)
        bsp_abort("bsp_get: The size was negative, which is illegal\n");

    if ( src == NULL )
        bsp_abort("bsp_get: Source address cannot be identified by NULL\n");

    bsplib::Rdma::Memslot src_slot_id = lookup_usable_reg( src, "bsp_get" );
    bsplib::Rdma::Memblock src_slot = s_rdma->slot( pid, src_slot_id );

    if ( size_t(offset + nbytes) > src_slot.size )
        bsp_abort("bsp_get: Reads %zu bytes beyond registered "
                  "range [%p,%p+%zu) at process %d\n",
                  offset + nbytes - src_slot.size,
                  src_slot.addr, src_slot.addr, src_slot.size, pid );

    s_rdma->get( pid, src_slot_id , offset, dst, nbytes );
}

void bsp_hpget( bsp_pid_t pid, const void * src, bsp_size_t offset,
        void * dst, bsp_size_t nbytes )
{
#ifdef PROFILE
    TicToc t( TicToc::HPGET, nbytes );
#endif

    if (nbytes == 0) // ignore any empty reads
        return;

    if (!s_spmd && !s_spmd->closed())
        bsp_abort("bsp_hpget: can only be called within SPMD section\n");

    if (pid < 0 || pid > s_spmd->nprocs())
        bsp_abort("bsp_hpget: The source process ID does not exist\n");

    if (offset < 0)
        bsp_abort("bsp_hpget: The offset was negative, which is illegal\n");
    
    if (nbytes< 0)
        bsp_abort("bsp_hpget: The size was negative, which is illegal\n");

    if ( src == NULL )
        bsp_abort("bsp_hpget: Source address cannot be identified by NULL\n");

    bsplib::Rdma::Memslot src_slot_id = lookup_usable_reg( src, "bsp_hpget" );
    bsplib::Rdma::Memblock src_slot = s_rdma->slot( pid, src_slot_id );

    if ( size_t(offset + nbytes) > src_slot.size )
        bsp_abort("bsp_get: Reads %zu bytes beyond registered "
                  "range [%p,%p+%zu) at process %d\n",
                  offset + nbytes - src_slot.size,
                  src_slot.addr, src_slot.addr, src_slot.size, pid );

    if (pid == s_spmd->pid() )
        memcpy( dst, src, nbytes);
    else
        s_rdma->hpget( pid, src_slot_id , offset, dst, nbytes );
}

void bsp_set_tagsize( bsp_size_t * tag_nbytes )
{
    if (!s_spmd && !s_spmd->closed())
        bsp_abort("bsp_set_tagsize: can only be called within SPMD section\n");

    if (tag_nbytes == NULL)
        bsp_abort("bsp_set_tagsize: NULL ptr as argument is now allowed\n");

    if (*tag_nbytes < 0 )
        bsp_abort("bsp_set_tagsize: Tag size may not be negative\n");

    *tag_nbytes = s_bsmp->set_tag_size( *tag_nbytes );
}

void bsp_send( bsp_pid_t pid, const void * tag, const void * payload,
        bsp_size_t payload_nbytes )
{
    if (!s_spmd && !s_spmd->closed())
        bsp_abort("bsp_send: can only be called within SPMD section\n");

    if (pid < 0 || pid > s_spmd->nprocs())
        bsp_abort("bsp_send: The source process ID does not exist\n");

    if (payload_nbytes < 0 )
        bsp_abort("bsp_send: Payload size may not be negative\n");

    if (payload_nbytes == bsp_size_t(-1))
        bsp_abort("bsp_send: Invalid payload size, becaues it is also "
                 "used as end marker\n");

#ifdef PROFILE
    TicToc t( TicToc::BSMP, s_bsmp->send_tag_size() + payload_nbytes );
#endif

    s_bsmp->send( pid, tag, payload, payload_nbytes );
}

void bsp_qsize( bsp_size_t * nmessages, bsp_size_t * accum_nbytes )
{
    if (!s_spmd && !s_spmd->closed())
        bsp_abort("bsp_qsize: can only be called within SPMD section\n");

    if ( nmessages == NULL || accum_nbytes == NULL )
        bsp_abort("bsp_qsize: both arguments may not be NULL\n");

    size_t bsp_size_max = std::numeric_limits<bsp_size_t>::max();
    if ( s_bsmp->n_total_messages() > bsp_size_max )
        bsp_abort("bsp_qsize: Integer overflow while querying number of "
                "messages in queue. There are %zu messages while the data"
                " type allows only %zu\n", s_bsmp->n_total_messages(),
                bsp_size_max );

   if ( s_bsmp->total_payload() > bsp_size_max )
        bsp_abort("bsp_qsize: Integer overflow while querying total amount of "
                "payload in queue. The total payload size is %zu  while the data"
                " type allows only %zu\n", s_bsmp->total_payload(),
                bsp_size_max );

   * nmessages = s_bsmp->n_total_messages();
   * accum_nbytes = s_bsmp->total_payload();
}

void bsp_get_tag( bsp_size_t * status, void * tag )
{
    if (!s_spmd && !s_spmd->closed())
        bsp_abort("bsp_get_tag: can only be called within SPMD section\n");
#ifdef PROFILE
    TicToc t( TicToc::BSMP, s_bsmp->recv_tag_size() );
#endif

    if ( status == NULL )
        bsp_abort("bsp_get_tag: first arguments may not be NULL\n");

    size_t tag_size = s_bsmp->recv_tag_size();
    if ( tag == NULL && tag_size  > 0 )
        bsp_abort("bsp_get_tag: tag must point to a block of available "
                "memory of at least %zu bytes long\n", tag_size );

    if ( s_bsmp->empty() ) {
        *status = bsp_size_t(-1);
    }
    else {
        *status = s_bsmp->payload_size();
        std::memcpy( tag, s_bsmp->tag(), tag_size );
    }
}

void bsp_move( void * payload, bsp_size_t reception_nbytes )
{
     if (!s_spmd && !s_spmd->closed())
        bsp_abort("bsp_move: can only be called within SPMD section\n");

     if (reception_nbytes < 0 )
         bsp_abort("bsp_move: size argument may not be negative\n");

     if ( payload == NULL && reception_nbytes > 0 )
        bsp_abort("bsp_move: payload may not be NULL if reception_nbytes is non-zero\n");

     if ( s_bsmp->empty() )
         bsp_abort("bsp_move: Message queue was empty\n");

#ifdef PROFILE
     TicToc t( TicToc::BSMP, s_bsmp->payload_size() );
#endif

     size_t size = std::min< size_t >( reception_nbytes, s_bsmp->payload_size() );
     std::memcpy( payload, s_bsmp->payload(), size );
     s_bsmp->pop();
}

bsp_size_t bsp_hpmove( void ** tag_ptr, void ** payload_ptr )
{
     if (!s_spmd && !s_spmd->closed())
        bsp_abort("bsp_hpmove: can only be called within SPMD section\n");

#ifdef PROFILE
     TicToc t( TicToc::BSMP, s_bsmp->payload_size() );
#endif

     if ( tag_ptr == NULL || payload_ptr == NULL )
        bsp_abort("bsp_hpmove: pointer arugments may not be NULL\n");

     if ( s_bsmp->empty() )
         return bsp_size_t(-1);

     *tag_ptr = s_bsmp->tag();
     *payload_ptr = s_bsmp->payload();
     size_t size = s_bsmp->payload_size();
     s_bsmp->pop();
     return size;
}

//// MulticoreBSP for C compatibility layer follows here

void mcbsp_begin( mcbsp_pid_t P )
{
   mcbsp_pid_t uint_max = std::numeric_limits<bsp_pid_t>::max();
   bsp_begin( std::min( uint_max, P ) );
}


void mcbsp_push_reg( void * address, mcbsp_size_t size )
{
    if (!s_spmd && !s_spmd->closed())
        bsp_abort("bsp_push_reg: can only be called within SPMD section\n");

    s_rdma->push_reg( address, size );
}

void mcbsp_pop_reg( void * address )
{
    bsp_pop_reg( address );
}

void mcbsp_put( mcbsp_pid_t pid, const void * src,
             const void * dst, mcbsp_size_t offset, mcbsp_size_t size )
{
#ifdef PROFILE
    TicToc t( TicToc::PUT, size );
#endif

    if (size == 0) // ignore any empty writes
        return;

    if (!s_spmd && !s_spmd->closed())
        bsp_abort("bsp_put: can only be called within SPMD section\n");

    if (pid > mcbsp_pid_t(s_spmd->nprocs()))
        bsp_abort("bsp_put: The destination process ID does not exist\n");

    if ( dst == NULL )
        bsp_abort("bsp_put: Destination address cannot be identified by NULL\n");

    bsplib::Rdma::Memslot dst_slot_id = lookup_usable_reg( dst, "bsp_put"); 

    bsplib::Rdma::Memblock dst_slot = s_rdma->slot( pid, dst_slot_id );

    if ( size_t(offset + size) > dst_slot.size )
        bsp_abort("bsp_put: Writes %zu bytes beyond registered "
                  "range [%p,%p+%zu) at process %d\n",
                  offset + size - dst_slot.size,
                  dst_slot.addr, dst_slot.addr, dst_slot.size, pid );

    s_rdma->put( src, pid, dst_slot_id, offset, size );

}
void mcbsp_hpput( mcbsp_pid_t pid, const void * src,
             const void * dst, mcbsp_size_t offset, mcbsp_size_t size )
{
#ifdef PROFILE
    TicToc t( TicToc::HPPUT, size );
#endif

    if (size == 0) // ignore any empty writes
        return;

    if (!s_spmd && !s_spmd->closed())
        bsp_abort("bsp_put: can only be called within SPMD section\n");

    if (pid > mcbsp_pid_t(s_spmd->nprocs()))
        bsp_abort("bsp_hpput: The destination process ID does not exist\n");
    if ( dst == NULL )
        bsp_abort("bsp_hpput: Destination address cannot be identified by NULL\n");

    bsplib::Rdma::Memslot dst_slot_id = lookup_usable_reg( dst, "bsp_hpput" );
    bsplib::Rdma::Memblock dst_slot = s_rdma->slot( pid, dst_slot_id );

    if ( size_t(offset + size ) > dst_slot.size )
        bsp_abort("bsp_hpput: Writes %zu bytes beyond registered "
                  "range [%p,%p+%zu) at process %d\n",
                  offset + size - dst_slot.size,
                  dst_slot.addr, dst_slot.addr, dst_slot.size, pid );

    if (pid == mcbsp_pid_t(s_spmd->pid()) )
        memcpy( const_cast<void *>(dst), src, size );
    else
        s_rdma->hpput( src, pid, dst_slot_id, offset, size );

}

void mcbsp_get( mcbsp_pid_t pid, const void * src,
    mcbsp_size_t offset, const void * dst, mcbsp_size_t size )
{
#ifdef PROFILE
    TicToc t( TicToc::GET, size );
#endif

    if (size == 0) // ignore any empty reads
        return;

    if (!s_spmd && !s_spmd->closed())
        bsp_abort("bsp_get: can only be called within SPMD section\n");

    if (pid > mcbsp_pid_t(s_spmd->nprocs()))
        bsp_abort("bsp_get: The source process ID does not exist\n");

    if ( src == NULL )
        bsp_abort("bsp_get: Source address cannot be identified by NULL\n");

    bsplib::Rdma::Memslot src_slot_id = lookup_usable_reg( src, "bsp_get" );
    bsplib::Rdma::Memblock src_slot = s_rdma->slot( pid, src_slot_id );

    if ( size_t(offset + size ) > src_slot.size )
        bsp_abort("bsp_get: Reads %zu bytes beyond registered "
                  "range [%p,%p+%zu) at process %d\n",
                  offset + size - src_slot.size,
                  src_slot.addr, src_slot.addr, src_slot.size, pid );

    s_rdma->get( pid, src_slot_id , offset, const_cast<void*>(dst), size );
}

void mcbsp_hpget( mcbsp_pid_t pid, const void * src, 
    mcbsp_size_t offset, const void * dst, mcbsp_size_t size )
{
#ifdef PROFILE
    TicToc t( TicToc::HPGET, size );
#endif

    if (size == 0) // ignore any empty reads
        return;

    if (!s_spmd && !s_spmd->closed())
        bsp_abort("bsp_hpget: can only be called within SPMD section\n");

    if (pid > mcbsp_pid_t(s_spmd->nprocs()))
        bsp_abort("bsp_hpget: The source process ID does not exist\n");

    if ( src == NULL )
        bsp_abort("bsp_hpget: Source address cannot be identified by NULL\n");

    bsplib::Rdma::Memslot src_slot_id = lookup_usable_reg( src, "bsp_hpget" );
    bsplib::Rdma::Memblock src_slot = s_rdma->slot( pid, src_slot_id );

    if ( size_t(offset + size ) > src_slot.size )
        bsp_abort("bsp_get: Reads %zu bytes beyond registered "
                  "range [%p,%p+%zu) at process %d\n",
                  offset + size - src_slot.size,
                  src_slot.addr, src_slot.addr, src_slot.size, pid );

    if (pid == mcbsp_pid_t(s_spmd->pid()) )
        memcpy( const_cast<void *>(dst), src, size );
    else
        s_rdma->hpget( pid, src_slot_id , offset, const_cast<void*>(dst), size );

}

void mcbsp_set_tagsize( mcbsp_size_t * size )
{
    if (!s_spmd && !s_spmd->closed())
        bsp_abort("bsp_set_tagsize: can only be called within SPMD section\n");

    if (size == NULL)
        bsp_abort("bsp_set_tagsize: NULL ptr as argument is now allowed\n");

    *size = s_bsmp->set_tag_size( *size );
}

void mcbsp_send( mcbsp_pid_t pid, const void * tag, 
             const void * payload, mcbsp_size_t size )
{
    if (!s_spmd && !s_spmd->closed())
        bsp_abort("bsp_send: can only be called within SPMD section\n");

    if (pid > mcbsp_pid_t(s_spmd->nprocs()))
        bsp_abort("bsp_send: The source process ID does not exist\n");

    if (size == mcbsp_size_t(-1))
        bsp_abort("bsp_send: Invalid payload size, becaues it is also "
                 "used as end marker\n");

#ifdef PROFILE
    TicToc t( TicToc::BSMP, s_bsmp->send_tag_size() + size );
#endif

    s_bsmp->send( pid, tag, payload, size );
}

void mcbsp_qsize( mcbsp_nprocs_t * packets, 
                             mcbsp_size_t * total_size )
{
    if (!s_spmd && !s_spmd->closed())
        bsp_abort("bsp_qsize: can only be called within SPMD section\n");

    if ( packets == NULL || total_size == NULL )
        bsp_abort("bsp_qsize: both arguments may not be NULL\n");

    size_t mcbsp_p_max = std::numeric_limits<mcbsp_nprocs_t>::max();
    size_t mcbsp_size_max = std::numeric_limits<mcbsp_size_t>::max();
    if ( s_bsmp->n_total_messages() > mcbsp_p_max )
        bsp_abort("bsp_qsize: Integer overflow while querying number of "
                "messages in queue. There are %zu messages while the data"
                " type allows only %zu\n", s_bsmp->n_total_messages(),
                mcbsp_p_max );

   if ( s_bsmp->total_payload() > mcbsp_size_max )
        bsp_abort("bsp_qsize: Integer overflow while querying total amount of "
                "payload in queue. The total payload size is %zu  while the data"
                " type allows only %zu\n", s_bsmp->total_payload(),
                mcbsp_size_max );

   * packets = s_bsmp->n_total_messages();
   * total_size = s_bsmp->total_payload();
}


void mcbsp_get_tag( mcbsp_size_t * status, void * tag )
{
    if (!s_spmd && !s_spmd->closed())
        bsp_abort("bsp_get_tag: can only be called within SPMD section\n");

#ifdef PROFILE
    TicToc t( TicToc::BSMP, s_bsmp->recv_tag_size() );
#endif

    if ( status == NULL )
        bsp_abort("bsp_get_tag: first arguments may not be NULL\n");

    size_t tag_size = s_bsmp->recv_tag_size();
    if ( tag == NULL && tag_size  > 0 )
        bsp_abort("bsp_get_tag: tag must point to a block of available "
                "memory of at least %zu bytes long\n", tag_size );

    if ( s_bsmp->empty() ) {
        *status = mcbsp_size_t(-1);
    }
    else {
        *status = s_bsmp->payload_size();
        std::memcpy( tag, s_bsmp->tag(), tag_size );
    }
}


void mcbsp_move( void * payload, mcbsp_size_t reception_bytes )
{
     if (!s_spmd && !s_spmd->closed())
        bsp_abort("bsp_move: can only be called within SPMD section\n");

     if ( payload == NULL && reception_bytes > 0 )
        bsp_abort("bsp_move: payload may not be NULL if size parameter is non-zero\n");

     if ( s_bsmp->empty() )
         bsp_abort("bsp_move: Message queue was empty\n");

#ifdef PROFILE
     TicToc t( TicToc::BSMP, s_bsmp->payload_size() );
#endif

     size_t size =
         std::min< size_t >( reception_bytes, s_bsmp->payload_size() );
     std::memcpy( payload, s_bsmp->payload(), size );
     s_bsmp->pop();
}
 
