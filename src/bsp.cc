#include "bsp.h"
#include "spmd.h"
#include "rdma.h"
#include "bsmp.h"
#include "exception.h"

#include <cstdarg>
#include <cstdio>
#include <cstdlib>
#include <cassert>
#include <limits>

#include <mpi.h>

static bsplib::Spmd * s_spmd = NULL;
static bsplib::Rdma * s_rdma = NULL;
static bsplib::Bsmp * s_bsmp = NULL;

static void bsp_finalize(void)
{
    int init = 0;
    MPI_Initialized(&init);
    if (init)
        MPI_Finalize();
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
    std::vfprintf(stderr, format, ap );

    int mpi_init = 0;
    MPI_Initialized( & mpi_init );
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
        atexit( bsp_finalize );
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

        if (!s_spmd->ended())
            bsp_abort("bsp_init: spmd function did not end with call to bsp_end\n");

        delete s_spmd;
        s_spmd = NULL;
        std::exit(0);
    }
}

void bsp_begin( bsp_pid_t maxprocs )
{
    int mpi_init = 0;
    MPI_Initialized(&mpi_init);
    if (!mpi_init) {
        MPI_Init(NULL, NULL);
        atexit( bsp_finalize );
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
        delete s_spmd;
        s_spmd = NULL;
        std::exit(0);
    }
    s_rdma = new bsplib::Rdma( s_spmd->comm(), INT_MAX );
    s_bsmp = new bsplib::Bsmp( s_spmd->comm(), INT_MAX );
}

void bsp_end() 
{
    if (!s_spmd)
        bsp_abort("bsp_end: called without calling bsp_begin() before\n");

    if (s_spmd->ended())
        bsp_abort("bsp_end: may not be called multiple times\n");

    if (s_spmd->end_sync())
        bsp_abort("bsp_end: Some other processes have called bsp_sync instead\n");

    assert( s_spmd->ended() );
    
    delete s_rdma;
    s_rdma = NULL;
    delete s_bsmp;
    s_bsmp = NULL;
}

bsp_pid_t bsp_nprocs()
{
    if (!s_spmd) {
        int mpi_init = 0;
        MPI_Initialized( & mpi_init );
        if (!mpi_init) {
            MPI_Init(NULL, NULL);
            atexit( bsp_finalize );
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
    if (!s_spmd && !s_spmd->ended())
        bsp_abort("bsp_pid: can only be called within SPMD section\n");

    return s_spmd->pid();
}

void bsp_sync()
{
    if (!s_spmd && !s_spmd->ended())
        bsp_abort("bsp_sync: can only be called within SPMD section\n");

    if (s_spmd->normal_sync())
        bsp_abort("bsp_sync: some processes called bsp_end instead\n");

    try {
        s_rdma->sync( );
    }
    catch( bsplib :: exception & e ) {
        bsp_abort("bsp_sync - RDMA error by %s\n", e.str().c_str() );
    }

    try {
        s_bsmp->sync( );
    }
    catch( bsplib  :: exception & e ) {
        bsp_abort("bsp_sync - BSMP error by %s\n", e.str().c_str() );
    }
}

void bsp_push_reg( const void * addr, bsp_size_t size )
{
    if (!s_spmd && !s_spmd->ended())
        bsp_abort("bsp_push_reg: can only be called within SPMD section\n");

    if (size < 0)
        bsp_abort("bsp_push_reg: memory size must be positive\n");

    s_rdma->push_reg( const_cast<void*>(addr), size );
}

void bsp_pop_reg( const void * addr )
{
    if (!s_spmd && !s_spmd->ended())
        bsp_abort("bsp_pop_reg: can only be called within SPMD section\n");

    bsplib::Rdma::MemslotID slot = s_rdma->lookup_reg( addr );
    if ( s_rdma->no_slot() == slot )
        bsp_abort("bsp_pop_reg: memory at address %p was not registered\n",
                addr );

    s_rdma->pop_reg( slot );
}

void bsp_put( bsp_pid_t pid, const void * src, void * dst,
        bsp_size_t offset, bsp_size_t nbytes )
{
    if (nbytes == 0) // ignore any empty writes
        return;

    if (!s_spmd && !s_spmd->ended())
        bsp_abort("bsp_put: can only be called within SPMD section\n");

    if (pid < 0 || pid > s_spmd->nprocs())
        bsp_abort("bsp_put: The destination process ID does not exist\n");

    if (offset < 0)
        bsp_abort("bsp_put: The offset was negative, which is illegal\n");
    
    if (nbytes< 0)
        bsp_abort("bsp_put: The size was negative, which is illegal\n");

    if ( dst == NULL )
        bsp_abort("bsp_put: Destination address cannot be identified by NULL\n");

    using bsplib :: Rdma;
    Rdma::MemslotID dst_slot_id = s_rdma->lookup_reg( dst );
    if ( dst_slot_id == s_rdma->no_slot() )
        bsp_abort("bsp_put: Destination address %p was not registered\n",
                dst );

    if ( dst_slot_id == s_rdma->null_slot() )
        bsp_abort("bsp_put: Destination may not be used because it was registered with NULL\n");


    Rdma::Memslot dst_slot = s_rdma->slot( pid, dst_slot_id );

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
    if (nbytes == 0) // ignore any empty writes
        return;

    if (!s_spmd && !s_spmd->ended())
        bsp_abort("bsp_put: can only be called within SPMD section\n");

    if (pid < 0 || pid > s_spmd->nprocs())
        bsp_abort("bsp_hpput: The destination process ID does not exist\n");

    if (offset < 0)
        bsp_abort("bsp_hpput: The offset was negative, which is illegal\n");
    
    if (nbytes< 0)
        bsp_abort("bsp_hpput: The size was negative, which is illegal\n");

    if ( dst == NULL )
        bsp_abort("bsp_hpput: Destination address cannot be identified by NULL\n");

    using bsplib :: Rdma;
    Rdma::MemslotID dst_slot_id = s_rdma->lookup_reg( dst );
    if ( dst_slot_id == s_rdma->no_slot() )
        bsp_abort("bsp_hpput: Destination address %p was not registered\n",
                dst );

    if ( dst_slot_id == s_rdma->null_slot() )
        bsp_abort("bsp_hpput: Destination may not be used as it was registered with NULL\n");

    Rdma::Memslot dst_slot = s_rdma->slot( pid, dst_slot_id );

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
    if (nbytes == 0) // ignore any empty reads
        return;

    if (!s_spmd && !s_spmd->ended())
        bsp_abort("bsp_get: can only be called within SPMD section\n");

    if (pid < 0 || pid > s_spmd->nprocs())
        bsp_abort("bsp_get: The source process ID does not exist\n");

    if (offset < 0)
        bsp_abort("bsp_get: The offset was negative, which is illegal\n");
    
    if (nbytes< 0)
        bsp_abort("bsp_get: The size was negative, which is illegal\n");

    if ( src == NULL )
        bsp_abort("bsp_get: Source address cannot be identified by NULL\n");

    using bsplib :: Rdma;
    Rdma::MemslotID src_slot_id = s_rdma->lookup_reg( src );
    if ( src_slot_id == s_rdma->no_slot() ) {
        if ( s_rdma->was_just_pushed( src ) )
          bsp_abort("bsp_get: Destination address %p was just registered "
                 " with a bsp_push_reg(), but it hasn't become effective yet, "
                " because bsp_sync() hasn't been performed yet\n", src );
        else
          bsp_abort("bsp_get: Destination address %p was not registered\n", src );
    }

    if ( src_slot_id == s_rdma->null_slot() )
        bsp_abort("bsp_get: Source may not be used as it was registered with NULL\n");

    Rdma::Memslot src_slot = s_rdma->slot( pid, src_slot_id );


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
    if (nbytes == 0) // ignore any empty reads
        return;

    if (!s_spmd && !s_spmd->ended())
        bsp_abort("bsp_hpget: can only be called within SPMD section\n");

    if (pid < 0 || pid > s_spmd->nprocs())
        bsp_abort("bsp_hpget: The source process ID does not exist\n");

    if (offset < 0)
        bsp_abort("bsp_hpget: The offset was negative, which is illegal\n");
    
    if (nbytes< 0)
        bsp_abort("bsp_hpget: The size was negative, which is illegal\n");

    if ( src == NULL )
        bsp_abort("bsp_hpget: Source address cannot be identified by NULL\n");

    using bsplib :: Rdma;
    Rdma::MemslotID src_slot_id = s_rdma->lookup_reg( src );
    if ( src_slot_id == s_rdma->no_slot() )
        bsp_abort("bsp_get: Destination address %p was not registered\n",
                src );

    if ( src_slot_id == s_rdma->null_slot() )
        bsp_abort("bsp_hpget: Source may not be used as it was registered with NULL\n");

    Rdma::Memslot src_slot = s_rdma->slot( pid, src_slot_id );


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
    if (!s_spmd && !s_spmd->ended())
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
    if (!s_spmd && !s_spmd->ended())
        bsp_abort("bsp_send: can only be called within SPMD section\n");

    if (pid < 0 || pid > s_spmd->nprocs())
        bsp_abort("bsp_send: The source process ID does not exist\n");

    if (payload_nbytes < 0 )
        bsp_abort("bsp_send: Payload size may not be negative\n");

    if (payload_nbytes == bsp_size_t(-1))
        bsp_abort("bsp_send: Invalid payload size, becaues it is also "
                 "used as end marker\n");

    s_bsmp->send( pid, tag, payload, payload_nbytes );
}

void bsp_qsize( bsp_size_t * nmessages, bsp_size_t * accum_nbytes )
{
    if (!s_spmd && !s_spmd->ended())
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
    if (!s_spmd && !s_spmd->ended())
        bsp_abort("bsp_get_tag: can only be called within SPMD section\n");

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
     if (!s_spmd && !s_spmd->ended())
        bsp_abort("bsp_move: can only be called within SPMD section\n");

     if (reception_nbytes < 0 )
         bsp_abort("bsp_move: size argument may not be negative\n");

     if ( payload == NULL && reception_nbytes > 0 )
        bsp_abort("bsp_move: payload may not be NULL if reception_nbytes is non-zero\n");

     if ( s_bsmp->empty() )
         bsp_abort("bsp_move: Mesaage queue was empty\n");

     size_t size = std::min< size_t >( reception_nbytes, s_bsmp->payload_size() );
     std::memcpy( payload, s_bsmp->payload(), size );
     s_bsmp->pop();
}

bsp_size_t bsp_hpmove( void ** tag_ptr, void ** payload_ptr )
{
     if (!s_spmd && !s_spmd->ended())
        bsp_abort("bsp_hpmove: can only be called within SPMD section\n");

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


