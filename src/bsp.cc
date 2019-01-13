#include "bsp.h"
#include "spmd.h"
#include "rdma.h"

#include <cstdarg>
#include <cstdio>
#include <cstdlib>
#include <cassert>

#include <mpi.h>

static bsplib::Spmd * s_spmd = NULL;
static bsplib::Rdma * s_rdma = NULL;

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
    catch( bsplib :: Rdma :: exception & e ) {
        bsp_abort("bsp_sync/%s\n", e.str().c_str() );
    }
}

void bsp_put( bsp_pid_t pid, const void * src, void * dst,
        bsp_size_t offset, bsp_size_t nbytes )
{
    if (!s_spmd && !s_spmd->ended())
        bsp_abort("bsp_put: can only be called within SPMD section\n");

    if (nbytes == 0) // ignore any empty writes
        return;

    if (pid < 0 || pid > s_spmd->nprocs())
        bsp_abort("bsp_put: The destination process ID does not exist\n");

    if (offset < 0)
        bsp_abort("bsp_put: The offset was negative, which is illegal\n");
    
    if (nbytes< 0)
        bsp_abort("bsp_put: The size was negative, which is illegal\n");

    using bsplib :: Rdma;
    Rdma::MemslotID dst_slot_id = s_rdma->lookup_reg( dst );
    if ( dst_slot_id == Rdma::MemslotID(-1) )
        bsp_abort("bsp_put: Destination address %p was not registered\n",
                dst );

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
    if (!s_spmd && !s_spmd->ended())
        bsp_abort("bsp_put: can only be called within SPMD section\n");

    if (nbytes == 0) // ignore any empty writes
        return;

    if (pid < 0 || pid > s_spmd->nprocs())
        bsp_abort("bsp_hpput: The destination process ID does not exist\n");

    if (offset < 0)
        bsp_abort("bsp_hpput: The offset was negative, which is illegal\n");
    
    if (nbytes< 0)
        bsp_abort("bsp_hpput: The size was negative, which is illegal\n");


    using bsplib :: Rdma;
    Rdma::MemslotID dst_slot_id = s_rdma->lookup_reg( dst );
    if ( dst_slot_id == Rdma::MemslotID(-1) )
        bsp_abort("bsp_hpput: Destination address %p was not registered\n",
                dst );

    Rdma::Memslot dst_slot = s_rdma->slot( pid, dst_slot_id );

    if ( size_t(offset + nbytes) > dst_slot.size )
        bsp_abort("bsp_hpput: Writes %zu bytes beyond registered "
                  "range [%p,%p+%zu) at process %d\n",
                  offset + nbytes - dst_slot.size,
                  dst_slot.addr, dst_slot.addr, dst_slot.size, pid );

    if (pid == s_spmd->pid() )
        memcpy( dst, src, nbytes);
    else
        s_rdma->put( src, pid, dst_slot_id, offset, nbytes );

}

void bsp_get( bsp_pid_t pid, const void * src, bsp_size_t offset,
        void * dst, bsp_size_t nbytes )
{
    if (!s_spmd && !s_spmd->ended())
        bsp_abort("bsp_get: can only be called within SPMD section\n");

    if (nbytes == 0) // ignore any empty reads
        return;

    if (pid < 0 || pid > s_spmd->nprocs())
        bsp_abort("bsp_get: The source process ID does not exist\n");

    if (offset < 0)
        bsp_abort("bsp_get: The offset was negative, which is illegal\n");
    
    if (nbytes< 0)
        bsp_abort("bsp_get: The size was negative, which is illegal\n");


    using bsplib :: Rdma;
    Rdma::MemslotID src_slot_id = s_rdma->lookup_reg( src );
    if ( src_slot_id == Rdma::MemslotID(-1) )
        bsp_abort("bsp_get: Destination address %p was not registered\n",
                src );

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
    if (!s_spmd && !s_spmd->ended())
        bsp_abort("bsp_hpget: can only be called within SPMD section\n");

    if (nbytes == 0) // ignore any empty reads
        return;

    if (pid < 0 || pid > s_spmd->nprocs())
        bsp_abort("bsp_hpget: The source process ID does not exist\n");

    if (offset < 0)
        bsp_abort("bsp_hpget: The offset was negative, which is illegal\n");
    
    if (nbytes< 0)
        bsp_abort("bsp_hpget: The size was negative, which is illegal\n");


    using bsplib :: Rdma;
    Rdma::MemslotID src_slot_id = s_rdma->lookup_reg( src );
    if ( src_slot_id == Rdma::MemslotID(-1) )
        bsp_abort("bsp_get: Destination address %p was not registered\n",
                src );

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

