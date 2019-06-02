#include "unbuf.h"
#include "exception.h"

#include <limits>
#include <algorithm>

namespace bsplib {

Unbuf::Unbuf( size_t max_msg_size, MPI_Comm comm)
   : m_pid(-1), m_nprocs(-1)
   , m_comm(MPI_COMM_NULL)
   , m_max_msg_size( max_msg_size )
{
    MPI_Comm_dup( comm, &m_comm );
    MPI_Comm_size( m_comm, &m_nprocs );
    MPI_Comm_rank( m_comm, &m_pid );
}

Unbuf::~Unbuf() {
    MPI_Comm_free( &m_comm );
    m_comm = MPI_COMM_NULL;
}

int Unbuf :: send( int dst_pid, const void * addr, size_t size )
{
    size_t int_max = std::numeric_limits<int>::max();
    if ( m_sends.size() > (int_max/2 - m_pid)/m_nprocs )
        throw exception("bsp_hpput") << "Too many HP messages";

    Entry entry = { dst_pid, 
                    static_cast<char *>(const_cast<void*>(addr)),
                    size,
                    (int( m_sends.size() ) * m_nprocs + m_pid)*2
                  };
    m_sends.push_back( entry );
    return entry.tag;
}


void Unbuf :: send( int recv_tag, int dst_pid, const void * addr, size_t size )
{
    Entry entry = { dst_pid, 
                    static_cast<char *>(const_cast<void*>(addr)),
                    size, 
                    recv_tag };
    m_sends.push_back( entry );
}

int Unbuf :: recv( int src_pid, void * addr, size_t size )
{
    size_t int_max = std::numeric_limits<int>::max();
    if ( m_recvs.size() > int_max )
        throw exception("bsp_hpput") << "Too many HP messages";
 
    Entry entry = { src_pid, 
                    static_cast<char *>(addr),
                    size,
                    int (2* m_recvs.size() + 1) };
    m_recvs.push_back( entry );
    return entry.tag;
}


void Unbuf :: recv( int send_tag, int src_pid, void * addr, size_t size )
{
    Entry entry = { src_pid, 
                    static_cast<char *>(addr),
                    size, 
                    send_tag };
    m_recvs.push_back( entry );
}

void Unbuf :: start( )
{
   m_reqs.resize( m_sends.size() + m_recvs.size() );
   size_t j = 0;
   for ( size_t i = 0 ; i < m_recvs.size(); ++j, ++i ) {
       Entry & u = m_recvs[i];
       size_t size = std::min( m_max_msg_size, u.size );
       char * addr = const_cast<char *>( u.addr );
       char * next_addr = addr + size;
       u.addr = next_addr;
       u.size -= size;
       const int tag = u.tag;
       MPI_Irecv( addr, int(size), MPI_BYTE, u.pid, tag, m_comm, &m_reqs[j]);
   }

   MPI_Barrier( m_comm ); // so we can use ready sends

   for ( size_t i = 0 ; i < m_sends.size(); ++j, ++i ) {
       Entry & u = m_sends[i];
       size_t size = std::min( m_max_msg_size, u.size );
       char * addr = const_cast<char *>( u.addr );
       char * next_addr = addr + size;
       u.addr = next_addr;
       u.size -= size;
       const int tag = u.tag;
       MPI_Irsend( addr, int(size), MPI_BYTE, u.pid, tag, m_comm, &m_reqs[j]);
   }

   m_ready.resize( m_reqs.size() );
}

void Unbuf :: wait( ) 
{
    while (true)
    {
        int outcount = MPI_UNDEFINED;
        MPI_Waitsome( int(m_reqs.size()), m_reqs.data(), &outcount, 
                m_ready.data(), MPI_STATUSES_IGNORE );
    
        if (outcount == MPI_UNDEFINED )
            break;

        for (int i = 0; i < outcount; ++i ) {
            size_t r = m_ready[i];

            if (r < m_recvs.size()) {
               Entry & u = m_recvs[r];
               size_t size = std::min( m_max_msg_size, u.size );
               char * addr = const_cast<char *>( u.addr );
               char * next_addr = addr + size;
               u.addr = next_addr;
               u.size -= size;
               const int tag = u.tag;
               
               if (size > 0) 
                MPI_Irecv( addr, int(size), MPI_BYTE, 
                           u.pid, tag, m_comm, &m_reqs[r]);
            }
            else
            {
               Entry & u = m_sends[r - m_recvs.size()];
               size_t size = std::min( m_max_msg_size, u.size );
               char * addr = const_cast<char *>( u.addr );
               char * next_addr = addr  + size;
               u.addr = next_addr;
               u.size -= size;
               const int tag = u.tag;

               if (size > 0)
                MPI_Isend(addr, int(size), MPI_BYTE, 
                          u.pid, tag, m_comm, &m_reqs[r]);
            }
        }
    }

    m_sends.clear();
    m_recvs.clear();
    m_reqs.clear();
    m_ready.clear();
}


}
