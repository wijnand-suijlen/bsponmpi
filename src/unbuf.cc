#include "unbuf.h"

namespace bsplib {

Unbuf::Unbuf( size_t max_msg_size, MPI_Comm comm)
   : m_comm(MPI_COMM_NULL)
   , m_max_msg_size( max_msg_size )
{
    MPI_Comm_dup( comm, &m_comm );
}

Unbuf::~Unbuf() {
    MPI_Comm_free( &m_comm );
    m_comm = MPI_COMM_NULL;
}


void Unbuf :: send( int dst_pid, const void * addr, size_t size )
{
    Entry entry = { dst_pid, 
                    static_cast<char *>(const_cast<void*>(addr)),
                    size };
    m_sends.push_back( entry );
}

void Unbuf :: recv( int src_pid, void * addr, size_t size )
{
    Entry entry = { src_pid, 
                    static_cast<char *>(addr),
                    size };
    m_recvs.push_back( entry );
}

void Unbuf :: start( )
{
    const int tag = 0;
   m_reqs.resize( m_sends.size() + m_recvs.size() );
   size_t j = 0;
   for ( size_t i = 0 ; i < m_recvs.size(); ++j, ++i ) {
       Entry u = m_recvs[i];
       size_t size = std::min( m_max_msg_size, u.size );
       char * addr = const_cast<char *>( u.addr );
       char * next_addr = addr + size;
       u.addr = next_addr;
       MPI_Irecv( addr, u.size, MPI_BYTE, u.pid, tag, m_comm, &m_reqs[j]);
   }

   MPI_Barrier( m_comm ); // so we can use ready sends

   for ( size_t i = 0 ; i < m_sends.size(); ++j, ++i ) {
       Entry u = m_sends[i];
       size_t size = std::min( m_max_msg_size, u.size );
       char * addr = const_cast<char *>( u.addr );
       char * next_addr = addr + size;
       u.addr = next_addr;
       MPI_Irsend( addr, u.size, MPI_BYTE, u.pid, tag, m_comm, &m_reqs[j]);
   }

   m_ready.resize( m_reqs.size() );
}

void Unbuf :: wait( ) 
{
    while (true)
    {
        int outcount = MPI_UNDEFINED;
        MPI_Waitsome( m_reqs.size(), m_reqs.data(), &outcount, 
                m_ready.data(), MPI_STATUSES_IGNORE );

        if (outcount == MPI_UNDEFINED )
            break;

        for (int i = 0; i < outcount; ++i ) {
            size_t r = m_ready[i];
            const int tag = 0;

            if (r < m_recvs.size()) {
               Entry u = m_recvs[r];
               size_t size = std::min( m_max_msg_size, u.size );
               char * addr = const_cast<char *>( u.addr );
               char * next_addr = addr + size;
               u.addr = next_addr;
               
               if (size > 0)
                MPI_Irecv( addr, u.size, MPI_BYTE, u.pid, tag, m_comm, &m_reqs[r]);
            }
            else
            {
               Entry u = m_sends[r - m_recvs.size()];
               size_t size = std::min( m_max_msg_size, u.size );
               char * addr = const_cast<char *>( u.addr );
               char * next_addr = addr  + size;
               u.addr = next_addr;

               if (size > 0)
                MPI_Isend(addr, u.size, MPI_BYTE, u.pid, tag, m_comm, &m_reqs[r]);
            }
        }
    }

    m_sends.clear();
    m_recvs.clear();
    m_reqs.clear();
    m_ready.clear();
}


}
