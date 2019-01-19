#include "bsmp.h"
#include "exception.h"

namespace bsplib {
Bsmp::Bsmp(MPI_Comm comm, size_t max_msg_size)
    : m_recv_tag_size(0)
    , m_send_tag_size(0)
    , m_next_tag_size(0)
    , m_total_n_messages(0)
    , m_total_payload(0)
    , m_a2a(comm, max_msg_size)
    , m_tag_buffer()
    , m_payloads()
    , m_payload_buffer()
{

}

void Bsmp::sync() 
{
    m_payloads.clear();
    m_tag_buffer.clear();
    m_payload_buffer.clear();

    const unsigned marker = 0;
    for (int p = 0; p < m_a2a.nprocs(); ++p ) {
        serial( m_a2a, p, marker );
        serial( m_a2a, p, m_next_tag_size );
    }

    m_a2a.exchange( );

    m_total_n_messages = 0; 
    m_total_payload = 0;
    for (int p = 0 ; p < m_a2a.nprocs(); ++p )
    {
        size_t size;
        while ( deserial( m_a2a, p, size ), size != marker ) {
            size_t payload_size = size - 1;

            void * tag = m_tag_buffer.append( m_send_tag_size );
            m_payloads.push_back(
                    std::make_pair( m_payload_buffer.size(), payload_size ) 
                    );
            void * payload = m_payload_buffer.append( payload_size );
          
            std::memcpy( tag, m_a2a.recv_top(p), m_send_tag_size );
            m_a2a.recv_pop(p, m_send_tag_size );
            std::memcpy( payload, m_a2a.recv_top(p), payload_size );
            m_a2a.recv_pop(p, payload_size );

            m_total_n_messages += 1;
            m_total_payload += payload_size;
        }

        size_t tag_size = size;
        if (tag_size != m_next_tag_size)
            throw exception("bsp_set_tagsize") 
                << ": Not all processes anncounced the same tag size\n";
    }
    
    // update tag size
    m_recv_tag_size = m_send_tag_size;
    m_send_tag_size = m_next_tag_size;
}



}
