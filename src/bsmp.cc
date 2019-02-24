#include "bsmp.h"
#include "exception.h"

#ifdef PROFILE
#include "tictoc.h"
#endif

namespace bsplib {
Bsmp::Bsmp(MPI_Comm comm, size_t max_msg_size, size_t small_exch_size)
    : m_set_tag_size_counter(0)
    , m_recv_tag_size(0)
    , m_send_tag_size(0)
    , m_next_tag_size(0)
    , m_total_n_messages(0)
    , m_total_payload(0)
    , m_send_empty(true)
    , m_a2a(comm, max_msg_size, small_exch_size )
    , m_tag_buffer()
    , m_payloads()
    , m_payload_buffer()
{

}

void Bsmp::sync( bool dummy ) 
{
    m_payloads.clear();
    m_tag_buffer.clear();
    m_payload_buffer.clear();

    if (!dummy) {
#ifdef PROFILE
    TicToc t( TicToc::BSMP);
#endif

    const unsigned marker = 0;
    for (int p = 0; p < m_a2a.nprocs(); ++p ) {
        serial( m_a2a, p, marker );
        serial( m_a2a, p, m_set_tag_size_counter );
        serial( m_a2a, p, m_next_tag_size );
    
#ifdef PROFILE
        t.add_bytes( m_a2a.send_size(p) );
#endif
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

        size_t counter = size_t(-1);
        deserial( m_a2a, p, counter );
        if (counter != m_set_tag_size_counter )
            throw exception("bsp_set_tagsize") 
                << "Not all processes set the tag size same number of times\n";

        size_t tag_size = size_t(-1);
        deserial( m_a2a, p, tag_size );
        if (tag_size != m_next_tag_size)
            throw exception("bsp_set_tagsize") 
                << "Not all processes announced the same tag size\n"
                << "  -> process " << m_a2a.pid() 
                                   << " got " << m_next_tag_size << '\n'
                << "  -> process " << p << " got " << tag_size << '\n';
    } 
    } // if (!dummy)
    
    m_send_empty = true;

    // update tag size
    m_recv_tag_size = m_send_tag_size;
    m_send_tag_size = m_next_tag_size;
    m_set_tag_size_counter = 0;
}



}
