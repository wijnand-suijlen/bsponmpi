#ifndef BSPONMPI_BSMP_H
#define BSPONMPI_BSMP_H

#include "aligned_buf.h"
#include "a2a.h"
#include "dllexport.h"

#ifdef PROFILE
#include "tictoc.h"
#endif

namespace bsplib {

class DLL_LOCAL Bsmp {
public:
    explicit Bsmp(MPI_Comm comm, size_t max_msg_size,
            size_t small_exch_size ) ;

    size_t set_tag_size( size_t new_size ) 
    {
        m_send_empty = false;
        m_set_tag_size_counter += 1;
        size_t old_tag_size = m_next_tag_size;
        m_next_tag_size = new_size;
        return old_tag_size;
    }

    size_t recv_tag_size() const 
    { return m_recv_tag_size; }

    size_t send_tag_size() const 
    { return m_send_tag_size; }

    void send( int pid, const void * tag, const void * payload, size_t nbytes )
    {
#ifdef PROFILE
        TicToc t(TicToc::BSMP);
        size_t start = m_a2a.send_size(pid);
#endif
        m_send_empty = false;
        serial( m_a2a, pid, nbytes + 1 );
        m_a2a.send( pid, tag, m_send_tag_size );
        m_a2a.send( pid, payload, nbytes );
#ifdef PROFILE
        size_t end = m_a2a.send_size(pid);
        t.add_bytes(end-start);
#endif
    }

    size_t n_total_messages() const 
    { return m_total_n_messages; }

    size_t total_payload() const
    { return m_total_payload; }

    bool empty() const
    { return m_payloads.empty(); }

    void pop()
    { 
        m_total_n_messages -= 1;
        m_total_payload -= m_payloads.back().second;
        return m_payloads.pop_back(); 
    }

    void * tag() 
    {   return m_tag_buffer.entry( m_payloads.size() - 1, m_recv_tag_size ); }

    void * payload()
    {   
        return static_cast<char *>(m_payload_buffer.data()) 
            + m_payloads.back().first;
    }

    size_t payload_size() const
    {   return m_payloads.back().second; }

    bool is_dummy() const 
    { return m_send_empty; }

    void sync(bool dummy) ;

private:
    size_t m_set_tag_size_counter;
    size_t m_recv_tag_size;
    size_t m_send_tag_size;
    size_t m_next_tag_size;
    size_t m_total_n_messages;
    size_t m_total_payload;
    bool m_send_empty;

    A2A m_a2a;
    AlignedBuf m_tag_buffer;
    std::vector< std::pair<size_t, size_t> > m_payloads;
    AlignedBuf m_payload_buffer;
};


}

#endif
