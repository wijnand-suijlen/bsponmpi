#include "rdma.h"
#include "exception.h"

#include <cstring>

namespace bsplib {

void Rdma::push_reg( void * addr, size_t size )
{
    Memslot slot = {addr, size};
    m_send_push_pop_comm_buf.m_pushed_slots.push_back( slot );
}

void Rdma::pop_reg( MemslotID slot )
{
    m_send_push_pop_comm_buf.m_popped_slots.push_back( slot );
}

void Rdma::put( const void * src,
        int dst_pid, MemslotID dst_slot, size_t dst_offset, size_t size )
{
    char * addr = static_cast<char *>( slot(dst_pid, dst_slot).addr );
    char * null = NULL;
    size_t dst_addr = addr - null + dst_offset;

    serial( m_second_exchange, dst_pid, size );
    serial( m_second_exchange, dst_pid, dst_addr );

    m_second_exchange.send( dst_pid, src, size );
}

void Rdma::hpput( const void * src,
        int dst_pid, MemslotID dst_slot, size_t dst_offset, size_t size )
{
    char * addr = static_cast<char *>( slot(dst_pid, dst_slot).addr );
    char * null = NULL;
    size_t dst_addr = addr - null + dst_offset;

    size_t src_addr = static_cast<const char *>(src) - null;

    int tag = m_unbuf.send( dst_pid, src, size );

    Action action = { Action::HPPUT, dst_pid, m_pid, dst_pid,
                      tag, src_addr, dst_addr, size };
    m_send_actions.push_back( action );
   
}

void Rdma::get( int src_pid, MemslotID src_slot, size_t src_offset,
        void * dst, size_t size )
{
    char * null = NULL;
    size_t dst_addr = static_cast<char *>(dst) - null;

    char * addr = static_cast<char *>( slot( src_pid, src_slot).addr );
    size_t src_addr = addr - null + src_offset;

    Action action = { Action::GET, src_pid, src_pid, m_pid, 
                      0, src_addr, dst_addr, size  };
    m_send_actions.push_back( action );
}

void Rdma::hpget( int src_pid, MemslotID src_slot, size_t src_offset,
        void * dst, size_t size )
{
    char * null = NULL;
    size_t dst_addr = static_cast<char *>(dst) - null;

    char * addr = static_cast<char *>( slot( src_pid, src_slot ).addr );
    size_t src_addr = addr - null + src_offset;

    int tag = m_unbuf.recv( src_pid, dst, size );
    Action action = { Action::HPGET, src_pid, src_pid, m_pid, 
                      tag, src_addr, dst_addr, size  };
    m_send_actions.push_back( action );

}

void Rdma::write_gets()
{
    for (int p = 0; p < m_nprocs; ++p){
        m_second_exchange.recv_pop( p, m_recv_actions.m_get_buffer_offset[p] );
        while ( m_second_exchange.recv_size( p ) > 0 ) {
            size_t size, addr;
            deserial( m_second_exchange, p, size );
            deserial( m_second_exchange, p, addr );
           
            char * dst_addr = NULL;
            dst_addr += addr;

            std::memcpy( dst_addr, m_second_exchange.recv_top(p), size );
            m_second_exchange.recv_pop(p, size );
        }
   }
}


void Rdma::write_puts()
{
    for (int p = 0; p < m_nprocs; ++p){
        m_second_exchange.recv_rewind( p );
        while ( m_second_exchange.recv_size( p ) > 0 ) {
            size_t size, addr;
            deserial( m_second_exchange, p, size );
            deserial( m_second_exchange, p, addr );
           
            char * dst_addr = NULL;
            dst_addr += addr;

            std::memcpy( dst_addr, m_second_exchange.recv_top(p), size );
            m_second_exchange.recv_pop(p, size );
        }
   }
}

void Rdma::sync()
{
    m_send_actions.set_get_buffer_offset( m_second_exchange );

    // first exchange: gets, unbuffered requests, push/pop registers
    m_send_push_pop_comm_buf.serialize( m_first_exchange );
    m_send_actions.serialize( m_first_exchange );
   
    m_first_exchange.exchange( );

    m_recv_push_pop_comm_buf.deserialize( m_first_exchange );
    m_recv_actions.deserialize( m_first_exchange );
   
    m_recv_actions.execute(*this);
    
    // start the unbuffered requests
    m_unbuf.start();

    // perform the buffered puts (which includes the buffered gets)
    m_second_exchange.exchange( );
    write_gets();
    write_puts() ;

    // wait for the unbuffered requests to finish
    m_unbuf.wait();

    // execute the push/pop registers
    m_recv_push_pop_comm_buf.execute(*this);

    // reset to initial state
    m_send_push_pop_comm_buf.clear();
    m_recv_push_pop_comm_buf.clear();
    m_send_actions.clear();
    m_recv_actions.clear();
}


void Rdma::ActionBuf::serialize( A2A & a2a ) 
{
    for (int p = 0; p < a2a.nprocs(); ++p ) {
        serial( a2a, p, m_get_buffer_offset[ p ] );
        serial( a2a, p, m_counts[ p ] );
    }

    for (size_t i = 0; i < m_actions.size(); ++i ) 
    {
        Action a = m_actions[i];
        serial( a2a, a.target_pid, unsigned(a.kind) );
        serial( a2a, a.target_pid, a.src_pid );
        serial( a2a, a.target_pid, a.dst_pid );
        serial( a2a, a.target_pid, a.tag );
        serial( a2a, a.target_pid, a.src_addr );
        serial( a2a, a.target_pid, a.dst_addr );
        serial( a2a, a.target_pid, a.size );
    }
}

void Rdma::ActionBuf::deserialize( A2A & a2a ) 
{
    size_t n = 0;
    for (int p = 0; p < a2a.nprocs(); ++p )
    {
        deserial( a2a, p, m_get_buffer_offset[ p ] );
        deserial( a2a, p, m_counts[ p ] );
        n += m_counts[p];
    }
    m_actions.resize( n );

    size_t i = 0;
    for (int p = 0; p < a2a.nprocs(); ++p )
    {
        for (size_t j = 0; j < m_counts[p]; ++j) {
            Action & a = m_actions[i++];
            a.target_pid = p;
            unsigned kind;
            deserial( a2a, p, kind );
            a.kind = Action::Kind(kind);
            deserial( a2a, p, a.src_pid );
            deserial( a2a, p, a.dst_pid );
            deserial( a2a, p, a.tag );
            deserial( a2a, p, a.src_addr );
            deserial( a2a, p, a.dst_addr );
            deserial( a2a, p, a.size );
        }
    }
    assert( i == m_actions.size() );
}

void Rdma::ActionBuf::execute( Rdma & rdma )
{
    for (size_t i = 0; i < m_actions.size(); ++i ) 
    {
        Action a = m_actions[i];
        switch( a.kind ) {
            case Action::GET : {
                serial( rdma.m_second_exchange, a.dst_pid, a.size );
                serial( rdma.m_second_exchange, a.dst_pid, a.dst_addr );
                char * src = NULL;
                src += a.src_addr;
                rdma.m_second_exchange.send( a.dst_pid, src, a.size );
                break;
            }

            case Action::HPPUT : {
                char * null = NULL;
                rdma.m_unbuf.recv( a.tag, a.src_pid, null+a.dst_addr, a.size );
                break;
            }

            case Action::HPGET : {
                char * null = NULL;
                rdma.m_unbuf.send( a.tag, a.dst_pid, null + a.src_addr, a.size );
                break;
            }
        }
    }
}


void Rdma::PushPopCommBuf::serialize( A2A & a2a ) 
{
    typedef UIntSerialize< size_t >    SizeSer;

    const size_t n_pushed_slots = m_pushed_slots.size();
    const size_t n_popped_slots = m_popped_slots.size();

    // broadcast pushed slots
    { SizeSer::Buffer nbuf;
      const int n_nbuf = SizeSer::write( n_pushed_slots, nbuf );
      for (int p = 0; p < a2a.nprocs(); ++p) {
        a2a.send( p, nbuf, n_nbuf );

        for (size_t s = 0; s < m_pushed_slots.size(); ++s) {
            Memslot slot = m_pushed_slots[s];
            char * null = NULL;
            size_t addr = static_cast<char *>(slot.addr) - null;
            size_t size = slot.size;
            
            serial( a2a, p, addr );
            serial( a2a, p, size );
        }
    } }

    // broadcast popped slots
    { SizeSer::Buffer nbuf;
      const int n_nbuf = SizeSer::write( n_popped_slots, nbuf );
      for (int p = 0; p < a2a.nprocs(); ++p) {
        a2a.send( p, nbuf, n_nbuf );

        for (size_t s = 0; s < m_popped_slots.size(); ++s) {
            MemslotID slotid = m_popped_slots[s];
            serial( a2a, p, slotid );
        }
    } }
}

void Rdma::PushPopCommBuf::deserialize( A2A & a2a ) 
{
    size_t n_pushed_slots;

    // read pushed slots
    for (int p = 0; p < a2a.nprocs(); ++p) {
        size_t n;
        deserial( a2a, p, n );

        if ( p == 0 ) {
            n_pushed_slots = n;
            m_pushed_slots.resize( a2a.nprocs() * n );
        } else if (n != n_pushed_slots )
            throw exception("bsp_push_reg") <<
                ": Not all processes registered the same number of memory "
                "blocks. For example, process 0 registered " << n_pushed_slots
                << " blocks, while process " << p << " registerd " << n;

        for (size_t s = 0; s < n_pushed_slots; ++s) {
            char * null = NULL;
            size_t addr, size;
            deserial( a2a, p, addr);
            deserial( a2a, p, size);

            Memslot slot = { null + addr, size };
            m_pushed_slots[ p + s * a2a.nprocs() ] = slot;
        }
    }

    size_t n_popped_slots;
    for (int p = 0; p < a2a.nprocs(); ++p) {
        size_t n;
        deserial( a2a, p, n );

        if ( p == 0 ) {
            n_popped_slots = n;
            m_popped_slots.resize( n, no_slot() );
        } else if (n != n_popped_slots )
            throw exception("bsp_push_reg") <<
                ": Not all processes registered the same number of memory "
                "blocks. For example, process 0 registered " << n_popped_slots
                << " blocks, while process " << p << " registerd " << n;

        for (size_t s = 0; s < n_popped_slots; ++s) {
            MemslotID slotid = no_slot();
            deserial( a2a, p, slotid );
            assert( slotid != no_slot() );

            if ( slotid != null_slot() ) {
                if ( m_popped_slots[s] == no_slot() ) {
                    m_popped_slots[s] = slotid;
                }
                else {
                    if (m_popped_slots[s] != slotid ) {
                        throw exception("bsp_pop_reg") <<
                            ": Processes 0 and " << p << " did not pop "
                            "the same memory registration";
                    }
                }
            }
        }
    }
    
}

void Rdma::PushPopCommBuf::execute( Rdma & rdma ) 
{
    /* do memory deregistrations */
    const int pid = rdma.m_pid;
    const int nprocs = rdma.m_nprocs;

    for (size_t s = 0; s < m_popped_slots.size(); ++s) {
        MemslotID slot = m_popped_slots[s];

        if ( slot != no_slot() ) {
            void * addr = rdma.m_used_slots[ pid + nprocs * slot ].addr;
            Reg::iterator reg_entry = rdma.m_register.find( addr );
            assert( reg_entry != rdma.m_register.end() );
            reg_entry->second.pop_back();
            if (reg_entry->second.empty())
                rdma.m_register.erase( reg_entry );
            
            for (int p = 0; p < nprocs; ++p) {
                rdma.m_used_slots[ slot * nprocs + p ].addr = NULL;
                rdma.m_used_slots[ slot * nprocs + p ].size = 0;
            }
            rdma.m_free_slots.push_back( slot );
        }
        else
        { // then all processes have popped NULL
          // Let's delete the first full block of NULL we can encounter

            Reg :: iterator i = rdma.m_register.find( NULL );
            assert( i != rdma.m_register.end() );

            std::list<MemslotID>::iterator j = i->second.begin();
            for ( ; j != i->second.end(); ++j) {

                slot = *j;

                bool all_nulls = true;
                for ( int p = 0; p < nprocs; ++p ) {
                    if ( rdma.m_used_slots[ slot * nprocs + p ].addr != NULL ) {
                        all_nulls = false;
                        break;
                    }
                }

                if ( all_nulls ) {
                    i->second.erase( j );
                    if (i->second.empty()) {
                        rdma.m_register.erase( i );
                    }
                    break;
                }
            }
            if ( j == i->second.end() ) 
                throw exception("bsp_pop_reg") <<
                    ": Tried to deregister NULL on all processes, but could not find a matching regisration (bsp_push_reg)";

        }
    }

    /* do memory registrations */
    for (size_t s = 0; s < m_pushed_slots.size()/nprocs; ++s) {

        MemslotID newslot = size_t(-1);
        if ( rdma.m_free_slots.empty() ) {
            newslot = rdma.m_used_slots.size() / nprocs;

            if ( rdma.m_used_slots.capacity() < rdma.m_used_slots.size() + nprocs )
                rdma.m_used_slots.reserve( std::max(
                            2*rdma.m_used_slots.capacity(),
                            rdma.m_used_slots.size() + nprocs ) );

            rdma.m_used_slots.resize( rdma.m_used_slots.size() + nprocs );
        }
        else {
            newslot = rdma.m_free_slots.back();
            rdma.m_free_slots.pop_back();
        }
        std::memcpy( &rdma.m_used_slots[ newslot * nprocs ], 
                     &m_pushed_slots[ s * nprocs], 
                      nprocs * sizeof(Memslot) );

        void * addr = m_pushed_slots[ s*nprocs + pid].addr;
        rdma.m_register[ addr ].push_back( newslot );
    }
}

bool Rdma::PushPopCommBuf::was_pushed( void * ptr ) const
{
    for (size_t i = 0; i < m_pushed_slots.size(); ++i)
        if ( m_pushed_slots[i].addr == ptr )
            return true;

    return false;
}


}
