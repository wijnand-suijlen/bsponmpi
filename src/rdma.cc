#include "rdma.h"
#include "exception.h"
#ifdef PROFILE
#include "tictoc.h"
#endif

#include <cstring>
#include <algorithm>

namespace bsplib {

Rdma :: Rdma( MPI_Comm comm, size_t max_msg_size, size_t small_exch_size,
      double alpha, double beta  )
    : m_first_exchange( comm, max_msg_size, small_exch_size, alpha, beta )
    , m_second_exchange( comm, max_msg_size, small_exch_size, alpha, beta )
    , m_pid(m_first_exchange.pid())
    , m_nprocs(m_first_exchange.nprocs())
    , m_register()
    , m_used_slots()
    , m_free_slots()
    , m_cached_slot( no_slot() )
    , m_send_actions( m_nprocs )
    , m_recv_actions( m_nprocs )
    , m_unbuf( max_msg_size, comm )
    , m_send_push_pop_comm_buf()
    , m_recv_push_pop_comm_buf()
#ifdef PROFILE
    , m_tictoc( TicToc :: TOTAL )
#endif
{}


void Rdma::push_reg( void * addr, size_t size )
{
    Memslot id = no_slot();
    if ( m_free_slots.empty() ) {
        id = m_used_slots.size() / m_nprocs;

        if ( m_used_slots.capacity() < m_used_slots.size() + m_nprocs )
            m_used_slots.reserve( std::max(
                        2*m_used_slots.capacity(),
                        m_used_slots.size() + m_nprocs ) );

        m_used_slots.resize( m_used_slots.size() + m_nprocs );
    }
    else {
        id = m_free_slots.back();
        m_free_slots.pop_back();
    }

    Memblock block = { addr, size, Memblock::PUSHED };
    slot( m_pid, id ) = block;
    m_register[ addr ].push_back( id );

    PushPopCommBuf::PushEntry entry = { block, id };
    m_send_push_pop_comm_buf.m_pushed_slots.push_back( entry );
}

void Rdma::pop_reg( Memslot id )
{
    assert( id == null_slot() || !( slot(m_pid, id ).status & Memblock::POPPED) );
    if ( id != null_slot() ) {
        slot( m_pid, id ).status |= Memblock::POPPED;
    }
    m_send_push_pop_comm_buf.m_popped_slots.push_back( id );
}

void Rdma::put( const void * src,
        int dst_pid, Memslot dst_slot, size_t dst_offset, size_t size )
{
#ifdef PROFILE
    TicToc t( TicToc::PUT );
    size_t t_a = m_second_exchange.send_size(dst_pid);
#endif

    assert( ! (slot( m_pid, dst_slot ).status & Memblock::PUSHED) );

    serial( m_second_exchange, dst_pid, size );
    serial( m_second_exchange, dst_pid, dst_slot );
    serial( m_second_exchange, dst_pid, dst_offset );

    m_second_exchange.send( dst_pid, src, size );
#ifdef PROFILE
    size_t t_b = m_second_exchange.send_size(dst_pid);
    t.add_bytes(t_b-t_a);
#endif
}

void Rdma::hpput( const void * src,
        int dst_pid, Memslot dst_slot, size_t dst_offset, size_t size )
{
#ifdef PROFILE
    TicToc t( TicToc::HPPUT );
#endif
    assert( ! (slot( m_pid, dst_slot ).status & Memblock::PUSHED) );
    
    Memslot src_slot = m_local_slots.size();
    m_local_slots.push_back( const_cast<void *>(src) );

    int tag = m_unbuf.send( dst_pid, src, size );

    Action action = { Action::HPPUT, dst_pid, m_pid, dst_pid,
                      tag, src_slot, dst_slot , dst_offset, size };
    m_send_actions.push_back( action );
}

void Rdma::get( int src_pid, Memslot src_slot, size_t src_offset,
        void * dst, size_t size )
{
#ifdef PROFILE
    TicToc t( TicToc::GET );
#endif
    assert( !( slot( m_pid, src_slot ).status & Memblock::PUSHED) );

    Memslot dst_slot = m_local_slots.size();
    m_local_slots.push_back( dst );

    Action action = { Action::GET, src_pid, src_pid, m_pid, 
                      0, src_slot, dst_slot, src_offset, size  };
    m_send_actions.push_back( action );
}

void Rdma::hpget( int src_pid, Memslot src_slot, size_t src_offset,
        void * dst, size_t size )
{
#ifdef PROFILE
    TicToc t( TicToc::HPGET );
#endif
    assert( !( slot( m_pid, src_slot ).status & Memblock::PUSHED) );
    
    Memslot dst_slot = m_local_slots.size();
    m_local_slots.push_back( dst );

    int tag = m_unbuf.recv( src_pid, dst, size );
    Action action = { Action::HPGET, src_pid, src_pid, m_pid, 
                      tag, src_slot, dst_slot, src_offset, size };
    m_send_actions.push_back( action );

}

void Rdma::write_gets()
{
#ifdef PROFILE
    TicToc t(TicToc::GET );
#endif
    for (int p = 0; p < m_nprocs; ++p){
        size_t start = m_recv_actions.m_get_buffer_offset[p];
        m_second_exchange.recv_pop( p, start );
        while ( m_second_exchange.recv_size( p ) > 0 ) {
            size_t size, slot;
            deserial( m_second_exchange, p, size );
            deserial( m_second_exchange, p, slot);
           
            char * dst_addr = static_cast<char *>(m_local_slots[ slot ]);

            std::memcpy( dst_addr, m_second_exchange.recv_top(p), size );
            m_second_exchange.recv_pop(p, size );
        }
#ifdef PROFILE
        size_t end = m_second_exchange.recv_pos(p);
        t.add_bytes( end - start );
#endif
   }
}


void Rdma::write_puts()
{
#ifdef PROFILE
    TicToc t(TicToc::PUT);
#endif
    for (int p = 0; p < m_nprocs; ++p){
        size_t end = m_recv_actions.m_get_buffer_offset[p];
        m_second_exchange.recv_rewind( p );
        while ( m_second_exchange.recv_pos( p ) < end ) {
            size_t size, dst_slot, dst_offset;
            deserial( m_second_exchange, p, size );
            deserial( m_second_exchange, p, dst_slot);
            deserial( m_second_exchange, p, dst_offset );
           
            char * dst_addr =
               static_cast<char *>(slot( m_pid, dst_slot ).addr );
            dst_addr += dst_offset;
            std::memcpy( dst_addr, m_second_exchange.recv_top(p), size );
            m_second_exchange.recv_pop(p, size );
        }
#ifdef PROFILE
        t.add_bytes( end );
#endif
   }
}

bool Rdma::sync(bool dummy_bsmp)
{
    m_send_actions.set_dummy_bsmp(dummy_bsmp);
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

    bool all_dummy_bsmp = m_recv_actions.get_dummy_bsmp();
    // reset to initial state
    m_send_push_pop_comm_buf.clear();
    m_recv_push_pop_comm_buf.clear();
    m_send_actions.clear();
    m_recv_actions.clear();
    m_local_slots.clear();
    return all_dummy_bsmp;
}


void Rdma::ActionBuf::serialize( A2A & a2a ) 
{
    for (int p = 0; p < a2a.nprocs(); ++p ) {
        serial( a2a, p, m_dummy_bsmp );
        serial( a2a, p, m_get_buffer_offset[ p ] );
        serial( a2a, p, m_counts[ p ] );
    }

    for (size_t i = 0; i < m_actions.size(); ++i ) 
    {
        Action a = m_actions[i];
#ifdef PROFILE
        TicToc t( a.kind==Action::GET? TicToc::GET:
                  a.kind==Action::HPPUT? TicToc::HPPUT:
                  TicToc::HPGET );
        size_t start = a2a.send_size(a.target_pid);
#endif
        serial( a2a, a.target_pid, unsigned(a.kind) );
        serial( a2a, a.target_pid, a.src_pid );
        serial( a2a, a.target_pid, a.dst_pid );
        serial( a2a, a.target_pid, a.tag );
        serial( a2a, a.target_pid, a.src_slot );
        serial( a2a, a.target_pid, a.dst_slot );
        serial( a2a, a.target_pid, a.offset );
        serial( a2a, a.target_pid, a.size );
#ifdef PROFILE
        size_t end = a2a.send_size(a.target_pid);
        t.add_bytes(end-start);
#endif
    }
}

void Rdma::ActionBuf::deserialize( A2A & a2a ) 
{
    size_t n = 0;
    m_dummy_bsmp = true;
    for (int p = 0; p < a2a.nprocs(); ++p )
    {
        unsigned dummy;
        deserial( a2a, p, dummy );
        deserial( a2a, p, m_get_buffer_offset[ p ] );
        deserial( a2a, p, m_counts[ p ] );
        n += m_counts[p];
        m_dummy_bsmp = m_dummy_bsmp && dummy;
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
#ifdef PROFILE
            TicToc t( a.kind==Action::GET? TicToc::GET:
                  a.kind==Action::HPPUT? TicToc::HPPUT:
                  TicToc::HPGET );
#endif
            deserial( a2a, p, a.src_pid );
            deserial( a2a, p, a.dst_pid );
            deserial( a2a, p, a.tag );
            deserial( a2a, p, a.src_slot );
            deserial( a2a, p, a.dst_slot );
            deserial( a2a, p, a.offset );
            deserial( a2a, p, a.size );
        }
    }
    assert( i == m_actions.size() );
}

void Rdma::ActionBuf::execute( Rdma & rdma )
{
    const int pid = rdma.m_pid;
    for (size_t i = 0; i < m_actions.size(); ++i ) 
    {
        Action a = m_actions[i];
#ifdef PROFILE
        TicToc t( a.kind==Action::GET? TicToc::GET:
                  a.kind==Action::HPPUT? TicToc::HPPUT:
                  TicToc::HPGET );
#endif
        switch( a.kind ) {
            case Action::GET : {
#ifdef PROFILE
                size_t start=rdma.m_second_exchange.send_size(a.dst_pid);
#endif
                char * addr = static_cast< char *>( 
                        rdma.slot( pid, a.src_slot ).addr );
                addr += a.offset;
                serial( rdma.m_second_exchange, a.dst_pid, a.size );
                serial( rdma.m_second_exchange, a.dst_pid, a.dst_slot );
                rdma.m_second_exchange.send( a.dst_pid, addr, a.size );
#ifdef PROFILE
                size_t end=rdma.m_second_exchange.send_size(a.dst_pid);
                t.add_bytes(end-start);
#endif
                break;
            }

            case Action::HPPUT : {
                char * addr = static_cast< char *>(
                        rdma.slot( pid, a.dst_slot ).addr );
                addr += a.offset;
                rdma.m_unbuf.recv( a.tag, a.src_pid, addr, a.size );
#ifdef PROFILE
                t.add_bytes(a.size);
#endif
                break;
            }

            case Action::HPGET : {
                char * addr = static_cast< char *>(
                        rdma.slot( pid, a.src_slot ).addr );
                addr += a.offset;
                rdma.m_unbuf.send( a.tag, a.dst_pid, addr, a.size );
#ifdef PROFILE
                t.add_bytes(a.size);
#endif
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
            PushEntry entry = m_pushed_slots[s];
            char * null = NULL;
            size_t addr = static_cast<char *>(entry.block.addr) - null;
            size_t size = entry.block.size;
            size_t id = entry.slot;
            
            serial( a2a, p, addr );
            serial( a2a, p, size );
            serial( a2a, p, id );
        }
    } }

    // broadcast popped slots
    { SizeSer::Buffer nbuf;
      const int n_nbuf = SizeSer::write( n_popped_slots, nbuf );
      for (int p = 0; p < a2a.nprocs(); ++p) {
        a2a.send( p, nbuf, n_nbuf );

        for (size_t s = 0; s < m_popped_slots.size(); ++s) {
            Memslot slotid = m_popped_slots[s];
            serial( a2a, p, slotid );
        }
    } }
}

void Rdma::PushPopCommBuf::deserialize( A2A & a2a ) 
{
    size_t n_pushed_slots = 0;

    // read pushed slots
    for (int p = 0; p < a2a.nprocs(); ++p) {
        size_t n;
        deserial( a2a, p, n );

        if ( p == 0 ) {
            n_pushed_slots = n;
            m_pushed_slots.resize( a2a.nprocs() * n );
        } else if (n != n_pushed_slots )
            throw exception("bsp_push_reg") <<
                "Not all processes registered the same number of memory "
                "blocks. For example, process 0 registered " << n_pushed_slots
                << " blocks, while process " << p << " registerd " << n;

        for (size_t s = 0; s < n_pushed_slots; ++s) {
            char * null = NULL;
            size_t addr, size, id;
            deserial( a2a, p, addr);
            deserial( a2a, p, size);
            deserial( a2a, p, id );

            PushEntry entry = { { null + addr, size, Memblock::NORMAL }, id };
            m_pushed_slots[ p + s * a2a.nprocs() ] = entry ;
        }
    }

    // read popped slots
    size_t n_popped_slots;
    for (int p = 0; p < a2a.nprocs(); ++p) {
        size_t n;
        deserial( a2a, p, n );

        if ( p == 0 ) {
            n_popped_slots = n;
            m_popped_slots.resize( n, no_slot() );
        } else if (n != n_popped_slots )
            throw exception("bsp_pop_reg") <<
                "Not all processes popped the same number of memory "
                "blocks. For example, process 0 registered " << n_popped_slots
                << " blocks, while process " << p << " registerd " << n;

        for (size_t s = 0; s < n_popped_slots; ++s) {
            Memslot slotid = no_slot();
            deserial( a2a, p, slotid );
            assert( slotid != no_slot() );

            if ( slotid != null_slot() ) {
                if ( m_popped_slots[s] == no_slot() ) {
                    m_popped_slots[s] = slotid;
                }
                else {
                    if (m_popped_slots[s] != slotid ) {
                        throw exception("bsp_pop_reg") <<
                            "Processes 0 and " << p << " did not pop "
                            "the same memory registration";
                    }
                }
            }
        }
    }
    
}

void Rdma::PushPopCommBuf::execute( Rdma & rdma ) 
{
    const int pid = rdma.m_pid;
    const int nprocs = rdma.m_nprocs;
 
    /* do memory registrations */
    for (size_t s = 0; s < m_pushed_slots.size()/nprocs; ++s) {
        Memslot id = m_pushed_slots[s*nprocs].slot;

        for ( int p = 0 ; p < nprocs; ++p ) {
            if ( m_pushed_slots[ s*nprocs + p ].slot != id )
                throw exception("bsp_push_reg") << "Internal error; inconsistent slot id";

            rdma.slot( p, id ).addr = m_pushed_slots[ s*nprocs + p ].block.addr;
            rdma.slot( p, id ).size = m_pushed_slots[ s*nprocs + p ].block.size;
            rdma.slot( p, id ).status &= ~Memblock::PUSHED; 
        }
    }

    /* do memory deregistrations */
    for (size_t s = 0; s < m_popped_slots.size(); ++s) {
        Memslot id = m_popped_slots[s];

        Reg::iterator reg_entry;
        RegStack :: iterator rs; 

        if ( id != no_slot() ) {
            assert( rdma.slot(pid,id).addr == NULL
                   || rdma.slot( pid, id ).status & Memblock::POPPED );
            void * addr = rdma.slot( pid, id).addr;
            reg_entry = rdma.m_register.find( addr );
            assert( reg_entry != rdma.m_register.end() );
            RegStack & stack = reg_entry->second;
            RegStack :: reverse_iterator j ; 
            for (j = stack.rbegin(); j != stack.rend(); ++j )
                if ( *j == id ) break;

            if ( j == stack.rend() )
                throw exception("bsp_pop_reg") << "Internal error: inconsistent slot id";

            rs = j.base();
            --rs;
        }
        else
        { // then all processes have popped NULL
          // Let's delete the first full block of NULL we can encounter

            reg_entry = rdma.m_register.find( NULL );
            if( reg_entry == rdma.m_register.end() )
                throw exception("bsp_pop_reg") <<
                    "Tried to deregister NULL on all processes, but NULL was never registered";

            RegStack & stack = reg_entry->second;

            for ( rs = stack.begin(); rs != stack.end(); ++rs ) {

                id = *rs;

                bool all_nulls = true;
                for ( int p = 0; p < nprocs; ++p ) {
                    if ( rdma.slot( p, id ).addr != NULL ) {
                        all_nulls = false;
                        break;
                    }
                }

                if ( all_nulls ) {
                   break;
                }
            }
            if ( rs == stack.end() ) 
                throw exception("bsp_pop_reg") <<
                    "Tried to deregister NULL on all processes, but could not find a matching regisration (bsp_push_reg)";
        }

        assert( *rs == id );
        reg_entry->second.erase( rs );

        if ( reg_entry->second.empty())
            rdma.m_register.erase( reg_entry );

        for (int p = 0; p < nprocs; ++p) {
            rdma.slot( p, id ).addr = NULL;
            rdma.slot( p, id ).size = 0;
            rdma.slot( p, id ).status = Memblock::NORMAL;
        }
        rdma.m_free_slots.push_back( id );
    }
}


}
