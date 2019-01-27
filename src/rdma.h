#ifndef BSPONMPI_RDMA_H
#define BSPONMPI_RDMA_H

#include "a2a.h"
#include "unbuf.h"
#include "uintserialize.h"
#include "dllexport.h"

#ifdef HAS_CXX11_UNORDERED_MAP
#include <unordered_map>
#else
#include <tr1/unordered_map>
#endif

#include <mpi.h>

namespace bsplib {

class DLL_LOCAL Rdma { 
public:
    struct Memblock {
        void * addr;
        size_t size;
        enum Status { NORMAL=0x0, PUSHED=0x1, POPPED=0x2 };
        unsigned status;
    };
    typedef size_t Memslot;

    Rdma( MPI_Comm comm , size_t max_msg_size  )
        : m_first_exchange( comm, max_msg_size )
        , m_second_exchange( comm, max_msg_size )
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
    {}

    void push_reg( void * addr, size_t size );
    void pop_reg( Memslot slot );

    static Memslot no_slot() { return Memslot(-1); }
    static Memslot null_slot() { return Memslot(-2); }

    Memslot lookup_reg( const void * addr, bool pushed, bool popped ) const
    {   
        if (addr == NULL ) return null_slot();
        
        if ( m_cached_slot != no_slot() 
                && slot( m_pid, m_cached_slot).addr == addr
                && slot( m_pid, m_cached_slot).status == Memblock::NORMAL )
            return m_cached_slot;

        Reg::const_iterator i = m_register.find( const_cast<void*>(addr) );
        if (i != m_register.end()) {
            for (RRSIt j = i->second.rbegin(); j != i->second.rend(); ++j) {
                unsigned status  = slot( m_pid, *j ).status ;
                if ( status == unsigned(Memblock::NORMAL))  {
                    m_cached_slot = *j;
                    return m_cached_slot;
                }
                if ( pushed == bool (status & Memblock::PUSHED)
                  &&  popped == bool (status & Memblock::POPPED)) {
                    return *j;
                }
                             
            }
        }
        return no_slot(); 
    }

    const Memblock & slot( int pid, Memslot slot ) const
    { return m_used_slots[ slot * m_nprocs + pid ]; }

    Memblock & slot( int pid, Memslot slot )
    { return m_used_slots[ slot * m_nprocs + pid ]; }

    void put( const void * src, int dst_pid, Memslot dst_slot, size_t dst_offset,
            size_t size );

    void hpput( const void * src, int dst_pid, Memslot dst_slot, size_t dst_offset,
            size_t size );

    void get( int src_pid, Memslot src_slot, size_t src_offset, void * dst,
            size_t size );

    void hpget( int src_pid, Memslot src_slot, size_t src_offset, void * dst,
            size_t size );

    bool sync( bool dummy_bsmp );

private:
    typedef std::list< Memslot > RegStack;
    typedef RegStack :: const_iterator RSIt;
    typedef RegStack :: const_reverse_iterator RRSIt;

#ifdef HAS_CXX11_UNORDERED_MAP
    typedef std::unordered_map< void * , RegStack > Reg;
#else
    typedef std::tr1::unordered_map< void *, RegStack > Reg;
#endif
    A2A m_first_exchange;
    A2A m_second_exchange;
    int m_pid, m_nprocs;
    Reg m_register;
    std::vector< Memblock > m_used_slots;
    std::vector< Memslot > m_free_slots;
    mutable Memslot m_cached_slot;

    struct Action { 
        enum Kind { GET, HPPUT, HPGET } kind; 
        int target_pid, src_pid, dst_pid, tag;
        size_t src_addr, dst_addr;
        size_t size;
    };

    struct ActionBuf {
        ActionBuf( int nprocs )
            : m_actions()
            , m_counts( nprocs )
            , m_get_buffer_offset( nprocs )
        {}

        std::vector< Action > m_actions;
        std::vector< size_t > m_counts;
        std::vector< size_t > m_get_buffer_offset;
        unsigned m_dummy_bsmp;

        void push_back( const Action & action )
        {   
            m_actions.push_back( action );
            m_counts[ action.target_pid ] += 1;
        }

        void set_get_buffer_offset( A2A & buf ) 
        {
            for ( int p = 0; p < buf.nprocs(); ++p ) 
                m_get_buffer_offset[ p ] = buf.send_size( p );
        }

        void set_dummy_bsmp( bool dummy )
        { m_dummy_bsmp = dummy; }

        bool get_dummy_bsmp()
        { return m_dummy_bsmp; }

        void serialize( A2A & a2a );
        void deserialize( A2A & a2a );
        void execute( Rdma & rdma );
      
        void clear() 
        {
            m_dummy_bsmp = false;
            m_actions.clear();
            std::fill( m_counts.begin(), m_counts.end(), 0);
            std::fill( m_get_buffer_offset.begin(), m_get_buffer_offset.end(), 0);
        }
    };

    ActionBuf m_send_actions;
    ActionBuf m_recv_actions;
        

    /// Unbuffered communications
    Unbuf m_unbuf;
    
    // Comm structure to exchange push & pop of registers 
    struct PushPopCommBuf {
        struct PushEntry { Memblock block; Memslot slot; };
        std::vector< PushEntry > m_pushed_slots;
        std::vector< Memslot >   m_popped_slots;
        
        void clear()
        { m_pushed_slots.clear(); m_popped_slots.clear(); }

        void serialize( A2A & a2a );
        void deserialize( A2A & a2a );
        void execute( Rdma & rdma );
    };
    PushPopCommBuf m_send_push_pop_comm_buf;
    PushPopCommBuf m_recv_push_pop_comm_buf;

    void write_gets();
    void write_puts();

};

}

#endif
