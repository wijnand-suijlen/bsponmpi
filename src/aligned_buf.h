#ifndef BSPONMPI_ALIGNED_BUF_H
#define BSPONMPI_ALIGNED_BUF_H

#if __cplusplus >= 201103L
#include <cstddef>
#endif

#include <cstdlib>
#include <cstring>
#include <algorithm>

namespace bsplib {

class AlignedBuf {
public:
    explicit AlignedBuf()
#if __cplusplus >= 201103L
        : m_alignment( alignof( std::max_align_t))
#else
        : m_alignment( sizeof(void*) )
#endif
        , m_ptr( NULL )
        , m_aligned_ptr( NULL )
        , m_capacity( 0 )
        , m_size( 0 )
        {}

    size_t capacity() const 
    { return m_capacity; }

    size_t size() const
    { return m_size; }

    void * data()
    { return m_aligned_ptr; }

    const void * data() const 
    { return m_aligned_ptr; }

    void * append( size_t extra_bytes ) {
        extra_bytes = align( extra_bytes );
        if ( m_size + extra_bytes > m_capacity ) {
            size_t newsize = std::max( 2 * m_capacity, m_size + extra_bytes );
            void * ptr = std::malloc( newsize + 2*m_alignment );
            void * aligned_ptr = align( ptr );
            std::memcpy( aligned_ptr, m_aligned_ptr, m_size );

            std::free( m_ptr );
            m_ptr = ptr;
            m_aligned_ptr = aligned_ptr; 
            m_capacity = newsize;
        }

        void * result = static_cast<char*>( m_aligned_ptr ) + m_size;
        m_size += extra_bytes;
        return result;
    }

    void * entry( size_t n, size_t elem_size ) {
        return static_cast<char *>(m_aligned_ptr) + n * align(elem_size); 
    }

   
    void clear() 
    {
        m_size = 0;
    }

private:
    AlignedBuf( const AlignedBuf & ); // copying prohibited
    AlignedBuf & operator=(const AlignedBuf & ); // assignment prohibited

    size_t align( size_t x ) const {
        return ( x + m_alignment - 1)/m_alignment * m_alignment;
    }
    char * align( char * x ) const {
        char * null = NULL;
        return null + align( size_t(x - null) );
    }

    void * align( void * x ) const
    { return align( static_cast<char *>(x) ); }

    const size_t m_alignment;
    void * m_ptr;
    void * m_aligned_ptr;
    size_t m_capacity;
    size_t m_size;
};

}

#endif

