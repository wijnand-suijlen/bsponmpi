#ifndef BSPONMPI_UINTSERIALIZE_H
#define BSPONMPI_UINTSERIALIZE_H

#include <climits>
#include <cassert>

namespace bsplib {

// Helper class to serialize and deserialize integers
// small values need less space
template <class UInt>
struct UIntSerialize{ 
    enum { value = (sizeof(UInt)*CHAR_BIT + CHAR_BIT - 2 )  / (CHAR_BIT-1) };

    typedef unsigned char Buffer[value];

    static int write( UInt x, Buffer & bytes ) {
        assert( x >= 0 );
        int i = 0;
        do {
            unsigned char m = (1 << (CHAR_BIT-1));
            unsigned char c = x & (m-1); // cut out least significant bits
            x >>= (CHAR_BIT-1);          // and save these in c
            c |= (x > 0) * m ;  // set highest signif. bit if continuing
            bytes[i++] = c;
        } while ( x != 0 );
        return i;
    }

    static int read( const unsigned char * cs, UInt & x ) {
        x = 0;
        int i = 0;
        unsigned char m = (1 << (CHAR_BIT-1));
        unsigned char c = 0;
        UInt power = 1;
        do {
            c = cs[i++];
            x += power * (c & (m-1));
            power <<= (CHAR_BIT-1);
        } while (c & m);

        return i;
    }
};

}

#endif
