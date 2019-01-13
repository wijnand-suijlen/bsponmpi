#include "uintserialize.h"

#include <cassert>
#include <limits>
#include <cstdio>
#include <cstring>

using namespace bsplib;

int main( int argc, char ** argv )
{
    (void) argc; (void) argv;

    UIntSerialize< unsigned int > intser;
    UIntSerialize< unsigned int > :: Buffer intbuf;
    UIntSerialize< unsigned int > :: Buffer readbuf;

    assert( sizeof(intbuf) == 5 );

    unsigned int y = 0;
    assert( 1 == intser.write( 0, intbuf ) );
    assert( intbuf[0] == '\0' );

    std::memset( readbuf, 0, sizeof(readbuf) );
    std::memcpy( readbuf, intbuf, 1 );
    assert( 1 == intser.read( readbuf, y ) );
    assert( y == 0 );

    assert( 1 == intser.write( 120, intbuf ) );
    assert( intbuf[0] == (unsigned char) 120 );

    std::memset( readbuf, 0, sizeof(readbuf) );
    std::memcpy( readbuf, intbuf, 1 );
    assert( 1 == intser.read( readbuf, y ) );
    assert( y == 120 );

    assert( 2 == intser.write( 256, intbuf ) );
    assert( intbuf[0] == (unsigned char) 128 );
    assert( intbuf[1] == (unsigned char) 2 );

    std::memset( readbuf, 0, sizeof(readbuf) );
    std::memcpy( readbuf, intbuf, 2 );
    assert( 2 == intser.read( readbuf, y ) );
    assert( y == 256 );

    unsigned uint_max = std::numeric_limits<unsigned>::max();
    assert( 5 == intser.write( uint_max, intbuf ) );
    assert( intbuf[0] == (unsigned char) 255 );
    assert( intbuf[1] == (unsigned char) 255 );
    assert( intbuf[2] == (unsigned char) 255 );
    assert( intbuf[3] == (unsigned char) 255 );
    assert( intbuf[4] == (unsigned char) 15 );
         
    std::memset( readbuf, 0, sizeof(readbuf) );
    std::memcpy( readbuf, intbuf, 5 );
    assert( 5 == intser.read( readbuf, y ) );
    assert( y == uint_max );

    return 0;
}
