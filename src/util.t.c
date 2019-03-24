#include "util.h"
#include <stdio.h>
#include <stdlib.h>

#define TEST_EQ( type, format, a, b ) \
   do { type _x = (a);  \
        type _y = (b);  \
        if ( _x != _y ) { \
            fprintf( stderr, "%s:%d\n"                       \
                     "    Test    " #a " == " #b " fails\n"     \
                     "    because " format " != " format "\n",  \
                     __FILE__, __LINE__, _x, _y );           \
            abort();                                    \
        } } while(0)


int main( int argc, char ** argv )
{
    (void) argc;
    (void) argv;

    TEST_EQ( int, "%d", int_log2( 0 ), INT_MIN );
    TEST_EQ( int, "%d", int_log2( 1 ), 0 );
    TEST_EQ( int, "%d", int_log2( 2 ), 1 );
    TEST_EQ( int, "%d", int_log2( 3 ), 1 );
    TEST_EQ( int, "%d", int_log2( 4 ), 2 );
    TEST_EQ( int, "%d", int_log2( 5 ), 2 );
    TEST_EQ( int, "%d", int_log2( 7 ), 2 );
    TEST_EQ( int, "%d", int_log2( 8 ), 3 );
    TEST_EQ( int, "%d", int_log2( 15 ), 3 );
    TEST_EQ( int, "%d", int_log2( 16 ), 4 );
    TEST_EQ( int, "%d", int_log2( 32 ), 5 );
    TEST_EQ( int, "%d", int_log2( 33 ), 5 );
    TEST_EQ( int, "%d", int_log2( 127 ), 6 );
    TEST_EQ( int, "%d", int_log2( 128 ), 7 );
    TEST_EQ( int, "%d", int_log2( 129 ), 7 );
    TEST_EQ( int, "%d", int_log2( 300 ), 8 );
    TEST_EQ( int, "%d", int_log2( 70000 ), 16 );
    TEST_EQ( int, "%d", int_log2( 20000000 ), 24 );
    TEST_EQ( int, "%d", int_log2( UINT32_MAX ) , 31 );
    TEST_EQ( int, "%d", uint64_log2( UINT64_MAX ) , 63 );

    TEST_EQ( int, "%d", int_log(50, 0), INT_MIN );
    TEST_EQ( int, "%d", int_log(3, 1), 0 );
    TEST_EQ( int, "%d", int_log(4, 4), 1 );
    TEST_EQ( int, "%d", int_log(5, 25), 2);
    TEST_EQ( int, "%d", int_log(UINT_MAX, UINT_MAX), 0 );
    TEST_EQ( int, "%d", int_log(3, 81), 4 );
    TEST_EQ( int, "%d", int_log(3, 80), 3 );
    TEST_EQ( int, "%d", int_log(3, 160), 4 );
    TEST_EQ( int, "%d", int_log(15, 1248234), 5 );

    TEST_EQ( int, "%d", int_pow(0, UINT_MAX), 0);
    TEST_EQ( int, "%d", int_pow(1, 10), 1 );
    TEST_EQ( int, "%d", int_pow(2, 31), 1<<31);
    TEST_EQ( int, "%d", int_pow(3, 4), 81);
    TEST_EQ( int, "%d", int_pow(5, 3), 125 );

    return 0;
}
