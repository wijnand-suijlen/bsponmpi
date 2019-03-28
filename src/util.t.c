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

void test_math_funcs(void)
{
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
}

void test_rng(void) {
    int i, n = 1000;
    size_t seed = 0;
    int mean = 1000;
    int sigma = 32;
    int match_1sigma=0;
    int match_2sigma=0;
    int match_3sigma=0;
    int * buckets = calloc( n, sizeof(int) );

    /*** test uniformity ***/

    for (i = 0; i < n*n; ++i ) {
        unsigned j = n * rand_next(&seed);
        buckets[ j ] += 1 ;
    }

    /* because we do so many tests, we can approximate
     * the number in each bucket is normally distributed
     * with the same variance as a binomial distribution
     * with 1000000 trials. So sigma =approx 32. Roughly
     * 68% of the buckets should have a value
     * between 968 and 1032 
     * 95% of the buckets should be between 936 and 1064
     * 99.7% of the buckets between 902 and 1096*/

    for (i = 0; i < n; ++i ) {
        if ( buckets[i] > mean-sigma && buckets[i] < mean + sigma) {
            match_1sigma ++;
        }
        if ( buckets[i] > mean-2*sigma && buckets[i] < mean + 2*sigma) {
            match_2sigma ++;
        }
        if ( buckets[i] > mean-3*sigma && buckets[i] < mean + 3*sigma) {
            match_3sigma ++;
        }
        
        if (buckets[i] == 0 ) {
            fprintf(stderr, "RNG has hole at i = %d\n", i);
            abort();
        }
    }

    if ( match_1sigma < 650 ) {
        fprintf(stderr, "RNG failed 1 sigma criterium\n");
        abort();
    }

    if ( match_2sigma < 920 ) {
        fprintf(stderr, "RNG failed 2 sigma criterium\n");
        abort();
    }

    if ( match_3sigma < 990 ) {
        fprintf(stderr, "RNG failed 3 sigma criterium\n");
        abort();
    }

    free(buckets);
}

int int_is_equal( const void * a, const void * b )
{
    const int * x = a, * y = b;
    return *x == *y;
}

size_t int_hash( const void * a )
{
    const int * x = a;
    return *x;
}

void test_hash_table(void)
{
    int i, x, n = 1000;
    hash_table_t table;
    int * entries;
    hash_table_create( &table, 10, int_is_equal, int_hash );

    entries = calloc( n, sizeof(int) );

    for ( i = 0; i < n; ++i ) {
        entries[i] = i;
        hash_table_new_item( &table, &entries[i] );
    }

    if ( table.n_buckets > 1024 ) {
        fprintf(stderr, "Hash table has grown too quickly\n");
        abort();
    }

    TEST_EQ( int, "%d", table.n_items ,  n );    

    x = 10;
    TEST_EQ( int, "%d", * (int *)
            hash_table_get_item( &table, &x ), 10 );

    for ( i = 0; i < n; ++i ) {
        TEST_EQ( void *, "%p", hash_table_new_item( &table, &entries[i] ),
                &entries[i] );
    }

    for ( i = 0; i < n; ++i ) {
        TEST_EQ( void *, "%p", hash_table_delete_item(&table,&entries[i]),
                &entries[i] );
    }

    TEST_EQ( int, "%d", table.n_items ,  0 );    

    free(entries);
}

int main( int argc, char ** argv )
{
    (void) argc;
    (void) argv;

    test_math_funcs();
    test_rng();
    test_hash_table();

    return 0;
}
