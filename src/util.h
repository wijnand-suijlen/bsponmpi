#ifndef BSPONMPI_UTIL_H
#define BSPONMPI_UTIL_H

#include <stdint.h>
#include <stddef.h>
#include <limits.h>

#include "dllexport.h"

#define MAX( a, b ) ( (a) < (b) ? (b) : (a) )
#define MIN( a, b ) ( (a) > (b) ? (b) : (a) )

#ifdef __cplusplus
extern "C" {
#endif

#ifdef UINT64_MAX
DLL_LOCAL int uint64_log2( uint64_t x );
#endif    
DLL_LOCAL int uint32_log2( uint32_t x );

#ifdef UINT64_MAX
  #if UINT_MAX == UINT64_MAX
    #define int_log2(x) uint64_log2(x)
  #else
    #define int_log2(x) uint32_log2(x)
  #endif

  #if ULONG_MAX == UINT64_MAX
    #define long_log2( x ) uint64_log2( x )
  #else
    #define long_log2( x ) uint32_log2( x )
  #endif
#else
  #define int_log2(x)  uint32_log2( x)
  #define long_log2(x) uint32_log2( x )
#endif

DLL_LOCAL unsigned int_log( unsigned base, unsigned x );

DLL_LOCAL unsigned int_pow( unsigned base, unsigned n );


#ifdef UINT64_MAX
DLL_LOCAL double rand_next_uint64( uint64_t * next );
#endif

DLL_LOCAL double rand_next_uint32( uint32_t * next );

#ifdef UINT64_MAX
  #if SIZE_MAX == UINT64_MAX
     #define rand_next( next )  rand_next_uint64( next )
  #elif SIZE_MAX == UINT32_MAX
     #define rand_next( next )  rand_next_uint32( next )
  #endif
#else
  #define rand_next( next )  rand_next_uint32( next )
#endif


typedef struct universal_hash_function {
    size_t mult, add;
    int log2_buckets;
} universal_hash_function_t;

DLL_LOCAL universal_hash_function_t new_universal_hash_function( size_t * seed,
       unsigned buckets );


typedef struct hash_table_bucket hash_table_bucket_t;
typedef struct hash_table {
    unsigned n_items;
    unsigned n_buckets;
    hash_table_bucket_t * buckets;
    hash_table_bucket_t * free_buckets;
    hash_table_bucket_t * allocs;
    int (*is_equal)(const void * a, const void * b);
    size_t (*hash)(const void * x );
    size_t seed;
    universal_hash_function_t ghash;

    double max_load;
    unsigned max_collisions;
} hash_table_t ;

DLL_LOCAL void hash_table_create( hash_table_t * table, unsigned initial_size,
       int (*is_equal)(const void * a, const void * b),
       size_t (*hash)(const void * x) ) ;
DLL_LOCAL void hash_table_destroy( hash_table_t * table );
DLL_LOCAL void hash_table_clear( hash_table_t * table );
DLL_LOCAL void * hash_table_new_item( hash_table_t * table, void * key );
DLL_LOCAL void * hash_table_get_item( const hash_table_t * table, void * key );
DLL_LOCAL void * hash_table_delete_item( hash_table_t * table, void * key );
DLL_LOCAL unsigned hash_table_get_size( const hash_table_t * table );



#ifdef __cplusplus
}
#endif

#endif

