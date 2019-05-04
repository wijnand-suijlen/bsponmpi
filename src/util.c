#include "util.h"
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

int uint32_log2( uint32_t x )
{
    int a0, a1, a2, a3, a4;
    if (x == 0) return INT_MIN;
    a4 = (x > 0xffff)?16:0;
    x >>= a4;
    a3 = (x > 0xff)?8:0;
    x >>= a3;
    a2 = (x > 0xf)?4:0;
    x >>= a2;
    a1 = (x > 0x3)?2:0;
    x >>= a1;
    a0 = (x > 0x1)?1:0;
    return a0 + a1 + a2 + a3 + a4;
}
       
#ifdef UINT64_MAX
int uint64_log2( uint64_t x )
{
    int a0, a1, a2, a3, a4, a5;
    if (x == 0) return INT_MIN;
    a5 = (x > 0xffffffff)?32:0;
    x >>= a5;
    a4 = (x > 0xffff)?16:0;
    x >>= a4;
    a3 = (x > 0xff)?8:0;
    x >>= a3;
    a2 = (x > 0xf)?4:0;
    x >>= a2;
    a1 = (x > 0x3)?2:0;
    x >>= a1;
    a0 = (x > 0x1)?1:0;
    return a0 + a1 + a2 + a3 + a4 + a5;
}
#endif



unsigned int_log( unsigned base, unsigned x )
{
    unsigned q = base;
    unsigned y = 1, result = 0;
    unsigned powers[ CHAR_BIT*sizeof(unsigned) ];

    /*
     We use that we want to find a binary representation
     of the log base q of x
    
     Find n s.t. q^n = x
    
     Or when n is in base 2
    
     q^(a_0 2^0 + ... + a_m 2^m ) = x 
    
       <=>
    
     q^a_0 * (q^2)^a_1 * (q^4)^a_2 * ... * (q^(2^m))^a_m = x

     So, first find the largest k s.t. q^(2^k) < x
    */

    unsigned k = 0;
    unsigned power = 1;
    unsigned next_power = q;

    if ( q == 0 || q == 1 )
        return INT_MAX;

    if ( q == 2 )
        return int_log2( x );

    if ( x == 0 )
        return INT_MIN;

    while ( next_power <= x ) {
        power = next_power;
        powers[ k ] = power;
        if ( next_power > UINT_MAX/next_power)
            break;
        next_power = next_power * next_power;
        k += 1;
    }

    /* Now we can find the logarithm */
    while (k > 0) {
        k--;
        if (y * powers[k] <= x) {
            y *= powers[k];
            result += (1 << k);
        }
    }
        
    return result;
}


unsigned int_pow( unsigned base, unsigned n )
{
    unsigned result = 1;
    unsigned power = base;
    while (n > 0) {
        if ( n & 0x1 )
            result *= power;
            
        n >>= 1;
        power = power * power;
    }
        
    return result; 
}

double rand_next_uint64( uint64_t * next )
{
    /* According to Wikipedia this comes from Knuth */
    uint64_t a = 6364136223846793005;
    uint64_t c = 1442695040888963407;
    uint64_t x = (*next) * a + c; 
    *next =  x;
    return x / (1.0 + UINT64_MAX) ;
}

double rand_next_uint32( uint32_t * next )
{
    /* Same as the POSIX manual of rand() */
    uint32_t a = 1103515245;
    uint32_t c = 12345;
    uint32_t x = (*next) * a + c; 
    const uint32_t m = (1u<<31) - 1;
    *next = x;
    return (x & m )/ (1.0 + m);
}


double rand_next( size_t * next ) 
{
#ifdef UINT64_MAX
  #if SIZE_MAX == UINT64_MAX
    return rand_next_uint64( (uint64_t *) next );
  #elif SIZE_MAX == UINT32_MAX
    return rand_next_uint32( (uint32_t *) next );
  #endif
#else
  return rand_next_uint32( (uint32_t *) next );
#endif
}

#define hash( universal_hash_function, integer ) \
    ( ( (size_t) universal_hash_function.mult * (integer) + \
        universal_hash_function.add ) >> ( CHAR_BIT * sizeof(size_t) \
            - universal_hash_function.log2_buckets ) )

universal_hash_function_t new_universal_hash_function( size_t * seed,
       unsigned buckets )
{
    universal_hash_function_t result;
    result.mult = SIZE_MAX * rand_next( seed );
    result.add  = SIZE_MAX * rand_next( seed );
    result.log2_buckets = int_log2( buckets );
    return result;
}

struct hash_table_bucket {
    void * data;
    hash_table_bucket_t * next;
};

unsigned hash_table_get_size( const hash_table_t * table )
{
    return table->n_items;
}

void hash_table_create( hash_table_t * table, unsigned initial_size,
       int (*is_equal)(const void * a, const void * b),
       size_t (*hash)(const void * x) ) 
{
    unsigned i;
    if (initial_size < 2 ) initial_size = 2;

    table->n_items = 0;
    table->n_buckets = 1u << (int_log2( initial_size-1) + 1 );
    table->allocs  
        = calloc( 2*table->n_buckets, sizeof(hash_table_bucket_t));
    if ( table->allocs == NULL) {
        fprintf( stderr, "hash_table_create: Insufficient memory\n");
        abort();
    }

    table->buckets = table->allocs + 1;
    table->free_buckets  = table->allocs + table->n_buckets + 1;
    
    /* initialize free list */
    for ( i = 0; i < initial_size - 1; ++i )
        table->free_buckets[i].next = &table->free_buckets[i+1];
    table->free_buckets[initial_size-1].next = NULL;

    table->is_equal = is_equal;
    table->hash = hash;
    table->seed = 0;
    table->ghash 
        = new_universal_hash_function( &table->seed, table->n_buckets );

    /* some configurable constants */
    table->max_load = 1.0;
    table->max_collisions = 5;
}

void hash_table_destroy( hash_table_t * table )
{
    while( table->allocs )
    {
        hash_table_bucket_t * next = table->allocs->next;
        free(table->allocs);
        table->allocs = next;
    }
    table->n_items = 0;
    table->n_buckets = 0;
    table->allocs = NULL;
    table->buckets = NULL;
    table->free_buckets = NULL;
    table->is_equal = NULL;
    table->hash = NULL;
}

void hash_table_clear( hash_table_t * table )
{
    unsigned i;
    /* return all blocks to the free list */
    for ( i = 0; i < table->n_buckets; ++i ) {
        hash_table_bucket_t *c, * b = table->buckets[i].next;
        while (b) {
           c = b->next;
           b->next = table->free_buckets;
           b->data = NULL;
           table->free_buckets = b;
           b = c; 
        }
        table->buckets[i].next = NULL;
        table->buckets[i].data = NULL;
    }
    table->n_items = 0;
}


static void hash_table_increase_buckets( hash_table_t * table )
{
    unsigned i, j, new_n_buckets ;
    hash_table_bucket_t * new_table ;
    universal_hash_function_t new_ghash;
    assert( table->n_buckets > 0 );
    
    new_table = calloc(table->n_buckets * 2 + 1, sizeof(hash_table_bucket_t));
    if (!new_table) {
        fprintf(stderr, "hash_table_increase_buckets: "
                "insufficient memory\n");
        abort();
    }
    new_table->next = table->allocs;
    table->allocs = new_table;
    new_table += 1;

    new_n_buckets = 2 * table->n_buckets;
    new_ghash = new_universal_hash_function(&table->seed, new_n_buckets);
    
    /* rehash existing items*/
    for ( i = 0; i < table->n_buckets ; ++i ) {
        hash_table_bucket_t * b = table->buckets[i].next;
        while (b) {
            hash_table_bucket_t * c = b->next;
            j = hash( new_ghash, (* table->hash )( b->data ) );
            
            b->next = new_table[j].next;
            new_table[j].next = b;

            b = c;
        }
    }

    /* return old buckets to free list */
    for ( i = 0; i < table->n_buckets - 1; ++i ) {
        table->buckets[i].next = &table->buckets[i+1];
        table->buckets[i].data = NULL;
    } 
    table->buckets[table->n_buckets-1].next = table->free_buckets;
    table->buckets[table->n_buckets-1].data = NULL;
    table->free_buckets = table->buckets;

    /* commit changes */
    table->buckets = new_table;
    table->n_buckets = new_n_buckets;
    table->ghash = new_ghash;
}


static void hash_table_increase_free_buckets( hash_table_t * table )
{
    unsigned i, n = table->n_buckets;
    hash_table_bucket_t * new_buckets = 
        calloc( n, sizeof(hash_table_bucket_t) );
    if (!new_buckets) {
        fprintf(stderr, "hash_table_increase_free_buckets: "
                "insufficient memory\n");
        abort();
    }

    new_buckets->next = table->allocs;
    table->allocs = new_buckets;

    /* assign to free list */
    for ( i = 1; i < n - 1; ++i ) {
        new_buckets[i].next = &new_buckets[i+1];
        new_buckets[i].data = NULL;
    } 
    new_buckets[n-1].next = table->free_buckets;
    new_buckets[n-1].data = NULL;
    table->free_buckets = &new_buckets[1];
}

void * hash_table_new_item( hash_table_t * table, void * key )
{
    int try;
    double load = (table->n_items + 1.0)/table->n_buckets;
    hash_table_bucket_t ** b;
    hash_table_bucket_t * new_bucket;

    if (load > table->max_load)
        hash_table_increase_buckets( table );

    for ( try = 0; try < 5; ++try ) {
        unsigned collisions = 0;
        unsigned i = hash( table->ghash, (* table->hash )( key ) );
        b = & table->buckets[i].next;
        
        while (*b) {
            if ( (*table->is_equal)(key, (*b)->data) )
                return (*b)->data;
            collisions += 1;
            b = &(*b)->next;
        }

        if ( collisions > table->max_collisions ) 
            hash_table_increase_buckets( table );
        else
            goto found_item;
    }

    fprintf(stderr, "hash_table_new_item: Beaten by the odds; "
            "too many rehashes\n");
    abort();

found_item:

    if ( table->free_buckets == NULL )
        hash_table_increase_free_buckets( table );

    /* pop a free bucket from the free list */
    new_bucket = table->free_buckets;
    table->free_buckets = table->free_buckets->next;
    new_bucket->next = NULL;
    new_bucket->data = NULL;

    /* append it to the end of the bucket list */
    *b = new_bucket;
    (*b)->data = key;
    table->n_items += 1;
    return (*b)->data;
}

void * hash_table_get_item( const hash_table_t * table, void * key )
{
    unsigned i = hash( table->ghash, (* table->hash )( key ) );
    const hash_table_bucket_t *b = table->buckets[i].next;

    while (b) {
        if ( (*table->is_equal)(key, b->data) )
            return b->data;
        b = b->next;
    }

    return NULL;
}

void * hash_table_delete_item( hash_table_t * table, void * key )
{
    unsigned i = hash( table->ghash, (* table->hash )( key ) );
    hash_table_bucket_t **b = & table->buckets[i].next;
   
    while (*b) {
        if ( (*table->is_equal)(key, (*b)->data) )
        {
            hash_table_bucket_t * c = *b;
            *b = (*b)->next;
            table->n_items -= 1;
            c->next = table->free_buckets;
            table->free_buckets = c;
            return c->data;
        }
        b = &(*b)->next;
    }
    return NULL;
}




