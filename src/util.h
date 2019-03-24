#ifndef BSPONMPI_UTIL_H
#define BSPONMPI_UTIL_H

#include <stdint.h>
#include <limits.h>

#define MAX( a, b ) ( (a) < (b) ? (b) : (a) )
#define MIN( a, b ) ( (a) > (b) ? (b) : (a) )

#ifdef __cplusplus
extern "C" {
#endif

#ifdef UINT64_MAX
int uint64_log2( uint64_t x );
#endif    
int uint32_log2( uint32_t x );

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

unsigned int_log( unsigned base, unsigned x );

unsigned int_pow( unsigned base, unsigned n );

#ifdef __cplusplus
}
#endif

#endif

