#include "util.h"

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



