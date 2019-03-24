#ifndef BSPONMPI_TEST_H
#define BSPONMPI_TEST_H

#include <stdio.h>
#include <stdlib.h>

#include "dllexport.h"

#define EXPECT_EQ( format, actual, expected ) \
    do { if ( (actual) != (expected) ) { \
        printf("TEST FAILED: Expected %s to be equal to %s = " \
                format " but got " format "\n"\
                "   " __FILE__ ":%d\n", \
               #actual, #expected, (expected), (actual), __LINE__ );  \
        abort(); } \
    } while(0)

#define EXPECT_OP( format, a, op, b) \
    do { if ( !((a) op (b)) ) { \
        printf("TEST FAILED: Expected %s %s %s but " \
                format " " #op " " format "\n"\
                "   " __FILE__ ":%d\n", \
                #a, #op, #b, (a), (b),  __LINE__ );  \
        abort(); } \
    } while(0)

#define EXPECT_IMPLIES( a, b) \
    do { if ( !(!(a) || (b)) ) { \
        printf("TEST FAILED: Expected %s implies %s but " \
                 "%s = %d and %s = %d\n"\
                "   " __FILE__ ":%d\n", \
                #a, #b, #a, (int) (a), #b, (int) (b),  __LINE__ );  \
        abort(); } \
    } while(0)


#define EXPECT( cond ) \
    do { if ( !(cond) ) { \
        printf("TEST FAILED: Expected " #cond " to be true\n" \
                "   " __FILE__ ":%d\n", \
                __LINE__ );  \
        abort(); }\
    } while(0)


DLL_PUBLIC void bsp_intern_expect_success(void);
DLL_PUBLIC void bsp_intern_expect_abort( const char * message );

#define TEST( name, result ) \
        void name(void); \
        int main( int argc, char ** argv ) \
        { (void) argc; (void) argv; \
            bsp_intern_expect_ ## result; \
            name(); \
            return 0; \
        }\
        void name() 



#endif
