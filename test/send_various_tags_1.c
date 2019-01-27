#include <bsp.h>
#include "test.h"

#include <stdint.h>

TEST( send_various_tags_1, success() )
{
    uint16_t t1;
    int32_t t2;
    double t3;
    bsp_pid_t s, P, t;
    bsp_size_t n, m;
    int i, j;
    char * cs = NULL;

    bsp_begin( bsp_nprocs() );

    s = bsp_pid();
    P = bsp_nprocs();

    for ( i = 0; i < s; ++i ) {
        bsp_send( i, NULL, NULL, 0 );
    }

    t = sizeof(t1);
    bsp_set_tagsize( &t );
    EXPECT_EQ("%d", t, 0 );
    bsp_sync();

    bsp_qsize( &n, &m );
    EXPECT_EQ( "%d", n, P-s-1 );
    EXPECT_EQ( "%d", m, 0 );

    t = sizeof(t2);
    bsp_set_tagsize( &t );
    EXPECT_EQ("%d", t, 2 );

    for ( i = 0; i < s; ++i ) {
        uint16_t x = s;
        bsp_send( i, &x, NULL, 0 );
    }

    for ( i = 0; i < P-s-1; ++i) {
        int status; 
        bsp_get_tag( &status, NULL );
        EXPECT_EQ( "%d", status, 0 );
        bsp_move( NULL, 0);
    }
    { int status; 
      bsp_get_tag( &status, NULL );
      EXPECT_EQ("%d", status, -1);
    }

    bsp_sync();
           
    bsp_qsize( &n, &m );
    EXPECT_EQ( "%d", n, P-s-1 );
    EXPECT_EQ( "%d", m, 0 );

    t = sizeof(t3);
    bsp_set_tagsize( &t );
    EXPECT_EQ("%d", t, 4 );

    for ( i = 0; i < P-s-1; ++i) {
        int status; uint16_t x;
        bsp_get_tag( &status, &x );
        EXPECT_EQ( "%d", status, 0 );
        EXPECT_EQ( "%d", (int) x, P-i-1 );
        bsp_move( NULL, 0 );
    }
    { int status; uint16_t x = 0;
      bsp_get_tag( &status, &x);
      EXPECT_EQ("%d", status, -1);
      EXPECT_EQ("%d", (int) x, 0);
    }

    cs = calloc( P-s, sizeof(char) );
    for ( i = 0; i < P-s; ++i )
        cs[i] = s + i;

    for ( i = 0; i < s; ++i ) {
        uint32_t x = s;
        bsp_send( i, &x, cs, P-s);
    }
    free(cs); cs = NULL;
    bsp_sync();

    bsp_qsize( &n, &m );
    EXPECT_EQ( "%d", n, P-s-1 );

    /* payload size
       s0: receives P-1 + P-2 + ... + 1 = P(P-1)/2
       s1: receives P-2 + P-3 + ... + 1 = (P-1)(P-2)/2
     */
    EXPECT_EQ( "%d", m, (P-s)*(P-s-1)/2 );

    t = 0;
    bsp_set_tagsize( &t );
    EXPECT_EQ("%d", t, 8 );

    for ( i = 0; i < P-s-1; ++i) {
        void *tp, *pp;
        int status; uint32_t x, *px;
        bsp_get_tag( &status, &x );
        EXPECT_EQ( "%d", status, i+1 );
        EXPECT_EQ( "%d", (int) x, P-i-1 );
        status = bsp_hpmove( &tp, &pp );
        px = tp; cs = pp;
        EXPECT_EQ( "%d", status, i+1 );
        EXPECT_EQ( "%d", (int) *px, P-i-1); 
        for ( j = 0; j < i+1 ; ++j ) {
            EXPECT_EQ( "%d", (int) cs[j],  P-i+j-1 );
        }
    }
    { int status; uint32_t x;
      bsp_get_tag( &status, &x);
      EXPECT_EQ("%d", status, -1);
    }

    bsp_end();
}
