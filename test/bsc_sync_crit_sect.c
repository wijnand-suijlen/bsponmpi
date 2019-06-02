#ifdef _MSC_VER
#define _CRT_SECURE_NO_WARNINGS
#endif
#include <bsp.h>
#include <bsc.h>
#include <stdio.h>
#include "test.h"

TEST( bsc_sync_crit_sect, success() )
{
    int rc;
    FILE * file;
    char filename[] = "bsc_sync_crit_sect_test_file";
    bsp_begin( bsp_nprocs() );

    /* this test writes a file in the current working directory
     * and checks whether it is the only one have done so */
    bsc_sync( 2*bsp_pid() );

    file = fopen( filename, "r" );
    EXPECT_OP( "%p", (void *) file, ==, NULL);

    file = fopen( filename, "w" );
    EXPECT_OP( "%p", (void *) file, !=, NULL);

    fprintf(file,"TEST FILE written by %d/%d\n", bsp_pid(), bsp_nprocs());
    rc = fclose(file);
    EXPECT_EQ( "%d", rc, 0 );

    bsc_sync( 2*bsp_pid() + 1); 

    rc = remove( filename );
    EXPECT_EQ( "%d", rc, 0 );

    bsc_sync( 2*bsp_nprocs() );

    bsp_end();
}
