#include <bsp.h>



int main(int argc, char ** argv) 
{
    (void) argc;
    (void) argv;
    bsp_abort("Dit is een test abort. Twee = %d\n", 2);


    return 0;
}
