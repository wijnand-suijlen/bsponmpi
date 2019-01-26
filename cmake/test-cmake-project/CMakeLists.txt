cmake_minimum_required(VERSION 2.8)
project(TestProject)

find_program(BSPRUN bsprun)

enable_testing()
add_executable(hello hello.c)
add_test(NAME test_hello 
        COMMAND ${BSPRUN} -n 5 $<TARGET_FILE:hello> )
