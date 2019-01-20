cmake_minimum_required(VERSION 3.13)
project(BSPonMPI2) 

find_package(MPI)
find_package(doxygen)

set(CMAKE_C_FLAGS "-Wall -Wextra -g -pedantic -std=c89")
set(CMAKE_CXX_FLAGS "-Wall -Wextra -g -pedantic -std=c++11" )

enable_testing()

function( add_unit_test NAME  )
    set( exe "unit_test_${NAME}")
    add_executable( ${exe} src/${NAME}.t.cc ${ARGN} )
    set( testname "unit_test_${NAME}" )
    add_test( NAME "${testname}" COMMAND ${exe} )
endfunction()

function( add_mpi_unit_test NAME P )
    set( exe "unit_test_${NAME}")
    add_executable( ${exe} src/${NAME}.t.cc ${ARGN} )
    target_link_libraries(${exe} ${MPI_C_LIBRARIES} )
    set( testname "unit_test_${NAME}_p${P}" )
    add_test( NAME "${testname}" COMMAND ${MPIEXEC} -n ${P} ./${exe} )
endfunction()
 
function( add_mpi_test NAME NR STATUS P )
    set( exe "test_${NAME}_${NR}")
    add_executable( ${exe} test/${NAME}.${NR}.${STATUS}.c )
    target_link_libraries( ${exe} bsponmpi )
    set( testname "${NAME}_${NR}_p${P}" )
    add_test( NAME "${testname}" COMMAND ${MPIEXEC} -n ${P} ./${exe} )
    if (STATUS STREQUAL fail)
       set_tests_properties( "${testname}" PROPERTIES WILL_FAIL TRUE)
    endif()
endfunction()

# THe library itself
include_directories( ${MPI_INCLUDE_PATH} )
include_directories( include )
add_library( bsponmpi SHARED
        src/bsp.cc src/spmd.cc src/a2a.cc src/unbuf.cc src/rdma.cc
        src/bsmp.cc  
        )
target_link_libraries( bsponmpi ${MPI_C_LIBRARIES})


# Unit tests
add_mpi_unit_test( a2a 3 src/a2a.cc )
add_mpi_unit_test( unbuf 4 src/unbuf.cc )
add_unit_test( uintserialize )

# API tests
add_mpi_test( nprocs 1 ok 2 )
add_mpi_test( nprocs 2 ok 4 )
#add_c_test( abort 1 fail )
add_mpi_test( abort 2 fail 2 )
add_mpi_test( begin 1 fail 2 )
add_mpi_test( sync 1 ok 3 )
add_mpi_test( sync 2 fail 4 )
add_mpi_test( init 1 ok 2 )
add_mpi_test( init 2 ok 8 )
add_mpi_test( paper_example_reverse 1 ok 5 ) 
add_mpi_test( paper_example_put_array 1 ok 7 ) 
add_mpi_test( paper_example_bsp_sum 1 ok 3 )

# Documentation
set(DOXYGEN_BRIEF_MEMBER_DESC YES)
set(DOXYGEN_JAVADOC_AUTOBRIEF YES)
set(DOXYGEN_GENERATE_HTML YES)
set(DOXYGEN_GENERATE_MAN  YES)
set(DOXYGEN_MAN_LINKS YES) # create individual entries for each primitive
doxygen_add_docs( manual include/bsp.h )
