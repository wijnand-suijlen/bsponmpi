cmake_minimum_required(VERSION 3.13)
project(BSPonMPI) 

set(CPACK_PACKAGE_NAME          "BSPonMPI")
set(CPACK_PACKAGE_VERSION_MAJOR "1" )
set(CPACK_PACKAGE_VERSION_MINOR "0" )
set(CPACK_PACKAGE_VERSION_PATCH "0")
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY 
        "BSPonMPI implements the BSPlib parallel programming framework on MPI")
set(CPACK_RESOURCE_FILE_LICENSE "${PROJECT_SOURCE_DIR}/LICENSE")
set(CPACK_STRIP_FILES TRUE)
set(CPACK_PACKAGE_CONTACT "wijnandsuijlen@gmail.com")
set(CPACK_RPM_PACKAGE_LICENSE "MIT License")
set(CPACK_RPM_PACKAGE_REQUIRES "mpi-devel")

include(CPack)

set(SOVERSION "${CPACK_PACKAGE_VERSION_MAJOR}.${CPACK_PACKAGE_VERSION_MINOR}.${CPACK_PACKAGE_VERSION_PATCH}")

find_package(Doxygen)
find_package(MPI)

if (NOT MPI_FOUND)
    message(SEND_ERROR "MPI was not found. Without it, BSPonMPI can't be built nor function")
endif()
string(REPLACE mpiexec mpirun MPIRUN_EXECUTABLE ${MPIEXEC_EXECUTABLE})

# do some compile dependent set-up
if (CMAKE_COMPILER_IS_GNUCC)
set(CMAKE_C_FLAGS_DEBUG 
        "${CMAKE_C_FLAGS_DEBUG} -Wall -Wextra -g  -pedantic -std=c89")

set(CMAKE_C_FLAGS_RELEASE 
        "${CMAKE_C_FLAGS_RELEASE} -Wall -Wextra -pedantic -std=c89")
endif()

try_compile(HAS_TR1_UNORDERED_MAP 
        ${CMAKE_BINARY_DIR}/cmake 
        ${CMAKE_SOURCE_DIR}/cmake/has_tr1_unordered_map.cc)

try_compile(HAS_CXX11_UNORDERED_MAP 
        ${CMAKE_BINARY_DIR}/cmake
        ${CMAKE_SOURCE_DIR}/cmake/has_cxx11_unordered_map.cc)

if (HAS_TR1_UNORDERED_MAP)
    add_definitions(-DHAS_TR1_UNORDERED_MAP)
endif()

if (HAS_CXX11_UNORDERED_MAP)
    add_definitions(-DHAS_CXX11_UNORDERED_MAP)
endif()

if (CMAKE_COMPILER_IS_GNUCXX)
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wall -Wextra -g -pedantic")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -Wall -Wextra -pedantic")
endif()

# Testing

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

install( TARGETS bsponmpi
         RUNTIME DESTINATION bin
         LIBRARY DESTINATION lib
         ARCHIVE DESTINATION lib
       )

# Compiler and bsprun frontends
set(BSPONMPI_LIB "-L${CMAKE_INSTALL_PREFIX}/lib -lbsponmpi")
set(BSPONMPI_INCLUDE "${CMAKE_INSTALL_PREFIX}/include")

configure_file( bspcc.in bspcc @ONLY )
configure_file( bspcc.in bspcxx @ONLY )
configure_file( bsprun.in bsprun @ONLY )
install( PROGRAMS ${CMAKE_BINARY_DIR}/bspcc
                  ${CMAKE_BINARY_DIR}/bspcxx
                  ${CMAKE_BINARY_DIR}/bsprun
         DESTINATION bin
       )


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
add_mpi_test( paper_example_all_gather_sparse_vec 1 ok 10 )

# Documentation
if (DOXYGEN_FOUND)
    message( STATUS "Doxygen was not found - Documentation will not be generated")
    set(DOXYGEN_BRIEF_MEMBER_DESC YES)
    set(DOXYGEN_JAVADOC_AUTOBRIEF YES)
    set(DOXYGEN_MACRO_EXPANSION YES)
    set(DOXYGEN_EXPAND_ONLY_PREDEF YES)
    set(DOXYGEN_EXPAND_AS_DEFINED "DLL_LOCAL;DLL_PUBLIC")
    set(DOXYGEN_GENERATE_HTML YES)
    set(DOXYGEN_GENERATE_MAN  YES)
    set(DOXYGEN_MAN_LINKS YES) # create individual entries for each primitive
    doxygen_add_docs( manual include/bsp.h )
    add_custom_target(docs ALL DEPENDS manual)

    install(FILES include/bsp.h 
            DESTINATION include )
    install(DIRECTORY ${CMAKE_BINARY_DIR}/man
            DESTINATION share)
    install(DIRECTORY ${CMAKE_BINARY_DIR}/html 
            DESTINATION share/doc/bsponmpi)
    install(FILES LICENSE README
            DESTINATION share/doc/bsponmpi)
    
endif()


