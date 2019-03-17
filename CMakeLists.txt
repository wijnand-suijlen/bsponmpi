cmake_minimum_required(VERSION 2.8.11)
project(BSPonMPI) 

set(CPACK_PACKAGE_NAME          "bsponmpi")
set(CPACK_PACKAGE_VERSION_MAJOR "1")
set(CPACK_PACKAGE_VERSION_MINOR "0")
set(CPACK_PACKAGE_VERSION_PATCH "1")
set(CPACK_PACKAGE_VERSION "${CPACK_PACKAGE_VERSION_MAJOR}.${CPACK_PACKAGE_VERSION_MINOR}.${CPACK_PACKAGE_VERSION_PATCH}" )
set(CPACK_PACKAGE_FILE_NAME     "${CPACK_PACKAGE_NAME}-${CPACK_PACKAGE_VERSION}-${CMAKE_SYSTEM_NAME}-${CMAKE_SYSTEM_PROCESSOR}")
set(CPACK_RESOURCE_FILE_LICENSE "${PROJECT_SOURCE_DIR}/LICENSE")
set(CPACK_STRIP_FILES TRUE)
set(CPACK_RPM_PACKAGE_RELEASE 1)
set(CPACK_PACKAGE_CONTACT "bsponmpi@suijlen.net")
set(CPACK_RPM_USER_BINARY_SPECFILE "${PROJECT_SOURCE_DIR}/bsponmpi.spec.in")
set(CPACK_RPM_EXCLUDE_FROM_AUTO_FILELIST "/usr;/usr/bin;/usr/lib;/usr/include;/usr/share;/usr/share/man;/usr/share/man/man1;/usr/share/man/man3;/usr/share/doc")
set(CPACK_RPM_SPEC_MORE_DEFINE "%define _rpmdir %_topdir/RPMS")
if (CMAKE_VERSION VERSION_LESS 3)
# Workaround for older CMake versions
    set(CPACK_RPM_SPEC_MORE_DEFINE "%define _rpmdir %_topdir")
else()
# Worarounds for newer CMake versions
    cmake_policy(SET CMP0042 NEW)
endif()

find_program(RPMBUILD rpmbuild)
find_program(DPKG     dpkg-deb)
set(generators)
if (RPMBUILD)
    list(APPEND generators "RPM")
    message(STATUS "Can build RPM package. Do: make package")
endif()
if(DPKG)
    list(APPEND generators "DEB")
    message(STATUS "Can build DEB package. Do: make package")
endif()
set(CPACK_GENERATOR "${generators}" )

include(CPack)
include(CheckSymbolExists)

string(TIMESTAMP BUILD_DATE)
set(SOVERSION "${CPACK_PACKAGE_VERSION}")

set(DOXYGEN_SKIP_DOT TRUE)
find_package(Doxygen)
find_package(MPI)

if (NOT MPI_FOUND)
    message(SEND_ERROR "MPI was not found. Without it, BSPonMPI can't be built nor function")
endif()
if (NOT MPIEXEC_EXECUTABLE) 
    # In previous versions of CMake, this used to be MPIEXEC
    set(MPIEXEC_EXECUTABLE "${MPIEXEC}")
endif()
string(REPLACE mpiexec mpirun mpirun_exec ${MPIEXEC_EXECUTABLE})
# fallback to mpiexec if mpirun does not exist
if (NOT EXISTS "${mpirun_exec}")
  set(mpirun_exec "${MPIEXEC_EXECUTABLE}")
endif()
set( MPIRUN_EXECUTABLE "${mpirun_exec}" CACHE FILEPATH "The path to the mpirun launcher (or equivalent)" )
set( MPIRUN_NPROCS_PARAM "${MPIEXEC_NUMPROC_FLAG};" CACHE STRING "Parameter to mpirun for setting the number of processes. Use a semicolon ; to denotes that the actual number of processes is the next separate command line parameter")
string(REPLACE ";" " " MPIRUN_NPROCS_PARAM_SCRIPT "${MPIRUN_NPROCS_PARAM}")

# Disable C++ bindings in OpenMPI
add_definitions(-DOMPI_SKIP_MPICXX)
# Disable C++ bindings in Platform MPI
add_definitions(-D_MPICC_H)

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

check_symbol_exists( clock_gettime "time.h" HAS_CLOCK_GETTIME )
if (HAS_CLOCK_GETTIME)
    add_definitions(-DHAS_CLOCK_GETTIME)
endif()

check_symbol_exists( _dupenv_s "stdlib.h" HAS_DUPENV_S_WIN32)
if (HAS_DUPENV_S_WIN32 AND WIN32)
    add_definitions(-DHAS_DUPENV_S_WIN32)
endif()

if (CMAKE_COMPILER_IS_GNUCXX)
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wall -Wextra -g -pedantic")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -Wall -Wextra -pedantic")
endif()

# On Apple and Win32 we want to build a shared lib by default. On other
# machines we prefer static
set(prefer_shared_libs OFF)
if (APPLE OR WIN32)
  set(prefer_shared_libs ON)
endif()
option(BUILD_SHARED_LIBS "Build BSPonMPI as a shared library; When you get strange linking problems, you might want to toggle this option" "${prefer_shared_libs}" ) 

# Macro's to facilitate Testing

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
    add_test( NAME "${testname}" COMMAND ${MPIRUN_EXECUTABLE} ${MPIRUN_NPROCS_PARAM}${P} $<TARGET_FILE:${exe}> )
endfunction()
 
function( add_mpi_test NAME P )
    foreach ( flavour msg rma rma_prof )
        set( exe "test_${NAME}_${flavour}")
        add_executable( ${exe} test/${NAME}.c )
        set_target_properties( ${exe} PROPERTIES 
                INCLUDE_DIRECTORIES 
                "${PROJECT_SOURCE_DIR}/include;${PROJECT_SOURCE_DIR}/test;${PROJECT_SOURCE_DIR}/src;${MPI_INCLUDE_PATH}" )
        target_link_libraries( ${exe} "bsponmpi_${flavour}" )
        set( testname "${NAME}_${flavour}_p${P}" )
        add_test( NAME "${testname}" COMMAND ${MPIRUN_EXECUTABLE} ${MPIRUN_NPROCS_PARAM}${P} $<TARGET_FILE:${exe}> )
    endforeach()
endfunction()

# The library itself
include_directories( ${MPI_INCLUDE_PATH} )
include_directories( include )
set(BSPONMPI_SOURCES
        src/bsp.cc src/spmd.cc src/a2a.cc src/unbuf.cc src/rdma.cc
        src/bsmp.cc src/tictoc.cc
   )

# library using message passing
function ( make_bsponmpi name )
    add_library( ${name} ${BSPONMPI_SOURCES})
    target_link_libraries( ${name} ${MPI_C_LIBRARIES})
    target_compile_definitions(${name} PRIVATE BUILDING_DLL ${ARGN})
    install( TARGETS ${name}
             RUNTIME DESTINATION bin
             LIBRARY DESTINATION lib
             ARCHIVE DESTINATION lib
           )
endfunction()

make_bsponmpi( bsponmpi_msg )
make_bsponmpi( bsponmpi_rma USE_ONESIDED )
make_bsponmpi( bsponmpi_msg_prof PROFILE )
make_bsponmpi( bsponmpi_rma_prof USE_ONESIDED PROFILE )

# Installation
install(FILES include/bsp.h 
         DESTINATION include )

# Compiler and bsprun frontends
set(BSPONMPI_LIB_DIR "${CMAKE_INSTALL_PREFIX}/lib")
set(BSPONMPI_INCLUDE "${CMAKE_INSTALL_PREFIX}/include")
find_program(BASH bash)

if (BASH)

configure_file( bspcc.in bspcc @ONLY )
configure_file( bspcc.in bspcxx @ONLY )
configure_file( bsprun.in bsprun @ONLY )
install( PROGRAMS ${CMAKE_BINARY_DIR}/bspcc
                  ${CMAKE_BINARY_DIR}/bspcxx
                  ${CMAKE_BINARY_DIR}/bsprun
         DESTINATION bin
       )
install( FILES man/man1/bsprun.1 man/man1/bspcc.1
         DESTINATION share/man/man1/ )

else()
  message(STATUS "Could not find bash; bspcc and bsprun frontends will not be installed" )
endif()


# Unit tests
add_mpi_unit_test( a2a 3 src/a2a.cc )
add_mpi_unit_test( a2a_rma 4 src/a2a.cc )
target_compile_definitions( unit_test_a2a_rma PRIVATE USE_ONESIDED)

add_mpi_unit_test( unbuf 4 src/unbuf.cc )
add_unit_test( uintserialize )

# API tests
add_mpi_test( abort_1 1 )
add_mpi_test( abort_2 2 )
add_mpi_test( begin_1 2 )
add_mpi_test( init_1 2 )
add_mpi_test( init_2 8 )
add_mpi_test( nprocs_1 2 )
add_mpi_test( nprocs_2 4 )
add_mpi_test( nprocs_3 5 )
add_mpi_test( paper_example_reverse 4 ) 
add_mpi_test( paper_example_put_array 3 ) 
add_mpi_test( paper_example_bsp_sum 3 )
add_mpi_test( paper_example_all_gather_sparse_vec 3 )
add_mpi_test( sync_1 3 )
add_mpi_test( sync_2 4 )
add_mpi_test( sync_3 3 )
add_mpi_test( push_reg_1 1 )
add_mpi_test( push_reg_2 2 )
add_mpi_test( push_reg_3 3 )
add_mpi_test( push_reg_4 2 )
add_mpi_test( push_reg_5 4 )
add_mpi_test( push_reg_6 3 )
add_mpi_test( push_reg_7 3 )
add_mpi_test( push_reg_8 2 )
add_mpi_test( push_reg_9 4 )
add_mpi_test( push_reg_10 4 )
add_mpi_test( push_reg_11 3 )
add_mpi_test( push_reg_12 2 )
add_mpi_test( push_reg_13 4 )
add_mpi_test( push_reg_14 2 )
add_mpi_test( push_reg_15 3 )
add_mpi_test( push_reg_16 3 )
add_mpi_test( push_reg_17 2 )
add_mpi_test( push_reg_18 4 )    
add_mpi_test( push_reg_19 2 )    
add_mpi_test( pop_reg_1 3 )
add_mpi_test( pop_reg_2 4 )
add_mpi_test( pop_reg_3 4 )
add_mpi_test( pop_reg_4 2 )
add_mpi_test( pop_reg_5 3 )
add_mpi_test( pop_reg_6 4 )
add_mpi_test( pop_reg_7 3 )    

add_mpi_test( put_1 3 )
add_mpi_test( put_2 4 )
add_mpi_test( put_3 4 )
add_mpi_test( put_4 3 )
add_mpi_test( put_5 4 )
add_mpi_test( put_6 4 )
add_mpi_test( put_7 1 )

add_mpi_test( get_1 3 )
add_mpi_test( get_2 4 )
add_mpi_test( get_3 4 )
add_mpi_test( get_4 3 )
add_mpi_test( get_5 4 )
add_mpi_test( get_6 4 )

add_mpi_test( putget_1 4 )
add_mpi_test( putget_2 5 )

add_mpi_test( hpput_1 3 )
add_mpi_test( hpget_1 3 )

add_mpi_test( set_tagsize_1 2 )
add_mpi_test( set_tagsize_2 3 )
add_mpi_test( set_tagsize_3 4 )
add_mpi_test( set_tagsize_4 3 )
add_mpi_test( set_tagsize_5 4 )
add_mpi_test( set_tagsize_spmd 2 )
add_mpi_test( set_tagsize_null 2 )
add_mpi_test( set_tagsize_negative_size 2 )

add_mpi_test( get_tag_spmd 2 )
add_mpi_test( get_tag_status_null 2 )
add_mpi_test( get_tag_tag_null 2 )

add_mpi_test( send_spmd 2 )
add_mpi_test( send_pid_negative 2 )
add_mpi_test( send_pid_greater_nprocs 2 )
add_mpi_test( send_size_negative 2 )

add_mpi_test( qsize_spmd 2 )
add_mpi_test( qsize_nmessages_null 2 )
add_mpi_test( qsize_accum_nbytes_null 2 )

add_mpi_test( default_tagsize_1 2 )

add_mpi_test( send_various_tags_1 4 )

add_mpi_test( hpmove_spmd 2 )
add_mpi_test( hpmove_payloadptr_null 2 )
add_mpi_test( hpmove_tagptr_null 2 )
add_mpi_test( hpmove_queue_empty 2 )

add_mpi_test( move_1 2 )
add_mpi_test( move_spmd 2 )
add_mpi_test( move_payload_null 2 )
add_mpi_test( move_queue_empty 2 )
add_mpi_test( move_nbytes_negative 2 )

add_mpi_test( pid_1 3)
add_mpi_test( pid_2 3)

option( ENABLE_LARGE_TESTS "Enable tests that require large amounts of memory and/or time" OFF)
if (ENABLE_LARGE_TESTS)
    message(STATUS "Tests that demand large amount of memory and/or time are enabled")
    add_mpi_test( hpput_one_int_max_size_msg 2 )

    #Test only relevant if bsp_size_t is un unassigned data type
    #add_mpi_test( send_size_endtag 2 )
endif()

if (BASH)
  configure_file( build_tests.sh.in build_tests.sh @ONLY)
  add_test( NAME build_tests COMMAND ${BASH} build_tests.sh )
endif()

# Documentation
if (NOT COMMAND doxygen_add_docs)
   
endif()
if (DOXYGEN_FOUND)
    configure_file(Doxyfile.in Doxyfile @ONLY)
    set(man3 ${PROJECT_BINARY_DIR}/man/man3)
    set(html ${PROJECT_BINARY_DIR}/html)
    add_custom_command(OUTPUT ${html}/index.html
            COMMAND ${DOXYGEN_EXECUTABLE} ${PROJECT_BINARY_DIR}/Doxyfile
            MAIN_DEPENDENCY ${PROJECT_SOURCE_DIR}/include/bsp.h
            WORKING_DIRECTORY ${PROJECT_BINARY_DIR}
            COMMENT "Invoking Doxygen to generate user documentation"
            )
    set(doc_output ${html}/index.html)

    install( FILES ${man3}/BSP_BSMP.3
                   ${man3}/BSP_DRMA.3
                   ${man3}/BSP_SPMD.3
                   ${man3}/BSP_TYPES.3
                   ${man3}/MCBSP.3
                   ${man3}/bsp_abort.3
                   ${man3}/bsp_begin.3
                   ${man3}/bsp_end.3
                   ${man3}/bsp_get.3
                   ${man3}/bsp_get_tag.3
                   ${man3}/bsp_hpget.3
                   ${man3}/bsp_hpmove.3
                   ${man3}/bsp_hpput.3
                   ${man3}/bsp_init.3
                   ${man3}/bsp_move.3
                   ${man3}/bsp_nprocs.3
                   ${man3}/bsp_pid.3
                   ${man3}/bsp_pid_t.3
                   ${man3}/bsp_pop_reg.3
                   ${man3}/bsp_push_reg.3
                   ${man3}/bsp_put.3
                   ${man3}/bsp_qsize.3
                   ${man3}/bsp_send.3
                   ${man3}/bsp_set_tagsize.3
                   ${man3}/bsp_size_t.3
                   ${man3}/bsp_sync.3
                   ${man3}/bsp_time.3
            DESTINATION share/man/man3/ )


    install(DIRECTORY ${CMAKE_BINARY_DIR}/html 
            DESTINATION share/doc/bsponmpi)
    install(FILES LICENSE README
            DESTINATION share/doc/bsponmpi)

    find_program(PDFLATEX pdflatex)
    find_program(MAKEINDEX makeindex)

    if (PDFLATEX AND MAKEINDEX)
       file(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/latex_test)
       file(COPY cmake/latex_test.tex 
               DESTINATION ${PROJECT_BINARY_DIR}/latex_test)
       execute_process(COMMAND ${PDFLATEX} latex_test
               WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/latex_test
               RESULT_VARIABLE LATEX_WORKS
               OUTPUT_FILE cmake-latex-test.log
               ERROR_FILE cmake-latex-test.log
               )


       if (LATEX_WORKS EQUAL 0)
       
           set(pdf "${PROJECT_BINARY_DIR}/latex/refman.pdf")
           add_custom_command( OUTPUT ${pdf}
                   COMMAND ${PDFLATEX} refman
                   COMMAND ${MAKEINDEX} refman.idx
                   COMMAND ${PDFLATEX} refman
                   COMMAND ${PDFLATEX} refman
                   COMMAND ${MAKEINDEX} refman.idx
                   COMMAND ${PDFLATEX} refman
                   WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/latex
                   MAIN_DEPENDENCY ${PROJECT_SOURCE_DIR}/include/bsp.h 
                   DEPENDS ${html}/index.html
                   COMMENT "Invoking Latex to generate PDF of manual"
                   )

           install(FILES ${pdf}
                  DESTINATION share/doc/bsponmpi
                  RENAME bsponmpi.pdf
                  )

           list(APPEND doc_output ${pdf})
        else()
            message(STATUS "The Latex installation misses some packages - A PDF of the reference documentation will not be generated")
        endif()
    else()
       message( STATUS "Latex and makeindex were not found - A PDF of the reference documentation will not be generated")
    endif()


    add_custom_target(docs ALL DEPENDS ${doc_output} )

else()
    message( STATUS "Doxygen was not found - Documentation will not be generated")
endif()


### post install check
configure_file( cmake/post-install-check.cmake.in post-install-check.cmake @ONLY)

if (BASH AND NOT MSVC)
  install(SCRIPT ${PROJECT_BINARY_DIR}/post-install-check.cmake)
endif()
