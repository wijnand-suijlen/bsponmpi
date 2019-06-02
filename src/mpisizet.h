#ifndef BSPONMPI_MPISIZET_H
#define BSPONMPI_MPISIZET_H

#include <mpi.h>

#ifdef __cplusplus
#include <cstddef>
#include <climits>
#else 
#include <stddef.h>
#include <limits.h>
#endif

#if   SIZE_MAX == UCHAR_MAX
  #define MY_MPI_SIZE_T MPI_UNSIGNED_CHAR
#elif SIZE_MAX == USHRT_MAX
  #define MY_MPI_SIZE_T MPI_UNSIGNED_SHORT
#elif SIZE_MAX == UINT_MAX
  #define MY_MPI_SIZE_T MPI_UNSIGNED
#elif SIZE_MAX == ULONG_MAX
  #define MY_MPI_SIZE_T MPI_UNSIGNED_LONG
#elif SIZE_MAX == ULLONG_MAX
  #define MY_MPI_SIZE_T MPI_UNSIGNED_LONG_LONG
#else
  #error "No MPI data type matches a size_t"
#endif

#endif
