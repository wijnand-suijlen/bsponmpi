#ifndef BSPONMPI_DLLEXPORT_H
#define BSPONMPI_DLLEXPORT_H


/* copied from https://gcc.gnu.org/wiki/Visibility */
#if defined _WIN32 || defined __CYGWIN__
  #ifdef BUILDING_DLL
    #ifdef __GNUC__
      #define DLL_PUBLIC __attribute__ ((dllexport))
    #elif defined __clang__
      #define DLL_PUBLIC __attribute__ ((dllexport))
    #else
      #define DLL_PUBLIC __declspec(dllexport)
    #endif
  #else
    #ifdef __GNUC__
      #define DLL_PUBLIC __attribute__ ((dllimport))
    #elif defined __clang__
      #define DLL_PUBLIC __attribute__ ((dllimport))
    #else
      #define DLL_PUBLIC __declspec(dllimport) 
    #endif
  #endif
  #define DLL_LOCAL
#else
  #if __GNUC__ >= 4 || defined __clang__
    #define DLL_PUBLIC __attribute__ ((visibility ("default")))
    #define DLL_LOCAL  __attribute__ ((visibility ("hidden")))
  #else
    #define DLL_PUBLIC
    #define DLL_LOCAL
  #endif
#endif

#endif
