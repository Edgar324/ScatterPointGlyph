/*
 * MATLAB Compiler: 5.1 (R2014a)
 * Date: Sun Jan 17 12:56:14 2016
 * Arguments: "-B" "macro_default" "-W" "lib:libNcut" "-T" "link:lib" "ncutW.m" 
 */

#ifndef __libNcut_h
#define __libNcut_h 1

#if defined(__cplusplus) && !defined(mclmcrrt_h) && defined(__linux__)
#  pragma implementation "mclmcrrt.h"
#endif
#include "mclmcrrt.h"
#ifdef __cplusplus
extern "C" {
#endif

#if defined(__SUNPRO_CC)
/* Solaris shared libraries use __global, rather than mapfiles
 * to define the API exported from a shared library. __global is
 * only necessary when building the library -- files including
 * this header file to use the library do not need the __global
 * declaration; hence the EXPORTING_<library> logic.
 */

#ifdef EXPORTING_libNcut
#define PUBLIC_libNcut_C_API __global
#else
#define PUBLIC_libNcut_C_API /* No import statement needed. */
#endif

#define LIB_libNcut_C_API PUBLIC_libNcut_C_API

#elif defined(_HPUX_SOURCE)

#ifdef EXPORTING_libNcut
#define PUBLIC_libNcut_C_API __declspec(dllexport)
#else
#define PUBLIC_libNcut_C_API __declspec(dllimport)
#endif

#define LIB_libNcut_C_API PUBLIC_libNcut_C_API


#else

#define LIB_libNcut_C_API

#endif

/* This symbol is defined in shared libraries. Define it here
 * (to nothing) in case this isn't a shared library. 
 */
#ifndef LIB_libNcut_C_API 
#define LIB_libNcut_C_API /* No special import/export declaration */
#endif

extern LIB_libNcut_C_API 
bool MW_CALL_CONV libNcutInitializeWithHandlers(
       mclOutputHandlerFcn error_handler, 
       mclOutputHandlerFcn print_handler);

extern LIB_libNcut_C_API 
bool MW_CALL_CONV libNcutInitialize(void);

extern LIB_libNcut_C_API 
void MW_CALL_CONV libNcutTerminate(void);



extern LIB_libNcut_C_API 
void MW_CALL_CONV libNcutPrintStackTrace(void);

extern LIB_libNcut_C_API 
bool MW_CALL_CONV mlxNcutW(int nlhs, mxArray *plhs[], int nrhs, mxArray *prhs[]);



extern LIB_libNcut_C_API bool MW_CALL_CONV mlfNcutW(int nargout, mxArray** NcutDiscrete, mxArray** NcutEigenvectors, mxArray** NcutEigenvalues, mxArray* W, mxArray* nbcluster);

#ifdef __cplusplus
}
#endif
#endif
