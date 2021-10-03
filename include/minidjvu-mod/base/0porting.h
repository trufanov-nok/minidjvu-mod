/*
 * 0porting.h - a portability header
 */

/*
 * "0porting" keeps several typedefs and MDJVU_FUNCTION/MDJVU_IMPLEMENT macros.
 * 
 * To compile a native Windows DLL, you need to put "__declspec(dllexport)" before
 * every exported function. To use a DLL, you have to put "__declspec(dllimport)"
 * before every function prototype you use. That's why we need macros here.
 * 
 * So, a function prototype
 * 
 *     MDJVU_FUNCTION void mdjvu_foo(void);
 * 
 * under Windows (when compiling the library) will expand into
 * 
 *     __declspec(dllexport) void mdjvu_foo(void);
 * 
 * and when compiling any other application - into
 * 
 *     __declspec(dllimport) void mdjvu_foo(void);
 * 
 * Under Linux this will lead to
 * 
 *     void mdjvu_foo(void);
 * 
 * in both cases.
 * 
 * (Also, under C++, there will be an `extern "C"' modifier).
 */


#ifndef MDJVU_USE_TIFFIO /* kluge not to typedef twice when using tiffio.h */
    #ifdef HAVE_STDINT_H
        #include <stdint.h>
        typedef int32_t int32;
        typedef uint32_t uint32;
        typedef int16_t int16;
        typedef uint16_t uint16;
    #else
        typedef int int32;
        typedef unsigned int uint32;
        typedef unsigned short uint16;
        typedef short int16;
    #endif
#endif

#ifndef HAVE_STDINT_H
    #define INT32_MAX 0x7FFFFFFF
#endif

#define MDJVU_INT32_FORMAT "%d"
#define MDJVU_INT16_FORMAT "%d"
#define MDJVU_UINT32_FORMAT "%u"
#define MDJVU_UINT16_FORMAT "%u"

/* MDJVU_FUNCTION and MDJVU_IMPLEMENT are prefixes of exported functions.
 * MDJVU_FUNCTION is for declarations, MDJVU_IMPLEMENT is for implementations.
 * So, it's like this:
 *
 *  // foo.h
 *  MDJVU_FUNCTION mdjvu_foo(void);
 *
 *  // foo.c
 *  MDJVU_IMPLEMENT mdjvu_foo(void)
 *  {
 *      ...
 *  }
 */

#ifdef _MSC_VER
#define _CRT_SECURE_NO_WARNINGS

#define __LITTLE_ENDIAN  1234
#define __BIG_ENDIAN  4321
#define __BYTE_ORDER __LITTLE_ENDIAN
#endif

#if defined(__cplusplus)
    #define MDJVU_C_EXPORT_PREFIX extern "C"
#else
    #define MDJVU_C_EXPORT_PREFIX
#endif

/* This Microsoft abomination of __declspec does not exist under mingw.
 */
#define MDJVU_SUPPRESS_DECLSPEC
#if (defined(windows) || defined(WIN32) || defined(_WIN32)) && !defined(MDJVU_SUPPRESS_DECLSPEC)
    #ifdef MINIDJVU_INCLUDED_FROM_INSIDE
        #define MDJVU_FUNCTION MDJVU_C_EXPORT_PREFIX __declspec(dllexport)
        #define MDJVU_IMPLEMENT __declspec(dllexport)
    #else
        #define MDJVU_FUNCTION MDJVU_C_EXPORT_PREFIX __declspec(dllimport)
        #define MDJVU_IMPLEMENT __declspec(dllimport)
    #endif
#else
    #define MDJVU_FUNCTION MDJVU_C_EXPORT_PREFIX
    #define MDJVU_IMPLEMENT
#endif

/* Convenience macros. */
#define MDJVU_MALLOC(T) ((T *) malloc(sizeof(T)))
#define MDJVU_MALLOCV(T,N) ((T *) malloc((N) * sizeof(T)))
#define MDJVU_CALLOC(T) ((T *) calloc(1, sizeof(T)))
#define MDJVU_CALLOCV(T,N) ((T *) calloc(N, sizeof(T)))
#define MDJVU_FREE(P) free(P)
#define MDJVU_FREEV(P) free(P)


/* Check that the portability typedefs work as expected.
 * If not, returns an error message.
 * Returns NULL if OK.
 */
MDJVU_FUNCTION const char *mdjvu_check_sanity(void);
