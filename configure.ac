AC_PREREQ([2.67])
AC_INIT([minidjvu-mod],[0.9m01],[trufanovan@gmail.com],,[https://github.com/trufanov-nok/minidjvu_mod/])
AC_CONFIG_AUX_DIR(config)
AC_CONFIG_HEADER(config.h:config/config.h.in)
AC_CONFIG_MACRO_DIR([m4])
AC_PREFIX_DEFAULT(/usr/local)

AM_INIT_AUTOMAKE([subdir-objects foreign -Wall])

AC_LANG([C++])

AM_GNU_GETTEXT([external])
AM_GNU_GETTEXT_VERSION([0.18.1])

# Checks for programs.
AC_PROG_CXX
AC_PROG_AWK
AC_PROG_CC
AC_PROG_CPP
AC_PROG_INSTALL
AC_PROG_LN_S
AC_PROG_MAKE_SET
AC_PATH_PROG(GZIP, gzip)
AC_PATH_PROG(RM, rm)
AM_PROG_AR

LT_INIT

PKG_INSTALLDIR

# Checks for libraries.

AC_CHECK_LIB(z, inflate)
AC_CHECK_LIB(jpeg, jpeg_destroy_decompress)
AC_CHECK_LIB(tiff, TIFFOpen)
AC_CHECK_LIB(jemalloc,malloc)
# Check for OpenMP
AC_OPENMP


# Checks for header files.
AC_CHECK_HEADERS([libintl.h locale.h stdint.h stdlib.h string.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_SYS_LARGEFILE
AC_HEADER_STDBOOL
AC_C_INLINE
AC_TYPE_SIZE_T
AC_TYPE_INT16_T
AC_TYPE_INT32_T
AC_TYPE_SIZE_T
AC_TYPE_UINT16_T
AC_TYPE_UINT32_T

# Checks for library functions.
#AC_FUNC_MALLOC
#AC_FUNC_REALLOC
AC_CHECK_FUNCS([memset pow setlocale strcspn strrchr strerror])


###############################################
# start # checks required for settings-reader #
###############################################

AC_CHECK_DECL([_WIN32],[have_os_win32=yes],[have_os_win32=no])
AM_CONDITIONAL([HAVE_OS_WIN32], [test "x${have_os_win32}" = "xyes"])

AC_CHECK_DECL([__APPLE__],[have_os_apple=yes],[have_os_apple=no])
AM_CONDITIONAL([HAVE_OS_APPLE], [test "x${have_os_apple}" = "xyes"])


# C++ 
AC_LANG(C++)
AC_CXX_BOOL
AC_CXX_EXCEPTIONS
AC_CXX_TYPENAME
AC_CXX_STDINCLUDES
AC_CXX_NAMESPACES
AC_CXX_MEMBER_TEMPLATES
AC_CXX_INTEL_ATOMIC_BUILTINS
AC_CXX_GCCTLS

AC_CHECK_HEADERS(wchar.h wctype.h inttypes.h)
AC_CHECK_TYPES(wchar_t)
AC_CHECK_TYPES(mbstate_t,,,[#include "wchar.h"])

# Search for PTHREADS (when not on windows)
have_pthread=no
AM_COND_IF([HAVE_OS_WIN32],,[AC_PATH_PTHREAD([have_pthread=yes])])

# ----------------------------------------
# Stuff added to config.h
# ----------------------------------------

# Fence
AH_TOP([
#ifndef CONFIG_H
#define CONFIG_H
/* config.h: begin */
])

# L18N Macros
AH_BOTTOM([

/* - Miscellaneous */
#define AUTOCONF 1
#if defined(__CYGWIN32__) || !defined(_WIN32)
# define UNIX 1
#endif

/* - WCHAR etc.*/
#if ! defined(HAVE_WCHAR_T)
#define HAS_WCHAR 0
#define HAS_WCTYPE 0
#define HAS_MBSTATE 0
#else
#define HAS_WCHAR 1
#if defined(HAVE_WCTYPE_H) && defined(HAVE_ISWSPACE)
#define HAS_WCTYPE 1
#endif
#if defined(HAVE_MBSTATE_T) && defined(HAVE_WCRTOMB)
#define HAS_MBSTATE 1
#endif
#endif

/* config.h: end */
#endif
])

#############################################
# end # checks required for settings-reader #
#############################################

# Trailer
AC_CONFIG_FILES([Makefile
                 doc/Makefile
                 doc/ru/Makefile
                 include/Makefile
		 po/Makefile.in])
AC_OUTPUT
