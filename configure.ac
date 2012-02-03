#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.59])
AC_INIT([OPM Polymer Module], [0.1], [atgeirr@sintef.no],dnl
        [opmpolymer], [https://public.ict.sintef.no/opm/hg/opm-polymer])

AM_INIT_AUTOMAKE([-Wall -Werror foreign subdir-objects])

m4_ifdef([AM_SILENT_RULES], [AM_SILENT_RULES([yes])])
m4_ifdef([LT_INIT], [LT_INIT])

AC_CONFIG_MACRO_DIR([m4])
AC_CONFIG_SRCDIR([opm/polymer/HelloPolymer.hpp])
AC_CONFIG_HEADERS([config.h])

# Checks for programs.
AC_PROG_CC
AC_PROG_CXX
AM_PROG_CC_C_O
AC_PROG_RANLIB
m4_ifdef([LT_INIT], [], [AC_PROG_LIBTOOL])

# AX_LAPACK requires a working F77 compiler or, rather, its runtime
# support libraries.
AC_PROG_F77

# F77 name mangling
AC_F77_WRAPPERS

# Checks for libraries.
AX_LAPACK
AX_BOOST_BASE([1.37])
AX_BOOST_SYSTEM
AX_BOOST_DATE_TIME
AX_BOOST_FILESYSTEM
AX_BOOST_UNIT_TEST_FRAMEWORK

AC_SEARCH_LIBS([process_grdecl], [opmcore],
               [opmcore_lib=yes],
               [opmcore_lib=no],
               [$BOOST_LDFLAGS $BOOST_SYSTEM_LIB
                $LAPACK_LIBS $BLAS_LIBS $LIBS $FLIBS])


# Checks for header files.
AC_CHECK_HEADERS([opm/core/grid.h], [],
[AC_MSG_ERROR(
  [cannot proceed without 'opm-core' header files.
   *** Update CPPFLAGS and try again. ***])])

AC_CHECK_HEADERS([suitesparse/umfpack.h],
                 [umfpack_header=yes],
                 [umfpack_header=no])


# Checks for typedefs, structures, and compiler characteristics.

# Checks for library functions.
AC_SEARCH_LIBS([umfpack_dl_solve], [umfpack],
               [umfpack_lib=yes], [umfpack_lib=no])

AC_SEARCH_LIBS([preprocess], [opmcore], [],
[AC_MSG_ERROR(
  [cannot proceed without CPGPreprocess library.
   *** Update LDFLAGS and try again. ***])],
[$BOOST_SYSTEM_LIB $LAPACK_LIBS $BLAS_LIBS])

AC_SEARCH_LIBS([csrmatrix_elm_index], [opmcore], [],
[AC_MSG_ERROR(
  [cannot proceed without opmpressure library.
   *** Update LDFLAGS and try again. ***])],
[$BOOST_SYSTEM_LIB $LAPACK_LIBS $BLAS_LIBS])

AM_CONDITIONAL([UMFPACK],
  [test "x$umfpack_header" != "xno" -a "x$umfpack_lib" != "xno"])

m4_ifdef([AM_COND_IF],
[AM_COND_IF([UMFPACK], [],
 [AC_MSG_NOTICE([Found no working installation of UMFPACK.
                 UMFPACK support is disabled.])])
])

AC_CONFIG_FILES([
  Makefile
  examples/Makefile
])

AC_OUTPUT
