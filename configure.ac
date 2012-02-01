#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.59])
AC_INIT([OPM Polymer Module], [0.1], [atgeirr@sintef.no],dnl
        [opmpolymer], [https://public.ict.sintef.no/opm/hg/opm-polymer])

AM_INIT_AUTOMAKE([-Wall -Werror foreign subdir-objects])

m4_ifdef([AM_SILENT_RULES], [AM_SILENT_RULES([yes])])
m4_ifdef([LT_INIT], [LT_INIT])

AC_CONFIG_SRCDIR([opm/polymer/HelloPolymer.hpp])
AC_CONFIG_HEADERS([config.h])

# Checks for programs.
AC_PROG_CXX
AM_PROG_CC_C_O
m4_ifdef([LT_INIT], [], [AC_PROG_LIBTOOL])

# Checks for libraries.

# Checks for header files.

# Checks for typedefs, structures, and compiler characteristics.

# Checks for library functions.

AC_CONFIG_FILES([
  Makefile
  examples/Makefile
])

AC_OUTPUT
