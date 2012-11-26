AC_DEFUN([_ERT_SOURCE_TEXT],
[
AC_LANG_PROGRAM(
[[
#include <stddef.h>
#include <ecl_util.h>
]],dnl
[[
int sz;
sz = ecl_util_get_sizeof_ctype(ECL_INT_TYPE);
]])[]dnl
])[]dnl

# ----------------------------------------------------------------------

AC_DEFUN([ERT],
[
AC_ARG_WITH([ert],
            [AS_HELP_STRING([--with-ert=<root>], [Use ERT libraries])],
            [], [with_ert=no])

use_ert=no

AS_IF([test x"${with_ert}" != x"no"],
[
  _ert_LDFLAGS_SAVE="${LDFLAGS}"
  _ert_LIBS_SAVE="${LIBS}"
  _ert_CPPFLAGS_SAVE="${CPPFLAGS}"
  _ert_CFLAGS_SAVE="${CFLAGS}"

  ERT_CPPFLAGS=
  ERT_LDFLAGS=
  ERT_LIBS="-lecl -lgeometry -lert_util -lpthread -lz -lgomp"
  AS_IF([test x"${with_ert}" != x"yes"],
    [ERT_LDFLAGS="-L${with_ert}/lib"
     ERT_CPPFLAGS="-I${with_ert}/include"], [:])[]dnl

  CFLAGS="-std=gnu99"
  CPPFLAGS="${ERT_CPPFLAGS} ${CPPFLAGS}"
  LDFLAGS="${ERT_LDFLAGS} ${LDFLAGS}"
  LIBS="${ERT_LIBS} ${LIBS}"

  AC_LINK_IFELSE([_ERT_SOURCE_TEXT],
    [use_ert=yes], [use_ert=no])

  LIBS="${_ert_LIBS_SAVE}"
  CPPFLAGS="${_ert_CPPFLAGS_SAVE}"
  LDFLAGS="${_ert_LDFLAGS_SAVE}"
  CFLAGS="${_ert_CFLAGS_SAVE}"

  AS_IF([test x"${use_ert}" = x"yes"],
   [AC_SUBST([ERT_CPPFLAGS])
    AC_SUBST([ERT_LDFLAGS])
    AC_SUBST([ERT_LIBS])
    AC_DEFINE([HAVE_ERT], [1],
              [Are the ERT libraries available for reading and writing ECLIPSE files.])],dnl
   [:])[]dnl
], [:])[]dnl

AM_CONDITIONAL([HAVE_ERT], [test x"${use_ert}" != x"no"])

# AC_MSG_ERROR(
# [**** ERT_CPPFLAGS = ${ERT_CPPFLAGS} ****
# **** ERT_LDFLAGS = ${ERT_LDFLAGS} ****
# **** ERT_LIBS = ${ERT_LIBS} ****
# ])
])
