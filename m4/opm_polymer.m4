dnl -*- autoconf -*-

dnl locate opm-polymer library itself; this macro is called by every module
dnl that depends on opm-polymer.
AC_DEFUN([OPM_POLYMER_CHECK_MODULE],
[
 OPM_CHECK_PKG_MODULE([opm-polymer],[1.0],[OPM module for polymer simulations])
])

dnl find all prerequisites of opm-polymer; nothing to do here since this
dnl is done by the CMake module and then stored in the -config file.
AC_DEFUN([OPM_POLYMER_CHECKS],[])
