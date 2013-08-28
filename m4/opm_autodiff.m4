dnl -*- autoconf -*-

dnl locate opm-autodiff library itself; this macro is called by every module
dnl that depends on opm-autodiff.
AC_DEFUN([OPM_AUTODIFF_CHECK_MODULE],
[
 OPM_CHECK_PKG_MODULE([opm-autodiff],[1.0],[Utilities for automatic differentiation and simulators based on AD])
])

dnl find all prerequisites of opm-autodiff; nothing to do here since this
dnl is done by the CMake module and then stored in the -config file.
AC_DEFUN([OPM_AUTODIFF_CHECKS],[])
