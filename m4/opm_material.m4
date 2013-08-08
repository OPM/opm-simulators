dnl -*- autoconf -*-

dnl locate opm-material library itself; this macro is called by every module
dnl that depends on opm-material
AC_DEFUN([OPM_MATERIAL_CHECK_MODULE],
[
 OPM_CHECK_PKG_MODULE([opm-material],[1.0],[Provides Thermodynamic Relations, Capillary Pressure Curves, etc.])
])

dnl find all prerequisites of opm-core; nothing to do here since this
dnl is done by the CMake module and then stored in the -config file.
AC_DEFUN([OPM_MATERIAL_CHECKS],[])
