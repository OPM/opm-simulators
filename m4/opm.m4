dnl -*- autoconf -*-

dnl OPM_PKG_CONFIG_MODULE (name, version, description)
dnl
dnl Common routine to include configuration module for an OPM project
AC_DEFUN([OPM_CHECK_PKG_MODULE],[
 dnl local variables representing parameters
 m4_pushdef([_opm_name], [$1])
 m4_pushdef([_opm_version], [$2])
 m4_pushdef([_opm_description], [$3])
 
 dnl macro-friendly version of the name; uppercase and with dashes
 dnl replaced with underscores
 m4_pushdef([_opm_module], [m4_translit(_opm_name,[-],[_])])
 m4_pushdef([_OPM_MODULE], [m4_toupper(_opm_module)])

 dnl if we are given the location as a parameter, look there first
 AC_ARG_WITH(_opm_name,
  AS_HELP_STRING([--with-_opm_name=PATH],[_opm_description directory]))

 AS_IF([test -n "$with_[]_opm_module"],[
  export PKG_CONFIG_PATH=$with_[]_opm_module:$PKG_CONFIG_PATH
 ])

 dnl let pkg-config do the heavy lifting of finding the .pc file
 PKG_CHECK_MODULES(_OPM_MODULE,[_opm_name = _opm_version],[
  AC_DEFINE(HAVE_[]_OPM_MODULE,[1],[_opm_description available])
 ])

 dnl TODO: here we could call PKG_CONFIG --variable if we need more

 dnl make flag available for Makefiles too
 AM_CONDITIONAL(HAVE_[]_OPM_MODULE, test x$HAVE_[]_OPM_MODULE = x1)

 dnl add our libraries to the global list of compiler and linker options
 DUNE_CPPFLAGS="$DUNE_CPPFLAGS $_OPM_MODULE[]_CFLAGS"
 DUNE_LIBS="$DUNE_LIBS $_OPM_MODULE[]_LIBS"

 # add this module to summary (if we are used with dunecontrol)
 ifdef([DUNE_MODULE_ADD_SUMMARY_ENTRY],[
  DUNE_MODULE_ADD_SUMMARY_ENTRY(_opm_name)
 ])

 dnl cleanup
 m4_popdef([_OPM_MODULE])
 m4_popdef([_opm_module])
 m4_popdef([_opm_description])
 m4_popdef([_opm_version])
 m4_popdef([_opm_name])
])
