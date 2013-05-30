# - Find Notay's Algebraic Multigrid Solver
#
# Set the path to the source directory of AGMG in the cache variable
# AGMG_ROOT.
#
# Note the difference between AGMG_DIR and AGMG_ROOT. The former will
# cause find_package to switch to config mode and search for a file
# named agmg-config.cmake, thereby bypassing this module altogether,
# whereas the latter communicates the location of the library to this
# module.
#
# When found, add the contents of AGMG_SOURCES to your own list of
# sources to compile and link for the target.
#
# Use define_fc_func from UseFortranWrappers to write FC_FUNC to your
# own config.h, and declare the function dagmg using this macro.

find_file (AGMG_SOURCES
  dagmg.f90
  PATHS ${AGMG_ROOT}
  PATH_SUFFIXES SRC
  DOC "Yvan Notay's Algebraic Multigrid Solver, Double Precision version"
  NO_DEFAULT_PATH
  )

# make sure that we can compile Fortran code
if (AGMG_SOURCES)
  enable_language (Fortran)
endif (AGMG_SOURCES)

# set value for config.h
if (AGMG_SOURCES)
  set (HAVE_AGMG 1)
else (AGMG_SOURCES)
  set (HAVE_AGMG 0)
endif (AGMG_SOURCES)

# handle REQUIRED and QUIET standard options
include (FindPackageHandleStandardArgs)
find_package_handle_standard_args (AGMG
  DEFAULT_MSG
  AGMG_SOURCES
  CMAKE_Fortran_COMPILER_SUPPORTS_F90
  )
