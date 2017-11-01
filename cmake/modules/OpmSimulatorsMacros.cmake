# .. cmake_module::
#
# This module's content is executed whenever a Dune module requires or
# suggests opm-simulators!
#

find_package(Eigen3)

if(EIGEN3_FOUND)
  set(HAVE_EIGEN "${EIGEN3_FOUND}")

  dune_register_package_flags(
    INCLUDE_DIRS "${EIGEN3_INCLUDE_DIRS}")
endif()

find_package(Boost
  COMPONENTS filesystem regex system date_time
  REQUIRED)
