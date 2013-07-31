# - Write version information into the source code
#
# Add an unconditional target to the Makefile which checks the current
# SHA of the source directory and write to a header file if and *only*
# if this has changed (thus we avoid unnecessary rebuilds). By having
# this in the Makefile, we get updated version information even though
# we haven't done any reconfiguring.
#
# The time it takes to probe the VCS for this information and write it
# to the miniature file in negligable.
#
# If the build type is Debug, then we only write a static version
# information as it gets tiresome to rebuild the project everytime one
# makes changes to any of the unit tests.

string (TOUPPER "${CMAKE_BUILD_TYPE}" cmake_build_type_upper_)
if (cmake_build_type_upper_ MATCHES DEBUG)
  file (WRITE "${PROJECT_BINARY_DIR}/project-version.h"
	"#define PROJECT_VERSION \"${${project}_LABEL} (debug)\"\n"
	)
else ()
  if (NOT GIT_FOUND)
	find_package (Git)
  endif ()

  add_custom_target (update-version ALL
	COMMAND ${CMAKE_COMMAND}
	-DCMAKE_HOME_DIRECTORY=${CMAKE_HOME_DIRECTORY}
	-DGIT_EXECUTABLE=${GIT_EXECUTABLE}
	-DPROJECT_SOURCE_DIR=${PROJECT_SOURCE_DIR}
	-DPROJECT_BINARY_DIR=${PROJECT_BINARY_DIR}
	-DPROJECT_LABEL=${${project}_LABEL}
	-P ${PROJECT_SOURCE_DIR}/cmake/Scripts/WriteVerSHA.cmake
	COMMENT "Updating version information"
	)

  # the target above gets built every time thanks to the "ALL" modifier,
  # but it must also be done before the main library so it can pick up
  # any changes it does.
  add_dependencies (${${project}_TARGET} update-version)
endif ()
