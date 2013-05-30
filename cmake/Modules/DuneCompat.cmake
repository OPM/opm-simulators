# - Dunecontrol compatibility
#
# Enables this build to be a part of a dunecontrol chain. The
# DUNE_CHECK_MODULES macro greps the top-level Makefile for the
# presence of the abs_top_srcdir variable (!) and uses that as
# the include directory of a module. Also, the contents are not
# checked so if the variable is not present, it generates an
# invalid command line (-I without argument) which causes the
# autoconf probe to fail. This module patches our Makefile (!!)
# so the necessary string will be there; in itself this string
# has no use for us, it is solemnly to satisfy the M4 scripts.

if (CMAKE_GENERATOR MATCHES "Unix Makefiles")
  # we need an up-to-date, patched Makefile. this is always checked when
  # the makefile is run, thus the need for a marker file to keep a
  # timestamp to see when it was last patched (by us)
  # amazingly, nothing depends on the generated Makefile, so this can be
  # run whenever in the build without trigging a compile of e.g. config.h
  add_custom_target (dune-compat ALL
	COMMAND ${CMAKE_COMMAND} -DCMAKE_HOME_DIRECTORY=${CMAKE_HOME_DIRECTORY} -P ${PROJECT_SOURCE_DIR}/cmake/Scripts/DuneCompat2.cmake
	COMMENT "Patching Makefile to be DUNE compatible"
	)
endif (CMAKE_GENERATOR MATCHES "Unix Makefiles")

# dunecontrol refuses to use a build tree as module directory unless
# there is a dune.module in it
if (NOT PROJECT_SOURCE_DIR STREQUAL PROJECT_BINARY_DIR)
  execute_process (COMMAND
	${CMAKE_COMMAND} -E copy_if_different ${PROJECT_SOURCE_DIR}/dune.module ${PROJECT_BINARY_DIR}/dune.module
	)
endif (NOT PROJECT_SOURCE_DIR STREQUAL PROJECT_BINARY_DIR)
