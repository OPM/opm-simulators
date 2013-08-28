# - Build satellites that are dependent of main library
#
# Enumerate all source code in a "satellite" directory such as tests/,
# compile each of them and optionally set them as a test for CTest to
# run. They will be linked to the main library created by the project.
#
# The following suffices must be defined for the opm prefix passed as
# parameter:
#
#	_LINKER_FLAGS   Necessary flags to link with this library
#	_TARGET         CMake target which creates the library
#	_LIBRARIES      Other dependencies that must also be linked

# Synopsis:
#	opm_compile_satellites (opm satellite excl_all test_regexp)
#
# Parameters:
#	opm             Prefix of the variable which contain information
#	                about the library these satellites depends on, e.g.
#	                pass "opm-core" if opm-core_TARGET is the name of
#	                the target the builds this library. Variables with
#	                suffixes _TARGET and _LIBRARIES must exist.
#
#	satellite       Prefix of variable which contain the names of the
#	                files, e.g. pass "tests" if the files are in the
#	                variable tests_SOURCES. Variables with suffixes
#	                _DATAFILES, _SOURCES and _DIR should exist. This
#	                name is also used as name of the target that builds
#	                all these files.
#
#	excl_all        EXCLUDE_FROM_ALL if these targets should not be built by
#	                default, otherwise empty string.
#
#	test_regexp     Regular expression which picks the name of a test
#	                out of the filename, or blank if no test should be
#	                setup.
#
# Example:
#	opm_compile_satellites (opm-core test "" "^test_([^/]*)$")
#
macro (opm_compile_satellites opm satellite excl_all test_regexp)
  # if we are going to build the tests always, then make sure that
  # the datafiles are present too
  if (NOT (${excl_all} MATCHES "EXCLUDE_FROM_ALL"))
	set (_incl_all "ALL")
  else (NOT (${excl_all} MATCHES "EXCLUDE_FROM_ALL"))
	set (_incl_all "")
  endif (NOT (${excl_all} MATCHES "EXCLUDE_FROM_ALL"))

  # if a set of datafiles has been setup, pull those in
  add_custom_target (${satellite} ${_incl_all})
  if (${satellite}_DATAFILES)
	add_dependencies (${satellite} ${${satellite}_DATAFILES})
  endif (${satellite}_DATAFILES)

  # compile each of these separately
  foreach (_sat_FILE IN LISTS ${satellite}_SOURCES)
	get_filename_component (_sat_NAME "${_sat_FILE}" NAME_WE)
	add_executable (${_sat_NAME} ${excl_all} ${_sat_FILE})
	add_dependencies (${satellite} ${_sat_NAME})
	set_target_properties (${_sat_NAME} PROPERTIES
	  LINK_FLAGS "${${opm}_LINKER_FLAGS_STR}"
	  )
	# are we building a test? luckily, the testing framework doesn't
	# require anything else, so we don't have to figure out where it
	# should go in the library list
	if (NOT "${test_regexp}" STREQUAL "")
	  set (_test_lib "${Boost_UNIT_TEST_FRAMEWORK_LIBRARY}")
	else (NOT "${test_regexp}" STREQUAL "")
	  set (_test_lib "")
	endif (NOT "${test_regexp}" STREQUAL "")
	target_link_libraries (${_sat_NAME} ${${opm}_TARGET} ${${opm}_LIBRARIES} ${_test_lib})
	strip_debug_symbols (${_sat_NAME} _sat_DEBUG)
	list (APPEND ${satellite}_DEBUG ${_sat_DEBUG})

	# variable with regular expression doubles as a flag for
	# whether tests should be setup or not
	if (NOT "${test_regexp}" STREQUAL "")
	  foreach (_regexp IN ITEMS ${test_regexp})
		if ("${_sat_NAME}" MATCHES "${_regexp}")
		  string (REGEX REPLACE "${_regexp}" "\\1" _sat_FANCY "${_sat_NAME}")
		endif ("${_sat_NAME}" MATCHES "${_regexp}")
	  endforeach (_regexp)
	  get_target_property (_sat_LOC ${_sat_NAME} LOCATION)
	  if (CMAKE_VERSION VERSION_LESS "2.8.4")
		add_test (
		  NAME ${_sat_FANCY}
		  COMMAND ${CMAKE_COMMAND} -E chdir "${PROJECT_BINARY_DIR}/${${satellite}_DIR}" ${_sat_LOC}
		  )
	  else (CMAKE_VERSION VERSION_LESS "2.8.4")
		add_test (${_sat_FANCY} ${_sat_LOC})
		# run the test in the directory where the data files are
		set_tests_properties (${_sat_FANCY} PROPERTIES
		  WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/${${satellite}_DIR}
		  )
	  endif (CMAKE_VERSION VERSION_LESS "2.8.4")
	endif(NOT "${test_regexp}" STREQUAL "")

	# if this program on the list of files that should be distributed?
	# we check by the name of the source file
	list (FIND ${satellite}_SOURCES_DIST "${_sat_FILE}" _is_util)
	if (NOT (_is_util EQUAL -1))
	  install (TARGETS ${_sat_NAME} RUNTIME
		DESTINATION bin${${opm}_VER_DIR}/
		)
	endif (NOT (_is_util EQUAL -1))
  endforeach (_sat_FILE)
endmacro (opm_compile_satellites opm prefix)

# Synopsis:
#	opm_data (satellite target dirname files)
#
# provides these output variables:
#
#	${satellite}_INPUT_FILES   List of all files that are copied
#	${satellite}_DATAFILES     Name of target which copies these files
#
# Example:
#
#	opm_data (tests datafiles "tests/")
#
macro (opm_data satellite target dirname)
  # even if there are no datafiles, create the directory so the
  # satellite programs have a homedir to run in
  execute_process (
	COMMAND ${CMAKE_COMMAND} -E make_directory ${dirname}
	)

  # if ever huge test datafiles are necessary, then change this
  # into "create_symlink" (on UNIX only, apparently)
  set (make_avail "copy")

  # provide datafiles as inputs for the tests, by copying them
  # to a tests/ directory in the output tree (if different)
  set (${satellite}_INPUT_FILES)
  if (NOT PROJECT_SOURCE_DIR STREQUAL PROJECT_BINARY_DIR)
	foreach (input_datafile IN LISTS ${satellite}_DATA)
	  file (RELATIVE_PATH rel_datafile "${PROJECT_SOURCE_DIR}" ${input_datafile})
	  set (output_datafile "${PROJECT_BINARY_DIR}/${rel_datafile}")
	  add_custom_command (
		OUTPUT ${output_datafile}
		COMMAND ${CMAKE_COMMAND}
		ARGS -E ${make_avail} ${input_datafile} ${output_datafile}
		DEPENDS ${input_datafile}
		VERBATIM
		)
	  list (APPEND ${satellite}_INPUT_FILES "${output_datafile}")
	endforeach (input_datafile)
  endif(NOT PROJECT_SOURCE_DIR STREQUAL PROJECT_BINARY_DIR)

  # setup a target which does all the copying
  set (${satellite}_DATAFILES "${target}")
  add_custom_target (${${satellite}_DATAFILES}
	DEPENDS ${${satellite}_INPUT_FILES}
	COMMENT "Making \"${satellite}\" data available in output tree"
	)
endmacro (opm_data satellite target dirname files)
