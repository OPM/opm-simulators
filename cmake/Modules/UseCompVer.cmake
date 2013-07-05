# - Get compiler version

# probe the GCC version, returns empty string if GCC is not compiler
function (get_gcc_version language ver_name)
  if(CMAKE_${language}_COMPILER_ID STREQUAL GNU)  
	# exec_program is deprecated, but execute_process does't work :-(
	exec_program (${CMAKE_${language}_COMPILER}
	  ARGS ${CMAKE_${language}_COMPILER_ARG1} -dumpversion
	  OUTPUT_VARIABLE _version
	  )
	set (${ver_name} ${_version} PARENT_SCOPE)
  else (CMAKE_${language}_COMPILER_ID STREQUAL GNU)
	set (${ver_name} "" PARENT_SCOPE)
  endif (CMAKE_${language}_COMPILER_ID STREQUAL GNU)
endfunction (get_gcc_version ver_name)

# less reliable, but includes the patch number
function (get_gcc_patch language ver_name)
  if(CMAKE_${language}_COMPILER_ID STREQUAL GNU)
	# exec_program is deprecated, but execute_process does't work :-(
	exec_program (${CMAKE_${language}_COMPILER}
	  ARGS ${CMAKE_${language}_COMPILER_ARG1} --version
	  OUTPUT_VARIABLE _version
	  )
	# split multi-line string into list
	if (WIN32)
	  string (REPLACE "\r\n" ";" _version "${_version}")
	else (WIN32)
	  string (REPLACE "\n" ";" _version "${_version}")
	endif (WIN32)
	# only keep first line
	list (GET _version 0 _version)
	# extract version number from it (this is the fragile part)
	string (REGEX REPLACE "^[^\\(]+(\\([^\\)]*\\))?[\ \t]*([0-9]+\\.[0-9]+\\.[0-9]+)(.*\\(.*\\))?" "\\2" _version "${_version}")
	# return this to the caller
	set (${ver_name} ${_version} PARENT_SCOPE)
  else (CMAKE_${language}_COMPILER_ID STREQUAL GNU)
	set (${ver_name} "" PARENT_SCOPE)
  endif (CMAKE_${language}_COMPILER_ID STREQUAL GNU)
endfunction (get_gcc_patch language ver_name)

function (compiler_info)
  if (CMAKE_COMPILER_IS_GNUCXX)
	get_gcc_patch (CXX version)
	message (STATUS "GNU C++ compiler version: ${version}")
  endif (CMAKE_COMPILER_IS_GNUCXX)
endfunction (compiler_info)

function (get_ld_version ver_name)
  # run linker to get the version number. interestingly, this option works
  # (for our purposes) on all major platforms (Linux, Mac OS X and Windows);
  # it returns the program version although it may have ended in error
  exec_program (${CMAKE_LINKER}
	ARGS "-v"
	OUTPUT_VARIABLE _version
	)

  # keep only first line, even on Mac OS X there is no line end
  list (GET _version 0 _version)

  # format of the version string is platform-specific
  if (NOT WIN32)
	if (APPLE)
	  string (REGEX REPLACE ".*, from Apple (.*\)" "\\1" _version "${_version}")
	else (APPLE)
	  # assuming some GNU toolchain now
	  string (REGEX REPLACE "GNU ([a-zA-Z0-9_]*) (version|\\(.*\\)) ([^\\ ]*).*" "\\1 \\3" _version "${_version}")
	endif (APPLE)
  endif (NOT WIN32)

  # return the string to the caller
  set (${ver_name} "${_version}" PARENT_SCOPE)
endfunction (get_ld_version ver_name)

function (linker_info)
  get_ld_version (version)
  message (STATUS "Linker: ${version}")
endfunction (linker_info)
