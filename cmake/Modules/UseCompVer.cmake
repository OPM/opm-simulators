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
