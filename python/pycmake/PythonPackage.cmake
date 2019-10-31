#  Copyright (C)  2016 Statoil ASA, Norway.
#
#  pymake is free software: you can redistribute it and/or modify it under the
#  terms of the GNU General Public License as published by the Free Software
#  Foundation, either version 3 of the License, or (at your option) any later
#  version.
#
#  pymake is distributed in the hope that it will be useful, but WITHOUT ANY
#  WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
#  A PARTICULAR PURPOSE.
#
#  See the GNU General Public License at <http://www.gnu.org/licenses/gpl.html>
#  for more details

# This module exports the following functions:
# * to_path_list(var path)
#       which takes a list of paths and constructs a valid PYTHONPATH string,
#       regardless  of platform (path1:path2 for unix, path1;path2 for windows)
#
# * add_python_package(package_name package_path python_files) where 
#       which takes a package name, a path and the series of python files that
#       makes that package. It exports the cmake target package_${package_name}
#       in copies all ${python_files} sources to python/${package_path}, and
#       sets up so you can install with `make install`.
#
# * add_python_test(testname python_test_file)
#       which sets up a test target (using pycmake_test_runner.py, distributed
#       with this module) and registeres it with ctest.
#
# * add_python_example(example testname test_file [args...])
#       which sets up an example program which will be run with the arguments
#       [args...] (can be empty) Useful to make sure some program runs
#       correctly with the given arguments, and which will report as a unit
#       test failure.

configure_file(${CMAKE_CURRENT_LIST_DIR}/pycmake_test_runner.py ${CMAKE_BINARY_DIR}/python/tests/pycmake_test_runner.py COPYONLY)

function(to_path_list var path1)
    if("${CMAKE_HOST_SYSTEM}" MATCHES ".*Windows.*")
        set(sep "\\;")
    else()
        set(sep ":")
    endif()
    set(result "${path1}") # First element doesn't require separator at all...
    foreach(path ${ARGN})
        set(result "${result}${sep}${path}") # .. but other elements do.
    endforeach()
    set(${var} "${result}" PARENT_SCOPE)
endfunction()

if (EXISTS "/etc/debian_version")
    set( PYTHON_PACKAGE_PATH "dist-packages")
else()
    set( PYTHON_PACKAGE_PATH "site-packages")
endif()
set(PYTHON_INSTALL_PREFIX "lib/python${PYTHON_VERSION_MAJOR}.${PYTHON_VERSION_MINOR}/${PYTHON_PACKAGE_PATH}" CACHE STRING "Subdirectory to install Python modules in")

function(add_python_package PACKAGE_NAME PACKAGE_PATH PYTHON_FILES)
    add_custom_target(package_${PACKAGE_NAME} ALL)

    foreach (file ${PYTHON_FILES})
        add_custom_command(TARGET package_${PACKAGE_NAME}
                COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_BINARY_DIR}/python/${PACKAGE_PATH}
                COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_SOURCE_DIR}/${file} ${CMAKE_BINARY_DIR}/python/${PACKAGE_PATH}
                )
    endforeach ()
    set_target_properties(package_${PACKAGE_NAME} PROPERTIES PACKAGE_INSTALL_PATH ${CMAKE_INSTALL_PREFIX}/${PYTHON_INSTALL_PREFIX}/${PACKAGE_PATH})
    install(FILES ${PYTHON_FILES} DESTINATION ${CMAKE_INSTALL_PREFIX}/${PYTHON_INSTALL_PREFIX}/${PACKAGE_PATH})
endfunction()

function(add_python_test TESTNAME PYTHON_TEST_FILE arg)
    configure_file(${PYTHON_TEST_FILE} ${PYTHON_TEST_FILE} COPYONLY)

    add_test(NAME ${TESTNAME}
            WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/python/tests
            COMMAND python pycmake_test_runner.py ${PYTHON_TEST_FILE} ${arg}
            )

    to_path_list(pythonpath "${CMAKE_BINARY_DIR}/python" "$ENV{PYTHONPATH}")
    set_tests_properties(${TESTNAME} PROPERTIES ENVIRONMENT "PYTHONPATH=${pythonpath}")
endfunction()

function(add_python_example TESTNAME PYTHON_TEST_FILE)
    configure_file(${PYTHON_TEST_FILE} ${PYTHON_TEST_FILE} COPYONLY)

    add_test(NAME ${TESTNAME}
            WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/examples
            COMMAND python ${PYTHON_TEST_FILE} ${ARGN}
            )
    to_path_list(pythonpath "${CMAKE_BINARY_DIR}/python" "$ENV{PYTHONPATH}")
    set_tests_properties(${TESTNAME} PROPERTIES ENVIRONMENT "PYTHONPATH=${pythonpath}")
endfunction()

