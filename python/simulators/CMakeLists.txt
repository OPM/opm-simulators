#  Note: If the new find_package(Python3) is used in the top level CMakeLists.txt, the variable
#    ${PYTHON_EXECUTABLE} is set there to ${Python3_EXECUTABLE}
#
# NOTE: The variable ${PYBIND11_SYSTEM} is set in python/CMakeLists.txt
#   to the value "SYSTEM" or unset, depending on the current version of Pybind11.
#   The value is then forwarded to target_include_directories(), see
#
#  https://cmake.org/cmake/help/latest/command/target_include_directories.html
#  https://pybind11.readthedocs.io/en/stable/compiling.html
#
pybind11_add_module(simulators ${PYBIND11_SYSTEM}
  $<TARGET_OBJECTS:moduleVersion>
  Pybind11Exporter.cpp
  PyBlackOilSimulator.cpp
  )

set(PYTHON_OPM_SIMULATORS_PACKAGE_PATH ${PROJECT_BINARY_DIR}/python/opm/simulators)
set_target_properties( simulators PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${PYTHON_OPM_SIMULATORS_PACKAGE_PATH} )

target_link_libraries( simulators PRIVATE opmsimulators )

execute_process(COMMAND ${PYTHON_EXECUTABLE} -c "
import site, sys
try:
    sys.stdout.write(site.getsitepackages()[-1])
except e:
    sys.stdout.write('')" OUTPUT_VARIABLE PYTHON_SITE_PACKAGES_PATH)

if (PYTHON_SITE_PACKAGES_PATH MATCHES ".*/dist-packages/?" AND
      CMAKE_INSTALL_PREFIX MATCHES "^/usr.*")
  # dist-packages is only used if we install below /usr and python's site packages
  # path matches dist-packages
  set(PYTHON_PACKAGE_PATH "dist-packages")
else()
  set(PYTHON_PACKAGE_PATH "site-packages")
endif()

if(OPM_ENABLE_PYTHON_TESTS)
  set(PYTHON_PATH ${PROJECT_BINARY_DIR}/python:${opm-common_DIR}/python:$ENV{PYTHONPATH})
  # NOTE: See comment in test_basic.py for the reason why we are
  #   splitting the python tests into multiple add_test() tests instead
  #   of having a single "python -m unittest" test call that will run all
  #   the tests in the "test" sub directory.
  add_test(NAME python_basic
      WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/python
      COMMAND ${CMAKE_COMMAND}
      -E env PYTHONPATH=${PYTHON_PATH} ${PYTHON_EXECUTABLE}
      -m unittest test/test_basic.py)
  add_test(NAME python_schedule
      WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/python
      COMMAND ${CMAKE_COMMAND}
      -E env PYTHONPATH=${PYTHON_PATH} ${PYTHON_EXECUTABLE}
      -m unittest test/test_schedule.py)
  add_test(NAME python_throw
      WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/python
      COMMAND ${CMAKE_COMMAND}
      -E env PYTHONPATH=${PYTHON_PATH} ${PYTHON_EXECUTABLE}
      -m unittest test/test_throw.py)
endif()

find_file(PYTHON_INSTALL_PY install.py
  PATHS ${opm-common_DIR} ${opm-common_PYTHON_COMMON_DIR}
  PATH_SUFFIXES python NO_DEFAULT_PATH REQUIRED)

# NOTE: instead of using file( COPY ...) which copies the files at configure time (not at build time)
#  we should add copying of the files at build time such that running "make" or "make all" will
#  update the files if they have been modified. We use the install.py script in opm-common, see also
#  CMakeLists.txt in opm-common
add_custom_target(copy_python ALL
    COMMAND ${PYTHON_EXECUTABLE} "${PYTHON_INSTALL_PY}"
        ${PROJECT_SOURCE_DIR}/python/opm ${PROJECT_BINARY_DIR}/python 0
    COMMAND ${PYTHON_EXECUTABLE} "${PYTHON_INSTALL_PY}"
        ${PROJECT_SOURCE_DIR}/python/test ${PROJECT_BINARY_DIR}/python 0
    COMMAND ${PYTHON_EXECUTABLE} "${PYTHON_INSTALL_PY}"
        ${PROJECT_SOURCE_DIR}/python/test_data ${PROJECT_BINARY_DIR}/python 0
)

# Since the installation of Python code is nonstandard it is protected by an
# extra cmake switch, OPM_INSTALL_PYTHON. If you prefer you can still invoke
# setup.py install manually - optionally with the generated script
# setup-install.sh - and completely bypass cmake in the installation phase.
if (OPM_INSTALL_PYTHON)
  include(PyInstallPrefix)  # from opm-common
  install(TARGETS simulators DESTINATION ${DEST_PREFIX}${CMAKE_INSTALL_PREFIX}/${PYTHON_INSTALL_PREFIX}/opm)
  install(
    CODE "execute_process(
      COMMAND ${PYTHON_EXECUTABLE}
              ${PYTHON_INSTALL_PY}
              ${PROJECT_BINARY_DIR}/python/opm
              ${DEST_PREFIX}${CMAKE_INSTALL_PREFIX}/${PYTHON_INSTALL_PREFIX} 1)")
endif()
