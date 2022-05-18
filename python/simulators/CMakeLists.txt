# NOTE: we assume that add_subdirectory( pybind11 ) is called from the
#    parent folder's CMakeLists.txt before this CMakeLists.txt is loaded.
#    Therefore, pybind11's CMakeLists.txt has already run
#    find_package(PYTHON) to define variables like ${PYTHON_EXECUTABLE}
#
# NOTE: The variable ${PYBIND11_SYSTEM} is set in python/CMakeLists.txt
#   to the value "SYSTEM" or unset, depending on the current version of Pybind11.
#   The value is then forwarded to target_include_directories(), see
#
#  https://cmake.org/cmake/help/latest/command/target_include_directories.html
#  https://pybind11.readthedocs.io/en/stable/compiling.html
#
pybind11_add_module(simulators ${PYBIND11_SYSTEM}
  PyBlackOilSimulator.cpp
  Pybind11Exporter.cpp)

set(PYTHON_OPM_SIMULATORS_PACKAGE_PATH ${PROJECT_BINARY_DIR}/python/opm2/simulators)
set_target_properties( simulators PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${PYTHON_OPM_SIMULATORS_PACKAGE_PATH} )

target_sources(simulators
  PRIVATE
  $<TARGET_OBJECTS:moduleVersion>
  $<TARGET_OBJECTS:flow_libblackoil>)

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
  if(Python3_EXECUTABLE AND NOT PYTHON_EXECUTABLE)
    set(PYTHON_EXECUTABLE ${Python3_EXECUTABLE})
  endif()
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
endif()

# NOTE: instead of using file( COPY ...) which copies the files at configure time (not at build time)
#  we should add copying of the files at build time such that running "make" or "make all" will
#  update the files if they have been modified. We use the install.py script in opm-common, see also
#  CMakeLists.txt in opm-common
add_custom_target(copy_python ALL
    COMMAND ${PYTHON_EXECUTABLE} ${opm-common_DIR}/python/install.py
      ${PROJECT_SOURCE_DIR}/python ${PROJECT_BINARY_DIR} 0)

# Since the installation of Python code is nonstandard it is protected by an
# extra cmake switch, OPM_INSTALL_PYTHON. If you prefer you can still invoke
# setup.py install manually - optionally with the generated script
# setup-install.sh - and completely bypass cmake in the installation phase.
if (OPM_INSTALL_PYTHON)
  set(PYTHON_INSTALL_PREFIX "lib/python${PYTHON_VERSION_MAJOR}.${PYTHON_VERSION_MINOR}/${PYTHON_PACKAGE_PATH}" CACHE STRING "Subdirectory to install Python modules in")
  install(TARGETS simulators DESTINATION ${DEST_PREFIX}${CMAKE_INSTALL_PREFIX}/${PYTHON_INSTALL_PREFIX}/opm)
  include(PyInstallPrefix)  # from opm-common
  install( CODE "execute_process(COMMAND ${PYTHON_EXECUTABLE} ${opm-common_DIR}/python/install.py ${PROJECT_BINARY_DIR}/python/opm2 ${DEST_PREFIX}${CMAKE_INSTALL_PREFIX}/${PYTHON_INSTALL_PREFIX} 1)")
endif()
