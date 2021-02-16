# NOTE: we assume that add_subdirectory( pybind11 ) is called from the
#       parent folder's CMakeLists.txt before this CMakeLists.txt is loaded.
#       Therefore, pybind11's CMakeLists.txt has already run
#       find_package(PYTHON) to define variables like
#       ${PYTHON_EXECUTABLE}
#
pybind11_add_module(simulators simulators.cpp SYSTEM)

set_target_properties( simulators PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/python/opm2 )

target_sources(simulators
  PRIVATE
  $<TARGET_OBJECTS:moduleVersion>
  $<TARGET_OBJECTS:flow_libblackoil>)

target_link_libraries( simulators PRIVATE opmsimulators ${opm-simulators_PRIVATE_LIBRARIES})

set(PYTHON_PACKAGE_PATH "site-packages")
set(PYTHON_INSTALL_PREFIX "lib/python${PYTHON_VERSION_MAJOR}.${PYTHON_VERSION_MINOR}/${PYTHON_PACKAGE_PATH}" CACHE STRING "Subdirectory to install Python modules in")

install(TARGETS simulators DESTINATION ${DEST_PREFIX}${CMAKE_INSTALL_PREFIX}/${PYTHON_INSTALL_PREFIX}/opm)

file( COPY ${PROJECT_SOURCE_DIR}/python/test
      DESTINATION ${PROJECT_BINARY_DIR}/python)
file( COPY ${PROJECT_SOURCE_DIR}/python/test_data
      DESTINATION ${PROJECT_BINARY_DIR}/python)

if(OPM_ENABLE_PYTHON_TESTS)
  add_test(NAME python_tests
      WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/python
      COMMAND ${CMAKE_COMMAND}
            -E env PYTHONPATH=${PROJECT_BINARY_DIR}/python:$ENV{PYTHONPATH}
            ${PYTHON_EXECUTABLE} -m unittest )
endif()
