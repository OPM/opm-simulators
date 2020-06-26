pybind11_add_module(simulators simulators.cpp SYSTEM)

set_target_properties( simulators PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/python/opm2 )

target_sources(simulators
  PRIVATE
    ../../opm/simulators/utils/moduleVersion.cpp
    ../../flow/flow_ebos_blackoil.cpp)

target_link_libraries( simulators PRIVATE opmsimulators )

set(PYTHON_PACKAGE_PATH "site-packages")
set(PYTHON_INSTALL_PREFIX "lib/python${PYTHON_VERSION_MAJOR}.${PYTHON_VERSION_MINOR}/${PYTHON_PACKAGE_PATH}" CACHE STRING "Subdirectory to install Python modules in")

install(TARGETS simulators DESTINATION ${DEST_PREFIX}${CMAKE_INSTALL_PREFIX}/${PYTHON_INSTALL_PREFIX}/opm)
