pybind11_add_module(simulators simulators.cpp)

set_target_properties( simulators PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/python/opm2 )

target_sources(simulators
  PRIVATE
    ../../flow/flow_ebos_blackoil.cpp)

target_link_libraries( simulators PRIVATE opmsimulators )

install(TARGETS simulators DESTINATION ${PYTHON_INSTALL_PREFIX}/simulators)
