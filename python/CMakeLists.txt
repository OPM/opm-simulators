find_package(pybind11 2.2 CONFIG)

if(NOT pybind11_FOUND)
  add_subdirectory( pybind11 )
  # SYSTEM is supported for embedded pybind source (>= 2.3.0)
  set(PYBIND11_SYSTEM "SYSTEM")
elseif(pybind11_VERSION VERSION_GREATER_EQUAL "2.3.0")
  # SYSTEM is only support for pybind >= 2.3.0 (Unfortunately querying the
  # version does not work. AT least we know the embedded one.
  set(PYBIND11_SYSTEM "SYSTEM")
endif()
add_subdirectory( simulators )
