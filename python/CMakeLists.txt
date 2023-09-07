find_package(pybind11 2.2 CONFIG)

if(NOT pybind11_FOUND)
  include(DownloadPyBind11)
elseif(pybind11_VERSION VERSION_GREATER_EQUAL "2.3.0")
  if(CMAKE_VERSION VERSION_LESS 3.12)
    # SYSTEM is only support for pybind >= 2.3.0 and cmake <= 3.12.0 (Unfortunately querying the
    # version does not work. AT least we know the embedded one.
    set(PYBIND11_SYSTEM "SYSTEM")
  endif()
endif()
add_subdirectory( simulators )
