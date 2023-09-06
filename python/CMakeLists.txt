find_package(pybind11 2.11 CONFIG)

if(NOT pybind11_FOUND)
  add_subdirectory( pybind11 )
  # SYSTEM is supported for embedded pybind source (>= 2.3.0)
  # pybind11_add_module SYSTEM argument was deprecated in pybind11 version 2.6
  # See: https://pybind11.readthedocs.io/en/stable/upgrade.html
  #  if(pybind11_VERSION VERSION_LESS "2.6.0")
  #      set(PYBIND11_SYSTEM "SYSTEM")
elseif(pybind11_VERSION VERSION_GREATER_EQUAL "2.3.0")
    if (pybind11_VERSION VERSION_LESS "2.6.0")
        if(CMAKE_VERSION VERSION_LESS 3.12)
            # SYSTEM is only support for pybind >= 2.3.0 and pybind < 2.6.0 and cmake <= 3.12.0 (Unfortunately querying the
            # version does not work. AT least we know the embedded one.
            set(PYBIND11_SYSTEM "SYSTEM")
        endif()
    endif()
endif()
add_subdirectory( simulators )
