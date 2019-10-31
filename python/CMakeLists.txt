find_package(PythonInterp REQUIRED)

set( CMAKE_CXX_FLAGS "-std=c++11" )

set( PYTHON_DIR ${CMAKE_SOURCE_DIR}/python )

list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/python/pycmake)

option( USE_RPATH "Embed RPATH in libraries and binaries" OFF)
if (USE_RPATH)
    SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}")
    set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
endif ()

add_subdirectory( pybind11 )

