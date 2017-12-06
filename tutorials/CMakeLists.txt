# Build system for opm-core tutorials
# Note: This is independent of the main opm-core build system.
# It assumes an opm-core library is installed at the system level.

project(opm-core-tutorials)

cmake_minimum_required(VERSION 2.8)

find_package(opm-core REQUIRED)

include_directories(${opm-core_INCLUDE_DIRS})
add_definitions(${opm-core_DEFINITIONS})

foreach(tut 1 2 3 4)
  add_executable(tutorial${tut} tutorial${tut}.cpp)
  target_link_libraries(tutorial${tut} ${opm-core_LIBRARIES})
endforeach()
