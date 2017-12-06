# Build system for opm-simulators tutorials
# Note: This is independent of the main opm-simulators build system.
# It assumes an opm-simulators library is installed at the system level.

project(opm-simulators-tutorials)

cmake_minimum_required(VERSION 2.8)

find_package(opm-core REQUIRED)
find_package(opm-simulators REQUIRED)

include_directories(${opm-core_INCLUDE_DIRS} ${opm-simulators_INCLUDE_DIRS})
add_definitions(${opm-core_DEFINITIONS} ${opm-simulators_DEFINITIONS})

foreach(tut 1 2 3 4)
  add_executable(tutorial${tut} tutorial${tut}.cpp)
  target_link_libraries(tutorial${tut} ${opm-core_LIBRARIES}
                                       ${opm-simulators_LIBRARIES})
endforeach()
