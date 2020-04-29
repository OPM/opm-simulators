###########################################################################
#                                                                         #
# Note: The bulk of the build system is located in the cmake/ directory.  #
#       This file only contains the specializations for this particular   #
#       project. Most likely you are interested in editing one of these   #
#       files instead:                                                    #
#                                                                         #
#       dune.module                              Name and version number  #
#       CMakeLists_files.cmake                   Path of source files     #
#       cmake/Modules/${project}-prereqs.cmake   Dependencies             #
#                                                                         #
###########################################################################

# Mandatory call to project

project(opm-simulators C CXX)

cmake_minimum_required (VERSION 2.8)


option(SIBLING_SEARCH "Search for other modules in sibling directories?" ON)
set( USE_OPENMP_DEFAULT OFF ) # Use of OpenMP is considered experimental
option(BUILD_FLOW "Build the production oriented flow simulator?" ON)
option(BUILD_FLOW_BLACKOIL_ONLY "Build the production oriented flow simulator only supporting the blackoil model?" OFF)
option(BUILD_FLOW_VARIANTS "Build the variants for flow by default?" OFF)
option(BUILD_EBOS "Build the research oriented ebos simulator?" ON)
option(BUILD_EBOS_EXTENSIONS "Build the variants for various extensions of ebos by default?" OFF)
option(BUILD_EBOS_DEBUG_EXTENSIONS "Build the ebos variants which are purely for debugging by default?" OFF)

option(ENABLE_3DPROPS_TESTING "Build and use the new experimental 3D properties" OFF)
if (ENABLE_3DPROPS_TESTING)
  add_definitions(-DENABLE_3DPROPS_TESTING)
endif()

if(SIBLING_SEARCH AND NOT opm-common_DIR)
  # guess the sibling dir
  get_filename_component(_leaf_dir_name ${PROJECT_BINARY_DIR} NAME)
  get_filename_component(_parent_full_dir ${PROJECT_BINARY_DIR} DIRECTORY)
  get_filename_component(_parent_dir_name ${_parent_full_dir} NAME)
  #Try if <module-name>/<build-dir> is used
  get_filename_component(_modules_dir ${_parent_full_dir} DIRECTORY)
  if(IS_DIRECTORY ${_modules_dir}/opm-common/${_leaf_dir_name})
    set(opm-common_DIR ${_modules_dir}/opm-common/${_leaf_dir_name})
  else()
    string(REPLACE ${PROJECT_NAME} opm-common _opm_common_leaf ${_leaf_dir_name})
    if(NOT _leaf_dir_name STREQUAL _opm_common_leaf
        AND IS_DIRECTORY ${_parent_full_dir}/${_opm_common_leaf})
      # We are using build directories named <prefix><module-name><postfix>
      set(opm-common_DIR ${_parent_full_dir}/${_opm_common_leaf})
    elseif(IS_DIRECTORY ${_parent_full_dir}/opm-common)
      # All modules are in a common build dir
      set(opm-common_DIR "${_parent_full_dir}/opm-common")
    endif()
  endif()
endif()
if(opm-common_DIR AND NOT IS_DIRECTORY ${opm-common_DIR})
  message(WARNING "Value ${opm-common_DIR} passed to variable"
    " opm-common_DIR is not a directory")
endif()

find_package(opm-common REQUIRED)

include(OpmInit)
OpmSetPolicies()

# not the same location as most of the other projects? this hook overrides
macro (dir_hook)
endmacro (dir_hook)

# project information is in dune.module. Read this file and set variables.
# we cannot generate dune.module since it is read by dunecontrol before
# the build starts, so it makes sense to keep the data there then.
include (OpmInit)

# Look for the opm-tests repository; if found the variable
# HAVE_OPM_TESTS will be set to true.
include(Findopm-tests)

# list of prerequisites for this particular project; this is in a
# separate file (in cmake/Modules sub-directory) because it is shared
# with the find module
include ("${project}-prereqs")

# Make sure we are using the same compiler underneath
# NVCC as for the rest. In the case that NVCC does not support
# that compiler it will error out. Unfortunately this will only
# work for CMake >= 3.8. We found no way to make FindCUDA.cmake error
# out. It seems to ignore CMAKE_NVCC_FLAGS and CMAKE. Additionally
# our way of specifying cuda source files never worked for CMake
# version < 3.8. Hence we deactivate cuda for these versions.
# We use "CMAKE_VERSION VERSION_GREATER 3.7.9" instead of
# CMAKE_VERSION VERSION_GREATER_EQUAL 3.8, because of backwarsd
# compatibility to cmake 3.6 and lower.
if(NOT CMAKE_DISABLE_FIND_PACKAGE_CUDA AND
    CMAKE_VERSION VERSION_GREATER 3.7.9)
  if(CMAKE_BUILD_TYPE)
    set(_flags_suffix "_${CMAKE_BUILD_TYPE}")
  endif()
  if(NOT DEFINED ENV{CUDAHOSTCXX} AND NOT DEFINED CMAKE_CUDA_HOST_COMPILER AND
      (NOT CMAKE_CUDA_FLAGS${_flags_suffix} OR NOT CMAKE_CUDA_FLAGS${_flags_suffix} MATCHES ".*-ccbin .*"))
    message(STATUS "Setting CUDA host compiler CMAKE_CUDA_HOST_COMPILER to ${CMAKE_CXX_COMPILER} to "
      "prevent incompatibilities. Note that this might report that there "
      "is not CUDA compiler if your system's CUDA compiler does not support "
      "${CMAKE_CXX_COMPILER}.")
    # check_language does not seem to care about ${CMAKE_CUDA_FLAGS} or $(CUDA_NVCC_FLAGS}.
    # Hence we set CMAKE_CUDA_HOST_COMPILER to our C++ compiler.
    # In check_language(CUDA) we will get an error if we in addition put
    # "-ccbin ${CMAKE_CXX_COMPILER}" into CMAKE_CUDA_FLAGS. It results
    # in "${NVCC} -ccbin=${CMAKE_CXX_COMPILER}   -ccbin ${CMAKE_CXX_COMPILER}"
    # which causes nvcc to abort
    set(CMAKE_CUDA_HOST_COMPILER ${CMAKE_CXX_COMPILER})
    set(ENV{CUDAHOSTCXX} ${CMAKE_CUDA_HOST_COMPILER}) # The only thing honored by check_language(CUDA)!
  endif()
  include(CheckLanguage)
  check_language(CUDA)
  if(CMAKE_CUDA_COMPILER)
    # OPTIONAL is ignored. Hence the magic above to check whether enabling CUDA works
    enable_language(CUDA OPTIONAL)
    # While the documentation says that it is deprecated, FindCUDA seems the
    # only easy way to determine the cublas and cusparse libraries.
    # Hence we call it unconditionally
    # The WellContributions kernel uses __shfl_down_sync, which was introduced in CUDA 9.0
    find_package(CUDA)
  endif()
  if(CUDA_FOUND AND CUDA_VERSION VERSION_LESS "9.0")
    set(CUDA_FOUND OFF)
    message(WARNING "Deactivating CUDA as we require version 9.0 or newer."
      " Found only CUDA version ${CUDA_VERSION}.")
  endif()
endif()

if(CUDA_FOUND)
  set(HAVE_CUDA 1)
  include_directories(${CUDA_INCLUDE_DIRS})
endif()


find_package(OpenCL)

if(OpenCL_FOUND)
  # the current OpenCL implementation relies on cl.hpp, not cl2.hpp
  # make sure it is available, otherwise disable OpenCL
  find_file(CL_HPP CL/cl.hpp HINTS ${OpenCL_INCLUDE_DIRS})
  if(CL_HPP)
    set(HAVE_OPENCL 1)
    include_directories(${OpenCL_INCLUDE_DIRS})
  else()
    message(WARNING " OpenCL was found, but this version of opm-simulators relies on CL/cl.hpp, which implements OpenCL 1.0, 1.1 and 1.2.\n Deactivating OpenCL")
    set(OpenCL_FOUND OFF)
    set(OPENCL_FOUND OFF)
  endif()
endif()

# read the list of components from this file (in the project directory);
# it should set various lists with the names of the files to include
include (CMakeLists_files.cmake)

macro (config_hook)
  opm_need_version_of ("dune-common")
  opm_need_version_of ("dune-istl")
  if(dune-fem_FOUND)
      opm_need_version_of ("dune-fem")
  endif()
  opm_need_version_of ("opm-models")
endmacro (config_hook)

macro (prereqs_hook)
endmacro (prereqs_hook)

macro (sources_hook)
endmacro (sources_hook)

macro (fortran_hook)
endmacro (fortran_hook)

macro (files_hook)
endmacro (files_hook)

macro (tests_hook)
endmacro (tests_hook)

# all setup common to the OPM library modules is done here
include (OpmLibMain)

if (HAVE_OPM_TESTS)
    include (${CMAKE_CURRENT_SOURCE_DIR}/compareECLFiles.cmake)
endif()

opm_set_test_driver(${CMAKE_CURRENT_SOURCE_DIR}/tests/run-parallel-unitTest.sh "")

opm_add_test(test_gatherconvergencereport
  DEPENDS "opmsimulators"
  LIBRARIES opmsimulators ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY}
  SOURCES
    tests/test_gatherconvergencereport.cpp
  CONDITION
    MPI_FOUND AND Boost_UNIT_TEST_FRAMEWORK_FOUND
  DRIVER_ARGS
    5 ${PROJECT_BINARY_DIR}
)

opm_add_test(test_gatherdeferredlogger
  DEPENDS "opmsimulators"
  LIBRARIES opmsimulators ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY}
  SOURCES
    tests/test_gatherdeferredlogger.cpp
  CONDITION
    MPI_FOUND AND Boost_UNIT_TEST_FRAMEWORK_FOUND
  DRIVER_ARGS
    5 ${PROJECT_BINARY_DIR}
)

include(OpmBashCompletion)

if (NOT BUILD_FLOW)
  set(FLOW_DEFAULT_ENABLE_IF "FALSE")
else()
  set(FLOW_DEFAULT_ENABLE_IF "TRUE")
endif()

add_library(moduleVersion OBJECT opm/simulators/utils/moduleVersion.cpp)

# Strictly we only depend on the update-version target,
# but this is not exposed in a super-build.
add_dependencies(moduleVersion opmsimulators)

# the production oriented general-purpose ECL simulator
opm_add_test(flow
  ONLY_COMPILE
  ALWAYS_ENABLE
  DEFAULT_ENABLE_IF ${FLOW_DEFAULT_ENABLE_IF}
  DEPENDS opmsimulators
  LIBRARIES opmsimulators
  SOURCES
  flow/flow.cpp
  flow/flow_ebos_blackoil.cpp
  flow/flow_ebos_gasoil.cpp
  flow/flow_ebos_oilwater.cpp
  flow/flow_ebos_polymer.cpp
  flow/flow_ebos_foam.cpp
  flow/flow_ebos_brine.cpp
  flow/flow_ebos_oilwater_brine.cpp
  flow/flow_ebos_solvent.cpp
  flow/flow_ebos_energy.cpp
  flow/flow_ebos_oilwater_polymer.cpp
  flow/flow_ebos_oilwater_polymer_injectivity.cpp
  $<TARGET_OBJECTS:moduleVersion>)


if (NOT BUILD_FLOW_BLACKOIL_ONLY)
  set(FLOW_BLACKOIL_ONLY_DEFAULT_ENABLE_IF "FALSE")
else()
  set(FLOW_BLACKOIL_ONLY_DEFAULT_ENABLE_IF "TRUE")
endif()

# the production oriented general-purpose ECL simulator
opm_add_test(flow_blackoil
  ONLY_COMPILE
  ALWAYS_ENABLE
  DEFAULT_ENABLE_IF ${FLOW_BLACKOIL_ONLY_DEFAULT_ENABLE_IF}
  DEPENDS opmsimulators
  LIBRARIES opmsimulators
  SOURCES
  flow/flow.cpp
  flow/flow_ebos_blackoil.cpp
  $<TARGET_OBJECTS:moduleVersion>)
target_compile_definitions(flow_blackoil PRIVATE "FLOW_BLACKOIL_ONLY")

if (NOT BUILD_FLOW_VARIANTS)
  set(FLOW_VARIANTS_DEFAULT_ENABLE_IF "FALSE")
else()
  set(FLOW_VARIANTS_DEFAULT_ENABLE_IF "TRUE")
endif()

# Variant versions of Flow.
opm_add_test(flow_blackoil_dunecpr
  ONLY_COMPILE
  DEFAULT_ENABLE_IF ${FLOW_VARIANTS_DEFAULT_ENABLE_IF}
  SOURCES
  flow/flow_blackoil_dunecpr.cpp
  $<TARGET_OBJECTS:moduleVersion>
  EXE_NAME flow_blackoil_dunecpr
  DEPENDS opmsimulators
  LIBRARIES opmsimulators)

opm_add_test(flow_blackoil_dunecpr_noflex
  ONLY_COMPILE
  DEFAULT_ENABLE_IF ${FLOW_VARIANTS_DEFAULT_ENABLE_IF}
  SOURCES flow/flow_blackoil_dunecpr_noflex.cpp
  EXE_NAME flow_blackoil_dunecpr_noflex
  DEPENDS "opmsimulators"
  LIBRARIES "opmsimulators")

opm_add_test(flow_blackoil_dunecpr_seq
  ONLY_COMPILE
  DEFAULT_ENABLE_IF ${FLOW_VARIANTS_DEFAULT_ENABLE_IF}
  SOURCES flow/flow_blackoil_dunecpr_seq.cpp
  EXE_NAME flow_blackoil_dunecpr_seq
  DEPENDS "opmsimulators"
  LIBRARIES "opmsimulators")


opm_add_test(flow_onephase
  ONLY_COMPILE
  DEFAULT_ENABLE_IF ${FLOW_VARIANTS_DEFAULT_ENABLE_IF}
  SOURCES
  flow/flow_onephase.cpp
  $<TARGET_OBJECTS:moduleVersion>
  EXE_NAME flow_onephase
  DEPENDS opmsimulators
  LIBRARIES opmsimulators)

opm_add_test(flow_onephase_energy
  ONLY_COMPILE
  DEFAULT_ENABLE_IF ${FLOW_VARIANTS_DEFAULT_ENABLE_IF}
  SOURCES
  flow/flow_onephase_energy.cpp
  $<TARGET_OBJECTS:moduleVersion>
  EXE_NAME flow_onephase_energy
  DEPENDS opmsimulators
  LIBRARIES opmsimulators)

if (BUILD_FLOW)
  install(TARGETS flow DESTINATION bin)
  opm_add_bash_completion(flow)

  add_test(NAME flow__version
    COMMAND flow --version)
  set_tests_properties(flow__version PROPERTIES
    PASS_REGULAR_EXPRESSION "${${project}_LABEL}")
endif()

if (NOT BUILD_EBOS)
  set(EBOS_DEFAULT_ENABLE_IF "FALSE")
else()
  set(EBOS_DEFAULT_ENABLE_IF "TRUE")
endif()

# the research oriented general-purpose ECL simulator ("ebos" == &ecl
# &black-&oil &simulator)
set(MEBOS_TARGETS "")
foreach(OBJ blackoil solvent polymer foam brine gasoil oilwater oilwaterpolymer thermal)
  add_library(ebos_lib${OBJ} OBJECT EXCLUDE_FROM_ALL ebos/ebos_${OBJ}.cc)
  list(APPEND MEBOS_TARGETS $<TARGET_OBJECTS:ebos_lib${OBJ}>)
endforeach()

opm_add_test(ebos
  ONLY_COMPILE
  DEFAULT_ENABLE_IF ${EBOS_DEFAULT_ENABLE_IF}
  ALWAYS_ENABLE
  EXE_NAME ebos
  LIBRARIES opmsimulators
  SOURCES ebos/ebos_main.cc $<TARGET_OBJECTS:ebos_libblackoil>)

if (BUILD_EBOS)
  install(TARGETS ebos DESTINATION bin)
  opm_add_bash_completion(ebos)
endif()

if (NOT BUILD_EBOS_EXTENSIONS)
  set(EBOS_EXTENSIONS_DEFAULT_ENABLE_IF "FALSE")
else()
  set(EBOS_EXTENSIONS_DEFAULT_ENABLE_IF "TRUE")
endif()

foreach(OBJ solvent polymer foam brine gasoil oilwater oilwaterpolymer thermal)
  opm_add_test(ebos_${OBJ}
    ONLY_COMPILE
    DEFAULT_ENABLE_IF ${EBOS_EXTENSIONS_DEFAULT_ENABLE_IF}
    SOURCES ebos/ebos_${OBJ}_main.cc $<TARGET_OBJECTS:ebos_lib${OBJ}>
    EXE_NAME ebos_${OBJ}
    LIBRARIES opmsimulators)
endforeach()

opm_add_test(mebos
  ONLY_COMPILE
  DEFAULT_ENABLE_IF ${EBOS_EXTENSIONS_DEFAULT_ENABLE_IF}
  SOURCES ebos/mebos_main.cc
          ${MEBOS_TARGETS}
  EXE_NAME mebos
  LIBRARIES opmsimulators)

if (NOT BUILD_EBOS_DEBUG_EXTENSIONS)
  set(EBOS_DEBUG_EXTENSIONS_DEFAULT_ENABLE_IF "FALSE")
else()
  set(EBOS_DEBUG_EXTENSIONS_DEFAULT_ENABLE_IF "TRUE")
endif()

opm_add_test(ebos_altidx
  ONLY_COMPILE
  DEFAULT_ENABLE_IF ${EBOS_DEBUG_EXTENSIONS_DEFAULT_ENABLE_IF}
  SOURCES ebos/ebos_altidx.cc
  EXE_NAME ebos_altidx
  DEPENDS opmsimulators
  LIBRARIES opmsimulators)

opm_add_test(ebos_plain
  ONLY_COMPILE
  DEFAULT_ENABLE_IF ${EBOS_DEBUG_EXTENSIONS_DEFAULT_ENABLE_IF}
  SOURCES ebos/ebos_plain.cc
  EXE_NAME ebos_plain
  DEPENDS opmsimulators
  LIBRARIES opmsimulators)

if (BUILD_EBOS_EXTENSIONS)
  foreach(TGT ebos_solvent ebos_polymer ebos_foam ebos_brine ebos_gasoil ebos_oilwater ebos_oilwaterpolymer ebos_thermal mebos)
    install(TARGETS ${TGT} DESTINATION bin)
    opm_add_bash_completion(${TGT})
  endforeach()
endif()

if (OPM_ENABLE_PYTHON)
  add_subdirectory(python)
endif()

add_custom_target(extra_test ${CMAKE_CTEST_COMMAND} -C ExtraTests)

# must link libraries after target 'flow' has been defined
if(CUDA_FOUND)
  target_link_libraries( opmsimulators PUBLIC ${CUDA_cusparse_LIBRARY} )
  target_link_libraries( opmsimulators PUBLIC ${CUDA_cublas_LIBRARY} )
endif()

if(OpenCL_FOUND)
  target_link_libraries( opmsimulators ${OpenCL_LIBRARIES} )
endif()

