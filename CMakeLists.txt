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

cmake_minimum_required (VERSION 3.10)


option(SIBLING_SEARCH "Search for other modules in sibling directories?" ON)
set( USE_OPENMP_DEFAULT OFF ) # Use of OpenMP is considered experimental
option(BUILD_FLOW "Build the production oriented flow simulator?" ON)
option(BUILD_FLOW_VARIANTS "Build the variants for flow by default?" OFF)
option(BUILD_FLOW_POLY_GRID "Build flow blackoil with polyhedral grid" OFF)
option(OPM_ENABLE_PYTHON "Enable python bindings?" OFF)
option(OPM_ENABLE_PYTHON_TESTS "Enable tests for the python bindings?" ON)
option(OPM_INSTALL_PYTHON "Install python bindings?" ON)
option(USE_CHOW_PATEL_ILU "Use the iterative ILU by Chow and Patel?" OFF)
option(USE_CHOW_PATEL_ILU_GPU "Run iterative ILU decomposition on GPU? Requires USE_CHOW_PATEL_ILU" OFF)
option(USE_CHOW_PATEL_ILU_GPU_PARALLEL "Try to use more parallelism on the GPU during the iterative ILU decomposition? Requires USE_CHOW_PATEL_ILU_GPU" OFF)
option(BUILD_FLOW_ALU_GRID "Build flow blackoil with alu grid" OFF)
option(USE_DAMARIS_LIB "Use the Damaris library for asynchronous I/O?" OFF)

# The following was copied from CMakeLists.txt in opm-common.
# TODO: factor out the common parts in opm-common and opm-simulator as a cmake module
if (OPM_ENABLE_PYTHON)
  # We need to be compatible with older CMake versions
  # that do not offer FindPython3
  # e.g. Ubuntu LTS 18.04 uses cmake 3.10
  if(${CMAKE_VERSION} VERSION_LESS "3.12.0")
    find_package(PythonInterp REQUIRED)
    if(PYTHON_VERSION_MAJOR LESS 3)
      message(SEND_ERROR "OPM requires version 3 of Python but only version ${PYTHON_VERSION_STRING} was found")
    endif()
    set(Python3_EXECUTABLE ${PYTHON_EXECUTABLE})
    set(Python3_LIBRARIES ${PYTHON_LIBRARIES})
    set(Python3_VERSION "${PYTHON_VERSION_STRING}")
    set(Python3_VERSION_MINOR ${PYTHON_VERSION_MINOR})
  else()
    # Be backwards compatible.
    if(PYTHON_EXECUTABLE AND NOT Python3_EXECUTABLE)
      set(Python3_EXECUTABLE ${PYTHON_EXECUTABLE})
    endif()
    find_package(Python3 REQUIRED COMPONENTS Interpreter Development)
  endif()
  if(Python3_VERSION_MINOR LESS 3)
    # Python native namespace packages requires python >= 3.3
    message(SEND_ERROR "OPM requires python >= 3.3 but only version ${Python3_VERSION} was found")
  endif()
  # Compatibility settings for PythonInterp and PythonLibs
  # used e.g. in FindCwrap, pybind11
  set(PYTHON_EXECUTABLE ${Python3_EXECUTABLE})
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
# CMAKE_VERSION VERSION_GREATER_EQUAL 3.8, because of backwards
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
  set(COMPILE_BDA_BRIDGE 1)
  include_directories(${CUDA_INCLUDE_DIRS})
endif()

find_package(OpenCL)

if(OpenCL_FOUND)
  # the current OpenCL implementation relies on cl2.hpp, not cl.hpp
  # make sure it is available, otherwise disable OpenCL
  find_file(CL2_HPP CL/cl2.hpp HINTS ${OpenCL_INCLUDE_DIRS})
  if(CL2_HPP)
    set(HAVE_OPENCL 1)
    set(COMPILE_BDA_BRIDGE 1)
    include_directories(${OpenCL_INCLUDE_DIRS})
    find_file(OPENCL_HPP CL/opencl.hpp HINTS ${OpenCL_INCLUDE_DIRS})
    if(OPENCL_HPP)
      set(HAVE_OPENCL_HPP 1)
    endif()
  else()
    message(WARNING " OpenCL was found, but this version of opm-simulators relies on CL/cl2.hpp, which implements OpenCL 1.0, 1.1 and 1.2.\n Deactivating OpenCL")
    set(OpenCL_FOUND OFF)
    set(OPENCL_FOUND OFF)
  endif()
  if(USE_CHOW_PATEL_ILU)
    add_compile_options(-DCHOW_PATEL=1)
    if(USE_CHOW_PATEL_ILU_GPU)
      add_compile_options(-DCHOW_PATEL_GPU=1)
      if(USE_CHOW_PATEL_ILU_GPU_PARALLEL)
        add_compile_options(-DCHOW_PATEL_GPU_PARALLEL=1)
      else()
        add_compile_options(-DCHOW_PATEL_GPU_PARALLEL=0)
      endif()
    else()
      add_compile_options(-DCHOW_PATEL_GPU=0)
      add_compile_options(-DCHOW_PATEL_GPU_PARALLEL=0)
    endif()
  endif()
else()
  if(USE_CHOW_PATEL_ILU)
    message(FATAL_ERROR " CHOW_PATEL_ILU only works for openclSolver, but OpenCL was not found")
  endif()
endif()

find_package(amgcl)
if(amgcl_FOUND)
  set(HAVE_AMGCL 1)
  set(COMPILE_BDA_BRIDGE 1)
  # Linking to target angcl::amgcl drags in OpenMP and -fopenmp as a compile
  # flag. With that nvcc fails as it does not that flag.
  # Hence we set AMGCL_INCLUDE_DIRS.
  get_property(AMGCL_INCLUDE_DIRS TARGET amgcl::amgcl PROPERTY INTERFACE_INCLUDE_DIRECTORIES)
  include_directories(SYSTEM ${AMGCL_INCLUDE_DIRS})
endif()

if(OpenCL_FOUND)
  find_package(VexCL)
  if(VexCL_FOUND)
    set(HAVE_VEXCL 1)
    set(COMPILE_BDA_BRIDGE 1)
    # generator expressions in vexcl do not seem to work and therefore
    # we cannot use the imported target. Hence we exract the needed info
    # from the targets
    get_property(VEXCL_INCLUDE_DIRS TARGET VexCL::Common PROPERTY INTERFACE_INCLUDE_DIRECTORIES)
    get_property(VEXCL_LINK_LIBRARIES TARGET VexCL::Common PROPERTY INTERFACE_LINK_LIBRARIES)
    get_property(VEXCL_COMPILE_DEFINITIONS TARGET VexCL::OpenCL PROPERTY INTERFACE_COMPILE_DEFINITIONS)
    set(VEXCL_LINK_LIBRARIES "${VEXCL_LINK_LIBRARIES};OpenCL::OpenCL")
    add_library(OPM::VexCL::OpenCL INTERFACE IMPORTED)
    set_target_properties(OPM::VexCL::OpenCL PROPERTIES
      INTERFACE_COMPILE_DEFINITIONS "${VEXCL_COMPILE_DEFINITIONS}"
      INTERFACE_LINK_LIBRARIES "${VEXCL_LINK_LIBRARIES}")
    target_include_directories(OPM::VexCL::OpenCL SYSTEM INTERFACE "${VEXCL_INCLUDE_DIRS}")
  endif()
endif()

find_package(rocalution)
if(ROCALUTION_FOUND)
  set(HAVE_ROCALUTION 1)
  set(COMPILE_BDA_BRIDGE 1)
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

  if(NOT fmt_FOUND)
    add_definitions(-DFMT_HEADER_ONLY)
    list(APPEND EXTRA_INCLUDES SYSTEM ${PROJECT_SOURCE_DIR}/external/fmtlib/include)
  endif()
  include_directories(${EXTRA_INCLUDES})

  include(UseDamaris)
endmacro (config_hook)

macro (prereqs_hook)
endmacro (prereqs_hook)

macro (sources_hook)
  if(OPENCL_FOUND)
    include(opencl-source-provider)
    list(APPEND opm-simulators_SOURCES ${PROJECT_BINARY_DIR}/clSources.cpp)
  endif()
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
    -n 4
    -b ${PROJECT_BINARY_DIR}
)

opm_add_test(test_gatherdeferredlogger
  DEPENDS "opmsimulators"
  LIBRARIES opmsimulators ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY}
  SOURCES
    tests/test_gatherdeferredlogger.cpp
  CONDITION
    MPI_FOUND AND Boost_UNIT_TEST_FRAMEWORK_FOUND
  DRIVER_ARGS
    -n 4
    -b ${PROJECT_BINARY_DIR}
)

opm_add_test(test_parallelwellinfo_mpi
  EXE_NAME
    test_parallelwellinfo
  CONDITION
    MPI_FOUND AND Boost_UNIT_TEST_FRAMEWORK_FOUND
  DRIVER_ARGS
    -n 4
    -b ${PROJECT_BINARY_DIR}
  NO_COMPILE
    )

opm_add_test(test_broadcast
  DEPENDS "opmsimulators"
  LIBRARIES opmsimulators ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY}
  SOURCES
    tests/test_broadcast.cpp
  CONDITION
    MPI_FOUND AND Boost_UNIT_TEST_FRAMEWORK_FOUND
  DRIVER_ARGS
    -n 4
    -b ${PROJECT_BINARY_DIR}
)

include(OpmBashCompletion)

if (NOT BUILD_FLOW)
  set(FLOW_DEFAULT_ENABLE_IF "FALSE")
else()
  set(FLOW_DEFAULT_ENABLE_IF "TRUE")
endif()

if (NOT BUILD_FLOW_VARIANTS)
  set(FLOW_VARIANTS_DEFAULT_ENABLE_IF "FALSE")
else()
  set(FLOW_VARIANTS_DEFAULT_ENABLE_IF "TRUE")
endif()

if (NOT BUILD_FLOW_POLY_GRID)
  set(FLOW_POLY_ONLY_DEFAULT_ENABLE_IF "FALSE")
else()
  set(FLOW_POLY_ONLY_DEFAULT_ENABLE_IF "TRUE")
endif()

add_library(moduleVersion OBJECT opm/simulators/utils/moduleVersion.cpp)
set_property(TARGET moduleVersion PROPERTY POSITION_INDEPENDENT_CODE ON)

# Strictly we only depend on the update-version target,
# but this is not exposed in a super-build.
add_dependencies(moduleVersion opmsimulators)

set(FLOW_MODELS blackoil brine energy extbo foam gasoil gaswater
                oilwater oilwater_brine gaswater_brine oilwater_polymer
                oilwater_polymer_injectivity micp polymer solvent
                gasoil_energy brine_saltprecipitation
                gaswater_saltprec_vapwat brine_precsalt_vapwat
                blackoil_legacyassembly gasoildiffuse)
set(FLOW_VARIANT_MODELS brine_energy onephase onephase_energy)

set(FLOW_TGTS)
foreach(OBJ ${COMMON_MODELS} ${FLOW_MODELS} ${FLOW_VARIANT_MODELS})
  add_library(flow_lib${OBJ} OBJECT flow/flow_ebos_${OBJ}.cpp)
  list(APPEND FLOW_TGTS $<TARGET_OBJECTS:flow_lib${OBJ}>)
  if(TARGET fmt::fmt)
    target_link_libraries(flow_lib${OBJ} fmt::fmt)
  endif()
  opm_add_test(flow_${OBJ}
               ONLY_COMPILE
               SOURCES
               flow/flow_${OBJ}.cpp
               $<TARGET_OBJECTS:moduleVersion>
               $<TARGET_OBJECTS:flow_lib${OBJ}>
               EXE_NAME flow_${OBJ}
               DEPENDS opmsimulators
               LIBRARIES opmsimulators)
endforeach()
set_property(TARGET flow_libblackoil PROPERTY POSITION_INDEPENDENT_CODE ON)

foreach(OBJ ${FLOW_VARIANT_MODELS})
  set_property(TARGET flow_lib${OBJ} PROPERTY EXCLUDE_FROM_ALL ${FLOW_VARIANTS_DEFAULT_ENABLE_IF})
endforeach()

add_library(flow_libblackoil_poly OBJECT flow/flow_ebos_blackoil_poly.cpp)
if(TARGET fmt::fmt)
  target_link_libraries(flow_libblackoil_poly fmt::fmt)
 endif()
set_property(TARGET flow_libblackoil_poly PROPERTY EXCLUDE_FROM_ALL ${FLOW_POLY_ONLY_DEFAULT_ENABLE_IF})

opm_add_test(flow
  ONLY_COMPILE
  ALWAYS_ENABLE
  DEFAULT_ENABLE_IF ${FLOW_DEFAULT_ENABLE_IF}
  DEPENDS opmsimulators
  LIBRARIES opmsimulators
  SOURCES
  flow/flow.cpp
  ${FLOW_TGTS}
  $<TARGET_OBJECTS:moduleVersion>
  )

opm_add_test(flow_poly
  ONLY_COMPILE
  ALWAYS_ENABLE
  DEFAULT_ENABLE_IF ${FLOW_POLY_ONLY_DEFAULT_ENABLE_IF}
  DEPENDS opmsimulators
  LIBRARIES opmsimulators
  SOURCES
  flow/flow_blackoil_poly.cpp
  $<TARGET_OBJECTS:flow_libblackoil_poly>
  $<TARGET_OBJECTS:moduleVersion>)
target_compile_definitions(flow_poly PRIVATE USE_POLYHEDRALGRID)

opm_add_test(flow_distribute_z
  ONLY_COMPILE
  ALWAYS_ENABLE
  DEFAULT_ENABLE_IF ${FLOW_DEFAULT_ENABLE_IF}
  DEPENDS opmsimulators
  LIBRARIES opmsimulators
  SOURCES
  flow/flow_distribute_z.cpp
  ${FLOW_TGTS}
  $<TARGET_OBJECTS:moduleVersion>
  )

if(dune-alugrid_FOUND)
  if (NOT BUILD_FLOW_ALU_GRID)
    set(FLOW_ALUGRID_ONLY_DEFAULT_ENABLE_IF "FALSE")
  else()
    set(FLOW_ALUGRID_ONLY_DEFAULT_ENABLE_IF "TRUE")
  endif()

  opm_add_test(flow_alugrid
    ONLY_COMPILE
    ALWAYS_ENABLE
    DEFAULT_ENABLE_IF ${FLOW_ALUGRID_ONLY_DEFAULT_ENABLE_IF}
    DEPENDS opmsimulators
    LIBRARIES opmsimulators
    SOURCES
    flow/flow_blackoil_alugrid.cpp
    $<TARGET_OBJECTS:moduleVersion>)
  target_compile_definitions(flow_alugrid PRIVATE USE_ALUGRID)
endif()

if (BUILD_FLOW)
  install(TARGETS flow DESTINATION bin)
  opm_add_bash_completion(flow)

  add_test(NAME flow__version
    COMMAND flow --version)
  set_tests_properties(flow__version PROPERTIES
    PASS_REGULAR_EXPRESSION "${${project}_LABEL}")
endif()

if (OPM_ENABLE_PYTHON)
  add_subdirectory(python)
endif()

add_custom_target(extra_test ${CMAKE_CTEST_COMMAND} -C ExtraTests)

# must link libraries after target 'opmsimulators' has been defined

if(CUDA_FOUND)
  target_link_libraries( opmsimulators PUBLIC ${CUDA_cusparse_LIBRARY} )
  target_link_libraries( opmsimulators PUBLIC ${CUDA_cublas_LIBRARY} )
endif()

if(OpenCL_FOUND)
  target_link_libraries( opmsimulators PUBLIC ${OpenCL_LIBRARIES} )
endif()

if(ROCALUTION_FOUND)
  target_include_directories(opmsimulators PUBLIC ${rocalution_INCLUDE_DIR}/rocalution)
endif()

if(VexCL_FOUND)
  target_link_libraries( opmsimulators PUBLIC OPM::VexCL::OpenCL )
endif()

if(Damaris_FOUND)
  target_link_libraries(opmsimulators PUBLIC damaris)
endif()

install(DIRECTORY doc/man1 DESTINATION ${CMAKE_INSTALL_MANDIR}
  FILES_MATCHING PATTERN "*.1")


