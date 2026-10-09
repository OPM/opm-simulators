# SPDX-License-Identifier: GPL-3.0-or-later
# Standalone RFC prototype. Not included by the production build.
cmake_minimum_required(VERSION 3.23)
project(opm_vulkan_experiment LANGUAGES C CXX)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

option(BUILD_OPM_ADAPTER "Build the experimental Flow executables against installed OPM" ON)
option(OPM_VULKAN_PACKAGE_ASSERTIONS "Installed OPM was built with WITH_NDEBUG=OFF" ON)
option(OPM_VULKAN_GPU_TESTS "Register tests that require a hardware Vulkan device" OFF)
include(CTest)

find_package(Vulkan REQUIRED)
find_package(Python3 REQUIRED COMPONENTS Interpreter)
find_program(GLSLANG glslangValidator REQUIRED)

set(shader "${CMAKE_CURRENT_SOURCE_DIR}/bicgstab.comp")
foreach(variant IN ITEMS bicgstab bicgstab_dilu)
  set(defines -DRESLAB_SPECIALIZE=1)
  if(variant STREQUAL "bicgstab_dilu")
    list(APPEND defines -DRESLAB_DILU=1)
  endif()
  add_custom_command(
    OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/${variant}_spv.hpp"
    COMMAND "${GLSLANG}" -V ${defines} "${shader}"
            -o "${CMAKE_CURRENT_BINARY_DIR}/${variant}.spv"
    COMMAND "${Python3_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/embed_shader.py"
            "${CMAKE_CURRENT_BINARY_DIR}/${variant}.spv"
            "${CMAKE_CURRENT_BINARY_DIR}/${variant}_spv.hpp"
    DEPENDS "${shader}" embed_shader.py
    VERBATIM)
endforeach()

add_library(opm_vulkan_experiment VulkanSolver.cpp
  "${CMAKE_CURRENT_BINARY_DIR}/bicgstab_spv.hpp"
  "${CMAKE_CURRENT_BINARY_DIR}/bicgstab_dilu_spv.hpp")
target_include_directories(opm_vulkan_experiment
  PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}" PRIVATE "${CMAKE_CURRENT_BINARY_DIR}")
target_link_libraries(opm_vulkan_experiment PUBLIC Vulkan::Vulkan)
target_compile_options(opm_vulkan_experiment PRIVATE -Wall -Wextra -Wno-missing-field-initializers)

if(BUILD_TESTING)
  add_executable(test_vulkan_native test_native.cpp)
  target_link_libraries(test_vulkan_native PRIVATE opm_vulkan_experiment)
  if(OPM_VULKAN_GPU_TESTS)
    foreach(mode IN ITEMS block-jacobi dilu dilu-sweeps cpr)
      add_test(NAME vulkan_${mode} COMMAND test_vulkan_native --${mode})
      set_tests_properties(vulkan_${mode} PROPERTIES LABELS "gpu;vulkan" RUN_SERIAL TRUE)
    endforeach()
    add_test(NAME vulkan_cpr_cooperative COMMAND test_vulkan_native --cpr)
    set_tests_properties(vulkan_cpr_cooperative PROPERTIES
      LABELS "gpu;vulkan" RUN_SERIAL TRUE
      ENVIRONMENT "RESLAB_VULKAN_COOPERATIVE=1;RESLAB_VULKAN_AGGREGATION_STRENGTH=0.25")
    add_test(NAME vulkan_cpr_fpf COMMAND test_vulkan_native --cpr)
    set_tests_properties(vulkan_cpr_fpf PROPERTIES
      LABELS "gpu;vulkan" RUN_SERIAL TRUE
      ENVIRONMENT "RESLAB_VULKAN_COOPERATIVE=1;RESLAB_VULKAN_AGGREGATION_STRENGTH=0.25;RESLAB_VULKAN_CPR_FPF=1;RESLAB_VULKAN_SWEEP_COUNT=3;RESLAB_VULKAN_PRESSURE_CYCLES=2")
  endif()
endif()

if(BUILD_OPM_ADAPTER)
  find_package(opm-simulators REQUIRED)
  foreach(target IN ITEMS flow_vulkan flow_vulkan_cpu)
    add_executable(${target} flow_vulkan.cpp moduleVersion.cpp)
    target_link_libraries(${target} PRIVATE opm_vulkan_experiment opmsimulators)
    # EnsureFinalized changes object layout with NDEBUG: match the installed ABI.
    if(OPM_VULKAN_PACKAGE_ASSERTIONS)
      target_compile_options(${target} PRIVATE -UNDEBUG)
    endif()
  endforeach()
  target_compile_definitions(flow_vulkan_cpu PRIVATE RESLAB_CPU_REFERENCE=1)
endif()
