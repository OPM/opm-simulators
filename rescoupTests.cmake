# Copyright 2026 Equinor
#
# This file is part of the Open Porous Media project (OPM).
#
# OPM is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# OPM is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with OPM.  If not, see <http://www.gnu.org/licenses/>.

# Reservoir coupling integration tests.
# Included from CMakeLists.txt inside if(MPI_FOUND).

# MPIEXEC_EXECUTABLE may point to a different MPI implementation than
# the one the binaries were compiled against. This happens when CMake's
# FindMPI locates MPIEXEC_EXECUTABLE via PATH independently of
# MPI_C_COMPILER (e.g., system OpenMPI found first when compiled against
# a custom MPICH).
#
# Resolution order:
#   1. If MPIEXEC_EXECUTABLE is in the same directory as MPI_C_COMPILER,
#      it was either explicitly set by the user (-DMPIEXEC_EXECUTABLE=...)
#      or correctly auto-detected — use it directly.
#   2. Otherwise, search for a matching launcher in MPI_C_COMPILER's bin
#      directory (workaround for PATH-based detection picking wrong impl).
#   3. Fall back to MPIEXEC_EXECUTABLE if nothing else is found.
get_filename_component(_rescoup_mpi_bin_dir "${MPI_C_COMPILER}" DIRECTORY)
get_filename_component(_rescoup_mpiexec_dir "${MPIEXEC_EXECUTABLE}" DIRECTORY)
if("${_rescoup_mpiexec_dir}" STREQUAL "${_rescoup_mpi_bin_dir}")
  set(_rescoup_mpiexec "${MPIEXEC_EXECUTABLE}")
else()
  find_program(_rescoup_mpiexec
    NAMES mpiexec.hydra mpiexec mpirun
    HINTS "${_rescoup_mpi_bin_dir}"
    NO_DEFAULT_PATH
  )
  if(NOT _rescoup_mpiexec)
    set(_rescoup_mpiexec "${MPIEXEC_EXECUTABLE}")
  endif()
endif()

# Select the simulator binary. When USE_DEV_SIMULATOR_IN_TESTS is set, use the
# faster-to-build development simulator (flow_blackoil) instead of the full flow
# target. This matches the pattern used by the other test cmake files.
if(USE_DEV_SIMULATOR_IN_TESTS)
  set(_rescoup_simulator flow_blackoil)
else()
  set(_rescoup_simulator flow)
endif()

# Slave parse error test.
add_custom_target(test_rc_slave_parsing_err
  COMMAND ${CMAKE_COMMAND} -E copy_directory
    ${CMAKE_CURRENT_SOURCE_DIR}/tests/rescoup/slave_parse_error
    ${CMAKE_CURRENT_BINARY_DIR}/tests/rescoup/slave_parse_error
  COMMAND ${CMAKE_COMMAND} -E copy_directory
    ${CMAKE_CURRENT_SOURCE_DIR}/tests/include
    ${CMAKE_CURRENT_BINARY_DIR}/tests/include
  COMMENT "Copying slave_parse_error test data to build tree"
  DEPENDS ${_rescoup_simulator}
)
add_test(NAME rc_slave_parsing_err
  COMMAND
    ${CMAKE_CURRENT_SOURCE_DIR}/tests/rescoup/slave_parse_error/run_ctest.sh
    $<TARGET_FILE:${_rescoup_simulator}>
    ${_rescoup_mpiexec}
  WORKING_DIRECTORY
    ${CMAKE_CURRENT_BINARY_DIR}/tests/rescoup/slave_parse_error
  CONFIGURATIONS Integration
)
set_tests_properties(rc_slave_parsing_err PROPERTIES TIMEOUT 30)

# Multi-slave parse error test: one slave fails to parse while another is
# healthy. Verifies the master aborts the healthy slave cleanly (exit 1) instead
# of deadlocking on a collective MPI_Comm_disconnect.
add_custom_target(test_rc_slave_parsing_err_2slaves
  COMMAND ${CMAKE_COMMAND} -E copy_directory
    ${CMAKE_CURRENT_SOURCE_DIR}/tests/rescoup/slave_parse_error_2slaves
    ${CMAKE_CURRENT_BINARY_DIR}/tests/rescoup/slave_parse_error_2slaves
  COMMAND ${CMAKE_COMMAND} -E copy_directory
    ${CMAKE_CURRENT_SOURCE_DIR}/tests/include
    ${CMAKE_CURRENT_BINARY_DIR}/tests/include
  COMMENT "Copying slave_parse_error_2slaves test data to build tree"
  DEPENDS ${_rescoup_simulator}
)
add_test(NAME rc_slave_parsing_err_2slaves
  COMMAND
    ${CMAKE_CURRENT_SOURCE_DIR}/tests/rescoup/slave_parse_error_2slaves/run_ctest.sh
    $<TARGET_FILE:${_rescoup_simulator}>
    ${_rescoup_mpiexec}
  WORKING_DIRECTORY
    ${CMAKE_CURRENT_BINARY_DIR}/tests/rescoup/slave_parse_error_2slaves
  CONFIGURATIONS Integration
)
set_tests_properties(rc_slave_parsing_err_2slaves PROPERTIES TIMEOUT 40)
