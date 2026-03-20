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

# Slave parse error test.
add_custom_target(test_rc_slave_parsing_err
  COMMAND ${CMAKE_COMMAND} -E copy_directory
    ${CMAKE_CURRENT_SOURCE_DIR}/tests/rescoup/slave_parse_error
    ${CMAKE_CURRENT_BINARY_DIR}/tests/rescoup/slave_parse_error
  COMMAND ${CMAKE_COMMAND} -E copy_directory
    ${CMAKE_CURRENT_SOURCE_DIR}/tests/include
    ${CMAKE_CURRENT_BINARY_DIR}/tests/include
  COMMENT "Copying slave_parse_error test data to build tree"
  DEPENDS flow
)
add_test(NAME rc_slave_parsing_err
  COMMAND
    ${CMAKE_CURRENT_SOURCE_DIR}/tests/rescoup/slave_parse_error/run_ctest.sh
    $<TARGET_FILE:flow>
    ${_rescoup_mpiexec}
  WORKING_DIRECTORY
    ${CMAKE_CURRENT_BINARY_DIR}/tests/rescoup/slave_parse_error
  CONFIGURATIONS Integration
)
set_tests_properties(rc_slave_parsing_err PROPERTIES TIMEOUT 30)
