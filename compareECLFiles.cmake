# This script manages the addition of tests.
# The tests are orchestrated by a shell script,
# configured using opm_set_test_driver()
# and then the appropriate helper macro is called to
# register the ctest entry through the opm_add_test macro.
# Information such as the binary to call and test tolerances
# are passed from the build system to the driver script through
# command line parameters. See the opm_add_test() documentation for
# details on the parameters passed to the macro.

# Define some paths
set(BASE_RESULT_PATH ${PROJECT_BINARY_DIR}/tests/results)

###########################################################################
# TEST: compareECLFiles
###########################################################################

# Input:
#   - casename: basename (no extension)
#
# Details:
#   - This test class compares output from a simulation to reference files.
macro (add_test_compareECLFiles casename filename simulator abs_tol rel_tol prefix dirprefix)
  if(${ARGC} GREATER 7)
    set(DIR ${ARGV7})
  else()
    set(DIR ${casename})
  endif()
  if(${ARGC} GREATER 8)
    set(TEST_ARGS ${OPM_DATA_ROOT}/${DIR}/${ARGV1} deckfilename=${OPM_DATA_ROOT}/${DIR}/${filename})
  else()
    set(TEST_ARGS ${OPM_DATA_ROOT}/${DIR}/${filename})
  endif()
  set(RESULT_PATH ${BASE_RESULT_PATH}${dirprefix}/${simulator}+${casename})
  opm_add_test(${prefix}_${simulator}+${filename} NO_COMPILE
               EXE_NAME ${simulator}
               DRIVER_ARGS ${OPM_DATA_ROOT}/${DIR} ${RESULT_PATH}
                           ${CMAKE_BINARY_DIR}/bin
                           ${filename}
                           ${abs_tol} ${rel_tol}
                           ${COMPARE_SUMMARY_COMMAND}
                           ${COMPARE_ECL_COMMAND}
               TEST_ARGS ${TEST_ARGS})
endmacro (add_test_compareECLFiles)

###########################################################################
# TEST: add_test_compare_restarted_simulation
###########################################################################

# Input:
#   - casename: basename (no extension)
#
# Details:
#   - This test class compares the output from a restarted simulation
#     to that of a non-restarted simulation.
macro (add_test_compare_restarted_simulation casename filename simulator abs_tol rel_tol)

  set(RESULT_PATH ${BASE_RESULT_PATH}/restart/${simulator}+${casename})
  opm_add_test(compareRestartedSim_${simulator}+${filename} NO_COMPILE
               EXE_NAME ${simulator}
               DRIVER_ARGS ${OPM_DATA_ROOT}/${casename} ${RESULT_PATH}
                           ${CMAKE_BINARY_DIR}/bin
                           ${filename}
                           ${abs_tol} ${rel_tol}
                           ${COMPARE_SUMMARY_COMMAND}
                           ${COMPARE_ECL_COMMAND}
               TEST_ARGS ${OPM_DATA_ROOT}/${casename}/${filename})
endmacro (add_test_compare_restarted_simulation)

###########################################################################
# TEST: add_test_compare_parallel_simulation
###########################################################################

# Input:
#   - casename: basename (no extension)
#
# Details:
#   - This test class compares the output from a parallel simulation
#     to the output from the serial instance of the same model.
macro (add_test_compare_parallel_simulation casename filename simulator abs_tol rel_tol)
  set(RESULT_PATH ${BASE_RESULT_PATH}/parallel/${simulator}+${casename})

  # Add test that runs flow_mpi and outputs the results to file
  opm_add_test(compareParallelSim_${simulator}+${filename} NO_COMPILE
               EXE_NAME ${simulator}
               DRIVER_ARGS ${OPM_DATA_ROOT}/${casename} ${RESULT_PATH}
                           ${CMAKE_BINARY_DIR}/bin
                           ${filename}
                           ${abs_tol} ${rel_tol}
                           ${COMPARE_SUMMARY_COMMAND}
                           ${COMPARE_ECL_COMMAND}
               TEST_ARGS ${OPM_DATA_ROOT}/${casename}/${filename})
endmacro (add_test_compare_parallel_simulation)

if(NOT TARGET test-suite)
  add_custom_target(test-suite)
endif()

# Regression tests
opm_set_test_driver(${PROJECT_SOURCE_DIR}/tests/run-regressionTest.sh "")

# Set absolute tolerance to be used passed to the macros in the following tests
set(abs_tol 2e-2)
set(rel_tol 1e-5)

add_test_compareECLFiles(spe1 SPE1CASE2 flow_ebos ${abs_tol} ${rel_tol} compareECLFiles "")
add_test_compareECLFiles(spe1_2p SPE1CASE2_2P flow_ebos ${abs_tol} ${rel_tol} compareECLFiles "" spe1)
add_test_compareECLFiles(spe1 SPE1CASE2 flow_legacy ${abs_tol} ${rel_tol} compareECLFiles "")
add_test_compareECLFiles(spe1_2p SPE1CASE2_2P flow_legacy ${abs_tol} ${rel_tol} compareECLFiles "" spe1)
add_test_compareECLFiles(spe1 SPE1CASE1 flow_sequential ${abs_tol} ${rel_tol} compareECLFiles "")
add_test_compareECLFiles(spe3 SPE3CASE1 flow_ebos ${abs_tol} ${rel_tol} compareECLFiles "")
add_test_compareECLFiles(spe3 SPE3CASE1 flow_legacy ${abs_tol} ${rel_tol} compareECLFiles "")
add_test_compareECLFiles(spe9 SPE9_CP_SHORT flow_ebos ${abs_tol} ${rel_tol} compareECLFiles "")
add_test_compareECLFiles(spe9 SPE9_CP_SHORT flow_legacy ${abs_tol} ${rel_tol} compareECLFiles "")
add_test_compareECLFiles(msw_2d_h 2D_H__ flow_multisegment ${abs_tol} ${rel_tol} compareECLFiles "")
add_test_compareECLFiles(polymer_simple2D 2D_THREEPHASE_POLY_HETER flow_polymer ${abs_tol} ${rel_tol} compareECLFiles "" polymer_simple2D run.param)
add_test_compareECLFiles(spe5 SPE5CASE1 flow_solvent ${abs_tol} ${rel_tol} compareECLFiles "")

# Restart tests
opm_set_test_driver(${PROJECT_SOURCE_DIR}/tests/run-restart-regressionTest.sh "")

# Cruder tolerances for the restarted tests
set(abs_tol_restart 2e-1)
set(rel_tol_restart 4e-5)
foreach(sim flow_legacy flow_ebos)
  add_test_compare_restarted_simulation(spe1 SPE1CASE2_ACTNUM ${sim} ${abs_tol_restart} ${rel_tol_restart})
  add_test_compare_restarted_simulation(spe9 SPE9_CP_SHORT ${sim} ${abs_tol_restart} ${rel_tol_restart})
endforeach()

# Init tests
opm_set_test_driver(${PROJECT_SOURCE_DIR}/tests/run-init-regressionTest.sh "")

foreach(sim flow_legacy flow_ebos)
  add_test_compareECLFiles(norne NORNE_ATW2013 ${sim} ${abs_tol} ${rel_tol} compareECLInitFiles /init)
endforeach()

# Parallel tests
if(MPI_FOUND)
  opm_set_test_driver(${PROJECT_SOURCE_DIR}/tests/run-parallel-regressionTest.sh "")

  # Different tolerances for these tests
  set(abs_tol_parallel 0.02)
  set(rel_tol_parallel 1e-5)

  foreach(sim flow_mpi flow_ebos)
    add_test_compare_parallel_simulation(spe1 SPE1CASE2 ${sim} ${abs_tol_parallel} ${rel_tol_parallel})
    add_test_compare_parallel_simulation(spe9 SPE9_CP_SHORT ${sim} ${abs_tol_parallel} ${rel_tol_parallel})
  endforeach()
  add_test_compare_parallel_simulation(spe3 SPE3CASE1 flow_mpi ${abs_tol_parallel} ${rel_tol_parallel})
endif()
