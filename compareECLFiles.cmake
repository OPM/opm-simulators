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
function(add_test_compareECLFiles)
  set(oneValueArgs CASENAME FILENAME SIMULATOR ABS_TOL REL_TOL DIR DIR_PREFIX PREFIX)
  set(multiValueArgs TEST_ARGS)
  cmake_parse_arguments(PARAM "$" "${oneValueArgs}" "${multiValueArgs}" ${ARGN} )
  if(NOT PARAM_DIR)
    set(PARAM_DIR ${PARAM_CASENAME})
  endif()
  if(NOT PARAM_PREFIX)
    set(PARAM_PREFIX compareECLFiles)
  endif()
  set(RESULT_PATH ${BASE_RESULT_PATH}${PARAM_DIR_PREFIX}/${PARAM_SIMULATOR}+${PARAM_CASENAME})
  set(TEST_ARGS ${OPM_DATA_ROOT}/${PARAM_DIR}/${PARAM_FILENAME} ${PARAM_TEST_ARGS})
  opm_add_test(${PARAM_PREFIX}_${PARAM_SIMULATOR}+${PARAM_FILENAME} NO_COMPILE
               EXE_NAME ${PARAM_SIMULATOR}
               DRIVER_ARGS ${OPM_DATA_ROOT}/${PARAM_DIR} ${RESULT_PATH}
                           ${CMAKE_BINARY_DIR}/bin
                           ${PARAM_FILENAME}
                           ${PARAM_ABS_TOL} ${PARAM_REL_TOL}
                           ${COMPARE_SUMMARY_COMMAND}
                           ${COMPARE_ECL_COMMAND}
               TEST_ARGS ${TEST_ARGS})
endfunction()

###########################################################################
# TEST: add_test_compare_restarted_simulation
###########################################################################

# Input:
#   - casename: basename (no extension)
#
# Details:
#   - This test class compares the output from a restarted simulation
#     to that of a non-restarted simulation.
function(add_test_compare_restarted_simulation)
  set(oneValueArgs CASENAME FILENAME SIMULATOR ABS_TOL REL_TOL)
  set(multiValueArgs TEST_ARGS)
  cmake_parse_arguments(PARAM "$" "${oneValueArgs}" "${multiValueArgs}" ${ARGN} )

  set(RESULT_PATH ${BASE_RESULT_PATH}/restart/${PARAM_SIMULATOR}+${PARAM_CASENAME})
  set(TEST_ARGS ${OPM_DATA_ROOT}/${PARAM_CASENAME}/${PARAM_FILENAME} ${PARAM_TEST_ARGS})

  opm_add_test(compareRestartedSim_${PARAM_SIMULATOR}+${PARAM_FILENAME} NO_COMPILE
               EXE_NAME ${PARAM_SIMULATOR}
               DRIVER_ARGS ${OPM_DATA_ROOT}/${PARAM_CASENAME} ${RESULT_PATH}
                           ${CMAKE_BINARY_DIR}/bin
                           ${PARAM_FILENAME}
                           ${PARAM_ABS_TOL} ${PARAM_REL_TOL}
                           ${COMPARE_SUMMARY_COMMAND}
                           ${COMPARE_ECL_COMMAND}
               TEST_ARGS ${TEST_ARGS})
endfunction()

###########################################################################
# TEST: add_test_compare_parallel_simulation
###########################################################################

# Input:
#   - casename: basename (no extension)
#
# Details:
#   - This test class compares the output from a parallel simulation
#     to the output from the serial instance of the same model.
function(add_test_compare_parallel_simulation)
  set(oneValueArgs CASENAME FILENAME SIMULATOR ABS_TOL REL_TOL)
  set(multiValueArgs TEST_ARGS)
  cmake_parse_arguments(PARAM "$" "${oneValueArgs}" "${multiValueArgs}" ${ARGN} )

  set(RESULT_PATH ${BASE_RESULT_PATH}/parallel/${PARAM_SIMULATOR}+${PARAM_CASENAME})
  set(TEST_ARGS ${OPM_DATA_ROOT}/${PARAM_CASENAME}/${PARAM_FILENAME} ${PARAM_TEST_ARGS})

  # Add test that runs flow_mpi and outputs the results to file
  opm_add_test(compareParallelSim_${PARAM_SIMULATOR}+${PARAM_FILENAME} NO_COMPILE
               EXE_NAME ${PARAM_SIMULATOR}
               DRIVER_ARGS ${OPM_DATA_ROOT}/${PARAM_CASENAME} ${RESULT_PATH}
                           ${CMAKE_BINARY_DIR}/bin
                           ${PARAM_FILENAME}
                           ${PARAM_ABS_TOL} ${PARAM_REL_TOL}
                           ${COMPARE_SUMMARY_COMMAND}
                           ${COMPARE_ECL_COMMAND}
               TEST_ARGS ${TEST_ARGS})
endfunction()

if(NOT TARGET test-suite)
  add_custom_target(test-suite)
endif()

# Regression tests
opm_set_test_driver(${PROJECT_SOURCE_DIR}/tests/run-regressionTest.sh "")

# Set absolute tolerance to be used passed to the macros in the following tests
set(abs_tol 2e-2)
set(rel_tol 1e-5)
set(coarse_rel_tol 1e-2)

foreach(SIM flow flow_ebos flow_legacy)
  add_test_compareECLFiles(CASENAME spe1
                           FILENAME SPE1CASE2
                           SIMULATOR ${SIM}
                           ABS_TOL ${abs_tol}
                           REL_TOL ${coarse_rel_tol})
endforeach()

add_test_compareECLFiles(CASENAME spe1_2p
                         FILENAME SPE1CASE2_2P
                         SIMULATOR flow
                         ABS_TOL ${abs_tol}
                         REL_TOL ${rel_tol}
                         DIR spe1)

add_test_compareECLFiles(CASENAME spe1_2p
                         FILENAME SPE1CASE2_2P
                         SIMULATOR flow_legacy
                         ABS_TOL ${abs_tol}
                         REL_TOL ${coarse_rel_tol}
                         DIR spe1)

add_test_compareECLFiles(CASENAME spe1
                         FILENAME SPE1CASE1
                         SIMULATOR flow_sequential
                         ABS_TOL ${abs_tol}
                         REL_TOL ${rel_tol})

foreach(SIM flow flow_ebos flow_legacy)
  add_test_compareECLFiles(CASENAME spe3
                           FILENAME SPE3CASE1
                           SIMULATOR ${SIM}
                           ABS_TOL ${abs_tol}
                           REL_TOL ${coarse_rel_tol}
                           TEST_ARGS tolerance_wells=1e-6 max_iter=20)
endforeach()

foreach(SIM flow flow_ebos flow_legacy)
  add_test_compareECLFiles(CASENAME spe9
                           FILENAME SPE9_CP_SHORT
                           SIMULATOR ${SIM}
                           ABS_TOL ${abs_tol}
                           REL_TOL ${rel_tol})
endforeach()

foreach(SIM flow flow_ebos)
  add_test_compareECLFiles(CASENAME spe9group
                           FILENAME SPE9_CP_GROUP
                           SIMULATOR ${SIM}
                           ABS_TOL ${abs_tol}
                           REL_TOL ${rel_tol})
endforeach()

add_test_compareECLFiles(CASENAME msw_2d_h
                         FILENAME 2D_H__
                         SIMULATOR flow_multisegment
                         ABS_TOL ${abs_tol}
                         REL_TOL ${rel_tol})

add_test_compareECLFiles(CASENAME polymer_simple2D
                         FILENAME 2D_THREEPHASE_POLY_HETER
                         SIMULATOR flow_polymer
                         ABS_TOL ${abs_tol}
                         REL_TOL ${rel_tol})

add_test_compareECLFiles(CASENAME spe5
                         FILENAME SPE5CASE1
                         SIMULATOR flow_solvent
                         ABS_TOL ${abs_tol}
                         REL_TOL ${coarse_rel_tol}
                         TEST_ARGS max_iter=13)

# Restart tests
opm_set_test_driver(${PROJECT_SOURCE_DIR}/tests/run-restart-regressionTest.sh "")

# Cruder tolerances for the restarted tests
set(abs_tol_restart 2e-1)
set(rel_tol_restart 4e-5)
foreach(sim flow flow_legacy flow_ebos)
  add_test_compare_restarted_simulation(CASENAME spe1
                                        FILENAME SPE1CASE2_ACTNUM
                                        SIMULATOR ${sim}
                                        ABS_TOL ${abs_tol_restart}
                                        REL_TOL ${rel_tol_restart})
  add_test_compare_restarted_simulation(CASENAME spe9
                                        FILENAME SPE9_CP_SHORT
                                        SIMULATOR ${sim}
                                        ABS_TOL ${abs_tol_restart}
                                        REL_TOL ${rel_tol_restart})
endforeach()

# Init tests
opm_set_test_driver(${PROJECT_SOURCE_DIR}/tests/run-init-regressionTest.sh "")

foreach(sim flow flow_legacy flow_ebos)
  add_test_compareECLFiles(CASENAME norne
                           FILENAME NORNE_ATW2013
                           SIMULATOR ${sim}
                           ABS_TOL ${abs_tol}
                           REL_TOL ${rel_tol}
                           PREFIX compareECLInitFiles
                           DIR_PREFIX /init)
endforeach()

# Parallel tests
if(MPI_FOUND)
  opm_set_test_driver(${PROJECT_SOURCE_DIR}/tests/run-parallel-regressionTest.sh "")

  # Different tolerances for these tests
  set(abs_tol_parallel 0.02)
  set(rel_tol_parallel 1e-5)

  foreach(sim flow flow_mpi flow_ebos)
    add_test_compare_parallel_simulation(CASENAME spe1
                                         FILENAME SPE1CASE2
                                         SIMULATOR ${sim}
                                         ABS_TOL ${abs_tol_parallel}
                                         REL_TOL ${rel_tol_parallel})
    add_test_compare_parallel_simulation(CASENAME spe9
                                         FILENAME SPE9_CP_SHORT
                                         SIMULATOR ${sim}
                                         ABS_TOL ${abs_tol_parallel}
                                         REL_TOL ${rel_tol_parallel})
  endforeach()
  add_test_compare_parallel_simulation(CASENAME spe3
                                       FILENAME SPE3CASE1
                                       SIMULATOR flow_mpi
                                       ABS_TOL ${abs_tol_parallel}
                                       REL_TOL ${rel_tol_parallel})
endif()
