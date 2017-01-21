# Set absolute tolerance to be used for testing
set(abs_tol 2e-2)
set(rel_tol 1e-5)

# Define some paths
set(BASE_RESULT_PATH ${PROJECT_BINARY_DIR}/tests/results)

###########################################################################
# TEST: compareECLFiles
###########################################################################

# Input:
#   - casename: basename (no extension)
#
macro (add_test_compareECLFiles casename filename simulator)

  set(RESULT_PATH ${BASE_RESULT_PATH}/${simulator}+${casename})
  opm_add_test(compareECLFiles_${simulator}+${filename} NO_COMPILE
               EXE_NAME ${simulator}
               DRIVER_ARGS ${OPM_DATA_ROOT}/${casename} ${RESULT_PATH}
                           ${CMAKE_BINARY_DIR}/bin
                           ${filename}
                           ${abs_tol} ${rel_tol}
                           ${COMPARE_SUMMARY_COMMAND}
                           ${COMPARE_ECL_COMMAND}
               TEST_ARGS ${OPM_DATA_ROOT}/${casename}/${filename}.DATA )
endmacro (add_test_compareECLFiles)

###########################################################################
# TEST: compareECLRestartFiles
###########################################################################

# Input:
#   - casename: basename (no extension)
#
macro (add_test_compareECLRestartFiles casename filename simulator)

  set(RESULT_PATH ${BASE_RESULT_PATH}/restart/${simulator}+${casename})
  opm_add_test(compareECLRestartFiles_${simulator}+${filename} NO_COMPILE
               EXE_NAME ${simulator}
               DRIVER_ARGS ${OPM_DATA_ROOT}/${casename} ${RESULT_PATH}
                           ${CMAKE_BINARY_DIR}/bin
                           ${filename}
                           ${abs_tol} ${rel_tol}
                           ${COMPARE_SUMMARY_COMMAND}
                           ${COMPARE_ECL_COMMAND}
               TEST_ARGS ${OPM_DATA_ROOT}/${casename}/${filename})
endmacro (add_test_compareECLRestartFiles)

###########################################################################
# TEST: compareECLInitFiles
###########################################################################

# Input:
#   - casename: basename (no extension)
#
macro (add_test_compareECLInitFiles casename filename simulator)

  set(RESULT_PATH ${BASE_RESULT_PATH}/init/${simulator}+${casename})
  opm_add_test(compareECLInitFiles_${simulator}+${filename} NO_COMPILE
               EXE_NAME ${simulator}
               DRIVER_ARGS ${OPM_DATA_ROOT}/${casename} ${RESULT_PATH}
                           ${CMAKE_BINARY_DIR}/bin
                           ${filename}
                           ${abs_tol} ${rel_tol}
                           ${COMPARE_SUMMARY_COMMAND}
                           ${COMPARE_ECL_COMMAND}
               TEST_ARGS ${OPM_DATA_ROOT}/${casename}/${filename})
endmacro (add_test_compareECLInitFiles)

###########################################################################
# TEST: parallelECLFiles
###########################################################################

# Input:
#   - casename: basename (no extension)
#
macro (add_test_parallelECLFiles casename filename simulator)
  set(abs_tol 0.20)
  set(rel_tol 4e-4)
  set(RESULT_PATH ${BASE_RESULT_PATH}/parallel/${simulator}+${casename})

  # Add test that runs flow_mpi and outputs the results to file
  opm_add_test(parallelECLFiles_${simulator}+${filename} NO_COMPILE
               EXE_NAME ${simulator}
               DRIVER_ARGS ${OPM_DATA_ROOT}/${casename} ${RESULT_PATH}
                           ${CMAKE_BINARY_DIR}/bin
                           ${filename}
                           ${abs_tol} ${rel_tol}
                           ${COMPARE_SUMMARY_COMMAND}
                           ${COMPARE_ECL_COMMAND}
               TEST_ARGS ${OPM_DATA_ROOT}/${casename}/${filename})
endmacro (add_test_parallelECLFiles)

if(NOT TARGET test-suite)
  add_custom_target(test-suite)
endif()

# Regression tests
opm_set_test_driver(${PROJECT_SOURCE_DIR}/tests/run-regressionTest.sh "")

add_test_compareECLFiles(spe1 SPE1CASE2 flow)
add_test_compareECLFiles(spe1 SPE1CASE1 flow_sequential)
add_test_compareECLFiles(spe3 SPE3CASE1 flow)
add_test_compareECLFiles(spe9 SPE9_CP_SHORT flow)

# Restart tests
opm_set_test_driver(${PROJECT_SOURCE_DIR}/tests/run-restart-regressionTest.sh "")

add_test_compareECLRestartFiles(spe1 SPE1CASE2_ACTNUM flow)
add_test_compareECLRestartFiles(spe9 SPE9_CP_SHORT flow)

# Init tests
opm_set_test_driver(${PROJECT_SOURCE_DIR}/tests/run-init-regressionTest.sh "")

add_test_compareECLInitFiles(norne NORNE_ATW2013 flow)

# Parallel tests
if(MPI_FOUND)
  opm_set_test_driver(${PROJECT_SOURCE_DIR}/tests/run-parallel-regressionTest.sh "")

  add_test_parallelECLFiles(spe1 SPE1CASE2 flow_mpi)
  add_test_parallelECLFiles(spe3 SPE3CASE1 flow_mpi)
  add_test_parallelECLFiles(spe9 SPE9_CP_SHORT flow_mpi)
endif()
