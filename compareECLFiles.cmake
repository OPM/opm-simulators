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
# TEST: runSim
###########################################################################

# Input:
#   - casename: basename (no extension)
#
# Details:
#   - This test class simply runs a simulation.
function(add_test_runSimulator)
  set(oneValueArgs CASENAME FILENAME SIMULATOR DIR DIR_PREFIX PROCS CONFIGURATION POST_COMMAND)
  set(multiValueArgs TEST_ARGS)
  cmake_parse_arguments(PARAM "$" "${oneValueArgs}" "${multiValueArgs}" ${ARGN} )
  if(NOT PARAM_DIR)
    set(PARAM_DIR ${PARAM_CASENAME})
  endif()
  set(RESULT_PATH ${BASE_RESULT_PATH}${PARAM_DIR_PREFIX}/${PARAM_SIMULATOR}+${PARAM_CASENAME})
  set(TEST_ARGS ${PARAM_TEST_ARGS})
  set(DRIVER_ARGS -i ${OPM_TESTS_ROOT}/${PARAM_DIR}
                  -r ${RESULT_PATH}
                  -b ${PROJECT_BINARY_DIR}/bin
                  -f ${PARAM_FILENAME})
 if(PARAM_PROCS)
   list(APPEND DRIVER_ARGS -n ${PARAM_PROCS})
 endif()
 if(PARAM_POST_COMMAND)
   list(APPEND DRIVER_ARGS -p "${PARAM_POST_COMMAND}")
 endif()
  opm_add_test(runSimulator/${PARAM_CASENAME} NO_COMPILE
               EXE_NAME ${PARAM_SIMULATOR}
               DRIVER_ARGS ${DRIVER_ARGS}
               TEST_ARGS ${TEST_ARGS}
               CONFIGURATION ${PARAM_CONFIGURATION})
 if(PARAM_PROCS)
    set_tests_properties(runSimulator/${PARAM_CASENAME} PROPERTIES PROCESSORS ${PARAM_PROCS})
  endif()
endfunction()

###########################################################################
# TEST: compareECLFiles
###########################################################################

# Input:
#   - casename: basename (no extension)
#
# Details:
#   - This test class compares output from a simulation to reference files.
function(add_test_compareECLFiles)
  set(oneValueArgs CASENAME FILENAME SIMULATOR ABS_TOL REL_TOL DIR DIR_PREFIX PREFIX RESTART_STEP RESTART_SCHED)
  set(multiValueArgs TEST_ARGS)
  cmake_parse_arguments(PARAM "$" "${oneValueArgs}" "${multiValueArgs}" ${ARGN} )
  if(NOT PARAM_DIR)
    set(PARAM_DIR ${PARAM_CASENAME})
  endif()
  if(NOT PARAM_PREFIX)
    set(PARAM_PREFIX compareECLFiles)
  endif()
  set(RESULT_PATH ${BASE_RESULT_PATH}${PARAM_DIR_PREFIX}/${PARAM_SIMULATOR}+${PARAM_CASENAME})
  set(TEST_ARGS ${PARAM_TEST_ARGS})
  set(DRIVER_ARGS -i ${OPM_TESTS_ROOT}/${PARAM_DIR}
                  -r ${RESULT_PATH}
                  -b ${PROJECT_BINARY_DIR}/bin
                  -f ${PARAM_FILENAME}
                  -a ${PARAM_ABS_TOL}
                  -t ${PARAM_REL_TOL}
                  -c ${COMPARE_ECL_COMMAND}
                  -d ${RST_DECK_COMMAND})
   if(PARAM_RESTART_STEP)
     list(APPEND DRIVER_ARGS -s ${PARAM_RESTART_STEP})
   endif()
  if(PARAM_RESTART_SCHED)
   list(APPEND DRIVER_ARGS -h ${PARAM_RESTART_SCHED})
  endif()
  opm_add_test(${PARAM_PREFIX}_${PARAM_SIMULATOR}+${PARAM_FILENAME} NO_COMPILE
               EXE_NAME ${PARAM_SIMULATOR}
               DRIVER_ARGS ${DRIVER_ARGS}
               TEST_ARGS ${TEST_ARGS})
  set_tests_properties(${PARAM_PREFIX}_${PARAM_SIMULATOR}+${PARAM_FILENAME} PROPERTIES
                        DIRNAME ${PARAM_DIR}
                        FILENAME ${PARAM_FILENAME}
                        SIMULATOR ${PARAM_SIMULATOR}
                        TESTNAME ${PARAM_CASENAME})
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
  set(oneValueArgs CASENAME FILENAME SIMULATOR TEST_NAME ABS_TOL REL_TOL DIR RESTART_STEP)
  set(multiValueArgs TEST_ARGS)
  cmake_parse_arguments(PARAM "$" "${oneValueArgs}" "${multiValueArgs}" ${ARGN} )
  if(NOT PARAM_DIR)
    set(PARAM_DIR ${PARAM_CASENAME})
  endif()
  if (PARAM_TEST_NAME)
    set(TEST_NAME ${PARAM_TEST_NAME})
  else()
    set(TEST_NAME compareRestartedSim_${PARAM_SIMULATOR}+${PARAM_FILENAME})
  endif()

  set(RESULT_PATH ${BASE_RESULT_PATH}/restart/${PARAM_SIMULATOR}+${PARAM_CASENAME})
  opm_add_test(${TEST_NAME} NO_COMPILE
               EXE_NAME ${PARAM_SIMULATOR}
               DRIVER_ARGS -i ${OPM_TESTS_ROOT}/${PARAM_DIR}
                           -r ${RESULT_PATH}
                           -b ${PROJECT_BINARY_DIR}/bin
                           -f ${PARAM_FILENAME}
                           -a ${PARAM_ABS_TOL}
                           -t ${PARAM_REL_TOL}
                           -c ${COMPARE_ECL_COMMAND}
                           -d ${RST_DECK_COMMAND}
                           -s ${PARAM_RESTART_STEP}
               TEST_ARGS ${PARAM_TEST_ARGS})
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
  set(oneValueArgs CASENAME FILENAME SIMULATOR ABS_TOL REL_TOL DIR MPI_PROCS)
  set(multiValueArgs TEST_ARGS)
  cmake_parse_arguments(PARAM "$" "${oneValueArgs}" "${multiValueArgs}" ${ARGN} )

  if(NOT PARAM_DIR)
    set(PARAM_DIR ${PARAM_CASENAME})
  endif()

  if(PARAM_MPI_PROCS)
    set(MPI_PROCS ${PARAM_MPI_PROCS})
  else()
    set(MPI_PROCS 4)
  endif()

  set(RESULT_PATH ${BASE_RESULT_PATH}/parallel/${PARAM_SIMULATOR}+${PARAM_CASENAME})
  set(TEST_ARGS ${OPM_TESTS_ROOT}/${PARAM_DIR}/${PARAM_FILENAME} ${PARAM_TEST_ARGS})
  set(DRIVER_ARGS -i ${OPM_TESTS_ROOT}/${PARAM_DIR}
                  -r ${RESULT_PATH}
                  -b ${PROJECT_BINARY_DIR}/bin
                  -f ${PARAM_FILENAME}
                  -a ${PARAM_ABS_TOL}
                  -t ${PARAM_REL_TOL}
                  -c ${COMPARE_ECL_COMMAND}
                  -n ${MPI_PROCS})

  # Add test that runs flow_mpi and outputs the results to file
  opm_add_test(compareParallelSim_${PARAM_SIMULATOR}+${PARAM_FILENAME} NO_COMPILE
               EXE_NAME ${PARAM_SIMULATOR}
               DRIVER_ARGS ${DRIVER_ARGS}
               TEST_ARGS ${TEST_ARGS})
  set_tests_properties(compareParallelSim_${PARAM_SIMULATOR}+${PARAM_FILENAME}
                       PROPERTIES PROCESSORS ${MPI_PROCS})
endfunction()


###########################################################################
# TEST: add_test_compare_parallel_restarted_simulation
###########################################################################

# Input:
#   - casename: basename (no extension)
#
# Details:
#   - This test class compares the output from a restarted parallel simulation
#     to that of a non-restarted parallel simulation.
function(add_test_compare_parallel_restarted_simulation)
  set(oneValueArgs CASENAME FILENAME SIMULATOR ABS_TOL REL_TOL DIR MPI_PROCS RESTART_STEP TEST_NAME)
  set(multiValueArgs TEST_ARGS)
  cmake_parse_arguments(PARAM "$" "${oneValueArgs}" "${multiValueArgs}" ${ARGN} )

  if(NOT PARAM_DIR)
    set(PARAM_DIR ${PARAM_CASENAME})
  endif()

  if (PARAM_TEST_NAME)
    set(TEST_NAME ${PARAM_TEST_NAME})
  else()
    set(TEST_NAME compareParallelRestartedSim_${PARAM_SIMULATOR}+${PARAM_FILENAME})
  endif()

  if(PARAM_MPI_PROCS)
    set(MPI_PROCS ${PARAM_MPI_PROCS})
  else()
    set(MPI_PROCS 4)
  endif()

  set(RESULT_PATH ${BASE_RESULT_PATH}/parallelRestart/${PARAM_SIMULATOR}+${PARAM_CASENAME})
  set(DRIVER_ARGS -i ${OPM_TESTS_ROOT}/${PARAM_DIR}
                  -r ${RESULT_PATH}
                  -b ${PROJECT_BINARY_DIR}/bin
                  -f ${PARAM_FILENAME}
                  -a ${PARAM_ABS_TOL}
                  -t ${PARAM_REL_TOL}
                  -c ${COMPARE_ECL_COMMAND}
                  -s ${PARAM_RESTART_STEP}
                  -d ${RST_DECK_COMMAND}
                  -n ${MPI_PROCS})

  opm_add_test(${TEST_NAME} NO_COMPILE
               EXE_NAME ${PARAM_SIMULATOR}
               DRIVER_ARGS ${DRIVER_ARGS}
               TEST_ARGS ${PARAM_TEST_ARGS})
  set_tests_properties(${TEST_NAME} PROPERTIES PROCESSORS ${MPI_PROCS})
endfunction()


###########################################################################
# TEST: add_test_split_comm
###########################################################################

# Input:
#   - casename: basename (no extension)
#
# Details:
#   - This test class compares the output from a parallel simulation
#     to that of a parallel simulation running with a custom communicator.
function(add_test_split_comm)
  set(oneValueArgs CASENAME FILENAME SIMULATOR ABS_TOL REL_TOL DIR MPI_PROCS)
  set(multiValueArgs TEST_ARGS)
  cmake_parse_arguments(PARAM "$" "${oneValueArgs}" "${multiValueArgs}" ${ARGN} )

  if(NOT PARAM_DIR)
    set(PARAM_DIR ${PARAM_CASENAME})
  endif()

  if(PARAM_MPI_PROCS)
    set(MPI_PROCS ${PARAM_MPI_PROCS})
  else()
    set(MPI_PROCS 4)
  endif()

  set(RESULT_PATH ${BASE_RESULT_PATH}/parallelSplitComm/${PARAM_SIMULATOR}+${PARAM_CASENAME})
  set(DRIVER_ARGS -i ${OPM_TESTS_ROOT}/${PARAM_DIR}
                  -r ${RESULT_PATH}
                  -b ${PROJECT_BINARY_DIR}/bin
                  -f ${PARAM_FILENAME}
                  -a ${PARAM_ABS_TOL}
                  -t ${PARAM_REL_TOL}
                  -c ${COMPARE_ECL_COMMAND}
                  -n ${MPI_PROCS})

  opm_add_test(compareParallelSplitComm_${PARAM_SIMULATOR}+${PARAM_FILENAME} NO_COMPILE
               EXE_NAME ${PARAM_SIMULATOR}
               DRIVER_ARGS ${DRIVER_ARGS}
               TEST_ARGS ${PARAM_TEST_ARGS})
  set_tests_properties(compareParallelSplitComm_${PARAM_SIMULATOR}+${PARAM_FILENAME}
                       PROPERTIES PROCESSORS ${MPI_PROCS})
endfunction()


###########################################################################

###########################################################################
# TEST: compareDamarisFiles
###########################################################################

# Input:
#   - casename: basename (no extension)
#
# Details:
#   - This test class compares output from a simulation to reference files.
function(add_test_compareDamarisFiles)
  set(oneValueArgs CASENAME FILENAME SIMULATOR ABS_TOL REL_TOL DIR DIR_PREFIX PREFIX MPI_PROCS)
  set(multiValueArgs TEST_ARGS)
  cmake_parse_arguments(PARAM "$" "${oneValueArgs}" "${multiValueArgs}" ${ARGN} )
  if(NOT PARAM_DIR)
    set(PARAM_DIR ${PARAM_CASENAME})
  endif()
  if(NOT PARAM_PREFIX)
    set(PARAM_PREFIX compareDamarisFiles)
  endif()
  if(PARAM_MPI_PROCS)
    set(MPI_PROCS ${PARAM_MPI_PROCS})
  else()
    set(MPI_PROCS 4)
  endif()
  set(RESULT_PATH ${BASE_RESULT_PATH}${PARAM_DIR_PREFIX}/${PARAM_SIMULATOR}+${PARAM_CASENAME})
  set(TEST_ARGS ${PARAM_TEST_ARGS})
  set(DRIVER_ARGS -i ${OPM_TESTS_ROOT}/${PARAM_DIR}
                  -r ${RESULT_PATH}
                  -b ${PROJECT_BINARY_DIR}/bin
                  -f ${PARAM_FILENAME}
                  -a ${PARAM_ABS_TOL}
                  -t ${PARAM_REL_TOL}
                  -c ${HDF5_DIFF_EXECUTABLE}
                  -n ${MPI_PROCS})
  set(TEST_NAME ${PARAM_PREFIX}_${PARAM_SIMULATOR}+${PARAM_FILENAME})
  opm_add_test(${TEST_NAME} NO_COMPILE
               EXE_NAME ${PARAM_SIMULATOR}
               DRIVER_ARGS ${DRIVER_ARGS}
               TEST_ARGS ${TEST_ARGS})
  set_tests_properties(${TEST_NAME} PROPERTIES PROCESSORS ${MPI_PROCS})
endfunction()


if(NOT TARGET test-suite)
  add_custom_target(test-suite)
endif()

# Simple execution tests
opm_set_test_driver(${PROJECT_SOURCE_DIR}/tests/run-test.sh "")
add_test_runSimulator(CASENAME norne
                      FILENAME NORNE_ATW2013
                      SIMULATOR flow
                      CONFIGURATION extra)

add_test_runSimulator(CASENAME norne_parallel
                      FILENAME NORNE_ATW2013
                      SIMULATOR flow
                      DIR norne
                      PROCS 4
                      CONFIGURATION extra)

# Tests that are run based on simulator results, but not necessarily direct comparison to reference results
add_test_runSimulator(CASENAME tuning_xxxmbe
                      FILENAME 01_TUNING_XXXMBE
                      SIMULATOR flow
											DIR tuning
                      TEST_ARGS --output-extra-convergence-info=iterations --enable-tuning=true
                      POST_COMMAND $<TARGET_FILE:test_tuning_xxxMBE>)

add_test_runSimulator(CASENAME notuning_xxxmbe
                      FILENAME 01_TUNING_XXXMBE
                      SIMULATOR flow
											DIR tuning
                      TEST_ARGS --output-extra-convergence-info=iterations --enable-tuning=false
                      POST_COMMAND $<TARGET_FILE:test_tuning_xxxMBE>)

set_tests_properties(runSimulator/notuning_xxxmbe PROPERTIES WILL_FAIL TRUE)

add_test_runSimulator(CASENAME tuning_tsinit_nextstep
                      FILENAME 02_TUNING_TSINIT_NEXTSTEP
                      SIMULATOR flow
											DIR tuning
                      TEST_ARGS --enable-tuning=true
                      POST_COMMAND $<TARGET_FILE:test_tuning_tsinit_nextstep>)


include (${CMAKE_CURRENT_SOURCE_DIR}/regressionTests.cmake)
include (${CMAKE_CURRENT_SOURCE_DIR}/restartTests.cmake)

# PORV test
opm_set_test_driver(${PROJECT_SOURCE_DIR}/tests/run-porv-acceptanceTest.sh "")
add_test_compareECLFiles(CASENAME norne
                         FILENAME NORNE_ATW2013
                         SIMULATOR flow
                         ABS_TOL 1e-5
                         REL_TOL 1e-8
                         PREFIX comparePORV
                         DIR_PREFIX /porv)

# Init tests
opm_set_test_driver(${PROJECT_SOURCE_DIR}/tests/run-init-regressionTest.sh "")

add_test_compareECLFiles(CASENAME norne
                         FILENAME NORNE_ATW2013
                         SIMULATOR flow
                         ABS_TOL ${abs_tol}
                         REL_TOL ${rel_tol}
                         PREFIX compareECLInitFiles
                         DIR_PREFIX /init)

# This is not a proper regression test; the test will load a norne case prepared
# for restart and run one single timestep - of length one day. The results are not
# verified in any way.
add_test(NAME NORNE_RESTART
         COMMAND flow --output-dir=${BASE_RESULT_PATH}/norne-restart ${OPM_TESTS_ROOT}/norne/NORNE_ATW2013_RESTART.DATA)


if(MPI_FOUND)
  include (${CMAKE_CURRENT_SOURCE_DIR}/parallelRestartTests.cmake)

  # Single test to verify that we treat custom communicators correctly.
  opm_set_test_driver(${PROJECT_SOURCE_DIR}/tests/run-split-comm-test.sh "")
  add_test_split_comm(CASENAME spe1
                      FILENAME SPE1CASE2
                      SIMULATOR flow
                      ABS_TOL 0.0
                      REL_TOL 0.0)

  # Single test for damaris
  if(Damaris_FOUND AND USE_DAMARIS_LIB)
      opm_set_test_driver(${PROJECT_SOURCE_DIR}/tests/run-damaris-regressionTest.sh "")
      add_test_compareDamarisFiles(CASENAME spe1_damaris
                      FILENAME SPE1CASE1
                      SIMULATOR flow
                      ABS_TOL ${abs_tol}
                      REL_TOL ${rel_tol}
                      DIR spe1)
    endif()

  include (${CMAKE_CURRENT_SOURCE_DIR}/parallelTests.cmake)
endif()


if(OPM_TESTS_ROOT)
  add_custom_target(update_data
                    COMMAND ${CMAKE_COMMAND} --build ${PROJECT_BINARY_DIR} --target check --parallel || exit 0
                    COMMAND ${CMAKE_COMMAND} -E env "REASON=Local\ data\ update" ${PROJECT_SOURCE_DIR}/tests/update_reference_data.sh ${OPM_TESTS_ROOT} ${PROJECT_BINARY_DIR}
                    COMMENT "Updating reference data")
endif()
