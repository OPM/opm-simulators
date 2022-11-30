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
  set(oneValueArgs CASENAME FILENAME SIMULATOR DIR DIR_PREFIX PROCS)
  set(multiValueArgs TEST_ARGS)
  cmake_parse_arguments(PARAM "$" "${oneValueArgs}" "${multiValueArgs}" ${ARGN} )
  if(NOT PARAM_DIR)
    set(PARAM_DIR ${PARAM_CASENAME})
  endif()
  set(RESULT_PATH ${BASE_RESULT_PATH}${PARAM_DIR_PREFIX}/${PARAM_SIMULATOR}+${PARAM_CASENAME})
  set(TEST_ARGS ${PARAM_TEST_ARGS})
  opm_add_test(runSimulator/${PARAM_CASENAME} NO_COMPILE
               EXE_NAME ${PARAM_SIMULATOR}
               DRIVER_ARGS -i ${OPM_TESTS_ROOT}/${PARAM_DIR}
                           -r ${RESULT_PATH}
                           -b ${PROJECT_BINARY_DIR}/bin
                           -f ${PARAM_FILENAME}
                           -n ${PARAM_PROCS}
               TEST_ARGS ${TEST_ARGS}
               CONFIGURATION extra)
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
  set(oneValueArgs CASENAME FILENAME SIMULATOR ABS_TOL REL_TOL DIR)
  set(multiValueArgs TEST_ARGS)
  cmake_parse_arguments(PARAM "$" "${oneValueArgs}" "${multiValueArgs}" ${ARGN} )

  if(NOT PARAM_DIR)
    set(PARAM_DIR ${PARAM_CASENAME})
  endif()

  set(RESULT_PATH ${BASE_RESULT_PATH}/parallel/${PARAM_SIMULATOR}+${PARAM_CASENAME})
  set(TEST_ARGS ${OPM_TESTS_ROOT}/${PARAM_DIR}/${PARAM_FILENAME} ${PARAM_TEST_ARGS})

  # Add test that runs flow_mpi and outputs the results to file
  opm_add_test(compareParallelSim_${PARAM_SIMULATOR}+${PARAM_FILENAME} NO_COMPILE
               EXE_NAME ${PARAM_SIMULATOR}
               DRIVER_ARGS -i ${OPM_TESTS_ROOT}/${PARAM_DIR}
                           -r ${RESULT_PATH}
                           -b ${PROJECT_BINARY_DIR}/bin
                           -f ${PARAM_FILENAME}
                           -a ${PARAM_ABS_TOL}
                           -t ${PARAM_REL_TOL}
                           -c ${COMPARE_ECL_COMMAND}
               TEST_ARGS ${TEST_ARGS})
  set_tests_properties(compareParallelSim_${PARAM_SIMULATOR}+${PARAM_FILENAME}
                       PROPERTIES RUN_SERIAL 1)
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
  set(oneValueArgs CASENAME FILENAME SIMULATOR ABS_TOL REL_TOL DIR MPI_PROCS RESTART_STEP)
  set(multiValueArgs TEST_ARGS)
  cmake_parse_arguments(PARAM "$" "${oneValueArgs}" "${multiValueArgs}" ${ARGN} )
  if(NOT PARAM_DIR)
    set(PARAM_DIR ${PARAM_CASENAME})
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
                  -d ${RST_DECK_COMMAND})
  if(PARAM_MPI_PROCS)
    list(APPEND DRIVER_ARGS -n ${PARAM_MPI_PROCS})
  endif()

  opm_add_test(compareParallelRestartedSim_${PARAM_SIMULATOR}+${PARAM_FILENAME} NO_COMPILE
               EXE_NAME ${PARAM_SIMULATOR}
               DRIVER_ARGS ${DRIVER_ARGS}
               TEST_ARGS ${PARAM_TEST_ARGS})
  set_tests_properties(compareParallelRestartedSim_${PARAM_SIMULATOR}+${PARAM_FILENAME}
                       PROPERTIES RUN_SERIAL 1)
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
  set(oneValueArgs CASENAME FILENAME SIMULATOR ABS_TOL REL_TOL DIR)
  set(multiValueArgs TEST_ARGS)
  cmake_parse_arguments(PARAM "$" "${oneValueArgs}" "${multiValueArgs}" ${ARGN} )
  if(NOT PARAM_DIR)
    set(PARAM_DIR ${PARAM_CASENAME})
  endif()
  set(RESULT_PATH ${BASE_RESULT_PATH}/parallelSplitComm/${PARAM_SIMULATOR}+${PARAM_CASENAME})
  set(DRIVER_ARGS -i ${OPM_TESTS_ROOT}/${PARAM_DIR}
                  -r ${RESULT_PATH}
                  -b ${PROJECT_BINARY_DIR}/bin
                  -f ${PARAM_FILENAME}
                  -a ${PARAM_ABS_TOL}
                  -t ${PARAM_REL_TOL}
                  -c ${COMPARE_ECL_COMMAND})
  if(PARAM_MPI_PROCS)
    list(APPEND DRIVER_ARGS -n ${PARAM_MPI_PROCS})
  endif()

  opm_add_test(compareParallelSplitComm_${PARAM_SIMULATOR}+${PARAM_FILENAME} NO_COMPILE
               EXE_NAME ${PARAM_SIMULATOR}
               DRIVER_ARGS ${DRIVER_ARGS}
               TEST_ARGS ${PARAM_TEST_ARGS})
  set_tests_properties(compareParallelSplitComm_${PARAM_SIMULATOR}+${PARAM_FILENAME}
                       PROPERTIES RUN_SERIAL 1)
endfunction()


###########################################################################

### Helper functions ###

macro(parse_json_common_params)
  string(JSON test_cases GET ${JSON_STRING} "test_cases")
  string(JSON simulator GET ${JSON_STRING} "simulator")
  string(JSON abs_tol GET ${JSON_STRING} "abs_tol")
  string(JSON rel_tol GET ${JSON_STRING} "rel_tol")
  string(JSON dir ERROR_VARIABLE err GET ${JSON_STRING} dir)
  if (err)
    # Common case - use directory hosting json file
    get_filename_component(dir ${json_file} DIRECTORY)
    get_filename_component(dir ${dir} NAME)
  endif()
endmacro()

macro(parse_json_test_definition idx)
  string(JSON casename GET ${test_cases} ${idx} casename)
  string(JSON filename GET ${test_cases} ${idx} filename)
  string(JSON loc_abs_tol ERROR_VARIABLE err GET ${test_cases} ${idx} abs_tol)
  if (err)
    set(loc_abs_tol ${abs_tol})
  endif()
  string(JSON loc_rel_tol ERROR_VARIABLE err GET ${test_cases} ${idx} rel_tol)
  if (err)
    set(loc_rel_tol ${rel_tol})
  endif()
  set(PARAMS CASENAME ${casename}
             SIMULATOR ${simulator}
             FILENAME ${filename}
             ABS_TOL ${loc_abs_tol}
             REL_TOL ${loc_rel_tol})
  list(APPEND PARAMS DIR ${dir})
  string(JSON restart_sched ERROR_VARIABLE err GET ${test_cases} ${idx} restart_sched)
  if (NOT err)
    list(APPEND PARAMS RESTART_SCHED ${restart_sched})
  endif()
  string(JSON restart_step ERROR_VARIABLE err GET ${test_cases} ${idx} restart_step)
  if (NOT err)
    list(APPEND PARAMS RESTART_STEP ${restart_step})
  endif()
  string(JSON test_args ERROR_VARIABLE err GET ${test_cases} ${idx} test_args)
  if (NOT err)
    list(APPEND PARAMS TEST_ARGS ${test_args})
  endif()
endmacro()

if(NOT TARGET test-suite)
  add_custom_target(test-suite)
endif()

# Simple execution tests
opm_set_test_driver(${PROJECT_SOURCE_DIR}/tests/run-test.sh "")
add_test_runSimulator(CASENAME norne
                      FILENAME NORNE_ATW2013
                      SIMULATOR flow
                      PROCS 1)

add_test_runSimulator(CASENAME norne_parallel
                      FILENAME NORNE_ATW2013
                      SIMULATOR flow
                      DIR norne
                      PROCS 4)

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

  include (${CMAKE_CURRENT_SOURCE_DIR}/parallelTests.cmake)
endif()


if(OPM_TESTS_ROOT)
  add_custom_target(update_data
                    COMMAND ${CMAKE_COMMAND} --build ${PROJECT_BINARY_DIR} --target check --parallel || exit 0
                    COMMAND ${CMAKE_COMMAND} -E env "REASON=Local\ data\ update" ${PROJECT_SOURCE_DIR}/tests/update_reference_data.sh ${OPM_TESTS_ROOT} ${PROJECT_BINARY_DIR}
                    COMMENT "Updating reference data")
endif()
