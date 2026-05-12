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

# Select which executable should be used for blackoil cases that currently
# register the generic flow simulator.
set(OPM_BLACKOIL_TEST_SIMULATOR "flow" CACHE STRING
    "Simulator executable used for blackoil test cases")
set_property(CACHE OPM_BLACKOIL_TEST_SIMULATOR PROPERTY STRINGS flow flow_blackoil)

set(OPM_BLACKOIL_ALLOWED_RUNSPEC_MODEL_KEYWORDS
  OIL GAS WATER DISGAS DISOIL VAPOIL)
set(OPM_BLACKOIL_DISALLOWED_RUNSPEC_MODEL_KEYWORDS
  COMPS DISGASW VAPWAT TEMP THERMAL BRINE PRECSALT POLYMER FOAM SOLVENT WSOLVENT
  CO2SOL H2SOL CO2STORE H2STORE MICP DIFFUSE MECH)

function(get_blackoil_reference_simulator simulator reference_simulator_var)
  if("${simulator}" STREQUAL "${OPM_BLACKOIL_TEST_SIMULATOR}")
    set(${reference_simulator_var} "flow" PARENT_SCOPE)
  else()
    set(${reference_simulator_var} "${simulator}" PARENT_SCOPE)
  endif()
endfunction()

function(resolve_test_deck_file deck_dir deck_name resolved_deck_var)
  set(deck_root "${OPM_TESTS_ROOT}/${deck_dir}")
  foreach(suffix "" ".DATA" ".data" ".BASE" ".base" ".INC" ".inc")
    set(candidate "${deck_root}/${deck_name}${suffix}")
    if(EXISTS "${candidate}")
      set(${resolved_deck_var} "${candidate}" PARENT_SCOPE)
      return()
    endif()
  endforeach()
  set(${resolved_deck_var} "" PARENT_SCOPE)
endfunction()

function(collect_runspec_blackoil_keywords deck_path inherited_section visited_var keywords_var)
  get_filename_component(abs_deck "${deck_path}" ABSOLUTE)
  set(visited ${${visited_var}})
  list(FIND visited "${abs_deck}" visited_index)
  if(NOT visited_index EQUAL -1)
    set(${visited_var} "${visited}" PARENT_SCOPE)
    set(${keywords_var} "" PARENT_SCOPE)
    return()
  endif()

  list(APPEND visited "${abs_deck}")
  if(NOT EXISTS "${abs_deck}")
    set(${visited_var} "${visited}" PARENT_SCOPE)
    set(${keywords_var} "" PARENT_SCOPE)
    return()
  endif()

  file(READ "${abs_deck}" deck_content)
  string(REGEX REPLACE "--[^\n]*" "" deck_content "${deck_content}")
  string(REGEX MATCHALL "'[^']*'|\"[^\"]*\"|/|[^ \t\r\n/]+" tokens "${deck_content}")

  set(current_section "${inherited_section}")
  set(local_keywords "")
  list(LENGTH tokens token_count)
  set(index 0)

  while(index LESS token_count)
    list(GET tokens ${index} token)
    if(token STREQUAL "/")
      math(EXPR index "${index} + 1")
      continue()
    endif()

    string(REGEX REPLACE "^[\"']|[\"']$" "" bare_token "${token}")
    string(TOUPPER "${bare_token}" upper_token)

    if(upper_token STREQUAL "RUNSPEC" OR
       upper_token STREQUAL "GRID" OR
       upper_token STREQUAL "EDIT" OR
       upper_token STREQUAL "PROPS" OR
       upper_token STREQUAL "REGIONS" OR
       upper_token STREQUAL "SOLUTION" OR
       upper_token STREQUAL "SUMMARY" OR
       upper_token STREQUAL "SCHEDULE")
      set(current_section "${upper_token}")
      math(EXPR index "${index} + 1")
      continue()
    endif()

    if(current_section STREQUAL "RUNSPEC" AND upper_token STREQUAL "INCLUDE")
      math(EXPR index "${index} + 1")
      while(index LESS token_count)
        list(GET tokens ${index} include_token)
        if(include_token STREQUAL "/")
          break()
        endif()

        string(REGEX REPLACE "^[\"']|[\"']$" "" include_name "${include_token}")
        if(NOT include_name MATCHES "[$<>]")
          get_filename_component(deck_dirname "${abs_deck}" DIRECTORY)
          if(IS_ABSOLUTE "${include_name}")
            set(include_path "${include_name}")
          else()
            set(include_path "${deck_dirname}/${include_name}")
          endif()
          if(EXISTS "${include_path}")
            collect_runspec_blackoil_keywords("${include_path}" "${current_section}" visited include_keywords)
            list(APPEND local_keywords ${include_keywords})
          endif()
        endif()

        math(EXPR index "${index} + 1")
      endwhile()
      math(EXPR index "${index} + 1")
      continue()
    endif()

    if(current_section STREQUAL "RUNSPEC")
      list(FIND OPM_BLACKOIL_ALLOWED_RUNSPEC_MODEL_KEYWORDS "${upper_token}" allowed_index)
      list(FIND OPM_BLACKOIL_DISALLOWED_RUNSPEC_MODEL_KEYWORDS "${upper_token}" disallowed_index)
      if(NOT "${allowed_index}" STREQUAL "-1" OR NOT "${disallowed_index}" STREQUAL "-1")
        list(APPEND local_keywords "${upper_token}")
      endif()
    endif()

    math(EXPR index "${index} + 1")
  endwhile()

  list(REMOVE_DUPLICATES local_keywords)
  set(${visited_var} "${visited}" PARENT_SCOPE)
  set(${keywords_var} "${local_keywords}" PARENT_SCOPE)
endfunction()

function(test_supports_blackoil_simulator deck_dir deck_name supports_var)
  resolve_test_deck_file("${deck_dir}" "${deck_name}" resolved_deck)
  if(NOT resolved_deck)
    set(${supports_var} FALSE PARENT_SCOPE)
    return()
  endif()

  set(visited "")
  collect_runspec_blackoil_keywords("${resolved_deck}" "" visited keywords)
  foreach(phase IN ITEMS OIL GAS WATER)
    list(FIND keywords "${phase}" phase_index)
    if(phase_index EQUAL -1)
      set(${supports_var} FALSE PARENT_SCOPE)
      return()
    endif()
  endforeach()

  foreach(disallowed_keyword IN LISTS OPM_BLACKOIL_DISALLOWED_RUNSPEC_MODEL_KEYWORDS)
    list(FIND keywords "${disallowed_keyword}" disallowed_index)
    if(NOT disallowed_index EQUAL -1)
      set(${supports_var} FALSE PARENT_SCOPE)
      return()
    endif()
  endforeach()

  set(${supports_var} TRUE PARENT_SCOPE)
endfunction()

function(resolve_blackoil_test_simulator simulator is_blackoil_data resolved_simulator_var)
  if(("${simulator}" STREQUAL "flow" OR "${simulator}" STREQUAL "flow_blackoil") AND NOT ${is_blackoil_data})
    set(${resolved_simulator_var} "flow" PARENT_SCOPE)
  else()
    set(${resolved_simulator_var} "${simulator}" PARENT_SCOPE)
  endif()
endfunction()

function(resolve_unlabeled_test_simulator simulator resolved_simulator_var)
  if("${simulator}" STREQUAL "flow" OR "${simulator}" STREQUAL "flow_blackoil")
    set(${resolved_simulator_var} "flow" PARENT_SCOPE)
  else()
    set(${resolved_simulator_var} "${simulator}" PARENT_SCOPE)
  endif()
endfunction()

function(label_blackoil_test test_name simulator is_blackoil)
  if(${is_blackoil})
    if("${simulator}" STREQUAL "flow" OR "${simulator}" STREQUAL "flow_blackoil")
      set_property(TEST ${test_name} APPEND PROPERTY LABELS ${ARGN})
    endif()
  endif()
endfunction()

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
  resolve_unlabeled_test_simulator("${PARAM_SIMULATOR}" TEST_SIMULATOR)
  set(RESULT_PATH ${BASE_RESULT_PATH}${PARAM_DIR_PREFIX}/${TEST_SIMULATOR}+${PARAM_CASENAME})
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
            EXE_NAME ${TEST_SIMULATOR}
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

  test_supports_blackoil_simulator("${PARAM_DIR}" "${PARAM_FILENAME}" IS_BLACKOIL_DATA)
  resolve_blackoil_test_simulator("${PARAM_SIMULATOR}" ${IS_BLACKOIL_DATA} TEST_SIMULATOR)
  get_blackoil_reference_simulator("${TEST_SIMULATOR}" REFERENCE_SIMULATOR)

  set(RESULT_PATH ${BASE_RESULT_PATH}${PARAM_DIR_PREFIX}/${TEST_SIMULATOR}+${PARAM_CASENAME})
  set(TEST_ARGS ${PARAM_TEST_ARGS})
  set(DRIVER_ARGS -i ${OPM_TESTS_ROOT}/${PARAM_DIR}
                  -r ${RESULT_PATH}
                  -b ${PROJECT_BINARY_DIR}/bin
                  -f ${PARAM_FILENAME}
                  -a ${PARAM_ABS_TOL}
                  -t ${PARAM_REL_TOL}
                  -c $<TARGET_FILE:compareECL>
                  -d $<TARGET_FILE:rst_deck>
                  -u ${REFERENCE_SIMULATOR})
  if(PARAM_RESTART_STEP)
    list(APPEND DRIVER_ARGS -s ${PARAM_RESTART_STEP})
  endif()
  if(PARAM_RESTART_SCHED STREQUAL "false" OR PARAM_RESTART_SCHED STREQUAL "true")
    list(APPEND DRIVER_ARGS -h ${PARAM_RESTART_SCHED})
  endif()
  opm_add_test(${PARAM_PREFIX}_${TEST_SIMULATOR}+${PARAM_FILENAME} NO_COMPILE
               EXE_NAME ${TEST_SIMULATOR}
               DRIVER_ARGS ${DRIVER_ARGS}
               TEST_ARGS ${TEST_ARGS})
  set_tests_properties(${PARAM_PREFIX}_${TEST_SIMULATOR}+${PARAM_FILENAME} PROPERTIES
                        DIRNAME ${PARAM_DIR}
                        FILENAME ${PARAM_FILENAME}
                        SIMULATOR ${TEST_SIMULATOR}
                        TESTNAME ${PARAM_CASENAME})

  if(PARAM_PREFIX STREQUAL "compareECLFiles")
    label_blackoil_test(${PARAM_PREFIX}_${TEST_SIMULATOR}+${PARAM_FILENAME} "${TEST_SIMULATOR}" ${IS_BLACKOIL_DATA} blackoil blackoil-regression)
  elseif(PARAM_PREFIX STREQUAL "compareECLInitFiles" OR PARAM_PREFIX STREQUAL "comparePORV")
    label_blackoil_test(${PARAM_PREFIX}_${TEST_SIMULATOR}+${PARAM_FILENAME} "${TEST_SIMULATOR}" ${IS_BLACKOIL_DATA} blackoil blackoil-reference)
  endif()
endfunction()

###########################################################################
# TEST: compareSeparateECLFiles
###########################################################################

# Input:
#   - casename: basename (no extension)
#   - filename1 (no extension)
#   - filename2 (no extension)
#
# Details:
#   - This test class compares two separate simulations
function(add_test_compareSeparateECLFiles)
  set(oneValueArgs CASENAME FILENAME1 FILENAME2 DIR1 DIR2 SIMULATOR ABS_TOL REL_TOL IGNORE_EXTRA_KW DIR_PREFIX MPI_PROCS)
  set(multiValueArgs TEST_ARGS)
  cmake_parse_arguments(PARAM "$" "${oneValueArgs}" "${multiValueArgs}" ${ARGN} )
  if(NOT PARAM_PREFIX)
    set(PARAM_PREFIX compareSeparateECLFiles)
  endif()
  if(PARAM_MPI_PROCS)
    set(MPI_PROCS ${PARAM_MPI_PROCS})
  else()
    set(MPI_PROCS 1)
  endif()
  if(NOT PARAM_DIR)
    set(PARAM_DIR ${PARAM_CASENAME})
  endif()

  test_supports_blackoil_simulator("${PARAM_DIR1}" "${PARAM_FILENAME1}" IS_BLACKOIL_DATA_1)
  test_supports_blackoil_simulator("${PARAM_DIR2}" "${PARAM_FILENAME2}" IS_BLACKOIL_DATA_2)
  if(IS_BLACKOIL_DATA_1 AND IS_BLACKOIL_DATA_2)
    set(IS_BLACKOIL_DATA TRUE)
  else()
    set(IS_BLACKOIL_DATA FALSE)
  endif()
  resolve_blackoil_test_simulator("${PARAM_SIMULATOR}" ${IS_BLACKOIL_DATA} TEST_SIMULATOR)

  set(RESULT_PATH ${BASE_RESULT_PATH}${PARAM_DIR_PREFIX}/${TEST_SIMULATOR}+${PARAM_CASENAME})
  set(TEST_ARGS ${PARAM_TEST_ARGS})
  set(DRIVER_ARGS -i ${OPM_TESTS_ROOT}/${PARAM_DIR1}
                  -j ${OPM_TESTS_ROOT}/${PARAM_DIR2}
                  -f ${PARAM_FILENAME1}
                  -g ${PARAM_FILENAME2}
                  -r ${RESULT_PATH}
                  -b ${PROJECT_BINARY_DIR}/bin
                  -a ${PARAM_ABS_TOL}
                  -t ${PARAM_REL_TOL}
                  -c $<TARGET_FILE:compareECL>
                  -n ${MPI_PROCS})
  if(PARAM_IGNORE_EXTRA_KW)
    list(APPEND DRIVER_ARGS -y ${PARAM_IGNORE_EXTRA_KW})
  endif()
  opm_add_test(${PARAM_PREFIX}_${TEST_SIMULATOR}+${PARAM_CASENAME} NO_COMPILE
               EXE_NAME ${TEST_SIMULATOR}
               DRIVER_ARGS ${DRIVER_ARGS}
               TEST_ARGS ${TEST_ARGS})
  set_tests_properties(${PARAM_PREFIX}_${TEST_SIMULATOR}+${PARAM_CASENAME} PROPERTIES
                        DIRNAME ${PARAM_DIR}
                        FILENAME1 ${PARAM_FILENAME1}
                        FILENAME2 ${PARAM_FILENAME2}
                        SIMULATOR ${TEST_SIMULATOR}
                        TESTNAME ${PARAM_CASENAME}
                        PROCESSORS ${MPI_PROCS})
  if(CMAKE_VERSION VERSION_GREATER_EQUAL 3.22)
    set_tests_properties(${PARAM_PREFIX}_${TEST_SIMULATOR}+${PARAM_CASENAME} PROPERTIES
                         ENVIRONMENT_MODIFICATION PYTHONPATH=path_list_append:${opm-common_DIR}/python)
  endif()

  if(IS_BLACKOIL_DATA)
    label_blackoil_test(${PARAM_PREFIX}_${TEST_SIMULATOR}+${PARAM_CASENAME} "${TEST_SIMULATOR}" TRUE blackoil blackoil-comparison)
  endif()
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

  test_supports_blackoil_simulator("${PARAM_DIR}" "${PARAM_FILENAME}" IS_BLACKOIL_DATA)
  resolve_blackoil_test_simulator("${PARAM_SIMULATOR}" ${IS_BLACKOIL_DATA} TEST_SIMULATOR)

  if (PARAM_TEST_NAME)
    set(TEST_NAME ${PARAM_TEST_NAME})
  else()
    set(TEST_NAME compareRestartedSim_${TEST_SIMULATOR}+${PARAM_FILENAME})
  endif()

  set(RESULT_PATH ${BASE_RESULT_PATH}/restart/${TEST_SIMULATOR}+${PARAM_CASENAME})
  opm_add_test(${TEST_NAME} NO_COMPILE
               EXE_NAME ${TEST_SIMULATOR}
               DRIVER_ARGS -i ${OPM_TESTS_ROOT}/${PARAM_DIR}
                           -r ${RESULT_PATH}
                           -b ${PROJECT_BINARY_DIR}/bin
                           -f ${PARAM_FILENAME}
                           -a ${PARAM_ABS_TOL}
                           -t ${PARAM_REL_TOL}
                           -c $<TARGET_FILE:compareECL>
                           -d $<TARGET_FILE:rst_deck>
                           -s ${PARAM_RESTART_STEP}
               TEST_ARGS ${PARAM_TEST_ARGS})

  label_blackoil_test(${TEST_NAME} "${TEST_SIMULATOR}" ${IS_BLACKOIL_DATA} blackoil blackoil-restart blackoil-consistency)
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
  set(oneValueArgs CASENAME FILENAME SIMULATOR ABS_TOL REL_TOL DIR POSTFIX MPI_PROCS)
  set(multiValueArgs TEST_ARGS)
  cmake_parse_arguments(PARAM "$" "${oneValueArgs}" "${multiValueArgs}" ${ARGN} )

  if(NOT PARAM_DIR)
    set(PARAM_DIR ${PARAM_CASENAME})
  endif()

  if(NOT PARAM_POSTFIX)
    set(PARAM_POSTFIX "")
  else()
    set(PARAM_POSTFIX +${PARAM_POSTFIX})
  endif()

  if(PARAM_MPI_PROCS)
    set(MPI_PROCS ${PARAM_MPI_PROCS})
  else()
    set(MPI_PROCS 4)
  endif()

  test_supports_blackoil_simulator("${PARAM_DIR}" "${PARAM_FILENAME}" IS_BLACKOIL_DATA)
  resolve_blackoil_test_simulator("${PARAM_SIMULATOR}" ${IS_BLACKOIL_DATA} TEST_SIMULATOR)

  if(MPIEXEC_MAX_NUMPROCS GREATER_EQUAL MPI_PROCS)
    # Local computer system has at least ${MPI_PROCS} CPUs. Register test.
    set(RESULT_PATH ${BASE_RESULT_PATH}/parallel/${TEST_SIMULATOR}+${PARAM_CASENAME})
    set(TEST_ARGS ${OPM_TESTS_ROOT}/${PARAM_DIR}/${PARAM_FILENAME} ${PARAM_TEST_ARGS})

    set(DRIVER_ARGS -i ${OPM_TESTS_ROOT}/${PARAM_DIR}
                    -r ${RESULT_PATH}
                    -b ${PROJECT_BINARY_DIR}/bin
                    -f ${PARAM_FILENAME}
                    -a ${PARAM_ABS_TOL}
                    -t ${PARAM_REL_TOL}
                    -c $<TARGET_FILE:compareECL>
                    -n ${MPI_PROCS})

    # Add test that runs flow_mpi and outputs the results to file
    opm_add_test(compareParallelSim_${TEST_SIMULATOR}+${PARAM_FILENAME}${PARAM_POSTFIX} NO_COMPILE
                 EXE_NAME ${TEST_SIMULATOR}
                 DRIVER_ARGS ${DRIVER_ARGS}
                 TEST_ARGS ${TEST_ARGS})
    set_tests_properties(compareParallelSim_${TEST_SIMULATOR}+${PARAM_FILENAME}${PARAM_POSTFIX}
                         PROPERTIES PROCESSORS ${MPI_PROCS})

    label_blackoil_test(compareParallelSim_${TEST_SIMULATOR}+${PARAM_FILENAME}${PARAM_POSTFIX} "${TEST_SIMULATOR}" ${IS_BLACKOIL_DATA} blackoil blackoil-parallel blackoil-consistency)
  endif()
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

  test_supports_blackoil_simulator("${PARAM_DIR}" "${PARAM_FILENAME}" IS_BLACKOIL_DATA)
  resolve_blackoil_test_simulator("${PARAM_SIMULATOR}" ${IS_BLACKOIL_DATA} TEST_SIMULATOR)

  if(MPIEXEC_MAX_NUMPROCS GREATER_EQUAL MPI_PROCS)
    # Local computer system has at least ${MPI_PROCS} CPUs. Register test.
    set(RESULT_PATH ${BASE_RESULT_PATH}/parallelRestart/${TEST_SIMULATOR}+${PARAM_CASENAME})
    set(DRIVER_ARGS -i ${OPM_TESTS_ROOT}/${PARAM_DIR}
                    -r ${RESULT_PATH}
                    -b ${PROJECT_BINARY_DIR}/bin
                    -f ${PARAM_FILENAME}
                    -a ${PARAM_ABS_TOL}
                    -t ${PARAM_REL_TOL}
                    -c $<TARGET_FILE:compareECL>
                    -s ${PARAM_RESTART_STEP}
                    -d $<TARGET_FILE:rst_deck>
                    -n ${MPI_PROCS})

    opm_add_test(${TEST_NAME} NO_COMPILE
                    EXE_NAME ${TEST_SIMULATOR}
                 DRIVER_ARGS ${DRIVER_ARGS}
                 TEST_ARGS ${PARAM_TEST_ARGS})
    set_tests_properties(${TEST_NAME} PROPERTIES PROCESSORS ${MPI_PROCS})

    label_blackoil_test(${TEST_NAME} "${TEST_SIMULATOR}" ${IS_BLACKOIL_DATA} blackoil blackoil-parallel blackoil-restart blackoil-consistency)
  endif()
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

  test_supports_blackoil_simulator("${PARAM_DIR}" "${PARAM_FILENAME}" IS_BLACKOIL_DATA)
  resolve_blackoil_test_simulator("${PARAM_SIMULATOR}" ${IS_BLACKOIL_DATA} TEST_SIMULATOR)

  set(RESULT_PATH ${BASE_RESULT_PATH}/parallelSplitComm/${TEST_SIMULATOR}+${PARAM_CASENAME})
  set(DRIVER_ARGS -i ${OPM_TESTS_ROOT}/${PARAM_DIR}
                  -r ${RESULT_PATH}
                  -b ${PROJECT_BINARY_DIR}/bin
                  -f ${PARAM_FILENAME}
                  -a ${PARAM_ABS_TOL}
                  -t ${PARAM_REL_TOL}
                  -c $<TARGET_FILE:compareECL>
                  -n ${MPI_PROCS})

  opm_add_test(compareParallelSplitComm_${TEST_SIMULATOR}+${PARAM_FILENAME} NO_COMPILE
               EXE_NAME ${TEST_SIMULATOR}
               DRIVER_ARGS ${DRIVER_ARGS}
               TEST_ARGS ${PARAM_TEST_ARGS})
  set_tests_properties(compareParallelSplitComm_${TEST_SIMULATOR}+${PARAM_FILENAME}
                       PROPERTIES PROCESSORS ${MPI_PROCS})

  label_blackoil_test(compareParallelSplitComm_${TEST_SIMULATOR}+${PARAM_FILENAME} "${TEST_SIMULATOR}" ${IS_BLACKOIL_DATA} blackoil blackoil-parallel blackoil-consistency)
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

  test_supports_blackoil_simulator("${PARAM_DIR}" "${PARAM_FILENAME}" IS_BLACKOIL_DATA)
  resolve_blackoil_test_simulator("${PARAM_SIMULATOR}" ${IS_BLACKOIL_DATA} TEST_SIMULATOR)
  get_blackoil_reference_simulator("${TEST_SIMULATOR}" REFERENCE_SIMULATOR)

  set(RESULT_PATH ${BASE_RESULT_PATH}${PARAM_DIR_PREFIX}/${TEST_SIMULATOR}+${PARAM_CASENAME})
  set(TEST_ARGS ${PARAM_TEST_ARGS})
  set(DRIVER_ARGS -i ${OPM_TESTS_ROOT}/${PARAM_DIR}
                  -r ${RESULT_PATH}
                  -b ${PROJECT_BINARY_DIR}/bin
                  -f ${PARAM_FILENAME}
                  -a ${PARAM_ABS_TOL}
                  -t ${PARAM_REL_TOL}
                  -c ${HDF5_DIFF_EXECUTABLE}
                  -n ${MPI_PROCS}
                  -u ${REFERENCE_SIMULATOR})
  set(TEST_NAME ${PARAM_PREFIX}_${TEST_SIMULATOR}+${PARAM_FILENAME})
  opm_add_test(${TEST_NAME} NO_COMPILE
               EXE_NAME ${TEST_SIMULATOR}
               DRIVER_ARGS ${DRIVER_ARGS}
               TEST_ARGS ${TEST_ARGS})
  set_tests_properties(${TEST_NAME} PROPERTIES PROCESSORS ${MPI_PROCS}
                                    DIRNAME ${PARAM_DIR}
                                    FILENAME ${PARAM_FILENAME}
                                    SIMULATOR ${TEST_SIMULATOR}
                                    TESTNAME ${PARAM_CASENAME})

  label_blackoil_test(${TEST_NAME} "${TEST_SIMULATOR}" ${IS_BLACKOIL_DATA} blackoil blackoil-reference)
endfunction()

###########################################################################

# Adds several tests cases with similar parameters
# cases Variable name of list with test cases
# prefix Prefix to use
# argn Parameters for cases
macro(add_multiple_tests cases prefix)
  foreach(case ${${cases}})
    string(TOLOWER ${case} test)
    add_test_compareECLFiles(
        CASENAME ${prefix}${test}
        FILENAME ${case}
        ${ARGN}
    )
  endforeach()
endmacro()

###########################################################################

# Adds several tests cases in a numerical range with similar parameters
# start Start of range
# end End fof range
# ftemplate File name template to use
# prefix Prefix to use
# argn Parameters for cases
macro(add_multiple_test_range start end ftemplate prefix)
    foreach(case RANGE ${start} ${end})
      add_test_compareECLFiles(
          CASENAME ${prefix}_${case}
          FILENAME ${ftemplate}${case}
          ${ARGN}
      )
    endforeach()
endmacro()

###########################################################################

if(NOT TARGET test-suite)
  add_custom_target(test-suite)
endif()

# Simple execution tests
opm_set_test_driver(${PROJECT_SOURCE_DIR}/tests/run-test.sh "")
add_test_runSimulator(CASENAME norne
                      FILENAME NORNE_ATW2013
                      SIMULATOR ${OPM_BLACKOIL_TEST_SIMULATOR}
                      CONFIGURATION extra)

add_test_runSimulator(CASENAME norne_parallel
                      FILENAME NORNE_ATW2013
                      SIMULATOR ${OPM_BLACKOIL_TEST_SIMULATOR}
                      DIR norne
                      PROCS 4
                      CONFIGURATION extra)

add_test_runSimulator(CASENAME spe1case1_carfin
                      FILENAME SPE1CASE1_CARFIN
                      SIMULATOR ${OPM_BLACKOIL_TEST_SIMULATOR}
                      DIR lgr
                      TEST_ARGS --parsing-strictness=low --enable-ecl-output=true --enable-vtk-output=true)

add_test_runSimulator(CASENAME spe1case1_carfin_gr
                      FILENAME SPE1CASE1_CARFIN_GR
                      SIMULATOR ${OPM_BLACKOIL_TEST_SIMULATOR}
                      DIR lgr
                      TEST_ARGS --parsing-strictness=low --enable-ecl-output=true --enable-vtk-output=true)

if(MPI_FOUND)
  add_test_runSimulator(CASENAME spe1case1_carfin_parallel
                        FILENAME SPE1CASE1_CARFIN
                        SIMULATOR ${OPM_BLACKOIL_TEST_SIMULATOR}
                        DIR lgr
                        PROCS 4
                        TEST_ARGS --parsing-strictness=low --enable-ecl-output=false --enable-vtk-output=true)
endif()

add_test_runSimulator(CASENAME dryrun
                      FILENAME CO2STORE_PRECSALT
                      SIMULATOR ${OPM_BLACKOIL_TEST_SIMULATOR}
                      DIR co2store
                      TEST_ARGS --enable-dry-run=true --enable-ecl-output=false --enable-vtk-output=true)

# Tests that are run based on simulator results, but not necessarily direct comparison to reference results
add_test_runSimulator(CASENAME tuning_trgmbe
                      FILENAME 01_TUNING_TRGMBE
                      SIMULATOR ${OPM_BLACKOIL_TEST_SIMULATOR}
											DIR tuning
                      TEST_ARGS --output-extra-convergence-info=iterations --enable-tuning=true
                      POST_COMMAND $<TARGET_FILE:test_tuning_trgmbe>)

add_test_runSimulator(CASENAME notuning_trgmbe
                      FILENAME 01_TUNING_TRGMBE
                      SIMULATOR ${OPM_BLACKOIL_TEST_SIMULATOR}
											DIR tuning
                      TEST_ARGS --output-extra-convergence-info=iterations --enable-tuning=false
                      POST_COMMAND $<TARGET_FILE:test_tuning_trgmbe>)

set_tests_properties(runSimulator/notuning_trgmbe PROPERTIES WILL_FAIL TRUE)

add_test_runSimulator(CASENAME tuning_tsinit_nextstep
                      FILENAME 02_TUNING_TSINIT_NEXTSTEP
                      SIMULATOR ${OPM_BLACKOIL_TEST_SIMULATOR}
											DIR tuning
                      TEST_ARGS --enable-tuning=true
                      POST_COMMAND $<TARGET_FILE:test_tuning_tsinit_nextstep>)

get_property(opm-common_EMBEDDED_PYTHON TARGET opmcommon PROPERTY EMBEDDED_PYTHON)
if (opm-common_EMBEDDED_PYTHON)
  include (${CMAKE_CURRENT_SOURCE_DIR}/pyactionActionXComparisons.cmake)
endif ()
include (${CMAKE_CURRENT_SOURCE_DIR}/regressionTests.cmake)
include (${CMAKE_CURRENT_SOURCE_DIR}/comparisonTests.cmake)
include (${CMAKE_CURRENT_SOURCE_DIR}/restartTests.cmake)

# PORV test
opm_set_test_driver(${PROJECT_SOURCE_DIR}/tests/run-porv-acceptanceTest.sh "")
add_test_compareECLFiles(CASENAME norne
                         FILENAME NORNE_ATW2013
                         SIMULATOR ${OPM_BLACKOIL_TEST_SIMULATOR}
                         ABS_TOL 1e-5
                         REL_TOL 1e-8
                         PREFIX comparePORV
                         DIR_PREFIX /porv)

# Init tests
opm_set_test_driver(${PROJECT_SOURCE_DIR}/tests/run-init-regressionTest.sh "")

add_test_compareECLFiles(CASENAME norne_init
                         FILENAME NORNE_ATW2013
                         SIMULATOR ${OPM_BLACKOIL_TEST_SIMULATOR}
                         ABS_TOL ${abs_tol}
                         REL_TOL ${rel_tol}
                         PREFIX compareECLInitFiles
                         DIR norne
                         DIR_PREFIX /init)

set(_operate_work_tests
  OPERATE_ENDPOINTS-01
  OPERATER_ENDPOINTS-01
)

add_multiple_tests(_operate_work_tests
  "operate_work_"
  SIMULATOR ${OPM_BLACKOIL_TEST_SIMULATOR}
  ABS_TOL ${abs_tol}
  REL_TOL ${rel_tol}
  PREFIX compareECLInitFiles
  DIR operate
  DIR_PREFIX /init
)

# This is not a proper regression test; the test will load a norne case prepared
# for restart and run one single timestep - of length one day. The results are not
# verified in any way.
resolve_unlabeled_test_simulator("${OPM_BLACKOIL_TEST_SIMULATOR}" NORNE_RESTART_SIMULATOR)
add_test(NAME NORNE_RESTART
         COMMAND ${NORNE_RESTART_SIMULATOR} --output-dir=${BASE_RESULT_PATH}/norne-restart ${OPM_TESTS_ROOT}/norne/NORNE_ATW2013_RESTART.DATA)


if(MPI_FOUND)
  include (${CMAKE_CURRENT_SOURCE_DIR}/parallelRestartTests.cmake)

  # Single test to verify that we treat custom communicators correctly.
  opm_set_test_driver(${PROJECT_SOURCE_DIR}/tests/run-split-comm-test.sh "")
  add_test_split_comm(CASENAME spe1
                      FILENAME SPE1CASE2
                      SIMULATOR ${OPM_BLACKOIL_TEST_SIMULATOR}
                      ABS_TOL 0.0
                      REL_TOL 0.0)

  # Single test for damaris
  if(Damaris_FOUND AND USE_DAMARIS_LIB)
      opm_set_test_driver(${PROJECT_SOURCE_DIR}/tests/run-damaris-regressionTest.sh "")
      add_test_compareDamarisFiles(CASENAME spe1_damaris
                      FILENAME SPE1CASE1
                      SIMULATOR ${OPM_BLACKOIL_TEST_SIMULATOR}
                      ABS_TOL ${abs_tol}
                      REL_TOL 0.01
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
