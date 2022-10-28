# Restart tests
opm_set_test_driver(${PROJECT_SOURCE_DIR}/tests/run-restart-regressionTest.sh "")

function(AddRestartTests directory)
  file(GLOB_RECURSE tests "${directory}/*.json")
  foreach(json_file ${tests})
    file(READ ${json_file} JSON_STRING)
    parse_json_common_params()
    string(JSON size LENGTH ${test_cases})
    math(EXPR COUNT "${size}-1")
    foreach(idx RANGE ${COUNT})
      parse_json_test_definition(${idx})
      add_test_compare_restarted_simulation(${PARAMS})
    endforeach()
  endforeach()
endfunction()

AddRestartTests("tests/definitions/restart")

# Basic restart tests which only compare the summary output, this test driver should
# only be used in situations where it is challenging to get agreement in the restart file.
opm_set_test_driver(${PROJECT_SOURCE_DIR}/tests/run-summary-restart-regressionTest.sh "")

add_test_compare_restarted_simulation(CASENAME spe1
                                      FILENAME SPE1CASE2_ACTNUM
                                      SIMULATOR flow
                                      TEST_NAME restart_spe1_summary
                                      ABS_TOL 2e-1
                                      REL_TOL 4e-4
                                      RESTART_STEP 6
                                      TEST_ARGS --sched-restart=false)
