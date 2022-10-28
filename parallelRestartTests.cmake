# Parallel tests
opm_set_test_driver(${PROJECT_SOURCE_DIR}/tests/run-parallel-restart-regressionTest.sh "")

function(AddParallelRestartTests directory)
  file(GLOB_RECURSE tests "${directory}/*.json")
  foreach(json_file ${tests})
    file(READ ${json_file} JSON_STRING)
    string(JSON size LENGTH ${JSON_STRING})
    math(EXPR COUNT "${size}-1")

    foreach(idx RANGE ${COUNT})
      parse_json_test_definition(${idx})
      add_test_compare_parallel_restarted_simulation(${PARAMS})
    endforeach()
  endforeach()
endfunction()

AddParallelRestartTests("tests/definitions/parallel_restart")
