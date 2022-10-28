# Regression tests
opm_set_test_driver(${PROJECT_SOURCE_DIR}/tests/run-regressionTest.sh "")

function(AddRegressionTests directory)
  file(GLOB_RECURSE tests "${directory}/*.json")
  foreach(json_file ${tests})
    file(READ ${json_file} JSON_STRING)
    string(JSON size LENGTH ${JSON_STRING})
    math(EXPR COUNT "${size}-1")

    foreach(idx RANGE ${COUNT})
      parse_json_test_definition(${idx})
      add_test_compareECLFiles(${PARAMS})
    endforeach()
  endforeach()
endfunction()

AddRegressionTests("tests/definitions/regression")

if (opm-common_EMBEDDED_PYTHON)
  add_test_compareECLFiles(CASENAME udq_pyaction
                           FILENAME PYACTION_WCONPROD
                           SIMULATOR flow
                           ABS_TOL 2e-2
                           REL_TOL 1e-5
                           DIR udq_actionx)
endif()
