# Parallel tests
opm_set_test_driver(${PROJECT_SOURCE_DIR}/tests/run-parallel-restart-regressionTest.sh "")

file(GLOB_RECURSE tests "${OPM_TESTS_ROOT}/parallelRestartTests.json")
  foreach(json_file ${tests})
  file(READ ${json_file} JSON_STRING)
  parse_json_common_params()
  string(JSON size LENGTH ${test_cases})
  math(EXPR COUNT "${size}-1")
  foreach(idx RANGE ${COUNT})
    parse_json_test_definition(${idx})
    add_test_compare_parallel_restarted_simulation(${PARAMS})
  endforeach()
endforeach()

add_test_compare_parallel_restarted_simulation(CASENAME ctaquifer_2d_oilwater
                                               FILENAME 2D_OW_CTAQUIFER
                                               SIMULATOR flow
                                               ABS_TOL ${abs_tol_restart}
                                               REL_TOL ${rel_tol_restart}
                                               RESTART_STEP 15
                                               DIR aquifer-oilwater
                                               TEST_ARGS --enable-tuning=true --linear-solver-reduction=1e-7 --tolerance-cnv=5e-6 --tolerance-mb=1e-6)
add_test_compare_parallel_restarted_simulation(CASENAME fetkovich_2d
                                               FILENAME 2D_FETKOVICHAQUIFER
                                               SIMULATOR flow
                                               ABS_TOL ${abs_tol_restart}
                                               REL_TOL ${rel_tol_restart}
                                               RESTART_STEP 30
                                               DIR aquifer-fetkovich
                                               TEST_ARGS --enable-tuning=true --linear-solver-reduction=1e-7 --tolerance-cnv=5e-6 --tolerance-mb=1e-6)

add_test_compare_parallel_restarted_simulation(CASENAME numerical_aquifer_3d_2aqu
                                               FILENAME 3D_2AQU_NUM
                                               SIMULATOR flow
                                               ABS_TOL 0.12
                                               REL_TOL 5.0e-2
                                               RESTART_STEP 3
                                               DIR aquifer-num
                                               TEST_ARGS --enable-tuning=true --tolerance-cnv=0.00003 --time-step-control=pid --linear-solver=cpr)

add_test_compare_parallel_restarted_simulation(CASENAME numerical_aquifer_3d_1aqu
                                               FILENAME 3D_1AQU_3CELLS
                                               SIMULATOR flow
                                               ABS_TOL 0.12
                                               REL_TOL 5.0e-2
                                               RESTART_STEP 3
                                               DIR aquifer-num
                                               TEST_ARGS --enable-tuning=true --tolerance-cnv=0.00003 --time-step-control=pid --linear-solver=cpr)
