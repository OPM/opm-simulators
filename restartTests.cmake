# Restart tests
opm_set_test_driver(${PROJECT_SOURCE_DIR}/tests/run-restart-regressionTest.sh "")
# Cruder tolerances for the restarted tests
set(abs_tol_restart 2e-1)
set(rel_tol_restart 4e-4)

file(GLOB_RECURSE tests "${OPM_TESTS_ROOT}/restartTests.json")
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

add_test_compare_restarted_simulation(CASENAME spe1
                                      FILENAME SPE1CASE2_ACTNUM
                                      SIMULATOR flow
                                      ABS_TOL ${abs_tol_restart}
                                      REL_TOL ${rel_tol_restart}
                                      RESTART_STEP 6
                                      TEST_ARGS --sched-restart=false)

add_test_compare_restarted_simulation(CASENAME spe9
                                      FILENAME SPE9_CP_SHORT
                                      SIMULATOR flow
                                      ABS_TOL ${abs_tol_restart}
                                      REL_TOL ${rel_tol_restart}
                                      RESTART_STEP 15
                                      TEST_ARGS --sched-restart=false  --tolerance-mb=1e-7)

add_test_compare_restarted_simulation(CASENAME ctaquifer_2d_oilwater
                                      FILENAME 2D_OW_CTAQUIFER
                                      SIMULATOR flow
                                      ABS_TOL ${abs_tol_restart}
                                      REL_TOL ${rel_tol_restart}
                                      DIR aquifer-oilwater
                                      RESTART_STEP 15
                                      TEST_ARGS --sched-restart=true)

add_test_compare_restarted_simulation(CASENAME fetkovich_2d
                                      FILENAME 2D_FETKOVICHAQUIFER
                                      SIMULATOR flow
                                      ABS_TOL ${abs_tol_restart}
                                      REL_TOL ${rel_tol_restart}
                                      RESTART_STEP 30
                                      DIR aquifer-fetkovich
                                      TEST_ARGS --sched-restart=true)

add_test_compare_restarted_simulation(CASENAME numerical_aquifer_3d_1aqu
                                      FILENAME 3D_1AQU_3CELLS
                                      SIMULATOR flow
                                      ABS_TOL 0.4
                                      REL_TOL 4.0e-3
                                      RESTART_STEP 3
                                      DIR aquifer-num
                                      TEST_ARGS --sched-restart=true --enable-tuning=true)

add_test_compare_restarted_simulation(CASENAME numerical_aquifer_3d_2aqu
                                      FILENAME 3D_2AQU_NUM
                                      SIMULATOR flow
                                      ABS_TOL 0.4
                                      REL_TOL 4.0e-3
                                      RESTART_STEP 3
                                      DIR aquifer-num
                                      TEST_ARGS --sched-restart=true --enable-tuning=true)

# Basic restart tests which only compare the summary output, this test driver should
# only be used in situations where it is challenging to get agreement in the restart file.
opm_set_test_driver(${PROJECT_SOURCE_DIR}/tests/run-summary-restart-regressionTest.sh "")

add_test_compare_restarted_simulation(CASENAME spe1
                                      FILENAME SPE1CASE2_ACTNUM
                                      SIMULATOR flow
                                      TEST_NAME restart_spe1_summary
                                      ABS_TOL ${abs_tol_restart}
                                      REL_TOL ${rel_tol_restart}
                                      RESTART_STEP 6
                                      TEST_ARGS --sched-restart=false)
