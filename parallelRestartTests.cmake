# Parallel tests
opm_set_test_driver(${PROJECT_SOURCE_DIR}/tests/run-parallel-restart-regressionTest.sh "")
add_test_compare_parallel_restarted_simulation(CASENAME spe1
                                               FILENAME SPE1CASE2_ACTNUM
                                               SIMULATOR flow
                                               ABS_TOL ${abs_tol_restart}
                                               REL_TOL ${rel_tol_restart}
                                               RESTART_STEP 6
                                               TEST_ARGS --sched-restart=false --enable-adaptive-time-stepping=false)

add_test_compare_parallel_restarted_simulation(CASENAME ctaquifer_2d_oilwater
                                               FILENAME 2D_OW_CTAQUIFER
                                               SIMULATOR flow
                                               ABS_TOL ${abs_tol_restart}
                                               REL_TOL ${rel_tol_restart}
                                               RESTART_STEP 15
                                               DIR aquifer-oilwater
                                               TEST_ARGS --enable-tuning=true --linear-solver-reduction=1e-7 --tolerance-cnv=5e-6 --tolerance-mb=1e-6 --enable-adaptive-time-stepping=false)

add_test_compare_parallel_restarted_simulation(CASENAME fetkovich_2d
                                               FILENAME 2D_FETKOVICHAQUIFER
                                               SIMULATOR flow
                                               ABS_TOL ${abs_tol_restart}
                                               REL_TOL ${rel_tol_restart}
                                               RESTART_STEP 30
                                               DIR aquifer-fetkovich
                                               TEST_ARGS --enable-tuning=true --linear-solver-reduction=1e-7 --tolerance-cnv=5e-6 --tolerance-mb=1e-6 --enable-adaptive-time-stepping=false)

add_test_compare_parallel_restarted_simulation(CASENAME numerical_aquifer_3d_2aqu
                                               FILENAME 3D_2AQU_NUM
                                               SIMULATOR flow
                                               ABS_TOL 0.12
                                               REL_TOL 5.0e-2
                                               RESTART_STEP 3
                                               DIR aquifer-num
                                               TEST_ARGS --enable-tuning=true --tolerance-cnv=0.00003 --time-step-control=pid --linear-solver=cpr_trueimpes --enable-adaptive-time-stepping=false)

add_test_compare_parallel_restarted_simulation(CASENAME numerical_aquifer_3d_1aqu
                                               FILENAME 3D_1AQU_3CELLS
                                               SIMULATOR flow
                                               ABS_TOL 0.12
                                               REL_TOL 5.0e-2
                                               RESTART_STEP 3
                                               DIR aquifer-num
                                               TEST_ARGS --enable-tuning=true --tolerance-cnv=0.00003 --time-step-control=pid --linear-solver=cpr_trueimpes --enable-adaptive-time-stepping=false)

add_test_compare_parallel_restarted_simulation(CASENAME aquflux_01
                                               FILENAME AQUFLUX-01
                                               SIMULATOR flow
                                               ABS_TOL ${abs_tol_restart}
                                               REL_TOL 5.0e-2
                                               RESTART_STEP 3
                                               DIR aquifers
                                               TEST_ARGS --solver-max-time-step-in-days=1)

add_test_compare_parallel_restarted_simulation(CASENAME aquflux_02
                                               FILENAME AQUFLUX-02
                                               SIMULATOR flow
                                               ABS_TOL ${abs_tol_restart}
                                               REL_TOL 5.0e-2
                                               RESTART_STEP 50
                                               DIR aquifers
                                               TEST_ARGS --solver-max-time-step-in-days=1)

# Serialized restart tests
if(HDF5_FOUND)
  opm_set_test_driver(${PROJECT_SOURCE_DIR}/tests/run-serialization-regressionTest.sh "")
  add_test_compare_parallel_restarted_simulation(CASENAME spe1_serialized
                                                 DIR spe1
                                                 FILENAME SPE1CASE1
                                                 SIMULATOR flow
                                                 TEST_NAME compareParallelSerializedSim_flow+spe1
                                                 ABS_TOL 2e-2
                                                 REL_TOL 1e-5
                                                 RESTART_STEP 94
                                                 TEST_ARGS --tolerance-mb=1e-7
                                                 MPI_PROCS 4)
endif()
