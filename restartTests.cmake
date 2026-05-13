# Restart tests
opm_set_test_driver(${PROJECT_SOURCE_DIR}/tests/run-restart-regressionTest.sh "")
# Cruder tolerances for the restarted tests
set(abs_tol_restart 2e-1)
set(rel_tol_restart 4e-4)

# This is not a proper regression test; the test will load a norne case prepared
# for restart and run one single timestep - of length one day. The results are not
# verified in any way.
if(USE_DEV_SIMULATOR_IN_TESTS)
  set(RESTART_SIM flow_blackoil)
else()
  set(RESTART_SIM flow)
endif()

add_test(
  NAME
    NORNE_RESTART
  COMMAND
    ${RESTART_SIM}
    --output-dir=${BASE_RESULT_PATH}/norne-restart
    ${OPM_TESTS_ROOT}/norne/NORNE_ATW2013_RESTART.DATA
)

add_test_compare_restarted_simulation(
  CASENAME
    spe1
  FILENAME
    SPE1CASE2_ACTNUM
  SIMULATOR
    flow
  DEV_SIMULATOR
    flow_blackoil
  ABS_TOL
    ${abs_tol_restart}
  REL_TOL
    ${rel_tol_restart}
  RESTART_STEP
    6
  TEST_ARGS
    --sched-restart=false
)

add_test_compare_restarted_simulation(
  CASENAME
    spe9
  FILENAME
    SPE9_CP_SHORT
  SIMULATOR
    flow
  DEV_SIMULATOR
    flow_blackoil
  ABS_TOL
    ${abs_tol_restart}
  REL_TOL
    ${rel_tol_restart}
  RESTART_STEP
    15
  TEST_ARGS
    --sched-restart=false
    --tolerance-mb=1e-7
)

add_test_compare_restarted_simulation(
  CASENAME
    ctaquifer_2d_oilwater
  FILENAME
    2D_OW_CTAQUIFER
  SIMULATOR
    flow
  DEV_SIMULATOR
    flow_oilwater
  ABS_TOL
    ${abs_tol_restart}
  REL_TOL
    ${rel_tol_restart}
  DIR
    aquifer-oilwater
  RESTART_STEP
    15
  TEST_ARGS
    --sched-restart=true
)

add_test_compare_restarted_simulation(
  CASENAME
    fetkovich_2d
  FILENAME
    2D_FETKOVICHAQUIFER
  SIMULATOR
    flow
  DEV_SIMULATOR
    flow_oilwater
  ABS_TOL
    ${abs_tol_restart}
  REL_TOL
    ${rel_tol_restart}
  RESTART_STEP
    30
  DIR
    aquifer-fetkovich
  TEST_ARGS
    --sched-restart=true
)

add_test_compare_restarted_simulation(
  CASENAME
    numerical_aquifer_3d_1aqu
  FILENAME
    3D_1AQU_3CELLS
  SIMULATOR
    flow
  DEV_SIMULATOR
    flow_blackoil
  ABS_TOL
    0.4
  REL_TOL
    4.0e-3
  RESTART_STEP
    3
  DIR
    aquifer-num
  TEST_ARGS
    --enable-tuning=true
    --relaxed-max-pv-fraction=0.0
    --enable-drift-compensation=false
)

add_test_compare_restarted_simulation(
  CASENAME
    numerical_aquifer_3d_2aqu
  FILENAME
    3D_2AQU_NUM
  SIMULATOR
    flow
  DEV_SIMULATOR
    flow_blackoil
  ABS_TOL
    0.4
  REL_TOL
    4.0e-3
  RESTART_STEP
    3
  DIR
    aquifer-num
  TEST_ARGS
    --sched-restart=true
    --enable-tuning=true
)

add_test_compare_restarted_simulation(
  CASENAME
    aquflux_01
  FILENAME
    AQUFLUX-01
  SIMULATOR
    flow
  DEV_SIMULATOR
    flow_blackoil
  ABS_TOL
    ${abs_tol_restart}
  REL_TOL
    3.0e-3
  RESTART_STEP
    3
  DIR
    aquifers
  TEST_ARGS
    --enable-tuning=true
)

add_test_compare_restarted_simulation(
  CASENAME
    aquflux_02
  FILENAME
    AQUFLUX-02
  SIMULATOR
    flow
  DEV_SIMULATOR
    flow_oilwater
  ABS_TOL
    ${abs_tol_restart}
  REL_TOL
    ${rel_tol_restart}
  RESTART_STEP
    50
  DIR
    aquifers
  TEST_ARGS
    --enable-tuning=true
)

add_test_compare_restarted_simulation(
  CASENAME
    network_01_restart
  FILENAME
    NETWORK-01-RESTART
  SIMULATOR
    flow
  DEV_SIMULATOR
    flow_blackoil
  ABS_TOL
    ${abs_tol_restart}
  REL_TOL
    ${rel_tol_restart}
  RESTART_STEP
    5
  DIR
    network
  TEST_ARGS
    --enable-tuning=true
    --local-well-solve-control-switching=true
)

add_test_compare_restarted_simulation(
  CASENAME
    network_01_reroute_restart
  FILENAME
    NETWORK-01-REROUTE-RESTART
  SIMULATOR
    flow
  DEV_SIMULATOR
    flow_blackoil
  ABS_TOL
    ${abs_tol_restart}
  REL_TOL
    ${rel_tol_restart}
  RESTART_STEP
    5
  DIR
    network
  TEST_ARGS
    --enable-tuning=true
    --local-well-solve-control-switching=true
)

add_test_compare_restarted_simulation(
  CASENAME
    spe1_temp
  FILENAME
    SPE1CASE2_TEMP
  SIMULATOR
    flow
  DEV_SIMULATOR
    flow_blackoil_temp
  ABS_TOL
    ${abs_tol_restart}
  REL_TOL
    5.0e-2
  RESTART_STEP
    3
  DIR
    spe1
  TEST_ARGS
    --solver-max-time-step-in-days=1
)

# Restart run in which a UDQ defining expression has exactly 128
# characters.  Verifies that we don't overflow the ZUDL character
# limit in restart files.
add_test_compare_restarted_simulation(
  CASENAME
    udq_reg_02
  FILENAME
    UDQ_REG-02
  SIMULATOR
    flow
  DEV_SIMULATOR
    flow_blackoil
  ABS_TOL
    ${abs_tol_restart}
  REL_TOL
    ${rel_tol_restart}
  RESTART_STEP
    2
  DIR
    udq_actionx
  TEST_ARGS
    --enable-tuning=true
)

# The dynamic MSW data is not written to /read from the restart file
# We therefore accept significant deviation in the results.
# Note also that we use --sched-restart=true since some necessary
# MSW info is still lacking in the restart file.
set(abs_tol_restart_msw 2e2)
set(rel_tol_restart_msw 1e-3)

add_test_compare_restarted_simulation(
  CASENAME
    msw_3d_hfa
  FILENAME
    3D_MSW
  SIMULATOR
    flow
  DEV_SIMULATOR
    flow_blackoil
  ABS_TOL
    ${abs_tol_restart_msw}
  REL_TOL
    ${rel_tol_restart_msw}
  RESTART_STEP
    10
  TEST_ARGS
    --enable-adaptive-time-stepping=false
    --sched-restart=true
    --tolerance-wells=1e-7
)


# Basic restart tests which only compare the summary output, this test driver should
# only be used in situations where it is challenging to get agreement in the restart file.
opm_set_test_driver(${PROJECT_SOURCE_DIR}/tests/run-summary-restart-regressionTest.sh "")

add_test_compare_restarted_simulation(
  CASENAME
    spe1_actnum
  FILENAME
    SPE1CASE2_ACTNUM
  SIMULATOR
    flow
  DEV_SIMULATOR
    flow_blackoil
  TEST_NAME
    restart_spe1_summary
  ABS_TOL
    ${abs_tol_restart}
  REL_TOL
    ${rel_tol_restart}
  RESTART_STEP
    6
  DIR
    spe1
  TEST_ARGS
    --sched-restart=false
)

# Serialized restart tests
if(HDF5_FOUND)
  opm_set_test_driver(${PROJECT_SOURCE_DIR}/tests/run-serialization-regressionTest.sh "")
  add_test_compare_restarted_simulation(
    CASENAME
      spe1_serialized
    DIR
      spe1
    FILENAME
      SPE1CASE1
    SIMULATOR
      flow
    DEV_SIMULATOR
      flow_blackoil
    TEST_NAME
      compareSerializedSim_${RESTART_SIM}+spe1
    ABS_TOL
      2e-2
    REL_TOL
      1e-5
    RESTART_STEP
      94
    TEST_ARGS
      --tolerance-mb=1e-7
  )
endif()
