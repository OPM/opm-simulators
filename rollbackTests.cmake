# Deterministic timestep-replay regression tests
#
# Run a case, record the accepted substep end times, then re-run with the
# hardcoded timestep controller replaying those times and compare the output.
# --truncate-time-step-to-float=true makes the recorded step sizes exactly
# reproducible on replay (see Simulator::setTimeStepSize).

set(timestep_replay_abs_tol 2e-14)
set(timestep_replay_rel_tol 2e-14)

opm_set_test_driver(${PROJECT_SOURCE_DIR}/tests/run-timestep-replay-regressionTest.sh "")

add_test_compareECLFiles(
  CASENAME
    spe1_timestep_replay
  FILENAME
    SPE1CASE1
  SIMULATOR
    flow
  DEV_SIMULATOR
    flow_blackoil
  PREFIX
    compareTimestepReplay
  ABS_TOL
    ${timestep_replay_abs_tol}
  REL_TOL
    ${timestep_replay_rel_tol}
  DIR
    spe1
  TEST_ARGS
    --enable-state-rollback=true
    --truncate-time-step-to-float=true
  TEST_ARGS_REPLAY
    --initial-time-step-in-days=11111111
)

add_test_compareECLFiles(
  CASENAME
    spe9_timestep_replay
  FILENAME
    SPE9_CP_SHORT
  SIMULATOR
    flow
  DEV_SIMULATOR
    flow_blackoil
  PREFIX
    compareTimestepReplay
  ABS_TOL
    ${timestep_replay_abs_tol}
  REL_TOL
    ${timestep_replay_rel_tol}
  DIR
    spe9
  TEST_ARGS
    --enable-state-rollback=true
    --truncate-time-step-to-float=true
  TEST_ARGS_REPLAY
    --initial-time-step-in-days=11111111
)

# Larger model2 case: forcing the full report step initially exercises failed
# steps and the per-step state rollback before the replay must still match.
add_test_compareECLFiles(
  CASENAME
    model2_base_timestep_replay
  FILENAME
    0_BASE_MODEL2
  SIMULATOR
    flow
  DEV_SIMULATOR
    flow_blackoil
  PREFIX
    compareTimestepReplay
  ABS_TOL
    ${timestep_replay_abs_tol}
  REL_TOL
    ${timestep_replay_rel_tol}
  DIR
    model2
  TEST_ARGS
    --enable-state-rollback=true
    --truncate-time-step-to-float=true
    --full-time-step-initially=true
    --cpr-reuse-setup=0
    --newton-max-iterations=8
    --tolerance-cnv-relaxed=1e-3
  TEST_ARGS_REPLAY
    --initial-time-step-in-days=11111111
)

# ACTIONX dynamic event trigger cases
add_test_compareECLFiles(
  CASENAME
    actionx_gconprod_timestep_replay
  FILENAME
    ACTIONX_GCONPROD
  SIMULATOR
    flow
  DEV_SIMULATOR
    flow_blackoil
  PREFIX
    compareTimestepReplay
  ABS_TOL
    ${timestep_replay_abs_tol}
  REL_TOL
    ${timestep_replay_rel_tol}
  DIR
    actionx
  TEST_ARGS
    --enable-state-rollback=true
    --truncate-time-step-to-float=true
  TEST_ARGS_REPLAY
    --initial-time-step-in-days=11111111
)

add_test_compareECLFiles(
  CASENAME
    actionx_udq_timestep_replay
  FILENAME
    ACTIONX_UDQ
  SIMULATOR
    flow
  DEV_SIMULATOR
    flow_blackoil
  PREFIX
    compareTimestepReplay
  ABS_TOL
    ${timestep_replay_abs_tol}
  REL_TOL
    ${timestep_replay_rel_tol}
  DIR
    actionx
  TEST_ARGS
    --enable-state-rollback=true
    --truncate-time-step-to-float=true
  TEST_ARGS_REPLAY
    --initial-time-step-in-days=11111111
)

# Relative permeability hysteresis (EHYSTR) case
add_test_compareECLFiles(
  CASENAME
    model2_hysteresis_timestep_replay
  FILENAME
    7_HYSTERESIS_MODEL2
  SIMULATOR
    flow
  DEV_SIMULATOR
    flow_blackoil
  PREFIX
    compareTimestepReplay
  ABS_TOL
    ${timestep_replay_abs_tol}
  REL_TOL
    ${timestep_replay_rel_tol}
  DIR
    model2
  TEST_ARGS
    --enable-state-rollback=true
    --truncate-time-step-to-float=true
    --cpr-reuse-setup=0
  TEST_ARGS_REPLAY
    --initial-time-step-in-days=11111111
)

# Rock compaction hysteresis (ROCKCOMP) case
add_test_compareECLFiles(
  CASENAME
    spe1_rock2dtr_timestep_replay
  FILENAME
    SPE1CASE2_ROCK2DTR
  SIMULATOR
    flow
  DEV_SIMULATOR
    flow_blackoil
  PREFIX
    compareTimestepReplay
  ABS_TOL
    ${timestep_replay_abs_tol}
  REL_TOL
    ${timestep_replay_rel_tol}
  DIR
    spe1
  TEST_ARGS
    --enable-state-rollback=true
    --truncate-time-step-to-float=true
    --full-time-step-initially=true
  TEST_ARGS_REPLAY
    --initial-time-step-in-days=11111111
)

# Well cycling (WCYCLE) timer case
add_test_compareECLFiles(
  CASENAME
    wcycle1_timestep_replay
  FILENAME
    WCYCLE-1
  SIMULATOR
    flow
  DEV_SIMULATOR
    flow_blackoil
  PREFIX
    compareTimestepReplay
  ABS_TOL
    ${timestep_replay_abs_tol}
  REL_TOL
    ${timestep_replay_rel_tol}
  DIR
    wcycle
  TEST_ARGS
    --enable-state-rollback=true
    --truncate-time-step-to-float=true
  TEST_ARGS_REPLAY
    --initial-time-step-in-days=11111111
)

# Multisegment Wells (MSW) case
add_test_compareECLFiles(
  CASENAME
    msw_model3_timestep_replay
  FILENAME
    2_MSWWELL_MODEL3
  SIMULATOR
    flow
  DEV_SIMULATOR
    flow_blackoil
  PREFIX
    compareTimestepReplay
  ABS_TOL
    ${timestep_replay_abs_tol}
  REL_TOL
    ${timestep_replay_rel_tol}
  DIR
    model3
  TEST_ARGS
    --enable-state-rollback=true
    --truncate-time-step-to-float=true
    --cpr-reuse-setup=0
  TEST_ARGS_REPLAY
    --initial-time-step-in-days=11111111
)

# Group hierarchy and dynamic rate allocation (MOD4_GRP) case
add_test_compareECLFiles(
  CASENAME
    mod4_grp_timestep_replay
  FILENAME
    MOD4_GRP
  SIMULATOR
    flow
  DEV_SIMULATOR
    flow_blackoil
  PREFIX
    compareTimestepReplay
  ABS_TOL
    ${timestep_replay_abs_tol}
  REL_TOL
    ${timestep_replay_rel_tol}
  DIR
    model4
  TEST_ARGS
    --enable-state-rollback=true
    --truncate-time-step-to-float=true
  TEST_ARGS_REPLAY
    --initial-time-step-in-days=11111111
)
