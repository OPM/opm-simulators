# Regression tests
opm_set_test_driver(${PROJECT_SOURCE_DIR}/tests/run-comparison.sh "")

# Use same tolerances as in regressionTests
set(abs_tol 2e-2)
set(rel_tol 1e-5)
set(coarse_rel_tol 1e-2)

# Tests comparing VFPTABLE oneliner vs multiliner lift curves (constant delta pressure)
add_test_compareSeparateECLFiles(
  CASENAME
    spe1_metric_vfp1_multiliner_vs_oneliner
  DIR1
    vfpprod_spe1_oneliner
  FILENAME1
    SPE1CASE1_METRIC_VFP1_MULTILINER
  DIR2
    vfpprod_spe1_oneliner
  FILENAME2
    SPE1CASE1_METRIC_VFP1_ONELINER
  SIMULATOR
    flow
  DEV_SIMULATOR
    flow_blackoil
  ABS_TOL
    ${abs_tol}
  REL_TOL
    ${rel_tol}
  IGNORE_EXTRA_KW
    BOTH
  MPI_PROCS
    1
)

add_test_compareSeparateECLFiles(
  CASENAME
    1dcompositional_flow_comp3_2p_vs_flowexp_comp3_2p
  DIR1
    compositional
  FILENAME1
    1D_COMP
  DIR2
    compositional
  FILENAME2
    1D_COMP
  SIMULATOR
    flow_comp3_2p
  REFERENCE_SIMULATOR
    flowexp_comp3_2p
  ABS_TOL
    ${abs_tol}
  REL_TOL
    ${rel_tol}
  TEST_ARGS
    --enable-adaptive-time-stepping=false
  MPI_PROCS
    1
)

add_test_compareSeparateECLFiles(
  CASENAME
    1dcompositional_prewell_flow_comp3_2p_vs_flowexp_comp3_2p
  DIR1
    compositional
  FILENAME1
    1D_COMP_PREWELL
  DIR2
    compositional
  FILENAME2
    1D_COMP_PREWELL
  SIMULATOR
    flow_comp3_2p
  REFERENCE_SIMULATOR
    flowexp_comp3_2p
  ABS_TOL
    ${abs_tol}
  REL_TOL
    ${rel_tol}
  TEST_ARGS
    --enable-adaptive-time-stepping=false
  MPI_PROCS
    1
)
