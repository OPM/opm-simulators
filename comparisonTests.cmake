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

# Dual-continuum runs: a corner-point (COORD/ZCORN) deck must give the same
# results as the block-centred (DX/DY/DZ/TOPS) deck describing the same grid.
foreach(model dp dk)
  string(TOUPPER ${model} MODEL)
  add_test_compareSeparateECLFiles(
    CASENAME
      dualporo_${model}_cornerpoint_vs_blockcentred
    DIR1
      dualporo
    FILENAME1
      ${MODEL}_BASIC
    DIR2
      dualporo
    FILENAME2
      ${MODEL}_CP1
    SIMULATOR
      flow_oilwater
    ABS_TOL
      ${abs_tol}
    REL_TOL
      ${rel_tol}
    MPI_PROCS
      1
  )
endforeach()

# Dual-continuum runs on a corner-point grid with inactive cells, MINPV and PINCH
# (model1 with two continua). With no flow between the continua, the cells of a continuum
# must behave as the single porosity model made of them. The time steps are fixed and the
# tolerances tight, so that the runs can be compared cell by cell.
opm_set_test_driver(${PROJECT_SOURCE_DIR}/tests/run-continuum-comparison.sh "")
set(continuum_test_args
  --solver-max-time-step-in-days=1
  --tolerance-cnv=1e-9
  --tolerance-cnv-relaxed=1e-9
  --tolerance-mb=1e-12
  --linear-solver-reduction=1e-12)
foreach(pinch "" _ALL)
  add_test_compareContinuumToSinglePorosity(
    CASENAME
      dualporo_model1_fracture_continuum${pinch}
    DIR1
      dualporo_model1
    FILENAME1
      M1_SP${pinch}
    DIR2
      dualporo_model1
    FILENAME2
      M1_DP${pinch}
    CONTINUUM
      fracture
    SIMULATOR
      flow_blackoil
    TEST_ARGS
      ${continuum_test_args}
  )
  add_test_compareContinuumToSinglePorosity(
    CASENAME
      dualporo_model1_matrix_continuum${pinch}
    DIR1
      dualporo_model1
    FILENAME1
      M1_SPM${pinch}
    DIR2
      dualporo_model1
    FILENAME2
      M1_DKM${pinch}
    CONTINUUM
      matrix
    SIMULATOR
      flow_blackoil
    TEST_ARGS
      ${continuum_test_args}
  )
endforeach()

# The matrix-fracture coupling together with inactive cells, MINPV and PINCH: a twin
# removed by MINPV takes its coupling with it.
opm_set_test_driver(${PROJECT_SOURCE_DIR}/tests/run-test.sh "")
foreach(case M1_DP_SIG M1_DK_ALL_SIG)
  add_test_runSimulator(
    CASENAME
      dualporo_model1_coupled_${case}
    DIR
      dualporo_model1
    FILENAME
      ${case}
    SIMULATOR
      flow_blackoil
  )
endforeach()
