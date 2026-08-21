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

# Refining part of the grid must not change the field totals.  The two decks are
# identical but for the CARFIN block, so any difference in FPR/FOIP/FGIP/FWIP/
# FRPV/FHPV means the region arrays or the in-place sums are not mapped onto the
# leaf grid correctly.  The tolerance sits an order of magnitude above the
# discretisation drift between the two grids (~1.3e-4 by the end of the run) and
# an order below the ~1.4e-2 error it is there to catch.
set(lgr_rel_tol 1e-3)
add_test_compareSeparateECLFiles(
  CASENAME
    lgr_field_totals_invariant_under_refinement
  DIR1
    lgr
  FILENAME1
    SPE1CASE1_CARFIN1-3DCORNERPOINT_XYZ
  DIR2
    lgr
  FILENAME2
    SPE1CASE1_CARFIN1-3DCORNERPOINT_XYZ_NOLGR
  SIMULATOR
    flow
  DEV_SIMULATOR
    flow_blackoil
  ABS_TOL
    ${abs_tol}
  REL_TOL
    ${lgr_rel_tol}
  IGNORE_EXTRA_KW
    BOTH
  MPI_PROCS
    1
  TEST_ARGS
    --parsing-strictness=low
)
