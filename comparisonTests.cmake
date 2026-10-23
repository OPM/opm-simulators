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

# The following tests are disabled until the decks CO2STORE_GW_MULTPV*.DATA
# are available in opm-tests.

# Test of the outlier cap in the relaxed CNV pore volume fraction test
# (--relaxed-pv-outlier-cap-multiplier). The pore volume fraction test is
# disabled by default and is therefore enabled explicitly here. A single cell
# with a large MULTPV mimics an aquifer. Without the cap, its pore volume
# hides a region of roughly 10% of the cells that violate the strict CNV
# tolerance. The relaxed CNV tolerance is then used there, and the solution
# differs from the reference by more than 2e-2 (relative) in FGIPL. With the
# cap the solution matches the reference case, where relaxed CNV is never
# used (XXXCNV = TRGCNV).
#add_test_compareSeparateECLFiles(
#  CASENAME
#    co2store_gw_multpv_relaxed_cnv
#  DIR1
#    co2store
#  FILENAME1
#    CO2STORE_GW_MULTPV
#  DIR2
#    co2store
#  FILENAME2
#    CO2STORE_GW_MULTPV_STRICT
#  SIMULATOR
#    flow
#  DEV_SIMULATOR
#    flow_gaswater_dissolution
#  ABS_TOL
#    ${abs_tol}
#  REL_TOL
#    1e-3
#  MPI_PROCS
#    1
#  TEST_ARGS
#    --enable-tuning=true
#    --relaxed-max-pv-fraction=0.03
#)

# Test of the pore volume floor in the per-cell CNV measure
# (--cnv-pv-floor-fraction, default 0.01). A layer of sliver cells with
# MULTPV = 1e-4 would otherwise dictate CNV convergence. With the floor the
# solution must still match the reference case, where the CNV measured with
# the cells' own pore volumes must satisfy the strict tolerance
# (XXXCNV = TRGCNV).
#add_test_compareSeparateECLFiles(
#  CASENAME
#    co2store_gw_multpv_sliver_cnv_pv_floor
#  DIR1
#    co2store
#  FILENAME1
#    CO2STORE_GW_MULTPV_SLIVER
#  DIR2
#    co2store
#  FILENAME2
#    CO2STORE_GW_MULTPV_SLIVER_STRICT
#  SIMULATOR
#    flow
#  DEV_SIMULATOR
#    flow_gaswater_dissolution
#  ABS_TOL
#    ${abs_tol}
#  REL_TOL
#    1e-3
#  MPI_PROCS
#    1
#  TEST_ARGS
#    --enable-tuning=true
#)
