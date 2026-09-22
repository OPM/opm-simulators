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


###########################################################################
# Numerical aquifers: grid vs aux mode compared against each other, not reference data
###########################################################################

opm_set_test_driver(${PROJECT_SOURCE_DIR}/tests/run-numerical-aquifer-mode-comparison.sh "")

# AQUNUM-02 requests BPR at aquifer cells, which exist only in grid mode; hence -x.
foreach(case AQUNUM-01 AQUNUM-02 AQUNUM-03 AQUNUM-04)
  string(TOLOWER ${case} test)
  set(aquifer_mode_extra_args "")
  if(${case} STREQUAL AQUNUM-02)
    set(aquifer_mode_extra_args -x)
  endif()

  opm_add_test(compareNumericalAquiferModes_flow+${test}
    EXE_TARGET
      flow
    DRIVER_ARGS
      -i ${OPM_TESTS_ROOT}/aquifers
      -f ${case}
      -r ${BASE_RESULT_PATH}/flow+aquifer_modes_${test}
      -a ${abs_tol}
      -t ${rel_tol}
      -c $<TARGET_FILE:compareECL>
      ${aquifer_mode_extra_args}
  )
endforeach()
