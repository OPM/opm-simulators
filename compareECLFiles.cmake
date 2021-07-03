# This script manages the addition of tests.
# The tests are orchestrated by a shell script,
# configured using opm_set_test_driver()
# and then the appropriate helper macro is called to
# register the ctest entry through the opm_add_test macro.
# Information such as the binary to call and test tolerances
# are passed from the build system to the driver script through
# command line parameters. See the opm_add_test() documentation for
# details on the parameters passed to the macro.

# Define some paths
set(BASE_RESULT_PATH ${PROJECT_BINARY_DIR}/tests/results)

###########################################################################
# TEST: runSim
###########################################################################

# Input:
#   - casename: basename (no extension)
#
# Details:
#   - This test class simply runs a simulation.
function(add_test_runSimulator)
  set(oneValueArgs CASENAME FILENAME SIMULATOR DIR DIR_PREFIX PROCS)
  set(multiValueArgs TEST_ARGS)
  cmake_parse_arguments(PARAM "$" "${oneValueArgs}" "${multiValueArgs}" ${ARGN} )
  if(NOT PARAM_DIR)
    set(PARAM_DIR ${PARAM_CASENAME})
  endif()
  set(RESULT_PATH ${BASE_RESULT_PATH}${PARAM_DIR_PREFIX}/${PARAM_SIMULATOR}+${PARAM_CASENAME})
  set(TEST_ARGS ${OPM_TESTS_ROOT}/${PARAM_DIR}/${PARAM_FILENAME} ${PARAM_TEST_ARGS})
  opm_add_test(runSimulator/${PARAM_CASENAME} NO_COMPILE
               EXE_NAME ${PARAM_SIMULATOR}
               DRIVER_ARGS ${OPM_TESTS_ROOT}/${PARAM_DIR}
                           ${RESULT_PATH}
                           ${PROJECT_BINARY_DIR}/bin
                           ${PARAM_FILENAME}
                           ${PARAM_PROCS}
               TEST_ARGS ${TEST_ARGS}
               CONFIGURATION extra)
endfunction()

###########################################################################
# TEST: compareECLFiles
###########################################################################

# Input:
#   - casename: basename (no extension)
#
# Details:
#   - This test class compares output from a simulation to reference files.
function(add_test_compareECLFiles)
  set(oneValueArgs CASENAME FILENAME SIMULATOR ABS_TOL REL_TOL DIR DIR_PREFIX PREFIX)
  set(multiValueArgs TEST_ARGS)
  cmake_parse_arguments(PARAM "$" "${oneValueArgs}" "${multiValueArgs}" ${ARGN} )
  if(NOT PARAM_DIR)
    set(PARAM_DIR ${PARAM_CASENAME})
  endif()
  if(NOT PARAM_PREFIX)
    set(PARAM_PREFIX compareECLFiles)
  endif()
  set(RESULT_PATH ${BASE_RESULT_PATH}${PARAM_DIR_PREFIX}/${PARAM_SIMULATOR}+${PARAM_CASENAME})
  set(TEST_ARGS ${OPM_TESTS_ROOT}/${PARAM_DIR}/${PARAM_FILENAME} ${PARAM_TEST_ARGS})
  opm_add_test(${PARAM_PREFIX}_${PARAM_SIMULATOR}+${PARAM_FILENAME} NO_COMPILE
               EXE_NAME ${PARAM_SIMULATOR}
               DRIVER_ARGS ${OPM_TESTS_ROOT}/${PARAM_DIR} ${RESULT_PATH}
                           ${PROJECT_BINARY_DIR}/bin
                           ${PARAM_FILENAME}
                           ${PARAM_ABS_TOL} ${PARAM_REL_TOL}
                           ${COMPARE_ECL_COMMAND}
               TEST_ARGS ${TEST_ARGS})
endfunction()

###########################################################################
# TEST: add_test_compare_restarted_simulation
###########################################################################

# Input:
#   - casename: basename (no extension)
#
# Details:
#   - This test class compares the output from a restarted simulation
#     to that of a non-restarted simulation.
function(add_test_compare_restarted_simulation)
  set(oneValueArgs CASENAME FILENAME SIMULATOR TEST_NAME ABS_TOL REL_TOL DIR)
  set(multiValueArgs TEST_ARGS)
  cmake_parse_arguments(PARAM "$" "${oneValueArgs}" "${multiValueArgs}" ${ARGN} )
  if(NOT PARAM_DIR)
    set(PARAM_DIR ${PARAM_CASENAME})
  endif()
  if (PARAM_TEST_NAME)
    set(TEST_NAME ${PARAM_TEST_NAME})
  else()
    set(TEST_NAME compareRestartedSim_${PARAM_SIMULATOR}+${PARAM_FILENAME})
  endif()

  set(RESULT_PATH ${BASE_RESULT_PATH}/restart/${PARAM_SIMULATOR}+${PARAM_CASENAME})
  opm_add_test(${TEST_NAME} NO_COMPILE
               EXE_NAME ${PARAM_SIMULATOR}
               DRIVER_ARGS ${OPM_TESTS_ROOT}/${PARAM_DIR} ${RESULT_PATH}
                           ${PROJECT_BINARY_DIR}/bin
                           ${PARAM_FILENAME}
                           ${PARAM_ABS_TOL} ${PARAM_REL_TOL}
                           ${COMPARE_ECL_COMMAND}
                           ${OPM_PACK_COMMAND}
               TEST_ARGS ${PARAM_TEST_ARGS})
endfunction()

###########################################################################
# TEST: add_test_compare_parallel_simulation
###########################################################################

# Input:
#   - casename: basename (no extension)
#
# Details:
#   - This test class compares the output from a parallel simulation
#     to the output from the serial instance of the same model.
function(add_test_compare_parallel_simulation)
  set(oneValueArgs CASENAME FILENAME SIMULATOR ABS_TOL REL_TOL DIR)
  set(multiValueArgs TEST_ARGS)
  cmake_parse_arguments(PARAM "$" "${oneValueArgs}" "${multiValueArgs}" ${ARGN} )

  if(NOT PARAM_DIR)
    set(PARAM_DIR ${PARAM_CASENAME})
  endif()

  set(RESULT_PATH ${BASE_RESULT_PATH}/parallel/${PARAM_SIMULATOR}+${PARAM_CASENAME})
  set(TEST_ARGS ${OPM_TESTS_ROOT}/${PARAM_DIR}/${PARAM_FILENAME} ${PARAM_TEST_ARGS})

  # Add test that runs flow_mpi and outputs the results to file
  opm_add_test(compareParallelSim_${PARAM_SIMULATOR}+${PARAM_FILENAME} NO_COMPILE
               EXE_NAME ${PARAM_SIMULATOR}
               DRIVER_ARGS ${OPM_TESTS_ROOT}/${PARAM_DIR} ${RESULT_PATH}
                           ${PROJECT_BINARY_DIR}/bin
                           ${PARAM_FILENAME}
                           ${PARAM_ABS_TOL} ${PARAM_REL_TOL}
                           ${COMPARE_ECL_COMMAND}
               TEST_ARGS ${TEST_ARGS})
  set_tests_properties(compareParallelSim_${PARAM_SIMULATOR}+${PARAM_FILENAME}
                       PROPERTIES RUN_SERIAL 1)
endfunction()


###########################################################################
# TEST: add_test_compare_parallel_restarted_simulation
###########################################################################

# Input:
#   - casename: basename (no extension)
#
# Details:
#   - This test class compares the output from a restarted parallel simulation
#     to that of a non-restarted parallel simulation.
function(add_test_compare_parallel_restarted_simulation)
  set(oneValueArgs CASENAME FILENAME SIMULATOR ABS_TOL REL_TOL DIR)
  set(multiValueArgs TEST_ARGS)
  cmake_parse_arguments(PARAM "$" "${oneValueArgs}" "${multiValueArgs}" ${ARGN} )
  if(NOT PARAM_DIR)
    set(PARAM_DIR ${PARAM_CASENAME})
  endif()
  set(RESULT_PATH ${BASE_RESULT_PATH}/parallelRestart/${PARAM_SIMULATOR}+${PARAM_CASENAME})

  opm_add_test(compareParallelRestartedSim_${PARAM_SIMULATOR}+${PARAM_FILENAME} NO_COMPILE
               EXE_NAME ${PARAM_SIMULATOR}
               DRIVER_ARGS ${OPM_TESTS_ROOT}/${PARAM_DIR} ${RESULT_PATH}
                           ${PROJECT_BINARY_DIR}/bin
                           ${PARAM_FILENAME}
                           ${PARAM_ABS_TOL} ${PARAM_REL_TOL}
                           ${COMPARE_ECL_COMMAND}
                           ${OPM_PACK_COMMAND}
               TEST_ARGS ${PARAM_TEST_ARGS})
  set_tests_properties(compareParallelRestartedSim_${PARAM_SIMULATOR}+${PARAM_FILENAME}
                       PROPERTIES RUN_SERIAL 1)
endfunction()

if(NOT TARGET test-suite)
  add_custom_target(test-suite)
endif()

# Simple execution tests
opm_set_test_driver(${PROJECT_SOURCE_DIR}/tests/run-test.sh "")
add_test_runSimulator(CASENAME norne
                      FILENAME NORNE_ATW2013
                      SIMULATOR flow
                      PROCS 1)

add_test_runSimulator(CASENAME norne_parallel
                      FILENAME NORNE_ATW2013
                      SIMULATOR flow
                      DIR norne
                      PROCS 4)

# Regression tests
opm_set_test_driver(${PROJECT_SOURCE_DIR}/tests/run-regressionTest.sh "")

# Set absolute tolerance to be used passed to the macros in the following tests
set(abs_tol 2e-2)
set(rel_tol 1e-5)
set(coarse_rel_tol 1e-2)

add_test_compareECLFiles(CASENAME spe12
                         FILENAME SPE1CASE2
                         SIMULATOR flow
                         ABS_TOL ${abs_tol}
                         REL_TOL ${coarse_rel_tol}
                         DIR spe1)

add_test_compareECLFiles(CASENAME spe1_2p
                         FILENAME SPE1CASE2_2P
                         SIMULATOR flow
                         ABS_TOL ${abs_tol}
                         REL_TOL ${rel_tol}
                         DIR spe1)

add_test_compareECLFiles(CASENAME spe1_oilgas
                         FILENAME SPE1CASE2_OILGAS
                         SIMULATOR flow
                         ABS_TOL ${abs_tol}
                         REL_TOL ${coarse_rel_tol}
                         DIR spe1)

add_test_compareECLFiles(CASENAME spe1
                         FILENAME SPE1CASE1
                         SIMULATOR flow
                         ABS_TOL ${abs_tol}
                         REL_TOL ${rel_tol})

add_test_compareECLFiles(CASENAME spe1_import
                         FILENAME SPE1CASE1_IMPORT
                         SIMULATOR flow
                         ABS_TOL ${abs_tol}
                         REL_TOL ${rel_tol}
                         DIR spe1)


add_test_compareECLFiles(CASENAME spe1_nowells
                         FILENAME SPE1CASE2_NOWELLS
                         SIMULATOR flow
                         ABS_TOL ${abs_tol}
                         REL_TOL ${rel_tol}
                         DIR spe1)

add_test_compareECLFiles(CASENAME spe1_thermal
                         FILENAME SPE1CASE2_THERMAL
                         SIMULATOR flow
                         ABS_TOL ${abs_tol}
                         REL_TOL ${rel_tol}
                         DIR spe1)

add_test_compareECLFiles(CASENAME spe1_thermal_watvisc
                         FILENAME SPE1CASE2_THERMAL_WATVISC
                         SIMULATOR flow
                         ABS_TOL ${abs_tol}
                         REL_TOL ${rel_tol}
                         DIR spe1)

add_test_compareECLFiles(CASENAME spe1_rockcomp
                         FILENAME SPE1CASE2_ROCK2DTR
                         SIMULATOR flow
                         ABS_TOL ${abs_tol}
                         REL_TOL ${rel_tol}
                         DIR spe1)

add_test_compareECLFiles(CASENAME spe1_brine
                         FILENAME SPE1CASE1_BRINE
                         SIMULATOR flow
                         ABS_TOL ${abs_tol}
                         REL_TOL ${rel_tol}
                         DIR spe1_brine)

add_test_compareECLFiles(CASENAME spe1_metric_vfp1
                         FILENAME SPE1CASE1_METRIC_VFP1
                         SIMULATOR flow
                         ABS_TOL ${abs_tol}
                         REL_TOL ${rel_tol}
                         DIR vfpprod_spe1)

if(BUILD_FLOW_VARIANTS)
  add_test_compareECLFiles(CASENAME spe1_water
                           FILENAME SPE1CASE1_WATER
                           SIMULATOR flow_onephase
                           ABS_TOL ${abs_tol}
                           REL_TOL ${rel_tol}
                           DIR spe1)

  add_test_compareECLFiles(CASENAME spe1_thermal_onephase
                           FILENAME SPE1CASE2_THERMAL_ONEPHASE
                           SIMULATOR flow_onephase_energy
                           ABS_TOL ${abs_tol}
                           REL_TOL ${rel_tol}
                           DIR spe1)
endif()

add_test_compareECLFiles(CASENAME spe1_spider
                           FILENAME SPIDER_CAKESLICE
                           SIMULATOR flow
                           ABS_TOL ${abs_tol}
                           REL_TOL ${rel_tol}
                           DIR radial_grid)

add_test_compareECLFiles(CASENAME ctaquifer_2d_oilwater
                         FILENAME 2D_OW_CTAQUIFER
                         SIMULATOR flow
                         ABS_TOL ${abs_tol}
                         REL_TOL ${rel_tol}
                         DIR aquifer-oilwater)

add_test_compareECLFiles(CASENAME fetkovich_2d
                         FILENAME 2D_FETKOVICHAQUIFER
                         SIMULATOR flow
                         ABS_TOL ${abs_tol}
                         REL_TOL ${rel_tol}
                         DIR aquifer-fetkovich)

add_test_compareECLFiles(CASENAME numerical_aquifer_3d_2aqu
                         FILENAME 3D_2AQU_NUM
                         SIMULATOR flow
                         ABS_TOL ${abs_tol}
                         REL_TOL ${rel_tol}
                         DIR aquifer-num
                         TEST_ARGS --relaxed-max-pv-fraction=0 --tolerance-cnv=0.00003 --time-step-control=pid --linsolver=cpr)

add_test_compareECLFiles(CASENAME numerical_aquifer_3d_1aqu
                         FILENAME 3D_1AQU_3CELLS
                         SIMULATOR flow
                         ABS_TOL ${abs_tol}
                         REL_TOL ${rel_tol}
                         DIR aquifer-num
                         TEST_ARGS --relaxed-max-pv-fraction=0 --tolerance-cnv=0.00003 --time-step-control=pid --linsolver=cpr)

add_test_compareECLFiles(CASENAME spe3
                         FILENAME SPE3CASE1
                         SIMULATOR flow
                         ABS_TOL ${abs_tol}
                         REL_TOL ${coarse_rel_tol}
                         TEST_ARGS --tolerance-wells=1e-6 --flow-newton-max-iterations=20)

add_test_compareECLFiles(CASENAME spe9
                         FILENAME SPE9_CP_SHORT
                         SIMULATOR flow
                         ABS_TOL ${abs_tol}
                         REL_TOL ${rel_tol})

add_test_compareECLFiles(CASENAME spe9group
                         FILENAME SPE9_CP_GROUP
                         SIMULATOR flow
                         ABS_TOL ${abs_tol}
                         REL_TOL ${rel_tol})

add_test_compareECLFiles(CASENAME spe9group_resv
                         FILENAME SPE9_CP_GROUP_RESV
                         SIMULATOR flow
                         ABS_TOL ${abs_tol}
                         REL_TOL ${rel_tol}
                         DIR spe9group)

add_test_compareECLFiles(CASENAME msw_2d_h
                         FILENAME 2D_H__
                         SIMULATOR flow
                         ABS_TOL ${abs_tol}
                         REL_TOL ${coarse_rel_tol})

add_test_compareECLFiles(CASENAME msw_3d_hfa
                         FILENAME 3D_MSW
                         SIMULATOR flow
                         ABS_TOL ${abs_tol}
                         REL_TOL ${rel_tol}
                         TEST_ARGS --tolerance-pressure-ms-wells=10)

add_test_compareECLFiles(CASENAME polymer_oilwater
                         FILENAME 2D_OILWATER_POLYMER
                         SIMULATOR flow
                         ABS_TOL ${abs_tol}
                         REL_TOL ${rel_tol}
                         TEST_ARGS --tolerance-mb=1.e-7)

add_test_compareECLFiles(CASENAME polymer_injectivity
                         FILENAME 2D_POLYMER_INJECTIVITY
                         SIMULATOR flow
                         ABS_TOL ${abs_tol}
                         REL_TOL ${rel_tol}
                         TEST_ARGS --tolerance-mb=1.e-7 --tolerance-wells=1.e-6)

add_test_compareECLFiles(CASENAME polymer_simple2D
                         FILENAME 2D_THREEPHASE_POLY_HETER
                         SIMULATOR flow
                         ABS_TOL ${abs_tol}
                         REL_TOL ${coarse_rel_tol}
                         TEST_ARGS --tolerance-mb=1.e-7)

add_test_compareECLFiles(CASENAME spe5
                         FILENAME SPE5CASE1
                         SIMULATOR flow
                         ABS_TOL ${abs_tol}
                         REL_TOL ${coarse_rel_tol}
                         TEST_ARGS --flow-newton-max-iterations=20)

add_test_compareECLFiles(CASENAME spe5_co2eor
                         FILENAME SPE5CASE1_DYN
                         SIMULATOR flow
                         ABS_TOL ${abs_tol}
                         REL_TOL ${coarse_rel_tol}
                         TEST_ARGS --flow-newton-max-iterations=20)

add_test_compareECLFiles(CASENAME wecon_wtest
                         FILENAME 3D_WECON
                         SIMULATOR flow
                         ABS_TOL ${abs_tol}
                         REL_TOL ${coarse_rel_tol})

add_test_compareECLFiles(CASENAME msw_model_1
                         FILENAME MSW_MODEL_1
                         SIMULATOR flow
                         ABS_TOL ${abs_tol}
                         REL_TOL ${rel_tol}
                         DIR model1
                         TEST_ARGS --solver-max-time-step-in-days=5.0)

add_test_compareECLFiles(CASENAME base_model_1
                         FILENAME BASE_MODEL_1
                         SIMULATOR flow
                         ABS_TOL ${abs_tol}
                         REL_TOL ${rel_tol}
                         DIR model1)

add_test_compareECLFiles(CASENAME faults_model_1
                         FILENAME FAULTS_MODEL_1
                         SIMULATOR flow
                         ABS_TOL ${abs_tol}
                         REL_TOL ${rel_tol}
                         DIR model1
                         TEST_ARGS --solver-max-time-step-in-days=5.0)

add_test_compareECLFiles(CASENAME base_model2
                         FILENAME 0_BASE_MODEL2
                         SIMULATOR flow
                         ABS_TOL ${abs_tol}
                         REL_TOL ${rel_tol}
                         DIR model2)

add_test_compareECLFiles(CASENAME 0a1_grpctl_stw_model2
                         FILENAME 0A1_GRCTRL_LRAT_ORAT_BASE_MODEL2_STW
                         SIMULATOR flow
                         ABS_TOL ${abs_tol}
                         REL_TOL ${rel_tol}
                         DIR model2)

add_test_compareECLFiles(CASENAME 0a1_grpctl_msw_model2
                         FILENAME 0A1_GRCTRL_LRAT_ORAT_BASE_MODEL2_MSW
                         SIMULATOR flow
                         ABS_TOL ${abs_tol}
                         REL_TOL ${rel_tol}
                         DIR model2
                         TEST_ARGS --solver-max-time-step-in-days=3)

add_test_compareECLFiles(CASENAME 0a2_grpctl_stw_model2
                         FILENAME 0A2_GRCTRL_LRAT_ORAT_GGR_BASE_MODEL2_STW
                         SIMULATOR flow
                         ABS_TOL ${abs_tol}
                         REL_TOL ${rel_tol}
                         DIR model2)

add_test_compareECLFiles(CASENAME 0a2_grpctl_msw_model2
                         FILENAME 0A2_GRCTRL_LRAT_ORAT_GGR_BASE_MODEL2_MSW
                         SIMULATOR flow
                         ABS_TOL ${abs_tol}
                         REL_TOL ${rel_tol}
                         DIR model2
                         TEST_ARGS --solver-max-time-step-in-days=3)

add_test_compareECLFiles(CASENAME 0a3_grpctl_stw_model2
                         FILENAME 0A3_GRCTRL_LRAT_LRAT_BASE_MODEL2_STW
                         SIMULATOR flow
                         ABS_TOL ${abs_tol}
                         REL_TOL ${rel_tol}
                         DIR model2)

add_test_compareECLFiles(CASENAME 0a3_grpctl_msw_model2
                         FILENAME 0A3_GRCTRL_LRAT_LRAT_BASE_MODEL2_MSW
                         SIMULATOR flow
                         ABS_TOL ${abs_tol}
                         REL_TOL ${rel_tol}
                         DIR model2
                         TEST_ARGS --solver-max-time-step-in-days=3)

add_test_compareECLFiles(CASENAME 0a4_grpctl_stw_model2
                         FILENAME 0A4_GRCTRL_LRAT_LRAT_GGR_BASE_MODEL2_STW
                         SIMULATOR flow
                         ABS_TOL ${abs_tol}
                         REL_TOL ${rel_tol}
                         DIR model2)

add_test_compareECLFiles(CASENAME 0a4_grpctl_msw_model2
                         FILENAME 0A4_GRCTRL_LRAT_LRAT_GGR_BASE_MODEL2_MSW
                         SIMULATOR flow
                         ABS_TOL ${abs_tol}
                         REL_TOL ${rel_tol}
                         DIR model2
                         TEST_ARGS --solver-max-time-step-in-days=3)

add_test_compareECLFiles(CASENAME multregt_model2
                         FILENAME 1_MULTREGT_MODEL2
                         SIMULATOR flow
                         ABS_TOL ${abs_tol}
                         REL_TOL ${rel_tol}
                         DIR model2)

add_test_compareECLFiles(CASENAME udq_actionx
                         FILENAME UDQ_ACTIONX
                         SIMULATOR flow
                         ABS_TOL ${abs_tol}
                         REL_TOL ${rel_tol}
                         DIR udq_actionx)

add_test_compareECLFiles(CASENAME udq_wconprod
                         FILENAME UDQ_WCONPROD
                         SIMULATOR flow
                         ABS_TOL ${abs_tol}
                         REL_TOL ${rel_tol}
                         DIR udq_actionx)

add_test_compareECLFiles(CASENAME actionx_m1
                         FILENAME ACTIONX_M1
                         SIMULATOR flow
                         ABS_TOL ${abs_tol}
                         REL_TOL ${rel_tol}
                         DIR udq_actionx
                         TEST_ARGS --solver-max-time-step-in-days=1)

add_test_compareECLFiles(CASENAME pinch_multz_all
                         FILENAME PINCH_MULTZ_ALL
                         SIMULATOR flow
                         ABS_TOL ${abs_tol}
                         REL_TOL ${rel_tol}
                         DIR pinch)

add_test_compareECLFiles(CASENAME pinch_multz_all_barrier
                         FILENAME PINCH_MULTZ_ALL_BARRIER
                         SIMULATOR flow
                         ABS_TOL ${abs_tol}
                         REL_TOL ${rel_tol}
                         DIR pinch)

add_test_compareECLFiles(CASENAME udq_uadd
                         FILENAME UDQ_M1
                         SIMULATOR flow
                         ABS_TOL ${abs_tol}
                         REL_TOL ${rel_tol}
                         DIR udq_actionx)

add_test_compareECLFiles(CASENAME udq_undefined
                         FILENAME UDQ_M2
                         SIMULATOR flow
                         ABS_TOL ${abs_tol}
                         REL_TOL ${rel_tol}
                         DIR udq_actionx)

add_test_compareECLFiles(CASENAME udq_in_actionx
                         FILENAME UDQ_M3
                         SIMULATOR flow
                         ABS_TOL ${abs_tol}
                         REL_TOL ${rel_tol}
                         DIR udq_actionx)

add_test_compareECLFiles(CASENAME co2store
                         FILENAME CO2STORE
                         SIMULATOR flow
                         ABS_TOL ${abs_tol}
                         REL_TOL ${rel_tol}
                         DIR co2store)

add_test_compareECLFiles(CASENAME co2store_diffusive
                         FILENAME CO2STORE_DIFFUSIVE
                         SIMULATOR flow
                         ABS_TOL ${abs_tol}
                         REL_TOL ${rel_tol}
                         DIR co2store)

add_test_compareECLFiles(CASENAME co2store_drsdtcon
                         FILENAME CO2STORE_DRSDTCON
                         SIMULATOR flow
                         ABS_TOL ${abs_tol}
                         REL_TOL ${rel_tol}
                         DIR co2store)

if (opm-common_EMBEDDED_PYTHON)
  add_test_compareECLFiles(CASENAME udq_pyaction
                           FILENAME PYACTION_WCONPROD
                           SIMULATOR flow
                           ABS_TOL ${abs_tol}
                           REL_TOL ${rel_tol}
                           DIR udq_actionx)
endif()

add_test_compareECLFiles(CASENAME multxyz_model2
			  FILENAME 2_MULTXYZ_MODEL2
			  SIMULATOR flow
			  ABS_TOL ${abs_tol}
			  REL_TOL ${rel_tol}
			  DIR model2)

add_test_compareECLFiles(CASENAME multflt_model2
			  FILENAME 3_MULTFLT_MODEL2
			  SIMULATOR flow
			  ABS_TOL ${abs_tol}
			  REL_TOL ${rel_tol}
			  DIR model2)

add_test_compareECLFiles(CASENAME multpvv_model2
			  FILENAME 4_MINPVV_MODEL2
			  SIMULATOR flow
			  ABS_TOL ${abs_tol}
			  REL_TOL ${rel_tol}
			  DIR model2)

add_test_compareECLFiles(CASENAME swatinit_model2
			  FILENAME 5_SWATINIT_MODEL2
			  SIMULATOR flow
			  ABS_TOL ${abs_tol}
			  REL_TOL ${rel_tol}
			  DIR model2)

add_test_compareECLFiles(CASENAME endscale_model2
			  FILENAME 6_ENDSCALE_MODEL2
			  SIMULATOR flow
			  ABS_TOL ${abs_tol}
			  REL_TOL ${rel_tol}
			  DIR model2)

add_test_compareECLFiles(CASENAME hysteresis_model2
			  FILENAME 7_HYSTERESIS_MODEL2
			  SIMULATOR flow
			  ABS_TOL ${abs_tol}
			  REL_TOL ${rel_tol}
			  DIR model2)

add_test_compareECLFiles(CASENAME multiply_tranxyz_model2
			  FILENAME 8_MULTIPLY_TRANXYZ_MODEL2
			  SIMULATOR flow
			  ABS_TOL ${abs_tol}
			  REL_TOL ${rel_tol}
			  DIR model2)

add_test_compareECLFiles(CASENAME editnnc_model2
			  FILENAME 9_EDITNNC_MODEL2
			  SIMULATOR flow
			  ABS_TOL ${abs_tol}
			  REL_TOL ${rel_tol}
			  DIR model2)

add_test_compareECLFiles(CASENAME 9_1a_grpctl_stw_model2
                         FILENAME 9_1A_DEPL_MAX_RATE_MIN_BHP_STW
                         SIMULATOR flow
                         ABS_TOL ${abs_tol}
                         REL_TOL ${rel_tol}
                         DIR model2)

add_test_compareECLFiles(CASENAME 9_1a_grpctl_msw_model2
                         FILENAME 9_1A_DEPL_MAX_RATE_MIN_BHP_MSW
                         SIMULATOR flow
                         ABS_TOL ${abs_tol}
                         REL_TOL ${rel_tol}
                         DIR model2)

add_test_compareECLFiles(CASENAME 9_1b_grpctl_stw_model2
                         FILENAME 9_1B_DEPL_MAX_RATE_MIN_THP_STW
                         SIMULATOR flow
                         ABS_TOL ${abs_tol}
                         REL_TOL ${rel_tol}
                         DIR model2)

add_test_compareECLFiles(CASENAME 9_1b_grpctl_msw_model2
                         FILENAME 9_1B_DEPL_MAX_RATE_MIN_THP_MSW
                         SIMULATOR flow
                         ABS_TOL ${abs_tol}
                         REL_TOL ${rel_tol}
                         DIR model2)

add_test_compareECLFiles(CASENAME 9_2a_grpctl_stw_model2
                         FILENAME 9_2A_DEPL_GCONPROD_1L_STW
                         SIMULATOR flow
                         ABS_TOL ${abs_tol}
                         REL_TOL ${rel_tol}
                         DIR model2)

add_test_compareECLFiles(CASENAME 9_2a_grpctl_msw_model2
                         FILENAME 9_2A_DEPL_GCONPROD_1L_MSW
                         SIMULATOR flow
                         ABS_TOL ${abs_tol}
                         REL_TOL ${rel_tol}
                         DIR model2)

add_test_compareECLFiles(CASENAME 9_2b_grpctl_stw_model2
                         FILENAME 9_2B_DEPL_GCONPROD_2L_STW
                         SIMULATOR flow
                         ABS_TOL ${abs_tol}
                         REL_TOL ${rel_tol}
                         DIR model2)

add_test_compareECLFiles(CASENAME 9_2b_grpctl_msw_model2
                         FILENAME 9_2B_DEPL_GCONPROD_2L_MSW
                         SIMULATOR flow
                         ABS_TOL ${abs_tol}
                         REL_TOL ${rel_tol}
                         DIR model2)

add_test_compareECLFiles(CASENAME 9_3a_grpctl_stw_model2
                         FILENAME 9_3A_GINJ_REIN-G_STW
                         SIMULATOR flow
                         ABS_TOL ${abs_tol}
                         REL_TOL ${rel_tol}
                         DIR model2)

add_test_compareECLFiles(CASENAME 9_3a_grpctl_msw_model2
                         FILENAME 9_3A_GINJ_REIN-G_MSW
                         SIMULATOR flow
                         ABS_TOL ${abs_tol}
                         REL_TOL ${rel_tol}
                         DIR model2)

add_test_compareECLFiles(CASENAME 9_3b_grpctl_stw_model2
                         FILENAME 9_3B_GINJ_GAS_EXPORT_STW
                         SIMULATOR flow
                         ABS_TOL ${abs_tol}
                         REL_TOL ${rel_tol}
                         DIR model2)

add_test_compareECLFiles(CASENAME 9_3b_grpctl_msw_model2
                         FILENAME 9_3B_GINJ_GAS_EXPORT_MSW
                         SIMULATOR flow
                         ABS_TOL ${abs_tol}
                         REL_TOL ${rel_tol}
                         DIR model2)

add_test_compareECLFiles(CASENAME 9_3c_grpctl_stw_model2
                         FILENAME 9_3C_GINJ_GAS_GCONSUMP_STW
                         SIMULATOR flow
                         ABS_TOL ${abs_tol}
                         REL_TOL ${rel_tol}
                         DIR model2)

add_test_compareECLFiles(CASENAME 9_3c_grpctl_msw_model2
                         FILENAME 9_3C_GINJ_GAS_GCONSUMP_MSW
                         SIMULATOR flow
                         ABS_TOL ${abs_tol}
                         REL_TOL ${rel_tol}
                         DIR model2)

add_test_compareECLFiles(CASENAME 9_3d_grpctl_stw_model2
                         FILENAME 9_3D_GINJ_GAS_MAX_EXPORT_STW
                         SIMULATOR flow
                         ABS_TOL ${abs_tol}
                         REL_TOL ${rel_tol}
                         DIR model2)

add_test_compareECLFiles(CASENAME 9_3d_grpctl_msw_model2
                         FILENAME 9_3D_GINJ_GAS_MAX_EXPORT_MSW
                         SIMULATOR flow
                         ABS_TOL ${abs_tol}
                         REL_TOL ${rel_tol}
                         DIR model2)

add_test_compareECLFiles(CASENAME 9_3e_grpctl_stw_model2
                         FILENAME 9_3E_GAS_MIN_EXPORT_STW
                         SIMULATOR flow
                         ABS_TOL ${abs_tol}
                         REL_TOL ${rel_tol}
                         DIR model2)

add_test_compareECLFiles(CASENAME 9_3e_grpctl_msw_model2
                         FILENAME 9_3E_GAS_MIN_EXPORT_MSW
                         SIMULATOR flow
                         ABS_TOL ${abs_tol}
                         REL_TOL ${rel_tol}
                         DIR model2)

add_test_compareECLFiles(CASENAME 9_4a_grpctl_stw_model2
                         FILENAME 9_4A_WINJ_MAXWRATES_MAXBHP_GCONPROD_1L_STW
                         SIMULATOR flow
                         ABS_TOL ${abs_tol}
                         REL_TOL ${rel_tol}
                         DIR model2)

add_test_compareECLFiles(CASENAME 9_4a_grpctl_msw_model2
                         FILENAME 9_4A_WINJ_MAXWRATES_MAXBHP_GCONPROD_1L_MSW
                         SIMULATOR flow
                         ABS_TOL ${abs_tol}
                         REL_TOL ${rel_tol}
                         DIR model2)

add_test_compareECLFiles(CASENAME 9_4b_grpctl_stw_model2
                         FILENAME 9_4B_WINJ_VREP-W_STW
                         SIMULATOR flow
                         ABS_TOL ${abs_tol}
                         REL_TOL ${rel_tol}
                         DIR model2)

add_test_compareECLFiles(CASENAME 9_4b_grpctl_msw_model2
                         FILENAME 9_4B_WINJ_VREP-W_MSW
                         SIMULATOR flow
                         ABS_TOL ${abs_tol}
                         REL_TOL ${rel_tol}
                         DIR model2)

add_test_compareECLFiles(CASENAME 9_4c_grpctl_stw_model2
                         FILENAME 9_4C_WINJ_GINJ_VREP-W_REIN-G_STW
                         SIMULATOR flow
                         ABS_TOL ${abs_tol}
                         REL_TOL ${rel_tol}
                         DIR model2)

add_test_compareECLFiles(CASENAME 9_4c_grpctl_msw_model2
                         FILENAME 9_4C_WINJ_GINJ_VREP-W_REIN-G_MSW
                         SIMULATOR flow
                         ABS_TOL ${abs_tol}
                         REL_TOL ${rel_tol}
                         DIR model2)

add_test_compareECLFiles(CASENAME 9_4d_grpctl_stw_model2
                         FILENAME 9_4D_WINJ_GINJ_GAS_EXPORT_STW
                         SIMULATOR flow
                         ABS_TOL ${abs_tol}
                         REL_TOL ${rel_tol}
                         DIR model2)

add_test_compareECLFiles(CASENAME 9_4d_grpctl_msw_model2
                         FILENAME 9_4D_WINJ_GINJ_GAS_EXPORT_MSW
                         SIMULATOR flow
                         ABS_TOL ${abs_tol}
                         REL_TOL ${rel_tol}
                         DIR model2)

add_test_compareECLFiles(CASENAME model4_group
                         FILENAME MOD4_GRP
                         SIMULATOR flow
                         ABS_TOL ${abs_tol}
                         REL_TOL ${rel_tol}
                         DIR model4)

add_test_compareECLFiles(CASENAME model4_udq_group
                         FILENAME MOD4_UDQ_ACTIONX
                         SIMULATOR flow
                         ABS_TOL ${abs_tol}
                         REL_TOL ${rel_tol}
                         DIR model4)

add_test_compareECLFiles(CASENAME model6_msw
                         FILENAME 1_MSW_MODEL6
                         SIMULATOR flow
                         ABS_TOL ${abs_tol}
                         REL_TOL ${rel_tol}
                         DIR model6)

add_test_compareECLFiles(CASENAME wsegsicd
			  FILENAME TEST_WSEGSICD
			  SIMULATOR flow
			  ABS_TOL ${abs_tol}
			  REL_TOL ${rel_tol})

add_test_compareECLFiles(CASENAME wsegaicd
			  FILENAME BASE_MSW_WSEGAICD
			  SIMULATOR flow
			  ABS_TOL ${abs_tol}
			  REL_TOL ${rel_tol})

add_test_compareECLFiles(CASENAME nnc
                         FILENAME NNC_AND_EDITNNC
                         SIMULATOR flow
                         ABS_TOL ${abs_tol}
                         REL_TOL ${rel_tol}
                         DIR editnnc)

add_test_compareECLFiles(CASENAME spe1_foam
                         FILENAME SPE1FOAM
                         SIMULATOR flow
                         ABS_TOL ${abs_tol}
                         REL_TOL ${rel_tol}
                         DIR spe1_foam)

add_test_compareECLFiles(CASENAME bc_lab
                         FILENAME BC_LAB
                         SIMULATOR flow
                         ABS_TOL ${abs_tol}
                         REL_TOL ${rel_tol}
                         DIR bc_lab)

add_test_compareECLFiles(CASENAME norne_reperf
                         FILENAME NORNE_ATW2013_B1H_RE-PERF
                         SIMULATOR flow
                         ABS_TOL ${abs_tol}
                         REL_TOL ${rel_tol}
                         DIR norne)

add_test_compareECLFiles(CASENAME compl_smry
                         FILENAME COMPL_SMRY
                         SIMULATOR flow
                         ABS_TOL ${abs_tol}
                         REL_TOL ${rel_tol}
                         DIR compl_smry)

add_test_compareECLFiles(CASENAME 3d_tran_operator
                         FILENAME 3D_TRAN_OPERATOR
                         SIMULATOR flow
                         ABS_TOL ${abs_tol}
                         REL_TOL ${rel_tol}
                         DIR parallel_fieldprops)

# Restart tests
opm_set_test_driver(${PROJECT_SOURCE_DIR}/tests/run-restart-regressionTest.sh "")

# Cruder tolerances for the restarted tests
set(abs_tol_restart 2e-1)
set(rel_tol_restart 4e-4)

add_test_compare_restarted_simulation(CASENAME spe1
                                      FILENAME SPE1CASE2_ACTNUM
                                      SIMULATOR flow
                                      ABS_TOL ${abs_tol_restart}
                                      REL_TOL ${rel_tol_restart}
                                      TEST_ARGS --sched-restart=false)
add_test_compare_restarted_simulation(CASENAME spe9
                                      FILENAME SPE9_CP_SHORT
                                      SIMULATOR flow
                                      ABS_TOL ${abs_tol_restart}
                                      REL_TOL ${rel_tol_restart}
                                      TEST_ARGS --sched-restart=false)

add_test_compare_restarted_simulation(CASENAME udq_actionx
                                      FILENAME UDQ_M1
                                      SIMULATOR flow
                                      ABS_TOL ${abs_tol_restart}
                                      REL_TOL ${rel_tol_restart}
                                      TEST_ARGS --sched-restart=false)

# The dynamic MSW data is not written to /read from the restart file
# We therefore accept significant deviation in the results.
# Note also that we use --sched-restart=true since some necessary
# MSW info is still lacking in the restart file.
set(abs_tol_restart_msw 2e2)
set(rel_tol_restart_msw 1e-3)

add_test_compare_restarted_simulation(CASENAME msw_3d_hfa
                                      FILENAME 3D_MSW
                                      SIMULATOR flow
                                      ABS_TOL ${abs_tol_restart_msw}
                                      REL_TOL ${rel_tol_restart_msw}
                                      TEST_ARGS --enable-adaptive-time-stepping=false --sched-restart=true)


# Basic restart tests which only compare the summary output, this test driver should
# only be used in situations where it is challenging to get agreement in the restart file.
opm_set_test_driver(${PROJECT_SOURCE_DIR}/tests/run-summary-restart-regressionTest.sh "")

add_test_compare_restarted_simulation(CASENAME spe1
                                      FILENAME SPE1CASE2_ACTNUM
                                      SIMULATOR flow
                                      TEST_NAME restart_spe1_summary
                                      ABS_TOL ${abs_tol_restart}
                                      REL_TOL ${rel_tol_restart}
                                      TEST_ARGS --sched-restart=false)

add_test_compare_restarted_simulation(CASENAME udq_actionx
                                      FILENAME UDQ_M1
                                      SIMULATOR flow
                                      TEST_NAME restart_udq_m1_summary
                                      ABS_TOL ${abs_tol_restart}
                                      REL_TOL ${rel_tol_restart}
                                      TEST_ARGS --sched-restart=false --solver-max-time-step-in-days=1.0)

# PORV test
opm_set_test_driver(${PROJECT_SOURCE_DIR}/tests/run-porv-acceptanceTest.sh "")
add_test_compareECLFiles(CASENAME norne
                         FILENAME NORNE_ATW2013
                         SIMULATOR flow
                         ABS_TOL 1e-5
                         REL_TOL 1e-8
                         PREFIX comparePORV
                         DIR_PREFIX /porv)

# Init tests
opm_set_test_driver(${PROJECT_SOURCE_DIR}/tests/run-init-regressionTest.sh "")

add_test_compareECLFiles(CASENAME norne
                         FILENAME NORNE_ATW2013
                         SIMULATOR flow
                         ABS_TOL ${abs_tol}
                         REL_TOL ${rel_tol}
                         PREFIX compareECLInitFiles
                         DIR_PREFIX /init)

# This is not a proper regression test; the test will load a norne case prepared
# for restart and run one single timestep - of length one day. The results are not
# verified in any way.
add_test(NAME NORNE_RESTART
         COMMAND flow --output-dir=${BASE_RESULT_PATH}/norne-restart ${OPM_TESTS_ROOT}/norne/NORNE_ATW2013_RESTART.DATA)


# Parallel tests
if(MPI_FOUND)
  opm_set_test_driver(${PROJECT_SOURCE_DIR}/tests/run-parallel-restart-regressionTest.sh "")
  add_test_compare_parallel_restarted_simulation(CASENAME spe1
                                                 FILENAME SPE1CASE2_ACTNUM
                                                 SIMULATOR flow
                                                 ABS_TOL ${abs_tol_restart}
                                                 REL_TOL ${rel_tol_restart}
                                                 TEST_ARGS --sched-restart=false)


  opm_set_test_driver(${PROJECT_SOURCE_DIR}/tests/run-parallel-regressionTest.sh "")

  # Different tolerances for these tests
  set(abs_tol_parallel 0.02)
  set(rel_tol_parallel 8e-5)
  set(coarse_rel_tol_parallel 1e-2)

  add_test_compare_parallel_simulation(CASENAME spe1
                                       FILENAME SPE1CASE2
                                       SIMULATOR flow
                                       ABS_TOL ${abs_tol_parallel}
                                       REL_TOL ${rel_tol_parallel}
                                       TEST_ARGS --linear-solver-reduction=1e-7 --tolerance-cnv=5e-6 --tolerance-mb=1e-8)

  add_test_compare_parallel_simulation(CASENAME spe9
                                       FILENAME SPE9_CP_SHORT
                                       SIMULATOR flow
                                       ABS_TOL ${abs_tol_parallel}
                                       REL_TOL ${rel_tol_parallel}
                                       TEST_ARGS --linear-solver-reduction=1e-7 --tolerance-cnv=5e-6 --tolerance-mb=1e-8)

  # A test for distributed standard wells. We load distribute only along the z-axis
  add_test_compare_parallel_simulation(CASENAME spe9
                                       FILENAME SPE9_CP_SHORT
                                       SIMULATOR flow_distribute_z
                                       ABS_TOL ${abs_tol_parallel}
                                       REL_TOL ${rel_tol_parallel}
                                       TEST_ARGS --linear-solver-reduction=1e-7 --tolerance-cnv=5e-6 --tolerance-mb=1e-8)

  add_test_compare_parallel_simulation(CASENAME spe9group
                                       FILENAME SPE9_CP_GROUP
                                       SIMULATOR flow
                                       ABS_TOL ${abs_tol_parallel}
                                       REL_TOL ${coarse_rel_tol_parallel}
                                       TEST_ARGS --linear-solver-reduction=1e-7 --tolerance-cnv=5e-6 --tolerance-mb=1e-8)

  add_test_compare_parallel_simulation(CASENAME spe3
                                       FILENAME SPE3CASE1
                                       SIMULATOR flow
                                       ABS_TOL ${abs_tol_parallel}
                                       REL_TOL ${coarse_rel_tol_parallel}
                                       TEST_ARGS --linear-solver-reduction=1e-7 --tolerance-cnv=5e-6 --tolerance-mb=1e-8)

  add_test_compare_parallel_simulation(CASENAME spe1_solvent
                                       FILENAME SPE1CASE2_SOLVENT
                                       SIMULATOR flow
                                       ABS_TOL ${abs_tol_parallel}
                                       REL_TOL ${coarse_rel_tol_parallel}
                                       TEST_ARGS --linear-solver-reduction=1e-7 --tolerance-cnv=5e-6 --tolerance-mb=1e-8)

  add_test_compare_parallel_simulation(CASENAME polymer_simple2D
                                       FILENAME 2D_THREEPHASE_POLY_HETER
                                       SIMULATOR flow
                                       ABS_TOL ${abs_tol}
                                       REL_TOL ${coarse_rel_tol}
                                       TEST_ARGS --linear-solver-reduction=1e-7 --tolerance-cnv=5e-6 --tolerance-mb=1e-8)

  add_test_compare_parallel_simulation(CASENAME spe1_foam
                                       FILENAME SPE1FOAM
                                       SIMULATOR flow
                                       ABS_TOL ${abs_tol}
                                       REL_TOL ${coarse_rel_tol_parallel})

  add_test_compare_parallel_simulation(CASENAME spe1_thermal
                                       FILENAME SPE1CASE2_THERMAL
                                       SIMULATOR flow
                                       ABS_TOL ${abs_tol}
                                       REL_TOL ${rel_tol}
                                       DIR spe1
                                       TEST_ARGS --linear-solver-reduction=1e-7 --tolerance-cnv=5e-6 --tolerance-mb=1e-8)

if(BUILD_FLOW_VARIANTS)
  add_test_compare_parallel_simulation(CASENAME spe1_thermal
                                       FILENAME SPE1CASE2_THERMAL_ONEPHASE
                                       SIMULATOR flow_onephase_energy
                                       ABS_TOL ${abs_tol}
                                       REL_TOL ${rel_tol}
                                       DIR spe1
                                       TEST_ARGS --linear-solver-reduction=1e-7 --tolerance-cnv=5e-6 --tolerance-mb=1e-8)

  add_test_compare_parallel_simulation(CASENAME spe1_water
                                       FILENAME SPE1CASE1_WATER
                                       SIMULATOR flow_onephase
                                       ABS_TOL ${abs_tol}
                                       REL_TOL ${rel_tol}
                                       DIR spe1
                                       TEST_ARGS --linear-solver-reduction=1e-7 --tolerance-cnv=5e-6 --tolerance-mb=1e-8)
endif()

  add_test_compare_parallel_simulation(CASENAME spe1_brine
                                       FILENAME SPE1CASE1_BRINE
                                       SIMULATOR flow
                                       ABS_TOL ${abs_tol_parallel}
                                       REL_TOL ${coarse_rel_tol_parallel}
                                       TEST_ARGS --linear-solver-reduction=1e-7 --tolerance-cnv=5e-6 --tolerance-mb=1e-6)

  add_test_compare_parallel_simulation(CASENAME fetkovich_2d
                                       FILENAME 2D_FETKOVICHAQUIFER
                                       SIMULATOR flow
                                       ABS_TOL ${abs_tol_parallel}
                                       REL_TOL ${rel_tol_parallel}
                                       DIR aquifer-fetkovich
                                       TEST_ARGS --linear-solver-reduction=1e-7 --tolerance-cnv=5e-6 --tolerance-mb=1e-6)

  add_test_compare_parallel_simulation(CASENAME ctaquifer_2d_oilwater
                                       FILENAME 2D_OW_CTAQUIFER
                                       SIMULATOR flow
                                       ABS_TOL ${abs_tol_parallel}
                                       REL_TOL ${rel_tol_parallel}
                                       DIR aquifer-oilwater
                                       TEST_ARGS --linear-solver-reduction=1e-7 --tolerance-cnv=5e-6 --tolerance-mb=1e-6)

  add_test_compare_parallel_simulation(CASENAME 3d_tran_operator
                                       FILENAME 3D_TRAN_OPERATOR
                                       SIMULATOR flow
                                       ABS_TOL ${abs_tol_parallel}
                                       REL_TOL ${rel_tol_parallel}
                                       DIR parallel_fieldprops
                                       TEST_ARGS --linear-solver-reduction=1e-7 --tolerance-cnv=5e-6 --tolerance-mb=1e-6)

  add_test_compare_parallel_simulation(CASENAME numerical_aquifer_3d_2aqu
                                       FILENAME 3D_2AQU_NUM
                                       SIMULATOR flow
                                       ABS_TOL 0.12
                                       REL_TOL ${coarse_rel_tol_parallel}
                                       DIR aquifer-num
                                       TEST_ARGS --relaxed-max-pv-fraction=0 --tolerance-cnv=0.00003 --time-step-control=pid --linsolver=cpr)

  add_test_compare_parallel_simulation(CASENAME numerical_aquifer_3d_1aqu
                                       FILENAME 3D_1AQU_3CELLS
                                       SIMULATOR flow
                                       ABS_TOL ${abs_tol_parallel}
                                       REL_TOL ${coarse_rel_tol_parallel}
                                       DIR aquifer-num
                                       TEST_ARGS --relaxed-max-pv-fraction=0 --tolerance-cnv=0.00003 --time-step-control=pid --linsolver=cpr)

  add_test_compare_parallel_simulation(CASENAME actionx_m1
                                       FILENAME ACTIONX_M1
                                       SIMULATOR flow
                                       ABS_TOL ${abs_tol_parallel}
                                       REL_TOL ${coarse_rel_tol_parallel}
                                       DIR udq_actionx
                                       TEST_ARGS --linear-solver-reduction=1e-7 --tolerance-cnv=5e-6 --tolerance-mb=1e-6)
endif()
