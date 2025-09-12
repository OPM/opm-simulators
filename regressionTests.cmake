# Regression tests
opm_set_test_driver(${PROJECT_SOURCE_DIR}/tests/run-regressionTest.sh "")

# Set absolute tolerance to be used passed to the macros in the following tests
set(abs_tol 2e-2)
set(rel_tol 1e-5)
set(coarse_rel_tol 1e-2)

# Adds several tests cases with similar parameters
# cases Variable name of list with test cases
# prefix Prefix to use
# argn Parameters for cases
macro(add_multiple_tests cases prefix)
  foreach(case ${${cases}})
    string(TOLOWER ${case} test)
    add_test_compareECLFiles(
        CASENAME ${prefix}${test}
        FILENAME ${case}
        ${ARGN}
    )
  endforeach()
endmacro()

# Adds several tests cases in a numerical range with similar parameters
# start Start of range
# end End fof range
# ftemplate File name template to use
# prefix Prefix to use
# argn Parameters for cases
macro(add_multiple_test_range start end ftemplate prefix)
    foreach(case RANGE ${start} ${end})
      add_test_compareECLFiles(
          CASENAME ${prefix}_${case}
          FILENAME ${ftemplate}${case}
          ${ARGN}
      )
    endforeach()
endmacro()

add_test_compareECLFiles(CASENAME spe1flowexp
                         FILENAME SPE1CASE2
                         SIMULATOR flowexp_blackoil
                         ABS_TOL ${abs_tol}
                         REL_TOL ${rel_tol}
                         DIR spe1)

add_test_compareECLFiles(CASENAME 1dcompositional
                         FILENAME 1D_COMP
                         SIMULATOR flowexp_comp
                         ABS_TOL ${abs_tol}
                         REL_TOL ${rel_tol}
                         DIR compositional)

add_test_compareECLFiles(CASENAME spe12
                         FILENAME SPE1CASE2
                         SIMULATOR flow
                         ABS_TOL ${abs_tol}
                         REL_TOL ${coarse_rel_tol}
                         RESTART_SCHED false
                         RESTART_STEP 60
                         DIR spe1)

set(_spe1_tests
  SPE1CASE1
  SPE1CASE1_IMPORT
  SPE1CASE1_WATER
  SPE1CASE2_KRNUM
  SPE1CASE2_NOWELLS
  SPE1CASE2_ROCK2DTR
  SPE1CASE2_THERMAL
  SPE1CASE2_THERMAL_ONEPHASE
  SPE1CASE2_THERMAL_WATVISC
  SPE1CASE2_2P
)

add_multiple_tests(
  _spe1_tests
  ""
  SIMULATOR flow
  ABS_TOL ${abs_tol}
  REL_TOL ${rel_tol}
  DIR spe1
)

set(_spe1_coarse_tests
  SPE1CASE2_GASWATER
  SPE1CASE2_GASWATER_MSW
  SPE1CASE2_OILGAS
)

add_multiple_tests(
  _spe1_coarse_tests
  ""
  SIMULATOR flow
  ABS_TOL ${abs_tol}
  REL_TOL ${coarse_rel_tol}
  DIR spe1
)

set(_spe1_brine_tests
  SPE1CASE1_BRINE
  SPE1CASE2_BRINE_GASWATER
)

add_multiple_tests(
  _spe1_brine_tests
  ""
  SIMULATOR flow
  ABS_TOL ${abs_tol}
  REL_TOL ${rel_tol}
  DIR spe1_brine
)

set(_spe1_precsalt_tests
  GASWATER_VAPWAT_PRECSALT
  SPE1CASE1_PRECSALT
)

add_multiple_tests(
  _spe1_precsalt_tests
  ""
  SIMULATOR flow
  ABS_TOL ${abs_tol}
  REL_TOL ${rel_tol}
  DIR spe1_precsalt
)

add_test_compareECLFiles(CASENAME gasoil_precsalt
                         FILENAME GASCONDENSATE_VAPWAT_PRECSALT_REGRESSION
                         SIMULATOR flow
                         ABS_TOL ${abs_tol}
                         REL_TOL ${rel_tol}
                         DIR spe1_precsalt
                         TEST_ARGS --solver-max-time-step-in-days=0.05)

set(_network_tuning_tests
  NETWORK-01
  NETWORK-01_STANDARD
)

add_multiple_tests(
  _network_tuning_tests
  ""
  SIMULATOR flow
  ABS_TOL ${abs_tol}
  REL_TOL ${rel_tol}
  DIR network
  TEST_ARGS --enable-tuning=true
)

set(_network_tuning_local_switch_tests
  NETWORK-01-REROUTE
  NETWORK-01-REROUTE_STD
)

add_multiple_tests(
  _network_tuning_local_switch_tests
  ""
  SIMULATOR flow
  ABS_TOL ${abs_tol}
  REL_TOL ${rel_tol}
  DIR network
  TEST_ARGS --enable-tuning=true --local-well-solve-control-switching=true
)

add_test_compareECLFiles(CASENAME network_01_wtest
                         FILENAME NETWORK-01-WTEST
                         SIMULATOR flow
                         ABS_TOL ${abs_tol}
                         REL_TOL ${rel_tol}
                         DIR network
                         TEST_ARGS --enable-tuning=true)

add_test_compareECLFiles(CASENAME spe1_metric_vfp1
                         FILENAME SPE1CASE1_METRIC_VFP1
                         SIMULATOR flow
                         ABS_TOL ${abs_tol}
                         REL_TOL ${rel_tol}
                         DIR vfpprod_spe1)

add_test_compareECLFiles(CASENAME spe1_spider
                         FILENAME SPIDER_CAKESLICE
                         SIMULATOR flow
                         ABS_TOL ${abs_tol}
                         REL_TOL ${rel_tol}
                         DIR radial_grid)

add_test_compareECLFiles(CASENAME spe1_radial
                         FILENAME RADIAL_CAKESLICE
                         SIMULATOR flow
                         ABS_TOL ${abs_tol}
                         REL_TOL ${rel_tol}
                         DIR radial_grid)

add_test_compareECLFiles(CASENAME jfunc_01
                         FILENAME JFUNC-01
                         SIMULATOR flow
                         ABS_TOL ${abs_tol}
                         REL_TOL ${rel_tol}
                         DIR jfunc
                         TEST_ARGS --enable-tuning=true)

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
                         TEST_ARGS --tolerance-cnv=0.00003 --time-step-control=pid --linear-solver=cpr_trueimpes)

add_test_compareECLFiles(CASENAME numerical_aquifer_3d_1aqu
                         FILENAME 3D_1AQU_3CELLS
                         SIMULATOR flow
                         ABS_TOL ${abs_tol}
                         REL_TOL ${rel_tol}
                         DIR aquifer-num
                         TEST_ARGS --tolerance-cnv=0.00003 --time-step-control=pid --linear-solver=cpr_trueimpes)

add_test_compareECLFiles(CASENAME aquflux_01
                         FILENAME AQUFLUX-01
                         SIMULATOR flow
                         ABS_TOL ${abs_tol}
                         REL_TOL ${rel_tol}
                         DIR aquifers
                         TEST_ARGS --enable-tuning=true --relaxed-max-pv-fraction=0)

add_test_compareECLFiles(CASENAME aquflux_02
                         FILENAME AQUFLUX-02
                         SIMULATOR flow
                         ABS_TOL ${abs_tol}
                         REL_TOL ${rel_tol}
                         DIR aquifers
                         TEST_ARGS --solver-max-time-step-in-days=1)

add_test_compareECLFiles(CASENAME spe3
                         FILENAME SPE3CASE1
                         SIMULATOR flow
                         ABS_TOL ${abs_tol}
                         REL_TOL ${coarse_rel_tol}
                         TEST_ARGS --tolerance-wells=1e-6)

add_test_compareECLFiles(CASENAME spe9
                         FILENAME SPE9_CP_SHORT
                         SIMULATOR flow
                         ABS_TOL ${abs_tol}
                         REL_TOL ${rel_tol}
                         RESTART_STEP 10)

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
                         TEST_ARGS --solver-max-time-step-in-days=10)

add_test_compareECLFiles(CASENAME polymer_injectivity
                         FILENAME 2D_POLYMER_INJECTIVITY
                         SIMULATOR flow
                         ABS_TOL ${abs_tol}
                         REL_TOL ${rel_tol}
                         TEST_ARGS --tolerance-wells=1.e-6 --linear-solver=ilu0)

add_test_compareECLFiles(CASENAME polymer_simple2D
                         FILENAME 2D_THREEPHASE_POLY_HETER
                         SIMULATOR flow
                         ABS_TOL ${abs_tol}
                         REL_TOL ${coarse_rel_tol}
                         TEST_ARGS --tolerance-mb-relaxed=1.e-7)

add_test_compareECLFiles(CASENAME spe5
                         FILENAME SPE5CASE1
                         SIMULATOR flow
                         ABS_TOL ${abs_tol}
                         REL_TOL ${coarse_rel_tol})

add_test_compareECLFiles(CASENAME spe5_co2eor
                         FILENAME SPE5CASE1_DYN
                         SIMULATOR flow
                         ABS_TOL ${abs_tol}
                         REL_TOL ${coarse_rel_tol})

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

add_test_compareECLFiles(CASENAME base_model2_welpi
                         FILENAME 0B_WELPI_MODEL2
                         SIMULATOR flow
                         ABS_TOL ${abs_tol}
                         REL_TOL ${rel_tol}
                         DIR model2
                         TEST_ARGS --enable-tuning=true)

add_test_compareECLFiles(CASENAME 0a1_grpctl_msw_model2
                         FILENAME 0A1_GRCTRL_LRAT_ORAT_BASE_MODEL2_MSW
                         SIMULATOR flow
                         ABS_TOL ${abs_tol}
                         REL_TOL ${rel_tol}
                         DIR model2
                         TEST_ARGS --solver-max-time-step-in-days=3)

add_test_compareECLFiles(CASENAME 0a2_grpctl_msw_model2
                         FILENAME 0A2_GRCTRL_LRAT_ORAT_GGR_BASE_MODEL2_MSW
                         SIMULATOR flow
                         ABS_TOL ${abs_tol}
                         REL_TOL ${rel_tol}
                         DIR model2
                         TEST_ARGS --solver-max-time-step-in-days=3)

add_test_compareECLFiles(CASENAME 0a3_grpctl_msw_model2
                         FILENAME 0A3_GRCTRL_LRAT_LRAT_BASE_MODEL2_MSW
                         SIMULATOR flow
                         ABS_TOL ${abs_tol}
                         REL_TOL ${rel_tol}
                         DIR model2
                         TEST_ARGS --solver-max-time-step-in-days=3)

add_test_compareECLFiles(CASENAME 0a4_grpctl_msw_model2
                         FILENAME 0A4_GRCTRL_LRAT_LRAT_GGR_BASE_MODEL2_MSW
                         SIMULATOR flow
                         ABS_TOL ${abs_tol}
                         REL_TOL ${rel_tol}
                         DIR model2
                         TEST_ARGS --solver-max-time-step-in-days=3)

add_test_compareECLFiles(CASENAME udq_actionx
                         FILENAME UDQ_ACTIONX
                         SIMULATOR flow
                         ABS_TOL ${abs_tol}
                         REL_TOL ${rel_tol}
                         DIR udq_actionx
                         TEST_ARGS --solver-max-time-step-in-days=5)

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
                         TEST_ARGS --solver-max-time-step-in-days=0.2)

add_test_compareECLFiles(CASENAME waghyst1
                         FILENAME WAGHYSTR-01
                         SIMULATOR flow
                         ABS_TOL ${abs_tol}
                         REL_TOL ${rel_tol}
                         DIR waghystr)

add_test_compareECLFiles(CASENAME waghyst2
                         FILENAME WAGHYSTR-02
                         SIMULATOR flow
                         ABS_TOL ${abs_tol}
                         REL_TOL ${rel_tol}
                         DIR waghystr)

add_test_compareECLFiles(CASENAME gpmaint11
                         FILENAME GPMAINT-11
                         SIMULATOR flow
                         ABS_TOL ${abs_tol}
                         REL_TOL ${rel_tol}
                         DIR gpmaint)

add_test_compareECLFiles(CASENAME 3dwecon9
                         FILENAME 3D_WECON_9
                         SIMULATOR flow
                         ABS_TOL ${abs_tol}
                         REL_TOL ${rel_tol}
                         DIR wecon_wtest)

set(_gconprod_cases
  T1L
  T1W
  T2G
  T2O
)

add_multiple_tests(
  _gconprod_cases
  gconprod_
  SIMULATOR flow
  ABS_TOL ${abs_tol}
  REL_TOL ${rel_tol}
  DIR gconprod
)

set(_pinch_cases
  PINCH_MULTZ_ALL
  PINCH_MULTZ-_ALL
  PINCH_MULTZ_ALL_BARRIER
  PINCH_MULTZ-_ALL_BARRIER
  PINCH10_NOPINCH
  T1A_GAP T1A_NOGAP T1A_NOPINCH
  T1A1_NOGAP
  T2A1_GAP
  T2A_NOPINCH T2A_GAP
  T1B_NOPINCH
  T1B1_GAP
  T1B2_GAP
  T1B3_GAP
  T1B4_GAP
  T1B5_GAP
  T1B6_GAP
  T1B7_GAP
  T1C_NOPINCH
  T1C1_NOGAP T1C1_GAP
  T1C2_GAP T1C2_NOGAP
  T1C3_GAP T1C3_NOGAP
  T1D_NOPINCH
  T1D1_GAP T1D1_NOGAP)

add_multiple_tests(
  _pinch_cases
  pinch_
  SIMULATOR flow
  ABS_TOL ${abs_tol}
  REL_TOL ${rel_tol}
  DIR pinch
)

set(_udt_cases
  UDT-1D-01B
  UDT-1D-01
  UDT-1D-02
  UDT-1D-03
)

add_multiple_tests(
  _udt_cases
  ""
  SIMULATOR flow
  ABS_TOL ${abs_tol}
  REL_TOL ${rel_tol}
  TEST_ARGS --enable-tuning=true
  DIR udt
)

add_multiple_test_range(
  1
  6
  EQUALREG-0
  equalreg_multy
  SIMULATOR flow
  ABS_TOL ${abs_tol}
  REL_TOL ${rel_tol}
  DIR mult
)

add_multiple_test_range(
  1
  6
  ACTIONX_WELL_TEMPL-0
  actionx_well_templ
  SIMULATOR flow
  ABS_TOL ${abs_tol}
  REL_TOL ${rel_tol}
  DIR actionx
)

add_multiple_test_range(
  0
  8
  WCYCLE-
  WCYCLE
  SIMULATOR flow
  ABS_TOL ${abs_tol}
  REL_TOL ${rel_tol}
  DIR wcycle
  TEST_ARGS --enable-tuning=true
)

add_test_compareECLFiles(CASENAME udq_uadd
                         FILENAME UDQ_M1
                         SIMULATOR flow
                         ABS_TOL ${abs_tol}
                         REL_TOL ${rel_tol}
                         DIR udq_actionx
                         TEST_ARGS --tolerance-wells=1.e-8)

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

add_test_compareECLFiles(CASENAME reg_smry_in_fld_udq
                         FILENAME UDQ_REG-01
                         SIMULATOR flow
                         ABS_TOL ${abs_tol}
                         REL_TOL ${rel_tol}
                         DIR udq_actionx
                         TEST_ARGS --enable-tuning=true --time-step-control=pid)

# UDQ ASSIGN for subsets of group level UDQs.  Updates triggered from
# ACTIONX blocks.
add_test_compareECLFiles(CASENAME group_udq
                         FILENAME UDQ_GRP-01
                         SIMULATOR flow
                         ABS_TOL ${abs_tol}
                         REL_TOL ${rel_tol}
                         DIR udq_actionx
                         TEST_ARGS --solver-max-time-step-in-days=0.25)

set(_actionx_tests
  ACTIONX_GCONINJE
  ACTIONX_GCONPROD
  ACTIONX_UDQ
  ACTIONX_WCONHIST
  ACTIONX_WCONINJH
  ACTIONX_WEFAC
  UDQ-01
)

add_multiple_tests(
  _actionx_tests
  ""
  SIMULATOR flow
  ABS_TOL ${abs_tol}
  REL_TOL ${rel_tol}
  DIR actionx
)

add_multiple_test_range(
  1
  5
  CSKIN-0
  cskin
  SIMULATOR flow
  ABS_TOL ${abs_tol}
  REL_TOL ${rel_tol}
  DIR cskin
)

set(_co2store_cases
  CO2STORE
  CO2STORE_DIFFUSIVE
  CO2STORE_DRSDTCON
  CO2STORE_ENERGY
  CO2STORE_GASWAT
  CO2STORE_GW
  CO2STORE_GW_DIRICHLET
)

add_multiple_tests(
  _co2store_cases
  ""
  SIMULATOR flow
  ABS_TOL ${abs_tol}
  REL_TOL ${rel_tol}
  DIR co2store
)

set(_h2store_cases
  H2STORE
  H2STORE_DIFFUSIVE
  H2STORE_ENERGY
  H2STORE_GASWAT
  H2STORE_GW
)

add_multiple_tests(
  _h2store_cases
  ""
  SIMULATOR flow
  ABS_TOL ${abs_tol}
  REL_TOL ${rel_tol}
  DIR h2store
  TEST_ARGS --tolerance-cnv-relaxed=0.01
)

add_test_compareECLFiles(CASENAME ppcwmax
                         FILENAME PPCWMAX-01
                         SIMULATOR flow
                         ABS_TOL ${abs_tol}
                         REL_TOL ${rel_tol}
                         DIR ppcwmax)

if (opm-common_EMBEDDED_PYTHON)
  add_test_compareECLFiles(CASENAME udq_pyaction
                           FILENAME PYACTION_WCONPROD
                           SIMULATOR flow
                           ABS_TOL ${abs_tol}
                           REL_TOL ${rel_tol}
                           DIR udq_actionx
                           TEST_ARGS --solver-max-time-step-in-days=10)
  if(CMAKE_VERSION VERSION_GREATER_EQUAL 3.22)
    set_tests_properties(compareECLFiles_flow+PYACTION_WCONPROD PROPERTIES
                          ENVIRONMENT_MODIFICATION PYTHONPATH=path_list_append:${opm-common_DIR}/python)
  endif()
endif()

set(_model2_tests
  0_BASE_MODEL2
  0_BASE_MODEL2_LET
  0A1_GRCTRL_LRAT_ORAT_BASE_MODEL2_STW
  0A2_GRCTRL_LRAT_ORAT_GGR_BASE_MODEL2_STW
  0A3_GRCTRL_LRAT_LRAT_BASE_MODEL2_STW
  0A4_GRCTRL_LRAT_LRAT_GGR_BASE_MODEL2_STW
  0C_BASE_FBHPDEF
  1_MULTREGT_MODEL2
  2_MULTXYZ_MODEL2
  3_A_MPI_MULTFLT_SCHED_MODEL2
  5_SWATINIT_MODEL2
  6_ENDSCALE_MODEL2
  7_HYSTERESIS_MODEL2
  8_MULTIPLY_TRANXYZ_MODEL2
  9_EDITNNC_MODEL2
  9_1A_DEPL_MAX_RATE_MIN_BHP_STW
  9_1A_DEPL_MAX_RATE_MIN_BHP_MSW
  9_1B_DEPL_MAX_RATE_MIN_THP_STW
  9_1B_DEPL_MAX_RATE_MIN_THP_MSW
  9_2A_DEPL_GCONPROD_1L_STW
  9_2A_DEPL_GCONPROD_1L_MSW
  9_2B_DEPL_GCONPROD_2L_STW
  9_2B_DEPL_GCONPROD_2L_MSW
  9_3A_GINJ_REIN-G_STW
  9_3A_GINJ_REIN-G_MSW
  9_3B_GINJ_GAS_EXPORT_STW
  9_3B_GINJ_GAS_EXPORT_MSW
  9_3C_GINJ_GAS_GCONSUMP_STW
  9_3C_GINJ_GAS_GCONSUMP_MSW
  9_3D_GINJ_GAS_MAX_EXPORT_MSW
  9_3E_GAS_MIN_EXPORT_STW
  9_3E_GAS_MIN_EXPORT_MSW
  9_4A_WINJ_MAXWRATES_MAXBHP_GCONPROD_1L_STW
  9_4A_WINJ_MAXWRATES_MAXBHP_GCONPROD_1L_MSW
  9_4B_WINJ_VREP-W_STW
  9_4B_WINJ_VREP-W_MSW
  9_4C_WINJ_GINJ_VREP-W_REIN-G_STW
  9_4C_WINJ_GINJ_VREP-W_REIN-G_MSW
  9_4D_WINJ_GINJ_GAS_EXPORT_STW
  9_4D_WINJ_GINJ_GAS_EXPORT_MSW
)

add_multiple_tests(
  _model2_tests
  ""
  SIMULATOR flow
  ABS_TOL ${abs_tol}
  REL_TOL ${rel_tol}
  DIR model2
)

add_test_compareECLFiles(CASENAME multflt_model2
			  FILENAME 3_MULTFLT_MODEL2
			  SIMULATOR flow
			  ABS_TOL ${abs_tol}
			  REL_TOL ${rel_tol}
			  DIR model2
			  TEST_ARGS --solver-max-time-step-in-days=10)

add_test_compareECLFiles(CASENAME multpvv_model2
			  FILENAME 4_MINPVV_MODEL2
			  SIMULATOR flow
			  ABS_TOL ${abs_tol}
			  REL_TOL ${rel_tol}
			  DIR model2
			  TEST_ARGS --solver-max-time-step-in-days=10)

add_test_compareECLFiles(CASENAME 9_3d_grpctl_stw_model2
                         FILENAME 9_3D_GINJ_GAS_MAX_EXPORT_STW
                         SIMULATOR flow
                         ABS_TOL ${abs_tol}
                         REL_TOL ${rel_tol}
                         DIR model2
                         TEST_ARGS --time-step-after-event-in-days=1)

add_test_compareECLFiles(CASENAME model4_udq_group
                         FILENAME MOD4_UDQ_ACTIONX
                         SIMULATOR flow
                         ABS_TOL ${abs_tol}
                         REL_TOL ${rel_tol}
                         DIR model4)

add_test_compareECLFiles(CASENAME model4_gefac
                         FILENAME MOD4_GRP_GEFAC
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

add_test_compareECLFiles(CASENAME wsegvalv
			  FILENAME BASE_MSW_WSEGVALV
			  SIMULATOR flow
			  ABS_TOL ${abs_tol}
			  REL_TOL ${rel_tol})

add_test_compareECLFiles(CASENAME wsegvalv_2d_vert
                         FILENAME  MSW-2D-VERT-02
                         SIMULATOR flow
                         ABS_TOL ${abs_tol}
                         REL_TOL ${rel_tol}
                         DIR msw)

add_test_compareECLFiles(CASENAME nnc
                         FILENAME NNC_AND_EDITNNC
                         SIMULATOR flow
                         ABS_TOL ${abs_tol}
                         REL_TOL ${rel_tol}
                         DIR editnnc
                         TEST_ARGS --enable-opm-rst-file=1)

add_test_compareECLFiles(CASENAME nonnc
                         FILENAME NONNC
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

add_test_compareECLFiles(CASENAME spe1_solvent_foam
                         FILENAME SPE1CASE2_SOLVENT_FOAM
                         SIMULATOR flow
                         ABS_TOL ${abs_tol}
                         REL_TOL ${rel_tol}
                         DIR spe1_solvent)

add_test_compareECLFiles(CASENAME spe1_gaswater_solvent
                         FILENAME SPE1CASE2_GASWATER_SOLVENT
                         SIMULATOR flow
                         ABS_TOL ${abs_tol}
                         REL_TOL ${rel_tol}
                         DIR spe1_solvent)

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
                         DIR norne
                         TEST_ARGS --tolerance-wells=1.e-8)

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
                         DIR parallel_fieldprops
                         TEST_ARGS --enable-tuning=true --relaxed-max-pv-fraction=0)

add_test_compareECLFiles(CASENAME micp
                         FILENAME MICP
                         SIMULATOR flow
                         ABS_TOL ${abs_tol}
                         REL_TOL ${rel_tol}
                         DIR micp)

add_test_compareECLFiles(CASENAME 0_base_model6
                         FILENAME 0_BASE_MODEL6
                         SIMULATOR flow
                         ABS_TOL ${abs_tol}
                         REL_TOL ${rel_tol}
                         DIR model6)

add_test_compareECLFiles(CASENAME 0a_aquct_model6
                         FILENAME 0A_AQUCT_MODEL6
                         SIMULATOR flow
                         ABS_TOL ${abs_tol}
                         REL_TOL ${rel_tol}
                         DIR model6)

add_test_compareECLFiles(CASENAME 0b_rocktab_model6
                         FILENAME 0B_ROCKTAB_MODEL6
                         SIMULATOR flow
                         ABS_TOL ${abs_tol}
                         REL_TOL ${rel_tol}
                         DIR model6)

add_test_compareECLFiles(CASENAME base_wt_tracer
                         FILENAME BASE_WT_TRACER
                         SIMULATOR flow
                         ABS_TOL ${abs_tol}
                         REL_TOL ${rel_tol}
                         DIR tracer
                         RESTART_STEP 1,3,7)

add_multiple_test_range(
  1
  3
  MIN_BHP_
  min_bhp
  SIMULATOR flow
  ABS_TOL ${abs_tol}
  REL_TOL ${rel_tol}
  DIR wtest/bhp_min
  TEST_ARGS --enable-tuning=true
)

add_test_compareECLFiles(CASENAME min_thp_1
                         FILENAME MIN_THP_1
                         SIMULATOR flow
                         ABS_TOL ${abs_tol}
                         REL_TOL ${rel_tol}
                         DIR wtest/thp_min
                         TEST_ARGS --enable-tuning=true)

add_test_compareECLFiles(CASENAME max_gor_1
                         FILENAME MAX_GOR_1
                         SIMULATOR flow
                         ABS_TOL ${abs_tol}
                         REL_TOL ${rel_tol}
                         DIR wtest/wecon_gor_max
                         TEST_ARGS --enable-tuning=true)

add_test_compareECLFiles(CASENAME min_gasrate_1
                         FILENAME MIN_GASRATE_1
                         SIMULATOR flow
                         ABS_TOL ${abs_tol}
                         REL_TOL ${rel_tol}
                         DIR wtest/wecon_qg_min
                         TEST_ARGS --enable-tuning=true)

add_test_compareECLFiles(CASENAME min_qoil_1
                         FILENAME MIN_QOIL_1
                         SIMULATOR flow
                         ABS_TOL ${abs_tol}
                         REL_TOL ${rel_tol}
                         DIR wtest/wecon_qo_min
                         TEST_ARGS --enable-tuning=true)

add_multiple_test_range(
  1
  4
  MAX_WATERCUT_
  max_watercut
  SIMULATOR flow
  ABS_TOL ${abs_tol}
  REL_TOL ${rel_tol}
  DIR wtest/wecon_wct_max
  TEST_ARGS --enable-tuning=true
)

add_test_compareECLFiles(CASENAME max_wgr_1
                         FILENAME MAX_WGR_1
                         SIMULATOR flow
                         ABS_TOL ${abs_tol}
                         REL_TOL ${rel_tol}
                         DIR wtest/wecon_wgr_max
                         TEST_ARGS --enable-tuning=true)

add_test_compareECLFiles(CASENAME rxft_smry
                         FILENAME TEST_RXFT
                         SIMULATOR flow
                         ABS_TOL ${abs_tol}
                         REL_TOL ${rel_tol}
                         DIR rxft_smry
                         TEST_ARGS --enable-tuning=true --linear-solver-reduction=1e-7)

add_test_compareECLFiles(CASENAME bo_diffusion
                         FILENAME BO_DIFFUSE_CASE1
                         SIMULATOR flow
                         ABS_TOL ${abs_tol}
                         REL_TOL ${rel_tol}
                         DIR diffusion )

add_test_compareECLFiles(CASENAME fpr_nonhc
                         FILENAME WATER2F
                         SIMULATOR flow
                         ABS_TOL ${abs_tol}
                         REL_TOL ${rel_tol}
                         DIR water-1ph)

add_test_compareECLFiles(CASENAME actionx_wpimult
                         FILENAME ACTIONX_WPIMULT
                         SIMULATOR flow
                         ABS_TOL ${abs_tol}
                         REL_TOL ${rel_tol}
                         DIR actionx
                         TEST_ARGS --enable-tuning=true)

add_test_compareECLFiles(CASENAME wvfpexp_02
                         FILENAME WVFPEXP-02
                         SIMULATOR flow
                         ABS_TOL ${abs_tol}
                         REL_TOL ${rel_tol}
                         DIR wvfpexp)

set(_krnum_tests
  KRNUM-02X
  KRNUM-02Y
  KRNUM-02Z
  KRNUM-03X
  KRNUM-03Y
  KRNUM-03Z
)

add_multiple_tests(
  _krnum_tests
  ""
  SIMULATOR flow
  ABS_TOL ${abs_tol}
  REL_TOL ${rel_tol}
  DIR krnum
)

set(_gridunit_tests
  M_FIELD_GRID_CM
  M_FIELD_GRID_FEET
  M_FIELD_GRID_METRES
  M_METRIC_GRID_CM
  M_METRIC_GRID_FEET
  M_METRIC_GRID_METRES
)

add_multiple_tests(
  _gridunit_tests
  ""
  SIMULATOR flow
  ABS_TOL ${abs_tol}
  REL_TOL ${rel_tol}
  DIR gridunit
)

add_test_compareECLFiles(CASENAME 01_wgrupcon
                         FILENAME 01-WGRUPCON
                         SIMULATOR flow
                         ABS_TOL ${abs_tol}
                         REL_TOL ${rel_tol}
                         DIR wgrupcon
                         TEST_ARGS --enable-tuning=true)

add_test_compareECLFiles(CASENAME 02_wgrupcon
                         FILENAME 02-WGRUPCON
                         SIMULATOR flow
                         ABS_TOL ${abs_tol}
                         REL_TOL ${rel_tol}
                         DIR wgrupcon
                         TEST_ARGS --enable-tuning=true)
add_test_compareECLFiles(CASENAME winjmult_stdw
                         FILENAME WINJMULT_STDW
                         SIMULATOR flow
                         ABS_TOL ${abs_tol}
                         REL_TOL ${rel_tol}
                         DIR winjmult
                         TEST_ARGS --enable-tuning=true)
add_test_compareECLFiles(CASENAME winjmult_msw
                         FILENAME WINJMULT_MSW
                         SIMULATOR flow
                         ABS_TOL ${abs_tol}
                         REL_TOL ${rel_tol}
                         DIR winjmult
                         TEST_ARGS --enable-tuning=true)
add_test_compareECLFiles(CASENAME winjdam_stdw
                         FILENAME WINJDAM_STDW
                         SIMULATOR flow
                         ABS_TOL ${abs_tol}
                         REL_TOL ${rel_tol}
                         DIR winjdam
                         TEST_ARGS --enable-tuning=true)
add_test_compareECLFiles(CASENAME winjdam_msw
                         FILENAME WINJDAM_MSW
                         SIMULATOR flow
                         ABS_TOL ${abs_tol}
                         REL_TOL ${rel_tol}
                         DIR winjdam
                         TEST_ARGS --enable-tuning=true)
add_test_compareECLFiles(CASENAME 01_vappars
                         FILENAME VAPPARS-01
                         SIMULATOR flow
                         ABS_TOL ${abs_tol}
                         REL_TOL ${rel_tol}
                         DIR vappars)

add_test_compareECLFiles(CASENAME 6_uda_model5_stdw
  FILENAME 6_UDA_MODEL5_STDW
  SIMULATOR flow
  ABS_TOL ${abs_tol}
  REL_TOL ${rel_tol}
  DIR model5
  RESTART_STEP 1
  TEST_ARGS --enable-tuning=true
)

set(_mult_cases
  MULTFLT-01
  MULTFLT-02
  MULTFLT-03
  MULTREGT-01
  MULTXYZ05
  MULTXYZ06
  MULTXYZ07
  MULTXYZ08
  MULTXYZ11
)

add_multiple_tests(
  _mult_cases
  ""
  SIMULATOR flow
  ABS_TOL ${abs_tol}
  REL_TOL ${rel_tol}
  DIR mult
)

add_test_compareECLFiles(CASENAME gsatprod
                         FILENAME GSATPROD2
                         SIMULATOR flow
                         ABS_TOL ${abs_tol}
                         REL_TOL ${rel_tol}
                         DIR satellite)
if(BUILD_FLOW_POLY_GRID)
  add_test_compareECLFiles(CASENAME spe12_polyhedralgrid
                           FILENAME SPE1CASE2
                           SIMULATOR flow_blackoil_polyhedralgrid
                           ABS_TOL ${abs_tol}
                           REL_TOL ${coarse_rel_tol}
                           DIR spe1)
endif()

if(dune-alugrid_FOUND AND BUILD_FLOW_ALU_GRID AND MPI_FOUND)
  add_test_compareECLFiles(CASENAME spe12_alugrid
                           FILENAME SPE1CASE2
                           SIMULATOR flow_blackoil_alugrid
                           ABS_TOL ${abs_tol}
                           REL_TOL ${coarse_rel_tol}
                           DIR spe1)
endif()

if(BUILD_FLOW_FLOAT_VARIANTS)
  add_test_compareECLFiles(CASENAME spe1_float
                           FILENAME SPE1CASE1
                           SIMULATOR flow_blackoil_float
                           ABS_TOL ${abs_tol}
                           REL_TOL ${rel_tol}
                           DIR spe1
                           TEST_ARGS --tolerance-mb=1e-6)
endif()
