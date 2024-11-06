# Regression tests
opm_set_test_driver(${PROJECT_SOURCE_DIR}/tests/run-regressionTest.sh "")

# Set absolute tolerance to be used passed to the macros in the following tests
set(abs_tol 2e-2)
set(rel_tol 1e-5)
set(coarse_rel_tol 1e-2)

add_test_compareECLFiles(CASENAME spe1flowexp
                         FILENAME SPE1CASE2
                         SIMULATOR flowexp_blackoil
                         ABS_TOL ${abs_tol}
                         REL_TOL ${rel_tol}
                         DIR spe1)

add_test_compareECLFiles(CASENAME spe12
                         FILENAME SPE1CASE2
                         SIMULATOR flow
                         ABS_TOL ${abs_tol}
                         REL_TOL ${coarse_rel_tol}
                         RESTART_SCHED false
                         RESTART_STEP 60
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

add_test_compareECLFiles(CASENAME spe1_gaswater
                         FILENAME SPE1CASE2_GASWATER
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

add_test_compareECLFiles(CASENAME spe1_precsalt
                         FILENAME SPE1CASE1_PRECSALT
                         SIMULATOR flow
                         ABS_TOL ${abs_tol}
                         REL_TOL ${rel_tol}
                         DIR spe1_precsalt)

add_test_compareECLFiles(CASENAME gas_precsalt
                         FILENAME GASWATER_VAPWAT_PRECSALT
                         SIMULATOR flow
                         ABS_TOL ${abs_tol}
                         REL_TOL ${rel_tol}
                         DIR spe1_precsalt)

add_test_compareECLFiles(CASENAME gasoil_precsalt
                         FILENAME GASCONDENSATE_VAPWAT_PRECSALT_REGRESSION
                         SIMULATOR flow
                         ABS_TOL ${abs_tol}
                         REL_TOL ${rel_tol}
                         DIR spe1_precsalt
                         TEST_ARGS --solver-max-time-step-in-days=0.05)

add_test_compareECLFiles(CASENAME spe1_brine_gaswater
                         FILENAME SPE1CASE2_BRINE_GASWATER
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

add_test_compareECLFiles(CASENAME spe1_water
                         FILENAME SPE1CASE1_WATER
                         SIMULATOR flow
                         ABS_TOL ${abs_tol}
                         REL_TOL ${rel_tol}
                         DIR spe1)

add_test_compareECLFiles(CASENAME spe1_thermal_onephase
                         FILENAME SPE1CASE2_THERMAL_ONEPHASE
                         SIMULATOR flow
                         ABS_TOL ${abs_tol}
                         REL_TOL ${rel_tol}
                         DIR spe1)

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

add_test_compareECLFiles(CASENAME pinch_nopinch_1x1x10
                         FILENAME PINCH10_NOPINCH
                         SIMULATOR flow
                         ABS_TOL ${abs_tol}
                         REL_TOL ${rel_tol}
                         DIR pinch)

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
                         TEST_ARGS --solver-max-time-step-in-days=1)

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
                         TEST_ARGS --tolerance-wells=1e-6 --newton-max-iterations=20)

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
                         TEST_ARGS --solver-max-time-step-in-days=10 --tolerance-mb=1.e-7)

add_test_compareECLFiles(CASENAME polymer_injectivity
                         FILENAME 2D_POLYMER_INJECTIVITY
                         SIMULATOR flow
                         ABS_TOL ${abs_tol}
                         REL_TOL ${rel_tol}
                         TEST_ARGS --tolerance-mb=1.e-7 --tolerance-wells=1.e-6 --linear-solver=ilu0)

add_test_compareECLFiles(CASENAME polymer_simple2D
                         FILENAME 2D_THREEPHASE_POLY_HETER
                         SIMULATOR flow
                         ABS_TOL ${abs_tol}
                         REL_TOL ${coarse_rel_tol}
                         TEST_ARGS --tolerance-mb=1.e-7 --tolerance-mb-relaxed=1.e-7)

add_test_compareECLFiles(CASENAME spe5
                         FILENAME SPE5CASE1
                         SIMULATOR flow
                         ABS_TOL ${abs_tol}
                         REL_TOL ${coarse_rel_tol}
                         TEST_ARGS --newton-max-iterations=20)

add_test_compareECLFiles(CASENAME spe5_co2eor
                         FILENAME SPE5CASE1_DYN
                         SIMULATOR flow
                         ABS_TOL ${abs_tol}
                         REL_TOL ${coarse_rel_tol}
                         TEST_ARGS --newton-max-iterations=20)

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

add_test_compareECLFiles(CASENAME pinch_multz_all
                         FILENAME PINCH_MULTZ_ALL
                         SIMULATOR flow
                         ABS_TOL ${abs_tol}
                         REL_TOL ${rel_tol}
                         DIR pinch)

add_test_compareECLFiles(CASENAME pinch_multzm_all
                         FILENAME PINCH_MULTZ-_ALL
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

add_test_compareECLFiles(CASENAME pinch_multzm_all_barrier
                         FILENAME PINCH_MULTZ-_ALL_BARRIER
                         SIMULATOR flow
                         ABS_TOL ${abs_tol}
                         REL_TOL ${rel_tol}
                         DIR pinch)

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

set(_pinch_cases
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

foreach(pinch_case ${_pinch_cases})
  string(TOLOWER ${pinch_case} pinch_test)
  add_test_compareECLFiles(CASENAME pinch_${pinch_test}
    FILENAME ${pinch_case}
    SIMULATOR flow
    ABS_TOL ${abs_tol}
    REL_TOL ${rel_tol}
    DIR pinch)
endforeach()

set(_udt_cases
  1D-01B
  1D-01
  1D-02
  1D-03
)

foreach(udt_case ${_udt_cases})
  string(TOLOWER ${udt_case} udt_test)
  add_test_compareECLFiles(CASENAME udt_${udt_test}
    FILENAME UDT-${udt_case}
    SIMULATOR flow
    ABS_TOL ${abs_tol}
    REL_TOL ${rel_tol}
    TEST_ARGS --enable-tuning=true
    DIR udt)
endforeach()

foreach(eqreg_case RANGE 1 6)
  add_test_compareECLFiles(CASENAME equalreg_multy_0${eqreg_case}
    FILENAME EQUALREG-0${eqreg_case}
    SIMULATOR flow
    ABS_TOL ${abs_tol}
    REL_TOL ${rel_tol}
    DIR mult
  )
endforeach()

foreach(templ_case RANGE 1 6)
  add_test_compareECLFiles(CASENAME actionx_well_templ_0${templ_case}
    FILENAME ACTIONX_WELL_TEMPL-0${templ_case}
    SIMULATOR flow
    ABS_TOL ${abs_tol}
    REL_TOL ${rel_tol}
    DIR actionx
  )
endforeach()

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
                         TEST_ARGS --enable-tuning=true)

add_test_compareECLFiles(CASENAME udq_undefined_2
                         FILENAME UDQ-01
                         SIMULATOR flow
                         ABS_TOL ${abs_tol}
                         REL_TOL ${rel_tol}
                         DIR actionx)

add_test_compareECLFiles(CASENAME co2store
                         FILENAME CO2STORE
                         SIMULATOR flow
                         ABS_TOL ${abs_tol}
                         REL_TOL ${rel_tol}
                         DIR co2store)

add_test_compareECLFiles(CASENAME co2store_gw
                         FILENAME CO2STORE_GW
                         SIMULATOR flow
                         ABS_TOL ${abs_tol}
                         REL_TOL ${rel_tol}
                         DIR co2store)

add_test_compareECLFiles(CASENAME co2store_gw_dirichlet
                         FILENAME CO2STORE_GW_DIRICHLET
                         SIMULATOR flow
                         ABS_TOL ${abs_tol}
                         REL_TOL ${rel_tol}
                         DIR co2store)

add_test_compareECLFiles(CASENAME co2store_gaswat
                         FILENAME CO2STORE_GASWAT
                         SIMULATOR flow
                         ABS_TOL ${abs_tol}
                         REL_TOL ${rel_tol}
                         DIR co2store)

add_test_compareECLFiles(CASENAME ppcwmax
                         FILENAME PPCWMAX-01
                         SIMULATOR flow
                         ABS_TOL ${abs_tol}
                         REL_TOL ${rel_tol}
                         DIR ppcwmax)

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

add_test_compareECLFiles(CASENAME co2store_energy
                         FILENAME CO2STORE_ENERGY
                         SIMULATOR flow
                         ABS_TOL ${abs_tol}
                         REL_TOL ${rel_tol}
                         DIR co2store)

add_test_compareECLFiles(CASENAME h2store
                         FILENAME H2STORE
                         SIMULATOR flow
                         ABS_TOL ${abs_tol}
                         REL_TOL ${rel_tol}
                         DIR h2store
                         TEST_ARGS --tolerance-cnv-relaxed=0.01)

add_test_compareECLFiles(CASENAME h2store_gw
                         FILENAME H2STORE_GW
                         SIMULATOR flow
                         ABS_TOL ${abs_tol}
                         REL_TOL ${rel_tol}
                         DIR h2store
                         TEST_ARGS --tolerance-cnv-relaxed=0.01)

add_test_compareECLFiles(CASENAME h2store_gaswat
                         FILENAME H2STORE_GASWAT
                         SIMULATOR flow
                         ABS_TOL ${abs_tol}
                         REL_TOL ${rel_tol}
                         DIR h2store
                         TEST_ARGS --tolerance-cnv-relaxed=0.01)

add_test_compareECLFiles(CASENAME h2store_diffusive
                         FILENAME H2STORE_DIFFUSIVE
                         SIMULATOR flow
                         ABS_TOL ${abs_tol}
                         REL_TOL ${rel_tol}
                         DIR h2store
                         TEST_ARGS --tolerance-cnv-relaxed=0.01)

add_test_compareECLFiles(CASENAME h2store_energy
                         FILENAME H2STORE_ENERGY
                         SIMULATOR flow
                         ABS_TOL ${abs_tol}
                         REL_TOL ${rel_tol}
                         DIR h2store
                         TEST_ARGS --tolerance-cnv-relaxed=0.01)

if (opm-common_EMBEDDED_PYTHON)
  add_test_compareECLFiles(CASENAME udq_pyaction
                           FILENAME PYACTION_WCONPROD
                           SIMULATOR flow
                           ABS_TOL ${abs_tol}
                           REL_TOL ${rel_tol}
                           DIR udq_actionx
                           TEST_ARGS --solver-max-time-step-in-days=10)
endif()

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
                         DIR parallel_fieldprops)


add_test_compareECLFiles(CASENAME actionx_gconprod
                         FILENAME ACTIONX_GCONPROD
                         SIMULATOR flow
                         ABS_TOL ${abs_tol}
                         REL_TOL ${rel_tol}
                         DIR actionx)


add_test_compareECLFiles(CASENAME actionx_wefac
                         FILENAME ACTIONX_WEFAC
                         SIMULATOR flow
                         ABS_TOL ${abs_tol}
                         REL_TOL ${rel_tol}
                         DIR actionx)

add_test_compareECLFiles(CASENAME actionx_udq
                         FILENAME ACTIONX_UDQ
                         SIMULATOR flow
                         ABS_TOL ${abs_tol}
                         REL_TOL ${rel_tol}
                         DIR actionx)

add_test_compareECLFiles(CASENAME micp
                         FILENAME MICP
                         SIMULATOR flow
                         ABS_TOL ${abs_tol}
                         REL_TOL ${rel_tol}
                         DIR micp
                         TEST_ARGS --linear-solver=ilu0)

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


add_test_compareECLFiles(CASENAME min_bhp_1
                         FILENAME MIN_BHP_1
                         SIMULATOR flow
                         ABS_TOL ${abs_tol}
                         REL_TOL ${rel_tol}
                         DIR wtest/bhp_min
                         TEST_ARGS --enable-tuning=true)

add_test_compareECLFiles(CASENAME min_bhp_2
                         FILENAME MIN_BHP_2
                         SIMULATOR flow
                         ABS_TOL ${abs_tol}
                         REL_TOL ${rel_tol}
                         DIR wtest/bhp_min
                         TEST_ARGS --enable-tuning=true)

add_test_compareECLFiles(CASENAME min_bhp_3
                         FILENAME MIN_BHP_3
                         SIMULATOR flow
                         ABS_TOL ${abs_tol}
                         REL_TOL ${rel_tol}
                         DIR wtest/bhp_min
                         TEST_ARGS --enable-tuning=true)

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

add_test_compareECLFiles(CASENAME max_watercut_1
                         FILENAME MAX_WATERCUT_1
                         SIMULATOR flow
                         ABS_TOL ${abs_tol}
                         REL_TOL ${rel_tol}
                         DIR wtest/wecon_wct_max
                         TEST_ARGS --enable-tuning=true)

add_test_compareECLFiles(CASENAME max_watercut_2
                         FILENAME MAX_WATERCUT_2
                         SIMULATOR flow
                         ABS_TOL ${abs_tol}
                         REL_TOL ${rel_tol}
                         DIR wtest/wecon_wct_max
                         TEST_ARGS --enable-tuning=true)

add_test_compareECLFiles(CASENAME max_watercut_3
                         FILENAME MAX_WATERCUT_3
                         SIMULATOR flow
                         ABS_TOL ${abs_tol}
                         REL_TOL ${rel_tol}
                         DIR wtest/wecon_wct_max
                         TEST_ARGS --enable-tuning=true)

add_test_compareECLFiles(CASENAME max_watercut_4
                         FILENAME MAX_WATERCUT_4
                         SIMULATOR flow
                         ABS_TOL ${abs_tol}
                         REL_TOL ${rel_tol}
                         DIR wtest/wecon_wct_max
                         TEST_ARGS --enable-tuning=true)

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
                         TEST_ARGS --enable-tuning=true)

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

add_test_compareECLFiles(CASENAME spe1case2_krnum
                         FILENAME SPE1CASE2_KRNUM
                         SIMULATOR flow
                         ABS_TOL ${abs_tol}
                         REL_TOL ${rel_tol}
                         DIR spe1)

add_test_compareECLFiles(CASENAME model_field_gridunit_cm
                         FILENAME M_FIELD_GRID_CM
                         SIMULATOR flow
                         ABS_TOL ${abs_tol}
                         REL_TOL ${rel_tol}
                         DIR gridunit)

add_test_compareECLFiles(CASENAME model_field_gridunit_feet
                         FILENAME M_FIELD_GRID_FEET
                         SIMULATOR flow
                         ABS_TOL ${abs_tol}
                         REL_TOL ${rel_tol}
                         DIR gridunit)

add_test_compareECLFiles(CASENAME model_field_gridunit_metres
                         FILENAME M_FIELD_GRID_METRES
                         SIMULATOR flow
                         ABS_TOL ${abs_tol}
                         REL_TOL ${rel_tol}
                         DIR gridunit)

add_test_compareECLFiles(CASENAME model_metric_gridunit_cm
                         FILENAME M_METRIC_GRID_CM
                         SIMULATOR flow
                         ABS_TOL ${abs_tol}
                         REL_TOL ${rel_tol}
                         DIR gridunit)

add_test_compareECLFiles(CASENAME model_metric_gridunit_feet
                         FILENAME M_METRIC_GRID_FEET
                         SIMULATOR flow
                         ABS_TOL ${abs_tol}
                         REL_TOL ${rel_tol}
                         DIR gridunit)

add_test_compareECLFiles(CASENAME model_metric_gridunit_metres
                         FILENAME M_METRIC_GRID_METRES
                         SIMULATOR flow
                         ABS_TOL ${abs_tol}
                         REL_TOL ${rel_tol}
                         DIR gridunit)

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
add_test_compareECLFiles(CASENAME 01_multflt
                         FILENAME MULTFLT-01
                         SIMULATOR flow
                         ABS_TOL ${abs_tol}
                         REL_TOL ${rel_tol}
                         DIR mult)
add_test_compareECLFiles(CASENAME 02_multflt
                         FILENAME MULTFLT-02
                         SIMULATOR flow
                         ABS_TOL ${abs_tol}
                         REL_TOL ${rel_tol}
                         DIR mult)
add_test_compareECLFiles(CASENAME 03_multflt
                         FILENAME MULTFLT-03
                         SIMULATOR flow
                         ABS_TOL ${abs_tol}
                         REL_TOL ${rel_tol}
                         DIR mult)
add_test_compareECLFiles(CASENAME 01_multregt
                         FILENAME MULTREGT-01
                         SIMULATOR flow
                         ABS_TOL ${abs_tol}
                         REL_TOL ${rel_tol}
                         DIR mult)
add_test_compareECLFiles(CASENAME 05_multxyz
                         FILENAME MULTXYZ05
                         SIMULATOR flow
                         ABS_TOL ${abs_tol}
                         REL_TOL ${rel_tol}
                         DIR mult)
add_test_compareECLFiles(CASENAME 06_multxyz
                         FILENAME MULTXYZ06
                         SIMULATOR flow
                         ABS_TOL ${abs_tol}
                         REL_TOL ${rel_tol}
                         DIR mult)
add_test_compareECLFiles(CASENAME 07_multxyz
                         FILENAME MULTXYZ07
                         SIMULATOR flow
                         ABS_TOL ${abs_tol}
                         REL_TOL ${rel_tol}
                         DIR mult)
add_test_compareECLFiles(CASENAME 08_multxyz
                         FILENAME MULTXYZ08
                         SIMULATOR flow
                         ABS_TOL ${abs_tol}
                         REL_TOL ${rel_tol}
                         DIR mult)
add_test_compareECLFiles(CASENAME 11_multxyz
                         FILENAME MULTXYZ11
                         SIMULATOR flow
                         ABS_TOL ${abs_tol}
                         REL_TOL ${rel_tol}
                         DIR mult)

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

# Add tests from json file
function(AddJsonRegressionTests json_file)
  message(STATUS "Adding regression tests from ${json_file}")
  file(READ ${json_file} JSON_STRING)
  parse_json_common_params()
  string(JSON size LENGTH ${test_cases})
  math(EXPR COUNT "${size}-1")
  foreach(idx RANGE ${COUNT})
    parse_json_test_definition(${idx})
    cmake_parse_arguments(CURTEST "$" "CASENAME" "" ${PARAMS})
    message(STATUS "  - ${CURTEST_CASENAME}")
    add_test_compareECLFiles(${PARAMS})
  endforeach()
endfunction()

if(CMAKE_VERSION VERSION_GREATER_EQUAL 3.19)
  set(JSON_FILES 
        "tests/definitions/regression/cskin.json"
        "tests/definitions/regression/editnnc.json"
        "tests/definitions/regression/krnum.json"
        "tests/definitions/regression/model2.json"
        "tests/definitions/regression/network.json"
     )
  foreach(JSON_FILE ${JSON_FILES})
    AddJsonRegressionTests(${JSON_FILE})
    set_property(DIRECTORY APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS ${JSON_FILE})
  endforeach()
endif()
