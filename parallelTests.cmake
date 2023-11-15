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
                                     TEST_ARGS --linear-solver-reduction=1e-7 --tolerance-cnv=5e-6 --tolerance-mb=1e-8 --ecl-enable-drift-compensation=false)

add_test_compare_parallel_simulation(CASENAME spe1_gaswater
                                     FILENAME SPE1CASE2_GASWATER
                                     SIMULATOR flow
                                     ABS_TOL ${abs_tol_parallel}
                                     REL_TOL ${rel_tol_parallel}
                                     DIR spe1
                                     TEST_ARGS --linear-solver-reduction=1e-7 --tolerance-cnv=5e-6 --tolerance-mb=1e-8)

add_test_compare_parallel_simulation(CASENAME spe9
                                     FILENAME SPE9_CP_SHORT
                                     SIMULATOR flow
                                     ABS_TOL ${abs_tol_parallel}
                                     REL_TOL ${rel_tol_parallel}
                                     TEST_ARGS --linear-solver-reduction=1e-7 --tolerance-cnv=5e-6 --tolerance-mb=1e-8 --ecl-enable-drift-compensation=false)

# A test for distributed standard wells. We load distribute only along the z-axis
add_test_compare_parallel_simulation(CASENAME spe9_dist_z
                                     FILENAME SPE9_CP_SHORT
                                     DIR spe9
                                     SIMULATOR flow_distribute_z
                                     ABS_TOL ${abs_tol_parallel}
                                     REL_TOL ${rel_tol_parallel}
                                     TEST_ARGS --linear-solver-reduction=1e-7 --tolerance-cnv=5e-6 --tolerance-mb=1e-8 --ecl-enable-drift-compensation=false)

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
                                     TEST_ARGS --linear-solver-reduction=1e-7 --tolerance-cnv=5e-6 --tolerance-mb=1e-8 --tolerance-wells=1e-7)

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
                                     REL_TOL ${coarse_rel_tol_parallel}}
                                     DIR spe1
                                     TEST_ARGS --linear-solver-reduction=1e-7 --tolerance-cnv=5e-6 --tolerance-mb=1e-8)

add_test_compare_parallel_simulation(CASENAME spe1_thermal_onephase
                                     FILENAME SPE1CASE2_THERMAL_ONEPHASE
                                     SIMULATOR flow
                                     ABS_TOL ${abs_tol}
                                     REL_TOL ${rel_tol}
                                     DIR spe1
                                     TEST_ARGS --linear-solver-reduction=1e-7 --tolerance-cnv=5e-6 --tolerance-mb=1e-8)

add_test_compare_parallel_simulation(CASENAME spe1_water
                                     FILENAME SPE1CASE1_WATER
                                     SIMULATOR flow
                                     ABS_TOL ${abs_tol}
                                     REL_TOL ${rel_tol}
                                     DIR spe1
                                     TEST_ARGS --linear-solver-reduction=1e-7 --tolerance-cnv=5e-6 --tolerance-mb=1e-8 --tolerance-wells=1e-7)

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
                                     ABS_TOL 0.17
                                     REL_TOL ${coarse_rel_tol_parallel}
                                     DIR aquifer-num
                                     TEST_ARGS --tolerance-cnv=0.000003 --time-step-control=pid --linear-solver=cpr_trueimpes)

add_test_compare_parallel_simulation(CASENAME aquflux_01
                                     FILENAME AQUFLUX-01
                                     SIMULATOR flow
                                     ABS_TOL ${abs_tol}
                                     REL_TOL ${coarse_rel_tol_parallel}
                                     DIR aquifers
                                     TEST_ARGS --enable-tuning=true)

add_test_compare_parallel_simulation(CASENAME aquflux_02
                                     FILENAME AQUFLUX-02
                                     SIMULATOR flow
                                     ABS_TOL ${abs_tol}
                                     REL_TOL ${coarse_rel_tol_parallel}
                                     DIR aquifers
                                     TEST_ARGS --enable-tuning=true)

add_test_compare_parallel_simulation(CASENAME network_balance_01
		                             FILENAME NETWORK-01
		                             SIMULATOR flow
		                             ABS_TOL 0.04
		                             REL_TOL 0.02
		                             DIR network
		                             TEST_ARGS --enable-tuning=true)

add_test_compare_parallel_simulation(CASENAME numerical_aquifer_3d_1aqu
                                     FILENAME 3D_1AQU_3CELLS
                                     SIMULATOR flow
                                     ABS_TOL ${abs_tol_parallel}
                                     REL_TOL ${coarse_rel_tol_parallel}
                                     DIR aquifer-num
                                     TEST_ARGS --tolerance-cnv=0.00003 --time-step-control=pid --linear-solver=cpr_trueimpes)

add_test_compare_parallel_simulation(CASENAME actionx_m1
                                     FILENAME ACTIONX_M1
                                     SIMULATOR flow
                                     ABS_TOL ${abs_tol_parallel}
                                     REL_TOL ${coarse_rel_tol_parallel}
                                     DIR udq_actionx
                                     TEST_ARGS --solver-max-time-step-in-days=0.2 --linear-solver-reduction=1e-7 --tolerance-cnv=5e-6 --tolerance-mb=1e-6)

add_test_compare_parallel_simulation(CASENAME winjmult_msw
                                     FILENAME WINJMULT_MSW
                                     SIMULATOR flow
                                     ABS_TOL ${abs_tol}
                                     REL_TOL ${rel_tol}
                                     DIR winjmult
                                     TEST_ARGS --ecl-enable-drift-compensation=false --enable-tuning=true --linear-solver-reduction=1e-7 --tolerance-cnv=5e-6 --tolerance-mb=1e-8)

add_test_compare_parallel_simulation(CASENAME winjdam_msw
                                     FILENAME WINJDAM_MSW
                                     SIMULATOR flow
                                     ABS_TOL ${abs_tol}
                                     REL_TOL ${rel_tol}
                                     DIR winjdam
                                     TEST_ARGS --ecl-enable-drift-compensation=false --enable-tuning=true --linear-solver-reduction=1e-7 --tolerance-cnv=5e-6 --tolerance-mb=1e-8)

add_test_compare_parallel_simulation(CASENAME 3_a_mpi_multflt_mod2
                                     FILENAME 3_A_MPI_MULTFLT_SCHED_MODEL2
                                     SIMULATOR flow
                                     ABS_TOL ${abs_tol_parallel}
                                     REL_TOL ${rel_tol_parallel}
      			       DIR model2
                                     TEST_ARGS --linear-solver-reduction=1e-7 --tolerance-cnv=5e-6 --tolerance-mb=1e-8 --ecl-enable-drift-compensation=false)

add_test_compare_parallel_simulation(CASENAME rxft
                                     FILENAME TEST_RXFT
                                     SIMULATOR flow
                                     ABS_TOL ${abs_tol_parallel}
                                     REL_TOL 1.0e-3
                                     DIR rxft_smry
                                     TEST_ARGS --enable-tuning=true --linear-solver-reduction=1e-7 --tolerance-cnv=5e-6 --tolerance-mb=1e-8 --ecl-enable-drift-compensation=false)
