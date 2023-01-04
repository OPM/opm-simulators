opm_set_test_driver(${CMAKE_CURRENT_SOURCE_DIR}/tests/run-vtu-test.sh --simulation)
opm_set_test_default_working_directory(${PROJECT_BINARY_DIR}/tests)

opm_add_test(art2dgf
    NO_COMPILE
    DRIVER_ARGS --plain
    TEST_ARGS "data/fracture-raw.art")

opm_add_test(waterair_pvs_ni
             NO_COMPILE
             TEST_ARGS --grid-global-refinements=1)

opm_add_test(lens_immiscible_vcfv_ad
             NO_COMPILE
             TEST_ARGS --end-time=3000)

opm_add_test(lens_immiscible_vcfv_fd
             NO_COMPILE
             TEST_ARGS --end-time=3000)

opm_add_test(lens_immiscible_ecfv_ad
             NO_COMPILE
             TEST_ARGS --end-time=3000)

opm_add_test(lens_immiscible_ecfv_ad_23
             NO_COMPILE
             TEST_ARGS --end-time=3000)

opm_add_test(lens_immiscible_ecfv_ad_trans
             NO_COMPILE
             TEST_ARGS --end-time=3000)

# this test is identical to the simulation of the lens problem that
# uses the element centered finite volume discretization in
# conjunction with automatic differentiation
# (lens_immiscible_ecfv_ad). The only difference is that it uses
# multiple compile units in order to ensure that eWoms code can be
# used within libraries that use the same type tag within multiple
# compile units.
opm_add_test(lens_immiscible_ecfv_ad_mcu
             ONLY_COMPILE
             SOURCES
                examples/lens_immiscible_ecfv_ad_cu1.cc
                examples/lens_immiscible_ecfv_ad_cu2.cc
                examples/lens_immiscible_ecfv_ad_main.cc
             LIBRARIES opmsimulators)


opm_add_test(finger_immiscible_ecfv
             NO_COMPILE
             CONDITION ${DUNE_ALUGRID_FOUND})

opm_add_test(finger_immiscible_vcfv
             NO_COMPILE
             CONDITION ${DUNE_ALUGRID_FOUND})

opm_add_test(finger_immiscible_ecfv_adaptive
             EXE_NAME finger_immiscible_ecfv
             CONDITION ${DUNE_ALUGRID_FOUND} AND ${DUNE_FEM_FOUND}
             NO_COMPILE
             TEST_ARGS --enable-grid-adaptation=true --end-time=25e3)

opm_add_test(reservoir_blackoil_vcfv
             NO_COMPILE
             TEST_ARGS --end-time=8750000)
opm_add_test(reservoir_blackoil_ecfv
             NO_COMPILE
             TEST_ARGS --end-time=8750000)
opm_add_test(reservoir_ncp_vcfv
             NO_COMPILE
             TEST_ARGS --end-time=8750000)
opm_add_test(reservoir_ncp_ecfv
             NO_COMPILE
             TEST_ARGS --end-time=8750000)

opm_add_test(fracture_discretefracture
             NO_COMPILE
             CONDITION ${DUNE_ALUGRID_FOUND}
             TEST_ARGS --end-time=400)

opm_add_test(test_propertysystem
             DRIVER_ARGS --plain)

opm_add_test(test_quadrature
             DRIVER_ARGS --plain)

# test for the parallelization of the element centered finite volume
# discretization (using the non-isothermal NCP model and the parallel
# AMG linear solver)
opm_add_test(co2injection_ncp_ni_ecfv_parallel
             EXE_NAME co2injection_ncp_ni_ecfv
             NO_COMPILE
             PROCESSORS 4
             CONDITION ${MPI_FOUND}
             DRIVER_ARGS --parallel-simulation=4)

# test for the parallelization of the vertex centered finite volume
# discretization (using BiCGSTAB + ILU0)
opm_add_test(obstacle_immiscible_parallel
             EXE_NAME obstacle_immiscible
             NO_COMPILE
             PROCESSORS 4
             CONDITION ${MPI_FOUND}
             DRIVER_ARGS --parallel-simulation=4
             TEST_ARGS --end-time=1 --initial-time-step-size=1)

# test for the parallel AMG linear solver using the vertex centered
# finite volume discretization
opm_add_test(lens_immiscible_vcfv_fd_parallel
             EXE_NAME lens_immiscible_vcfv_fd
             NO_COMPILE
             PROCESSORS 4
             CONDITION ${MPI_FOUND}
             DRIVER_ARGS --parallel-simulation=4
             TEST_ARGS --end-time=250 --initial-time-step-size=250)

opm_add_test(lens_immiscible_vcfv_ad_parallel
             EXE_NAME lens_immiscible_vcfv_ad
             NO_COMPILE
             PROCESSORS 4
             CONDITION ${MPI_FOUND}
             DRIVER_ARGS --parallel-simulation=4
             TEST_ARGS --end-time=250 --initial-time-step-size=250)

opm_add_test(lens_immiscible_ecfv_ad_parallel
             EXE_NAME lens_immiscible_ecfv_ad
             NO_COMPILE
             PROCESSORS 4
             CONDITION ${MPI_FOUND}
             DRIVER_ARGS --parallel-simulation=4
             TEST_ARGS --end-time=250 --initial-time-step-size=250)

opm_add_test(obstacle_immiscible_parameters
             EXE_NAME obstacle_immiscible
             NO_COMPILE
             DEPENDS obstacle_immiscible
             DRIVER_ARGS --parameters)

opm_add_test(obstacle_pvs_restart
             EXE_NAME obstacle_pvs
             NO_COMPILE
             DEPENDS obstacle_pvs
             DRIVER_ARGS --restart
             TEST_ARGS --pvs-verbosity=2 --end-time=30000)

opm_add_test(test_tasklets
             DRIVER_ARGS --plain)

opm_add_test(test_mpiutil
             PROCESSORS 4
             CONDITION ${MPI_FOUND} AND Boost_UNIT_TEST_FRAMEWORK_FOUND
             DRIVER_ARGS --parallel-program=4)
