opm_set_test_driver(${CMAKE_CURRENT_SOURCE_DIR}/tests/run-parallel-unitTest.sh "")

opm_add_test(test_gatherconvergencereport
  DEPENDS "opmsimulators"
  LIBRARIES opmsimulators ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY}
  SOURCES
    tests/test_gatherconvergencereport.cpp
  CONDITION
    MPI_FOUND AND Boost_UNIT_TEST_FRAMEWORK_FOUND
  DRIVER_ARGS
    -n 4
    -b ${PROJECT_BINARY_DIR}
  PROCESSORS
    4
)

opm_add_test(test_gatherdeferredlogger
  DEPENDS "opmsimulators"
  LIBRARIES opmsimulators ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY}
  SOURCES
    tests/test_gatherdeferredlogger.cpp
  CONDITION
    MPI_FOUND AND Boost_UNIT_TEST_FRAMEWORK_FOUND
  DRIVER_ARGS
    -n 4
    -b ${PROJECT_BINARY_DIR}
  PROCESSORS
    4
)

opm_add_test(test_parallelwellinfo_mpi
  EXE_NAME
    test_parallelwellinfo
  CONDITION
    MPI_FOUND AND Boost_UNIT_TEST_FRAMEWORK_FOUND
  DRIVER_ARGS
    -n 4
    -b ${PROJECT_BINARY_DIR}
  NO_COMPILE
  PROCESSORS
    4
    )

foreach(NPROC 2 3 4)
  opm_add_test(test_parallel_wbp_sourcevalues_np${NPROC}
    EXE_NAME
      test_parallel_wbp_sourcevalues
    CONDITION
      MPI_FOUND AND Boost_UNIT_TEST_FRAMEWORK_FOUND
    DRIVER_ARGS
      -n ${NPROC}
      -b ${PROJECT_BINARY_DIR}
    NO_COMPILE
    PROCESSORS
      ${NPROC}
  )
endforeach()

opm_add_test(test_parallel_wbp_calculation
  SOURCES
    tests/test_parallel_wbp_calculation.cpp
  LIBRARIES
    opmsimulators ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY}
  CONDITION
    MPI_FOUND AND Boost_UNIT_TEST_FRAMEWORK_FOUND
  ONLY_COMPILE
)

opm_add_test(test_parallel_wbp_calculation_create
  EXE_NAME
    test_parallel_wbp_calculation
  CONDITION
    MPI_FOUND AND Boost_UNIT_TEST_FRAMEWORK_FOUND
  DRIVER_ARGS
    -n 2
    -b ${PROJECT_BINARY_DIR}
  TEST_ARGS
    --run_test=Create
  NO_COMPILE
  PROCESSORS
    2
)

opm_add_test(test_parallel_wbp_calculation_well_openconns
  EXE_NAME
    test_parallel_wbp_calculation
  CONDITION
    MPI_FOUND AND Boost_UNIT_TEST_FRAMEWORK_FOUND
  DRIVER_ARGS
    -n 2
    -b ${PROJECT_BINARY_DIR}
  TEST_ARGS
    --run_test=TopOfFormation_Well_OpenConns
  NO_COMPILE
  PROCESSORS
    2
)

foreach(NPROC 2 3 4)
  opm_add_test(test_rftcontainer_np${NPROC}
    EXE_NAME
      test_rftcontainer
    CONDITION
      MPI_FOUND AND Boost_UNIT_TEST_FRAMEWORK_FOUND
    DRIVER_ARGS
      -n ${NPROC}
      -b ${PROJECT_BINARY_DIR}
    NO_COMPILE
    PROCESSORS
      ${NPROC}
  )
endforeach()

foreach(NPROC 2 3 4)
  opm_add_test(test_parallel_region_phase_pvaverage_np${NPROC}
    EXE_NAME
      test_region_phase_pvaverage
    CONDITION
      MPI_FOUND AND Boost_UNIT_TEST_FRAMEWORK_FOUND
    DRIVER_ARGS
      -n ${NPROC}
      -b ${PROJECT_BINARY_DIR}
    TEST_ARGS
      --run_test=Parallel/*
    NO_COMPILE
    PROCESSORS
      ${NPROC}
  )
endforeach()

foreach(NPROC 2 3 4)
  opm_add_test(test_parallel_satfunc_consistency_checks_np${NPROC}
    EXE_NAME
      test_SatfuncConsistencyChecks_parallel
    CONDITION
      MPI_FOUND AND Boost_UNIT_TEST_FRAMEWORK_FOUND
    DRIVER_ARGS
      -n ${NPROC}
      -b ${PROJECT_BINARY_DIR}
    NO_COMPILE
    PROCESSORS
      ${NPROC}
  )
endforeach()

opm_add_test(test_broadcast
  DEPENDS "opmsimulators"
  LIBRARIES opmsimulators ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY}
  SOURCES
    tests/test_broadcast.cpp
  CONDITION
    MPI_FOUND AND Boost_UNIT_TEST_FRAMEWORK_FOUND
  DRIVER_ARGS
    -n 4
    -b ${PROJECT_BINARY_DIR}
  PROCESSORS
    4
)

opm_add_test(test_HDF5File_Parallel
  DEPENDS "opmsimulators"
  LIBRARIES opmsimulators ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY}
  SOURCES
    tests/test_HDF5File_Parallel.cpp
  CONDITION
    HDF5_FOUND AND MPI_FOUND AND Boost_UNIT_TEST_FRAMEWORK_FOUND
  DRIVER_ARGS
    -n 4
    -b ${PROJECT_BINARY_DIR}
  PROCESSORS
    4
)

opm_add_test(test_HDF5Serializer_Parallel
  DEPENDS "opmsimulators"
  LIBRARIES opmsimulators ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY}
  SOURCES
    tests/test_HDF5Serializer_Parallel.cpp
  CONDITION
    HDF5_FOUND AND MPI_FOUND AND Boost_UNIT_TEST_FRAMEWORK_FOUND
  DRIVER_ARGS
    -n 4
    -b ${PROJECT_BINARY_DIR}
  PROCESSORS
    4
)

opm_add_test(test_rstconv_parallel
  EXE_NAME
    test_rstconv
  CONDITION
    MPI_FOUND AND Boost_UNIT_TEST_FRAMEWORK_FOUND
  DRIVER_ARGS
    -n 4
    -b ${PROJECT_BINARY_DIR}
  NO_COMPILE
  PROCESSORS
    4
)
