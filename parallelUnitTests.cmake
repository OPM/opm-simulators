opm_set_test_driver(${CMAKE_CURRENT_SOURCE_DIR}/tests/run-parallel-unitTest.sh "")

if(NOT TARGET Boost::unit_test_framework)
  return()
endif()

opm_add_test(test_gatherconvergencereport
  DEPENDS
    opmsimulators
  LIBRARIES
    opmsimulators
    Boost::unit_test_framework
  SOURCES
    tests/test_gatherconvergencereport.cpp
  DRIVER_ARGS
    -n 4
  PROCESSORS
    4
)

opm_add_test(test_gatherdeferredlogger
  DEPENDS
    opmsimulators
  LIBRARIES
    opmsimulators
    Boost::unit_test_framework
  SOURCES
    tests/test_gatherdeferredlogger.cpp
  DRIVER_ARGS
    -n 4
  PROCESSORS
    4
)

opm_add_test(test_parallelwellinfo_mpi
  EXE_TARGET
    test_parallelwellinfo
  DRIVER_ARGS
    -n 4
  PROCESSORS
    4
)

foreach(NPROC 2 3 4)
  opm_add_test(test_parallel_wbp_sourcevalues_np${NPROC}
    EXE_TARGET
      test_parallel_wbp_sourcevalues
    DRIVER_ARGS
      -n ${NPROC}
    PROCESSORS
      ${NPROC}
  )
endforeach()

opm_add_executable(
  TARGET
    test_parallel_wbp_calculation
  SOURCES
    tests/test_parallel_wbp_calculation.cpp
  LIBRARIES
    opmcommon
    opmsimulators
    Boost::unit_test_framework
)

opm_add_test(test_parallel_wbp_calculation_create
  EXE_TARGET
    test_parallel_wbp_calculation
  DRIVER_ARGS
    -n 2
  TEST_ARGS
    --run_test=Create
  PROCESSORS
    2
)

opm_add_test(test_parallel_wbp_calculation_well_openconns
  EXE_TARGET
    test_parallel_wbp_calculation
  DRIVER_ARGS
    -n 2
  TEST_ARGS
    --run_test=TopOfFormation_Well_OpenConns
  PROCESSORS
    2
)

foreach(NPROC 2 3 4)
  opm_add_test(test_rftcontainer_np${NPROC}
    EXE_TARGET
      test_rftcontainer
    DRIVER_ARGS
      -n ${NPROC}
    PROCESSORS
      ${NPROC}
  )
endforeach()

foreach(NPROC 2 3 4)
  opm_add_test(test_parallel_region_phase_pvaverage_np${NPROC}
    EXE_TARGET
      test_region_phase_pvaverage
    DRIVER_ARGS
      -n ${NPROC}
    TEST_ARGS
      --run_test=Parallel/*
    PROCESSORS
      ${NPROC}
  )
endforeach()

foreach(NPROC 2 3 4)
  opm_add_test(test_parallel_satfunc_consistency_checks_np${NPROC}
    EXE_TARGET
      test_SatfuncConsistencyChecks_parallel
    DRIVER_ARGS
      -n ${NPROC}
    PROCESSORS
      ${NPROC}
  )
endforeach()

opm_add_test(test_broadcast
  DEPENDS
    opmsimulators
  LIBRARIES
    opmsimulators
    Boost::unit_test_framework
  SOURCES
    tests/test_broadcast.cpp
  DRIVER_ARGS
    -n 4
  PROCESSORS
    4
)

opm_add_test(test_HDF5File_Parallel
  DEPENDS
    opmsimulators
  LIBRARIES
    opmsimulators
    Boost::unit_test_framework
  SOURCES
    tests/test_HDF5File_Parallel.cpp
  CONDITION
    HDF5_FOUND
  DRIVER_ARGS
    -n 4
  PROCESSORS
    4
)

opm_add_test(test_HDF5Serializer_Parallel
  DEPENDS
    opmsimulators
  LIBRARIES
    opmcommon
    opmsimulators
    Boost::unit_test_framework
  SOURCES
    tests/test_HDF5Serializer_Parallel.cpp
  CONDITION
    HDF5_FOUND
  DRIVER_ARGS
    -n 4
  PROCESSORS
    4
)

opm_add_test(test_rstconv_parallel
  EXE_TARGET
    test_rstconv
  DRIVER_ARGS
    -n 4
  PROCESSORS
    4
)
