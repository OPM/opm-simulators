opm_set_test_driver(${CMAKE_CURRENT_SOURCE_DIR}/tests/run-parallel-unitTest.sh "")

if(NOT TARGET Boost::unit_test_framework)
  return()
endif()

opm_add_test(test_ParallelFieldProps_np2
  EXE_TARGET
    test_ParallelFieldProps
  DRIVER_ARGS
    -n 2
  PROCESSORS
    2
)

opm_add_test(test_compequil_invalid_eqlnum_np2
  EXE_TARGET
    test_compequil
  DRIVER_ARGS
    -n 2
  TEST_ARGS
    --run_test=InvalidEqlnumFailsOnAllRanks,MismatchedEqlnumSizeFailsOnAllRanks
  PROCESSORS
    2
)

opm_add_test(test_compequil_solvd_np2
  EXE_TARGET
    test_compequil_solvd
  DRIVER_ARGS
    -n 2
  PROCESSORS
    2
)

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

opm_add_test(test_wellsolverecovery
  DEPENDS
    opmsimulators
  LIBRARIES
    opmsimulators
    Boost::unit_test_framework
  SOURCES
    tests/test_wellsolverecovery.cpp
  DRIVER_ARGS
    -n 2
  PROCESSORS
    2
)

opm_add_test(test_parallel_exception_sync
  DEPENDS
    opmsimulators
  LIBRARIES
    opmsimulators
    Boost::unit_test_framework
  SOURCES
    tests/test_parallel_exception_sync.cpp
  DRIVER_ARGS
    -n 4
  PROCESSORS
    4
)

# One rank fails alone, and must abort the run with a message. This entry
# runs the binary of the entry above (EXE_TARGET, so no SOURCES), only its
# test case LoneFailureAborts. An aborted run never exits with 0, so the test
# passes when its output matches PASS_REGULAR_EXPRESSION instead, whatever
# the exit code. It fails, whatever else the output says, when it matches
# FAIL_REGULAR_EXPRESSION: the message of a rank that came back from its
# wait instead of being ended by the abort.
opm_add_test(test_parallel_exception_sync_abort
  EXE_TARGET
    test_parallel_exception_sync
  DRIVER_ARGS
    -n 4
  TEST_ARGS
    --run_test=LoneFailureAborts
  PROCESSORS
    4
)
set_tests_properties(parallel_exception_sync_abort
  PROPERTIES
  PASS_REGULAR_EXPRESSION
    "did not report a failure within 1 s"
  FAIL_REGULAR_EXPRESSION
    "The run should have been aborted"
)

# Rank 0 waits in the barrier for a failed time step while the other ranks
# take the maximum at the end of a parallel try-catch region. Rank 0 must
# abort with its message, as above, and no rank may come back from its wait:
# the barrier and the maximum must not be matched with each other.
opm_add_test(test_parallel_exception_sync_barrier_vs_region
  EXE_TARGET
    test_parallel_exception_sync
  DRIVER_ARGS
    -n 4
  TEST_ARGS
    --run_test=BarrierMeetsRegionEnd
  PROCESSORS
    4
)
set_tests_properties(parallel_exception_sync_barrier_vs_region
  PROPERTIES
  PASS_REGULAR_EXPRESSION
    "did not report a failure within 1 s .when a time step failed"
  FAIL_REGULAR_EXPRESSION
    "The run should have been aborted"
)

opm_add_test(test_parallelwellinfo_mpi
  EXE_TARGET
    test_parallelwellinfo
  DRIVER_ARGS
    -n 4
  PROCESSORS
    4
)

opm_add_test(test_blackoil_amg_np2
  EXE_TARGET
    test_blackoil_amg
  DRIVER_ARGS
    -n 2
  TEST_ARGS
    --run_test=IsolatedRowsHaveConsistentInitialSolution,UmfpackRequiresSequentialSolve
  PROCESSORS
    2
)

opm_add_test(test_blackoil_amg_transposed_cpr_np2
  EXE_TARGET
    test_blackoil_amg
  DRIVER_ARGS
    -n 2
  TEST_ARGS
    --run_test=TransposedCprKeepsCopiesConsistent
  PROCESSORS
    2
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

opm_add_executable(
  TARGET
    test_parallel_regionsetvariabledescriptor
  SOURCES
    tests/test_parallel_regionsetvariabledescriptor.cpp
  LIBRARIES
    opmcommon
    opmsimulators
    Boost::unit_test_framework
)

foreach(NPROC 2 3 4)
  opm_add_test(test_parallel_regionsetvariabledescriptor_np${NPROC}
    EXE_TARGET
      test_parallel_regionsetvariabledescriptor
    DRIVER_ARGS
      -n ${NPROC}
    PROCESSORS
      ${NPROC}
  )
endforeach()

opm_add_executable(
  TARGET
    test_parallel_regionvariablevalues
  SOURCES
    tests/test_parallel_regionvariablevalues.cpp
  LIBRARIES
    opmcommon
    opmsimulators
    Boost::unit_test_framework
)

foreach(NPROC 2 3 4)
  opm_add_test(test_parallel_regionvariablevalues_np${NPROC}
    EXE_TARGET
      test_parallel_regionvariablevalues
    DRIVER_ARGS
      -n ${NPROC}
    PROCESSORS
      ${NPROC}
  )
endforeach()

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
  opm_add_test(test_parallel_numaquifer_equilibration_check_np${NPROC}
    EXE_TARGET
      test_NumericalAquiferEquilibrationCheck
    DRIVER_ARGS
      -n ${NPROC}
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
