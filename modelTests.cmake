opm_set_test_driver(${PROJECT_SOURCE_DIR}/tests/run-vtu-test.sh "--simulation")

opm_add_test(art2dgf
  EXE_TARGET
    art2dgf
  DRIVER_ARGS
    --plain
  TEST_ARGS
    data/fracture-raw.art
  WORKING_DIRECTORY
    ${PROJECT_BINARY_DIR}/tests
)

foreach(tgt lens_immiscible_ecfv_ad
            lens_immiscible_ecfv_ad_23
            lens_immiscible_ecfv_ad_trans
            lens_immiscible_vcfv_ad
            lens_immiscible_vcfv_fd)
  opm_add_test(${tgt}
    EXE_TARGET
      ${tgt}
    TEST_ARGS
      --end-time=3000
    WORKING_DIRECTORY
      ${PROJECT_BINARY_DIR}/tests
  )
endforeach()

opm_add_test(waterair_pvs_ni
  EXE_TARGET
    waterair_pvs_ni
  TEST_ARGS
    --grid-global-refinements=1
  WORKING_DIRECTORY
    ${PROJECT_BINARY_DIR}/tests
)

set(PLAIN_TGT
  co2injection_flash_ecfv
  co2injection_flash_ni_ecfv
  co2injection_flash_ni_vcfv
  co2injection_flash_vcfv
  co2injection_immiscible_ecfv
  co2injection_immiscible_ni_ecfv
  co2injection_immiscible_ni_vcfv
  co2injection_immiscible_vcfv
  co2injection_ncp_ecfv
  co2injection_ncp_ni_vcfv
  co2injection_ncp_vcfv
  co2injection_pvs_ecfv
  co2injection_pvs_ni_vcfv
  co2injection_pvs_vcfv
  co2injection_ncp_ni_ecfv
  co2injection_pvs_ni_ecfv
  co2_ptflash_ecfv
  cuvette_pvs
  diffusion_flash
  diffusion_ncp
  diffusion_pvs
  groundwater_immiscible
  infiltration_pvs
  lens_richards_ecfv
  lens_richards_vcfv
  obstacle_immiscible
  obstacle_ncp
  obstacle_pvs
  outflow_pvs
  powerinjection_darcy_ad
  powerinjection_darcy_fd
  powerinjection_forchheimer_ad
  powerinjection_forchheimer_fd
  tutorial1
)

if(dune-alugrid_FOUND)
  list(APPEND PLAIN_TGT
    finger_immiscible_ecfv
    finger_immiscible_vcfv
  )
endif()

foreach(tgt ${PLAIN_TGT})
  opm_add_test(${tgt}
    EXE_TARGET
      ${tgt}
    WORKING_DIRECTORY
      ${PROJECT_BINARY_DIR}/tests
  )
endforeach()

foreach(tgt reservoir_blackoil_ecfv
            reservoir_blackoil_vcfv
            reservoir_ncp_ecfv
            reservoir_ncp_vcfv)
  opm_add_test(${tgt}
    EXE_TARGET
      ${tgt}
    TEST_ARGS
      --end-time=8750000
    WORKING_DIRECTORY
      ${PROJECT_BINARY_DIR}/tests
  )
endforeach()

if(dune-alugrid_FOUND)
  opm_add_test(fracture_discretefracture
    EXE_TARGET
      fracture_discretefracture
    TEST_ARGS
      --end-time=400
    WORKING_DIRECTORY
      ${PROJECT_BINARY_DIR}/tests
  )
endif()

if(dune-alugrid_FOUND AND dune-fem_FOUND)
  opm_add_test(finger_immiscible_ecfv_adaptive
    EXE_TARGET
      finger_immiscible_ecfv
    TEST_ARGS
      --enable-grid-adaptation=true
      --end-time=25e3
      --enable-async-vtk-output=false
    WORKING_DIRECTORY
      ${PROJECT_BINARY_DIR}/tests
  )
endif()

opm_add_test(obstacle_immiscible_parameters
  EXE_TARGET
    obstacle_immiscible
  DRIVER_ARGS
    --parameters
  WORKING_DIRECTORY
    ${PROJECT_BINARY_DIR}/tests
)

opm_add_test(obstacle_pvs_restart
  EXE_TARGET
    obstacle_pvs
  TEST_ARGS
    --pvs-verbosity=2
    --end-time=30000
  DRIVER_ARGS
    --restart
  WORKING_DIRECTORY
    ${PROJECT_BINARY_DIR}/tests
)
