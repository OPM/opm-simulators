# Regression tests
opm_set_test_driver(${PROJECT_SOURCE_DIR}/tests/run-regressionTest.sh "")

# Set absolute tolerance to be used passed to the macros in the following tests
set(abs_tol 2e-2)
set(rel_tol 1e-5)
set(coarse_rel_tol 1e-2)

if (opm-common_EMBEDDED_PYTHON)
  add_test_compareECLFiles(CASENAME udq_pyaction
                           FILENAME PYACTION_WCONPROD
                           SIMULATOR flow
                           ABS_TOL ${abs_tol}
                           REL_TOL ${rel_tol}
                           DIR udq_actionx
                           TEST_ARGS --solver-max-time-step-in-days=10)
endif()

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

add_test_compareECLFiles(CASENAME micp
                         FILENAME MICP
                         SIMULATOR flow
                         ABS_TOL ${abs_tol}
                         REL_TOL ${rel_tol}
                         DIR micp
                         TEST_ARGS --linear-solver=ilu0)

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

add_test_compareECLFiles(CASENAME wvfpexp_02
                         FILENAME WVFPEXP-02
                         SIMULATOR flow
                         ABS_TOL ${abs_tol}
                         REL_TOL ${rel_tol}
                         DIR wvfpexp)

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
        "tests/definitions/regression/actionx.json"
        "tests/definitions/regression/aquifer_num.json"
        "tests/definitions/regression/aquifer-fetkovich.json"
        "tests/definitions/regression/aquifer-oilwater.json"
        "tests/definitions/regression/aquifers.json"
        "tests/definitions/regression/co2store.json"
        "tests/definitions/regression/cskin.json"
        "tests/definitions/regression/editnnc.json"
        "tests/definitions/regression/flowexp_blackoil.json"
        "tests/definitions/regression/gridunit.json"
        "tests/definitions/regression/gpmaint.json"
        "tests/definitions/regression/h2store.json"
        "tests/definitions/regression/jfunc.json"
        "tests/definitions/regression/krnum.json"
        "tests/definitions/regression/max_watercut.json"
        "tests/definitions/regression/model1.json"
        "tests/definitions/regression/model2.json"
        "tests/definitions/regression/model4.json"
        "tests/definitions/regression/model6.json"
        "tests/definitions/regression/msw_2d_h.json"
        "tests/definitions/regression/msw_3d_hfa.json"
        "tests/definitions/regression/mult.json"
        "tests/definitions/regression/network.json"
        "tests/definitions/regression/pinch.json"
        "tests/definitions/regression/polymer_injectivity.json"
        "tests/definitions/regression/polymer_oilwater.json"
        "tests/definitions/regression/polymer_simple2D.json"
        "tests/definitions/regression/ppcwmax.json"
        "tests/definitions/regression/radial_grid.json"
        "tests/definitions/regression/spe1.json"
        "tests/definitions/regression/spe1_brine.json"
        "tests/definitions/regression/spe1_precsalt.json"
        "tests/definitions/regression/spe1_solvent.json"
        "tests/definitions/regression/spe3.json"
        "tests/definitions/regression/spe5.json"
        "tests/definitions/regression/spe5_co2eor.json"
        "tests/definitions/regression/spe9.json"
        "tests/definitions/regression/spe9group.json"
        "tests/definitions/regression/tracer.json"
        "tests/definitions/regression/udq_actionx.json"
        "tests/definitions/regression/udt.json"
        "tests/definitions/regression/vappars.json"
        "tests/definitions/regression/vfpprod_spe1.json"
        "tests/definitions/regression/waghystr.json"
        "tests/definitions/regression/wecon_wtest.json"
        "tests/definitions/regression/wgrupcon.json"
        "tests/definitions/regression/winjdam.json"
        "tests/definitions/regression/winjmult.json"
        "tests/definitions/regression/wsegaicd.json"
        "tests/definitions/regression/wsegsicd.json"
        "tests/definitions/regression/wsegvalv.json"
        "tests/definitions/regression/wtest_bhp.json"
        "tests/definitions/regression/wtest_thp.json"
     )
  foreach(JSON_FILE ${JSON_FILES})
    AddJsonRegressionTests(${JSON_FILE})
    set_property(DIRECTORY APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS ${JSON_FILE})
  endforeach()
endif()
