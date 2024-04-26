# Regression tests
opm_set_test_driver(${PROJECT_SOURCE_DIR}/tests/run-comparison.sh "")

# Set absolute tolerance to be used passed to the macros in the following tests
set(abs_tol 1e+5)
set(rel_tol 2e-5)
set(coarse_rel_tol 1e-2)
set(PYTHON_PATH ${PROJECT_BINARY_DIR}/python:${opm-common_DIR}/python:$ENV{PYTHONPATH})

add_test_compareSeparateECLFiles(CASENAME pyaction_gconprod_insert_kw
                                 DIR1 pyaction
                                 FILENAME1 PYACTION_GCONPROD_INSERT_KW
                                 DIR2 actionx
                                 FILENAME2 ACTIONX_GCONPROD
                                 SIMULATOR flow
                                 ABS_TOL ${abs_tol}
                                 REL_TOL ${rel_tol}
                                 IGNORE_EXTRA_KW BOTH
                                 ENVIRONMENT "PYTHONPATH=${PYTHON_PATH}")

add_test_compareSeparateECLFiles(CASENAME pyaction_mult+_insert_kw
                                 DIR1 pyaction
                                 FILENAME1 PYACTION_MULT+_INSERT_KW
                                 DIR2 actionx
                                 FILENAME2 ACTIONX_MULT+
                                 SIMULATOR flow
                                 ABS_TOL ${abs_tol}
                                 REL_TOL ${rel_tol}
                                 IGNORE_EXTRA_KW BOTH
                                 ENVIRONMENT "PYTHONPATH=${PYTHON_PATH}")

add_test_compareSeparateECLFiles(CASENAME pyaction_multx+_insert_kw
                                 DIR1 pyaction
                                 FILENAME1 PYACTION_MULTX+_INSERT_KW
                                 DIR2 actionx
                                 FILENAME2 ACTIONX_MULTX+
                                 SIMULATOR flow
                                 ABS_TOL ${abs_tol}
                                 REL_TOL ${rel_tol}
                                 IGNORE_EXTRA_KW BOTH
                                 ENVIRONMENT "PYTHONPATH=${PYTHON_PATH}")

add_test_compareSeparateECLFiles(CASENAME pyaction_multx-_insert_kw
                                 DIR1 pyaction
                                 FILENAME1 PYACTION_MULTX-_INSERT_KW
                                 DIR2 actionx
                                 FILENAME2 ACTIONX_MULTX-
                                 SIMULATOR flow
                                 ABS_TOL ${abs_tol}
                                 REL_TOL ${rel_tol}
                                 IGNORE_EXTRA_KW BOTH
                                 ENVIRONMENT "PYTHONPATH=${PYTHON_PATH}")

add_test_compareSeparateECLFiles(CASENAME pyaction_next_insert_kw
                                 DIR1 pyaction
                                 FILENAME1 PYACTION_NEXT_INSERT_KW
                                 DIR2 actionx
                                 FILENAME2 ACTIONX_NEXT
                                 SIMULATOR flow
                                 ABS_TOL ${abs_tol}
                                 REL_TOL ${rel_tol}
                                 IGNORE_EXTRA_KW BOTH
                                 ENVIRONMENT "PYTHONPATH=${PYTHON_PATH}")

add_test_compareSeparateECLFiles(CASENAME pyaction_wconprod_insert_kw
                                 DIR1 pyaction
                                 FILENAME1 PYACTION_WCONPROD_INSERT_KW
                                 DIR2 actionx
                                 FILENAME2 ACTIONX_WCONPROD
                                 SIMULATOR flow
                                 ABS_TOL ${abs_tol}
                                 REL_TOL ${rel_tol}
                                 IGNORE_EXTRA_KW BOTH
                                 ENVIRONMENT "PYTHONPATH=${PYTHON_PATH}")

add_test_compareSeparateECLFiles(CASENAME pyaction_wefac_insert_kw
                                 DIR1 pyaction
                                 FILENAME1 PYACTION_WEFAC_INSERT_KW
                                 DIR2 actionx
                                 FILENAME2 ACTIONX_WEFAC
                                 SIMULATOR flow
                                 ABS_TOL ${abs_tol}
                                 REL_TOL ${rel_tol}
                                 IGNORE_EXTRA_KW BOTH
                                 ENVIRONMENT "PYTHONPATH=${PYTHON_PATH}")
