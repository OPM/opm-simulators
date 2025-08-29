# Regression tests
opm_set_test_driver(${PROJECT_SOURCE_DIR}/tests/run-comparison.sh "")

# Set absolute tolerance to be used passed to the macros in the following tests
set(abs_tol 1e+5)
set(rel_tol 2e-5)
set(coarse_rel_tol 1e-2)
add_test_compareSeparateECLFiles(CASENAME pyaction_gconprod_insert_kw
                                 DIR1 pyaction
                                 FILENAME1 PYACTION_GCONPROD_INSERT_KW
                                 DIR2 actionx
                                 FILENAME2 ACTIONX_GCONPROD
                                 SIMULATOR flow
                                 ABS_TOL ${abs_tol}
                                 REL_TOL ${rel_tol}
                                 IGNORE_EXTRA_KW BOTH)

add_test_compareSeparateECLFiles(CASENAME pyaction_gconsump_insert_kw
                                 DIR1 pyaction
                                 FILENAME1 PYACTION_GCONSUMP_INSERT_KW
                                 DIR2 actionx
                                 FILENAME2 ACTIONX_GCONSUMP
                                 SIMULATOR flow
                                 ABS_TOL ${abs_tol}
                                 REL_TOL ${rel_tol}
                                 IGNORE_EXTRA_KW BOTH)

add_test_compareSeparateECLFiles(CASENAME actionx_gefac
                                 DIR1 model4
                                 FILENAME1 MOD4_GRP_GEFAC
                                 DIR2 actionx
                                 FILENAME2 ACTIONX_GEFAC
                                 SIMULATOR flow
                                 ABS_TOL ${abs_tol}
                                 REL_TOL ${rel_tol}
                                 IGNORE_EXTRA_KW BOTH)

add_test_compareSeparateECLFiles(CASENAME pyaction_gefac_insert_kw
                                 DIR1 pyaction
                                 FILENAME1 PYACTION_GEFAC_INSERT_KW
                                 DIR2 actionx
                                 FILENAME2 ACTIONX_GEFAC
                                 SIMULATOR flow
                                 ABS_TOL ${abs_tol}
                                 REL_TOL ${rel_tol}
                                 IGNORE_EXTRA_KW BOTH)

add_test_compareSeparateECLFiles(CASENAME pyaction_gruptree_insert_kw
                                 DIR1 pyaction
                                 FILENAME1 PYACTION_GRUPTREE_INSERT_KW
                                 DIR2 actionx
                                 FILENAME2 ACTIONX_GRUPTREE
                                 SIMULATOR flow
                                 ABS_TOL ${abs_tol}
                                 REL_TOL ${rel_tol}
                                 IGNORE_EXTRA_KW BOTH)

add_test_compareSeparateECLFiles(CASENAME pyaction_include_insert_kw
                                 DIR1 pyaction
                                 FILENAME1 PYACTION_INCLUDE_INSERT_KW
                                 DIR2 actionx
                                 FILENAME2 ACTIONX_INCLUDE
                                 SIMULATOR flow
                                 ABS_TOL ${abs_tol}
                                 REL_TOL ${rel_tol}
                                 IGNORE_EXTRA_KW BOTH)

add_test_compareSeparateECLFiles(CASENAME pyaction_mult+_insert_kw
                                 DIR1 pyaction
                                 FILENAME1 PYACTION_MULT+_INSERT_KW
                                 DIR2 actionx
                                 FILENAME2 ACTIONX_MULT+
                                 SIMULATOR flow
                                 ABS_TOL ${abs_tol}
                                 REL_TOL ${rel_tol}
                                 IGNORE_EXTRA_KW BOTH)

add_test_compareSeparateECLFiles(CASENAME pyaction_multx+_insert_kw
                                 DIR1 pyaction
                                 FILENAME1 PYACTION_MULTX+_INSERT_KW
                                 DIR2 actionx
                                 FILENAME2 ACTIONX_MULTX+
                                 SIMULATOR flow
                                 ABS_TOL ${abs_tol}
                                 REL_TOL ${rel_tol}
                                 IGNORE_EXTRA_KW BOTH)

add_test_compareSeparateECLFiles(CASENAME pyaction_multx-_insert_kw
                                 DIR1 pyaction
                                 FILENAME1 PYACTION_MULTX-_INSERT_KW
                                 DIR2 actionx
                                 FILENAME2 ACTIONX_MULTX-
                                 SIMULATOR flow
                                 ABS_TOL ${abs_tol}
                                 REL_TOL ${rel_tol}
                                 IGNORE_EXTRA_KW BOTH)

add_test_compareSeparateECLFiles(CASENAME pyaction_next_insert_kw
                                 DIR1 pyaction
                                 FILENAME1 PYACTION_NEXT_INSERT_KW
                                 DIR2 actionx
                                 FILENAME2 ACTIONX_NEXT
                                 SIMULATOR flow
                                 ABS_TOL ${abs_tol}
                                 REL_TOL ${rel_tol}
                                 IGNORE_EXTRA_KW BOTH)

add_test_compareSeparateECLFiles(CASENAME pyaction_wconprod_insert_kw
                                 DIR1 pyaction
                                 FILENAME1 PYACTION_WCONPROD_INSERT_KW
                                 DIR2 actionx
                                 FILENAME2 ACTIONX_WCONPROD
                                 SIMULATOR flow
                                 ABS_TOL ${abs_tol}
                                 REL_TOL ${rel_tol}
                                 IGNORE_EXTRA_KW BOTH)

add_test_compareSeparateECLFiles(CASENAME pyaction_wefac_insert_kw
                                 DIR1 pyaction
                                 FILENAME1 PYACTION_WEFAC_INSERT_KW
                                 DIR2 actionx
                                 FILENAME2 ACTIONX_WEFAC
                                 SIMULATOR flow
                                 ABS_TOL ${abs_tol}
                                 REL_TOL ${rel_tol}
                                 IGNORE_EXTRA_KW BOTH)

add_test_compareSeparateECLFiles(CASENAME pyaction_welpi_insert_kw
                                 DIR1 pyaction
                                 FILENAME1 PYACTION_WELPI_INSERT_KW
                                 DIR2 actionx
                                 FILENAME2 ACTIONX_WELPI
                                 SIMULATOR flow
                                 ABS_TOL ${abs_tol}
                                 REL_TOL ${rel_tol}
                                 IGNORE_EXTRA_KW BOTH)

add_test_compareSeparateECLFiles(CASENAME pyaction_wpimult_insert_kw
                                 DIR1 pyaction
                                 FILENAME1 PYACTION_WPIMULT_INSERT_KW
                                 DIR2 actionx
                                 FILENAME2 ACTIONX_WPIMULT
                                 SIMULATOR flow
                                 ABS_TOL ${abs_tol}
                                 REL_TOL ${rel_tol}
                                 IGNORE_EXTRA_KW BOTH)

add_test_compareSeparateECLFiles(CASENAME pyaction_wsegvalv_insert_kw
                                 DIR1 pyaction
                                 FILENAME1 PYACTION_WSEGVALV_INSERT_KW
                                 DIR2 actionx
                                 FILENAME2 ACTIONX_WSEGVALV
                                 SIMULATOR flow
                                 ABS_TOL ${abs_tol}
                                 REL_TOL ${rel_tol}
                                 IGNORE_EXTRA_KW BOTH)

add_test_compareSeparateECLFiles(CASENAME actionx_wlist
                                 DIR1 udq_actionx
                                 FILENAME1 ACTIONX_M1
                                 DIR2 actionx
                                 FILENAME2 ACTIONX_WLIST
                                 SIMULATOR flow
                                 ABS_TOL ${abs_tol}
                                 REL_TOL ${rel_tol}
                                 IGNORE_EXTRA_KW BOTH)

add_test_compareSeparateECLFiles(CASENAME pyaction_wlist_insert_kw
                                 DIR1 pyaction
                                 FILENAME1 PYACTION_WLIST_INSERT_KW
                                 DIR2 actionx
                                 FILENAME2 ACTIONX_WLIST
                                 SIMULATOR flow
                                 ABS_TOL ${abs_tol}
                                 REL_TOL ${rel_tol}
                                 IGNORE_EXTRA_KW BOTH)

add_test_compareSeparateECLFiles(CASENAME pyaction_wtest_insert_kw
                                 DIR1 pyaction
                                 FILENAME1 PYACTION_WTEST_INSERT_KW
                                 DIR2 actionx
                                 FILENAME2 ACTIONX_WTEST
                                 SIMULATOR flow
                                 ABS_TOL ${abs_tol}
                                 REL_TOL ${rel_tol}
                                 IGNORE_EXTRA_KW BOTH)

add_test_compareSeparateECLFiles(CASENAME pyaction_WTMULT_insert_kw
                                 DIR1 pyaction
                                 FILENAME1 PYACTION_WTMULT_INSERT_KW
                                 DIR2 actionx
                                 FILENAME2 ACTIONX_WTMULT
                                 SIMULATOR flow
                                 ABS_TOL ${abs_tol}
                                 REL_TOL ${rel_tol}
                                 IGNORE_EXTRA_KW BOTH)

if(MPI_FOUND)
  add_test_compareSeparateECLFiles(CASENAME pyaction_gconsump_insert_kw_4_procs
                                   DIR1 pyaction
                                   FILENAME1 PYACTION_GCONSUMP_INSERT_KW
                                   DIR2 actionx
                                   FILENAME2 ACTIONX_GCONSUMP
                                   SIMULATOR flow
                                   ABS_TOL ${abs_tol}
                                   REL_TOL ${rel_tol}
                                   IGNORE_EXTRA_KW BOTH
                                   MPI_PROCS 4)


  add_test_compareSeparateECLFiles(CASENAME actionx_gefac_4_procs
                                   DIR1 model4
                                   FILENAME1 MOD4_GRP_GEFAC
                                   DIR2 actionx
                                   FILENAME2 ACTIONX_GEFAC
                                   SIMULATOR flow
                                   ABS_TOL ${abs_tol}
                                   REL_TOL ${rel_tol}
                                   IGNORE_EXTRA_KW BOTH
                                   MPI_PROCS 4)

  add_test_compareSeparateECLFiles(CASENAME pyaction_gefac_insert_kw_4_procs
                                   DIR1 pyaction
                                   FILENAME1 PYACTION_GEFAC_INSERT_KW
                                   DIR2 actionx
                                   FILENAME2 ACTIONX_GEFAC
                                   SIMULATOR flow
                                   ABS_TOL ${abs_tol}
                                   REL_TOL ${rel_tol}
                                   IGNORE_EXTRA_KW BOTH
                                   MPI_PROCS 4)

  add_test_compareSeparateECLFiles(CASENAME pyaction_multx+_insert_kw_4_procs
                                   DIR1 pyaction
                                   FILENAME1 PYACTION_MULTX+_INSERT_KW
                                   DIR2 actionx
                                   FILENAME2 ACTIONX_MULTX+
                                   SIMULATOR flow
                                   ABS_TOL ${abs_tol}
                                   REL_TOL ${rel_tol}
                                   IGNORE_EXTRA_KW BOTH
                                   MPI_PROCS 4)

  add_test_compareSeparateECLFiles(CASENAME pyaction_mult+_insert_kw_4_procs
                                   DIR1 pyaction
                                   FILENAME1 PYACTION_MULT+_INSERT_KW
                                   DIR2 actionx
                                   FILENAME2 ACTIONX_MULT+
                                   SIMULATOR flow
                                   ABS_TOL ${abs_tol}
                                   REL_TOL ${rel_tol}
                                   IGNORE_EXTRA_KW BOTH
                                   MPI_PROCS 4)

  add_test_compareSeparateECLFiles(CASENAME pyaction_multx-_insert_kw_4_procs
                                   DIR1 pyaction
                                   FILENAME1 PYACTION_MULTX-_INSERT_KW
                                   DIR2 actionx
                                   FILENAME2 ACTIONX_MULTX-
                                   SIMULATOR flow
                                   ABS_TOL ${abs_tol}
                                   REL_TOL ${rel_tol}
                                   IGNORE_EXTRA_KW BOTH
                                   MPI_PROCS 4)

  add_test_compareSeparateECLFiles(CASENAME pyaction_next_insert_kw_4_procs
                                   DIR1 pyaction
                                   FILENAME1 PYACTION_NEXT_INSERT_KW
                                   DIR2 actionx
                                   FILENAME2 ACTIONX_NEXT
                                   SIMULATOR flow
                                   ABS_TOL ${abs_tol}
                                   REL_TOL ${rel_tol}
                                   IGNORE_EXTRA_KW BOTH
                                   MPI_PROCS 4)

  add_test_compareSeparateECLFiles(CASENAME pyaction_wefac_insert_kw_4_procs
                                   DIR1 pyaction
                                   FILENAME1 PYACTION_WEFAC_INSERT_KW
                                   DIR2 actionx
                                   FILENAME2 ACTIONX_WEFAC
                                   SIMULATOR flow
                                   ABS_TOL ${abs_tol}
                                   REL_TOL ${rel_tol}
                                   IGNORE_EXTRA_KW BOTH
                                   MPI_PROCS 4)

  add_test_compareSeparateECLFiles(CASENAME pyaction_welpi_insert_kw_4_procs
                                 DIR1 pyaction
                                 FILENAME1 PYACTION_WELPI_INSERT_KW
                                 DIR2 actionx
                                 FILENAME2 ACTIONX_WELPI
                                 SIMULATOR flow
                                 ABS_TOL ${abs_tol}
                                 REL_TOL ${rel_tol}
                                 IGNORE_EXTRA_KW BOTH,
                                 MPI_PROCS 4)

  add_test_compareSeparateECLFiles(CASENAME pyaction_wpimult_insert_kw_4_procs
                                 DIR1 pyaction
                                 FILENAME1 PYACTION_WPIMULT_INSERT_KW
                                 DIR2 actionx
                                 FILENAME2 ACTIONX_WPIMULT
                                 SIMULATOR flow
                                 ABS_TOL ${abs_tol}
                                 REL_TOL ${rel_tol}
                                 IGNORE_EXTRA_KW BOTH
                                 MPI_PROCS 4)

  add_test_compareSeparateECLFiles(CASENAME pyaction_wsegvalv_insert_kw_4_procs
                                   DIR1 pyaction
                                   FILENAME1 PYACTION_WSEGVALV_INSERT_KW
                                   DIR2 actionx
                                   FILENAME2 ACTIONX_WSEGVALV
                                   SIMULATOR flow
                                   ABS_TOL ${abs_tol}
                                   REL_TOL ${rel_tol}
                                   IGNORE_EXTRA_KW BOTH
                                   MPI_PROCS 4)

  add_test_compareSeparateECLFiles(CASENAME pyaction_WTMULT_insert_kw_4_procs
                                   DIR1 pyaction
                                   FILENAME1 PYACTION_WTMULT_INSERT_KW
                                   DIR2 actionx
                                   FILENAME2 ACTIONX_WTMULT
                                   SIMULATOR flow
                                   ABS_TOL ${abs_tol}
                                   REL_TOL ${rel_tol}
                                   IGNORE_EXTRA_KW BOTH
                                   MPI_PROCS 4)

  add_test_compareSeparateECLFiles(CASENAME actionx_wlist_4_procs
                                   DIR1 udq_actionx
                                   FILENAME1 ACTIONX_M1
                                   DIR2 actionx
                                   FILENAME2 ACTIONX_WLIST
                                   SIMULATOR flow
                                   ABS_TOL ${abs_tol}
                                   REL_TOL ${rel_tol}
                                   IGNORE_EXTRA_KW BOTH
                                   MPI_PROCS 4)

  add_test_compareSeparateECLFiles(CASENAME pyaction_wlist_insert_kw_4_procs
                                   DIR1 pyaction
                                   FILENAME1 PYACTION_WLIST_INSERT_KW
                                   DIR2 actionx
                                   FILENAME2 ACTIONX_WLIST
                                   SIMULATOR flow
                                   ABS_TOL ${abs_tol}
                                   REL_TOL ${rel_tol}
                                   IGNORE_EXTRA_KW BOTH
                                   MPI_PROCS 4)

  add_test_compareSeparateECLFiles(CASENAME pyaction_wtest_insert_kw_4_procs
                                   DIR1 pyaction
                                   FILENAME1 PYACTION_WTEST_INSERT_KW
                                   DIR2 actionx
                                   FILENAME2 ACTIONX_WTEST
                                   SIMULATOR flow
                                   ABS_TOL ${abs_tol}
                                   REL_TOL ${rel_tol}
                                   IGNORE_EXTRA_KW BOTH
                                   MPI_PROCS 4)
endif()
