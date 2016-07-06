# Set absolute tolerance to be used for testing
set(abs_tol 2e-2)
set(rel_tol 1e-5)
# Define some paths
set(RESULT_PATH ${PROJECT_BINARY_DIR}/tests/results)
# Create directory to store upscaling results in

###########################################################################
# TEST: compareECLFiles 
###########################################################################

# Input: 
#   - casename: basename (no extension) 
#
macro (add_test_compareECLFiles casename filename)
  
  # Add test that runs flow and outputs the results to file
  opm_add_test(compareECLFiles_${filename} NO_COMPILE
               EXE_NAME flow    
	       DRIVER_ARGS ${OPM_DATA_ROOT}/${casename} ${RESULT_PATH}
                           ${CMAKE_BINARY_DIR}/bin
                           ${filename}
                           ${abs_tol} ${rel_tol}
			   ${SUMMARY_REGRESSION_TEST_COMMAND} 
			   ${RESTART_REGRESSION_TEST_COMMAND} 
			   ${INIT_REGRESSION_TEST_COMMAND} 
               TEST_ARGS ${OPM_DATA_ROOT}/${casename}/${filename}.DATA )
endmacro (add_test_compareECLFiles)


if(NOT TARGET test-suite)
  add_custom_target(test-suite)
endif()


opm_set_test_driver(${PROJECT_SOURCE_DIR}/tests/run-regressionTest.sh "")

add_test_compareECLFiles(spe9 SPE9_CP)

