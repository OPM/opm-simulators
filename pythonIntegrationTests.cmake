file(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/tests/integration/HybridNewton)
execute_process(
  COMMAND "${Python3_EXECUTABLE}" -c "import tensorflow"
  RESULT_VARIABLE PYTHON_HAS_TENSORFLOW
  OUTPUT_QUIET
  ERROR_QUIET
)
if(PYTHON_HAS_TENSORFLOW EQUAL 0)
  add_test(NAME
            HybridNewton
           COMMAND
            ${Python3_EXECUTABLE}
            -m unittest discover
            ${PROJECT_SOURCE_DIR}/python/integration_tests/ml/hybrid_newton
          WORKING_DIRECTORY
            ${PROJECT_BINARY_DIR}/tests/integration/HybridNewton
    )
  if(CMAKE_VERSION VERSION_GREATER_EQUAL 3.22)
    set(modifications
        "PYTHONPATH=path_list_append:${opm-common_DIR}/python"
        "PYTHONPATH=path_list_append:${CMAKE_CURRENT_BINARY_DIR}/python"
        "FLOW_BINARY=set:$<TARGET_FILE:flow>"
    )
    set_tests_properties(HybridNewton PROPERTIES
                         ENVIRONMENT_MODIFICATION
                           "${modifications}"
    )
  endif()
endif()
