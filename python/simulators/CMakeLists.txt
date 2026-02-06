set(PYTHON_OPM_SIMULATORS_PACKAGE_PATH ${PROJECT_BINARY_DIR}/python/opm/simulators)

# Set the path to the input docstrings.json file and the output .hpp files
set(PYTHON_DOCSTRINGS_FILE "${PROJECT_SOURCE_DIR}/python/docstrings_simulators.json")
set(PYTHON_DOCSTRINGS_GENERATED_HPP "${PROJECT_BINARY_DIR}/python/PyBlackOilSimulatorDoc.hpp")
set(PYTHON_GW_DOCSTRINGS_GENERATED_HPP "${PROJECT_BINARY_DIR}/python/PyGasWaterSimulatorDoc.hpp")
set(PYTHON_OP_DOCSTRINGS_GENERATED_HPP "${PROJECT_BINARY_DIR}/python/PyOnePhaseSimulatorDoc.hpp")

# Command to run the Python script

find_file(PYTHON_GENERATE_DOCSTRINGS_PY generate_docstring_hpp.py
  PATHS ${opm-common_DIR} ${opm-common_PYTHON_COMMON_DIR}
  PATH_SUFFIXES python NO_DEFAULT_PATH REQUIRED)

add_custom_command(
    OUTPUT
      ${PYTHON_DOCSTRINGS_GENERATED_HPP}
    COMMAND
      ${CMAKE_COMMAND} -E env PYTHONPATH=${CMAKE_SOURCE_DIR}
      $<TARGET_FILE:Python3::Interpreter>
      ${PYTHON_GENERATE_DOCSTRINGS_PY}
      ${PYTHON_DOCSTRINGS_FILE}
      ${PYTHON_DOCSTRINGS_GENERATED_HPP}
      PYBLACKOILSIMULATORDOC_HPP
      "Opm::Pybind::DocStrings"
      "BlackOil"
    DEPENDS
      ${PYTHON_DOCSTRINGS_FILE}
    COMMENT
      "Generating PyBlackOilSimulatorDoc.hpp from JSON file"
)
add_custom_command(
    OUTPUT
      ${PYTHON_GW_DOCSTRINGS_GENERATED_HPP}
    COMMAND
      ${CMAKE_COMMAND} -E env PYTHONPATH=${CMAKE_SOURCE_DIR}
      $<TARGET_FILE:Python3::Interpreter>
      ${PYTHON_GENERATE_DOCSTRINGS_PY}
      ${PYTHON_DOCSTRINGS_FILE}
      ${PYTHON_GW_DOCSTRINGS_GENERATED_HPP}
      PYGASWATERSIMULATORDOC_HPP
      "Opm::Pybind::DocStrings"
      "GasWater"
    DEPENDS
      ${PYTHON_DOCSTRINGS_FILE}
    COMMENT
      "Generating PyGasWaterSimulatorDoc.hpp from JSON file"
)
add_custom_command(
    OUTPUT
      ${PYTHON_OP_DOCSTRINGS_GENERATED_HPP}
    COMMAND
      ${CMAKE_COMMAND} -E env PYTHONPATH=${CMAKE_SOURCE_DIR}
      $<TARGET_FILE:Python3::Interpreter>
      ${PYTHON_GENERATE_DOCSTRINGS_PY}
      ${PYTHON_DOCSTRINGS_FILE}
      ${PYTHON_OP_DOCSTRINGS_GENERATED_HPP}
      PYONEPHASESIMULATORDOC_HPP
      "Opm::Pybind::DocStrings"
      "OnePhase"
    DEPENDS
      ${PYTHON_DOCSTRINGS_FILE}
    COMMENT
      "Generating PyOnePhaseSimulatorDoc.hpp from JSON file"
)
# NOTE: The variable ${PYBIND11_SYSTEM} is set in python/CMakeLists.txt
#   to the value "SYSTEM" or unset, depending on the current version of Pybind11.
#   The value is then forwarded to target_include_directories(), see
#
#  https://cmake.org/cmake/help/latest/command/target_include_directories.html
#  https://pybind11.readthedocs.io/en/stable/compiling.html
#
pybind11_add_module(BlackOil ${PYBIND11_SYSTEM}
  $<TARGET_OBJECTS:moduleVersion>
  PyBlackOilSimulator.cpp
  ${PYTHON_DOCSTRINGS_GENERATED_HPP}  # Include the generated .hpp as a source file
)
opm_add_target_options(TARGET BlackOil)

pybind11_add_module(GasWater ${PYBIND11_SYSTEM}
  $<TARGET_OBJECTS:moduleVersion>
  PyGasWaterSimulator.cpp
  ${PYTHON_GW_DOCSTRINGS_GENERATED_HPP}  # Include the generated .hpp as a source file
)
opm_add_target_options(TARGET GasWater)

pybind11_add_module(OnePhase ${PYBIND11_SYSTEM}
  $<TARGET_OBJECTS:moduleVersion>
  PyOnePhaseSimulator.cpp
  ${PYTHON_OP_DOCSTRINGS_GENERATED_HPP}  # Include the generated .hpp as a source file
)
opm_add_target_options(TARGET OnePhase)

# Create a convenience target to build all Python simulator modules
add_custom_target(python_simulator_modules
  DEPENDS BlackOil GasWater OnePhase
  COMMENT "Building all Python simulator modules (BlackOil, GasWater, etc.)"
)

# Create a convenience target to build everything needed for Python tests
add_custom_target(python_tests
  DEPENDS python_simulator_modules copy_python
  COMMENT "Building Python test dependencies (modules and test files)"
)

find_file(PYTHON_INSTALL_PY install.py
  PATHS ${opm-common_DIR} ${opm-common_PYTHON_COMMON_DIR}
  PATH_SUFFIXES python NO_DEFAULT_PATH REQUIRED)

# NOTE: instead of using file( COPY ...) which copies the files at configure time (not at build time)
#  we should add copying of the files at build time such that running "make" or "make all" will
#  update the files if they have been modified. We use the install.py script in opm-common, see also
#  CMakeLists.txt in opm-common
add_custom_target(copy_python ALL
    COMMAND
      $<TARGET_FILE:Python3::Interpreter>
      "${PYTHON_INSTALL_PY}"
      ${PROJECT_SOURCE_DIR}/python/opm
      ${PROJECT_BINARY_DIR}/python
      0
    COMMAND
      $<TARGET_FILE:Python3::Interpreter>
      "${PYTHON_INSTALL_PY}"
      ${PROJECT_SOURCE_DIR}/python/test
      ${PROJECT_BINARY_DIR}/python
      0
    COMMAND
      $<TARGET_FILE:Python3::Interpreter>
      "${PYTHON_INSTALL_PY}"
      ${PROJECT_SOURCE_DIR}/python/test_data
      ${PROJECT_BINARY_DIR}/python
      0
)

foreach(target_name IN ITEMS BlackOil GasWater OnePhase)
  set_target_properties(
    ${target_name} PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${PYTHON_OPM_SIMULATORS_PACKAGE_PATH}
  )
  target_link_libraries( ${target_name} PRIVATE opmsimulators )
  # Add the binary (build) directory to the include directories for the target

    # Add the build directory where the generated hpp file will be
  #  to the include directories for the target
  target_include_directories(${target_name} PRIVATE ${PROJECT_BINARY_DIR}/python)
  add_dependencies(${target_name} copy_python)

  # Since the installation of Python code is nonstandard it is protected by an
  # extra cmake switch, OPM_INSTALL_PYTHON. If you prefer you can still invoke
  # setup.py install manually - optionally with the generated script
  # setup-install.sh - and completely bypass cmake in the installation phase.
  if (OPM_INSTALL_PYTHON)
    include(PyInstallPrefix)  # from opm-common
    install(
      TARGETS
        ${target_name}
      DESTINATION
        ${DEST_PREFIX}${CMAKE_INSTALL_PREFIX}/${PYTHON_INSTALL_PREFIX}/opm
    )
    install(
      CODE "execute_process(
            COMMAND
              $<TARGET_FILE:Python3::Interpreter>
              ${PYTHON_INSTALL_PY}
              ${PROJECT_BINARY_DIR}/python/opm
              ${DEST_PREFIX}${CMAKE_INSTALL_PREFIX}/${PYTHON_INSTALL_PREFIX}
              1
            )"
    )
  endif()
endforeach()

include(PyInstallPrefix)
if(OPM_ENABLE_PYTHON_TESTS)
  set(cases basic fluidstate_variables primary_variables schedule throw)
  if(MPI_FOUND)
    list(APPEND cases mpi)
  endif()

  # Create a fixture test that builds Python modules before running Python tests
  add_test(
    NAME
      python_modules_build
    COMMAND
      ${CMAKE_COMMAND}
      --build ${CMAKE_BINARY_DIR}
      --target python_tests
  )
  set_tests_properties(python_modules_build
    PROPERTIES
    FIXTURES_SETUP
      python_modules_fixture
  )

  # NOTE: See comment in test_basic.py for the reason why we are
  #   splitting the python tests into multiple add_test() tests instead
  #   of having a single "python -m unittest" test call that will run all
  #   the tests in the "test" sub directory.
  foreach(case_name IN ITEMS ${cases})
    add_test(
      NAME
        python_${case_name}
      WORKING_DIRECTORY
        ${PROJECT_BINARY_DIR}/python
      COMMAND
        $<TARGET_FILE:Python3::Interpreter>
        -m unittest
        test/test_${case_name}.py
    )
    set(env_mod
      PYTHONPATH=path_list_prepend:${PROJECT_BINARY_DIR}/python
      PYTHONPATH=path_list_prepend:${opm-common_DIR}/python
    )
    set_tests_properties(python_${case_name}
      PROPERTIES
      ENVIRONMENT_MODIFICATION
        "${env_mod}"
      FIXTURES_REQUIRED
        python_modules_fixture
    )
  endforeach()
endif()

file(
  COPY
    ${PROJECT_SOURCE_DIR}/python/opm/CMakeLists.txt
    ${PROJECT_SOURCE_DIR}/python/requirements.txt
    ${PROJECT_SOURCE_DIR}/python/README.md
    ${PROJECT_SOURCE_DIR}/python/MANIFEST.in
    ${PROJECT_SOURCE_DIR}/python/pyproject.toml
  DESTINATION
    ${PROJECT_BINARY_DIR}/python
)

set(opm-simulators_PYTHON_PACKAGE_VERSION ${OPM_PYTHON_PACKAGE_VERSION_TAG})

# Generate versioned setup.py
configure_file(${PROJECT_SOURCE_DIR}/python/setup.py.in
                ${PROJECT_BINARY_DIR}/python/setup.py.tmp @ONLY)
file(GENERATE
  OUTPUT
    ${PROJECT_BINARY_DIR}/python/setup.py
  INPUT
    ${PROJECT_BINARY_DIR}/python/setup.py.tmp
)
