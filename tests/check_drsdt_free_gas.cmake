# Runs a deck that is undersaturated everywhere and asserts that no free gas
# appears.  Invoked from CTest with -DFLOW=, -DDECK= and -DOUTDIR=.

execute_process(
  COMMAND ${CMAKE_COMMAND} -E rm -rf ${OUTDIR}
)
execute_process(
  COMMAND ${CMAKE_COMMAND} -E make_directory ${OUTDIR}
)
execute_process(
  COMMAND ${FLOW} ${DECK} --output-dir=${OUTDIR} --enable-vtk-output=false
  RESULT_VARIABLE run_status
  OUTPUT_QUIET
)
if(NOT run_status EQUAL 0)
  message(FATAL_ERROR "flow_blackoil failed with status ${run_status}")
endif()

file(GLOB prt "${OUTDIR}/*.PRT")
if(NOT prt)
  message(FATAL_ERROR "no PRT written to ${OUTDIR}")
endif()

file(STRINGS ${prt} balance REGEX "CURRENTLY IN PLACE")
if(NOT balance)
  message(FATAL_ERROR "no fluid-in-place report in ${prt}")
endif()

# ":CURRENTLY IN PLACE : <oil> : <water> : <free gas> <dissolved gas> <total> :"
# -- the free gas column must stay at zero for an undersaturated model.
foreach(line IN LISTS balance)
  string(REPLACE ":" ";" fields "${line}")
  list(GET fields 4 gas)
  string(STRIP "${gas}" gas)
  string(REGEX MATCH "^[0-9.eE+-]+" free_gas "${gas}")
  if(free_gas GREATER 1)
    message(FATAL_ERROR
      "free gas ${free_gas} in an undersaturated model; the DRSDT limiter "
      "drove the dissolved gas out of solution.\n  ${line}")
  endif()
endforeach()
