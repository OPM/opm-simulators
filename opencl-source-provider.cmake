set(BDA_DIR opm/simulators/linalg/bda)
set(KERNELS_DIR ${BDA_DIR}/opencl/kernels)

option(DEBUG_OPENCL_KERNELS_INTEL "Run ocloc to check kernel (works only on Intel)" OFF)
if(DEBUG_OPENCL_KERNELS_INTEL)
  set(DEBUG_OPENCL_DIR ${KERNELS_DIR}/.debug)

  execute_process(
    COMMAND ${CMAKE_COMMAND} -E make_directory ${DEBUG_OPENCL_DIR}
    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
  )
endif()

set(CL_SRC_FILE ${PROJECT_BINARY_DIR}/clSources.cpp)
file(WRITE ${CL_SRC_FILE} "// This file is auto-generated. Do not edit!\n\n")
file(APPEND ${CL_SRC_FILE} "#include <config.h>\n\n")
file(APPEND ${CL_SRC_FILE} "#include <${BDA_DIR}/opencl/openclKernels.hpp>\n\n")
file(APPEND ${CL_SRC_FILE} "namespace Opm\{\n\n")
file(APPEND ${CL_SRC_FILE} "namespace Accelerator\{\n\n")

file(GLOB CL_LIST "${KERNELS_DIR}/*.cl")

if(USE_CHOW_PATEL_ILU)
  list(REMOVE_ITEM CL_LIST "${PROJECT_SOURCE_DIR}/${KERNELS_DIR}/ILU_apply1_fm.cl")
  list(REMOVE_ITEM CL_LIST "${PROJECT_SOURCE_DIR}/${KERNELS_DIR}/ILU_apply2_fm.cl")
else()
  list(REMOVE_ITEM CL_LIST "${PROJECT_SOURCE_DIR}/${KERNELS_DIR}/ILU_apply1.cl")
  list(REMOVE_ITEM CL_LIST "${PROJECT_SOURCE_DIR}/${KERNELS_DIR}/ILU_apply2.cl")
endif()

foreach(CL ${CL_LIST})
  get_filename_component(FNAME ${CL} NAME_WE)

  file(APPEND ${CL_SRC_FILE} "const std::string OpenclKernels::${FNAME}_str = R\"\( \n")
  file(READ "${CL}" CL_CONTENT)
  file(APPEND ${CL_SRC_FILE} "${CL_CONTENT}")
  file(APPEND ${CL_SRC_FILE} "\)\"; \n\n")

  if(DEBUG_OPENCL_KERNELS_INTEL)
    execute_process(
      COMMAND ocloc -file ${CL} -device kbl -out_dir ${DEBUG_OPENCL_DIR}
      WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
      )
  endif()
endforeach()

file(APPEND ${CL_SRC_FILE} "\}\n")
file(APPEND ${CL_SRC_FILE} "\}")

if(DEBUG_OPENCL_KERNELS_INTEL)
  file(REMOVE_RECURSE ${DEBUG_DIR})
endif()
