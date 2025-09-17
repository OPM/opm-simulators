#include <hip/hip_runtime.h>

#include <opm/models/discretization/common/tpfalinearizer_gpu.hh>

#include <opm/common/utility/gpuDecorators.hpp>

#include <opm/simulators/linalg/gpuistl_hip/MiniMatrix.hpp>
#include <opm/simulators/linalg/gpuistl_hip/MiniVector.hpp>
#include <opm/simulators/linalg/gpuistl_hip/GpuVector.hpp>
#include <opm/simulators/linalg/gpuistl_hip/GpuSparseMatrix.hpp>

__global__ void testKernel() {
    printf("Hello from GPU!\n");
}

void kernelWrapper() {
    testKernel<<<1, 1>>>();
    hipDeviceSynchronize();
}
