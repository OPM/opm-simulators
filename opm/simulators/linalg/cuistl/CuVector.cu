#include <cublas.h>
#include <cuda.h>
#include <opm/simulators/linalg/cuistl/CuVector.hpp>
#include <opm/simulators/linalg/cuistl/cublas_safe_call.hpp>
#include <opm/simulators/linalg/cuistl/cuda_safe_call.hpp>

namespace Opm::cuistl
{

template <class T>
CuVector<T>::CuVector(const int numberOfElements)
    : cuBlasHandle(CuBlasHandle::getInstance())
{
    CUDA_SAFE_CALL(cudaMalloc(&dataOnDevice, sizeof(T) * numberOfElements));
}

template <class T>
CuVector<T>::CuVector(const T* dataOnHost, const int numberOfElements)
    : CuVector(numberOfElements)
{
    CUDA_SAFE_CALL(cudaMemcpy(dataOnDevice, dataOnHost, cudaMemcpyHostToDevice));
}

template <class T>
CuVector<T>::~CuVector()
{
    cudaFree(data);
}

template <typename T>
const T*
CuVector<T>::data() const
{
    return dataOnDevice;
}

template <typename T>
T*
CuVector<T>::data()
{
    return dataOnDevice;
}

template <class T>
CuVector<T>&
CuVector<T>::operator*=(const T& scalar)
{
    OPM_CUBLAS_SAFE_CALL(cublasDscal(cuBlasHandle.get(), numberOfElements, scalar, data(), 1));
    return *this;
}
class CuVector<double>;
} // namespace Opm::cuistl
