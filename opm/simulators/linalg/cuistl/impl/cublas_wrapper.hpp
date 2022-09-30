#ifndef OPM_CUBLASWRAPPER_HEADER_INCLUDED
#define OPM_CUBLASWRAPPER_HEADER_INCLUDED
#include <cublas_v2.h>
#include <opm/common/ErrorMacros.hpp>

namespace Opm::cuistl::impl
{

inline cublasStatus_t
cublasScal(cublasHandle_t handle,
           int n,
           const double* alpha, /* host or device pointer */
           double* x,
           int incx)
{
    return cublasDscal(handle,
                       n,
                       alpha, /* host or device pointer */
                       x,
                       incx);
}

inline cublasStatus_t
cublasScal(cublasHandle_t handle,
           int n,
           const float* alpha, /* host or device pointer */
           float* x,
           int incx)
{
    return cublasSscal(handle,
                       n,
                       alpha, /* host or device pointer */
                       x,
                       incx);
}

inline cublasStatus_t
cublasScal(cublasHandle_t handle,
           int n,
           const int* alpha, /* host or device pointer */
           int* x,
           int incx)
{
    OPM_THROW(std::runtime_error, "cublasScal multiplication for integer vectors is not implemented yet.");
}
inline cublasStatus_t
cublasAxpy(cublasHandle_t handle,
           int n,
           const double* alpha, /* host or device pointer */
           const double* x,
           int incx,
           double* y,
           int incy)
{
    return cublasDaxpy(handle,
                       n,
                       alpha, /* host or device pointer */
                       x,
                       incx,
                       y,
                       incy);
}

inline cublasStatus_t
cublasAxpy(cublasHandle_t handle,
           int n,
           const float* alpha, /* host or device pointer */
           const float* x,
           int incx,
           float* y,
           int incy)
{
    return cublasSaxpy(handle,
                       n,
                       alpha, /* host or device pointer */
                       x,
                       incx,
                       y,
                       incy);
}

inline cublasStatus_t
cublasAxpy(cublasHandle_t handle,
           int n,
           const int* alpha, /* host or device pointer */
           const int* x,
           int incx,
           int* y,
           int incy)
{
    OPM_THROW(std::runtime_error, "axpy multiplication for integer vectors is not implemented yet.");
}

inline cublasStatus_t
cublasDot(cublasHandle_t handle, int n, const double* x, int incx, const double* y, int incy, double* result)
{
    return cublasDdot(handle, n, x, incx, y, incy, result);
}

inline cublasStatus_t
cublasDot(cublasHandle_t handle, int n, const float* x, int incx, const float* y, int incy, float* result)
{
    return cublasSdot(handle, n, x, incx, y, incy, result);
}

inline cublasStatus_t
cublasDot(cublasHandle_t handle, int n, const int* x, int incx, const int* y, int incy, int* result)
{
    OPM_THROW(std::runtime_error, "inner product for integer vectors is not implemented yet.");
}

inline cublasStatus_t
cublasNrm2(cublasHandle_t handle, int n, const double* x, int incx, double* result)
{
    return cublasDnrm2(handle, n, x, incx, result);
}


inline cublasStatus_t
cublasNrm2(cublasHandle_t handle, int n, const float* x, int incx, float* result)
{
    return cublasSnrm2(handle, n, x, incx, result);
}

inline cublasStatus_t
cublasNrm2(cublasHandle_t handle, int n, const int* x, int incx, int* result)
{
    OPM_THROW(std::runtime_error, "norm2 for integer vectors is not implemented yet.");
}

} // namespace Opm::cuistl::impl
#endif
