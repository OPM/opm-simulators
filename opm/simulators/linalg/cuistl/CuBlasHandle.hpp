#ifndef CUBLASHANDLE_HPP
#define CUBLASHANDLE_HPP
#include <cublas_v2.h>
#include <memory>

namespace Opm::cuistl
{
class CuBlasHandle
{
public:
    // This should not be copyable.
    CuBlasHandle(const CuBlasHandle&) = delete;
    CuBlasHandle& operator=(const CuBlasHandle&) = delete;

    ~CuBlasHandle();

    cublasHandle_t get();

    static CuBlasHandle& getInstance();

private:
    CuBlasHandle();
    cublasHandle_t handle;
};

using CuBlasHandlePtr = std::unique_ptr<CuBlasHandle>;
} // namespace Opm::cuistl
#endif // CuBlasHandle_HPP
