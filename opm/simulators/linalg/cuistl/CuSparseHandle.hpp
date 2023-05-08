#ifndef CUSPARSEHANDLE_HPP
#define CUSPARSEHANDLE_HPP
#include <cusparse.h>
#include <memory>

namespace Opm::cuistl
{
class CuSparseHandle
{
public:
    // This should not be copyable.
    CuSparseHandle(const CuSparseHandle&) = delete;
    CuSparseHandle& operator=(const CuSparseHandle&) = delete;

    ~CuSparseHandle();

    cusparseHandle_t get();

    static CuSparseHandle& getInstance();

private:
    CuSparseHandle();
    cusparseHandle_t handle;
};

using CuSparseHandlePtr = std::unique_ptr<CuSparseHandle>;
} // namespace Opm::cuistl
#endif // CUSPARSEHANDLE_HPP
