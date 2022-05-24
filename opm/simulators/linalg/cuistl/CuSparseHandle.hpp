#ifndef CUSPARSEHANDLE_HPP
#define CUSPARSEHANDLE_HPP
#include <cusparse.h>
#include <memory>

namespace Opm::cuistl {
class CuSparseHandle
{
public:
    CuSparseHandle();

    ~CuSparseHandle();

    cusparseHandle_t get();
    const cusparseHandle_t get() const;
private:
    cusparseHandle_t handle;
};

using CuSparseHandlePtr = std::unique_ptr<CuSparseHandle>;
}
#endif // CUSPARSEHANDLE_HPP
