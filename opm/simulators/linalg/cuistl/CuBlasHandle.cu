#include <opm/simulators/linalg/cuistl/CuBlasHandle.hpp>
#include <opm/simulators/linalg/cuistl/cusparse_safe_call.hpp>
namespace Opm::cuistl
{


CuSparseHandle::CuSparseHandle()
{
    OPM_CUSPARSE_SAFE_CALL(cublasCreate(&handle));
}

CuSparseHandle::~CuSparseHandle()
{
    OPM_CUSPARSE_SAFE_CALL(cublasDestroy(handle));
}

cusparseHandle_t
CuSparseHandle::get()
{
    return handle;
}

CuSparseHandle&
CuSparseHandle::getInstance()
{
    static CuSparseHandle instance;
    return instance;
}

} // namespace Opm::cuistl
