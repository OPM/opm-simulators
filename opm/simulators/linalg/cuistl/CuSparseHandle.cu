#include <opm/simulators/linalg/cuistl/CuSparseHandle.hpp>
#include <opm/simulators/linalg/cuistl/cusparse_safe_call.hpp>
namespace Opm::cuistl {


CuSparseHandle::CuSparseHandle()
{
    OPM_CUSPARSE_SAFE_CALL(cusparseCreate(&handle));
}

CuSparseHandle::~CuSparseHandle()
{
    OPM_CUSPARSE_SAFE_CALL(cusparseDestroy(handle));
}

cusparseHandle_t CuSparseHandle::get()
{
    return handle;
}

}
