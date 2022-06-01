#ifndef CU_MATRIX_DESCRIPTION_HPP
#define CU_MATRIX_DESCRIPTION_HPP
#include <opm/simulators/linalg/cuistl/CuSparseResource.hpp>
#include <opm/simulators/linalg/cuistl/cusparse_safe_call.hpp>

namespace Opm::cuistl
{
using CuSparseMatrixDescription = CuSparseResource<cusparseMatDescr_t>;
using CuSparseMatrixDescriptionPtr = std::shared_ptr<CuSparseResource<cusparseMatDescr_t>>;
CuSparseMatrixDescriptionPtr
createMatrixDescription()
{
    auto description = std::make_shared<CuSparseMatrixDescription>();

    // Note: We always want to use zero base indexing.
    OPM_CUSPARSE_SAFE_CALL(cusparseSetMatIndexBase(description->get(), CUSPARSE_INDEX_BASE_ZERO));
    OPM_CUSPARSE_SAFE_CALL(cusparseSetMatType(description->get(), CUSPARSE_MATRIX_TYPE_GENERAL));

    return description;
}

CuSparseMatrixDescriptionPtr
createLowerDiagonalDescription()
{
    auto description = createMatrixDescription();
    OPM_CUSPARSE_SAFE_CALL(cusparseSetMatFillMode(description->get(), CUSPARSE_FILL_MODE_LOWER));
    OPM_CUSPARSE_SAFE_CALL(cusparseSetMatDiagType(description->get(), CUSPARSE_DIAG_TYPE_UNIT));
    return description;
}

CuSparseMatrixDescriptionPtr
createUpperDiagonalDescription()
{
    auto description = createMatrixDescription();
    OPM_CUSPARSE_SAFE_CALL(cusparseSetMatFillMode(description->get(), CUSPARSE_FILL_MODE_UPPER));
    OPM_CUSPARSE_SAFE_CALL(cusparseSetMatDiagType(description->get(), CUSPARSE_DIAG_TYPE_UNIT));

    return description;
}

} // namespace Opm::cuistl

#endif // CU_MATRIX_DESCRIPTION_HPP
