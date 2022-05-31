#ifndef CU_MATRIX_DESCRIPTION_HPP
#define CU_MATRIX_DESCRIPTION_HPP
#include <opm/simulators/linalg/cuistl/CuSparseResource.hpp>
#include <opm/simulators/linalg/cuistl/cusparse_safe_call.hpp>

namespace Opm::cuistl {
using CuSparseMatrixDescription = CuSparseResource<cusparseMatDescr_t>;
using CuSparseMatrixDescriptionPtr = std::shared_ptr<CuSparseResource<cusparseMatDescr_t>>;
CuSparseMatrixDescriptionPtr createMatrixDescription() {
    auto description = std::make_shared<CuSparseMatrixDescription>();

    // Note: We always want to use zero base indexing.
    CUSPARSE_SAFE_CALL(cusparseSetMatIndexBase(description.get(), CUSPARSE_INDEX_BASE_ZERO));
    CUSPARSE_SAFE_CALL(cusparseSetMatType(description.get(), CUSPARSE_MATRIX_TYPE_GENERAL));

    return description;
}

CuSparseMatrixDescriptionPtr createLowerDiagonalDescription() {
    auto description = createMatrixDescription();
    CUSPARSE_SAFE_CALL(cusparseSetMatFillMode(description.get(), CUSPARSE_FILL_MODE_LOWER));
    CUSPARSE_SAFE_CALL(cusparseSetMatDiagType(description.get(), CUSPARSE_DIAG_TYPE_UNIT));
    return description;
}

CuSparseMatrixDescriptionPtr createUpperDiagonalDescription() {
    auto description = createMatrixDescription();
    CUSPARSE_SAFE_CALL(cusparseSetMatFillMode(description.get(), CUSPARSE_FILL_MODE_UPPER));
    CUSPARSE_SAFE_CALL(cusparseSetMatDiagType(description.get(), CUSPARSE_DIAG_TYPE_UNIT));

    return description;
}

}

#endif // CU_MATRIX_DESCRIPTION_HPP
