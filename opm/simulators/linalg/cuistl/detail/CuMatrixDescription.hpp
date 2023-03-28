/*
  Copyright 2022-2023 SINTEF AS

  This file is part of the Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.
*/
#ifndef CU_MATRIX_DESCRIPTION_HPP
#define CU_MATRIX_DESCRIPTION_HPP
#include <opm/simulators/linalg/cuistl/detail/CuSparseResource.hpp>
#include <opm/simulators/linalg/cuistl/detail/cusparse_safe_call.hpp>

namespace Opm::cuistl::detail
{

/**
 * CuSparseMatrixDescription holder. This is internal information needed for most calls to the CuSparse API.
 */
using CuSparseMatrixDescription = CuSparseResource<cusparseMatDescr_t>;

/**
 * Pointer to CuSparseMatrixDescription holder. This is internal information needed for most calls to the CuSparse API.
 */
using CuSparseMatrixDescriptionPtr = std::shared_ptr<CuSparseResource<cusparseMatDescr_t>>;

/**
 * @brief createMatrixDescription creates a default matrix description
 * @return a matrix description to a general sparse matrix with zero based indexing.
 */
inline CuSparseMatrixDescriptionPtr
createMatrixDescription()
{
    auto description = std::make_shared<CuSparseMatrixDescription>();

    // Note: We always want to use zero base indexing.
    OPM_CUSPARSE_SAFE_CALL(cusparseSetMatType(description->get(), CUSPARSE_MATRIX_TYPE_GENERAL));
    OPM_CUSPARSE_SAFE_CALL(cusparseSetMatIndexBase(description->get(), CUSPARSE_INDEX_BASE_ZERO));

    return description;
}

/**
 * @brief createLowerDiagonalDescription creates a lower diagonal matrix description
 * @return a lower diagonal matrix description overlapped with options from ::Opm::cuistl::detail::createMatrixDescription()
 *
 * @note This will assume it has a unit diagonal
 */
inline CuSparseMatrixDescriptionPtr
createLowerDiagonalDescription()
{
    auto description = createMatrixDescription();
    OPM_CUSPARSE_SAFE_CALL(cusparseSetMatFillMode(description->get(), CUSPARSE_FILL_MODE_LOWER));
    OPM_CUSPARSE_SAFE_CALL(cusparseSetMatDiagType(description->get(), CUSPARSE_DIAG_TYPE_UNIT));
    return description;
}

/**
 * @brief createUpperDiagonalDescription creates an upper diagonal matrix description
 * @return an upper diagonal matrix description overlapped with options from ::Opm::cuistl::detail::createMatrixDescription()
 *
 * @note This will assume it has a non-unit diagonal.
 */
inline CuSparseMatrixDescriptionPtr
createUpperDiagonalDescription()
{
    auto description = createMatrixDescription();
    OPM_CUSPARSE_SAFE_CALL(cusparseSetMatFillMode(description->get(), CUSPARSE_FILL_MODE_UPPER));
    OPM_CUSPARSE_SAFE_CALL(cusparseSetMatDiagType(description->get(), CUSPARSE_DIAG_TYPE_NON_UNIT));

    return description;
}

} // namespace Opm::cuistl::detail

#endif // CU_MATRIX_DESCRIPTION_HPP
