/*
  Copyright 2025 Equinor ASA

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

#ifndef OPM_HYPRE_MATRIX_UTILS_HPP
#define OPM_HYPRE_MATRIX_UTILS_HPP

// Fix conflict with Dune's matrixmarket.hh header - HYPRE defines MM_MAX_LINE_LENGTH as macro
// but Dune expects it as enum value
#ifdef MM_MAX_LINE_LENGTH
#undef MM_MAX_LINE_LENGTH
#endif

#include <opm/simulators/linalg/gpuistl/hypreinterface/HypreErrorHandling.hpp>

#include <HYPRE.h>
#include <_hypre_utilities.h>

#include <algorithm>
#include <vector>

namespace Opm::gpuistl::HypreInterface
{
/**
 * @brief Get matrix values from Hypre matrix
 *
 * Retrieves the values from a Hypre matrix for verification or debugging purposes.
 *
 * @param hypre_matrix The Hypre matrix to read values from
 * @param ncols Array containing number of columns per row
 * @param rows Array containing row indices
 * @param cols Array containing column indices
 * @param use_gpu_backend Whether the Hypre backend is using GPU acceleration
 * @return Vector containing the retrieved matrix values
 * @throws HypreError if the operation fails
 * @note For GPU backend, uses HYPRE_ParCSRMatrixGetRow since HYPRE_IJMatrixGetValues not supported on GPU
 */
inline std::vector<HYPRE_Real>
getMatrixValues(HYPRE_IJMatrix hypre_matrix,
                const std::vector<HYPRE_Int>& ncols,
                const std::vector<HYPRE_BigInt>& rows,
                const std::vector<HYPRE_BigInt>& cols,
                bool use_gpu_backend = false)
{
    const auto N = rows.size();
    std::vector<HYPRE_Real> values(cols.size());

    if (use_gpu_backend) {
        // HYPRE_IJMatrixGetValues not supported on GPU, use HYPRE_ParCSRMatrixGetRow instead
        void* object;
        OPM_HYPRE_SAFE_CALL(HYPRE_IJMatrixGetObject(hypre_matrix, &object));
        auto par_csr_matrix = static_cast<HYPRE_ParCSRMatrix>(object);
        HYPRE_Int max_row_size = *std::max_element(ncols.begin(), ncols.end());
        HYPRE_Complex* row_values_host = hypre_CTAlloc(HYPRE_Complex, max_row_size, HYPRE_MEMORY_HOST);

        size_t value_idx = 0;
        for (size_t row_idx = 0; row_idx < N; ++row_idx) {
            HYPRE_Int row_ncols;
            HYPRE_BigInt* row_cols_device;
            HYPRE_Complex* row_values_device;
            OPM_HYPRE_SAFE_CALL(HYPRE_ParCSRMatrixGetRow(
                par_csr_matrix, rows[row_idx], &row_ncols, &row_cols_device, &row_values_device));
            if (row_ncols > 0) {
                hypre_TMemcpy(row_values_host,
                              row_values_device,
                              HYPRE_Complex,
                              row_ncols,
                              HYPRE_MEMORY_HOST,
                              HYPRE_MEMORY_DEVICE);
                for (HYPRE_Int j = 0; j < row_ncols; ++j, ++value_idx) {
                    values[value_idx] = static_cast<HYPRE_Real>(row_values_host[j]);
                }
            }
            OPM_HYPRE_SAFE_CALL(HYPRE_ParCSRMatrixRestoreRow(
                par_csr_matrix, rows[row_idx], &row_ncols, &row_cols_device, &row_values_device));
        }
        hypre_TFree(row_values_host, HYPRE_MEMORY_HOST);
    } else {
        OPM_HYPRE_SAFE_CALL(HYPRE_IJMatrixGetValues(hypre_matrix,
                                                    N,
                                                    const_cast<HYPRE_Int*>(ncols.data()),
                                                    const_cast<HYPRE_BigInt*>(rows.data()),
                                                    const_cast<HYPRE_BigInt*>(cols.data()),
                                                    values.data()));
    }
    return values;
}
} // namespace Opm::gpuistl::HypreInterface

#endif // OPM_HYPRE_MATRIX_UTILS_HPP
