/*
  Copyright 2024 Equinor ASA

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

#ifndef OPM_ROCSPARSEBILU0_HPP
#define OPM_ROCSPARSEBILU0_HPP

#include <opm/simulators/linalg/gpubridge/BlockedMatrix.hpp>
#include <opm/simulators/linalg/gpubridge/WellContributions.hpp>

#include <opm/simulators/linalg/gpubridge/rocm/rocsparsePreconditioner.hpp>

#include <rocblas/rocblas.h>
#include <rocsparse/rocsparse.h>

#include <hip/hip_version.h>

namespace Opm::Accelerator {

/// This class implements a Blocked ILU0 preconditioner
/// The decomposition is done on GPU, using exact decomposition, or ChowPatel decomposition
/// The preconditioner is applied via two exact triangular solves
template <class Scalar, unsigned int block_size>
class rocsparseBILU0 : public rocsparsePreconditioner<Scalar, block_size>
{
    typedef rocsparsePreconditioner<Scalar, block_size> Base;

    using Base::N;
    using Base::Nb;
    using Base::nnz;
    using Base::nnzb;
    using Base::verbosity;

private:

    rocsparse_mat_descr descr_M, descr_L, descr_U;
    rocsparse_mat_info ilu_info;
#if HIP_VERSION >= 50400000
    rocsparse_mat_info spmv_info;
#endif

    rocsparse_int *d_Mrows, *d_Mcols;
    Scalar *d_Mvals, *d_t;
    void *d_buffer; // buffer space, used by rocsparse ilu0 analysis

    std::size_t d_bufferSize_M=0, d_bufferSize_L=0, d_bufferSize_U=0, d_bufferSize=0;

public:

    rocsparseBILU0(int verbosity_);
    ~rocsparseBILU0();

    /// Initialize GPU and allocate memory
    /// \param[in] matrix     matrix A
    /// \param[in] jacMatrix  matrix for preconditioner
    /// \param[in] d_Arows Array of row indices
    /// \param[in] d_Acols Array of column indices
    bool initialize(std::shared_ptr<BlockedMatrix<Scalar>> matrix,
                    std::shared_ptr<BlockedMatrix<Scalar>> jacMatrix,
                    rocsparse_int *d_Arows,
                    rocsparse_int *d_Acols) override;

    /// Analysis, extract parallelism if specified
    bool analyze_matrix();

    /// Analysis, extract parallelism if specified
    /// \param[in] mat     matrix A
    bool analyze_matrix(BlockedMatrix<Scalar> *mat) override;

    /// Analysis, extract parallelism if specified
    /// \param[in] mat     matrix A
    /// \param[in] jacMat  matrix for preconditioner, analyze this as well
    bool analyze_matrix(BlockedMatrix<Scalar> *mat,
                        BlockedMatrix<Scalar> *jacMat) override;

    /// ILU decomposition
    /// \param[in] mat     matrix A to decompose
    bool create_preconditioner(BlockedMatrix<Scalar> *mat) override;

    /// ILU decomposition
    /// \param[in] mat     matrix A
    /// \param[in] jacMat  matrix for preconditioner, decompose this one if used
    bool create_preconditioner(BlockedMatrix<Scalar> *mat,
                               BlockedMatrix<Scalar> *jacMat) override;

    /// Apply preconditioner, x = prec(y)
    /// via Lz = y
    /// and Ux = z
    /// \param[in]  y  Input y vector
    /// \param[out] x  Output x vector
    /// \param wellContribs Well contributions
    void apply(const Scalar& y,
               Scalar& x,
               WellContributions<Scalar>& wellContribs) override;

    /// Copy matrix A values to GPU
    /// \param[in]  mVals  Input values
    void copy_system_to_gpu(Scalar *mVals) override;

    /// Copy matrix A values to GPU
    /// \param[in]  mVals  Input values
    /// \param[in]  mRows Array of matrix row indices
    /// \param[in]  mCols Array of matrix column indices
    /// \param[in]  reuse True to reuse old matrix
    void copy_values_to_gpu(Scalar *mVals, int *mRows, int *mCols, bool reuse);

    /// Update GPU values after a new assembly is done
    /// \param[in] b     New b vector
    void update_system_on_gpu(Scalar*, Scalar* b) override;

};
} // namespace Opm

#endif
