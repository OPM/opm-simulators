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

#include <mutex>
#include <vector>

#include <opm/simulators/linalg/bda/BlockedMatrix.hpp>

#include <opm/simulators/linalg/bda/rocm/rocsparsePreconditioner.hpp>

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
    double *d_Mvals, *d_t;
    void *d_buffer; // buffer space, used by rocsparse ilu0 analysis
    
public:

    rocsparseBILU0(int verbosity_);
    
    /// Initialize GPU and allocate memory
    /// \param[in] matrix     matrix A
    /// \param[in] jacMatrix  matrix for preconditioner
    bool initialize(std::shared_ptr<BlockedMatrix<Scalar>> matrix, std::shared_ptr<BlockedMatrix<Scalar>> jacMatrix, rocsparse_int *d_Arows, rocsparse_int *d_Acols);

    // analysis, extract parallelism if specified
    bool analyze_matrix(BlockedMatrix<Scalar> *mat);
    bool analyze_matrix(BlockedMatrix<Scalar> *mat, BlockedMatrix<Scalar> *jacMat);

    // ilu_decomposition
    bool create_preconditioner(BlockedMatrix<Scalar> *mat) override;
    bool create_preconditioner(BlockedMatrix<Scalar> *mat, BlockedMatrix<Scalar> *jacMat) override;

    // apply preconditioner, x = prec(y)
    // via Lz = y
    // and Ux = z
    void apply(double& y, double& x) override;

#if HAVE_OPENCL
    // apply preconditioner, x = prec(y)
    void apply(const cl::Buffer& y, cl::Buffer& x) {}
#endif

    void copy_system_to_gpu(double *mVals);
    void update_system(double *vals, double *b);
    void update_system_on_gpu(double *b);
};
} // namespace Opm

#endif

