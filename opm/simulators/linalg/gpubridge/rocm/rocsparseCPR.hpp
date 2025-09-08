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

#ifndef OPM_ROCSPARSECPR_HPP
#define OPM_ROCSPARSECPR_HPP

#include <mutex>

#include <opm/simulators/linalg/gpubridge/rocm/rocsparseBILU0.hpp>
#include <opm/simulators/linalg/gpubridge/Matrix.hpp>
#include <opm/simulators/linalg/gpubridge/CprCreation.hpp>
#include <opm/simulators/linalg/gpubridge/rocm/rocsparseMatrix.hpp>
#include <opm/simulators/linalg/gpubridge/rocm/rocsparsePreconditioner.hpp>
#include <opm/simulators/linalg/gpubridge/WellContributions.hpp>

#include <opm/simulators/linalg/gpubridge/rocm/rocsparseSolverBackend.hpp>

namespace Opm::Accelerator {

template<class Scalar> class BlockedMatrix;

/// This class implements a Constrained Pressure Residual (CPR) preconditioner
template <class Scalar, unsigned int block_size>
class rocsparseCPR : public rocsparsePreconditioner<Scalar, block_size>, public CprCreation<Scalar, block_size>
{
    typedef rocsparsePreconditioner<Scalar, block_size> Base;

    using Base::N;
    using Base::Nb;
    using Base::nnz;
    using Base::nnzb;
    using Base::verbosity;

private:
    std::vector<RocmMatrix<Scalar>> d_Amatrices, d_Rmatrices; // scalar matrices that represent the AMG hierarchy
    
    std::vector<RocmVector<int>> d_PcolIndices; // prolongation does not need a full matrix, only store colIndices
    std::vector<RocmVector<Scalar>> d_invDiags; // inverse of diagonal of Amatrices
    std::vector<RocmVector<Scalar>> d_t, d_f; // intermediate vectors used during amg cycle
    std::vector<RocmVector<Scalar>> d_u; // intermediate vectors used during amg cycle
    std::vector<RocmVector<Scalar>> d_rs;      // use before extracting the pressure
    std::vector<RocmVector<Scalar>> d_weights; // the quasiimpes weights, used to extract pressure
    std::unique_ptr<RocmMatrix<Scalar>> d_mat;   // stores blocked matrix
    std::shared_ptr<BlockedMatrix<Scalar>> h_mat;

    std::vector<RocmVector<Scalar>> d_coarse_y, d_coarse_x, d_yinput; // stores the scalar vectors
    std::once_flag rocm_buffers_allocated;  // only allocate OpenCL Buffers once

    std::unique_ptr<rocsparseBILU0<Scalar, block_size> > bilu0;                    // Blocked ILU0 preconditioner

    std::vector<std::shared_ptr<rocsparseBILU0<Scalar,1>>> smootherilu0s;

    Scalar *d_tmp, *d_yamg, *d_ywell;
    bool createfirsttime;

    #if HIP_VERSION >= 50400000
    rocsparse_mat_info spmv_info;
    #endif

    // Initialize and allocate matrices and vectors
    void init_rocm_buffers();
    void init_rocm_buffers_fineMatrix();

    // Deallocate matrices and vectors
    void clear_rocm_buffers();

    // Copy matrices and vectors to GPU
    void rocm_upload();

    // apply pressure correction to vector
    void apply_amg(const Scalar& y,
                   Scalar& x,
                   WellContributions<Scalar>& wellContribs);

    // Apply the AMG preconditioner
    void amg_cycle_gpu(const int level,
                       Scalar &y,
                       Scalar &x,
                       WellContributions<Scalar>& wellContribs);
    
public:

    rocsparseCPR(int verbosity);

    /// Initialize GPU and allocate memory
    /// \param[in] matrix     matrix A
    /// \param[in] jacMatrix  matrix for preconditioner
    /// \param[in] d_Arows Array of matrix row indices
    /// \param[in] d_Acols Array of matrix column indices
    bool initialize(std::shared_ptr<BlockedMatrix<Scalar>> matrix,
                    std::shared_ptr<BlockedMatrix<Scalar>> jacMatrix,
                    rocsparse_int *d_Arows,
                    rocsparse_int *d_Acols) override;
    

    /// Analysis, extract parallelism if specified
    /// \param[in] mat     matrix A
    bool analyze_matrix(BlockedMatrix<Scalar> *mat) override;
    
    /// Analysis, extract parallelism if specified
    /// \param[in] mat     matrix A
    /// \param[in] jacMat  matrix for preconditioner, analyze this as well
    bool analyze_matrix(BlockedMatrix<Scalar> *mat,
                        BlockedMatrix<Scalar> *jacMat) override;

    /// Create AMG preconditioner and perform ILU decomposition
    /// \param[in] mat     matrix A
    bool create_preconditioner(BlockedMatrix<Scalar> *mat) override;
    
    /// Create AMG preconditioner and perform ILU decomposition
    /// \param[in] mat     matrix A
    /// \param[in] jacMat  matrix for preconditioner, decompose this one if used
    bool create_preconditioner(BlockedMatrix<Scalar> *mat,
                               BlockedMatrix<Scalar> *jacMat) override;
    
    /// Apply preconditioner, x = prec(y)
    /// applies blocked ilu0
    /// also applies amg for pressure component
    /// \param[in]  y  Input y vector
    /// \param[out] x  Output x vector
    /// \param wellContribs Well contributions
    void apply(const Scalar& y,
               Scalar& x,
               WellContributions<Scalar>& wellContribs) override;
    
    /// Copy matrix A values to GPU
    /// \param[in]  b Input values
    void copy_system_to_gpu(Scalar *b) override;

    /// Update linear system to GPU
    /// \param[in] vals           Matrix values
    /// \param[in] b              input vector, contains N values
    void update_system_on_gpu(Scalar* vals, Scalar* b) override;
    
};

} // namespace Opm

#endif

