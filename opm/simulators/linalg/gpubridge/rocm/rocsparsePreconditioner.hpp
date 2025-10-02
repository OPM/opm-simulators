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

#ifndef OPM_ROCSPARSEPRECONDITIONER_HEADER_INCLUDED
#define OPM_ROCSPARSEPRECONDITIONER_HEADER_INCLUDED

#include <opm/simulators/linalg/gpubridge/Preconditioner.hpp>

#include <rocsparse/rocsparse.h>
#include <rocblas/rocblas.h>

namespace Opm::Accelerator {

template<class Scalar> class BlockedMatrix;

template <class Scalar, unsigned int block_size>
class rocsparsePreconditioner : public Preconditioner<Scalar, block_size>
{

protected:
    rocsparse_handle handle;
    rocblas_handle blas_handle;
    rocsparse_direction dir = rocsparse_direction_row;
    rocsparse_operation operation = rocsparse_operation_none;
    rocsparse_mat_descr descr_L, descr_U;

    hipStream_t stream;

    rocsparsePreconditioner(int verbosity_) :
    Preconditioner<Scalar, block_size>(verbosity_)
    {};

public:

    int nnzbs_prec = 0; // number of nnz blocks in preconditioner matrix M
    bool useJacMatrix = false;
    std::shared_ptr<BlockedMatrix<Scalar>> jacMat{}; // matrix for preconditioner

    virtual ~rocsparsePreconditioner() = default;

    static std::unique_ptr<rocsparsePreconditioner<Scalar, block_size>> create(PreconditionerType type,
                                                                               int verbosity);

    virtual bool initialize(std::shared_ptr<BlockedMatrix<Scalar>> matrix,
                            std::shared_ptr<BlockedMatrix<Scalar>> jacMatrix,
                            rocsparse_int* d_Arows,
                            rocsparse_int* d_Acols) = 0;

    virtual void copy_system_to_gpu(Scalar* b) = 0;

    /// Update linear system to GPU
    /// \param[in] vals           Matrix values
    /// \param[in] b              input vector, contains N values
    virtual void update_system_on_gpu(Scalar* vals, Scalar* b)=0;

    void set_matrix_analysis(rocsparse_mat_descr descr_L,
                             rocsparse_mat_descr descr_U);

    void set_context(rocsparse_handle handle,
                     rocblas_handle blas_handle,
                     rocsparse_direction dir,
                     rocsparse_operation operation,
                     hipStream_t stream);

    void setJacMat(const BlockedMatrix<Scalar>& jacMat);
};
} //namespace Opm

#endif
