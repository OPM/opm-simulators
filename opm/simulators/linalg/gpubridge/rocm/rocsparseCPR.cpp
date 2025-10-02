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

#include <config.h>
#include <opm/common/TimingMacros.hpp>
#include <opm/common/OpmLog/OpmLog.hpp>
#include <opm/common/ErrorMacros.hpp>
#include <dune/common/timer.hh>

#include <dune/common/shared_ptr.hh>

#include <opm/simulators/linalg/PreconditionerFactory.hpp>
#include <opm/simulators/linalg/PropertyTree.hpp>

#include <opm/simulators/linalg/gpubridge/GpuBridge.hpp>
#include <opm/simulators/linalg/gpubridge/BlockedMatrix.hpp>
#include <opm/simulators/linalg/gpubridge/rocm/rocsparseCPR.hpp>
#include <opm/simulators/linalg/gpubridge/rocm/hipKernels.hpp>

#include <opm/simulators/linalg/gpubridge/Misc.hpp>

#include <opm/simulators/linalg/gpubridge/rocm/rocsparseWellContributions.hpp>

#include <type_traits>

namespace Opm::Accelerator {

using Opm::OpmLog;
using Dune::Timer;

template <class Scalar, unsigned int block_size>
rocsparseCPR<Scalar, block_size>::rocsparseCPR(int verbosity_) :
    rocsparsePreconditioner<Scalar, block_size>(verbosity_)
{
    bilu0 = std::make_unique<rocsparseBILU0<Scalar, block_size> >(verbosity_);
    createfirsttime = true;
}

template <class Scalar, unsigned int block_size>
bool rocsparseCPR<Scalar, block_size>::
initialize(std::shared_ptr<BlockedMatrix<Scalar>> matrix,
           std::shared_ptr<BlockedMatrix<Scalar>> jacMatrix,
           rocsparse_int *d_Arows,
           rocsparse_int *d_Acols)
{
    h_mat = matrix;

    this->Nb = matrix->Nb;
    this->nnzb = matrix->nnzbs;
    this->N = Nb * block_size;
    this->nnz = nnzb * block_size * block_size;

    bilu0->set_context(this->handle, this->blas_handle, this->dir, this->operation, this->stream);
    bilu0->setJacMat(*this->jacMat.get());
    bilu0->initialize(matrix,jacMatrix,d_Arows,d_Acols);
    return true;
}

template <class Scalar, unsigned int block_size>
void rocsparseCPR<Scalar, block_size>::
copy_system_to_gpu(Scalar *b)
{
    bilu0->copy_system_to_gpu(b);
}

template <class Scalar, unsigned int block_size>
void rocsparseCPR<Scalar, block_size>::
update_system_on_gpu(Scalar* h_vals, Scalar *vals)
{
    bilu0->update_system_on_gpu(h_vals, vals);
    this->mat->nnzValues = h_vals;
}

template <class Scalar, unsigned int block_size>
bool rocsparseCPR<Scalar, block_size>::
analyze_matrix(BlockedMatrix<Scalar> *mat_)
{
    this->Nb = mat_->Nb;
    this->nnzb = mat_->nnzbs;
    this->N = Nb * block_size;
    this->nnz = nnzb * block_size * block_size;

    bool success = bilu0->analyze_matrix(mat_);

    this->mat = mat_;

    return success;
}

template <class Scalar, unsigned int block_size>
bool rocsparseCPR<Scalar, block_size>::
analyze_matrix(BlockedMatrix<Scalar> *mat_,
               BlockedMatrix<Scalar> *jacMat_)
{
    this->Nb = mat_->Nb;
    this->nnzb = mat_->nnzbs;
    this->N = Nb * block_size;
    this->nnz = nnzb * block_size * block_size;

    bool success = bilu0->analyze_matrix(mat_, jacMat_);
    this->mat = mat_;

    return success;
}

template <class Scalar, unsigned int block_size>
bool rocsparseCPR<Scalar, block_size>::
create_preconditioner([[maybe_unused]] BlockedMatrix<Scalar> *mat_,
                      BlockedMatrix<Scalar> *jacMat_)
{
    return create_preconditioner(jacMat_);
}

template <class Scalar, unsigned int block_size>
bool rocsparseCPR<Scalar, block_size>::
create_preconditioner(BlockedMatrix<Scalar> *mat_) {
    bool result = bilu0->create_preconditioner(mat_);
    bool amgcreated = this->create_preconditioner_amg(this->mat);

    if(amgcreated)
    {
        smootherilu0s.clear();
    }

    for (int i=0; i<this->num_levels-1; i++) {
        std::shared_ptr<rocsparseBILU0<Scalar,1>> smoother;

        //convert Matrix to BlockMatrix
        auto levelMat = std::make_shared<BlockedMatrix<Scalar>>(this->Amatrices[i].N, this->Amatrices[i].nnzs, 1, &(this->Amatrices[i].nnzValues[0]), &(this->Amatrices[i].colIndices[0]), &(this->Amatrices[i].rowPointers[0]));

        if(amgcreated)
        {
            smoother = std::make_shared<rocsparseBILU0<Scalar,1>>(verbosity);
            smootherilu0s.push_back(smoother);
            smoother->set_context(this->handle, this->blas_handle, this->dir, this->operation, this->stream);

            smoother->initialize(levelMat,levelMat,0,0);
            smoother->copy_values_to_gpu(&(this->Amatrices[i].nnzValues[0]), &(this->Amatrices[i].rowPointers[0]), &(this->Amatrices[i].colIndices[0]), false);
            smoother->analyze_matrix(levelMat.get());
         }
         else {
             smoother = smootherilu0s[i];

             smoother->copy_values_to_gpu(&(this->Amatrices[i].nnzValues[0]), 0, 0, true);
             smoother->analyze_matrix();
         }

         smoother->create_preconditioner(levelMat.get());
     }

     auto init_func = std::bind(&rocsparseCPR::init_rocm_buffers_fineMatrix, this);
     std::call_once(rocm_buffers_allocated, init_func);

     if(amgcreated) {
         if(!createfirsttime) clear_rocm_buffers();
         init_rocm_buffers();
         HIP_CHECK(hipMemcpyAsync(d_weights.data()->nnzValues, this->weights.data(), sizeof(Scalar) * this->N, hipMemcpyHostToDevice, this->stream));
         createfirsttime = false;
     }

     // upload matrices and vectors to GPU
     rocm_upload();

     return result;
}

template <class Scalar, unsigned int block_size>
void rocsparseCPR<Scalar, block_size>::
init_rocm_buffers_fineMatrix()
{
    d_mat = std::make_unique<RocmMatrix<Scalar>>(h_mat.get()->Nb, h_mat.get()->Nb, h_mat.get()->nnzbs, block_size);
    d_weights.emplace_back(this->N);
}

template <class Scalar, unsigned int block_size>
void rocsparseCPR<Scalar, block_size>::
init_rocm_buffers()
{
    d_Amatrices.reserve(this->num_levels);
    d_Rmatrices.reserve(this->num_levels - 1);
    d_invDiags.reserve(this->num_levels - 1);
    d_f.reserve(this->num_levels - 1);
    d_u.reserve(this->num_levels - 1);
    d_PcolIndices.reserve(this->num_levels - 1);
    d_invDiags.reserve(this->num_levels - 1);
    d_t.reserve(this->num_levels - 1);

    for (Matrix<Scalar>& m : this->Amatrices) {
        d_Amatrices.emplace_back(m.N, m.M, m.nnzs, 1);
    }

    for (Matrix<Scalar>& m : this->Rmatrices) {
        d_Rmatrices.emplace_back(m.N, m.M, m.nnzs, 1);
        d_f.emplace_back(m.N);
        d_u.emplace_back(m.N);
        d_PcolIndices.emplace_back(m.M);
        d_invDiags.emplace_back(m.M);
        d_t.emplace_back(m.M);
    }

    HIP_CHECK(hipMalloc((void**)&d_tmp, sizeof(Scalar) * this->N));
    HIP_CHECK(hipMalloc((void**)&d_yamg, sizeof(Scalar) * this->N));
    HIP_CHECK(hipMalloc((void**)&d_ywell, sizeof(Scalar) * this->N));

    d_rs.emplace_back(this->N);

    d_coarse_y.emplace_back(this->Nb);
    d_coarse_x.emplace_back(this->Nb);
    d_yinput.emplace_back(this->Nb);
}

template <class Scalar, unsigned int block_size>
void rocsparseCPR<Scalar, block_size>::
clear_rocm_buffers()
{
    d_Amatrices.clear();
    d_Rmatrices.clear();
    d_invDiags.clear();
    d_f.clear();
    d_u.clear();
    d_PcolIndices.clear();
    d_t.clear();

    HIP_CHECK(hipFree(d_tmp));
    HIP_CHECK(hipFree(d_yamg));
    HIP_CHECK(hipFree(d_ywell));

    d_rs.clear();

    d_coarse_y.clear();
    d_coarse_x.clear();
    d_yinput.clear();
}

template <class Scalar, unsigned int block_size>
void rocsparseCPR<Scalar, block_size>::
rocm_upload()
{
     d_mat->upload(h_mat.get(), this->stream);

     for (unsigned int i = 0; i < this->Rmatrices.size(); ++i) {
         d_Amatrices[i].upload(&this->Amatrices[i], this->stream);
         d_Rmatrices[i].upload(&this->Rmatrices[i], this->stream);

         HIP_CHECK(hipMemcpyAsync(d_invDiags[i].nnzValues, this->invDiags[i].data(), sizeof(Scalar) * this->Amatrices[i].N, hipMemcpyHostToDevice, this->stream));
         HIP_CHECK(hipMemcpyAsync(d_PcolIndices[i].nnzValues, this->PcolIndices[i].data(), sizeof(int) * this->Amatrices[i].N, hipMemcpyHostToDevice, this->stream));
     }
}

template <class Scalar, unsigned int block_size>
void rocsparseCPR<Scalar, block_size>::
amg_cycle_gpu(const int level,
              Scalar &y,
              Scalar &x,
              [[maybe_unused]] WellContributions<Scalar>& wellContribs)
{
    Scalar one  = 1.0;
    Scalar zero = 0.0;
    RocmMatrix<Scalar> *A = &d_Amatrices[level];
    RocmMatrix<Scalar> *R = &d_Rmatrices[level];

    std::shared_ptr<rocsparseBILU0<Scalar,1>> S;

    int Ncur = A->Nb;

    rocsparse_mat_descr descr_R;
    ROCSPARSE_CHECK(rocsparse_create_mat_descr(&descr_R));

    if (level == this->num_levels - 1) {
        // solve coarsest level
        std::vector<Scalar> h_y(Ncur), h_x(Ncur, 0);

        HIP_CHECK(hipMemcpyAsync(h_y.data(), &y, sizeof(Scalar) * Ncur, hipMemcpyDeviceToHost, this->stream));

        // The if constexpr is needed to make the code compile
        // since the umfpack member is an 'int' with float Scalar.
        // We will never get here with float Scalar as we throw earlier.
        // Solve coarsest level using umfpack
        if constexpr (std::is_same_v<Scalar,double>) {
            this->umfpack.apply(h_x.data(), h_y.data());
        }

        HIP_CHECK(hipMemcpyAsync(&x, h_x.data(), sizeof(Scalar) * Ncur, hipMemcpyHostToDevice, this->stream));

        return;
    }
    else
        S = smootherilu0s[level];

    RocmVector<Scalar>& t = d_t[level]; // this is y hierarchy ==> this vector needs to be the updated
    RocmVector<Scalar>& f = d_f[level]; // this is x=lhs hierarchy starting from the first level
    RocmVector<Scalar>& u = d_u[level]; // u was 0-initialized earlier : this is y=rhs from the first level

    // presmooth
    for (unsigned i = 0; i < this->num_pre_smooth_steps; ++i){
        S->apply(y, x, wellContribs);

        if constexpr (std::is_same_v<Scalar,float>) {
            ROCBLAS_CHECK(rocblas_saxpy(this->blas_handle, Ncur, &one, &x, 1, t.nnzValues, 1)); //t is output, x input
        } else {
            ROCBLAS_CHECK(rocblas_daxpy(this->blas_handle, Ncur, &one, &x, 1, t.nnzValues, 1)); //t is output, x input
        }

        HipKernels<Scalar>::residual(A->nnzValues, A->colIndices, A->rowPointers, &x, &y, &y, Ncur, 1, this->stream);
    }

    // move to coarser level
    ROCSPARSE_CHECK(rocsparse_create_mat_info(&spmv_info));
    if constexpr (std::is_same_v<Scalar,float>) {
        ROCSPARSE_CHECK(rocsparse_scsrmv_analysis(this->handle, this->operation, R->Nb, R->Mb, R->nnzbs, descr_R, R->nnzValues, R->rowPointers, R->colIndices, spmv_info));
    } else {
        ROCSPARSE_CHECK(rocsparse_dcsrmv_analysis(this->handle, this->operation, R->Nb, R->Mb, R->nnzbs, descr_R, R->nnzValues, R->rowPointers, R->colIndices, spmv_info));
    }

    if constexpr (std::is_same_v<Scalar,float>) {
            ROCSPARSE_CHECK(rocsparse_scsrmv(this->handle, this->operation,
                                                R->Nb, R->Mb, R->nnzbs, &one, descr_R,
                                                R->nnzValues, R->rowPointers, R->colIndices, spmv_info,
                                                &y/*t.nnzValues*/, &zero, f.nnzValues));
    } else {
            ROCSPARSE_CHECK(rocsparse_dcsrmv(this->handle, this->operation,
                                                R->Nb, R->Mb, R->nnzbs, &one, descr_R,
                                                R->nnzValues, R->rowPointers, R->colIndices, spmv_info,
                                                &y/*t.nnzValues*/, &zero, f.nnzValues));
    }

    ROCSPARSE_CHECK(rocsparse_destroy_mat_info(spmv_info));

    amg_cycle_gpu(level + 1, *f.nnzValues, *u.nnzValues, wellContribs);//f is input = y, u is output = x!

    HIP_CHECK(hipMemsetAsync(&x, 0, sizeof(Scalar) * Ncur, this->stream));
    if(level < this->num_levels - 2)
    {
        HipKernels<Scalar>::prolongate_vector(d_t[level+1].nnzValues, &x, d_PcolIndices[level].nnzValues, Ncur, this->stream);
    }
    else {
        HipKernels<Scalar>::prolongate_vector(u.nnzValues, &x, d_PcolIndices[level].nnzValues, Ncur, this->stream);
    }

    // accumulate update
    if constexpr (std::is_same_v<Scalar,float>) {
        ROCBLAS_CHECK(rocblas_saxpy(this->blas_handle, Ncur, &one, &x, 1, t.nnzValues, 1));
    } else {
        ROCBLAS_CHECK(rocblas_daxpy(this->blas_handle, Ncur, &one, &x, 1, t.nnzValues, 1));
    }

    // postsmooth
    for (unsigned i = 0; i < this->num_post_smooth_steps; ++i){
        // update defect
        HipKernels<Scalar>::residual(A->nnzValues, A->colIndices, A->rowPointers, &x, &y, &y, Ncur, 1, this->stream); // params: x@in, y@in, y@out = (A.., x, rhs, out, Nb, block_size, stream) -> computes out = rhs - A * x <==> y -= Ax

        //smoother applier
        S->apply(y, x, wellContribs);

        // accumulate update
        if constexpr (std::is_same_v<Scalar,float>) {
            ROCBLAS_CHECK(rocblas_saxpy(this->blas_handle, Ncur, &one, &x, 1, t.nnzValues, 1));
        } else {
            ROCBLAS_CHECK(rocblas_daxpy(this->blas_handle, Ncur, &one, &x, 1, t.nnzValues, 1));
        }
    }
}

// x = prec(y)
template <class Scalar, unsigned int block_size>
void rocsparseCPR<Scalar, block_size>::
apply_amg(const Scalar& y,
          Scalar& x,
          [[maybe_unused]] WellContributions<Scalar>& wellContribs)
{
    HIP_CHECK(hipMemsetAsync(d_coarse_x.data()->nnzValues, 0, sizeof(Scalar) * this->Nb, this->stream));

    for (unsigned int i = 0; i < d_t.size(); ++i) {
        HIP_CHECK(hipMemsetAsync(d_t[i].nnzValues, 0, sizeof(Scalar) * d_t[i].size, this->stream));
    }

    HipKernels<Scalar>::full_to_pressure_restriction(&y, d_weights.data()->nnzValues, d_coarse_y.data()->nnzValues, Nb, this->stream);

    //NOTE: input b=y is is overwritten, so duplicate into d_yinput
    HIP_CHECK(hipMemcpyAsync(d_yinput.data()->nnzValues, d_coarse_y.data()->nnzValues, sizeof(Scalar) * d_coarse_y.data()->size, hipMemcpyDeviceToDevice, this->stream));

    amg_cycle_gpu(0, *(d_coarse_y.data()->nnzValues), *(d_coarse_x.data()->nnzValues), wellContribs);

    RocmMatrix<Scalar> *A = &d_Amatrices[0];

    // update defect
    HipKernels<Scalar>::residual(A->nnzValues, A->colIndices, A->rowPointers, d_t[0].nnzValues, d_yinput.data()->nnzValues, d_yinput.data()->nnzValues, d_yinput.data()->size, 1, this->stream); // params: x@in, y@in, y@out = (A.., x, rhs, out, Nb, block_size, stream)

    //NOTE: this initialization is required to get mathcing result with spmv in applyscaleadd of the well operator
    HIP_CHECK(hipMemsetAsync(&x, 0,  sizeof(Scalar) *N, this->stream));
    HIP_CHECK(hipStreamSynchronize(this->stream));

    HipKernels<Scalar>::add_coarse_pressure_correction(d_t[0].nnzValues, &x, this->pressure_idx, Nb, this->stream);
}

template <class Scalar, unsigned int block_size>
void rocsparseCPR<Scalar, block_size>::
apply(const Scalar& y,
      Scalar& x,
      [[maybe_unused]] WellContributions<Scalar>& wellContribs)
{
    Scalar one  = 1.0;
    Scalar mone = -1.0;
    rocsparse_mat_descr descr_mat;

    HIP_CHECK(hipMemsetAsync(&x, 0,  sizeof(Scalar) *N, this->stream));
    HIP_CHECK(hipMemcpyAsync(d_yamg, &y, sizeof(Scalar) * this->N, hipMemcpyDeviceToDevice, this->stream));
    HIP_CHECK(hipStreamSynchronize(this->stream));

    apply_amg(*d_yamg, x, wellContribs);

    ROCSPARSE_CHECK(rocsparse_create_mat_descr(&descr_mat));

#if HIP_VERSION >= 50400000
    ROCSPARSE_CHECK(rocsparse_create_mat_info(&spmv_info));
#endif

#if HIP_VERSION >= 60000000
    if constexpr (std::is_same_v<Scalar,float>) {
        ROCSPARSE_CHECK(rocsparse_sbsrmv_analysis(this->handle, this->dir, this->operation, Nb, Nb, N, descr_mat, d_mat->nnzValues, d_mat->rowPointers, d_mat->colIndices, block_size, spmv_info));
    } else {
        ROCSPARSE_CHECK(rocsparse_dbsrmv_analysis(this->handle, this->dir, this->operation, Nb, Nb, N, descr_mat, d_mat->nnzValues, d_mat->rowPointers, d_mat->colIndices, block_size, spmv_info));
    }
#elif HIP_VERSION >= 50400000
    if constexpr (std::is_same_v<Scalar,float>) {
        ROCSPARSE_CHECK(rocsparse_sbsrmv_ex_analysis(this->handle, this->dir, this->operation, Nb, Nb, N, descr_mat, d_mat->nnzValues, d_mat->rowPointers, d_mat->colIndices, block_size, spmv_info));
    } else {
        ROCSPARSE_CHECK(rocsparse_dbsrmv_ex_analysis(this->handle, this->dir, this->operation, Nb, Nb, N, descr_mat, d_mat->nnzValues, d_mat->rowPointers, d_mat->colIndices, block_size, spmv_info));
    }
#endif

    //NOTE: below is the correct y -= AX but it is curretly not working so I am debugging spmv = Ax above alone!
#if HIP_VERSION >= 60000000
    if constexpr (std::is_same_v<Scalar,float>) {
            ROCSPARSE_CHECK(rocsparse_sbsrmv(this->handle, this->dir, this->operation,
                                                Nb, Nb, d_mat->nnzbs /*N*/, &mone, descr_mat,
                                                d_mat->nnzValues, d_mat->rowPointers, d_mat->colIndices, block_size, spmv_info, &x, &one, d_yamg));
    } else {
            ROCSPARSE_CHECK(rocsparse_dbsrmv(this->handle, this->dir, this->operation,
                                                Nb, Nb, d_mat->nnzbs /*N*/, &mone, descr_mat,
                                                d_mat->nnzValues, d_mat->rowPointers, d_mat->colIndices, block_size, spmv_info, &x, &one, d_yamg));
    }
#elif HIP_VERSION >= 50400000
    if constexpr (std::is_same_v<Scalar,float>) {
            ROCSPARSE_CHECK(rocsparse_sbsrmv_ex(this->handle, this->dir, this->operation,
                                                Nb, Nb, d_mat->nnzbs /*N*/, &mone, descr_mat,
                                                d_mat->nnzValues, d_mat->rowPointers, d_mat->colIndices, block_size, spmv_info, &x, &one, d_yamg));
    } else {
            ROCSPARSE_CHECK(rocsparse_dbsrmv_ex(this->handle, this->dir, this->operation,
                                                Nb, Nb, d_mat->nnzbs /*N*/, &mone, descr_mat,
                                                d_mat->nnzValues, d_mat->rowPointers, d_mat->colIndices, block_size, spmv_info, &x, &one, d_yamg));
    }
#else
    if constexpr (std::is_same_v<Scalar,float>) {
            ROCSPARSE_CHECK(rocsparse_sbsrmv(this->handle, this->dir, this->operation,
                                                Nb, Nb, d_mat->nnzbs /*N*/, &mone, descr_mat,
                                                d_mat->nnzValues, d_mat->rowPointers, d_mat->colIndices, block_size, &x, &one, d_yamg));
    } else {
            ROCSPARSE_CHECK(rocsparse_dbsrmv(this->handle, this->dir, this->operation,
                                                Nb, Nb, d_mat->nnzbs /*N*/, &mone, descr_mat,
                                                d_mat->nnzValues, d_mat->rowPointers, d_mat->colIndices, block_size, &x, &one, d_yamg));
    }
#endif

#if HIP_VERSION >= 50400000
    ROCSPARSE_CHECK(rocsparse_destroy_mat_info(spmv_info));
#endif

    ROCSPARSE_CHECK(rocsparse_destroy_mat_descr(descr_mat));

    HIP_CHECK(hipMemsetAsync(d_ywell, 0,  sizeof(Scalar) *N, this->stream));
    // NOTE: we need to also apply the wells to the defect!!!
    // apply wellContributions
    if (wellContribs.getNumWells() > 0) {
        static_cast<WellContributionsRocsparse<Scalar>&>(wellContribs).apply(&x, d_ywell);
    }

    // add contrib of wells to y: Ax = Ax + alpha * scaleAddRes_ <==> Ax.axpy( alpha, scaleAddRes_ )
    if constexpr (std::is_same_v<Scalar,float>) {
        ROCBLAS_CHECK(rocblas_saxpy(this->blas_handle, N, &mone, d_ywell, 1, d_yamg, 1));
    } else {
        ROCBLAS_CHECK(rocblas_daxpy(this->blas_handle, N, &mone, d_ywell, 1, d_yamg, 1));
    }

    HIP_CHECK(hipMemsetAsync(d_tmp, 0, N * sizeof(Scalar), this->stream));

    bilu0->apply(*d_yamg, *d_tmp, wellContribs);

    HIP_CHECK(hipStreamSynchronize(this->stream));

    // Accumulate update:  *levelContext.update += *levelContext.lhs;
    if constexpr (std::is_same_v<Scalar,float>) {
        ROCBLAS_CHECK(rocblas_saxpy(this->blas_handle, N, &one, d_tmp, 1, &x, 1));
    } else {
        ROCBLAS_CHECK(rocblas_daxpy(this->blas_handle, N, &one, d_tmp, 1, &x, 1));
    }
}

#define INSTANTIATE_TYPE(T)           \
    template class rocsparseCPR<T,1>; \
    template class rocsparseCPR<T,2>; \
    template class rocsparseCPR<T,3>; \
    template class rocsparseCPR<T,4>; \
    template class rocsparseCPR<T,5>; \
    template class rocsparseCPR<T,6>;

INSTANTIATE_TYPE(double)

#if FLOW_INSTANTIATE_FLOAT
INSTANTIATE_TYPE(float)
#endif

} // namespace Opm
