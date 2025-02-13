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

#include <type_traits>

namespace Opm::Accelerator {

using Opm::OpmLog;
using Dune::Timer;

template <class Scalar, unsigned int block_size>
rocsparseCPR<Scalar, block_size>::rocsparseCPR(int verbosity_) :
    rocsparsePreconditioner<Scalar, block_size>(verbosity_)
{
    bilu0 = std::make_unique<rocsparseBILU0<Scalar, block_size> >(verbosity_);
}

template <class Scalar, unsigned int block_size>
bool rocsparseCPR<Scalar, block_size>::
initialize(std::shared_ptr<BlockedMatrix<Scalar>> matrix,
           std::shared_ptr<BlockedMatrix<Scalar>> jacMatrix,
           rocsparse_int *d_Arows,
           rocsparse_int *d_Acols)
{
    this->Nb = matrix->Nb;
    this->nnzb = matrix->nnzbs;
    this->N = Nb * block_size;
    this->nnz = nnzb * block_size * block_size;

    bilu0->set_context(this->handle, this->dir, this->operation, this->stream);
    bilu0->setJacMat(*this->jacMat.get());
    bilu0->initialize(matrix,jacMatrix,d_Arows,d_Acols);
    return true;
}

template <class Scalar, unsigned int block_size>
void rocsparseCPR<Scalar, block_size>::
copy_system_to_gpu(Scalar *b) {
    bilu0->copy_system_to_gpu(b);
}

template <class Scalar, unsigned int block_size>
void rocsparseCPR<Scalar, block_size>::
update_system_on_gpu(Scalar *vals) {
    bilu0->update_system_on_gpu(vals);
}

template <class Scalar, unsigned int block_size>
bool rocsparseCPR<Scalar, block_size>::
analyze_matrix(BlockedMatrix<Scalar> *mat_) {
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
create_preconditioner(BlockedMatrix<Scalar> *mat_,
                      BlockedMatrix<Scalar> *jacMat_)
{
    Dune::Timer t_bilu0;
    bool result = bilu0->create_preconditioner(mat_, jacMat_);
    if (verbosity >= 3) {
        std::ostringstream out;
        out << "rocsparseCPR create_preconditioner bilu0(): " << t_bilu0.stop() << " s";
        OpmLog::info(out.str());
    }
    
    Dune::Timer t_amg;
    this->create_preconditioner_amg(this->mat);

    if (verbosity >= 3) {
        std::ostringstream out;
        out << "rocsparseCPR create_preconditioner_amg(): " << t_amg.stop() << " s";
        OpmLog::info(out.str());
    }
    
    auto init_func = std::bind(&rocsparseCPR::init_rocm_buffers, this);
    std::call_once(rocm_buffers_allocated, init_func);

    // upload matrices and vectors to GPU
    rocm_upload();

    return result;
}

template <class Scalar, unsigned int block_size>
bool rocsparseCPR<Scalar, block_size>::
create_preconditioner(BlockedMatrix<Scalar> *mat_) {
    Dune::Timer t_bilu0;
    bool result = bilu0->create_preconditioner(mat_);
    if (verbosity >= 3) {
        std::ostringstream out;
        out << "rocsparseCPR create_preconditioner bilu0(): " << t_bilu0.stop() << " s";
        OpmLog::info(out.str());
    }

    Dune::Timer t_amg;
    this->create_preconditioner_amg(this->mat); 
    
    if (verbosity >= 3) {
        std::ostringstream out;
        out << "rocsparseCPR create_preconditioner_amg(): " << t_amg.stop() << " s";
        OpmLog::info(out.str());
    }
    
    auto init_func = std::bind(&rocsparseCPR::init_rocm_buffers, this);
    std::call_once(rocm_buffers_allocated, init_func);

    // upload matrices and vectors to GPU
    rocm_upload();
    
    return result;
}

template <class Scalar, unsigned int block_size>
void rocsparseCPR<Scalar, block_size>::
init_rocm_buffers() {
    d_Amatrices.reserve(this->num_levels);
    d_Rmatrices.reserve(this->num_levels - 1);
    d_invDiags.reserve(this->num_levels - 1);
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
    
    d_weights.emplace_back(this->N);
    d_rs.emplace_back(this->N);
    d_mat = std::make_unique<RocmMatrix<Scalar>>(this->Nb, this->Nb, this->nnzb, block_size);
    
    d_coarse_y.emplace_back(this->Nb);
    d_coarse_x.emplace_back(this->Nb);
}

template <class Scalar, unsigned int block_size>
void rocsparseCPR<Scalar, block_size>::
rocm_upload() {
     d_mat->upload(this->mat, this->stream);
     
     HIP_CHECK(hipMemcpyAsync(d_weights.data()->nnzValues, this->weights.data(), sizeof(Scalar) * this->N, hipMemcpyHostToDevice, this->stream));
    
    for (unsigned int i = 0; i < this->Rmatrices.size(); ++i) {
        d_Amatrices[i].upload(&this->Amatrices[i], this->stream);
        
        HIP_CHECK(hipMemcpyAsync(d_invDiags[i].nnzValues, this->invDiags[i].data(), sizeof(Scalar) * this->Amatrices[i].N, hipMemcpyHostToDevice, this->stream));
        HIP_CHECK(hipMemcpyAsync(d_PcolIndices[i].nnzValues, this->PcolIndices[i].data(), sizeof(int) * this->Amatrices[i].N, hipMemcpyHostToDevice, this->stream));
    }
    
    for (unsigned int i = 0; i < this->Rmatrices.size(); ++i) {
        d_Rmatrices[i].upload(&this->Rmatrices[i], this->stream);
    }
}

template <class Scalar, unsigned int block_size>
void rocsparseCPR<Scalar, block_size>::
amg_cycle_gpu(const int level,
              Scalar &y,
              Scalar &x)
{
    RocmMatrix<Scalar> *A = &d_Amatrices[level];
    RocmMatrix<Scalar> *R = &d_Rmatrices[level];
    int Ncur = A->Nb;
    
    rocsparse_mat_info spmv_info;
    rocsparse_mat_descr descr_R;
    ROCSPARSE_CHECK(rocsparse_create_mat_info(&spmv_info));
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
    
    int Nnext = d_Amatrices[level+1].Nb;

    RocmVector<Scalar>& t = d_t[level];
    RocmVector<Scalar>& f = d_f[level];
    RocmVector<Scalar>& u = d_u[level]; // u was 0-initialized earlier

    // presmooth
    Scalar jacobi_damping = 0.65; // default value in amgcl: 0.72
    for (unsigned i = 0; i < this->num_pre_smooth_steps; ++i){
        HipKernels<Scalar>::residual(A->nnzValues, A->colIndices, A->rowPointers, &x, &y, t.nnzValues, Ncur, 1, this->stream);
        HipKernels<Scalar>::vmul(jacobi_damping, d_invDiags[level].nnzValues, t.nnzValues, &x, Ncur, this->stream);
    }

    // move to coarser level
    HipKernels<Scalar>::residual(A->nnzValues, A->colIndices, A->rowPointers, &x, &y, t.nnzValues, Ncur, 1, this->stream);

// TODO: understand why rocsparse spmv library function does not here. 
//     ROCSPARSE_CHECK(rocsparse_dbsrmv(this->handle, this->dir, this->operation,
//                                          R->Nb, R->Mb, R->nnzbs, &one, descr_R,
//                                          R->nnzValues, R->rowPointers, R->colIndices, 1,
//                                          t.nnzValues, &zero, f.nnzValues));
    HipKernels<Scalar>::spmv(R->nnzValues, R->colIndices, R->rowPointers, t.nnzValues, f.nnzValues, Nnext, 1, this->stream);

    amg_cycle_gpu(level + 1, *f.nnzValues, *u.nnzValues);
    HipKernels<Scalar>::prolongate_vector(u.nnzValues, &x, d_PcolIndices[level].nnzValues, Ncur, this->stream);

    // postsmooth
    for (unsigned i = 0; i < this->num_post_smooth_steps; ++i){
        HipKernels<Scalar>::residual(A->nnzValues, A->colIndices, A->rowPointers, &x, &y, t.nnzValues, Ncur, 1, this->stream);
        HipKernels<Scalar>::vmul(jacobi_damping, d_invDiags[level].nnzValues, t.nnzValues, &x, Ncur, this->stream);
    }
}

// x = prec(y)
template <class Scalar, unsigned int block_size>
void rocsparseCPR<Scalar, block_size>::
apply_amg(const Scalar& y, 
          Scalar& x)
{
    HIP_CHECK(hipMemsetAsync(d_coarse_x.data()->nnzValues, 0, sizeof(Scalar) * this->Nb, this->stream));
    
    for (unsigned int i = 0; i < d_u.size(); ++i) {
        d_u[i].upload(this->Rmatrices[i].nnzValues.data(), this->stream);
    }
    
    HipKernels<Scalar>::residual(d_mat->nnzValues, d_mat->colIndices, d_mat->rowPointers, &x, &y, d_rs.data()->nnzValues, this->Nb, block_size, this->stream);
    
    HipKernels<Scalar>::full_to_pressure_restriction(d_rs.data()->nnzValues, d_weights.data()->nnzValues, d_coarse_y.data()->nnzValues, Nb, this->stream);
    
    amg_cycle_gpu(0, *(d_coarse_y.data()->nnzValues), *(d_coarse_x.data()->nnzValues));

    HipKernels<Scalar>::add_coarse_pressure_correction(d_coarse_x.data()->nnzValues, &x, this->pressure_idx, Nb, this->stream);
}

template <class Scalar, unsigned int block_size>
void rocsparseCPR<Scalar, block_size>::
apply(const Scalar &y,
      Scalar& x)
{
    Dune::Timer t_bilu0;

    bilu0->apply(y, x);

    if (verbosity >= 3) {
        HIP_CHECK(hipStreamSynchronize(this->stream));
        std::ostringstream out;
        out << "rocsparseCPR apply bilu0(): " << t_bilu0.stop() << " s";
        OpmLog::info(out.str());
    }

    Dune::Timer t_amg;
    apply_amg(y, x);
    if (verbosity >= 3) {
        std::ostringstream out;
        out << "rocsparseCPR apply amg(): " << t_amg.stop() << " s";
        OpmLog::info(out.str());
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
