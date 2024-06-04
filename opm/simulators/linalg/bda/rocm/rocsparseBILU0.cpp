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
#include <algorithm>

#include <opm/common/OpmLog/OpmLog.hpp>
#include <opm/common/ErrorMacros.hpp>
#include <dune/common/timer.hh>

#include <opm/simulators/linalg/bda/rocm/rocsparseBILU0.hpp> 
#include <opm/simulators/linalg/bda/Reorder.hpp>

#include <sstream>

#include <iostream> //Razvan
#include <hip/hip_runtime_api.h>

#define HIP_CHECK(STAT)                                  \
    do {                                                 \
        const hipError_t stat = (STAT);                  \
        if(stat != hipSuccess)                           \
        {                                                \
            std::ostringstream oss;                      \
            oss << "rocsparseBILU0::hip ";       \
            oss << "error: " << hipGetErrorString(stat); \
            OPM_THROW(std::logic_error, oss.str());      \
        }                                                \
    } while(0)

#define ROCSPARSE_CHECK(STAT)                            \
    do {                                                 \
        const rocsparse_status stat = (STAT);            \
        if(stat != rocsparse_status_success)             \
        {                                                \
            std::ostringstream oss;                      \
            oss << "rocsparseBILU0::rocsparse "; \
            oss << "error: " << stat;                    \
            OPM_THROW(std::logic_error, oss.str());      \
        }                                                \
    } while(0)

#define ROCBLAS_CHECK(STAT)                              \
    do {                                                 \
        const rocblas_status stat = (STAT);              \
        if(stat != rocblas_status_success)               \
        {                                                \
            std::ostringstream oss;                      \
            oss << "rocsparseBILU0::rocblas ";   \
            oss << "error: " << stat;                    \
            OPM_THROW(std::logic_error, oss.str());      \
        }                                                \
    } while(0)
    
#if HAVE_OPENMP
#include <thread>
#include <omp.h>
extern std::shared_ptr<std::thread> copyThread;
#endif //HAVE_OPENMP

namespace Opm::Accelerator {

using Opm::OpmLog;
using Dune::Timer;

template <class Scalar, unsigned int block_size>
rocsparseBILU0<Scalar, block_size>::
rocsparseBILU0(int verbosity_) : 
    Base(verbosity_)
{
}

template <class Scalar, unsigned int block_size>
bool rocsparseBILU0<Scalar, block_size>::
initialize(std::shared_ptr<BlockedMatrix<Scalar>> matrix, std::shared_ptr<BlockedMatrix<Scalar>> jacMatrix, rocsparse_int *d_Arows, rocsparse_int *d_Acols) { 
    this->Nb = matrix->Nb;
    this->N = Nb * block_size;
    this->nnzb = matrix->nnzbs;
    this->nnz = nnzb * block_size * block_size;
    this->nnzbs_prec = this->nnzb;

    if (jacMatrix) {
        this->useJacMatrix = true;
        this->nnzbs_prec = jacMatrix->nnzbs;
        this->jacMat = jacMatrix;
    }

     HIP_CHECK(hipMalloc((void**)&d_t, sizeof(double) * this->N));

    if (this->useJacMatrix) {
        HIP_CHECK(hipMalloc((void**)&d_Mrows, sizeof(rocsparse_int) * (Nb + 1)));
        HIP_CHECK(hipMalloc((void**)&d_Mcols, sizeof(rocsparse_int) * this->nnzbs_prec));
        HIP_CHECK(hipMalloc((void**)&d_Mvals, sizeof(double) * this->nnzbs_prec * block_size * block_size));
    } else { // preconditioner matrix is same
        HIP_CHECK(hipMalloc((void**)&d_Mvals, sizeof(double) * this->nnzbs_prec * block_size * block_size));
        d_Mcols = d_Acols;
        d_Mrows = d_Arows;
    }
    
    return true;
} // end initialize()

template <class Scalar, unsigned int block_size>
bool rocsparseBILU0<Scalar, block_size>::
analyze_matrix(BlockedMatrix<Scalar> *mat) {
    return analyze_matrix(mat, &(*this->jacMat));
}


template <class Scalar, unsigned int block_size>
bool rocsparseBILU0<Scalar, block_size>::
analyze_matrix(BlockedMatrix<Scalar> *mat, BlockedMatrix<Scalar> *jacMat) {
    std::size_t d_bufferSize_M, d_bufferSize_L, d_bufferSize_U, d_bufferSize;
    Timer t;
    
#if HIP_VERSION >= 50400000
    ROCSPARSE_CHECK(rocsparse_create_mat_info(&spmv_info));
#endif

    ROCSPARSE_CHECK(rocsparse_create_mat_descr(&descr_M));

    ROCSPARSE_CHECK(rocsparse_create_mat_descr(&descr_L));
    ROCSPARSE_CHECK(rocsparse_set_mat_fill_mode(descr_L, rocsparse_fill_mode_lower));
    ROCSPARSE_CHECK(rocsparse_set_mat_diag_type(descr_L, rocsparse_diag_type_unit));

    ROCSPARSE_CHECK(rocsparse_create_mat_descr(&descr_U));
    ROCSPARSE_CHECK(rocsparse_set_mat_fill_mode(descr_U, rocsparse_fill_mode_upper));
    ROCSPARSE_CHECK(rocsparse_set_mat_diag_type(descr_U, rocsparse_diag_type_non_unit));
    
    ROCSPARSE_CHECK(rocsparse_dbsrilu0_buffer_size(this->handle, this->dir, Nb, this->nnzbs_prec, descr_M, d_Mvals, d_Mrows, d_Mcols, block_size, ilu_info, &d_bufferSize_M));
    
    ROCSPARSE_CHECK(rocsparse_dbsrsv_buffer_size(this->handle, this->dir, this->operation, Nb, this->nnzbs_prec,
                               descr_L, d_Mvals, d_Mrows, d_Mcols, block_size, ilu_info, &d_bufferSize_L));

    ROCSPARSE_CHECK(rocsparse_dbsrsv_buffer_size(this->handle, this->dir, this->operation, Nb, this->nnzbs_prec,
                               descr_U, d_Mvals, d_Mrows, d_Mcols, block_size, ilu_info, &d_bufferSize_U));

    d_bufferSize = std::max(d_bufferSize_M, std::max(d_bufferSize_L, d_bufferSize_U));

    HIP_CHECK(hipMalloc((void**)&d_buffer, d_bufferSize));

    // analysis of ilu LU decomposition
    ROCSPARSE_CHECK(rocsparse_dbsrilu0_analysis(this->handle, this->dir, \
                               Nb, this->nnzbs_prec, descr_M, d_Mvals, d_Mrows, d_Mcols, \
                               block_size, ilu_info, rocsparse_analysis_policy_reuse, rocsparse_solve_policy_auto, d_buffer));

    int zero_position = 0;
    rocsparse_status status = rocsparse_bsrilu0_zero_pivot(this->handle, ilu_info, &zero_position);
    if (rocsparse_status_success != status) {
        printf("L has structural and/or numerical zero at L(%d,%d)\n", zero_position, zero_position);
        return false;
    }

    // analysis of ilu apply
    ROCSPARSE_CHECK(rocsparse_dbsrsv_analysis(this->handle, this->dir, this->operation, \
                             Nb, this->nnzbs_prec, descr_L, d_Mvals, d_Mrows, d_Mcols, \
                             block_size, ilu_info, rocsparse_analysis_policy_reuse, rocsparse_solve_policy_auto, d_buffer));
    ROCSPARSE_CHECK(rocsparse_dbsrsv_analysis(this->handle, this->dir, this->operation, \
                             Nb, this->nnzbs_prec, descr_U, d_Mvals, d_Mrows, d_Mcols, \
                             block_size, ilu_info, rocsparse_analysis_policy_reuse, rocsparse_solve_policy_auto, d_buffer));

    if (verbosity >= 3) {
    	HIP_CHECK(hipStreamSynchronize(this->stream));
        std::ostringstream out;
        out << "rocsparseBILU0::analyze_matrix(): " << t.stop() << " s";
        OpmLog::info(out.str());
    }

    return true;
}

template <class Scalar, unsigned int block_size>
bool rocsparseBILU0<Scalar, block_size>::
create_preconditioner(BlockedMatrix<Scalar> *mat) {
    return create_preconditioner(mat, &*this->jacMat);
}

template <class Scalar, unsigned int block_size>
bool rocsparseBILU0<Scalar, block_size>::
create_preconditioner(BlockedMatrix<Scalar> *mat, BlockedMatrix<Scalar> *jacMat) {
    Timer t;
    bool result = true;
    
    ROCSPARSE_CHECK(rocsparse_dbsrilu0(this->handle, this->dir, Nb, this->nnzbs_prec, descr_M, d_Mvals, d_Mrows, d_Mcols, block_size, ilu_info, rocsparse_solve_policy_auto, d_buffer));

    // Check for zero pivot
    int zero_position = 0;
    rocsparse_status status = rocsparse_bsrilu0_zero_pivot(this->handle, ilu_info, &zero_position);
    if(rocsparse_status_success != status)
    {
        printf("L has structural and/or numerical zero at L(%d,%d)\n", zero_position, zero_position);
        return false;
    }

    if (verbosity >= 3) {
        HIP_CHECK(hipStreamSynchronize(this->stream));
        std::ostringstream out;
        out << "rocsparseBILU0::create_preconditioner(): " << t.stop() << " s";
        OpmLog::info(out.str());
    }
    return result;
} // end create_preconditioner()

template <class Scalar, unsigned int block_size>
void rocsparseBILU0<Scalar, block_size>::
copy_system_to_gpu(double *d_Avals) {
    Timer t;
    
    if (this->useJacMatrix) {
#if HAVE_OPENMP
        if (omp_get_max_threads() > 1) {
           copyThread->join();
        }
#endif
        HIP_CHECK(hipMemcpyAsync(d_Mrows, this->jacMat->rowPointers, sizeof(rocsparse_int) * (Nb + 1), hipMemcpyHostToDevice, this->stream));
        HIP_CHECK(hipMemcpyAsync(d_Mcols, this->jacMat->colIndices, sizeof(rocsparse_int) * this->nnzbs_prec, hipMemcpyHostToDevice, this->stream));
        HIP_CHECK(hipMemcpyAsync(d_Mvals, this->jacMat->nnzValues, sizeof(double) * this->nnzbs_prec * block_size * block_size, hipMemcpyHostToDevice, this->stream));
    } else {
        HIP_CHECK(hipMemcpyAsync(d_Mvals, d_Avals, sizeof(double) * nnz, hipMemcpyDeviceToDevice, this->stream));
    }

    if (verbosity >= 3) {
        HIP_CHECK(hipStreamSynchronize(this->stream));
        std::ostringstream out;
        out << "rocsparseBILU0::copy_system_to_gpu(): " << t.stop() << " s";
        OpmLog::info(out.str());
    }
} // end copy_system_to_gpu()

// don't copy rowpointers and colindices, they stay the same
template <class Scalar, unsigned int block_size>
void rocsparseBILU0<Scalar, block_size>::
update_system_on_gpu(double *d_Avals) {
    Timer t;

    if (this->useJacMatrix) {
#if HAVE_OPENMP
        if (omp_get_max_threads() > 1) {
           copyThread->join();
        }
#endif
        HIP_CHECK(hipMemcpyAsync(d_Mvals, this->jacMat->nnzValues, sizeof(double) * this->nnzbs_prec * block_size * block_size, hipMemcpyHostToDevice, this->stream));
    } else {
        HIP_CHECK(hipMemcpyAsync(d_Mvals, d_Avals, sizeof(double) * nnz, hipMemcpyDeviceToDevice, this->stream));
    }

    if (verbosity >= 3) {
        HIP_CHECK(hipStreamSynchronize(this->stream));
        std::ostringstream out;
        out << "rocsparseSolver::update_system_on_gpu(): " << t.stop() << " s";
        OpmLog::info(out.str());
    }
} // end update_system_on_gpu()

template <class Scalar, unsigned int block_size>
void rocsparseBILU0<Scalar, block_size>::
apply(double& y, double& x) {
    double zero = 0.0;
    double one  = 1.0;

    Timer t_apply;

    ROCSPARSE_CHECK(rocsparse_dbsrsv_solve(this->handle, this->dir, \
                              this->operation, Nb, this->nnzbs_prec, &one, \
                              descr_L, d_Mvals, d_Mrows, d_Mcols, block_size, ilu_info, &y, d_t, rocsparse_solve_policy_auto, d_buffer));

    ROCSPARSE_CHECK(rocsparse_dbsrsv_solve(this->handle, this->dir, \
                              this->operation, Nb, this->nnzbs_prec, &one, \
                              descr_U, d_Mvals, d_Mrows, d_Mcols, block_size, ilu_info, d_t, &x, rocsparse_solve_policy_auto, d_buffer));
        
    if (verbosity >= 3) {
        std::ostringstream out;
        HIP_CHECK(hipStreamSynchronize(this->stream));
        out << "rocsparseBILU0 apply: " << t_apply.stop() << " s";
        OpmLog::info(out.str());
    }
}

#define INSTANCE_TYPE(T)                \
    template class rocsparseBILU0<T,1>; \
    template class rocsparseBILU0<T,2>; \
    template class rocsparseBILU0<T,3>; \
    template class rocsparseBILU0<T,4>; \
    template class rocsparseBILU0<T,5>; \
    template class rocsparseBILU0<T,6>;

INSTANCE_TYPE(double)
} // namespace Opm
