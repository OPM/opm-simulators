/*
  Copyright 2023 Equinor ASA

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
#include <sstream>
#include <stdexcept>

#include <opm/common/OpmLog/OpmLog.hpp>
#include <opm/common/ErrorMacros.hpp>
#include <dune/common/timer.hh>

// WellContributions are included via the solver
// MultisegmentWellContribution includes the cuda runtime if found by CMake
// this leads to inclusion of both amd_hip_vector_types.h and vector_types.h
// which both define vector types like uchar2, short3 and double4.
// Restore the value (if defined) afterwards.
#ifdef HAVE_CUDA
#define HIP_HAVE_CUDA_DEFINED HAVE_CUDA
#endif

#undef HAVE_CUDA

#include <opm/simulators/linalg/bda/rocsparseSolverBackend.hpp>
#include <opm/simulators/linalg/bda/rocsparseWellContributions.hpp>

#include <opm/simulators/linalg/bda/BdaResult.hpp>

#include <hip/hip_runtime_api.h>
#include <hip/hip_version.h>

#ifdef HIP_HAVE_CUDA_DEFINED
#define HAVE_CUDA HIP_HAVE_CUDA_DEFINED
#undef HIP_HAVE_CUDA_DEFINED
#endif

#define HIP_CHECK(STAT)                                  \
    do {                                                 \
        const hipError_t stat = (STAT);                  \
        if(stat != hipSuccess)                           \
        {                                                \
            std::ostringstream oss;                      \
            oss << "rocsparseSolverBackend::hip ";       \
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
            oss << "rocsparseSolverBackend::rocsparse "; \
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
            oss << "rocsparseSolverBackend::rocblas ";   \
            oss << "error: " << stat;                    \
            OPM_THROW(std::logic_error, oss.str());      \
        }                                                \
    } while(0)

#include <cstddef>

namespace Opm
{
namespace Accelerator
{

using Opm::OpmLog;
using Dune::Timer;

template <unsigned int block_size>
rocsparseSolverBackend<block_size>::rocsparseSolverBackend(int verbosity_, int maxit_, double tolerance_, unsigned int platformID_, unsigned int deviceID_) : BdaSolver<block_size>(verbosity_, maxit_, tolerance_, platformID_, deviceID_) {
    int numDevices = 0;
    HIP_CHECK(hipGetDeviceCount(&numDevices));
    if (static_cast<int>(deviceID) >= numDevices) {
        OPM_THROW(std::runtime_error, "Error chosen too high HIP device ID");
    }
    HIP_CHECK(hipSetDevice(deviceID));

    ROCSPARSE_CHECK(rocsparse_create_handle(&handle));
    ROCBLAS_CHECK(rocblas_create_handle(&blas_handle));

    ROCSPARSE_CHECK(rocsparse_get_version(handle, &ver));
    ROCSPARSE_CHECK(rocsparse_get_git_rev(handle, rev));

    std::ostringstream out;
    out << "rocSPARSE version: " << ver / 100000 << "." << ver / 100 % 1000 << "."
        << ver % 100 << "-" << rev << "\n";
    OpmLog::info(out.str());

    HIP_CHECK(hipStreamCreate(&stream));
    ROCSPARSE_CHECK(rocsparse_set_stream(handle, stream));
    ROCBLAS_CHECK(rocblas_set_stream(blas_handle, stream));
}


template <unsigned int block_size>
rocsparseSolverBackend<block_size>::~rocsparseSolverBackend() {
    hipError_t hipstatus = hipStreamSynchronize(stream);
    if(hipstatus != hipSuccess){
        OpmLog::error("Could not synchronize with hipStream");
    }
    hipstatus = hipStreamDestroy(stream);
    if(hipstatus != hipSuccess){
        OpmLog::error("Could not destroy hipStream");
    }
    rocsparse_status status1 = rocsparse_destroy_handle(handle);
    if(status1 != rocsparse_status_success){
        OpmLog::error("Could not destroy rocsparse handle");
    }
    rocblas_status status2 = rocblas_destroy_handle(blas_handle);
    if(status2 != rocblas_status_success){
        OpmLog::error("Could not destroy rocblas handle");
    }
}


template <unsigned int block_size>
void rocsparseSolverBackend<block_size>::gpu_pbicgstab([[maybe_unused]] WellContributions& wellContribs,
                                                       BdaResult& res)
{
    float it = 0.5;
    double rho, rhop, beta, alpha, nalpha, omega, nomega, tmp1, tmp2;
    double norm, norm_0;
    double zero = 0.0;
    double one  = 1.0;
    double mone = -1.0;

    Timer t_total, t_prec(false), t_spmv(false), t_well(false), t_rest(false);

    // set stream here, the WellContributions object is destroyed every linear solve
    // the number of wells can change every linear solve
    if(wellContribs.getNumWells() > 0){
        static_cast<WellContributionsRocsparse&>(wellContribs).setStream(stream);
    }

// HIP_VERSION is defined as (HIP_VERSION_MAJOR * 10000000 + HIP_VERSION_MINOR * 100000 + HIP_VERSION_PATCH)
#if HIP_VERSION >= 50400000
    ROCSPARSE_CHECK(rocsparse_dbsrmv_ex(handle, dir, operation,
                                        Nb, Nb, nnzb, &one, descr_M,
                                        d_Avals, d_Arows, d_Acols, block_size,
                                        spmv_info, d_x, &zero, d_r));
#else
    ROCSPARSE_CHECK(rocsparse_dbsrmv(handle, dir, operation,
                                        Nb, Nb, nnzb, &one, descr_M,
                                        d_Avals, d_Arows, d_Acols, block_size,
                                        d_x, &zero, d_r));
#endif
    ROCBLAS_CHECK(rocblas_dscal(blas_handle, N, &mone, d_r, 1));
    ROCBLAS_CHECK(rocblas_daxpy(blas_handle, N, &one, d_b, 1, d_r, 1));
    ROCBLAS_CHECK(rocblas_dcopy(blas_handle, N, d_r, 1, d_rw, 1));
    ROCBLAS_CHECK(rocblas_dcopy(blas_handle, N, d_r, 1, d_p, 1));
    ROCBLAS_CHECK(rocblas_dnrm2(blas_handle, N, d_r, 1, &norm_0));

    if (verbosity >= 2) {
        std::ostringstream out;
        out << std::scientific << "rocsparseSolver initial norm: " << norm_0;
        OpmLog::info(out.str());
    }

    if (verbosity >= 3) {
        t_rest.start();
    }
    for (it = 0.5; it < maxit; it += 0.5) {
        rhop = rho;
        ROCBLAS_CHECK(rocblas_ddot(blas_handle, N, d_rw, 1, d_r, 1, &rho));

        if (it > 1) {
            beta = (rho / rhop) * (alpha / omega);
            nomega = -omega;
            ROCBLAS_CHECK(rocblas_daxpy(blas_handle, N, &nomega, d_v, 1, d_p, 1));
            ROCBLAS_CHECK(rocblas_dscal(blas_handle, N, &beta, d_p, 1));
            ROCBLAS_CHECK(rocblas_daxpy(blas_handle, N, &one, d_r, 1, d_p, 1));
        }
        if (verbosity >= 3) {
            HIP_CHECK(hipStreamSynchronize(stream));
            t_rest.stop();
            t_prec.start();
        }

        // apply ilu0
        ROCSPARSE_CHECK(rocsparse_dbsrsv_solve(handle, dir, \
                              operation, Nb, nnzbs_prec, &one, \
                              descr_L, d_Mvals, d_Mrows, d_Mcols, block_size, ilu_info, d_p, d_t, rocsparse_solve_policy_auto, d_buffer));
        ROCSPARSE_CHECK(rocsparse_dbsrsv_solve(handle, dir, \
                              operation, Nb, nnzbs_prec, &one, \
                              descr_U, d_Mvals, d_Mrows, d_Mcols, block_size, ilu_info, d_t, d_pw, rocsparse_solve_policy_auto, d_buffer));
        if (verbosity >= 3) {
            HIP_CHECK(hipStreamSynchronize(stream));
            t_prec.stop();
            t_spmv.start();
        }

        // spmv
#if HIP_VERSION >= 50400000
        ROCSPARSE_CHECK(rocsparse_dbsrmv_ex(handle, dir, operation,
                                            Nb, Nb, nnzb, &one, descr_M,
                                            d_Avals, d_Arows, d_Acols, block_size,
                                            spmv_info, d_pw, &zero, d_v));
#else
        ROCSPARSE_CHECK(rocsparse_dbsrmv(handle, dir, operation,
                                            Nb, Nb, nnzb, &one, descr_M,
                                            d_Avals, d_Arows, d_Acols, block_size,
                                            d_pw, &zero, d_v));
#endif
        if (verbosity >= 3) {
            HIP_CHECK(hipStreamSynchronize(stream));
            t_spmv.stop();
            t_well.start();
        }

        // apply wellContributions
        if(wellContribs.getNumWells() > 0){
            static_cast<WellContributionsRocsparse&>(wellContribs).apply(d_pw, d_v);
        }
        if (verbosity >= 3) {
            HIP_CHECK(hipStreamSynchronize(stream));
            t_well.stop();
            t_rest.start();
        }

        ROCBLAS_CHECK(rocblas_ddot(blas_handle, N, d_rw, 1, d_v, 1, &tmp1));
        alpha = rho / tmp1;
        nalpha = -alpha;
        ROCBLAS_CHECK(rocblas_daxpy(blas_handle, N, &nalpha, d_v, 1, d_r, 1));
        ROCBLAS_CHECK(rocblas_daxpy(blas_handle, N, &alpha, d_pw, 1, d_x, 1));
        ROCBLAS_CHECK(rocblas_dnrm2(blas_handle, N, d_r, 1, &norm));
        if (verbosity >= 3) {
            HIP_CHECK(hipStreamSynchronize(stream));
            t_rest.stop();
        }

        if (norm < tolerance * norm_0) {
            break;
        }

        it += 0.5;

        // apply ilu0
        if (verbosity >= 3) {
            t_prec.start();
        }
        ROCSPARSE_CHECK(rocsparse_dbsrsv_solve(handle, dir, \
                              operation, Nb, nnzbs_prec, &one, \
                              descr_L, d_Mvals, d_Mrows, d_Mcols, block_size, ilu_info, d_r, d_t, rocsparse_solve_policy_auto, d_buffer));
        ROCSPARSE_CHECK(rocsparse_dbsrsv_solve(handle, dir, \
                              operation, Nb, nnzbs_prec, &one, \
                              descr_U, d_Mvals, d_Mrows, d_Mcols, block_size, ilu_info, d_t, d_s, rocsparse_solve_policy_auto, d_buffer));
        if (verbosity >= 3) {
            HIP_CHECK(hipStreamSynchronize(stream));
            t_prec.stop();
            t_spmv.start();
        }

        // spmv
#if HIP_VERSION >= 50400000
        ROCSPARSE_CHECK(rocsparse_dbsrmv_ex(handle, dir, operation,
                                            Nb, Nb, nnzb, &one, descr_M,
                                            d_Avals, d_Arows, d_Acols, block_size,
                                            spmv_info, d_s, &zero, d_t));
#else
        ROCSPARSE_CHECK(rocsparse_dbsrmv(handle, dir, operation,
                                            Nb, Nb, nnzb, &one, descr_M,
                                            d_Avals, d_Arows, d_Acols, block_size,
                                            d_s, &zero, d_t));
#endif
        if(verbosity >= 3){
            HIP_CHECK(hipStreamSynchronize(stream));
            t_spmv.stop();
            t_well.start();
        }

        // apply wellContributions
        if(wellContribs.getNumWells() > 0){
            static_cast<WellContributionsRocsparse&>(wellContribs).apply(d_s, d_t);
        }
        if (verbosity >= 3) {
            HIP_CHECK(hipStreamSynchronize(stream));
            t_well.stop();
            t_rest.start();
        }

        ROCBLAS_CHECK(rocblas_ddot(blas_handle, N, d_t, 1, d_r, 1, &tmp1));
        ROCBLAS_CHECK(rocblas_ddot(blas_handle, N, d_t, 1, d_t, 1, &tmp2));
        omega = tmp1 / tmp2;
        nomega = -omega;
        ROCBLAS_CHECK(rocblas_daxpy(blas_handle, N, &omega, d_s, 1, d_x, 1));
        ROCBLAS_CHECK(rocblas_daxpy(blas_handle, N, &nomega, d_t, 1, d_r, 1));

        ROCBLAS_CHECK(rocblas_dnrm2(blas_handle, N, d_r, 1, &norm));
        if (verbosity >= 3) {
            HIP_CHECK(hipStreamSynchronize(stream));
            t_rest.stop();
        }

        if (norm < tolerance * norm_0) {
            break;
        }

        if (verbosity > 1) {
            std::ostringstream out;
            out << "it: " << it << std::scientific << ", norm: " << norm;
            OpmLog::info(out.str());
        }
    }

    res.iterations = std::min(it, (float)maxit);
    res.reduction = norm / norm_0;
    res.conv_rate  = static_cast<double>(pow(res.reduction, 1.0 / it));
    res.elapsed = t_total.stop();
    res.converged = (it != (maxit + 0.5));

    if (verbosity >= 1) {
        std::ostringstream out;
        out << "=== converged: " << res.converged << ", conv_rate: " << res.conv_rate << ", time: " << res.elapsed << \
            ", time per iteration: " << res.elapsed / it << ", iterations: " << it;
        OpmLog::info(out.str());
    }
    if (verbosity >= 3) {
        std::ostringstream out;
        out << "rocsparseSolver::prec_apply:  " << t_prec.elapsed() << " s\n";
        out << "rocsparseSolver::spmv:        " << t_spmv.elapsed() << " s\n";
        out << "rocsparseSolver::well:        " << t_well.elapsed() << " s\n";
        out << "rocsparseSolver::rest:        " << t_rest.elapsed() << " s\n";
        out << "rocsparseSolver::total_solve: " << res.elapsed << " s\n";
        OpmLog::info(out.str());
    }
}


template <unsigned int block_size>
void rocsparseSolverBackend<block_size>::initialize(std::shared_ptr<BlockedMatrix> matrix, std::shared_ptr<BlockedMatrix> jacMatrix) {
    this->Nb = matrix->Nb;
    this->N = Nb * block_size;
    this->nnzb = matrix->nnzbs;
    this->nnz = nnzb * block_size * block_size;
    nnzbs_prec = nnzb;

    if (jacMatrix) {
        useJacMatrix = true;
        nnzbs_prec = jacMatrix->nnzbs;
    }

    std::ostringstream out;
    out << "Initializing GPU, matrix size: " << Nb << " blockrows, nnzb: " << nnzb << "\n";
    if (useJacMatrix) {
        out << "Blocks in ILU matrix: " << jacMatrix->nnzbs << "\n";
    }
    out << "Maxit: " << maxit << std::scientific << ", tolerance: " << tolerance << "\n";
    out << "PlatformID: " << platformID << ", deviceID: " << deviceID << "\n";
    OpmLog::info(out.str());
    out.str("");
    out.clear();

    mat = matrix;
    jacMat = jacMatrix;

    HIP_CHECK(hipMalloc((void**)&d_r, sizeof(double) * N));
    HIP_CHECK(hipMalloc((void**)&d_rw, sizeof(double) * N));
    HIP_CHECK(hipMalloc((void**)&d_p, sizeof(double) * N));
    HIP_CHECK(hipMalloc((void**)&d_pw, sizeof(double) * N));
    HIP_CHECK(hipMalloc((void**)&d_s, sizeof(double) * N));
    HIP_CHECK(hipMalloc((void**)&d_t, sizeof(double) * N));
    HIP_CHECK(hipMalloc((void**)&d_v, sizeof(double) * N));

    HIP_CHECK(hipMalloc((void**)&d_Arows, sizeof(rocsparse_int) * (Nb + 1)));
    HIP_CHECK(hipMalloc((void**)&d_Acols, sizeof(rocsparse_int) * nnzb));
    HIP_CHECK(hipMalloc((void**)&d_Avals, sizeof(double) * nnz));
    HIP_CHECK(hipMalloc((void**)&d_x, sizeof(double) * N));
    HIP_CHECK(hipMalloc((void**)&d_b, sizeof(double) * N));

    if (useJacMatrix) {
        HIP_CHECK(hipMalloc((void**)&d_Mrows, sizeof(rocsparse_int) * (Nb + 1)));
        HIP_CHECK(hipMalloc((void**)&d_Mcols, sizeof(rocsparse_int) * nnzbs_prec));
        HIP_CHECK(hipMalloc((void**)&d_Mvals, sizeof(double) * nnzbs_prec * block_size * block_size));
    } else { // preconditioner matrix is same
        HIP_CHECK(hipMalloc((void**)&d_Mvals, sizeof(double) * nnzbs_prec * block_size * block_size));
        d_Mcols = d_Acols;
        d_Mrows = d_Arows;
    }

    initialized = true;
} // end initialize()

template <unsigned int block_size>
void rocsparseSolverBackend<block_size>::copy_system_to_gpu(double *b) {
    Timer t;

    HIP_CHECK(hipMemcpyAsync(d_Arows, mat->rowPointers, sizeof(rocsparse_int) * (Nb + 1), hipMemcpyHostToDevice, stream));
    HIP_CHECK(hipMemcpyAsync(d_Acols, mat->colIndices, sizeof(rocsparse_int) * nnzb, hipMemcpyHostToDevice, stream));
    HIP_CHECK(hipMemcpyAsync(d_Avals, mat->nnzValues, sizeof(double) * nnz, hipMemcpyHostToDevice, stream));
    if (useJacMatrix) {
        HIP_CHECK(hipMemcpyAsync(d_Mrows, jacMat->rowPointers, sizeof(rocsparse_int) * (Nb + 1), hipMemcpyHostToDevice, stream));
        HIP_CHECK(hipMemcpyAsync(d_Mcols, jacMat->colIndices, sizeof(rocsparse_int) * nnzbs_prec, hipMemcpyHostToDevice, stream));
        HIP_CHECK(hipMemcpyAsync(d_Mvals, jacMat->nnzValues, sizeof(double) * nnzbs_prec * block_size * block_size, hipMemcpyHostToDevice, stream));
    } else {
        HIP_CHECK(hipMemcpyAsync(d_Mvals, d_Avals, sizeof(double) * nnz, hipMemcpyDeviceToDevice, stream));
    }
    HIP_CHECK(hipMemsetAsync(d_x, 0, sizeof(double) * N, stream));
    HIP_CHECK(hipMemcpyAsync(d_b, b, sizeof(double) * N, hipMemcpyHostToDevice, stream));

    if (verbosity >= 3) {
        HIP_CHECK(hipStreamSynchronize(stream));
        std::ostringstream out;
        out << "rocsparseSolver::copy_system_to_gpu(): " << t.stop() << " s";
        OpmLog::info(out.str());
    }
} // end copy_system_to_gpu()

// don't copy rowpointers and colindices, they stay the same
template <unsigned int block_size>
void rocsparseSolverBackend<block_size>::update_system_on_gpu(double *b) {
    Timer t;

    HIP_CHECK(hipMemcpyAsync(d_Avals, mat->nnzValues, sizeof(double) * nnz, hipMemcpyHostToDevice, stream));
    if (useJacMatrix) {
        HIP_CHECK(hipMemcpyAsync(d_Mvals, jacMat->nnzValues, sizeof(double) * nnzbs_prec * block_size * block_size, hipMemcpyHostToDevice, stream));
    } else {
        HIP_CHECK(hipMemcpyAsync(d_Mvals, d_Avals, sizeof(double) * nnz, hipMemcpyDeviceToDevice, stream));
    }
    HIP_CHECK(hipMemsetAsync(d_x, 0, sizeof(double) * N, stream));
    HIP_CHECK(hipMemcpyAsync(d_b, b, sizeof(double) * N, hipMemcpyHostToDevice, stream));

    if (verbosity >= 3) {
        HIP_CHECK(hipStreamSynchronize(stream));
        std::ostringstream out;
        out << "rocsparseSolver::update_system_on_gpu(): " << t.stop() << " s";
        OpmLog::info(out.str());
    }
} // end update_system_on_gpu()

template <unsigned int block_size>
bool rocsparseSolverBackend<block_size>::analyze_matrix() {
    std::size_t d_bufferSize_M, d_bufferSize_L, d_bufferSize_U, d_bufferSize;
    Timer t;

    ROCSPARSE_CHECK(rocsparse_set_pointer_mode(handle, rocsparse_pointer_mode_host));

    ROCSPARSE_CHECK(rocsparse_create_mat_info(&ilu_info));
#if HIP_VERSION >= 50400000
    ROCSPARSE_CHECK(rocsparse_create_mat_info(&spmv_info));
#endif

    ROCSPARSE_CHECK(rocsparse_create_mat_descr(&descr_A));
    ROCSPARSE_CHECK(rocsparse_create_mat_descr(&descr_M));

    ROCSPARSE_CHECK(rocsparse_create_mat_descr(&descr_L));
    ROCSPARSE_CHECK(rocsparse_set_mat_fill_mode(descr_L, rocsparse_fill_mode_lower));
    ROCSPARSE_CHECK(rocsparse_set_mat_diag_type(descr_L, rocsparse_diag_type_unit));

    ROCSPARSE_CHECK(rocsparse_create_mat_descr(&descr_U));
    ROCSPARSE_CHECK(rocsparse_set_mat_fill_mode(descr_U, rocsparse_fill_mode_upper));
    ROCSPARSE_CHECK(rocsparse_set_mat_diag_type(descr_U, rocsparse_diag_type_non_unit));

    ROCSPARSE_CHECK(rocsparse_dbsrilu0_buffer_size(handle, dir, Nb, nnzbs_prec,
                                 descr_M, d_Mvals, d_Mrows, d_Mcols, block_size, ilu_info, &d_bufferSize_M));
    ROCSPARSE_CHECK(rocsparse_dbsrsv_buffer_size(handle, dir, operation, Nb, nnzbs_prec,
                               descr_L, d_Mvals, d_Mrows, d_Mcols, block_size, ilu_info, &d_bufferSize_L));
    ROCSPARSE_CHECK(rocsparse_dbsrsv_buffer_size(handle, dir, operation, Nb, nnzbs_prec,
                               descr_U, d_Mvals, d_Mrows, d_Mcols, block_size, ilu_info, &d_bufferSize_U));

    d_bufferSize = std::max(d_bufferSize_M, std::max(d_bufferSize_L, d_bufferSize_U));

    HIP_CHECK(hipMalloc((void**)&d_buffer, d_bufferSize));

    // analysis of ilu LU decomposition
    ROCSPARSE_CHECK(rocsparse_dbsrilu0_analysis(handle, dir, \
                               Nb, nnzbs_prec, descr_M, d_Mvals, d_Mrows, d_Mcols, \
                               block_size, ilu_info, rocsparse_analysis_policy_reuse, rocsparse_solve_policy_auto, d_buffer));

    int zero_position = 0;
    rocsparse_status status = rocsparse_bsrilu0_zero_pivot(handle, ilu_info, &zero_position);
    if (rocsparse_status_success != status) {
        printf("L has structural and/or numerical zero at L(%d,%d)\n", zero_position, zero_position);
        return false;
    }

    // analysis of ilu apply
    ROCSPARSE_CHECK(rocsparse_dbsrsv_analysis(handle, dir, operation, \
                             Nb, nnzbs_prec, descr_L, d_Mvals, d_Mrows, d_Mcols, \
                             block_size, ilu_info, rocsparse_analysis_policy_reuse, rocsparse_solve_policy_auto, d_buffer));
    ROCSPARSE_CHECK(rocsparse_dbsrsv_analysis(handle, dir, operation, \
                             Nb, nnzbs_prec, descr_U, d_Mvals, d_Mrows, d_Mcols, \
                             block_size, ilu_info, rocsparse_analysis_policy_reuse, rocsparse_solve_policy_auto, d_buffer));

#if HIP_VERSION >= 50400000
    ROCSPARSE_CHECK(rocsparse_dbsrmv_ex_analysis(handle, dir, operation,
        Nb, Nb, nnzb,
        descr_A, d_Avals, d_Arows, d_Acols,
        block_size, spmv_info));
#endif

    if (verbosity >= 3) {
        HIP_CHECK(hipStreamSynchronize(stream));
        std::ostringstream out;
        out << "rocsparseSolver::analyze_matrix(): " << t.stop() << " s";
        OpmLog::info(out.str());
    }

    analysis_done = true;

    return true;
} // end analyze_matrix()


template <unsigned int block_size>
bool rocsparseSolverBackend<block_size>::create_preconditioner() {
    Timer t;

    bool result = true;
    ROCSPARSE_CHECK(rocsparse_dbsrilu0(handle, dir, Nb, nnzbs_prec, descr_M,
                    d_Mvals, d_Mrows, d_Mcols, block_size, ilu_info, rocsparse_solve_policy_auto, d_buffer));

    // Check for zero pivot
    int zero_position = 0;
    rocsparse_status status = rocsparse_bsrilu0_zero_pivot(handle, ilu_info, &zero_position);
    if(rocsparse_status_success != status)
    {
        printf("L has structural and/or numerical zero at L(%d,%d)\n", zero_position, zero_position);
        return false;
    }

    if (verbosity >= 3) {
        HIP_CHECK(hipStreamSynchronize(stream));
        std::ostringstream out;
        out << "rocsparseSolver::create_preconditioner(): " << t.stop() << " s";
        OpmLog::info(out.str());
    }
    return result;
} // end create_preconditioner()


template <unsigned int block_size>
void rocsparseSolverBackend<block_size>::solve_system(WellContributions &wellContribs, BdaResult &res) {
    Timer t;

    // actually solve
    gpu_pbicgstab(wellContribs, res);

    if (verbosity >= 3) {
        HIP_CHECK(hipStreamSynchronize(stream));
        std::ostringstream out;
        out << "rocsparseSolver::solve_system(): " << t.stop() << " s";
        OpmLog::info(out.str());
    }

} // end solve_system()


// copy result to host memory
// caller must be sure that x is a valid array
template <unsigned int block_size>
void rocsparseSolverBackend<block_size>::get_result(double *x) {
    Timer t;

    HIP_CHECK(hipMemcpyAsync(x, d_x, sizeof(double) * N, hipMemcpyDeviceToHost, stream));
    HIP_CHECK(hipStreamSynchronize(stream)); // always wait, caller might want to use x immediately

    if (verbosity >= 3) {
        std::ostringstream out;
        out << "rocsparseSolver::get_result(): " << t.stop() << " s";
        OpmLog::info(out.str());
    }
} // end get_result()


template <unsigned int block_size>
SolverStatus rocsparseSolverBackend<block_size>::solve_system(std::shared_ptr<BlockedMatrix> matrix,
                                                              double *b,
                                                              std::shared_ptr<BlockedMatrix> jacMatrix,
                                                              WellContributions& wellContribs,
                                                              BdaResult &res)
{
    if (initialized == false) {
        initialize(matrix, jacMatrix);
        copy_system_to_gpu(b);
        if (analysis_done == false) {
            if (!analyze_matrix()) {
                return SolverStatus::BDA_SOLVER_ANALYSIS_FAILED;
            }
        }
        if (!create_preconditioner()) {
            return SolverStatus::BDA_SOLVER_CREATE_PRECONDITIONER_FAILED;
        }
    } else {
        update_system_on_gpu(b);
        if (!create_preconditioner()) {
            return SolverStatus::BDA_SOLVER_CREATE_PRECONDITIONER_FAILED;
        }
    }
    solve_system(wellContribs, res);

    return SolverStatus::BDA_SOLVER_SUCCESS;
}


#define INSTANTIATE_BDA_FUNCTIONS(n)                                        \
template rocsparseSolverBackend<n>::rocsparseSolverBackend(                 \
    int, int, double, unsigned int, unsigned int);

INSTANTIATE_BDA_FUNCTIONS(1);
INSTANTIATE_BDA_FUNCTIONS(2);
INSTANTIATE_BDA_FUNCTIONS(3);
INSTANTIATE_BDA_FUNCTIONS(4);
INSTANTIATE_BDA_FUNCTIONS(5);
INSTANTIATE_BDA_FUNCTIONS(6);

#undef INSTANTIATE_BDA_FUNCTIONS

} // namespace Accelerator
} // namespace Opm
