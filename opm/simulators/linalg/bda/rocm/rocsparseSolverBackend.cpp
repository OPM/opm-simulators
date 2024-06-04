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

#include <opm/simulators/linalg/bda/rocm/rocsparseSolverBackend.hpp>
#include <opm/simulators/linalg/bda/rocm/rocsparseWellContributions.hpp>

#include <opm/simulators/linalg/bda/BdaResult.hpp>

#include <opm/simulators/linalg/bda/Preconditioner.hpp>

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

namespace Opm::Accelerator {

using Dune::Timer;

template<class Scalar, unsigned int block_size>
rocsparseSolverBackend<Scalar,block_size>::
rocsparseSolverBackend(int verbosity_, int maxit_, Scalar tolerance_,
                       unsigned int platformID_, unsigned int deviceID_, std::string linsolver)
    : Base(verbosity_, maxit_, tolerance_, platformID_, deviceID_)
{
    int numDevices = 0;
    bool use_cpr, use_isai;
    
    if (linsolver.compare("ilu0") == 0) {
        use_cpr = false;
        use_isai = false;
    } else if (linsolver.compare("cpr_quasiimpes") == 0) {
        use_cpr = true; 
        use_isai = false;
    } else if (linsolver.compare("isai") == 0) {
        OPM_THROW(std::logic_error, "Error rocsparseSolver does not support --linerar-solver=isai");
    } else if (linsolver.compare("cpr_trueimpes") == 0) {
        OPM_THROW(std::logic_error, "Error rocsparseSolver does not support --linerar-solver=cpr_trueimpes");
    } else {
        OPM_THROW(std::logic_error, "Error unknown value for argument --linear-solver, " + linsolver);
    }
    
    HIP_CHECK(hipGetDeviceCount(&numDevices));
    if (static_cast<int>(deviceID) >= numDevices) {
        OPM_THROW(std::runtime_error, "Invalid HIP device ID");
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
    
    using PreconditionerType = typename Opm::Accelerator::PreconditionerType;
    if (use_cpr) {
        prec = rocsparsePreconditioner<Scalar, block_size>::create(PreconditionerType::CPR, verbosity);
    } else {
        prec = rocsparsePreconditioner<Scalar, block_size>::create(PreconditionerType::BILU0, verbosity);
    }
    
    prec->set_context(handle, dir, operation, stream);
}

template<class Scalar, unsigned int block_size>
rocsparseSolverBackend<Scalar,block_size>::~rocsparseSolverBackend()
{
    hipError_t hipstatus = hipStreamSynchronize(stream);
    if (hipstatus != hipSuccess) {
        OpmLog::error("Could not synchronize with hipStream");
    }
    hipstatus = hipStreamDestroy(stream);
    if (hipstatus != hipSuccess) {
        OpmLog::error("Could not destroy hipStream");
    }
    rocsparse_status status1 = rocsparse_destroy_handle(handle);
    if (status1 != rocsparse_status_success) {
        OpmLog::error("Could not destroy rocsparse handle");
    }
    rocblas_status status2 = rocblas_destroy_handle(blas_handle);
    if (status2 != rocblas_status_success) {
        OpmLog::error("Could not destroy rocblas handle");
    }
}

template<class Scalar, unsigned int block_size>
void rocsparseSolverBackend<Scalar,block_size>::
gpu_pbicgstab([[maybe_unused]] WellContributions<Scalar>& wellContribs,
              BdaResult& res)
{
    float it = 0.5;
    Scalar rho, rhop, beta, alpha, nalpha, omega, nomega, tmp1, tmp2;
    Scalar norm, norm_0;
    Scalar zero = 0.0;
    Scalar one  = 1.0;
    Scalar mone = -1.0;

    Timer t_total, t_prec(false), t_spmv(false), t_well(false), t_rest(false);

    // set stream here, the WellContributions object is destroyed every linear solve
    // the number of wells can change every linear solve
    if (wellContribs.getNumWells() > 0) {
        static_cast<WellContributionsRocsparse<Scalar>&>(wellContribs).setStream(stream);
    }

// HIP_VERSION is defined as (HIP_VERSION_MAJOR * 10000000 + HIP_VERSION_MINOR * 100000 + HIP_VERSION_PATCH)
#if HIP_VERSION >= 60000000
    ROCSPARSE_CHECK(rocsparse_dbsrmv(handle, dir, operation,
                                     Nb, Nb, nnzb, &one, descr_A,
                                     d_Avals, d_Arows, d_Acols, block_size,
                                     spmv_info, d_x, &zero, d_r));
#elif HIP_VERSION >= 50400000
    ROCSPARSE_CHECK(rocsparse_dbsrmv_ex(handle, dir, operation,
                                        Nb, Nb, nnzb, &one, descr_A,
                                        d_Avals, d_Arows, d_Acols, block_size,
                                        spmv_info, d_x, &zero, d_r));
#else
    ROCSPARSE_CHECK(rocsparse_dbsrmv(handle, dir, operation,
                                        Nb, Nb, nnzb, &one, descr_A,
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
        prec->apply(*d_p, *d_pw);
        
        if (verbosity >= 3) {
            HIP_CHECK(hipStreamSynchronize(stream));
            t_prec.stop();
            t_spmv.start();
        }

        // spmv
#if HIP_VERSION >= 60000000
        ROCSPARSE_CHECK(rocsparse_dbsrmv(handle, dir, operation,
                                          Nb, Nb, nnzb, &one, descr_A,
                                          d_Avals, d_Arows, d_Acols, block_size,
                                          spmv_info, d_pw, &zero, d_v));
#elif HIP_VERSION >= 50400000
        ROCSPARSE_CHECK(rocsparse_dbsrmv_ex(handle, dir, operation,
                                            Nb, Nb, nnzb, &one, descr_A,
                                            d_Avals, d_Arows, d_Acols, block_size,
                                            spmv_info, d_pw, &zero, d_v));
#else
        ROCSPARSE_CHECK(rocsparse_dbsrmv(handle, dir, operation,
                                            Nb, Nb, nnzb, &one, descr_A,
                                            d_Avals, d_Arows, d_Acols, block_size,
                                            d_pw, &zero, d_v));
#endif
        if (verbosity >= 3) {
            HIP_CHECK(hipStreamSynchronize(stream));
            t_spmv.stop();
            t_well.start();
        }

        // apply wellContributions
        if (wellContribs.getNumWells() > 0) {
            static_cast<WellContributionsRocsparse<Scalar>&>(wellContribs).apply(d_pw, d_v);
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
        
        prec->apply(*d_r, *d_s);
       
        if (verbosity >= 3) {
            HIP_CHECK(hipStreamSynchronize(stream));
            t_prec.stop();
            t_spmv.start();
        }

        // spmv
#if HIP_VERSION >= 60000000
        ROCSPARSE_CHECK(rocsparse_dbsrmv(handle, dir, operation,
                                         Nb, Nb, nnzb, &one, descr_A,
                                         d_Avals, d_Arows, d_Acols, block_size,
                                         spmv_info, d_s, &zero, d_t));
#elif HIP_VERSION >= 50400000
        ROCSPARSE_CHECK(rocsparse_dbsrmv_ex(handle, dir, operation,
                                            Nb, Nb, nnzb, &one, descr_A,
                                            d_Avals, d_Arows, d_Acols, block_size,
                                            spmv_info, d_s, &zero, d_t));
#else
        ROCSPARSE_CHECK(rocsparse_dbsrmv(handle, dir, operation,
                                            Nb, Nb, nnzb, &one, descr_A,
                                            d_Avals, d_Arows, d_Acols, block_size,
                                            d_s, &zero, d_t));
#endif
        if (verbosity >= 3) {
            HIP_CHECK(hipStreamSynchronize(stream));
            t_spmv.stop();
            t_well.start();
        }

        // apply wellContributions
        if (wellContribs.getNumWells() > 0) {
            static_cast<WellContributionsRocsparse<Scalar>&>(wellContribs).apply(d_s, d_t);
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
        out << "=== converged: " << res.converged
            << ", conv_rate: " << res.conv_rate
            << ", time: " << res.elapsed << \
            ", time per iteration: " << res.elapsed / it
            << ", iterations: " << it;
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

template<class Scalar, unsigned int block_size>
void rocsparseSolverBackend<Scalar,block_size>::
initialize(std::shared_ptr<BlockedMatrix<Scalar>> matrix,
           std::shared_ptr<BlockedMatrix<Scalar>> jacMatrix)
{
    this->Nb = matrix->Nb;
    this->N = this->Nb * block_size;
    this->nnzb = matrix->nnzbs;
    this->nnz = this->nnzb * block_size * block_size;

    if (jacMatrix) {
        prec->useJacMatrix = true;
        prec->jacMat = jacMatrix;
        prec->nnzbs_prec = jacMatrix->nnzbs;
    } else {
        prec->useJacMatrix = false;
        prec->jacMat = matrix;
        prec->nnzbs_prec = this->nnzb;
    }

    std::ostringstream out;
    out << "Initializing GPU, matrix size: "
        << Nb << " blockrows, nnzb: " << nnzb << "\n";
    if (useJacMatrix) {
        out << "Blocks in ILU matrix: " << jacMatrix->nnzbs << "\n";
    }
    out << "Maxit: " << maxit
        << std::scientific << ", tolerance: " << tolerance << "\n"
        << "PlatformID: " << platformID << ", deviceID: " << deviceID << "\n";
    OpmLog::info(out.str());
    out.str("");
    out.clear();

    mat = matrix;

    HIP_CHECK(hipMalloc((void**)&d_r, sizeof(Scalar) * N));
    HIP_CHECK(hipMalloc((void**)&d_rw, sizeof(Scalar) * N));
    HIP_CHECK(hipMalloc((void**)&d_p, sizeof(Scalar) * N));
    HIP_CHECK(hipMalloc((void**)&d_pw, sizeof(Scalar) * N));
    HIP_CHECK(hipMalloc((void**)&d_s, sizeof(Scalar) * N));
    HIP_CHECK(hipMalloc((void**)&d_t, sizeof(Scalar) * N));
    HIP_CHECK(hipMalloc((void**)&d_v, sizeof(Scalar) * N));

    HIP_CHECK(hipMalloc((void**)&d_Arows, sizeof(rocsparse_int) * (Nb + 1)));
    HIP_CHECK(hipMalloc((void**)&d_Acols, sizeof(rocsparse_int) * nnzb));
    HIP_CHECK(hipMalloc((void**)&d_Avals, sizeof(Scalar) * nnz));
    HIP_CHECK(hipMalloc((void**)&d_x, sizeof(Scalar) * N));
    HIP_CHECK(hipMalloc((void**)&d_b, sizeof(Scalar) * N));

    prec->initialize(matrix, jacMatrix, d_Arows, d_Acols);//TODO: do we need all parameters?

    initialized = true;
} // end initialize()

template<class Scalar, unsigned int block_size>
void rocsparseSolverBackend<Scalar,block_size>::
copy_system_to_gpu(Scalar *b)
{
    Timer t;

    HIP_CHECK(hipMemcpyAsync(d_Arows, mat->rowPointers,
                             sizeof(rocsparse_int) * (Nb + 1), 
                             hipMemcpyHostToDevice, stream));
    HIP_CHECK(hipMemcpyAsync(d_Acols, mat->colIndices, 
                             sizeof(rocsparse_int) * nnzb,
                             hipMemcpyHostToDevice, stream));
    HIP_CHECK(hipMemcpyAsync(d_Avals, mat->nnzValues,
                             sizeof(Scalar) * nnz,
                             hipMemcpyHostToDevice, stream));
    HIP_CHECK(hipMemsetAsync(d_x, 0, N * sizeof(Scalar), stream));
    HIP_CHECK(hipMemcpyAsync(d_b, b, N * sizeof(Scalar) * N,
                             hipMemcpyHostToDevice, stream));
    
    prec->copy_system_to_gpu(d_Avals); 

    if (verbosity >= 3) {
        HIP_CHECK(hipStreamSynchronize(stream));
        
        c_copy += t.stop();
        std::ostringstream out;
        out << "-----rocsparseSolver::copy_system_to_gpu(): " << t.elapsed() << " s\n";
        out << "---rocsparseSolver::cum copy: " << c_copy << " s";
        OpmLog::info(out.str());
    }
} // end copy_system_to_gpu()

// don't copy rowpointers and colindices, they stay the same
template<class Scalar, unsigned int block_size>
void rocsparseSolverBackend<Scalar,block_size>::
update_system_on_gpu(Scalar* b)
{
    Timer t;

    HIP_CHECK(hipMemcpyAsync(d_Avals, mat->nnzValues, sizeof(Scalar) * nnz,
                             hipMemcpyHostToDevice, stream));
    HIP_CHECK(hipMemsetAsync(d_x, 0, N * sizeof(Scalar), stream));
    HIP_CHECK(hipMemcpyAsync(d_b, b, N* sizeof(Scalar),
                             hipMemcpyHostToDevice, stream));
    
    prec->update_system_on_gpu(d_Avals);
    
    if (verbosity >= 3) {
        HIP_CHECK(hipStreamSynchronize(stream));

        c_copy += t.stop();
        std::ostringstream out;
        out << "-----rocsparseSolver::update_system_on_gpu(): " << t.elapsed() << " s\n";
        out << "---rocsparseSolver::cum copy: " << c_copy << " s";
        OpmLog::info(out.str());
    }
} // end update_system_on_gpu()

template<class Scalar, unsigned int block_size>
bool rocsparseSolverBackend<Scalar,block_size>::
analyze_matrix()
{
    Timer t;

    ROCSPARSE_CHECK(rocsparse_set_pointer_mode(handle, rocsparse_pointer_mode_host));

#if HIP_VERSION >= 50400000
    ROCSPARSE_CHECK(rocsparse_create_mat_info(&spmv_info));
#endif

    ROCSPARSE_CHECK(rocsparse_create_mat_descr(&descr_A));

#if HIP_VERSION >= 60000000
    ROCSPARSE_CHECK(rocsparse_dbsrmv_analysis(handle, dir, operation,
                                              Nb, Nb, nnzb,
                                              descr_A, d_Avals, d_Arows, d_Acols,
                                              block_size, spmv_info));
#elif HIP_VERSION >= 50400000
    ROCSPARSE_CHECK(rocsparse_dbsrmv_ex_analysis(handle, dir, operation,
        Nb, Nb, nnzb,
        descr_A, d_Avals, d_Arows, d_Acols,
        block_size, spmv_info));
#endif

    if(!prec->analyze_matrix(&*mat)) {
        std::ostringstream out;
        out << "Warning: rocsparseSolver matrix analysis failed!";
        OpmLog::info(out.str());
        return false;
    }
    
    if (verbosity >= 3) {
        HIP_CHECK(hipStreamSynchronize(stream));
        std::ostringstream out;
        out << "rocsparseSolver::analyze_matrix(): " << t.stop() << " s";
        OpmLog::info(out.str());
    }

    analysis_done = true;

    return true;
} // end analyze_matrix()

template<class Scalar, unsigned int block_size>
bool rocsparseSolverBackend<Scalar,block_size>::
create_preconditioner()
{
    return prec->create_preconditioner(&*mat);
} // end create_preconditioner()

template<class Scalar, unsigned int block_size>
void rocsparseSolverBackend<Scalar,block_size>::
solve_system(WellContributions<Scalar>& wellContribs, BdaResult& res)
{
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
template<class Scalar, unsigned int block_size>
void rocsparseSolverBackend<Scalar,block_size>::
get_result(Scalar* x)
{
    Timer t;

    HIP_CHECK(hipMemcpyAsync(x, d_x, sizeof(Scalar) * N,
                             hipMemcpyDeviceToHost, stream));
    HIP_CHECK(hipStreamSynchronize(stream)); // always wait, caller might want to use x immediately

    if (verbosity >= 3) {
        std::ostringstream out;
        out << "rocsparseSolver::get_result(): " << t.stop() << " s";
        OpmLog::info(out.str());
    }
} // end get_result()

template<class Scalar, unsigned int block_size>
SolverStatus rocsparseSolverBackend<Scalar,block_size>::
solve_system(std::shared_ptr<BlockedMatrix<Scalar>> matrix,
             Scalar* b,
             std::shared_ptr<BlockedMatrix<Scalar>> jacMatrix,
             WellContributions<Scalar>& wellContribs,
             BdaResult& res)
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

#define INSTANTIATE_TYPE(T)                     \
    template class rocsparseSolverBackend<T,1>; \
    template class rocsparseSolverBackend<T,2>; \
    template class rocsparseSolverBackend<T,3>; \
    template class rocsparseSolverBackend<T,4>; \
    template class rocsparseSolverBackend<T,5>; \
    template class rocsparseSolverBackend<T,6>;

INSTANTIATE_TYPE(double)

} // namespace Opm::Accelerator
