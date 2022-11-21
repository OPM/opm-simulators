/*
  Copyright 2019 Equinor ASA

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

#include <cuda_runtime.h>
#include <sstream>

#include <opm/common/OpmLog/OpmLog.hpp>
#include <dune/common/timer.hh>

#include <opm/simulators/linalg/bda/cuda/cusparseSolverBackend.hpp>
#include <opm/simulators/linalg/bda/cuda/cuWellContributions.hpp>
#include <opm/simulators/linalg/bda/BdaResult.hpp>
#include <opm/simulators/linalg/bda/cuda/cuda_header.hpp>

#include "cublas_v2.h"
#include "cusparse_v2.h"
// For more information about cusparse, check https://docs.nvidia.com/cuda/cusparse/index.html

// iff true, the nonzeroes of the matrix are copied row-by-row into a contiguous, pinned memory array, then a single GPU memcpy is done
// otherwise, the nonzeroes of the matrix are assumed to be in a contiguous array, and a single GPU memcpy is enough
#define COPY_ROW_BY_ROW 0

namespace Opm
{
namespace Accelerator
{

using Opm::OpmLog;
using Dune::Timer;

const cusparseSolvePolicy_t policy = CUSPARSE_SOLVE_POLICY_USE_LEVEL;
const cusparseOperation_t operation  = CUSPARSE_OPERATION_NON_TRANSPOSE;
const cusparseDirection_t order = CUSPARSE_DIRECTION_ROW;


template <unsigned int block_size>
cusparseSolverBackend<block_size>::cusparseSolverBackend(int verbosity_, int maxit_, double tolerance_, unsigned int deviceID_) : BdaSolver<block_size>(verbosity_, maxit_, tolerance_, deviceID_) {

    // initialize CUDA device, stream and libraries
    cudaSetDevice(deviceID);
    cudaCheckLastError("Could not get device");
    struct cudaDeviceProp props;
    cudaGetDeviceProperties(&props, deviceID);
    cudaCheckLastError("Could not get device properties");
    std::ostringstream out;
    out << "Name GPU: " << props.name << ", Compute Capability: " << props.major << "." << props.minor;
    OpmLog::info(out.str());

    cudaStreamCreate(&stream);
    cudaCheckLastError("Could not create stream");

    cublasCreate(&cublasHandle);
    cudaCheckLastError("Could not create cublasHandle");
    cusparseCreate(&cusparseHandle);
    cudaCheckLastError("Could not create cusparseHandle");

    cublasSetStream(cublasHandle, stream);
    cudaCheckLastError("Could not set stream to cublas");
    cusparseSetStream(cusparseHandle, stream);
    cudaCheckLastError("Could not set stream to cusparse");
}

template <unsigned int block_size>
cusparseSolverBackend<block_size>::~cusparseSolverBackend() {
    finalize();
}

template <unsigned int block_size>
void cusparseSolverBackend<block_size>::gpu_pbicgstab(WellContributions& wellContribs, BdaResult& res) {
    Timer t_total, t_prec(false), t_spmv(false), t_well(false), t_rest(false);
    int n = N;
    double rho = 1.0, rhop;
    double alpha, nalpha, beta;
    double omega, nomega, tmp1, tmp2;
    double norm, norm_0;
    double zero = 0.0;
    double one  = 1.0;
    double mone = -1.0;
    float it;

    if (wellContribs.getNumWells() > 0) {
        static_cast<WellContributionsCuda&>(wellContribs).setCudaStream(stream);
    }

    cusparseDbsrmv(cusparseHandle, order, operation, Nb, Nb, nnzb, &one, descr_M, d_bVals, d_bRows, d_bCols, block_size, d_x, &zero, d_r);

    cublasDscal(cublasHandle, n, &mone, d_r, 1);
    cublasDaxpy(cublasHandle, n, &one, d_b, 1, d_r, 1);
    cublasDcopy(cublasHandle, n, d_r, 1, d_rw, 1);
    cublasDcopy(cublasHandle, n, d_r, 1, d_p, 1);
    cublasDnrm2(cublasHandle, n, d_r, 1, &norm_0);

    if (verbosity > 1) {
        std::ostringstream out;
        out << std::scientific << "cusparseSolver initial norm: " << norm_0;
        OpmLog::info(out.str());
    }

    for (it = 0.5; it < maxit; it += 0.5) {
        rhop = rho;
        cublasDdot(cublasHandle, n, d_rw, 1, d_r, 1, &rho);

        if (it > 1) {
            beta = (rho / rhop) * (alpha / omega);
            nomega = -omega;
            cublasDaxpy(cublasHandle, n, &nomega, d_v, 1, d_p, 1);
            cublasDscal(cublasHandle, n, &beta, d_p, 1);
            cublasDaxpy(cublasHandle, n, &one, d_r, 1, d_p, 1);
        }

        // apply ilu0
        cusparseDbsrsv2_solve(cusparseHandle, order, \
                              operation, Nb, nnzbs_prec, &one, \
                              descr_L, d_mVals, d_mRows, d_mCols, block_size, info_L, d_p, d_t, policy, d_buffer);
        cusparseDbsrsv2_solve(cusparseHandle, order, \
                              operation, Nb, nnzbs_prec, &one, \
                              descr_U, d_mVals, d_mRows, d_mCols, block_size, info_U, d_t, d_pw, policy, d_buffer);

        // spmv
        cusparseDbsrmv(cusparseHandle, order, \
                       operation, Nb, Nb, nnzb, \
                       &one, descr_M, d_bVals, d_bRows, d_bCols, block_size, d_pw, &zero, d_v);

        // apply wellContributions
        if (wellContribs.getNumWells() > 0) {
            static_cast<WellContributionsCuda&>(wellContribs).apply(d_pw, d_v);
        }

        cublasDdot(cublasHandle, n, d_rw, 1, d_v, 1, &tmp1);
        alpha = rho / tmp1;
        nalpha = -alpha;
        cublasDaxpy(cublasHandle, n, &nalpha, d_v, 1, d_r, 1);
        cublasDaxpy(cublasHandle, n, &alpha, d_pw, 1, d_x, 1);
        cublasDnrm2(cublasHandle, n, d_r, 1, &norm);

        if (norm < tolerance * norm_0) {
            break;
        }

        it += 0.5;

        // apply ilu0
        cusparseDbsrsv2_solve(cusparseHandle, order, \
                              operation, Nb, nnzbs_prec, &one, \
                              descr_L, d_mVals, d_mRows, d_mCols, block_size, info_L, d_r, d_t, policy, d_buffer);
        cusparseDbsrsv2_solve(cusparseHandle, order, \
                              operation, Nb, nnzbs_prec, &one, \
                              descr_U, d_mVals, d_mRows, d_mCols, block_size, info_U, d_t, d_s, policy, d_buffer);

        // spmv
        cusparseDbsrmv(cusparseHandle, order, \
                       operation, Nb, Nb, nnzb, &one, descr_M, \
                       d_bVals, d_bRows, d_bCols, block_size, d_s, &zero, d_t);

        // apply wellContributions
        if (wellContribs.getNumWells() > 0) {
            static_cast<WellContributionsCuda&>(wellContribs).apply(d_s, d_t);
        }

        cublasDdot(cublasHandle, n, d_t, 1, d_r, 1, &tmp1);
        cublasDdot(cublasHandle, n, d_t, 1, d_t, 1, &tmp2);
        omega = tmp1 / tmp2;
        nomega = -omega;
        cublasDaxpy(cublasHandle, n, &omega, d_s, 1, d_x, 1);
        cublasDaxpy(cublasHandle, n, &nomega, d_t, 1, d_r, 1);

        cublasDnrm2(cublasHandle, n, d_r, 1, &norm);


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

    if (verbosity > 0) {
        std::ostringstream out;
        out << "=== converged: " << res.converged << ", conv_rate: " << res.conv_rate << ", time: " << res.elapsed << \
            ", time per iteration: " << res.elapsed / it << ", iterations: " << it;
        OpmLog::info(out.str());
    }
}


template <unsigned int block_size>
void cusparseSolverBackend<block_size>::initialize(std::shared_ptr<BlockedMatrix> matrix, std::shared_ptr<BlockedMatrix> jacMatrix) {
    this->Nb = matrix->Nb;
    this->N = Nb * block_size;
    this->nnzb = matrix->nnzbs;
    this->nnz = nnzb * block_size * block_size;

    if (jacMatrix) {
        useJacMatrix = true;
        nnzbs_prec = jacMatrix->nnzbs;
    } else {
        nnzbs_prec = nnzb;
    }

    std::ostringstream out;
    out << "Initializing GPU, matrix size: " << Nb << " blockrows, nnz: " << nnzb << " blocks\n";
    if (useJacMatrix) {
        out << "Blocks in ILU matrix: " << nnzbs_prec << "\n";
    }
    out << "Maxit: " << maxit << std::scientific << ", tolerance: " << tolerance << "\n";
    OpmLog::info(out.str());

    cudaMalloc((void**)&d_x, sizeof(double) * N);
    cudaMalloc((void**)&d_b, sizeof(double) * N);
    cudaMalloc((void**)&d_r, sizeof(double) * N);
    cudaMalloc((void**)&d_rw, sizeof(double) * N);
    cudaMalloc((void**)&d_p, sizeof(double) * N);
    cudaMalloc((void**)&d_pw, sizeof(double) * N);
    cudaMalloc((void**)&d_s, sizeof(double) * N);
    cudaMalloc((void**)&d_t, sizeof(double) * N);
    cudaMalloc((void**)&d_v, sizeof(double) * N);
    cudaMalloc((void**)&d_bVals, sizeof(double) * nnz);
    cudaMalloc((void**)&d_bCols, sizeof(int) * nnzb);
    cudaMalloc((void**)&d_bRows, sizeof(int) * (Nb + 1));
    if (useJacMatrix) {
        cudaMalloc((void**)&d_mVals, sizeof(double) * nnzbs_prec * block_size * block_size);
        cudaMalloc((void**)&d_mCols, sizeof(int) * nnzbs_prec);
        cudaMalloc((void**)&d_mRows, sizeof(int) * (Nb + 1));
    } else {
        cudaMalloc((void**)&d_mVals, sizeof(double) * nnz);
        d_mCols = d_bCols;
        d_mRows = d_bRows;
    }
    cudaCheckLastError("Could not allocate enough memory on GPU");

#if COPY_ROW_BY_ROW
    cudaMallocHost((void**)&vals_contiguous, sizeof(double) * nnz);
    cudaCheckLastError("Could not allocate pinned memory");
#endif

    initialized = true;
} // end initialize()

template <unsigned int block_size>
void cusparseSolverBackend<block_size>::finalize() {
    if (initialized) {
        cudaFree(d_x);
        cudaFree(d_b);
        cudaFree(d_r);
        cudaFree(d_rw);
        cudaFree(d_p);
        cudaFree(d_pw);
        cudaFree(d_s);
        cudaFree(d_t);
        cudaFree(d_v);
        cudaFree(d_mVals);
        if (useJacMatrix) {
            cudaFree(d_mCols);
            cudaFree(d_mRows);
        }
        cudaFree(d_bVals);
        cudaFree(d_bCols);
        cudaFree(d_bRows);
        cudaFree(d_buffer);
        cusparseDestroyBsrilu02Info(info_M);
        cusparseDestroyBsrsv2Info(info_L);
        cusparseDestroyBsrsv2Info(info_U);
        cusparseDestroyMatDescr(descr_B);
        cusparseDestroyMatDescr(descr_M);
        cusparseDestroyMatDescr(descr_L);
        cusparseDestroyMatDescr(descr_U);
        cusparseDestroy(cusparseHandle);
        cublasDestroy(cublasHandle);
#if COPY_ROW_BY_ROW
        cudaFreeHost(vals_contiguous);
#endif
        cudaStreamDestroy(stream);
    }
} // end finalize()


template <unsigned int block_size>
void cusparseSolverBackend<block_size>::copy_system_to_gpu(std::shared_ptr<BlockedMatrix> matrix, double *b, std::shared_ptr<BlockedMatrix> jacMatrix) {
    Timer t;

#if COPY_ROW_BY_ROW
    int sum = 0;
    for (int i = 0; i < Nb; ++i) {
        int size_row = matrix->rowPointers[i + 1] - matrix->rowPointers[i];
        memcpy(vals_contiguous + sum, matrix->nnzValues + sum, size_row * sizeof(double) * block_size * block_size);
        sum += size_row * block_size * block_size;
    }
    cudaMemcpyAsync(d_bVals, vals_contiguous, nnz * sizeof(double), cudaMemcpyHostToDevice, stream);
#else
    cudaMemcpyAsync(d_bVals, matrix->nnzValues, nnz * sizeof(double), cudaMemcpyHostToDevice, stream);
    if (useJacMatrix) {
        cudaMemcpyAsync(d_mVals, jacMatrix->nnzValues, nnzbs_prec * block_size * block_size * sizeof(double), cudaMemcpyHostToDevice, stream);
    } else {
        cudaMemcpyAsync(d_mVals, d_bVals, nnz  * sizeof(double), cudaMemcpyDeviceToDevice, stream);
    }
#endif

    cudaMemcpyAsync(d_bCols, matrix->colIndices, nnzb * sizeof(int), cudaMemcpyHostToDevice, stream);
    cudaMemcpyAsync(d_bRows, matrix->rowPointers, (Nb + 1) * sizeof(int), cudaMemcpyHostToDevice, stream);
    if (useJacMatrix) {
        cudaMemcpyAsync(d_mCols, jacMatrix->colIndices, nnzbs_prec * sizeof(int), cudaMemcpyHostToDevice, stream);
        cudaMemcpyAsync(d_mRows, jacMatrix->rowPointers, (Nb + 1) * sizeof(int), cudaMemcpyHostToDevice, stream);
    }
    cudaMemcpyAsync(d_b, b, N * sizeof(double), cudaMemcpyHostToDevice, stream);
    cudaMemsetAsync(d_x, 0, sizeof(double) * N, stream);

    if (verbosity > 2) {
        cudaStreamSynchronize(stream);
        std::ostringstream out;
        out << "cusparseSolver::copy_system_to_gpu(): " << t.stop() << " s";
        OpmLog::info(out.str());
    }
} // end copy_system_to_gpu()


// don't copy rowpointers and colindices, they stay the same
template <unsigned int block_size>
void cusparseSolverBackend<block_size>::update_system_on_gpu(std::shared_ptr<BlockedMatrix> matrix, double *b, std::shared_ptr<BlockedMatrix> jacMatrix) {
    Timer t;

#if COPY_ROW_BY_ROW
    int sum = 0;
    for (int i = 0; i < Nb; ++i) {
        int size_row = matrix->rowPointers[i + 1] - matrix->rowPointers[i];
        memcpy(vals_contiguous + sum, matrix->nnzValues + sum, size_row * sizeof(double) * block_size * block_size);
        sum += size_row * block_size * block_size;
    }
    cudaMemcpyAsync(d_bVals, vals_contiguous, nnz * sizeof(double), cudaMemcpyHostToDevice, stream);
#else
    cudaMemcpyAsync(d_bVals, matrix->nnzValues, nnz * sizeof(double), cudaMemcpyHostToDevice, stream);
    if (useJacMatrix) {
        cudaMemcpyAsync(d_mVals, jacMatrix->nnzValues, nnzbs_prec * block_size * block_size * sizeof(double), cudaMemcpyHostToDevice, stream);
    } else {
        cudaMemcpyAsync(d_mVals, d_bVals, nnz  * sizeof(double), cudaMemcpyDeviceToDevice, stream);
    }
#endif

    cudaMemcpyAsync(d_b, b, N * sizeof(double), cudaMemcpyHostToDevice, stream);
    cudaMemsetAsync(d_x, 0, sizeof(double) * N, stream);

    if (verbosity > 2) {
        cudaStreamSynchronize(stream);
        std::ostringstream out;
        out << "cusparseSolver::update_system_on_gpu(): " << t.stop() << " s";
        OpmLog::info(out.str());
    }
} // end update_system_on_gpu()


template <unsigned int block_size>
bool cusparseSolverBackend<block_size>::analyse_matrix() {

    int d_bufferSize_M, d_bufferSize_L, d_bufferSize_U, d_bufferSize;
    Timer t;

    cusparseCreateMatDescr(&descr_B);
    cusparseCreateMatDescr(&descr_M);
    cusparseSetMatType(descr_B, CUSPARSE_MATRIX_TYPE_GENERAL);
    cusparseSetMatType(descr_M, CUSPARSE_MATRIX_TYPE_GENERAL);
    const cusparseIndexBase_t base_type = CUSPARSE_INDEX_BASE_ZERO;     // matrices from Flow are base0

    cusparseSetMatIndexBase(descr_B, base_type);
    cusparseSetMatIndexBase(descr_M, base_type);

    cusparseCreateMatDescr(&descr_L);
    cusparseSetMatIndexBase(descr_L, base_type);
    cusparseSetMatType(descr_L, CUSPARSE_MATRIX_TYPE_GENERAL);
    cusparseSetMatFillMode(descr_L, CUSPARSE_FILL_MODE_LOWER);
    cusparseSetMatDiagType(descr_L, CUSPARSE_DIAG_TYPE_UNIT);

    cusparseCreateMatDescr(&descr_U);
    cusparseSetMatIndexBase(descr_U, base_type);
    cusparseSetMatType(descr_U, CUSPARSE_MATRIX_TYPE_GENERAL);
    cusparseSetMatFillMode(descr_U, CUSPARSE_FILL_MODE_UPPER);
    cusparseSetMatDiagType(descr_U, CUSPARSE_DIAG_TYPE_NON_UNIT);
    cudaCheckLastError("Could not initialize matrix descriptions");

    cusparseCreateBsrilu02Info(&info_M);
    cusparseCreateBsrsv2Info(&info_L);
    cusparseCreateBsrsv2Info(&info_U);
    cudaCheckLastError("Could not create analysis info");

    cusparseDbsrilu02_bufferSize(cusparseHandle, order, Nb, nnzbs_prec,
                                 descr_M, d_mVals, d_mRows, d_mCols, block_size, info_M, &d_bufferSize_M);
    cusparseDbsrsv2_bufferSize(cusparseHandle, order, operation, Nb, nnzbs_prec,
                               descr_L, d_mVals, d_mRows, d_mCols, block_size, info_L, &d_bufferSize_L);
    cusparseDbsrsv2_bufferSize(cusparseHandle, order, operation, Nb, nnzbs_prec,
                               descr_U, d_mVals, d_mRows, d_mCols, block_size, info_U, &d_bufferSize_U);
    cudaCheckLastError();
    d_bufferSize = std::max(d_bufferSize_M, std::max(d_bufferSize_L, d_bufferSize_U));

    cudaMalloc((void**)&d_buffer, d_bufferSize);

    // analysis of ilu LU decomposition
    cusparseDbsrilu02_analysis(cusparseHandle, order, \
                               Nb, nnzbs_prec, descr_B, d_mVals, d_mRows, d_mCols, \
                               block_size, info_M, policy, d_buffer);

    int structural_zero;
    cusparseStatus_t status = cusparseXbsrilu02_zeroPivot(cusparseHandle, info_M, &structural_zero);
    if (CUSPARSE_STATUS_ZERO_PIVOT == status) {
        return false;
    }

    // analysis of ilu apply
    cusparseDbsrsv2_analysis(cusparseHandle, order, operation, \
                             Nb, nnzbs_prec, descr_L, d_mVals, d_mRows, d_mCols, \
                             block_size, info_L, policy, d_buffer);

    cusparseDbsrsv2_analysis(cusparseHandle, order, operation, \
                             Nb, nnzbs_prec, descr_U, d_mVals, d_mRows, d_mCols, \
                             block_size, info_U, policy, d_buffer);
    cudaCheckLastError("Could not analyse level information");

    if (verbosity > 2) {
        cudaStreamSynchronize(stream);
        std::ostringstream out;
        out << "cusparseSolver::analyse_matrix(): " << t.stop() << " s";
        OpmLog::info(out.str());
    }

    analysis_done = true;

    return true;
} // end analyse_matrix()

template <unsigned int block_size>
bool cusparseSolverBackend<block_size>::create_preconditioner() {
    Timer t;

    cusparseDbsrilu02(cusparseHandle, order, \
                      Nb, nnzbs_prec, descr_M, d_mVals, d_mRows, d_mCols, \
                      block_size, info_M, policy, d_buffer);
    cudaCheckLastError("Could not perform ilu decomposition");

    int structural_zero;
    // cusparseXbsrilu02_zeroPivot() calls cudaDeviceSynchronize()
    cusparseStatus_t status = cusparseXbsrilu02_zeroPivot(cusparseHandle, info_M, &structural_zero);
    if (CUSPARSE_STATUS_ZERO_PIVOT == status) {
        return false;
    }

    if (verbosity > 2) {
        cudaStreamSynchronize(stream);
        std::ostringstream out;
        out << "cusparseSolver::create_preconditioner(): " << t.stop() << " s";
        OpmLog::info(out.str());
    }
    return true;
} // end create_preconditioner()


template <unsigned int block_size>
void cusparseSolverBackend<block_size>::solve_system(WellContributions& wellContribs, BdaResult &res) {
    // actually solve
    gpu_pbicgstab(wellContribs, res);
    cudaStreamSynchronize(stream);
    cudaCheckLastError("Something went wrong during the GPU solve");
} // end solve_system()


// copy result to host memory
// caller must be sure that x is a valid array
template <unsigned int block_size>
void cusparseSolverBackend<block_size>::get_result(double *x) {
    Timer t;

    cudaMemcpyAsync(x, d_x, N * sizeof(double), cudaMemcpyDeviceToHost, stream);
    cudaStreamSynchronize(stream);

    if (verbosity > 2) {
        std::ostringstream out;
        out << "cusparseSolver::get_result(): " << t.stop() << " s";
        OpmLog::info(out.str());
    }
} // end get_result()



template <unsigned int block_size>
SolverStatus cusparseSolverBackend<block_size>::solve_system(std::shared_ptr<BlockedMatrix> matrix,
                                                              double *b,
                                                              std::shared_ptr<BlockedMatrix> jacMatrix,
                                                              WellContributions& wellContribs,
                                                              BdaResult &res)
{
    if (initialized == false) {
        initialize(matrix, jacMatrix);
        copy_system_to_gpu(matrix, b, jacMatrix);
    } else {
        update_system_on_gpu(matrix, b, jacMatrix);
    }
    if (analysis_done == false) {
        if (!analyse_matrix()) {
            return SolverStatus::BDA_SOLVER_ANALYSIS_FAILED;
        }
    }
    if (create_preconditioner()) {
        solve_system(wellContribs, res);
    } else {
        return SolverStatus::BDA_SOLVER_CREATE_PRECONDITIONER_FAILED;
    }
    return SolverStatus::BDA_SOLVER_SUCCESS;
}


#define INSTANTIATE_BDA_FUNCTIONS(n)                                                       \
template cusparseSolverBackend<n>::cusparseSolverBackend(int, int, double, unsigned int);  \

INSTANTIATE_BDA_FUNCTIONS(1);
INSTANTIATE_BDA_FUNCTIONS(2);
INSTANTIATE_BDA_FUNCTIONS(3);
INSTANTIATE_BDA_FUNCTIONS(4);
INSTANTIATE_BDA_FUNCTIONS(5);
INSTANTIATE_BDA_FUNCTIONS(6);

#undef INSTANTIATE_BDA_FUNCTIONS

} // namespace Accelerator
} // namespace Opm
