

#ifndef __NVCC__
    #error "Cannot compile for cusparse: NVIDIA compiler not found"
#endif

#include <cstdio>
#include <cstdlib>
#include <cuda_runtime.h>
#include <iostream>
#include <sys/time.h>

#include <opm/bda/cusparseSolverBackend.hpp>
#include <opm/bda/BdaResult.hpp>
#include <opm/bda/cuda_header.h>

#include "cublas_v2.h"
#include "cusparse_v2.h"

// print initial, intermediate and final norms, and used iterations
#define VERBOSE_BACKEND 0

// print more detailed timers of various solve elements and backend functions
#define PRINT_TIMERS_BACKEND 0

namespace Opm
{

    const cusparseSolvePolicy_t policy_L = CUSPARSE_SOLVE_POLICY_USE_LEVEL;
    const cusparseSolvePolicy_t policy_U = CUSPARSE_SOLVE_POLICY_USE_LEVEL;
    const cusparseSolvePolicy_t policy_M = CUSPARSE_SOLVE_POLICY_USE_LEVEL;
    //const cusparseSolvePolicy_t policy_L = CUSPARSE_SOLVE_POLICY_NO_LEVEL;
    //const cusparseSolvePolicy_t policy_U = CUSPARSE_SOLVE_POLICY_NO_LEVEL;
    //const cusparseSolvePolicy_t policy_M = CUSPARSE_SOLVE_POLICY_NO_LEVEL;
    const cusparseOperation_t trans_L  = CUSPARSE_OPERATION_NON_TRANSPOSE;
    const cusparseOperation_t trans_U  = CUSPARSE_OPERATION_NON_TRANSPOSE;
    const cusparseOperation_t trans_M  = CUSPARSE_OPERATION_NON_TRANSPOSE;
    const cusparseDirection_t dir = CUSPARSE_DIRECTION_ROW;
    //const cusparseDirection_t dir = CUSPARSE_DIRECTION_COLUMN;

    const int BLOCK_SIZE = 3;

    double second (void){
        struct timeval tv;
        gettimeofday(&tv, nullptr);
        return (double)tv.tv_sec + (double)tv.tv_usec / 1000000.0;
    }

    cusparseSolverBackend::cusparseSolverBackend(int maxit_, double tolerance_) : maxit(maxit_), tolerance(tolerance_), minit(0){
    }

    cusparseSolverBackend::~cusparseSolverBackend(){
        finalize();
    }

    // return true iff converged
    bool cusparseSolverBackend::gpu_pbicgstab(BdaResult& res){
        double t_total1, t_total2;
        int n = matrixN;
        double rho, rhop, beta, alpha, negalpha, omega, negomega, temp, temp2;
        double nrmr, nrmr0;
        rho = 1.0;
        double zero = 0.0;
        double one  = 1.0;
        double mone = -1.0;
        int i = 0;
        int j = 0;

        //compute initial residual r0=b-Ax0 (using initial guess in x)
        t_total1 = second();

        cudaSafeCall(cusparseDbsrmv(cusparseHandle, dir, trans_M, mb, mb, nnzb, &one, descr_M, d_bVals, d_bRows, d_bCols, blockDim, d_x, &zero, d_r));

        cudaSafeCall(cublasDscal(cublasHandle, n, &mone, d_r, 1));
        cudaSafeCall(cublasDaxpy(cublasHandle, n, &one, d_f, 1, d_r, 1));
        //copy residual r into r^{\hat} and p
        cudaSafeCall(cublasDcopy(cublasHandle, n, d_r, 1, d_rw, 1));
        cudaSafeCall(cublasDcopy(cublasHandle, n, d_r, 1, d_p, 1)); 
        cudaSafeCall(cublasDnrm2(cublasHandle, n, d_r, 1, &nrmr0));

#if VERBOSE_BACKEND
        printf("Initial norm: %.5e\n", nrmr0);
        printf("Tolerance: %.0e, nnzb: %d\n", tol, nnzb);
#endif

        for (i = 0; i < maxit; ++i){
            rhop = rho;
            cudaSafeCall(cublasDdot(cublasHandle, n, d_rw, 1, d_r, 1, &rho));

            if (i > 0){
                beta = (rho/rhop) * (alpha/omega);
                negomega = -omega;
                cudaSafeCall(cublasDaxpy(cublasHandle, n, &negomega, d_v, 1, d_p, 1));
                cudaSafeCall(cublasDscal(cublasHandle, n, &beta, d_p, 1));
                cudaSafeCall(cublasDaxpy(cublasHandle, n, &one, d_r, 1, d_p, 1));
            }

            //preconditioning step (lower and upper triangular solve)
            cudaSafeCall(cusparseDbsrsv2_solve(cusparseHandle,dir,\
                trans_L,mb,nnzb,&one,\
                descr_L,d_mVals,d_mRows,d_mCols,blockDim,info_L,d_p,d_t, policy_L, pBuffer));
            cudaSafeCall(cusparseDbsrsv2_solve(cusparseHandle,dir,\
                trans_U,mb,nnzb,&one,\
                descr_U,d_mVals,d_mRows,d_mCols,blockDim,info_U,d_t,d_pw, policy_U, pBuffer));

            //matrix-vector multiplication
            cudaSafeCall(cusparseDbsrmv(cusparseHandle, dir, \
                trans_M, mb, mb, nnzb, \
                &one, descr_M, d_bVals, d_bRows, d_bCols, blockDim, d_pw, &zero, d_v));

            cudaSafeCall(cublasDdot(cublasHandle, n, d_rw, 1, d_v, 1, &temp));
            alpha= rho / temp;
            negalpha = -(alpha);
            cudaSafeCall(cublasDaxpy(cublasHandle, n, &negalpha, d_v, 1, d_r, 1));
            cudaSafeCall(cublasDaxpy(cublasHandle, n, &alpha, d_pw, 1, d_x, 1));
            cudaSafeCall(cublasDnrm2(cublasHandle, n, d_r, 1, &nrmr));

            if (nrmr < tolerance*nrmr0 && i > minit){
                j = 5;
                break;
            }

            //preconditioning step (lower and upper triangular solve)
            cudaSafeCall(cusparseDbsrsv2_solve(cusparseHandle,dir,\
                trans_L,mb,nnzb,&one,\
                descr_L,d_mVals,d_mRows,d_mCols,blockDim,info_L,d_r,d_t, policy_L, pBuffer));
            cudaSafeCall(cusparseDbsrsv2_solve(cusparseHandle,dir,\
                trans_U,mb,nnzb,&one,\
                descr_U,d_mVals,d_mRows,d_mCols,blockDim,info_U,d_t,d_s, policy_U, pBuffer));

            //matrix-vector multiplication
            cudaSafeCall(cusparseDbsrmv(cusparseHandle, dir, \
                trans_M, mb, mb, nnzb, &one, descr_M, \
                d_bVals, d_bRows, d_bCols, blockDim, d_s, &zero, d_t));

            cudaSafeCall(cublasDdot(cublasHandle, n, d_t, 1, d_r, 1, &temp));
            cudaSafeCall(cublasDdot(cublasHandle, n, d_t, 1, d_t, 1, &temp2));
            omega = temp / temp2;
            negomega = -(omega);
            cudaSafeCall(cublasDaxpy(cublasHandle, n, &omega, d_s, 1, d_x, 1));
            cudaSafeCall(cublasDaxpy(cublasHandle, n, &negomega, d_t, 1, d_r, 1));

            cudaSafeCall(cublasDnrm2(cublasHandle, n, d_r, 1, &nrmr));


            if (nrmr < tolerance*nrmr0 && i > minit){
                i++;
                j=0;
                break;
            }
#if VERBOSE_BACKEND
            if(i % 1 == 0){
                printf("it: %d.0, norm: %.5e\n", i, nrmr);
            }
#endif
        }

        t_total2 = second();
#if PRINT_TIMERS_BACKEND
        printf("Total solve time: %.6f s\n", t_total2-t_total1);
#endif
#if VERBOSE_BACKEND
        printf("Iterations: %d.%d\n", i, j);
        printf("Final norm: %.5e\n", nrmr);
#endif
        res.iterations = i + (j == 5);
        res.reduction = nrmr/nrmr0;
        res.conv_rate  = static_cast<double>(pow(res.reduction,1.0/i));
        res.elapsed = t_total2-t_total1;
        res.converged = (i != maxit);
        return res.converged;
    }


    void cusparseSolverBackend::initialize(int N, int nnz, int dim){
        this->matrixN = N;
        this->nnz = nnz;
        this->blockDim = dim;
        mb = (N + dim - 1) / dim;
        printf("Initializing GPU, N: %d, nnz: %d, mb: %d\n", N, nnz, mb);
        printf("Minit: %d, maxit: %d, tol: %.1e\n", minit, maxit, tolerance);
        
        int deviceID = 0;
        cudaSafeCall(cudaSetDevice(deviceID));
        struct cudaDeviceProp props;
        cudaSafeCall(cudaGetDeviceProperties(&props, deviceID));
        std::cout << "Name: " << props.name << "\n";
        printf("CC: %d.%d\n", props.major, props.minor);

        cudaSafeCall(cudaStreamCreate(&stream));

        /* initialize cublas */
        cudaSafeCall(cublasCreate(&cublasHandle));

        /* initialize cusparse */
        cudaSafeCall(cusparseCreate(&cusparseHandle));

        /* allocate device memory for csr matrix and vectors */
        cudaSafeCall(cudaMalloc ((void**)&d_x, sizeof(d_x[0]) * N));
        cudaSafeCall(cudaMalloc ((void**)&d_f, sizeof(d_f[0]) * N));
        cudaSafeCall(cudaMalloc ((void**)&d_r, sizeof(d_r[0]) * N));
        cudaSafeCall(cudaMalloc ((void**)&d_rw,sizeof(d_rw[0])* N));
        cudaSafeCall(cudaMalloc ((void**)&d_p, sizeof(d_p[0]) * N));
        cudaSafeCall(cudaMalloc ((void**)&d_pw,sizeof(d_pw[0])* N));
        cudaSafeCall(cudaMalloc ((void**)&d_s, sizeof(d_s[0]) * N));
        cudaSafeCall(cudaMalloc ((void**)&d_t, sizeof(d_t[0]) * N));
        cudaSafeCall(cudaMalloc ((void**)&d_v, sizeof(d_v[0]) * N));
        cudaSafeCall(cudaMalloc ((void**)&d_bVals, sizeof(d_bVals[0]) * nnz));
        cudaSafeCall(cudaMalloc ((void**)&d_bCols, sizeof(d_bCols[0]) * nnz));
        cudaSafeCall(cudaMalloc ((void**)&d_bRows, sizeof(d_bRows[0]) * (mb+1)));
        cudaSafeCall(cudaMalloc ((void**)&d_mVals, sizeof(d_mVals[0]) * nnz));

        cudaSafeCall(cusparseSetStream(cusparseHandle, stream));
        cudaSafeCall(cudaMallocHost((void**)&x, sizeof(double) * N));
        cudaSafeCall(cudaMallocHost((void**)&v, sizeof(double) * N));
        cudaSafeCall(cudaMallocHost((void**)&tx, sizeof(double) * N));

        initialized = true;
    } // end initialize()

    void cusparseSolverBackend::finalize(){
        cudaSafeCall(cudaFree(d_x));
        cudaSafeCall(cudaFree(d_f));
        cudaSafeCall(cudaFree(d_r));
        cudaSafeCall(cudaFree(d_rw));
        cudaSafeCall(cudaFree(d_p));
        cudaSafeCall(cudaFree(d_pw));
        cudaSafeCall(cudaFree(d_s));
        cudaSafeCall(cudaFree(d_t));
        cudaSafeCall(cudaFree(d_v));
        cudaSafeCall(cudaFree(d_mVals));
        cudaSafeCall(cudaFree(d_bVals));
        cudaSafeCall(cudaFree(d_bCols));
        cudaSafeCall(cudaFree(d_bRows));
        cudaSafeCall(cudaFree(pBuffer));
        cudaSafeCall(cusparseDestroyBsrilu02Info(info_M));
        cudaSafeCall(cusparseDestroyBsrsv2Info(info_L));
        cudaSafeCall(cusparseDestroyBsrsv2Info(info_U));
        cudaSafeCall(cusparseDestroyMatDescr(descra));
        cudaSafeCall(cusparseDestroyMatDescr(descr_M));
        cudaSafeCall(cusparseDestroyMatDescr(descr_L));
        cudaSafeCall(cusparseDestroyMatDescr(descr_U));
        cudaSafeCall(cusparseDestroy(cusparseHandle));
        cudaSafeCall(cublasDestroy(cublasHandle));
        cudaHostUnregister(vals);
        cudaHostUnregister(cols);
        cudaHostUnregister(rows);
        cudaSafeCall(cudaStreamDestroy(stream));
        cudaSafeCall(cudaFreeHost(x));
        cudaSafeCall(cudaFreeHost(v));
        cudaSafeCall(cudaFreeHost(tx));
    } // end finalize()


    void cusparseSolverBackend::copy_system_to_gpu(double *vals, int *rows, int *cols, double *f){

        /* copy the csr matrix and vectors into device memory */
#if PRINT_TIMERS_BACKEND
        double t1, t2;
        t1 = second();
#endif

        // information: https://docs.nvidia.com/cuda/cuda-runtime-api/group__CUDART__MEMORY.html#group__CUDART__MEMORY_1ge8d5c17670f16ac4fc8fcb4181cb490c
        // possible flags for cudaHostRegister: cudaHostRegisterDefault, cudaHostRegisterPortable, cudaHostRegisterMapped, cudaHostRegisterIoMemory
        cudaSafeCall(cudaHostRegister(vals, nnz * sizeof(double), cudaHostRegisterDefault));
        cudaSafeCall(cudaHostRegister(cols, nnz * sizeof(int), cudaHostRegisterDefault));
        cudaSafeCall(cudaHostRegister(rows, (matrixN/BLOCK_SIZE+1) * sizeof(int), cudaHostRegisterDefault));
        cudaSafeCall(cudaMemcpyAsync(d_bVals, vals, (size_t)(nnz * sizeof(double)), cudaMemcpyHostToDevice, stream));
        cudaSafeCall(cudaMemcpyAsync(d_bCols, cols, (size_t)(nnz * sizeof(int)), cudaMemcpyHostToDevice, stream));
        cudaSafeCall(cudaMemcpyAsync(d_bRows, rows, (size_t)((matrixN/BLOCK_SIZE+1) * sizeof(int)), cudaMemcpyHostToDevice, stream));
        cudaSafeCall(cudaMemcpyAsync(d_f, f, (size_t)(matrixN * sizeof(d_f[0])), cudaMemcpyHostToDevice, stream));
        /* clean memory */
        cudaSafeCall(cudaMemsetAsync((void *)d_x, 0, sizeof(d_x[0]) * matrixN, stream));

        this->vals = vals;
        this->cols = cols;
        this->rows = rows;

#if PRINT_TIMERS_BACKEND
        t2 = second();
        printf("copy_system_to_gpu(): %f s\n", t2-t1);
#endif
    } // end copy_system_to_gpu()


    // don't copy rowpointers and colindices, they stay the same
    void cusparseSolverBackend::update_system_on_gpu(double *vals, double *f){

        /* copy the csr matrix and vectors into device memory */
#if PRINT_TIMERS_BACKEND
        double t1, t2;
        t1 = second();
#endif

        cudaSafeCall(cudaMemcpyAsync(d_bVals, vals, (size_t)(nnz * sizeof(vals[0])), cudaMemcpyHostToDevice, stream));
        cudaSafeCall(cudaMemcpyAsync(d_f, f, (size_t)(matrixN * sizeof(d_f[0])), cudaMemcpyHostToDevice, stream));
        /* clean memory */
        cudaSafeCall(cudaMemsetAsync((void *)d_x, 0, sizeof(d_x[0]) * matrixN, stream));

#if PRINT_TIMERS_BACKEND
        t2 = second();
        printf("update_system_on_gpu(): %f s\n", t2-t1);
#endif
    } // end update_system_on_gpu()


    void cusparseSolverBackend::reset_prec_on_gpu(){
        cudaSafeCall(cudaMemcpyAsync(d_mVals, d_bVals, (size_t)(nnz  * sizeof(d_mVals[0])), cudaMemcpyDeviceToDevice, stream));
    }


    void cusparseSolverBackend::analyse_matrix(){

        int pBufferSize_M, pBufferSize_L, pBufferSize_U, pBufferSize;

        // step 1: create a descriptor which contains
        cusparseCreateMatDescr(&descra);
        cusparseCreateMatDescr(&descr_M);
        cusparseSetMatType(descra,CUSPARSE_MATRIX_TYPE_GENERAL);
        cusparseSetMatType(descr_M, CUSPARSE_MATRIX_TYPE_GENERAL);
        cusparseIndexBase_t base_type;

        if(base){
            base_type = CUSPARSE_INDEX_BASE_ONE;
        }else{
            base_type = CUSPARSE_INDEX_BASE_ZERO;
        }

        cusparseSetMatIndexBase(descra, base_type);
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
        cudaCheckError();

        // step 2: create a empty info structure
        // we need one info for csrilu02 and two info's for csrsv2
        cusparseCreateBsrilu02Info(&info_M);
        cusparseCreateBsrsv2Info(&info_L);
        cusparseCreateBsrsv2Info(&info_U);
        cudaCheckError();

        cudaSafeCall(cudaMemcpyAsync(d_bRows, rows, sizeof(int)*(mb+1), cudaMemcpyHostToDevice, stream));
        cudaSafeCall(cudaMemcpyAsync(d_bCols, cols, sizeof(int)*nnz, cudaMemcpyHostToDevice, stream));
        cudaSafeCall(cudaMemcpyAsync(d_bVals, vals, sizeof(double)*nnz, cudaMemcpyHostToDevice, stream));

        // step 3: query how much memory used in csrilu02 and csrsv2, and allocate the buffer
        nnzb = nnz/BLOCK_SIZE/BLOCK_SIZE;
        cusparseDbsrilu02_bufferSize(cusparseHandle, dir, mb, nnzb,
            descr_M, d_bVals, d_bRows, d_bCols, blockDim, info_M, &pBufferSize_M);
        cusparseDbsrsv2_bufferSize(cusparseHandle, dir, trans_L, mb, nnzb,
            descr_L, d_bVals, d_bRows, d_bCols, blockDim, info_L, &pBufferSize_L);
        cusparseDbsrsv2_bufferSize(cusparseHandle, dir, trans_U, mb, nnzb,
            descr_U, d_bVals, d_bRows, d_bCols, blockDim, info_U, &pBufferSize_U);
        cudaCheckError();
        pBufferSize = std::max(pBufferSize_M, std::max(pBufferSize_L, pBufferSize_U));
        
        // pBuffer returned by cudaMalloc is automatically aligned to 128 bytes.
        cudaSafeCall(cudaMalloc((void**)&pBuffer, pBufferSize));

        // analysis of ilu LU decomposition
        cudaSafeCall(cusparseDbsrilu02_analysis(cusparseHandle, dir, \
            mb, nnzb, descra, d_bVals, d_bRows, d_bCols, \
            blockDim, info_M, policy_M, pBuffer));
        int structural_zero;

        status1 = cusparseXbsrilu02_zeroPivot(cusparseHandle, info_M, &structural_zero);
        if (CUSPARSE_STATUS_ZERO_PIVOT == status1){
            printf("WARNING block U(%d,%d) is not invertible\n", structural_zero, structural_zero);
            fprintf(stderr, "WARNING block U(%d,%d) is not invertible\n", structural_zero, structural_zero);
        }

        /* analyse the lower and upper triangular factors */
        cudaSafeCall(cusparseDbsrsv2_analysis(cusparseHandle, dir, trans_L, \
            mb, nnzb, descr_L, d_bVals, d_bRows, d_bCols, \
            blockDim, info_L, policy_L, pBuffer));
        cudaSafeCall(cudaDeviceSynchronize());

        cudaSafeCall(cusparseDbsrsv2_analysis(cusparseHandle, dir, trans_U, \
            mb, nnzb, descr_U, d_bVals, d_bRows, d_bCols, \
            blockDim, info_U, policy_U, pBuffer));
        cudaSafeCall(cudaDeviceSynchronize());

    } // end analyse_matrix()

    bool cusparseSolverBackend::create_preconditioner(){
        /* compute the lower and upper triangular factors using CUSPARSE csrilu0 routine (on the GPU) */
#if PRINT_TIMERS_BACKEND
        double t1, t2;
        t1 = second();
#endif
        d_mCols = d_bCols;
        d_mRows = d_bRows;
        cudaSafeCall(cusparseDbsrilu02(cusparseHandle, dir, \
            mb, nnzb, descr_M, d_mVals, d_mRows, d_mCols, \
            blockDim, info_M, policy_M, pBuffer));
        int numerical_zero;
        status1 = cusparseXbsrilu02_zeroPivot(cusparseHandle, info_M, &numerical_zero);
        if (CUSPARSE_STATUS_ZERO_PIVOT == status1){
           printf("WARNING block U(%d,%d) is not invertible\n", numerical_zero, numerical_zero);
           fprintf(stderr, "WARNING block U(%d,%d) is not invertible\n", numerical_zero, numerical_zero);
           return false;
        }
        cudaSafeCall(cudaDeviceSynchronize());
#if PRINT_TIMERS_BACKEND
        t2 = second();
        printf("Decomp time: %.6f s\n", t2-t1);
#endif
        return true;
    } // end create_preconditioner()


    // return true iff converged
    bool cusparseSolverBackend::solve_system(BdaResult &res){
        // actually solve
        bool converged = gpu_pbicgstab(res);
        cudaSafeCall(cudaDeviceSynchronize());
        return converged;
    } // end solve_system()


    double* cusparseSolverBackend::post_process(){
        /* copy the result into host memory */
#if PRINT_TIMERS_BACKEND
        double t1, t2;
        t1 = second();
#endif

        cudaSafeCall(cudaMemcpyAsync(tx, d_x, (size_t)(matrixN * sizeof(tx[0])), cudaMemcpyDeviceToHost, stream));
        cudaSafeCall(cudaStreamSynchronize(stream));

#if PRINT_TIMERS_BACKEND
        t2 = second();
        printf("Copy result back to CPU: %.6f s\n", t2-t1);
#endif

        return tx;
    } // end post_process()


    bool cusparseSolverBackend::isInitialized(){
        return initialized;
    }

}


