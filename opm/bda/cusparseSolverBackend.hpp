/*
 * Copyright 1993-2015 NVIDIA Corporation.  All rights reserved.
 *
 * Please refer to the NVIDIA end user license agreement (EULA) associated
 * with this source code for terms and conditions that govern your use of
 * this software. Any use, reproduction, disclosure, or distribution of
 * this software and related documentation outside the terms of the EULA
 * is strictly prohibited.
 *
 */

#ifndef OPM_cuSPARSESOLVER_BACKEND_HEADER_INCLUDED
#define OPM_cuSPARSESOLVER_BACKEND_HEADER_INCLUDED


#include "cublas_v2.h"
#include "cusparse_v2.h"

#include "opm/bda/BdaResult.hpp"

namespace Opm
{

class cusparseSolverBackend{

private:

    int minit;
    int maxit;
    double tolerance;

    cublasHandle_t cublasHandle = 0;
    cusparseHandle_t cusparseHandle = 0;
    cudaStream_t stream = 0;
    cusparseMatDescr_t descra, descr_M, descr_L, descr_U;
    bsrilu02Info_t info_M = 0;
    bsrsv2Info_t info_L, info_U;
    cusparseStatus_t status1;
    // a: csr matrix, b: bsr matrix, m: preconditioner
    double *d_bVals, *d_mVals;
    int    *d_bCols, *d_bRows;
    int    *d_mCols, *d_mRows;
    double *d_x, *d_f, *d_r, *d_rw, *d_p;
    double *d_pw, *d_s, *d_t, *d_v;
    double *vals;
    int    *cols, *rows;
    double *x, *tx, *f, *v, *xg, *bg;
    void *pBuffer;
    int matrixN, nnz, mb, nnzb;
    int *h_nnzTotal = &nnzb;
    int base = 0;
    int blockDim;

    bool initialized = false;

public:

    cusparseSolverBackend(int maxit, double tolerance);

    ~cusparseSolverBackend();

    // return true iff converged
    bool gpu_pbicgstab(BdaResult& res);

    void initialize(int N, int nnz, int dim);

    void finalize();

    void copy_system_to_gpu(double *vals, int *rows, int *cols, double *f);

    // don't copy rowpointers and colindices, they stay the same
    void update_system_on_gpu(double *vals, double *f);

    void reset_prec_on_gpu();

    void analyse_matrix();

    bool create_preconditioner();

    // return true iff converged
    bool solve_system(BdaResult &res);

    double* post_process();

    bool isInitialized();

}; // end class cusparseSolverBackend

}

#endif

