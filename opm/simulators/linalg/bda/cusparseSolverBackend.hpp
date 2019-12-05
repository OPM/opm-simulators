/*
  Copyright 2019 Big Data Accelerate

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

#ifndef OPM_cuSPARSESOLVER_BACKEND_HEADER_INCLUDED
#define OPM_cuSPARSESOLVER_BACKEND_HEADER_INCLUDED


#include "cublas_v2.h"
#include "cusparse_v2.h"

#include "opm/simulators/linalg/bda/BdaResult.hpp"

namespace Opm
{

class cusparseSolverBackend{

private:

    int minit;
    int maxit;
    double tolerance;

    cublasHandle_t cublasHandle;
    cusparseHandle_t cusparseHandle;
    cudaStream_t stream;
    cusparseMatDescr_t descr_B, descr_M, descr_L, descr_U;
    bsrilu02Info_t info_M;
    bsrsv2Info_t info_L, info_U;
    // b: bsr matrix, m: preconditioner
    double *d_bVals, *d_mVals;
    int *d_bCols, *d_mCols;
    int *d_bRows, *d_mRows;
    double *d_x, *d_b, *d_r, *d_rw, *d_p;
    double *d_pw, *d_s, *d_t, *d_v;
    double *vals;
    int *cols, *rows;
    double *x, *b;
    void *d_buffer;
    int N, Nb, nnz, nnzb;

    int BLOCK_SIZE;

    bool initialized = false;

    // verbosity
    // 0: print nothing during solves, only when initializing
    // 1: print number of iterations and final norm
    // 2: also print norm each iteration
    // 3: also print timings of different backend functions

    int verbosity = 0;

public:

    cusparseSolverBackend(int linear_solver_verbosity, int maxit, double tolerance);

    ~cusparseSolverBackend();

    // return true iff converged
    bool gpu_pbicgstab(BdaResult& res);

    void initialize(int N, int nnz, int dim);

    void finalize();

    void copy_system_to_gpu(double *vals, int *rows, int *cols, double *b);

    // don't copy rowpointers and colindices, they stay the same
    void update_system_on_gpu(double *vals, double *b);

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

