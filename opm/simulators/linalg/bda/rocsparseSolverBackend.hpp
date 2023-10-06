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

#ifndef OPM_ROCSPARSESOLVER_BACKEND_HEADER_INCLUDED
#define OPM_ROCSPARSESOLVER_BACKEND_HEADER_INCLUDED

#include <memory>

#include <opm/simulators/linalg/bda/BdaResult.hpp>
#include <opm/simulators/linalg/bda/BdaSolver.hpp>
#include <opm/simulators/linalg/bda/WellContributions.hpp>

#include <rocblas/rocblas.h>
#include <rocsparse/rocsparse.h>

#include <hip/hip_version.h>

namespace Opm
{
namespace Accelerator
{

/// This class implements a rocsparse-based ilu0-bicgstab solver on GPU
template <unsigned int block_size>
class rocsparseSolverBackend : public BdaSolver<block_size>
{
    typedef BdaSolver<block_size> Base;

    using Base::N;
    using Base::Nb;
    using Base::nnz;
    using Base::nnzb;
    using Base::verbosity;
    using Base::platformID;
    using Base::deviceID;
    using Base::maxit;
    using Base::tolerance;
    using Base::initialized;

private:

    bool useJacMatrix = false;

    bool analysis_done = false;
    std::shared_ptr<BlockedMatrix> mat = nullptr;                 // original matrix
    std::shared_ptr<BlockedMatrix> jacMat = nullptr;              // matrix for preconditioner
    int nnzbs_prec = 0;    // number of nnz blocks in preconditioner matrix M

    rocsparse_direction dir = rocsparse_direction_row;
    rocsparse_operation operation = rocsparse_operation_none;
    rocsparse_handle handle;
    rocblas_handle blas_handle;
    rocsparse_mat_descr descr_A, descr_M, descr_L, descr_U;
    rocsparse_mat_info ilu_info;
#if HIP_VERSION >= 50400000
    rocsparse_mat_info spmv_info;
#endif
    hipStream_t stream;

    rocsparse_int *d_Arows, *d_Mrows;
    rocsparse_int *d_Acols, *d_Mcols;
    double *d_Avals, *d_Mvals;
    double *d_x, *d_b, *d_r, *d_rw, *d_p;     // vectors, used during linear solve
    double *d_pw, *d_s, *d_t, *d_v;
    void *d_buffer; // buffer space, used by rocsparse ilu0 analysis
    int  ver;
    char rev[64];


    /// Solve linear system using ilu0-bicgstab
    /// \param[in] wellContribs   WellContributions, to apply them separately, instead of adding them to matrix A
    /// \param[inout] res         summary of solver result
    void gpu_pbicgstab(WellContributions& wellContribs, BdaResult& res);

    /// Initialize GPU and allocate memory
    /// \param[in] matrix     matrix A
    /// \param[in] jacMatrix  matrix for preconditioner
    void initialize(std::shared_ptr<BlockedMatrix> matrix, std::shared_ptr<BlockedMatrix> jacMatrix);

    /// Copy linear system to GPU
    /// \param[in] b              input vector, contains N values
    void copy_system_to_gpu(double *b);

    /// Update linear system to GPU
    /// \param[in] b              input vector, contains N values
    void update_system_on_gpu(double *b);

    /// Analyze sparsity pattern to extract parallelism
    /// \return true iff analysis was successful
    bool analyze_matrix();

    /// Perform ilu0-decomposition
    /// \return true iff decomposition was successful
    bool create_preconditioner();

    /// Solve linear system
    /// \param[in] wellContribs   WellContributions, to apply them separately, instead of adding them to matrix A
    /// \param[inout] res         summary of solver result
    void solve_system(WellContributions &wellContribs, BdaResult &res);

public:
    /// Construct a openclSolver
    /// \param[in] linear_solver_verbosity    verbosity of openclSolver
    /// \param[in] maxit                      maximum number of iterations for openclSolver
    /// \param[in] tolerance                  required relative tolerance for openclSolver
    /// \param[in] platformID                 the OpenCL platform to be used
    /// \param[in] deviceID                   the device to be used
    rocsparseSolverBackend(int linear_solver_verbosity, int maxit, double tolerance, unsigned int platformID, unsigned int deviceID);

    /// For the CPR coarse solver
    // rocsparseSolverBackend(int linear_solver_verbosity, int maxit, double tolerance, ILUReorder opencl_ilu_reorder);

    /// Destroy a openclSolver, and free memory
    ~rocsparseSolverBackend();

    /// Solve linear system, A*x = b, matrix A must be in blocked-CSR format
    /// \param[in] matrix         matrix A
    /// \param[in] b              input vector, contains N values
    /// \param[in] jacMatrix      matrix for preconditioner
    /// \param[in] wellContribs   WellContributions, to apply them separately, instead of adding them to matrix A
    /// \param[inout] res         summary of solver result
    /// \return                   status code
    SolverStatus solve_system(std::shared_ptr<BlockedMatrix> matrix, double *b,
        std::shared_ptr<BlockedMatrix> jacMatrix, WellContributions& wellContribs, BdaResult &res) override;

    /// Solve scalar linear system, for example a coarse system of an AMG preconditioner
    /// Data is already on the GPU
    // SolverStatus solve_system(BdaResult &res);

    /// Get result after linear solve, and peform postprocessing if necessary
    /// \param[inout] x          resulting x vector, caller must guarantee that x points to a valid array
    void get_result(double *x) override;

}; // end class rocsparseSolverBackend

} // namespace Accelerator
} // namespace Opm

#endif


