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

#include <opm/simulators/linalg/gpubridge/GpuResult.hpp>
#include <opm/simulators/linalg/gpubridge/GpuSolver.hpp>
#include <opm/simulators/linalg/gpubridge/WellContributions.hpp>

#include <opm/simulators/linalg/gpubridge/rocm/rocsparsePreconditioner.hpp>

#include <rocblas/rocblas.h>
#include <rocsparse/rocsparse.h>

#include <hip/hip_version.h>

namespace Opm::Accelerator {

/// This class implements a rocsparse-based ilu0-bicgstab solver on GPU
template<class Scalar, unsigned int block_size>
class rocsparseSolverBackend : public GpuSolver<Scalar,block_size>
{
    using Base = GpuSolver<Scalar,block_size>;

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
    double c_copy = 0.0; // cummulative timer measuring the total time it takes to transfer the data to the GPU

    bool useJacMatrix = false;

    bool analysis_done = false;
    std::shared_ptr<BlockedMatrix<Scalar>> mat{};                 // original matrix
    std::shared_ptr<BlockedMatrix<Scalar>> jacMat{};                 // jacobi matrix

    rocsparse_direction dir = rocsparse_direction_row;
    rocsparse_operation operation = rocsparse_operation_none;
    rocsparse_handle handle;
    rocblas_handle blas_handle;
    rocsparse_mat_descr descr_A;
#if HIP_VERSION >= 50400000
    rocsparse_mat_info spmv_info;
#endif
    hipStream_t stream;

    rocsparse_int *d_Arows, *d_Acols;
    Scalar *d_Avals;
    Scalar *d_x, *d_b, *d_r, *d_rw, *d_p; // vectors, used during linear solve
    Scalar *d_pw, *d_s, *d_t, *d_v;
    int  ver;
    char rev[64];

    std::unique_ptr<rocsparsePreconditioner<Scalar, block_size> > prec; // can perform blocked ILU0 and AMG on pressure component

    /// Solve linear system using ilu0-bicgstab
    /// \param[in] wellContribs   WellContributions, to apply them separately, instead of adding them to matrix A
    /// \param[inout] res         summary of solver result
    void gpu_pbicgstab(WellContributions<Scalar>& wellContribs, GpuResult& res);

    /// Initialize GPU and allocate memory
    /// \param[in] matrix     matrix A
    /// \param[in] jacMatrix  matrix for preconditioner
    void initialize(std::shared_ptr<BlockedMatrix<Scalar>> matrix,
                    std::shared_ptr<BlockedMatrix<Scalar>> jacMatrix);

    /// Copy linear system to GPU
    /// \param[in] b              input vector, contains N values
    void copy_system_to_gpu(Scalar* b);

    /// Update linear system to GPU
    /// \param[in] b              input vector, contains N values
    void update_system_on_gpu(Scalar* vals, Scalar* b);

    /// Analyze sparsity pattern to extract parallelism
    /// \return true iff analysis was successful
    bool analyze_matrix();

    /// Perform ilu0-decomposition
    /// \return true iff decomposition was successful
    bool create_preconditioner();

    /// Solve linear system
    /// \param[in] wellContribs   WellContributions, to apply them separately, instead of adding them to matrix A
    /// \param[inout] res         summary of solver result
    void solve_system(WellContributions<Scalar>& wellContribs, GpuResult& res);

public:
    /// Construct a rocsparseSolver
    /// \param[in] linear_solver_verbosity    verbosity of rocsparseSolver
    /// \param[in] maxit                      maximum number of iterations for rocsparseSolver
    /// \param[in] tolerance                  required relative tolerance for rocsparseSolver
    /// \param[in] platformID                 the OpenCL platform to be used
    /// \param[in] deviceID                   the device to be used
    /// \param[in] linsolver                  indicating the preconditioner, equal to the --linear-solver cmdline argument
    rocsparseSolverBackend(int linear_solver_verbosity,
                           int maxit,
                           Scalar tolerance,
                           unsigned int platformID,
                           unsigned int deviceID,
                           std::string linsolver);

    /// For the CPR coarse solver
    rocsparseSolverBackend(int linear_solver_verbosity, int maxit, Scalar tolerance, bool opencl_ilu_reorder);

    /// Destroy a openclSolver, and free memory
    ~rocsparseSolverBackend();

    /// Solve linear system, A*x = b, matrix A must be in blocked-CSR format
    /// \param[in] matrix         matrix A
    /// \param[in] b              input vector, contains N values
    /// \param[in] jacMatrix      matrix for preconditioner
    /// \param[in] wellContribs   WellContributions, to apply them separately, instead of adding them to matrix A
    /// \param[inout] res         summary of solver result
    /// \return                   status code
    SolverStatus solve_system(std::shared_ptr<BlockedMatrix<Scalar>> matrix,
                              Scalar* b,
                              std::shared_ptr<BlockedMatrix<Scalar>> jacMatrix,
                              WellContributions<Scalar>& wellContribs,
                              GpuResult& res) override;

    /// Get result after linear solve, and peform postprocessing if necessary
    /// \param[inout] x          resulting x vector, caller must guarantee that x points to a valid array
    void get_result(Scalar* x) override;

}; // end class rocsparseSolverBackend

} // namespace Opm::Accelerator

#endif
