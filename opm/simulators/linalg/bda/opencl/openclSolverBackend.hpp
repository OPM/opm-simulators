/*
  Copyright 2020 Equinor ASA

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

#ifndef OPM_OPENCLSOLVER_BACKEND_HEADER_INCLUDED
#define OPM_OPENCLSOLVER_BACKEND_HEADER_INCLUDED

#include <opm/simulators/linalg/bda/opencl/opencl.hpp>
#include <opm/simulators/linalg/bda/BdaResult.hpp>
#include <opm/simulators/linalg/bda/BdaSolver.hpp>
#include <opm/simulators/linalg/bda/WellContributions.hpp>

#include <opm/simulators/linalg/bda/opencl/Preconditioner.hpp>

namespace Opm
{
namespace Accelerator
{

/// This class implements a opencl-based ilu0-bicgstab solver on GPU
template <unsigned int block_size>
class openclSolverBackend : public BdaSolver<block_size>
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
    double *h_b = nullptr;                // b vector, on host
    std::vector<double> vals_contiguous;  // only used if COPY_ROW_BY_ROW is true in openclSolverBackend.cpp

    // OpenCL variables must be reusable, they are initialized in initialize()
    cl::Buffer d_Avals, d_Acols, d_Arows;        // matrix in BSR format on GPU
    cl::Buffer d_x, d_b, d_rb, d_r, d_rw, d_p;   // vectors, used during linear solve
    cl::Buffer d_pw, d_s, d_t, d_v;              // vectors, used during linear solve
    cl::Buffer d_tmp;                            // used as tmp GPU buffer for dot() and norm()

    std::vector<cl::Device> devices;

    bool useJacMatrix = false;

    std::unique_ptr<Preconditioner<block_size> > prec;
                                                                  // can perform blocked ILU0 and AMG on pressure component
    bool is_root;                                                 // allow for nested solvers, the root solver is called by BdaBridge
    bool analysis_done = false;
    std::shared_ptr<BlockedMatrix> mat = nullptr;                 // original matrix
    std::shared_ptr<BlockedMatrix> jacMat = nullptr;              // matrix for preconditioner
    bool opencl_ilu_parallel;                                     // parallelize ILU operations (with level_scheduling)
    std::vector<cl::Event> events;
    cl_int err;

    /// Solve linear system using ilu0-bicgstab
    /// \param[in] wellContribs   WellContributions, to apply them separately, instead of adding them to matrix A
    /// \param[inout] res         summary of solver result
    void gpu_pbicgstab(WellContributions& wellContribs, BdaResult& res);

    /// Initialize GPU and allocate memory
    /// \param[in] matrix     matrix A
    /// \param[in] jacMatrix  matrix for preconditioner
    void initialize(std::shared_ptr<BlockedMatrix> matrix, std::shared_ptr<BlockedMatrix> jacMatrix);

    /// Copy linear system to GPU
    void copy_system_to_gpu();

    /// Reassign pointers, in case the addresses of the Dune variables have changed
    /// \param[in] vals           array of nonzeroes, each block is stored row-wise and contiguous, contains nnz values
    /// \param[in] b              input vector b, contains N values
    void update_system(double *vals, double *b);

    /// Update linear system on GPU, don't copy rowpointers and colindices, they stay the same
    void update_system_on_gpu();

    /// Analyze sparsity pattern to extract parallelism
    /// \return true iff analysis was successful
    bool analyze_matrix();

    /// Create the preconditioner, only done once per linear solve
    /// \return true iff decomposition was successful
    bool create_preconditioner();

    /// Solve linear system
    /// \param[in] wellContribs   WellContributions, to apply them separately, instead of adding them to matrix A
    ///                           could be empty
    /// \param[inout] res         summary of solver result
    void solve_system(WellContributions &wellContribs, BdaResult &res);

public:
    std::shared_ptr<cl::Context> context;
    std::shared_ptr<cl::CommandQueue> queue;

    /// Construct a openclSolver
    /// \param[in] linear_solver_verbosity    verbosity of openclSolver
    /// \param[in] maxit                      maximum number of iterations for openclSolver
    /// \param[in] tolerance                  required relative tolerance for openclSolver
    /// \param[in] platformID                 the OpenCL platform to be used
    /// \param[in] deviceID                   the device to be used
    /// \param[in] opencl_ilu_parallel        whether to parallelize the ILU decomposition and application in OpenCL with level_scheduling
    /// \param[in] linsolver                  indicating the preconditioner, equal to the --linear-solver cmdline argument
    ///                                       only ilu0, cpr_quasiimpes and isai are supported
    openclSolverBackend(int linear_solver_verbosity, int maxit, double tolerance, unsigned int platformID, unsigned int deviceID,
        bool opencl_ilu_parallel, std::string linsolver);

    /// For the CPR coarse solver
    openclSolverBackend(int linear_solver_verbosity, int maxit, double tolerance, bool opencl_ilu_parallel);

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

    /// Set OpenCL objects
    /// This class either creates them based on platformID and deviceID or receives them through this function
    /// \param[in] context   the opencl context to be used
    /// \param[in] queue     the opencl queue to be used
    void setOpencl(std::shared_ptr<cl::Context>& context, std::shared_ptr<cl::CommandQueue>& queue);
    
}; // end class openclSolverBackend

} // namespace Accelerator
} // namespace Opm

#endif


