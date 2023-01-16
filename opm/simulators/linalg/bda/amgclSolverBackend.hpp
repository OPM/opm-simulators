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

#ifndef OPM_AMGCLSOLVER_BACKEND_HEADER_INCLUDED
#define OPM_AMGCLSOLVER_BACKEND_HEADER_INCLUDED

#include <opm/simulators/linalg/bda/BdaResult.hpp>
#include <opm/simulators/linalg/bda/BdaSolver.hpp>
#include <opm/simulators/linalg/bda/WellContributions.hpp>

#include <boost/property_tree/ptree.hpp>

#include <amgcl/amg.hpp>
#include <amgcl/backend/builtin.hpp>
#include <amgcl/adapter/crs_tuple.hpp>
#include <amgcl/make_block_solver.hpp>
#include <amgcl/relaxation/as_preconditioner.hpp>
#include <amgcl/relaxation/ilu0.hpp>
#include <amgcl/solver/bicgstab.hpp>
#include <amgcl/preconditioner/runtime.hpp>
#include <amgcl/value_type/static_matrix.hpp>

#include <memory>
#include <mutex>
#include <vector>

namespace Opm
{
namespace Accelerator
{

/// This class does not implement a solver, but converts the BCSR format to normal CSR and uses amgcl for solving
/// Note amgcl also implements blocked solvers, but looks like it needs unblocked input data
template <unsigned int block_size>
class amgclSolverBackend : public BdaSolver<block_size>
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

    typedef amgcl::static_matrix<double, block_size, block_size> dmat_type; // matrix value type in double precision
    typedef amgcl::static_matrix<double, block_size, 1> dvec_type; // the corresponding vector value type
    typedef amgcl::backend::builtin<dmat_type> CPU_Backend;

    typedef amgcl::make_solver<amgcl::runtime::preconditioner<CPU_Backend>, amgcl::runtime::solver::wrapper<CPU_Backend> > CPU_Solver;

private:

    // amgcl can use different backends, this lets the user choose
    enum Amgcl_backend_type {
        cpu,
        cuda,
        vexcl
    };

    // store matrix in CSR format
    std::vector<unsigned> A_rows, A_cols;
    std::vector<double> A_vals, rhs;
    std::vector<double> x;
    std::once_flag print_info;
    Amgcl_backend_type backend_type = cpu;

    boost::property_tree::ptree prm;         // amgcl parameters
    int iters = 0;
    double error = 0.0;

#if HAVE_CUDA
    std::once_flag cuda_initialize;
    void solve_cuda(double *b);
#endif

#if HAVE_VEXCL
    std::once_flag vexcl_initialize;
#endif
    /// Initialize host memory and determine amgcl parameters
    /// \param[in] Nb               number of blockrows
    /// \param[in] nnzbs            number of blocks
    void initialize(int Nb, int nnzbs);

    /// Convert the BCSR sparsity pattern to a CSR one
    /// \param[in] rows           array of rowPointers, contains N/dim+1 values
    /// \param[in] cols           array of columnIndices, contains nnz values
    void convert_sparsity_pattern(int *rows, int *cols);

    /// Convert the BCSR nonzero data to a CSR format
    /// \param[in] vals           array of nonzeroes, each block is stored row-wise and contiguous, contains nnz values
    /// \param[in] rows           array of rowPointers, contains N/dim+1 values
    void convert_data(double *vals, int *rows);

    /// Solve linear system
    /// \param[in] b              pointer to b vector
    /// \param[inout] res         summary of solver result
    void solve_system(double *b, BdaResult &res);

public:
    /// Construct a openclSolver
    /// \param[in] linear_solver_verbosity    verbosity of openclSolver
    /// \param[in] maxit                      maximum number of iterations for openclSolver
    /// \param[in] tolerance                  required relative tolerance for openclSolver
    /// \param[in] platformID                 the OpenCL platform to be used
    /// \param[in] deviceID                   the device to be used
    amgclSolverBackend(int linear_solver_verbosity, int maxit, double tolerance, unsigned int platformID, unsigned int deviceID);

    /// Destroy a openclSolver, and free memory
    ~amgclSolverBackend();

    /// Solve linear system, A*x = b, matrix A must be in blocked-CSR format
    /// \param[in] matrix         matrix A
    /// \param[in] b              input vector, contains N values
    /// \param[in] jacMatrix      matrix for preconditioner
    /// \param[in] wellContribs   WellContributions, to apply them separately, instead of adding them to matrix A
    /// \param[inout] res         summary of solver result
    /// \return                   status code
    SolverStatus solve_system(std::shared_ptr<BlockedMatrix> matrix, double *b,
        std::shared_ptr<BlockedMatrix> jacMatrix, WellContributions& wellContribs, BdaResult &res) override;
    
    /// Get result after linear solve, and peform postprocessing if necessary
    /// \param[inout] x          resulting x vector, caller must guarantee that x points to a valid array
    void get_result(double *x) override;

}; // end class amgclSolverBackend

} // namespace Accelerator
} // namespace Opm

#endif


