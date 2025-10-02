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

#include <opm/simulators/linalg/gpubridge/GpuResult.hpp>
#include <opm/simulators/linalg/gpubridge/GpuSolver.hpp>
#include <opm/simulators/linalg/gpubridge/WellContributions.hpp>

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
#include <type_traits>
#include <vector>

namespace Opm::Accelerator {

/// This class does not implement a solver, but converts the BCSR format to normal CSR and uses amgcl for solving
/// Note amgcl also implements blocked solvers, but looks like it needs unblocked input data
template<class Scalar, unsigned int block_size>
class amgclSolverBackend : public GpuSolver<Scalar,block_size>
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

    using dmat_type = amgcl::static_matrix<Scalar, block_size, block_size>; // matrix value type in double precision
    using dvec_type = amgcl::static_matrix<Scalar, block_size, 1>; // the corresponding vector value type
    using CPU_Backend = std::conditional_t<block_size == 1,
                                           amgcl::backend::builtin<Scalar>,
                                           amgcl::backend::builtin<dmat_type>>;

    using CPU_Solver = amgcl::make_solver<amgcl::runtime::preconditioner<CPU_Backend>,
                                          amgcl::runtime::solver::wrapper<CPU_Backend>>;

private:
    // amgcl can use different backends, this lets the user choose
    enum Amgcl_backend_type {
        cpu,
        cuda,
        vexcl
    };

    // store matrix in CSR format
    std::vector<unsigned> A_rows, A_cols;
    std::vector<Scalar> A_vals, rhs;
    std::vector<Scalar> x;
    std::once_flag print_info;
    Amgcl_backend_type backend_type = cpu;

    boost::property_tree::ptree prm;         // amgcl parameters
    int iters = 0;
    Scalar error = 0.0;

#if HAVE_CUDA
    std::once_flag cuda_initialize;
    void solve_cuda(Scalar* b);
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
    void convert_sparsity_pattern(const int *rows, const int *cols);

    /// Convert the BCSR nonzero data to a CSR format
    /// \param[in] vals           array of nonzeroes, each block is stored row-wise and contiguous, contains nnz values
    /// \param[in] rows           array of rowPointers, contains N/dim+1 values
    void convert_data(const Scalar* vals, const int* rows);

    /// Solve linear system
    /// \param[in] b              pointer to b vector
    /// \param[inout] res         summary of solver result
    void solve_system(Scalar* b, GpuResult& res);

public:
    /// Construct an amgcl solver
    /// \param[in] linear_solver_verbosity    verbosity of amgclSolver
    /// \param[in] maxit                      maximum number of iterations for amgclSolver
    /// \param[in] tolerance                  required relative tolerance for amgclSolver
    /// \param[in] platformID                 the OpenCL platform to be used
    /// \param[in] deviceID                   the device to be used
    amgclSolverBackend(int linear_solver_verbosity, int maxit,
                       Scalar tolerance, unsigned int platformID,
                       unsigned int deviceID);

    /// Destroy a openclSolver, and free memory
    ~amgclSolverBackend() override;

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

}; // end class amgclSolverBackend

} // namespace Opm::Accelerator

#endif
