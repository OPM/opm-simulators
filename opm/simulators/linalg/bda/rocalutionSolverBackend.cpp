/*
  Copyright 2022 Equinor ASA

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
#include <cmath>
#include <sstream>
#include <fmt/format.h>

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

#include <opm/simulators/linalg/bda/rocalutionSolverBackend.hpp>

#include <rocalution.hpp>
#include <base/matrix_formats_ind.hpp> // check if blocks are interpreted as row-major or column-major

#ifdef HIP_HAVE_CUDA_DEFINED
#define HAVE_CUDA HIP_HAVE_CUDA_DEFINED
#undef HIP_HAVE_CUDA_DEFINED
#endif

namespace Opm
{
namespace Accelerator
{

using Opm::OpmLog;
using Dune::Timer;

template <unsigned int block_size>
rocalutionSolverBackend<block_size>::rocalutionSolverBackend(int verbosity_, int maxit_, double tolerance_) : BdaSolver<block_size>(verbosity_, maxit_, tolerance_) {
    rocalution::init_rocalution();
    rocalution::info_rocalution();
    roc_solver = std::make_unique<rocalution::BiCGStab<rocalution::LocalMatrix<double>, rocalution::LocalVector<double>, double> >();
    roc_prec = std::make_unique<rocalution::ILU<rocalution::LocalMatrix<double>, rocalution::LocalVector<double>, double> >();
    roc_solver->Verbose(0);
    roc_solver->Init(/*abs_tol=*/1e-15, tolerance, /*divergence_tol=*/1e3, maxit);
}


template <unsigned int block_size>
rocalutionSolverBackend<block_size>::~rocalutionSolverBackend() {
    // normally, these rocalution variables are destroyed after the destructor automatically,
    // but sometimes it segfaults, both with test_rocalutionSolver and with an actual case
    // release both variables here to prevent that segfault
    roc_prec.release();
    roc_solver.release();
    rocalution::stop_rocalution();
}


template <unsigned int block_size>
void rocalutionSolverBackend<block_size>::initialize(BlockedMatrix *matrix) {
    this->Nb = matrix->Nb;
    this->N = Nb * block_size;
    this->nnzb = matrix->nnzbs;
    this->nnz = nnzb * block_size * block_size;

    std::ostringstream out;
    out << fmt::format("Initializing rocalution, matrix size: {} blockrows, nnzb: {}\n", Nb, nnzb);
    out << fmt::format("Maxit: {}, tolerance: {:e}\n", maxit, tolerance);
    OpmLog::info(out.str());

    h_x.resize(Nb * block_size);

    initialized = true;
} // end initialize()


template <unsigned int block_size>
void rocalutionSolverBackend<block_size>::convert_matrix(BlockedMatrix *matrix) {
    Timer t;

    for(int i = 0; i < Nb+1; ++i){
        tmp_rowpointers[i] = matrix->rowPointers[i];
    }
    for(int i = 0; i < nnzb; ++i){
        tmp_colindices[i] = matrix->colIndices[i];
    }

    // convert values inside block from row major to col major
    // this is the same as transposing a block
    // when compiling rocm from scratch, it is possible to choose the direction, making this transposing unnecessary
    // BCSR_IND_BASE == 0: rocalution expects column-major
    // BCSR_IND_BASE == 1: rocalution expects row-major
    if (BCSR_IND_BASE == 0) {
        for(int i = 0; i < nnzb; ++i){
            tmp_nnzvalues[i * block_size * block_size + 0] = matrix->nnzValues[i * block_size * block_size + 0];
            tmp_nnzvalues[i * block_size * block_size + 1] = matrix->nnzValues[i * block_size * block_size + 3];
            tmp_nnzvalues[i * block_size * block_size + 2] = matrix->nnzValues[i * block_size * block_size + 6];
            tmp_nnzvalues[i * block_size * block_size + 3] = matrix->nnzValues[i * block_size * block_size + 1];
            tmp_nnzvalues[i * block_size * block_size + 4] = matrix->nnzValues[i * block_size * block_size + 4];
            tmp_nnzvalues[i * block_size * block_size + 5] = matrix->nnzValues[i * block_size * block_size + 7];
            tmp_nnzvalues[i * block_size * block_size + 6] = matrix->nnzValues[i * block_size * block_size + 2];
            tmp_nnzvalues[i * block_size * block_size + 7] = matrix->nnzValues[i * block_size * block_size + 5];
            tmp_nnzvalues[i * block_size * block_size + 8] = matrix->nnzValues[i * block_size * block_size + 8];
        }
    }
    if (verbosity >= 3) {
        std::ostringstream out;
        out << "rocalutionSolver::convert_matrix(): " << t.stop() << " s";
        OpmLog::info(out.str());
    }
}


// copy result to host memory
// caller must be sure that x is a valid array
template <unsigned int block_size>
void rocalutionSolverBackend<block_size>::get_result(double *x) {
    Timer t;

    std::copy(h_x.begin(), h_x.end(), x);

    if (verbosity >= 3) {
        std::ostringstream out;
        out << "rocalutionSolver::get_result(): " << t.stop() << " s";
        OpmLog::info(out.str());
    }
} // end get_result()


template <unsigned int block_size>
SolverStatus rocalutionSolverBackend<block_size>::solve_system(std::shared_ptr<BlockedMatrix> matrix,
                                                           double *b,
                                                           [[maybe_unused]] std::shared_ptr<BlockedMatrix> jacMatrix,
                                                           [[maybe_unused]] WellContributions& wellContribs,
                                                           BdaResult &res)
{
    if (initialized == false) {
        initialize(matrix.get());
    }

    tmp_rowpointers = new int[Nb+1];
    tmp_colindices = new int[nnzb];
    tmp_nnzvalues = new double[nnzb*block_size*block_size];

    convert_matrix(matrix.get());

    rocalution::LocalVector<double> roc_x;
    rocalution::LocalVector<double> roc_rhs;
    rocalution::LocalMatrix<double> roc_mat;

    // this also transfers ownership to the allocated memory to rocalution
    // and sets the tmp_* pointers to nullptr
    roc_mat.SetDataPtrBCSR(
        &tmp_rowpointers,
        &tmp_colindices,
        &tmp_nnzvalues,
        "matrix A", nnzb, Nb, Nb, block_size);

    roc_mat.MoveToAccelerator();
    roc_x.MoveToAccelerator();
    roc_rhs.MoveToAccelerator();

    roc_x.Allocate("x", roc_mat.GetN());
    roc_rhs.Allocate("rhs", roc_mat.GetN());

    // initialize vectors
    roc_rhs.CopyFromData(b);
    roc_x.Zeros();

    roc_solver->Clear();
    roc_solver->SetOperator(roc_mat);
    roc_solver->SetPreconditioner(*roc_prec);

    // the implementation of ILU::ReBuildNumeric() does not exist at the time of writing
    // so it just calls ILU::Build() everytime
    roc_solver->ReBuildNumeric();

    double norm_0 = roc_rhs.Norm(); // since the initial guess is a vector with 0s, initial error is norm(b)

    // actually solve
    Dune::Timer t_solve;
    roc_solver->Solve(roc_rhs, &roc_x);

    // roc_solver->GetSolverStatus() returns:
    // 0, if no criteria has been reached yet
    // 1, if absolute tolerance has been reached
    // 2, if relative tolerance has been reached
    // 3, if divergence tolerance has been reached
    // 4, if maximum number of iteration has been reached

    res.elapsed = t_solve.stop();
    res.iterations = roc_solver->GetIterationCount();
    res.reduction = roc_solver->GetCurrentResidual() / norm_0;
    res.conv_rate  = static_cast<double>(pow(res.reduction, 1.0 / res.iterations));
    res.converged = (roc_solver->GetSolverStatus() == 2);


    // copy solution vector to host vector
    // if roc_x could be reused, this should be removed here
    // and roc_x should be directly copied into x in get_result()
    roc_x.MoveToHost();
    roc_x.CopyToData(h_x.data());

    if (verbosity >= 1) {
        std::ostringstream out;
        out << "=== converged: " << res.converged << ", conv_rate: " << res.conv_rate << ", time: " << res.elapsed << \
            ", time per iteration: " << res.elapsed / res.iterations << ", iterations: " << res.iterations;
        OpmLog::info(out.str());
    }

    return SolverStatus::BDA_SOLVER_SUCCESS;
}


#define INSTANTIATE_BDA_FUNCTIONS(n) \
template rocalutionSolverBackend<n>::rocalutionSolverBackend(int, int, double);

INSTANTIATE_BDA_FUNCTIONS(1);
INSTANTIATE_BDA_FUNCTIONS(2);
INSTANTIATE_BDA_FUNCTIONS(3);
INSTANTIATE_BDA_FUNCTIONS(4);
INSTANTIATE_BDA_FUNCTIONS(5);
INSTANTIATE_BDA_FUNCTIONS(6);

#undef INSTANTIATE_BDA_FUNCTIONS

} // namespace Accelerator
} // namespace Opm
