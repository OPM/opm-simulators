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

#ifndef OPM_ROCALUTIONSOLVER_BACKEND_HEADER_INCLUDED
#define OPM_ROCALUTIONSOLVER_BACKEND_HEADER_INCLUDED

#include <opm/simulators/linalg/bda/BdaResult.hpp>
#include <opm/simulators/linalg/bda/BdaSolver.hpp>
#include <opm/simulators/linalg/bda/WellContributions.hpp>

namespace rocalution {
template<class Matrix, class Vector, class Scalar> class BiCGStab;
template<class Matrix, class Vector, class Scalar> class ILU;
template<class Scalar> class LocalMatrix;
template<class Scalar> class LocalVector;
}

namespace Opm
{
namespace Accelerator
{

/// This class implements a rocalution based linear solver solver on GPU
/// It uses ilu0-bicgstab
template <unsigned int block_size>
class rocalutionSolverBackend : public BdaSolver<block_size>
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
    std::vector<double> h_x; // store solution vector on host
    int *tmp_rowpointers;    // store matrix on host, this pointer is given to and freed by rocalution
    int *tmp_colindices;     // store matrix on host, this pointer is given to and freed by rocalution
    double *tmp_nnzvalues;   // store matrix on host, this pointer is given to and freed by rocalution

    std::unique_ptr<rocalution::ILU<rocalution::LocalMatrix<double>, rocalution::LocalVector<double>, double> > roc_prec;
    std::unique_ptr<rocalution::BiCGStab<rocalution::LocalMatrix<double>, rocalution::LocalVector<double>, double> > roc_solver;

    /// Initialize sizes and allocate memory
    /// \param[in] matrix     matrix A
    void initialize(BlockedMatrix *matrix);

    /// Convert matrix to rocalution format
    /// copy matrix to raw pointers, which are given to and freed by rocalution
    /// \param[in] matrix     matrix A
    void convert_matrix(BlockedMatrix *matrix);

public:

    /// Construct a rocalutionSolver
    /// also initialize rocalution library and rocalution variables
    /// \param[in] linear_solver_verbosity    verbosity of rocalutionSolver
    /// \param[in] maxit                      maximum number of iterations for rocalutionSolver
    /// \param[in] tolerance                  required relative tolerance for rocalutionSolver
    rocalutionSolverBackend(int linear_solver_verbosity, int maxit, double tolerance);

    /// Destroy a rocalutionSolver, and free memory
    ~rocalutionSolverBackend();

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
    
}; // end class rocalutionSolverBackend

} // namespace Accelerator
} // namespace Opm

#endif


