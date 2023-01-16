/*
  Copyright 2019 Equinor ASA

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

#ifndef OPM_BDASOLVER_BACKEND_HEADER_INCLUDED
#define OPM_BDASOLVER_BACKEND_HEADER_INCLUDED


#include <opm/simulators/linalg/bda/BdaResult.hpp>
#include <opm/simulators/linalg/bda/BlockedMatrix.hpp>

#include <memory>
#include <string>

namespace Opm {

class WellContributions;

namespace Accelerator {
    enum class SolverStatus {
        BDA_SOLVER_SUCCESS,
        BDA_SOLVER_ANALYSIS_FAILED,
        BDA_SOLVER_CREATE_PRECONDITIONER_FAILED,
        BDA_SOLVER_UNKNOWN_ERROR
    };

    /// This class serves to simplify choosing between different backend solvers, such as cusparseSolver and openclSolver
    /// This class is abstract, no instantiations can of it can be made, only of its children
    template <unsigned int block_size>
    class BdaSolver
    {

    protected:

        // verbosity
        // 0: print nothing during solves, only when initializing
        // 1: print number of iterations and final norm
        // 2: also print norm each iteration
        // 3: also print timings of different backend functions

        int verbosity = 0;

        int maxit = 200;
        double tolerance = 1e-2;

        int N;           // number of rows
        int Nb;          // number of blocked rows (Nb*block_size == N)
        int nnz;         // number of nonzeroes (scalars)
        int nnzb;        // number of nonzero blocks (nnzb*block_size*block_size == nnz)

        unsigned int platformID = 0; // ID of OpenCL platform to be used, only used by openclSolver now
        unsigned int deviceID = 0;   // ID of the device to be used

        bool initialized = false;

    public:
        /// Construct a BdaSolver
        /// \param[in] linear_solver_verbosity    verbosity of solver
        /// \param[in] maxit                      maximum number of iterations for solver
        /// \param[in] tolerance                  required relative tolerance for solver
        /// \param[in] platformID                 the OpenCL platform to be used, only used in openclSolver
        /// \param[in] deviceID                   the device to be used
        BdaSolver(int linear_solver_verbosity, int max_it, double tolerance_) : verbosity(linear_solver_verbosity), maxit(max_it), tolerance(tolerance_) {};
        BdaSolver(int linear_solver_verbosity, int max_it, double tolerance_, unsigned int deviceID_) : verbosity(linear_solver_verbosity), maxit(max_it), tolerance(tolerance_), deviceID(deviceID_) {};
        BdaSolver(int linear_solver_verbosity, int max_it, double tolerance_, unsigned int platformID_, unsigned int deviceID_) : verbosity(linear_solver_verbosity), maxit(max_it), tolerance(tolerance_), platformID(platformID_), deviceID(deviceID_) {};

        /// Define virtual destructor, so that the derivedclass destructor will be called
        virtual ~BdaSolver() {};

        /// Define as pure virtual functions, so derivedclass must implement them
        virtual SolverStatus solve_system(std::shared_ptr<BlockedMatrix> matrix, double *b,
            std::shared_ptr<BlockedMatrix> jacMatrix, WellContributions& wellContribs, BdaResult &res) = 0;

        virtual void get_result(double *x) = 0;

    }; // end class BdaSolver

} // namespace Accelerator
} // namespace Opm

#endif
