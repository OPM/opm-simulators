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

#ifndef BDABRIDGE_HEADER_INCLUDED
#define BDABRIDGE_HEADER_INCLUDED

#include "dune/istl/solver.hh" // for struct InverseOperatorResult

#include "dune/istl/bcrsmatrix.hh"
#include <opm/simulators/linalg/matrixblock.hh>

#include <opm/simulators/linalg/bda/WellContributions.hpp>

#if HAVE_CUDA
#include <opm/simulators/linalg/bda/cusparseSolverBackend.hpp>
#endif

#if HAVE_OPENCL
#include <opm/simulators/linalg/bda/openclSolverBackend.hpp>
#endif

namespace Opm
{

typedef Dune::InverseOperatorResult InverseOperatorResult;

/// BdaBridge acts as interface between opm-simulators with the BdaSolvers
template <class BridgeMatrix, class BridgeVector, int block_size>
class BdaBridge
{
private:
    std::unique_ptr<bda::BdaSolver<block_size> > backend;
    bool use_gpu = false;

public:
    /// Construct a BdaBridge
    /// \param[in] gpu_mode                   to select if a gpu solver is used, is passed via command-line: '--gpu-mode=[none|cusparse|opencl]'
    /// \param[in] linear_solver_verbosity    verbosity of BdaSolver
    /// \param[in] maxit                      maximum number of iterations for BdaSolver
    /// \param[in] tolerance                  required relative tolerance for BdaSolver
    /// \param[in] platformID                 the OpenCL platform ID to be used
    /// \param[in] deviceID                   the device ID to be used by the cusparse- and openclSolvers, too high values could cause runtime errors
    BdaBridge(std::string gpu_mode, int linear_solver_verbosity, int maxit, double tolerance, unsigned int platformID, unsigned int deviceID);


    /// Solve linear system, A*x = b
    /// \warning Values of A might get overwritten!
    /// \param[in] mat          matrix A, should be of type Dune::BCRSMatrix
    /// \param[in] b            vector b, should be of type Dune::BlockVector
    /// \param[in] wellContribs contains all WellContributions, to apply them separately, instead of adding them to matrix A
    /// \param[inout] result    summary of solver result
    void solve_system(BridgeMatrix *mat, BridgeVector &b, WellContributions& wellContribs, InverseOperatorResult &result);

    /// Get the resulting x vector
    /// \param[inout] x    vector x, should be of type Dune::BlockVector
    void get_result(BridgeVector &x);

    /// Return whether the BdaBridge will use the GPU or not
    /// return whether the BdaBridge will use the GPU or not
    bool getUseGpu(){
        return use_gpu;
    }

}; // end class BdaBridge

}

#endif
