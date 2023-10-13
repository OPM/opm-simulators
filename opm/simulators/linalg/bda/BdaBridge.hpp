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

#include <opm/simulators/linalg/bda/BdaSolver.hpp>

namespace Opm
{

class WellContributions;

typedef Dune::InverseOperatorResult InverseOperatorResult;

/// BdaBridge acts as interface between opm-simulators with the BdaSolvers
template <class BridgeMatrix, class BridgeVector, int block_size>
class BdaBridge
{
private:
    int verbosity = 0;
    bool use_gpu = false;
    std::string accelerator_mode;
    std::unique_ptr<Opm::Accelerator::BdaSolver<block_size> > backend;
    std::shared_ptr<Opm::Accelerator::BlockedMatrix> matrix;  // 'stores' matrix, actually points to h_rows, h_cols and the received BridgeMatrix for the nonzeroes
    std::shared_ptr<Opm::Accelerator::BlockedMatrix> jacMatrix;  // 'stores' preconditioner matrix, actually points to h_rows, h_cols and the received BridgeMatrix for the nonzeroes
    std::vector<int> h_rows, h_cols;  // store the sparsity pattern of the matrix
    std::vector<int> h_jacRows, h_jacCols;  // store the sparsity pattern of the jacMatrix
    std::vector<typename BridgeMatrix::size_type> diagIndices;   // contains offsets of the diagonal blocks wrt start of the row, used for replaceZeroDiagonal()
    std::vector<typename BridgeMatrix::size_type> jacDiagIndices;   // same but for jacMatrix

public:
    /// Construct a BdaBridge
    /// \param[in] accelerator_mode           to select if an accelerated solver is used, is passed via command-line: '--accelerator-mode=[none|cusparse|opencl|amgcl|rocalution|rocsparse]'
    /// \param[in] linear_solver_verbosity    verbosity of BdaSolver
    /// \param[in] maxit                      maximum number of iterations for BdaSolver
    /// \param[in] tolerance                  required relative tolerance for BdaSolver
    /// \param[in] platformID                 the OpenCL platform ID to be used
    /// \param[in] deviceID                   the device ID to be used by the cusparse- and openclSolvers, too high values could cause runtime errors
    /// \param[in] opencl_ilu_parallel        whether to parallelize the ILU decomposition and application in OpenCL with level_scheduling
    /// \param[in] linsolver                  indicating the preconditioner, equal to the --linear-solver cmdline argument
    BdaBridge(std::string accelerator_mode, int linear_solver_verbosity, int maxit, double tolerance,
        unsigned int platformID, unsigned int deviceID, bool opencl_ilu_parallel, std::string linsolver);


    /// Solve linear system, A*x = b
    /// \warning Values of A might get overwritten!
    /// \param[in] bridgeMat       matrix A, should be of type Dune::BCRSMatrix
    /// \param[in] jacMat          matrix A, but modified for the preconditioner, should be of type Dune::BCRSMatrix
    /// \param[in] numJacobiBlocks number of jacobi blocks in jacMat
    /// \param[in] b               vector b, should be of type Dune::BlockVector
    /// \param[in] wellContribs    contains all WellContributions, to apply them separately, instead of adding them to matrix A
    /// \param[inout] result       summary of solver result
    void solve_system(BridgeMatrix *bridgeMat, BridgeMatrix *jacMat, int numJacobiBlocks, BridgeVector &b, WellContributions& wellContribs, InverseOperatorResult &result);

    /// Get the resulting x vector
    /// \param[inout] x    vector x, should be of type Dune::BlockVector
    void get_result(BridgeVector &x);

    /// Return whether the BdaBridge will use the GPU or not
    /// return whether the BdaBridge will use the GPU or not
    bool getUseGpu(){
        return use_gpu;
    }

    /// Store sparsity pattern into vectors
    /// \param[in] mat       input matrix, probably BCRSMatrix
    /// \param[out] h_rows   rowpointers
    /// \param[out] h_cols   columnindices
    static void copySparsityPatternFromISTL(const BridgeMatrix& mat, std::vector<int>& h_rows, std::vector<int>& h_cols);

    /// Initialize the WellContributions object with opencl context and queue
    /// those must be set before calling BlackOilWellModel::getWellContributions() in ISTL
    /// \param[in] wellContribs   container to hold all WellContributions
    /// \param[in] N              number of rows in scalar vector that wellContribs will be applied on
    void initWellContributions(WellContributions& wellContribs, unsigned N);

    /// Return the selected accelerator mode, this is input via the command-line
    std::string getAccleratorName(){
        return accelerator_mode;
    }

}; // end class BdaBridge

}

#endif
