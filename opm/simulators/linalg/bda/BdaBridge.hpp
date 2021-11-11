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
#include <opm/simulators/linalg/bda/ILUReorder.hpp>

#if HAVE_FPGA
#include <opm/simulators/linalg/bda/FPGASolverBackend.hpp>
#endif

namespace Opm
{

class WellContributions;

typedef Dune::InverseOperatorResult InverseOperatorResult;
using Opm::Accelerator::ILUReorder;

/// BdaBridge acts as interface between opm-simulators with the BdaSolvers
template <class BridgeMatrix, class BridgeVector, int block_size>
class BdaBridge
{
private:
    int verbosity = 0;
    bool use_gpu = false;
    bool use_fpga = false;
    std::string accelerator_mode;
    std::unique_ptr<Opm::Accelerator::BdaSolver<block_size> > backend;

public:
    /// Construct a BdaBridge
    /// \param[in] accelerator_mode           to select if an accelerated solver is used, is passed via command-line: '--accelerator-mode=[none|cusparse|opencl|fpga]'
    /// \param[in] fpga_bitstream             FPGA programming bitstream file name, is passed via command-line: '--fpga-bitstream=[<filename>]'
    /// \param[in] linear_solver_verbosity    verbosity of BdaSolver
    /// \param[in] maxit                      maximum number of iterations for BdaSolver
    /// \param[in] tolerance                  required relative tolerance for BdaSolver
    /// \param[in] platformID                 the OpenCL platform ID to be used
    /// \param[in] deviceID                   the device ID to be used by the cusparse- and openclSolvers, too high values could cause runtime errors
    /// \param[in] opencl_ilu_reorder         select either level_scheduling or graph_coloring, see ILUReorder.hpp for explanation
    BdaBridge(std::string accelerator_mode, std::string fpga_bitstream, int linear_solver_verbosity, int maxit, double tolerance, unsigned int platformID, unsigned int deviceID, std::string opencl_ilu_reorder);


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

    /// Initialize the WellContributions object with opencl context and queue
    /// those must be set before calling BlackOilWellModel::getWellContributions() in ISTL
    /// \param[in] wellContribs   container to hold all WellContributions
    void initWellContributions(WellContributions& wellContribs);

    /// Return whether the BdaBridge will use the FPGA or not
    /// return whether the BdaBridge will use the FPGA or not
    bool getUseFpga(){
        return use_fpga;
    }

    /// Return the selected accelerator mode, this is input via the command-line
    std::string getAccleratorName(){
        return accelerator_mode;
    }

}; // end class BdaBridge

}

#endif
