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

#include <config.h>

#include "dune/istl/bcrsmatrix.hh"
#include <opm/simulators/linalg/matrixblock.hh>

#include <opm/common/OpmLog/OpmLog.hpp>
#include <opm/common/ErrorMacros.hpp>
#include <opm/material/common/Unused.hpp>

#include <opm/simulators/linalg/bda/BdaBridge.hpp>
#include <opm/simulators/linalg/bda/BdaResult.hpp>

#if HAVE_CUDA
#include <opm/simulators/linalg/bda/cusparseSolverBackend.hpp>
#endif

#if HAVE_OPENCL
#include <opm/simulators/linalg/bda/openclSolverBackend.hpp>
#include <opm/simulators/linalg/bda/openclWellContributions.hpp>
#endif

#if HAVE_FPGA
#include <opm/simulators/linalg/bda/FPGASolverBackend.hpp>
#endif

#if HAVE_AMGCL
#include <opm/simulators/linalg/bda/amgclSolverBackend.hpp>
#endif

typedef Dune::InverseOperatorResult InverseOperatorResult;

namespace Opm
{

    using Opm::Accelerator::BdaResult;
    using Opm::Accelerator::BdaSolver;
    using Opm::Accelerator::SolverStatus;
    using Opm::Accelerator::ILUReorder;

template <class BridgeMatrix, class BridgeVector, int block_size>
BdaBridge<BridgeMatrix, BridgeVector, block_size>::BdaBridge(std::string accelerator_mode_,
                                                             [[maybe_unused]] std::string fpga_bitstream,
                                                             int linear_solver_verbosity, int maxit,
                                                             double tolerance,
                                                             [[maybe_unused]] unsigned int platformID,
                                                             unsigned int deviceID,
                                                             [[maybe_unused]] std::string opencl_ilu_reorder,
                                                             [[maybe_unused]] std::string linsolver)
: verbosity(linear_solver_verbosity), accelerator_mode(accelerator_mode_)
{
    if (accelerator_mode.compare("cusparse") == 0) {
#if HAVE_CUDA
        use_gpu = true;
        backend.reset(new Opm::Accelerator::cusparseSolverBackend<block_size>(linear_solver_verbosity, maxit, tolerance, deviceID));
#else
        OPM_THROW(std::logic_error, "Error cusparseSolver was chosen, but CUDA was not found by CMake");
#endif
    } else if (accelerator_mode.compare("opencl") == 0) {
#if HAVE_OPENCL
        use_gpu = true;
        ILUReorder ilu_reorder;
        if (opencl_ilu_reorder == "") {
            ilu_reorder = Opm::Accelerator::ILUReorder::GRAPH_COLORING;  // default when not selected by user
        } else if (opencl_ilu_reorder == "level_scheduling") {
            ilu_reorder = Opm::Accelerator::ILUReorder::LEVEL_SCHEDULING;
        } else if (opencl_ilu_reorder == "graph_coloring") {
            ilu_reorder = Opm::Accelerator::ILUReorder::GRAPH_COLORING;
        } else if (opencl_ilu_reorder == "none") {
            ilu_reorder = Opm::Accelerator::ILUReorder::NONE;
        } else {
            OPM_THROW(std::logic_error, "Error invalid argument for --opencl-ilu-reorder, usage: '--opencl-ilu-reorder=[level_scheduling|graph_coloring]'");
        }
        backend.reset(new Opm::Accelerator::openclSolverBackend<block_size>(linear_solver_verbosity, maxit, tolerance, platformID, deviceID, ilu_reorder, linsolver));
#else
        OPM_THROW(std::logic_error, "Error openclSolver was chosen, but OpenCL was not found by CMake");
#endif
    } else if (accelerator_mode.compare("fpga") == 0) {
#if HAVE_FPGA
        use_fpga = true;
        ILUReorder ilu_reorder;
        if (opencl_ilu_reorder == "") {
            ilu_reorder = Opm::Accelerator::ILUReorder::LEVEL_SCHEDULING;  // default when not selected by user
        } else if (opencl_ilu_reorder == "level_scheduling") {
            ilu_reorder = Opm::Accelerator::ILUReorder::LEVEL_SCHEDULING;
        } else if (opencl_ilu_reorder == "graph_coloring") {
            ilu_reorder = Opm::Accelerator::ILUReorder::GRAPH_COLORING;
        } else {
            OPM_THROW(std::logic_error, "Error invalid argument for --opencl-ilu-reorder, usage: '--opencl-ilu-reorder=[level_scheduling|graph_coloring]'");
        }
        backend.reset(new Opm::Accelerator::FpgaSolverBackend<block_size>(fpga_bitstream, linear_solver_verbosity, maxit, tolerance, ilu_reorder));
#else
        OPM_THROW(std::logic_error, "Error fpgaSolver was chosen, but FPGA was not enabled by CMake");
#endif
    } else if (accelerator_mode.compare("amgcl") == 0) {
#if HAVE_AMGCL
        use_gpu = true; // should be replaced by a 'use_bridge' boolean
        backend.reset(new Opm::Accelerator::amgclSolverBackend<block_size>(linear_solver_verbosity, maxit, tolerance, platformID, deviceID));
#else
        OPM_THROW(std::logic_error, "Error amgclSolver was chosen, but amgcl was not found by CMake");
#endif
    } else if (accelerator_mode.compare("none") == 0) {
        use_gpu = false;
        use_fpga = false;
    } else {
        OPM_THROW(std::logic_error, "Error unknown value for parameter 'AcceleratorMode', should be passed like '--accelerator-mode=[none|cusparse|opencl|fpga|amgcl]");
    }
}



template <class BridgeMatrix>
int checkZeroDiagonal(BridgeMatrix& mat) {
    static std::vector<typename BridgeMatrix::size_type> diag_indices;   // contains offsets of the diagonal nnzs
    int numZeros = 0;
    const int dim = 3;                    // might be replaced with mat[0][0].N() or BridgeMatrix::block_type::size()
    const double zero_replace = 1e-15;
    if (diag_indices.empty()) {
        int N = mat.N();
        diag_indices.reserve(N);
        for (typename BridgeMatrix::iterator r = mat.begin(); r != mat.end(); ++r) {
            auto diag = r->find(r.index());  // diag is an iterator
            assert(diag.index() == r.index());
            for (int rr = 0; rr < dim; ++rr) {
                auto& val = (*diag)[rr][rr]; // reference to easily change the value
                if (val == 0.0) {             // could be replaced by '< 1e-30' or similar
                    val = zero_replace;
                    ++numZeros;
                }
            }
            diag_indices.emplace_back(diag.offset());
        }
    }else{
        for (typename BridgeMatrix::iterator r = mat.begin(); r != mat.end(); ++r) {
            typename BridgeMatrix::size_type offset = diag_indices[r.index()];
            auto& diag_block = r->getptr()[offset]; // diag_block is a reference to MatrixBlock, located on column r of row r
            for (int rr = 0; rr < dim; ++rr) {
                auto& val = diag_block[rr][rr];
                if (val == 0.0) {                     // could be replaced by '< 1e-30' or similar
                    val = zero_replace;
                    ++numZeros;
                }
            }
        }
    }

    return numZeros;
}


// iterate sparsity pattern from Matrix and put colIndices and rowPointers in arrays
// sparsity pattern should stay the same
// this could be removed if Dune::BCRSMatrix features an API call that returns colIndices and rowPointers
template <class BridgeMatrix, class BridgeVector, int block_size>
void BdaBridge<BridgeMatrix, BridgeVector, block_size>::getSparsityPattern(const BridgeMatrix& mat, std::vector<int> &h_rows, std::vector<int> &h_cols) {

    h_rows.clear();
    h_cols.clear();

    // convert colIndices and rowPointers
    h_rows.emplace_back(0);
    for (typename BridgeMatrix::const_iterator r = mat.begin(); r != mat.end(); ++r) {
        for (auto c = r->begin(); c != r->end(); ++c) {
            h_cols.emplace_back(c.index());
        }
        h_rows.emplace_back(h_cols.size());
    }

    // h_rows and h_cols could be changed to 'unsigned int', but cusparse expects 'int'
    if (static_cast<unsigned int>(h_rows[mat.N()]) != mat.nonzeroes()) {
        OPM_THROW(std::logic_error, "Error size of rows do not sum to number of nonzeroes in BdaBridge::getSparsityPattern()");
    }
} // end getSparsityPattern()


template <class BridgeMatrix, class BridgeVector, int block_size>
void BdaBridge<BridgeMatrix, BridgeVector, block_size>::solve_system([[maybe_unused]] BridgeMatrix* mat,
                                                                     [[maybe_unused]] BridgeVector& b,
                                                                     [[maybe_unused]] WellContributions& wellContribs,
                                                                     [[maybe_unused]] InverseOperatorResult& res)
{

    if (use_gpu || use_fpga) {
        BdaResult result;
        result.converged = false;
        static std::vector<int> h_rows;
        static std::vector<int> h_cols;
        const int dim = (*mat)[0][0].N();
        const int Nb = mat->N();
        const int N = Nb * dim;
        const int nnzb = (h_rows.empty()) ? mat->nonzeroes() : h_rows.back();
        const int nnz = nnzb * dim * dim;

        if (dim != 3) {
            OpmLog::warning("BdaSolver only accepts blocksize = 3 at this time, will use Dune for the remainder of the program");
            use_gpu = use_fpga = false;
            return;
        }

        if (h_rows.capacity() == 0) {
            h_rows.reserve(Nb+1);
            h_cols.reserve(nnzb);
            getSparsityPattern(*mat, h_rows, h_cols);
        }

        Dune::Timer t_zeros;
        int numZeros = checkZeroDiagonal(*mat);
        if (verbosity >= 2) {
            std::ostringstream out;
            out << "Checking zeros took: " << t_zeros.stop() << " s, found " << numZeros << " zeros";
            OpmLog::info(out.str());
        }


        /////////////////////////
        // actually solve

        // assume that underlying data (nonzeroes) from mat (Dune::BCRSMatrix) are contiguous, if this is not the case, the chosen BdaSolver is expected to perform undefined behaviour
        SolverStatus status = backend->solve_system(N, nnz, dim, static_cast<double*>(&(((*mat)[0][0][0][0]))), h_rows.data(), h_cols.data(), static_cast<double*>(&(b[0][0])), wellContribs, result);
        switch(status) {
        case SolverStatus::BDA_SOLVER_SUCCESS:
            //OpmLog::info("BdaSolver converged");
            break;
        case SolverStatus::BDA_SOLVER_ANALYSIS_FAILED:
            OpmLog::warning("BdaSolver could not analyse level information of matrix, perhaps there is still a 0.0 on the diagonal of a block on the diagonal");
            break;
        case SolverStatus::BDA_SOLVER_CREATE_PRECONDITIONER_FAILED:
            OpmLog::warning("BdaSolver could not create preconditioner, perhaps there is still a 0.0 on the diagonal of a block on the diagonal");
            break;
        default:
            OpmLog::warning("BdaSolver returned unknown status code");
        }

        res.iterations = result.iterations;
        res.reduction = result.reduction;
        res.converged = result.converged;
        res.conv_rate = result.conv_rate;
        res.elapsed = result.elapsed;
    } else {
        res.converged = false;
    }
}


template <class BridgeMatrix, class BridgeVector, int block_size>
void BdaBridge<BridgeMatrix, BridgeVector, block_size>::get_result([[maybe_unused]] BridgeVector& x) {
    if (use_gpu || use_fpga) {
        backend->get_result(static_cast<double*>(&(x[0][0])));
    }
}

template <class BridgeMatrix, class BridgeVector, int block_size>
void BdaBridge<BridgeMatrix, BridgeVector, block_size>::initWellContributions([[maybe_unused]] WellContributions& wellContribs) {
    if(accelerator_mode.compare("opencl") == 0){
#if HAVE_OPENCL
        const auto openclBackend = static_cast<const Opm::Accelerator::openclSolverBackend<block_size>*>(backend.get());
        static_cast<WellContributionsOCL&>(wellContribs).setOpenCLEnv(openclBackend->context.get(), openclBackend->queue.get());
#else
        OPM_THROW(std::logic_error, "Error openclSolver was chosen, but OpenCL was not found by CMake");
#endif
    }
}

// the tests use Dune::FieldMatrix, Flow uses Opm::MatrixBlock
#define INSTANTIATE_BDA_FUNCTIONS(n)                                                                                           \
template class BdaBridge<Dune::BCRSMatrix<Opm::MatrixBlock<double, n, n>, std::allocator<Opm::MatrixBlock<double, n, n> > >,   \
Dune::BlockVector<Dune::FieldVector<double, n>, std::allocator<Dune::FieldVector<double, n> > >,                               \
n>;                                                                                                                            \
                                                                                                                               \
template class BdaBridge<Dune::BCRSMatrix<Dune::FieldMatrix<double, n, n>, std::allocator<Dune::FieldMatrix<double, n, n> > >, \
Dune::BlockVector<Dune::FieldVector<double, n>, std::allocator<Dune::FieldVector<double, n> > >,                               \
n>;


INSTANTIATE_BDA_FUNCTIONS(1);
INSTANTIATE_BDA_FUNCTIONS(2);
INSTANTIATE_BDA_FUNCTIONS(3);
INSTANTIATE_BDA_FUNCTIONS(4);
INSTANTIATE_BDA_FUNCTIONS(5);
INSTANTIATE_BDA_FUNCTIONS(6);

#undef INSTANTIATE_BDA_FUNCTIONS

} // namespace Opm
