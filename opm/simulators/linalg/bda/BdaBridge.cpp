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
#include <opm/common/TimingMacros.hpp>
#include "dune/istl/bcrsmatrix.hh"
#include <opm/simulators/linalg/matrixblock.hh>

#include <opm/common/OpmLog/OpmLog.hpp>
#include <opm/common/ErrorMacros.hpp>

#include <opm/simulators/linalg/bda/BdaBridge.hpp>
#include <opm/simulators/linalg/bda/BdaResult.hpp>
#include <opm/simulators/linalg/bda/WellContributions.hpp>

#if HAVE_CUDA
#include <opm/simulators/linalg/bda/cuda/cusparseSolverBackend.hpp>
#endif

#if HAVE_OPENCL
#include <opm/simulators/linalg/bda/opencl/openclSolverBackend.hpp>
#include <opm/simulators/linalg/bda/opencl/openclWellContributions.hpp>
#endif

#if HAVE_AMGCL
#include <opm/simulators/linalg/bda/amgclSolverBackend.hpp>
#endif

#if HAVE_ROCALUTION
#include <opm/simulators/linalg/bda/rocalutionSolverBackend.hpp>
#endif

#if HAVE_ROCSPARSE
#include <opm/simulators/linalg/bda/rocsparseSolverBackend.hpp>
#endif

typedef Dune::InverseOperatorResult InverseOperatorResult;

namespace Opm
{

    using Opm::Accelerator::BdaResult;
    using Opm::Accelerator::BdaSolver;
    using Opm::Accelerator::SolverStatus;

template <class BridgeMatrix, class BridgeVector, int block_size>
BdaBridge<BridgeMatrix, BridgeVector, block_size>::BdaBridge(std::string accelerator_mode_,
                                                             int linear_solver_verbosity,
                                                             [[maybe_unused]] int maxit,
                                                             [[maybe_unused]] double tolerance,
                                                             [[maybe_unused]] unsigned int platformID,
                                                             [[maybe_unused]] unsigned int deviceID,
                                                             [[maybe_unused]] bool opencl_ilu_parallel,
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
        backend.reset(new Opm::Accelerator::openclSolverBackend<block_size>(linear_solver_verbosity, maxit, tolerance, platformID, deviceID, opencl_ilu_parallel, linsolver));
#else
        OPM_THROW(std::logic_error, "Error openclSolver was chosen, but OpenCL was not found by CMake");
#endif
    } else if (accelerator_mode.compare("amgcl") == 0) {
#if HAVE_AMGCL
        use_gpu = true; // should be replaced by a 'use_bridge' boolean
        backend.reset(new Opm::Accelerator::amgclSolverBackend<block_size>(linear_solver_verbosity, maxit, tolerance, platformID, deviceID));
#else
        OPM_THROW(std::logic_error, "Error amgclSolver was chosen, but amgcl was not found by CMake");
#endif
    } else if (accelerator_mode.compare("rocalution") == 0) {
#if HAVE_ROCALUTION
        use_gpu = true; // should be replaced by a 'use_bridge' boolean
        backend.reset(new Opm::Accelerator::rocalutionSolverBackend<block_size>(linear_solver_verbosity, maxit, tolerance));
#else
        OPM_THROW(std::logic_error, "Error rocalutionSolver was chosen, but rocalution was not found by CMake");
#endif
    } else if (accelerator_mode.compare("rocsparse") == 0) {
#if HAVE_ROCSPARSE
        use_gpu = true; // should be replaced by a 'use_bridge' boolean
        backend.reset(new Opm::Accelerator::rocsparseSolverBackend<block_size>(linear_solver_verbosity, maxit, tolerance, platformID, deviceID));
#else
        OPM_THROW(std::logic_error, "Error rocsparseSolver was chosen, but rocsparse/rocblas was not found by CMake");
#endif
    } else if (accelerator_mode.compare("none") == 0) {
        use_gpu = false;
    } else {
        OPM_THROW(std::logic_error, "Error unknown value for parameter 'AcceleratorMode', should be passed like '--accelerator-mode=[none|cusparse|opencl|amgcl|rocalution|rocsparse]");
    }
}



template <class BridgeMatrix>
int replaceZeroDiagonal(BridgeMatrix& mat, std::vector<typename BridgeMatrix::size_type>& diag_indices) {
    int numZeros = 0;
    const int dim = mat[0][0].N();                    // might be replaced with BridgeMatrix::block_type::size()
    const double zero_replace = 1e-15;
    if (diag_indices.empty()) {
        int Nb = mat.N();
        diag_indices.reserve(Nb);
        for (typename BridgeMatrix::iterator r = mat.begin(); r != mat.end(); ++r) {
            auto diag = r->find(r.index());  // diag is an iterator
            assert(diag.index() == r.index()); // every row must have a diagonal block
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
void BdaBridge<BridgeMatrix, BridgeVector, block_size>::copySparsityPatternFromISTL(const BridgeMatrix& mat, std::vector<int> &h_rows, std::vector<int> &h_cols) {

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
        OPM_THROW(std::logic_error, "Error size of rows do not sum to number of nonzeroes in BdaBridge::copySparsityPatternFromISTL()");
    }
}


// check if the nnz values of the matrix are in contiguous memory
// this is done by checking if the distance between the last value of the last block of row 0 and
// the first value of the first row of row 1 is equal to 1
// if the matrix only has 1 row, it is always contiguous
template <class BridgeMatrix>
void checkMemoryContiguous(const BridgeMatrix& mat) {
    auto block_size = mat[0][0].N();
    auto row = mat.begin();
    auto last_of_row0 = row->begin();

    // last_of_row0 points to last block, not to row->end()
    for(auto tmp = row->begin(); tmp != row->end(); ++tmp) {
        last_of_row0 = tmp;
    }

    bool isContiguous = mat.N() < 2 || std::distance(&((*last_of_row0)[block_size-1][block_size-1]), &(*mat[1].begin())[0][0]) == 1;

    if (!isContiguous) {
        OPM_THROW(std::logic_error, "Error memory of Matrix looks not contiguous");
    }
}


template <class BridgeMatrix, class BridgeVector, int block_size>
void BdaBridge<BridgeMatrix, BridgeVector, block_size>::solve_system(BridgeMatrix* bridgeMat,
                                                                     BridgeMatrix* jacMat,
                                                                     int numJacobiBlocks,
                                                                     BridgeVector& b,
                                                                     WellContributions& wellContribs,
                                                                     InverseOperatorResult& res)
{
    if (use_gpu) {
        BdaResult result;
        result.converged = false;
        const int dim = (*bridgeMat)[0][0].N();
        const int Nb = bridgeMat->N();
        const int nnzb = bridgeMat->nonzeroes();

        if (dim != 3) {
            OpmLog::warning("BdaSolver only accepts blocksize = 3 at this time, will use Dune for the remainder of the program");
            use_gpu = false;
            return;
        }

        if (!matrix) {
            h_rows.reserve(Nb+1);
            h_cols.reserve(nnzb);
            copySparsityPatternFromISTL(*bridgeMat, h_rows, h_cols);
            checkMemoryContiguous(*bridgeMat);
            matrix = std::make_unique<Opm::Accelerator::BlockedMatrix>(Nb, nnzb, block_size, static_cast<double*>(&(((*bridgeMat)[0][0][0][0]))), h_cols.data(), h_rows.data());
        }

        Dune::Timer t_zeros;
        int numZeros = replaceZeroDiagonal(*bridgeMat, diagIndices);
        if (verbosity >= 2) {
            std::ostringstream out;
            out << "Checking zeros took: " << t_zeros.stop() << " s, found " << numZeros << " zeros";
            OpmLog::info(out.str());
        }

        if (numJacobiBlocks >= 2) {
            const int jacNnzb = (h_jacRows.empty()) ? jacMat->nonzeroes() : h_jacRows.back();

            if (!jacMatrix) {
                h_jacRows.reserve(Nb+1);
                h_jacCols.reserve(jacNnzb);
                copySparsityPatternFromISTL(*jacMat, h_jacRows, h_jacCols);
                checkMemoryContiguous(*jacMat);
                jacMatrix = std::make_unique<Opm::Accelerator::BlockedMatrix>(Nb, jacNnzb, block_size, static_cast<double*>(&(((*jacMat)[0][0][0][0]))), h_jacCols.data(), h_jacRows.data());
            }

            Dune::Timer t_zeros2;
            int jacNumZeros = replaceZeroDiagonal(*jacMat, jacDiagIndices);
            if (verbosity >= 2) {
                std::ostringstream out;
                out << "Checking zeros for jacMat took: " << t_zeros2.stop() << " s, found " << jacNumZeros << " zeros";
                OpmLog::info(out.str());
            }
        }

        /////////////////////////
        // actually solve
        // assume that underlying data (nonzeroes) from b (Dune::BlockVector) are contiguous, if this is not the case, the chosen BdaSolver is expected to perform undefined behaviour
        SolverStatus status = backend->solve_system(matrix, static_cast<double*>(&(b[0][0])), jacMatrix, wellContribs, result);

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
    if (use_gpu) {
        backend->get_result(static_cast<double*>(&(x[0][0])));
    }
}

template <class BridgeMatrix, class BridgeVector, int block_size>
void BdaBridge<BridgeMatrix, BridgeVector, block_size>::initWellContributions([[maybe_unused]] WellContributions& wellContribs,
                                                                              [[maybe_unused]] unsigned N) {
    if(accelerator_mode.compare("opencl") == 0){
#if HAVE_OPENCL
        const auto openclBackend = static_cast<const Opm::Accelerator::openclSolverBackend<block_size>*>(backend.get());
        static_cast<WellContributionsOCL&>(wellContribs).setOpenCLEnv(openclBackend->context.get(), openclBackend->queue.get());
#else
        OPM_THROW(std::logic_error, "Error openclSolver was chosen, but OpenCL was not found by CMake");
#endif
    }
    wellContribs.setVectorSize(N);
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
