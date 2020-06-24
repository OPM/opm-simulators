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
#include <memory>
#include <sstream>

#include <opm/common/OpmLog/OpmLog.hpp>
#include <opm/common/ErrorMacros.hpp>
#include <opm/material/common/Unused.hpp>

#include <opm/simulators/linalg/bda/BdaBridge.hpp>
#include <opm/simulators/linalg/bda/BdaSolverStatus.hpp>
#include <opm/simulators/linalg/bda/BdaResult.hpp>

#define PRINT_TIMERS_BRIDGE 0

typedef Dune::InverseOperatorResult InverseOperatorResult;

namespace Opm
{

    using bda::BdaResult;
    using bda::BdaSolver;
    using bda::BdaSolverStatus;

template <class BridgeMatrix, class BridgeVector, int block_size>
BdaBridge<BridgeMatrix, BridgeVector, block_size>::BdaBridge(std::string gpu_mode, int linear_solver_verbosity, int maxit, double tolerance)
{
    if (gpu_mode.compare("cusparse") == 0) {
#if HAVE_CUDA
        use_gpu = true;
        backend.reset(new bda::cusparseSolverBackend<block_size>(linear_solver_verbosity, maxit, tolerance));
#else
        OPM_THROW(std::logic_error, "Error cusparseSolver was chosen, but CUDA was not found by CMake");
#endif
    } else if (gpu_mode.compare("opencl") == 0) {
#if HAVE_OPENCL
        use_gpu = true;
        backend.reset(new bda::openclSolverBackend<block_size>(linear_solver_verbosity, maxit, tolerance));
#else
        OPM_THROW(std::logic_error, "Error openclSolver was chosen, but OpenCL was not found by CMake");
#endif
    } else if (gpu_mode.compare("none") == 0) {
        use_gpu = false;
    } else {
        OPM_THROW(std::logic_error, "Error unknown value for parameter 'GpuMode', should be passed like '--gpu-mode=[none|cusparse|opencl]");
    }
}



template <class BridgeMatrix>
int checkZeroDiagonal(BridgeMatrix& mat) {
    static std::vector<typename BridgeMatrix::size_type> diag_indices;   // contains offsets of the diagonal nnzs
    int numZeros = 0;
    const int dim = 3;                    // might be replaced with mat[0][0].N() or BridgeMatrix::block_type::size()
    const double zero_replace = 1e-15;
    if (diag_indices.size() == 0) {
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
// sparsity pattern should stay the same due to matrix-add-well-contributions
// this could be removed if Dune::BCRSMatrix features an API call that returns colIndices and rowPointers
template <class BridgeMatrix>
void getSparsityPattern(BridgeMatrix& mat, std::vector<int> &h_rows, std::vector<int> &h_cols) {
    int sum_nnzs = 0;

    // convert colIndices and rowPointers
    if (h_rows.size() == 0) {
        h_rows.emplace_back(0);
            for (typename BridgeMatrix::const_iterator r = mat.begin(); r != mat.end(); ++r) {
                int size_row = 0;
                for (auto c = r->begin(); c != r->end(); ++c) {
                    h_cols.emplace_back(c.index());
                    size_row++;
                }
                sum_nnzs += size_row;
                h_rows.emplace_back(sum_nnzs);
            }

        // h_rows and h_cols could be changed to 'unsigned int', but cusparse expects 'int'
        if (static_cast<unsigned int>(h_rows[mat.N()]) != mat.nonzeroes()) {
            OPM_THROW(std::logic_error, "Error size of rows do not sum to number of nonzeroes in BdaBridge::getSparsityPattern()");
        }
    }
} // end getSparsityPattern()


template <class BridgeMatrix, class BridgeVector, int block_size>
void BdaBridge<BridgeMatrix, BridgeVector, block_size>::solve_system(BridgeMatrix *mat OPM_UNUSED, BridgeVector &b OPM_UNUSED, WellContributions& wellContribs OPM_UNUSED, InverseOperatorResult &res OPM_UNUSED)
{

    if (use_gpu) {
        BdaResult result;
        result.converged = false;
        static std::vector<int> h_rows;
        static std::vector<int> h_cols;
        const int dim = (*mat)[0][0].N();
        const int N = mat->N()*dim;
        const int nnz = (h_rows.empty()) ? mat->nonzeroes()*dim*dim : h_rows.back()*dim*dim;

        if (dim != 3) {
            OpmLog::warning("cusparseSolver only accepts blocksize = 3 at this time, will use Dune for the remainder of the program");
            use_gpu = false;
            return;
        }

        if (h_rows.capacity() == 0) {
            h_rows.reserve(N+1);
            h_cols.reserve(nnz);
#if PRINT_TIMERS_BRIDGE
            Dune::Timer t;
#endif
            getSparsityPattern(*mat, h_rows, h_cols);
#if PRINT_TIMERS_BRIDGE
            std::ostringstream out;
            out << "getSparsityPattern() took: " << t.stop() << " s";
            OpmLog::info(out.str());
#endif
        }

#if PRINT_TIMERS_BRIDGE
        Dune::Timer t_zeros;
        int numZeros = checkZeroDiagonal(*mat);
        std::ostringstream out;
        out << "Checking zeros took: " << t_zeros.stop() << " s, found " << numZeros << " zeros";
        OpmLog::info(out.str());
#else
        checkZeroDiagonal(*mat);
#endif


        /////////////////////////
        // actually solve

        typedef BdaSolverStatus::Status Status;
        // assume that underlying data (nonzeroes) from mat (Dune::BCRSMatrix) are contiguous, if this is not the case, cusparseSolver is expected to perform undefined behaviour
        Status status = backend->solve_system(N, nnz, dim, static_cast<double*>(&(((*mat)[0][0][0][0]))), h_rows.data(), h_cols.data(), static_cast<double*>(&(b[0][0])), wellContribs, result);
        switch(status) {
        case Status::BDA_SOLVER_SUCCESS:
            //OpmLog::info("BdaSolver converged");
            break;
        case Status::BDA_SOLVER_ANALYSIS_FAILED:
            OpmLog::warning("BdaSolver could not analyse level information of matrix, perhaps there is still a 0.0 on the diagonal of a block on the diagonal");
            break;
        case Status::BDA_SOLVER_CREATE_PRECONDITIONER_FAILED:
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
    }else{
        res.converged = false;
    }
}


template <class BridgeMatrix, class BridgeVector, int block_size>
void BdaBridge<BridgeMatrix, BridgeVector, block_size>::get_result(BridgeVector &x OPM_UNUSED) {
    if (use_gpu) {
        backend->get_result(static_cast<double*>(&(x[0][0])));
    }
}

#define INSTANTIATE_BDA_FUNCTIONS(n)                                                                                        \
template BdaBridge<Dune::BCRSMatrix<Opm::MatrixBlock<double, n, n>, std::allocator<Opm::MatrixBlock<double, n, n> > >,      \
Dune::BlockVector<Dune::FieldVector<double, n>, std::allocator<Dune::FieldVector<double, n> > >,                            \
n>::BdaBridge                                                                                                               \
(std::string gpu_mode_, int linear_solver_verbosity, int maxit, double tolerance);                                          \
                                                                                                                            \
template void BdaBridge<Dune::BCRSMatrix<Opm::MatrixBlock<double, n, n>, std::allocator<Opm::MatrixBlock<double, n, n> > >, \
Dune::BlockVector<Dune::FieldVector<double, n>, std::allocator<Dune::FieldVector<double, n> > >,                            \
n>::solve_system                                                                                                            \
(Dune::BCRSMatrix<Opm::MatrixBlock<double, n, n>, std::allocator<Opm::MatrixBlock<double, n, n> > >*,                       \
    Dune::BlockVector<Dune::FieldVector<double, n>, std::allocator<Dune::FieldVector<double, n> > >&,                       \
    WellContributions&, InverseOperatorResult&);                                                                            \
                                                                                                                            \
template void BdaBridge<Dune::BCRSMatrix<Opm::MatrixBlock<double, n, n>, std::allocator<Opm::MatrixBlock<double, n, n> > >, \
Dune::BlockVector<Dune::FieldVector<double, n>, std::allocator<Dune::FieldVector<double, n> > >,                            \
n>::get_result                                                                                                              \
(Dune::BlockVector<Dune::FieldVector<double, n>, std::allocator<Dune::FieldVector<double, n> > >&);                         \

INSTANTIATE_BDA_FUNCTIONS(1);
INSTANTIATE_BDA_FUNCTIONS(2);
INSTANTIATE_BDA_FUNCTIONS(3);
INSTANTIATE_BDA_FUNCTIONS(4);

#undef INSTANTIATE_BDA_FUNCTIONS

} // namespace Opm


