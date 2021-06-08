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

#include <config.h>
#include <sstream>

#include <opm/common/OpmLog/OpmLog.hpp>
#include <opm/common/ErrorMacros.hpp>
#include <dune/common/timer.hh>

#include <opm/simulators/linalg/bda/BdaResult.hpp>
#include <opm/simulators/linalg/bda/amgclSolverBackend.hpp>

namespace bda
{

using Opm::OpmLog;
using Dune::Timer;

template <unsigned int block_size>
amgclSolverBackend<block_size>::amgclSolverBackend(int verbosity_, int maxit_, double tolerance_, unsigned int platformID_, unsigned int deviceID_) : BdaSolver<block_size>(verbosity_, maxit_, tolerance_, platformID_, deviceID_) {}


template <unsigned int block_size>
amgclSolverBackend<block_size>::~amgclSolverBackend() {}


template <unsigned int block_size>
void amgclSolverBackend<block_size>::initialize(int N_, int nnz_, int dim, double *vals, int *rows, int *cols) {
    this->N = N_;
    this->nnz = nnz_;
    this->nnzb = nnz_ / block_size / block_size;

    Nb = (N + dim - 1) / dim;
    std::ostringstream out;
    out << "Initializing amgclSolverBackend, matrix size: " << N << " blocks, nnzb: " << nnzb << "\n";
    out << "Maxit: " << maxit << std::scientific << ", tolerance: " << tolerance << "\n";
    out << "PlatformID: " << platformID << ", deviceID: " << deviceID << "\n";
    OpmLog::info(out.str());
    out.str("");
    out.clear();

    A_vals.resize(nnz);
    A_cols.resize(nnz);
    A_rows.resize(N + 1);
    rhs.resize(N);
    x.resize(N);

#if AMGCL_CUDA
    cusparseCreate(&bprm.cusparse_handle);
#endif

    initialized = true;
} // end initialize()



template <unsigned int block_size>
void amgclSolverBackend<block_size>::convert_sparsity_pattern(int *rows, int *cols) {
    Timer t;
    const unsigned int bs = block_size;
    int idx = 0; // indicates the unblocked write index

    A_rows[0] = 0;
    for (int row = 0; row < Nb; ++row) {
        int rowStart = rows[row];
        int rowEnd = rows[row+1];
        for (unsigned r = 0; r < bs; ++r) {
            for (int ij = rowStart; ij < rowEnd; ++ij) {
                for (unsigned c = 0; c < bs; ++c) {
                    A_cols[idx] = cols[ij] * bs + c;
                    idx++;
                }
            }
            A_rows[row*bs + r + 1] = idx;
        }
    }

    if (verbosity >= 3) {
        std::ostringstream out;
        out << "amgclSolverBackend::convert_sparsity_pattern(): " << t.stop() << " s";
        OpmLog::info(out.str());
    }
} // end convert_sparsity_pattern()

template <unsigned int block_size>
void amgclSolverBackend<block_size>::convert_data(double *vals, int *rows) {
    Timer t;
    const unsigned int bs = block_size;
    int idx = 0; // indicates the unblocked write index

    for (int row = 0; row < Nb; ++row) {
        int rowStart = rows[row];
        int rowEnd = rows[row+1];
        for (unsigned r = 0; r < bs; ++r) {
            for (int ij = rowStart; ij < rowEnd; ++ij) {
                for (unsigned c = 0; c < bs; ++c) {
                    A_vals[idx] = vals[ij*bs*bs + r*bs + c];
                    idx++;
                }
            }
        }
    }

    if (verbosity >= 3) {
        std::ostringstream out;
        out << "amgclSolverBackend::convert_data(): " << t.stop() << " s";
        OpmLog::info(out.str());
    }
} // end convert_data()


template <unsigned int block_size>
void amgclSolverBackend<block_size>::solve_system(double *b, WellContributions &wellContribs, BdaResult &res) {
    Timer t;
    int iters = 0;
    double error = 0.0;

    try {
        // create matrix object
        auto A = std::tie(N, A_rows, A_cols, A_vals);

        // set parameters
        typename Solver::params prm;
        prm.precond.damping = 0.9;
        prm.solver.maxiter = maxit;
        prm.solver.tol = tolerance;
        prm.solver.verbose = (verbosity >= 2);

        // create solver
        Solver solve(A, prm, bprm);

        // print info (once)
        std::call_once(print_info, [&](){
#if AMGCL_CUDA
            cudaDeviceProp prop;
            cudaGetDeviceProperties(&prop, deviceID);
            std::cout << prop.name << std::endl;
#endif
            // print solver structure
            std::cout << solve << std::endl;
        });

#if AMGCL_CUDA
        thrust::device_vector<double> B(b, b + N);
        thrust::device_vector<double> X(N, 0.0);

        // actually solve
        std::tie(iters, error) = solve(B, X);

        thrust::copy(X.begin(), X.end(), x.begin());
#else
        // reset x vector
        std::fill(x.begin(), x.end(), 0.0);

        // create blocked vectors
        auto b_ptr = reinterpret_cast<dvec_type*>(b);
        auto x_ptr = reinterpret_cast<dvec_type*>(x.data());
        auto B = amgcl::make_iterator_range(b_ptr, b_ptr + N / block_size);
        auto X = amgcl::make_iterator_range(x_ptr, x_ptr + N / block_size);

        // actually solve
        std::tie(iters, error) = solve(A, B, X);
#endif

    } catch (const std::exception& ex) {
        std::cerr << "Caught exception: " << ex.what() << std::endl;
        throw ex;
    }

    double time_elapsed = t.stop();

    res.iterations = iters;
    res.reduction = 0.0;
    res.elapsed = time_elapsed;
    res.converged = (iters != maxit);

    if (verbosity >= 1) {
        std::ostringstream out;
        out << "=== converged: " << res.converged << ", time: " << res.elapsed << \
            ", time per iteration: " << res.elapsed / iters << ", iterations: " << iters;
        OpmLog::info(out.str());
    }
    if (verbosity >= 3) {
        std::ostringstream out;
        out << "amgclSolverBackend::solve_system(): " << time_elapsed << " s";
        OpmLog::info(out.str());
    }

} // end solve_system()


// copy result to host memory
// caller must be sure that x is a valid array
template <unsigned int block_size>
void amgclSolverBackend<block_size>::get_result(double *x_) {
    Timer t;

    std::copy(x.begin(), x.end(), x_);

    if (verbosity >= 3) {
        std::ostringstream out;
        out << "amgclSolverBackend::get_result(): " << t.stop() << " s";
        OpmLog::info(out.str());
    }
} // end get_result()


template <unsigned int block_size>
SolverStatus amgclSolverBackend<block_size>::solve_system(int N_, int nnz_, int dim, double *vals, int *rows, int *cols, double *b, WellContributions& wellContribs, BdaResult &res) {
    if (initialized == false) {
        initialize(N_, nnz_, dim, vals, rows, cols);
        convert_sparsity_pattern(rows, cols);
    }
    convert_data(vals, rows);
    solve_system(b, wellContribs, res);
    return SolverStatus::BDA_SOLVER_SUCCESS;
}


#define INSTANTIATE_BDA_FUNCTIONS(n)                                                                \
template amgclSolverBackend<n>::amgclSolverBackend(int, int, double, unsigned int, unsigned int);   \

INSTANTIATE_BDA_FUNCTIONS(1);
INSTANTIATE_BDA_FUNCTIONS(2);
INSTANTIATE_BDA_FUNCTIONS(3);
INSTANTIATE_BDA_FUNCTIONS(4);

#undef INSTANTIATE_BDA_FUNCTIONS

} // namespace bda

