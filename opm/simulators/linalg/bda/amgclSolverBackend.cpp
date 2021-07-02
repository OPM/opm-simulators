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

#include <boost/property_tree/json_parser.hpp>


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
    out << "DeviceID: " << deviceID << "\n";
    OpmLog::info(out.str());
    out.str("");
    out.clear();

    A_vals.resize(nnz);
    A_cols.resize(nnz);
    A_rows.resize(N + 1);
    rhs.resize(N);
    x.resize(N);

    // try to read amgcl parameters via json file
    std::string filename = "amgcl_options.json";
    std::ifstream file(filename);
    std::string backend_type_string;

    if (file.is_open()) { // if file exists, read parameters from file
        boost::property_tree::read_json(file, prm);

        backend_type_string = prm.get<std::string>("backend_type"); // defaults to cpu if not specified        
        out << "Using parameters from " << filename << ":\n";
    } else { // otherwise use default parameters, same as Dune
        prm.put("backend_type", "cpu");
        prm.put("precond.class", "relaxation");
        prm.put("precond.type", "ilu0");
        prm.put("precond.damping", 0.9);
        prm.put("solver.type", "bicgstab");
        prm.put("solver.tol", tolerance);
        prm.put("solver.maxiter", maxit);
        prm.put("solver.verbose", verbosity >= 2);
        out << "Using default amgcl parameters:\n";
    }

    boost::property_tree::write_json(out, prm); // print amgcl parameters
    prm.erase("backend_type");                  // delete custom parameter, otherwise amgcl prints a warning

    if (backend_type_string == "cpu") {
        backend_type = Amgcl_backend_type::cpu;
    } else if (backend_type_string == "cuda") {
        backend_type = Amgcl_backend_type::cuda;
    } else if (backend_type_string == "vexcl") {
        backend_type = Amgcl_backend_type::vexcl;
    } else {
        OPM_THROW(std::logic_error, "Error unknown value for amgcl parameter 'backend_type'");
    }

    if (backend_type == Amgcl_backend_type::cuda) {
#if HAVE_CUDA
        cudaDeviceProp prop;
        cudaGetDeviceProperties(&prop, deviceID);
        out << prop.name << std::endl;
        cusparseCreate(&CUDA_bprm.cusparse_handle);
#else
        OPM_THROW(std::logic_error, "Error amgcl is trying to use CUDA, but CUDA was not found by CMake");
#endif
    }
    if (backend_type == Amgcl_backend_type::vexcl) {
#if !HAVE_VEXCL
        OPM_THROW(std::logic_error, "Error amgcl is trying to use VexCL, but VexCL was not found by CMake");
#endif
    }
    OpmLog::info(out.str());

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
        if (backend_type == Amgcl_backend_type::cuda) { // use CUDA
#if HAVE_CUDA
            // create matrix object
            auto A = std::tie(N, A_rows, A_cols, A_vals);

            // create solver and construct preconditioner
            // don't reuse this unless the preconditioner can be reused
            CUDA_Solver solve(A, prm, CUDA_bprm);

            // print solver structure (once)
            std::call_once(print_info, [&](){
                std::ostringstream out;
                out << solve << std::endl;
                OpmLog::info(out.str());
            });

            thrust::device_vector<double> B(b, b + N);
            thrust::device_vector<double> X(N, 0.0);

            // actually solve
            std::tie(iters, error) = solve(B, X);

            thrust::copy(X.begin(), X.end(), x.begin());
#endif
        } else if (backend_type == Amgcl_backend_type::cpu) { // use builtin backend (CPU)
            // create matrix object
            auto Atmp = std::tie(N, A_rows, A_cols, A_vals);
            auto A = amgcl::adapter::block_matrix<dmat_type>(Atmp);

            // create solver and construct preconditioner
            // don't reuse this unless the preconditioner can be reused
            CPU_Solver solve(A, prm);

            // print solver structure (once)
            std::call_once(print_info, [&](){
                std::ostringstream out;
                out << solve << std::endl;
                OpmLog::info(out.str());
            });

            // reset x vector
            std::fill(x.begin(), x.end(), 0.0);

            // create blocked vectors
            auto b_ptr = reinterpret_cast<dvec_type*>(b);
            auto x_ptr = reinterpret_cast<dvec_type*>(x.data());
            auto B = amgcl::make_iterator_range(b_ptr, b_ptr + N / block_size);
            auto X = amgcl::make_iterator_range(x_ptr, x_ptr + N / block_size);

            // actually solve
            std::tie(iters, error) = solve(B, X);
        } else if (backend_type == Amgcl_backend_type::vexcl) {
#if HAVE_VEXCL
            if constexpr(block_size != 1){
                vex::Context ctx(vex::Filter::Env && vex::Filter::Count(1));
                std::cout << ctx << std::endl;

                // allow vexcl to use blocked matrices
                vex::scoped_program_header h1(ctx, amgcl::backend::vexcl_static_matrix_declaration<double, block_size>());

                typedef amgcl::backend::vexcl<dmat_type> Backend;

                typedef amgcl::make_block_solver<
                    amgcl::relaxation::as_preconditioner<Backend, amgcl::relaxation::ilu0>,
                    amgcl::solver::bicgstab<Backend>
                    > Solver;

                typename Backend::params bprm;
                bprm.q = ctx;  // set vexcl context

                auto A = std::tie(N, A_rows, A_cols, A_vals);

                Solver solve(A, prm, bprm);

                auto b_ptr = reinterpret_cast<dvec_type*>(b);
                auto x_ptr = reinterpret_cast<dvec_type*>(x.data());
                vex::vector<dvec_type> B(ctx, N / block_size, b_ptr);
                vex::vector<dvec_type> X(ctx, N / block_size, x_ptr);

                std::tie(iters, error) = solve(B, X);

                vex::copy(X, x_ptr);
            }
#endif
        }
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

