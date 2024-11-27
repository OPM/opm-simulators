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

#include <opm/simulators/linalg/gpubridge/amgclSolverBackend.hpp>

#include <opm/simulators/linalg/gpubridge/GpuResult.hpp>

#include <opm/common/OpmLog/OpmLog.hpp>
#include <opm/common/ErrorMacros.hpp>
#include <dune/common/timer.hh>

#include <boost/property_tree/json_parser.hpp>

#if HAVE_VEXCL
#include <amgcl/backend/vexcl.hpp>
#include <amgcl/backend/vexcl_static_matrix.hpp>
#endif

#include <algorithm>
#include <fstream>
#include <ios>
#include <ostream>
#include <memory>
#include <mutex>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

namespace Opm::Accelerator {

using Dune::Timer;

template<class Scalar, unsigned int block_size>
amgclSolverBackend<Scalar,block_size>::
amgclSolverBackend(const int          verbosity_,
                   const int          maxit_,
                   const Scalar       tolerance_,
                   const unsigned int platformID_,
                   const unsigned int deviceID_)
    : Base(verbosity_, maxit_, tolerance_, platformID_, deviceID_)
{}

template<class Scalar, unsigned int block_size>
amgclSolverBackend<Scalar,block_size>::~amgclSolverBackend()
{}

template<class Scalar, unsigned int block_size>
void amgclSolverBackend<Scalar,block_size>::initialize(int Nb_, int nnzbs)
{
    this->Nb = Nb_;
    this->N = Nb * block_size;
    this->nnzb = nnzbs;
    this->nnz = nnzbs * block_size * block_size;

    std::ostringstream out;
    out << "Initializing amgclSolverBackend, matrix size: " << Nb
        << " blockrows, nnzb: " << nnzb << " blocks\n";
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
        try {
            boost::property_tree::read_json(file, prm);
        } catch (boost::property_tree::json_parser::json_parser_error& e) {
            OPM_THROW(std::logic_error, "Error cannot parse json file '" + filename + "'");
        }

        // the prm.get reads data from the file, with default values if not specified
        // the prm.put puts the data in the property_tree, so it gets printed
        backend_type_string = prm.get("backend_type", "cpu");
        prm.put("backend_type", backend_type_string);
        std::string t1 = prm.get("precond.class", "relaxation");
        prm.put("precond.class", t1);
        t1 = prm.get("precond.type", "ilu0");
        prm.put("precond.type", t1);
        double t2 = prm.get("precond.damping", 0.9);
        prm.put("precond.damping", t2);
        t1 = prm.get("solver.type", "bicgstab");
        prm.put("solver.type", t1);
        t2 = prm.get("solver.tol", tolerance);
        prm.put("solver.tol", t2);
        int t3 = prm.get("solver.maxiter", maxit);
        prm.put("solver.maxiter", t3);
        bool t4 = prm.get("solver.verbose", verbosity >= 2);
        prm.put("solver.verbose", t4);
        out << "Using parameters from " << filename
            << " (with default values for omitted parameters):\n";
    } else { // otherwise use default parameters, same as Dune
        prm.put("backend_type", "cpu"); // put it in the tree so it gets printed
        prm.put("precond.class", "relaxation");
        prm.put("precond.type", "ilu0");
        prm.put("precond.damping", 0.9);
        prm.put("solver.type", "bicgstab");
        prm.put("solver.tol", tolerance);
        prm.put("solver.maxiter", maxit);
        prm.put("solver.verbose", verbosity >= 2);
        backend_type_string = prm.get("backend_type", "cpu");
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
        OPM_THROW(std::logic_error,
                  "Error unknown value for amgcl parameter 'backend_type', use [cpu|cuda|vexcl]");
    }

    if (backend_type == Amgcl_backend_type::cuda) {
#if !HAVE_CUDA
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

template<class Scalar, unsigned int block_size>
void amgclSolverBackend<Scalar,block_size>::
convert_sparsity_pattern(int* rows, int* cols)
{
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

template<class Scalar, unsigned int block_size>
void amgclSolverBackend<Scalar,block_size>::
convert_data(Scalar* vals, int* rows)
{
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

#if HAVE_VEXCL
void initialize_vexcl(std::vector<cl::CommandQueue>& ctx,
                      unsigned int platformID, unsigned int deviceID)
{
    std::vector<cl::Platform> platforms;
    std::vector<cl::Device> devices;
    cl::Platform::get(&platforms);

    if (platforms.size() <= platformID) {
        OPM_THROW(std::logic_error, "Error chosen too high OpenCL platform ID");
    }

    std::string platform_name, device_name;
    platforms[platformID].getInfo(CL_PLATFORM_NAME, &platform_name);
    platforms[platformID].getDevices(CL_DEVICE_TYPE_ALL, &devices);

    if (devices.size() <= deviceID){
        OPM_THROW(std::logic_error, "Error chosen too high OpenCL device ID");
    }

    devices[deviceID].getInfo(CL_DEVICE_NAME, &device_name);

    cl::Context c(devices[deviceID]);
    cl::CommandQueue q(c, devices[deviceID]);
    ctx.push_back(q);

    std::ostringstream out;
    out << "Using VexCL on " << device_name << " (" << platform_name << ")\n";
    OpmLog::info(out.str());
}

template <typename vexcl_matrix_type, typename vexcl_vector_type,
          unsigned int block_size, typename Scalar, typename AIJInfo>
void solve_vexcl(const AIJInfo& A,
                 const boost::property_tree::ptree prm,
                 const std::vector<cl::CommandQueue>& ctx,
                 Scalar* b,
                 std::vector<Scalar>& x,
                 const int N,
                 int& iters,
                 Scalar& error)
{
    using Backend = amgcl::backend::vexcl<vexcl_matrix_type>;
    using Solver = amgcl::make_solver<amgcl::runtime::preconditioner<Backend>,
                                      amgcl::runtime::solver::wrapper<Backend>>;

    typename Solver::backend_params bprm;
    bprm.q = ctx;  // set vexcl context

    Solver solve(A, prm, bprm); // create solver

    auto b_ptr = reinterpret_cast<vexcl_vector_type*>(b);
    auto x_ptr = reinterpret_cast<vexcl_vector_type*>(x.data());
    vex::vector<vexcl_vector_type> B(ctx, N / block_size, b_ptr);
    vex::vector<vexcl_vector_type> X(ctx, N / block_size, x_ptr);

    std::tie(iters, error) = solve(B, X); // actually perform solve

    vex::copy(X, x_ptr);
}
#endif

template<class Scalar, unsigned int block_size>
void amgclSolverBackend<Scalar,block_size>::
solve_system(Scalar* b, GpuResult& res)
{
    Timer t;

    try {
        if (backend_type == Amgcl_backend_type::cuda) { // use CUDA
#if HAVE_CUDA
            solve_cuda(b);
#endif
        } else if (backend_type == Amgcl_backend_type::cpu) { // use builtin backend (CPU)
            // create matrix object
            auto Atmp = std::tie(N, A_rows, A_cols, A_vals);
            auto print = [this](const auto& solve)
            {
                // print solver structure (once)
                std::call_once(print_info, [&](){
                    std::ostringstream out;
                    out << solve << std::endl;
                    OpmLog::info(out.str());
                });
            };
            if constexpr (block_size == 1) {
                // create solver and construct preconditioner
                // don't reuse this unless the preconditioner can be reused
                CPU_Solver solve(Atmp, prm);

                print(solve);

                // reset x vector
                std::fill(x.begin(), x.end(), 0.0);

                std::vector<Scalar> b_(b, b + N);

                // create numa vectors
                typename CPU_Backend::params bprm;
                auto B = CPU_Backend::copy_vector(b_, bprm);
                auto X = CPU_Backend::copy_vector(x, bprm);

                // actually solve
                std::tie(iters, error) = solve(*B, *X);

                std::copy(&(*X)[0], &(*X)[0] + N, x.begin());
            } else {
                // create matrix object
                auto A = amgcl::adapter::block_matrix<dmat_type>(Atmp);

                // create solver and construct preconditioner
                // don't reuse this unless the preconditioner can be reused
                CPU_Solver solve(A, prm);

                // print solver structure (once)
                print(solve);

                // reset x vector
                std::fill(x.begin(), x.end(), 0.0);

                // create blocked vectors
                auto b_ptr = reinterpret_cast<dvec_type*>(b);
                auto x_ptr = reinterpret_cast<dvec_type*>(x.data());
                auto B = amgcl::make_iterator_range(b_ptr, b_ptr + N / block_size);
                auto X = amgcl::make_iterator_range(x_ptr, x_ptr + N / block_size);

                // actually solve
                std::tie(iters, error) = solve(B, X);
            }
        } else if (backend_type == Amgcl_backend_type::vexcl) {
#if HAVE_VEXCL
            static std::vector<cl::CommandQueue> ctx; // using CommandQueue directly instead of vex::Context
            std::call_once(vexcl_initialize, [&](){
                initialize_vexcl(ctx, platformID, deviceID);
            });
            if constexpr(block_size == 1){
                auto A = std::tie(N, A_rows, A_cols, A_vals);

                solve_vexcl<Scalar, Scalar, block_size>(A, prm, ctx, b, x, N, iters, error);
            } else {
                // allow vexcl to use blocked matrices
                vex::scoped_program_header h1(ctx,
                                              amgcl::backend::vexcl_static_matrix_declaration<Scalar, block_size>());

                auto Atmp = std::tie(N, A_rows, A_cols, A_vals);
                auto A = amgcl::adapter::block_matrix<dmat_type>(Atmp);

                solve_vexcl<dmat_type, dvec_type, block_size>(A, prm, ctx, b, x, N, iters, error);
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
        out << "=== converged: " << res.converged << ", time: " << res.elapsed
            << ", time per iteration: " << res.elapsed / iters << ", iterations: " << iters;
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
template<class Scalar, unsigned int block_size>
void amgclSolverBackend<Scalar,block_size>::get_result(Scalar* x_)
{
    Timer t;

    std::copy(x.begin(), x.end(), x_);

    if (verbosity >= 3) {
        std::ostringstream out;
        out << "amgclSolverBackend::get_result(): " << t.stop() << " s";
        OpmLog::info(out.str());
    }
} // end get_result()

template<class Scalar, unsigned int block_size>
SolverStatus amgclSolverBackend<Scalar,block_size>::
solve_system(std::shared_ptr<BlockedMatrix<Scalar>> matrix,
             Scalar* b,
             [[maybe_unused]] std::shared_ptr<BlockedMatrix<Scalar>> jacMatrix,
             [[maybe_unused]] WellContributions<Scalar>& wellContribs,
             GpuResult& res)
{
    if (initialized == false) {
        initialize(matrix->Nb, matrix->nnzbs);
        convert_sparsity_pattern(matrix->rowPointers, matrix->colIndices);
    }
    convert_data(matrix->nnzValues, matrix->rowPointers);
    solve_system(b, res);
    return SolverStatus::GPU_SOLVER_SUCCESS;
}

#define INSTANTIATE_TYPE(T)                 \
    template class amgclSolverBackend<T,1>; \
    template class amgclSolverBackend<T,2>; \
    template class amgclSolverBackend<T,3>; \
    template class amgclSolverBackend<T,4>; \
    template class amgclSolverBackend<T,5>; \
    template class amgclSolverBackend<T,6>;

INSTANTIATE_TYPE(double)

} // namespace Opm::Accelerator
