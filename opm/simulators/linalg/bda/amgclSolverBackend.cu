/*
  Copyright 2021 Equinor ASA

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

#include <amgcl/backend/cuda.hpp>
#include <amgcl/relaxation/cusparse_ilu0.hpp>
#include <opm/simulators/linalg/bda/amgclSolverBackend.hpp>

/// This file is only compiled when both amgcl and CUDA are found by CMake

namespace Opm
{
namespace Accelerator
{

using Opm::OpmLog;


template <unsigned int block_size>
void amgclSolverBackend<block_size>::solve_cuda(double *b) {
    typedef amgcl::backend::cuda<double> CUDA_Backend;
    typedef amgcl::make_solver<amgcl::runtime::preconditioner<CUDA_Backend>, amgcl::runtime::solver::wrapper<CUDA_Backend> > CUDA_Solver;

    static typename CUDA_Backend::params CUDA_bprm; // amgcl backend parameters, only used for cusparseHandle

    // initialize cusparse handle for amgcl, cannot merge this call_once with 'print solver structure'
    std::call_once(cuda_initialize, [&](){
        cudaDeviceProp prop;
        cudaGetDeviceProperties(&prop, deviceID);
        std::ostringstream out;
        out << prop.name << std::endl;
        OpmLog::info(out.str());
        cusparseCreate(&CUDA_bprm.cusparse_handle);
    });

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
}


#define INSTANTIATE_BDA_FUNCTIONS(n)                      \
template void amgclSolverBackend<n>::solve_cuda(double*); \

INSTANTIATE_BDA_FUNCTIONS(1);
INSTANTIATE_BDA_FUNCTIONS(2);
INSTANTIATE_BDA_FUNCTIONS(3);
INSTANTIATE_BDA_FUNCTIONS(4);
INSTANTIATE_BDA_FUNCTIONS(5);
INSTANTIATE_BDA_FUNCTIONS(6);

#undef INSTANTIATE_BDA_FUNCTIONS

} // namespace Accelerator
} // namespace Opm

