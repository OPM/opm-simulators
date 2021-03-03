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

#ifndef OPENCL_HPP
#define OPENCL_HPP

#include <string>

#include <opm/simulators/linalg/bda/opencl.hpp>

namespace bda
{

using spmv_kernel_type = cl::make_kernel<cl::Buffer&, cl::Buffer&, cl::Buffer&, const unsigned int,
                                         cl::Buffer&, cl::Buffer&, const unsigned int, cl::LocalSpaceArg>;
using ilu_apply1_kernel_type = cl::make_kernel<cl::Buffer&, cl::Buffer&, cl::Buffer&, cl::Buffer&, cl::Buffer&,
                                               cl::Buffer&, cl::Buffer&, const unsigned int, const unsigned int, cl::LocalSpaceArg>;
using ilu_apply2_kernel_type = cl::make_kernel<cl::Buffer&, cl::Buffer&, cl::Buffer&, cl::Buffer&, cl::Buffer&,
                                               cl::Buffer&, cl::Buffer&, const unsigned int, const unsigned int, cl::LocalSpaceArg>;
using stdwell_apply_kernel_type = cl::make_kernel<cl::Buffer&, cl::Buffer&, cl::Buffer&, cl::Buffer&,
                                                  cl::Buffer&, cl::Buffer&, cl::Buffer&, cl::Buffer&,
                                                  const unsigned int, const unsigned int, cl::Buffer&,
                                                  cl::LocalSpaceArg, cl::LocalSpaceArg, cl::LocalSpaceArg>;
using stdwell_apply_no_reorder_kernel_type = cl::make_kernel<cl::Buffer&, cl::Buffer&, cl::Buffer&, cl::Buffer&,
                                                             cl::Buffer&, cl::Buffer&, cl::Buffer&,
                                                             const unsigned int, const unsigned int, cl::Buffer&,
                                                             cl::LocalSpaceArg, cl::LocalSpaceArg, cl::LocalSpaceArg>;
using ilu_decomp_kernel_type = cl::make_kernel<const unsigned int, const unsigned int, cl::Buffer&, cl::Buffer&,
                                               cl::Buffer&, cl::Buffer&, cl::Buffer&, const int, cl::LocalSpaceArg>;

    std::string get_axpy_string();

    // returns partial sums, instead of the final dot product
    std::string get_dot_1_string();

    // returns partial sums, instead of the final norm
    // the square root must be computed on CPU
    std::string get_norm_string();

    // p = (p - omega * v) * beta + r
    std::string get_custom_string();

    // b = mat * x
    // algorithm based on:
    // Optimization of Block Sparse Matrix-Vector Multiplication on Shared-MemoryParallel Architectures,
    // Ryan Eberhardt, Mark Hoemmen, 2016, https://doi.org/10.1109/IPDPSW.2016.42
    std::string get_spmv_blocked_string();

    // ILU apply part 1: forward substitution
    // solves L*x=y where L is a lower triangular sparse blocked matrix
    std::string get_ILU_apply1_string(bool full_matrix);

    // ILU apply part 2: backward substitution
    // solves U*x=y where L is a lower triangular sparse blocked matrix
    std::string get_ILU_apply2_string(bool full_matrix);

    /// Generate string with the stdwell_apply kernels
    /// If reorder is true, the B/Ccols do not correspond with the x/y vector
    /// the x/y vector is reordered, use toOrder to address that
    /// \param[in] reorder   whether the matrix is reordered or not
    std::string get_stdwell_apply_string(bool reorder);

    std::string get_ilu_decomp_string();

} // end namespace bda

#endif
