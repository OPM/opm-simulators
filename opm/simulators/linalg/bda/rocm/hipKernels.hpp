/*
  Copyright 2024 Equinor ASA

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

#ifndef HIPKERNELS_HPP
#define HIPKERNELS_HPP

#include <string>
#include <memory>
#include <cstddef>

#include <hip/hip_runtime_api.h>
#include <hip/hip_version.h>

// #include <opm/simulators/linalg/bda/opencl/opencl.hpp>

namespace Opm
{

class HipKernels
{
private:
    static int verbosity;
    static double* tmp;     // used as tmp CPU buffer for dot() and norm()
    static bool initialized;
    static std::size_t preferred_workgroup_size_multiple; // stores CL_KERNEL_PREFERRED_WORK_GROUP_SIZE_MULTIPLE
    
    HipKernels(){}; // disable instantiation
    
public:
    static void init(int verbosity);
    static void full_to_pressure_restriction(const double* fine_y, double* weights, double* coarse_y, int Nb, hipStream_t stream);
    static void add_coarse_pressure_correction(double* coarse_x, double* fine_x, int pressure_idx, int Nb, hipStream_t stream);
    static void vmul(const double alpha, double* in1, double* in2, double* out, int N, hipStream_t stream);
    static void prolongate_vector(const double* in, double* out, const int* cols, int N, hipStream_t stream);
    static void residual(double* vals, int* cols, int* rows, double* x, const double* rhs, double* out, int Nb, unsigned int block_size, hipStream_t stream);
    static void spmv(double* vals, int* cols, int* rows, double* x, double* y, int Nb, unsigned int block_size, hipStream_t stream);
};

} // namespace Opm

#endif
