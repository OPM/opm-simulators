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
#include <memory>
#include <cstddef>

#include <opm/simulators/linalg/bda/opencl/opencl.hpp>

namespace Opm
{
namespace Accelerator
{

using spmv_blocked_kernel_type = cl::KernelFunctor<cl::Buffer&, cl::Buffer&, cl::Buffer&, const unsigned int,
                                         const cl::Buffer&, cl::Buffer&, const unsigned int, cl::LocalSpaceArg>;
using spmv_kernel_type = cl::KernelFunctor<cl::Buffer&, cl::Buffer&, cl::Buffer&, const unsigned int,
                                         const cl::Buffer&, cl::Buffer&, cl::LocalSpaceArg>;
using residual_blocked_kernel_type = cl::KernelFunctor<cl::Buffer&, cl::Buffer&, cl::Buffer&, const unsigned int,
                                         cl::Buffer&, const cl::Buffer&, cl::Buffer&, const unsigned int, cl::LocalSpaceArg>;
using residual_kernel_type = cl::KernelFunctor<cl::Buffer&, cl::Buffer&, cl::Buffer&, const unsigned int,
                                         cl::Buffer&, const cl::Buffer&, cl::Buffer&, cl::LocalSpaceArg>;
using ilu_apply1_kernel_type = cl::KernelFunctor<cl::Buffer&, cl::Buffer&, cl::Buffer&, cl::Buffer&, cl::Buffer&, const cl::Buffer&,
                                               cl::Buffer&, cl::Buffer&, const unsigned int, const unsigned int, cl::LocalSpaceArg>;
using ilu_apply2_kernel_type = cl::KernelFunctor<cl::Buffer&, cl::Buffer&, cl::Buffer&, cl::Buffer&, cl::Buffer&, cl::Buffer&,
                                               cl::Buffer&, cl::Buffer&, const unsigned int, const unsigned int, cl::LocalSpaceArg>;
using stdwell_apply_kernel_type = cl::KernelFunctor<cl::Buffer&, cl::Buffer&, cl::Buffer&, cl::Buffer&,
                                                             cl::Buffer&, cl::Buffer&, cl::Buffer&,
                                                             const unsigned int, const unsigned int, cl::Buffer&,
                                                             cl::LocalSpaceArg, cl::LocalSpaceArg, cl::LocalSpaceArg>;
using ilu_decomp_kernel_type = cl::KernelFunctor<const unsigned int, const unsigned int, cl::Buffer&, cl::Buffer&, cl::Buffer&,
                                               cl::Buffer&, cl::Buffer&, cl::Buffer&, const int, cl::LocalSpaceArg>;
using isaiL_kernel_type = cl::KernelFunctor<cl::Buffer&, cl::Buffer&, cl::Buffer&, cl::Buffer&, cl::Buffer&, cl::Buffer&,
                                  cl::Buffer&, cl::Buffer&, cl::Buffer&, const unsigned int>;
using isaiU_kernel_type = cl::KernelFunctor<cl::Buffer&, cl::Buffer&, cl::Buffer&, cl::Buffer&, cl::Buffer&, cl::Buffer&,
                                  cl::Buffer&, cl::Buffer&, cl::Buffer&, cl::Buffer&, cl::Buffer&, const unsigned int>;

class OpenclKernels
{
private:
    static int verbosity;
    static cl::CommandQueue *queue;
    static std::vector<double> tmp;     // used as tmp CPU buffer for dot() and norm()
    static bool initialized;
    static std::size_t preferred_workgroup_size_multiple; // stores CL_KERNEL_PREFERRED_WORK_GROUP_SIZE_MULTIPLE

    static std::unique_ptr<cl::KernelFunctor<cl::Buffer&, cl::Buffer&, cl::Buffer&, const unsigned int, cl::LocalSpaceArg> > dot_k;
    static std::unique_ptr<cl::KernelFunctor<cl::Buffer&, cl::Buffer&, const unsigned int, cl::LocalSpaceArg> > norm_k;
    static std::unique_ptr<cl::KernelFunctor<cl::Buffer&, const double, cl::Buffer&, const unsigned int> > axpy_k;
    static std::unique_ptr<cl::KernelFunctor<cl::Buffer&, const double, const unsigned int> > scale_k;
    static std::unique_ptr<cl::KernelFunctor<const double, cl::Buffer&, cl::Buffer&, cl::Buffer&, const unsigned int> > vmul_k;
    static std::unique_ptr<cl::KernelFunctor<cl::Buffer&, cl::Buffer&, cl::Buffer&, const double, const double, const unsigned int> > custom_k;
    static std::unique_ptr<cl::KernelFunctor<const cl::Buffer&, cl::Buffer&, cl::Buffer&, const unsigned int> > full_to_pressure_restriction_k;
    static std::unique_ptr<cl::KernelFunctor<cl::Buffer&, cl::Buffer&, const unsigned int, const unsigned int> > add_coarse_pressure_correction_k;
    static std::unique_ptr<cl::KernelFunctor<const cl::Buffer&, cl::Buffer&, const cl::Buffer&, const unsigned int> > prolongate_vector_k;
    static std::unique_ptr<spmv_blocked_kernel_type> spmv_blocked_k;
    static std::unique_ptr<spmv_blocked_kernel_type> spmv_blocked_add_k;
    static std::unique_ptr<spmv_kernel_type> spmv_k;
    static std::unique_ptr<spmv_kernel_type> spmv_noreset_k;
    static std::unique_ptr<residual_blocked_kernel_type> residual_blocked_k;
    static std::unique_ptr<residual_kernel_type> residual_k;
    static std::unique_ptr<ilu_apply1_kernel_type> ILU_apply1_k;
    static std::unique_ptr<ilu_apply2_kernel_type> ILU_apply2_k;
    static std::unique_ptr<stdwell_apply_kernel_type> stdwell_apply_k;
    static std::unique_ptr<ilu_decomp_kernel_type> ilu_decomp_k;
    static std::unique_ptr<isaiL_kernel_type> isaiL_k;
    static std::unique_ptr<isaiU_kernel_type> isaiU_k;

    OpenclKernels(){}; // disable instantiation

public:
    static const std::string axpy_str;
    static const std::string scale_str;
    static const std::string vmul_str;
    static const std::string dot_1_str;
    static const std::string norm_str;
    static const std::string custom_str;
    static const std::string full_to_pressure_restriction_str;
    static const std::string add_coarse_pressure_correction_str;
    static const std::string prolongate_vector_str;
    static const std::string spmv_blocked_str;
    static const std::string spmv_blocked_add_str;
    static const std::string spmv_str;
    static const std::string spmv_noreset_str;
    static const std::string residual_blocked_str;
    static const std::string residual_str;
#if CHOW_PATEL
    static const std::string ILU_apply1_str;
    static const std::string ILU_apply2_str;
#else
    static const std::string ILU_apply1_fm_str;
    static const std::string ILU_apply2_fm_str;
#endif
    static const std::string stdwell_apply_str;
    static const std::string ILU_decomp_str;
    static const std::string isaiL_str;
    static const std::string isaiU_str;

    static void init(cl::Context *context, cl::CommandQueue *queue, std::vector<cl::Device>& devices, int verbosity);

    static double dot(cl::Buffer& in1, cl::Buffer& in2, cl::Buffer& out, int N);
    static double norm(cl::Buffer& in, cl::Buffer& out, int N);
    static void axpy(cl::Buffer& in, const double a, cl::Buffer& out, int N);
    static void scale(cl::Buffer& in, const double a, int N);
    static void vmul(const double alpha, cl::Buffer& in1, cl::Buffer& in2, cl::Buffer& out, int N);
    static void custom(cl::Buffer& p, cl::Buffer& v, cl::Buffer& r, const double omega, const double beta, int N);
    static void full_to_pressure_restriction(const cl::Buffer& fine_y, cl::Buffer& weights, cl::Buffer& coarse_y, int Nb);
    static void add_coarse_pressure_correction(cl::Buffer& coarse_x, cl::Buffer& fine_x, int pressure_idx, int Nb);
    static void prolongate_vector(const cl::Buffer& in, cl::Buffer& out, const cl::Buffer& cols, int N);
    static void spmv(cl::Buffer& vals, cl::Buffer& cols, cl::Buffer& rows, const cl::Buffer& x, cl::Buffer& b, int Nb, unsigned int block_size, bool reset = true, bool add = false);
    static void residual(cl::Buffer& vals, cl::Buffer& cols, cl::Buffer& rows, cl::Buffer& x, const cl::Buffer& rhs, cl::Buffer& out, int Nb, unsigned int block_size);

    static void ILU_apply1(cl::Buffer& rowIndices, cl::Buffer& vals, cl::Buffer& cols, cl::Buffer& rows, cl::Buffer& diagIndex,
        const cl::Buffer& y, cl::Buffer& x, cl::Buffer& rowsPerColor, int color, int Nb, unsigned int block_size);

    static void ILU_apply2(cl::Buffer& rowIndices, cl::Buffer& vals, cl::Buffer& cols, cl::Buffer& rows, cl::Buffer& diagIndex,
        cl::Buffer& invDiagVals, cl::Buffer& x, cl::Buffer& rowsPerColor, int color, int Nb, unsigned int block_size);

    static void ILU_decomp(int firstRow, int lastRow, cl::Buffer& rowIndices, cl::Buffer& vals, cl::Buffer& cols, cl::Buffer& rows,
        cl::Buffer& diagIndex, cl::Buffer& invDiagVals, int Nb, unsigned int block_size);

    static void apply_stdwells(cl::Buffer& d_Cnnzs_ocl, cl::Buffer &d_Dnnzs_ocl, cl::Buffer &d_Bnnzs_ocl,
        cl::Buffer &d_Ccols_ocl, cl::Buffer &d_Bcols_ocl, cl::Buffer &d_x, cl::Buffer &d_y,
        int dim, int dim_wells, cl::Buffer &d_val_pointers_ocl, int num_std_wells);

    static void isaiL(cl::Buffer& diagIndex, cl::Buffer& colPointers, cl::Buffer& mapping, cl::Buffer& nvc,
            cl::Buffer& luIdxs, cl::Buffer& xxIdxs, cl::Buffer& dxIdxs, cl::Buffer& LUvals, cl::Buffer& invLvals, unsigned int Nb);

    static void isaiU(cl::Buffer& diagIndex, cl::Buffer& colPointers, cl::Buffer& rowIndices, cl::Buffer& mapping,
            cl::Buffer& nvc, cl::Buffer& luIdxs, cl::Buffer& xxIdxs, cl::Buffer& dxIdxs, cl::Buffer& LUvals,
            cl::Buffer& invDiagVals, cl::Buffer& invUvals, unsigned int Nb);
};

} // namespace Accelerator
} // namespace Opm

#endif
