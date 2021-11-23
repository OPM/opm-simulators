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

#include <opm/simulators/linalg/bda/opencl.hpp>

namespace Opm
{
namespace Accelerator
{

using spmv_blocked_kernel_type = cl::KernelFunctor<cl::Buffer&, cl::Buffer&, cl::Buffer&, const unsigned int,
                                         cl::Buffer&, cl::Buffer&, const unsigned int, cl::LocalSpaceArg>;
using spmv_kernel_type = cl::KernelFunctor<cl::Buffer&, cl::Buffer&, cl::Buffer&, const unsigned int,
                                         cl::Buffer&, cl::Buffer&, cl::LocalSpaceArg>;
using residual_blocked_kernel_type = cl::KernelFunctor<cl::Buffer&, cl::Buffer&, cl::Buffer&, const unsigned int,
                                         cl::Buffer&, const cl::Buffer&, cl::Buffer&, const unsigned int, cl::LocalSpaceArg>;
using residual_kernel_type = cl::KernelFunctor<cl::Buffer&, cl::Buffer&, cl::Buffer&, const unsigned int,
                                         cl::Buffer&, const cl::Buffer&, cl::Buffer&, cl::LocalSpaceArg>;
using ilu_apply1_kernel_type = cl::KernelFunctor<cl::Buffer&, cl::Buffer&, cl::Buffer&, cl::Buffer&, const cl::Buffer&,
                                               cl::Buffer&, cl::Buffer&, const unsigned int, const unsigned int, cl::LocalSpaceArg>;
using ilu_apply2_kernel_type = cl::KernelFunctor<cl::Buffer&, cl::Buffer&, cl::Buffer&, cl::Buffer&, cl::Buffer&,
                                               cl::Buffer&, cl::Buffer&, const unsigned int, const unsigned int, cl::LocalSpaceArg>;
using stdwell_apply_kernel_type = cl::KernelFunctor<cl::Buffer&, cl::Buffer&, cl::Buffer&, cl::Buffer&,
                                                  cl::Buffer&, cl::Buffer&, cl::Buffer&, cl::Buffer&,
                                                  const unsigned int, const unsigned int, cl::Buffer&,
                                                  cl::LocalSpaceArg, cl::LocalSpaceArg, cl::LocalSpaceArg>;
using stdwell_apply_no_reorder_kernel_type = cl::KernelFunctor<cl::Buffer&, cl::Buffer&, cl::Buffer&, cl::Buffer&,
                                                             cl::Buffer&, cl::Buffer&, cl::Buffer&,
                                                             const unsigned int, const unsigned int, cl::Buffer&,
                                                             cl::LocalSpaceArg, cl::LocalSpaceArg, cl::LocalSpaceArg>;
using ilu_decomp_kernel_type = cl::KernelFunctor<const unsigned int, const unsigned int, cl::Buffer&, cl::Buffer&,
                                               cl::Buffer&, cl::Buffer&, cl::Buffer&, const int, cl::LocalSpaceArg>;

class OpenclKernels
{
private:
    static int verbosity;
    static cl::CommandQueue *queue;
    static std::vector<double> tmp;     // used as tmp CPU buffer for dot() and norm()
    static bool initialized;

    enum matrix_operation {
        spmv_op,
        residual_op
    };

    static std::unique_ptr<cl::KernelFunctor<cl::Buffer&, cl::Buffer&, cl::Buffer&, const unsigned int, cl::LocalSpaceArg> > dot_k;
    static std::unique_ptr<cl::KernelFunctor<cl::Buffer&, cl::Buffer&, const unsigned int, cl::LocalSpaceArg> > norm_k;
    static std::unique_ptr<cl::KernelFunctor<cl::Buffer&, const double, cl::Buffer&, const unsigned int> > axpy_k;
    static std::unique_ptr<cl::KernelFunctor<cl::Buffer&, const double, const unsigned int> > scale_k;
    static std::unique_ptr<cl::KernelFunctor<const double, cl::Buffer&, cl::Buffer&, cl::Buffer&, const unsigned int> > vmul_k;
    static std::unique_ptr<cl::KernelFunctor<cl::Buffer&, cl::Buffer&, cl::Buffer&, const double, const double, const unsigned int> > custom_k;
    static std::unique_ptr<cl::KernelFunctor<const cl::Buffer&, cl::Buffer&, cl::Buffer&, const unsigned int> > full_to_pressure_restriction_k;
    static std::unique_ptr<cl::KernelFunctor<cl::Buffer&, cl::Buffer&, const unsigned int, const unsigned int> > add_coarse_pressure_correction_k;
    static std::unique_ptr<spmv_blocked_kernel_type> spmv_blocked_k;
    static std::unique_ptr<spmv_kernel_type> spmv_k;
    static std::unique_ptr<spmv_kernel_type> spmv_noreset_k;
    static std::unique_ptr<residual_blocked_kernel_type> residual_blocked_k;
    static std::unique_ptr<residual_kernel_type> residual_k;
    static std::unique_ptr<ilu_apply1_kernel_type> ILU_apply1_k;
    static std::unique_ptr<ilu_apply2_kernel_type> ILU_apply2_k;
    static std::unique_ptr<stdwell_apply_kernel_type> stdwell_apply_k;
    static std::unique_ptr<stdwell_apply_no_reorder_kernel_type> stdwell_apply_no_reorder_k;
    static std::unique_ptr<ilu_decomp_kernel_type> ilu_decomp_k;

    /// Generate string with axpy kernel
    /// a = a + alpha * b
    static std::string get_axpy_source();

    /// Generate string with scale kernel
    /// a = a * alpha
    static std::string get_scale_source();

    /// multiply vector with another vector and a scalar, element-wise
    /// add result to a third vector
    static std::string get_vmul_source();

    /// returns partial sums, instead of the final dot product
    /// partial sums are added on CPU
    static std::string get_dot_1_source();

    /// returns partial sums, instead of the final norm
    /// the square root must be computed on CPU
    static std::string get_norm_source();

    /// Generate string with custom kernel
    /// This kernel combines some ilubicgstab vector operations into 1
    /// p = (p - omega * v) * beta + r
    static std::string get_custom_source();

    /// Transform blocked vector to scalar vector using pressure-weights
    static std::string get_full_to_pressure_restriction_source();

    /// Add the coarse pressure solution back to the finer, complete solution
    static std::string get_add_coarse_pressure_correction_source();

    /// b = mat * x
    /// algorithm based on:
    /// Optimization of Block Sparse Matrix-Vector Multiplication on Shared-MemoryParallel Architectures,
    /// Ryan Eberhardt, Mark Hoemmen, 2016, https://doi.org/10.1109/IPDPSW.2016.42
    /// or
    /// res = rhs - (mat * x)
    static std::string get_blocked_matrix_operation_source(matrix_operation op);
    static std::string get_matrix_operation_source(matrix_operation op, bool spmv_reset = true);

    /// ILU apply part 1: forward substitution
    /// solves L*x=y where L is a lower triangular sparse blocked matrix
    /// this L can be it's own BSR matrix (if full_matrix is false),
    /// or it can be inside a normal, square matrix, in that case diagIndex indicates where the rows of L end
    /// \param[in] full_matrix   whether the kernel should operate on a full (square) matrix or not
    static std::string get_ILU_apply1_source(bool full_matrix);

    /// ILU apply part 2: backward substitution
    /// solves U*x=y where U is an upper triangular sparse blocked matrix
    /// this U can be it's own BSR matrix (if full_matrix is false),
    /// or it can be inside a normal, square matrix, in that case diagIndex indicates where the rows of U start
    /// \param[in] full_matrix   whether the kernel should operate on a full (square) matrix or not
    static std::string get_ILU_apply2_source(bool full_matrix);

    /// Generate string with the stdwell_apply kernels
    /// If reorder is true, the B/Ccols do not correspond with the x/y vector
    /// the x/y vector is reordered, use toOrder to address that
    /// \param[in] reorder   whether the matrix is reordered or not
    static std::string get_stdwell_apply_source(bool reorder);

    /// Generate string with the exact ilu decomposition kernel
    /// The kernel takes a full BSR matrix and performs inplace ILU decomposition
    static std::string get_ilu_decomp_source();

    OpenclKernels(){}; // disable instantiation

public:
    static void init(cl::Context *context, cl::CommandQueue *queue, std::vector<cl::Device>& devices, int verbosity);

    static double dot(cl::Buffer& in1, cl::Buffer& in2, cl::Buffer& out, int N);
    static double norm(cl::Buffer& in, cl::Buffer& out, int N);
    static void axpy(cl::Buffer& in, const double a, cl::Buffer& out, int N);
    static void scale(cl::Buffer& in, const double a, int N);
    static void vmul(const double alpha, cl::Buffer& in1, cl::Buffer& in2, cl::Buffer& out, int N);
    static void custom(cl::Buffer& p, cl::Buffer& v, cl::Buffer& r, const double omega, const double beta, int N);
    static void full_to_pressure_restriction(const cl::Buffer& fine_y, cl::Buffer& weights, cl::Buffer& coarse_y, int Nb);
    static void add_coarse_pressure_correction(cl::Buffer& coarse_x, cl::Buffer& fine_x, int pressure_idx, int Nb);
    static void spmv(cl::Buffer& vals, cl::Buffer& cols, cl::Buffer& rows, cl::Buffer& x, cl::Buffer& b, int Nb, unsigned int block_size, bool reset = true);
    static void residual(cl::Buffer& vals, cl::Buffer& cols, cl::Buffer& rows, cl::Buffer& x, const cl::Buffer& rhs, cl::Buffer& out, int Nb, unsigned int block_size);

    static void ILU_apply1(cl::Buffer& vals, cl::Buffer& cols, cl::Buffer& rows, cl::Buffer& diagIndex,
        const cl::Buffer& y, cl::Buffer& x, cl::Buffer& rowsPerColor, int color, int Nb, unsigned int block_size);

    static void ILU_apply2(cl::Buffer& vals, cl::Buffer& cols, cl::Buffer& rows, cl::Buffer& diagIndex,
        cl::Buffer& invDiagVals, cl::Buffer& x, cl::Buffer& rowsPerColor, int color, int Nb, unsigned int block_size);

    static void ILU_decomp(int firstRow, int lastRow, cl::Buffer& vals, cl::Buffer& cols, cl::Buffer& rows,
        cl::Buffer& diagIndex, cl::Buffer& invDiagVals, int Nb, unsigned int block_size);

    static void apply_stdwells_reorder(cl::Buffer& d_Cnnzs_ocl, cl::Buffer &d_Dnnzs_ocl, cl::Buffer &d_Bnnzs_ocl,
        cl::Buffer &d_Ccols_ocl, cl::Buffer &d_Bcols_ocl, cl::Buffer &d_x, cl::Buffer &d_y,
        cl::Buffer &d_toOrder, int dim, int dim_wells, cl::Buffer &d_val_pointers_ocl, int num_std_wells);

    static void apply_stdwells_no_reorder(cl::Buffer& d_Cnnzs_ocl, cl::Buffer &d_Dnnzs_ocl, cl::Buffer &d_Bnnzs_ocl,
        cl::Buffer &d_Ccols_ocl, cl::Buffer &d_Bcols_ocl, cl::Buffer &d_x, cl::Buffer &d_y,
        int dim, int dim_wells, cl::Buffer &d_val_pointers_ocl, int num_std_wells);
};

} // namespace Accelerator
} // namespace Opm

#endif
