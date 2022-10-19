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
#include <cmath>
#include <sstream>

#include <opm/common/OpmLog/OpmLog.hpp>
#include <opm/common/ErrorMacros.hpp>
#include <dune/common/timer.hh>

#include <opm/simulators/linalg/bda/opencl/openclKernels.hpp>
#include <opm/simulators/linalg/bda/opencl/ChowPatelIlu.hpp>  // defines CHOW_PATEL

namespace Opm
{
namespace Accelerator
{

using Opm::OpmLog;
using Dune::Timer;

// define static variables and kernels
int OpenclKernels::verbosity;
cl::CommandQueue *OpenclKernels::queue;
std::vector<double> OpenclKernels::tmp;
bool OpenclKernels::initialized = false;
std::size_t OpenclKernels::preferred_workgroup_size_multiple = 0;

std::unique_ptr<cl::KernelFunctor<cl::Buffer&, cl::Buffer&, cl::Buffer&, const unsigned int, cl::LocalSpaceArg> > OpenclKernels::dot_k;
std::unique_ptr<cl::KernelFunctor<cl::Buffer&, cl::Buffer&, const unsigned int, cl::LocalSpaceArg> > OpenclKernels::norm_k;
std::unique_ptr<cl::KernelFunctor<cl::Buffer&, const double, cl::Buffer&, const unsigned int> > OpenclKernels::axpy_k;
std::unique_ptr<cl::KernelFunctor<cl::Buffer&, const double, const unsigned int> > OpenclKernels::scale_k;
std::unique_ptr<cl::KernelFunctor<const double, cl::Buffer&, cl::Buffer&, cl::Buffer&, const unsigned int> > OpenclKernels::vmul_k;
std::unique_ptr<cl::KernelFunctor<cl::Buffer&, cl::Buffer&, cl::Buffer&, const double, const double, const unsigned int> > OpenclKernels::custom_k;
std::unique_ptr<cl::KernelFunctor<const cl::Buffer&, cl::Buffer&, cl::Buffer&, const unsigned int> > OpenclKernels::full_to_pressure_restriction_k;
std::unique_ptr<cl::KernelFunctor<cl::Buffer&, cl::Buffer&, const unsigned int, const unsigned int> > OpenclKernels::add_coarse_pressure_correction_k;
std::unique_ptr<cl::KernelFunctor<const cl::Buffer&, cl::Buffer&, const cl::Buffer&, const unsigned int> > OpenclKernels::prolongate_vector_k;
std::unique_ptr<spmv_blocked_kernel_type> OpenclKernels::spmv_blocked_k;
std::unique_ptr<spmv_blocked_kernel_type> OpenclKernels::spmv_blocked_add_k;
std::unique_ptr<spmv_kernel_type> OpenclKernels::spmv_k;
std::unique_ptr<spmv_kernel_type> OpenclKernels::spmv_noreset_k;
std::unique_ptr<residual_blocked_kernel_type> OpenclKernels::residual_blocked_k;
std::unique_ptr<residual_kernel_type> OpenclKernels::residual_k;
std::unique_ptr<ilu_apply1_kernel_type> OpenclKernels::ILU_apply1_k;
std::unique_ptr<ilu_apply2_kernel_type> OpenclKernels::ILU_apply2_k;
std::unique_ptr<stdwell_apply_kernel_type> OpenclKernels::stdwell_apply_k;
std::unique_ptr<ilu_decomp_kernel_type> OpenclKernels::ilu_decomp_k;
std::unique_ptr<isaiL_kernel_type> OpenclKernels::isaiL_k;
std::unique_ptr<isaiU_kernel_type> OpenclKernels::isaiU_k;

// divide A by B, and round up: return (int)ceil(A/B)
unsigned int ceilDivision(const unsigned int A, const unsigned int B)
{
    return A / B + (A % B > 0);
}

void OpenclKernels::init(cl::Context *context, cl::CommandQueue *queue_, std::vector<cl::Device>& devices, int verbosity_)
{
    if (initialized) {
        OpmLog::debug("Warning OpenclKernels is already initialized");
        return;
    }

    queue = queue_;
    verbosity = verbosity_;

    cl::Program::Sources sources;
    sources.emplace_back(axpy_str);
    sources.emplace_back(scale_str);
    sources.emplace_back(vmul_str);
    sources.emplace_back(dot_1_str);
    sources.emplace_back(norm_str);
    sources.emplace_back(custom_str);
    sources.emplace_back(full_to_pressure_restriction_str);
    sources.emplace_back(add_coarse_pressure_correction_str);
    sources.emplace_back(prolongate_vector_str);
    sources.emplace_back(spmv_blocked_str);
    sources.emplace_back(spmv_blocked_add_str);
    sources.emplace_back(spmv_str);
    sources.emplace_back(spmv_noreset_str);
    sources.emplace_back(residual_blocked_str);
    sources.emplace_back(residual_str);
#if CHOW_PATEL
    sources.emplace_back(ILU_apply1_str);
    sources.emplace_back(ILU_apply2_str);
#else
    sources.emplace_back(ILU_apply1_fm_str);
    sources.emplace_back(ILU_apply2_fm_str);
#endif
    sources.emplace_back(stdwell_apply_str);
    sources.emplace_back(ILU_decomp_str);
    sources.emplace_back(isaiL_str);
    sources.emplace_back(isaiU_str);

    cl::Program program = cl::Program(*context, sources);
    program.build(devices);

    // queue.enqueueNDRangeKernel() is a blocking/synchronous call, at least for NVIDIA
    // cl::KernelFunctor<> myKernel(); myKernel(args, arg1, arg2); is also blocking

    // actually creating the kernels
    dot_k.reset(new cl::KernelFunctor<cl::Buffer&, cl::Buffer&, cl::Buffer&, const unsigned int, cl::LocalSpaceArg>(cl::Kernel(program, "dot_1")));
    norm_k.reset(new cl::KernelFunctor<cl::Buffer&, cl::Buffer&, const unsigned int, cl::LocalSpaceArg>(cl::Kernel(program, "norm")));
    axpy_k.reset(new cl::KernelFunctor<cl::Buffer&, const double, cl::Buffer&, const unsigned int>(cl::Kernel(program, "axpy")));
    scale_k.reset(new cl::KernelFunctor<cl::Buffer&, const double, const unsigned int>(cl::Kernel(program, "scale")));
    vmul_k.reset(new cl::KernelFunctor<const double, cl::Buffer&, cl::Buffer&, cl::Buffer&, const unsigned int>(cl::Kernel(program, "vmul")));
    custom_k.reset(new cl::KernelFunctor<cl::Buffer&, cl::Buffer&, cl::Buffer&, const double, const double, const unsigned int>(cl::Kernel(program, "custom")));
    full_to_pressure_restriction_k.reset(new cl::KernelFunctor<const cl::Buffer&, cl::Buffer&, cl::Buffer&, const unsigned int>(cl::Kernel(program, "full_to_pressure_restriction")));
    add_coarse_pressure_correction_k.reset(new cl::KernelFunctor<cl::Buffer&, cl::Buffer&, const unsigned int, const unsigned int>(cl::Kernel(program, "add_coarse_pressure_correction")));
    prolongate_vector_k.reset(new cl::KernelFunctor<const cl::Buffer&, cl::Buffer&, const cl::Buffer&, const unsigned int>(cl::Kernel(program, "prolongate_vector")));
    spmv_blocked_k.reset(new spmv_blocked_kernel_type(cl::Kernel(program, "spmv_blocked")));
    spmv_blocked_add_k.reset(new spmv_blocked_kernel_type(cl::Kernel(program, "spmv_blocked_add")));
    spmv_k.reset(new spmv_kernel_type(cl::Kernel(program, "spmv")));
    spmv_noreset_k.reset(new spmv_kernel_type(cl::Kernel(program, "spmv_noreset")));
    residual_blocked_k.reset(new residual_blocked_kernel_type(cl::Kernel(program, "residual_blocked")));
    residual_k.reset(new residual_kernel_type(cl::Kernel(program, "residual")));
    ILU_apply1_k.reset(new ilu_apply1_kernel_type(cl::Kernel(program, "ILU_apply1")));
    ILU_apply2_k.reset(new ilu_apply2_kernel_type(cl::Kernel(program, "ILU_apply2")));
    stdwell_apply_k.reset(new stdwell_apply_kernel_type(cl::Kernel(program, "stdwell_apply")));
    ilu_decomp_k.reset(new ilu_decomp_kernel_type(cl::Kernel(program, "ilu_decomp")));
    isaiL_k.reset(new isaiL_kernel_type(cl::Kernel(program, "isaiL")));
    isaiU_k.reset(new isaiU_kernel_type(cl::Kernel(program, "isaiU")));

    // testing shows all kernels have the same preferred_workgroup_size_multiple
    // 32 for NVIDIA
    // 64 for AMD
    cl::Kernel(program, "ILU_apply1").getWorkGroupInfo(devices[0], CL_KERNEL_PREFERRED_WORK_GROUP_SIZE_MULTIPLE, &preferred_workgroup_size_multiple);

    initialized = true;
} // end get_opencl_kernels()

double OpenclKernels::dot(cl::Buffer& in1, cl::Buffer& in2, cl::Buffer& out, int N)
{
    const unsigned int work_group_size = 256;
    const unsigned int num_work_groups = ceilDivision(N, work_group_size);
    const unsigned int total_work_items = num_work_groups * work_group_size;
    const unsigned int lmem_per_work_group = sizeof(double) * work_group_size;
    Timer t_dot;
    tmp.resize(num_work_groups);

    cl::Event event = (*dot_k)(cl::EnqueueArgs(*queue, cl::NDRange(total_work_items), cl::NDRange(work_group_size)), in1, in2, out, N, cl::Local(lmem_per_work_group));

    queue->enqueueReadBuffer(out, CL_TRUE, 0, sizeof(double) * num_work_groups, tmp.data());

    double gpu_sum = 0.0;
    for (unsigned int i = 0; i < num_work_groups; ++i) {
        gpu_sum += tmp[i];
    }

    if (verbosity >= 4) {
        event.wait();
        std::ostringstream oss;
        oss << std::scientific << "OpenclKernels dot() time: " << t_dot.stop() << " s";
        OpmLog::info(oss.str());
    }

    return gpu_sum;
}

double OpenclKernels::norm(cl::Buffer& in, cl::Buffer& out, int N)
{
    const unsigned int work_group_size = 256;
    const unsigned int num_work_groups = ceilDivision(N, work_group_size);
    const unsigned int total_work_items = num_work_groups * work_group_size;
    const unsigned int lmem_per_work_group = sizeof(double) * work_group_size;
    Timer t_norm;
    tmp.resize(num_work_groups);

    cl::Event event = (*norm_k)(cl::EnqueueArgs(*queue, cl::NDRange(total_work_items), cl::NDRange(work_group_size)), in, out, N, cl::Local(lmem_per_work_group));

    queue->enqueueReadBuffer(out, CL_TRUE, 0, sizeof(double) * num_work_groups, tmp.data());

    double gpu_norm = 0.0;
    for (unsigned int i = 0; i < num_work_groups; ++i) {
        gpu_norm += tmp[i];
    }
    gpu_norm = sqrt(gpu_norm);

    if (verbosity >= 4) {
        event.wait();
        std::ostringstream oss;
        oss << std::scientific << "OpenclKernels norm() time: " << t_norm.stop() << " s";
        OpmLog::info(oss.str());
    }

    return gpu_norm;
}

void OpenclKernels::axpy(cl::Buffer& in, const double a, cl::Buffer& out, int N)
{
    const unsigned int work_group_size = 32;
    const unsigned int num_work_groups = ceilDivision(N, work_group_size);
    const unsigned int total_work_items = num_work_groups * work_group_size;
    Timer t_axpy;

    cl::Event event = (*axpy_k)(cl::EnqueueArgs(*queue, cl::NDRange(total_work_items), cl::NDRange(work_group_size)), in, a, out, N);

    if (verbosity >= 4) {
        event.wait();
        std::ostringstream oss;
        oss << std::scientific << "OpenclKernels axpy() time: " << t_axpy.stop() << " s";
        OpmLog::info(oss.str());
    }
}

void OpenclKernels::scale(cl::Buffer& in, const double a, int N)
{
    const unsigned int work_group_size = 32;
    const unsigned int num_work_groups = ceilDivision(N, work_group_size);
    const unsigned int total_work_items = num_work_groups * work_group_size;
    Timer t_scale;

    cl::Event event = (*scale_k)(cl::EnqueueArgs(*queue, cl::NDRange(total_work_items), cl::NDRange(work_group_size)), in, a, N);

    if (verbosity >= 4) {
        event.wait();
        std::ostringstream oss;
        oss << std::scientific << "OpenclKernels scale() time: " << t_scale.stop() << " s";
        OpmLog::info(oss.str());
    }
}

void OpenclKernels::vmul(const double alpha, cl::Buffer& in1, cl::Buffer& in2, cl::Buffer& out, int N)
{
    const unsigned int work_group_size = 32;
    const unsigned int num_work_groups = ceilDivision(N, work_group_size);
    const unsigned int total_work_items = num_work_groups * work_group_size;
    Timer t_vmul;

    cl::Event event = (*vmul_k)(cl::EnqueueArgs(*queue, cl::NDRange(total_work_items), cl::NDRange(work_group_size)), alpha, in1, in2, out, N);

    if (verbosity >= 4) {
        event.wait();
        std::ostringstream oss;
        oss << std::scientific << "OpenclKernels vmul() time: " << t_vmul.stop() << " s";
        OpmLog::info(oss.str());
    }
}

void OpenclKernels::custom(cl::Buffer& p, cl::Buffer& v, cl::Buffer& r,
                           const double omega, const double beta, int N)
{
    const unsigned int work_group_size = 32;
    const unsigned int num_work_groups = ceilDivision(N, work_group_size);
    const unsigned int total_work_items = num_work_groups * work_group_size;
    Timer t_custom;

    cl::Event event = (*custom_k)(cl::EnqueueArgs(*queue, cl::NDRange(total_work_items), cl::NDRange(work_group_size)), p, v, r, omega, beta, N);

    if (verbosity >= 4) {
        event.wait();
        std::ostringstream oss;
        oss << std::scientific << "OpenclKernels custom() time: " << t_custom.stop() << " s";
        OpmLog::info(oss.str());
    }
}

void OpenclKernels::full_to_pressure_restriction(const cl::Buffer& fine_y, cl::Buffer& weights, cl::Buffer& coarse_y, int Nb)
{
    const unsigned int work_group_size = 32;
    const unsigned int num_work_groups = ceilDivision(Nb, work_group_size);
    const unsigned int total_work_items = num_work_groups * work_group_size;
    Timer t;

    cl::Event event = (*full_to_pressure_restriction_k)(cl::EnqueueArgs(*queue, cl::NDRange(total_work_items), cl::NDRange(work_group_size)), fine_y, weights, coarse_y, Nb);

    if (verbosity >= 4) {
        event.wait();
        std::ostringstream oss;
        oss << std::scientific << "OpenclKernels full_to_pressure_restriction() time: " << t.stop() << " s";
        OpmLog::info(oss.str());
    }
}

void OpenclKernels::add_coarse_pressure_correction(cl::Buffer& coarse_x, cl::Buffer& fine_x, int pressure_idx, int Nb)
{
    const unsigned int work_group_size = 32;
    const unsigned int num_work_groups = ceilDivision(Nb, work_group_size);
    const unsigned int total_work_items = num_work_groups * work_group_size;
    Timer t;

    cl::Event event = (*add_coarse_pressure_correction_k)(cl::EnqueueArgs(*queue, cl::NDRange(total_work_items), cl::NDRange(work_group_size)), coarse_x, fine_x, pressure_idx, Nb);

    if (verbosity >= 4) {
        event.wait();
        std::ostringstream oss;
        oss << std::scientific << "OpenclKernels add_coarse_pressure_correction() time: " << t.stop() << " s";
        OpmLog::info(oss.str());
    }
}

void OpenclKernels::prolongate_vector(const cl::Buffer& in, cl::Buffer& out, const cl::Buffer& cols, int N)
{
    const unsigned int work_group_size = 32;
    const unsigned int num_work_groups = ceilDivision(N, work_group_size);
    const unsigned int total_work_items = num_work_groups * work_group_size;
    Timer t;

    cl::Event event = (*prolongate_vector_k)(cl::EnqueueArgs(*queue, cl::NDRange(total_work_items), cl::NDRange(work_group_size)), in, out, cols, N);

    if (verbosity >= 4) {
        event.wait();
        std::ostringstream oss;
        oss << std::scientific << "OpenclKernels prolongate_vector() time: " << t.stop() << " s";
        OpmLog::info(oss.str());
    }
}

void OpenclKernels::spmv(cl::Buffer& vals, cl::Buffer& cols, cl::Buffer& rows,
                         const cl::Buffer& x, cl::Buffer& b, int Nb,
                         unsigned int block_size, bool reset, bool add)
{
    const unsigned int work_group_size = 32;
    const unsigned int num_work_groups = ceilDivision(Nb, work_group_size);
    const unsigned int total_work_items = num_work_groups * work_group_size;
    const unsigned int lmem_per_work_group = sizeof(double) * work_group_size;
    Timer t_spmv;
    cl::Event event;

    if (block_size > 1) {
        if (add) {
            event = (*spmv_blocked_add_k)(cl::EnqueueArgs(*queue, cl::NDRange(total_work_items), cl::NDRange(work_group_size)),
                        vals, cols, rows, Nb, x, b, block_size, cl::Local(lmem_per_work_group));
        } else {
            event = (*spmv_blocked_k)(cl::EnqueueArgs(*queue, cl::NDRange(total_work_items), cl::NDRange(work_group_size)),
                        vals, cols, rows, Nb, x, b, block_size, cl::Local(lmem_per_work_group));
        }
    } else {
        if (reset) {
            event = (*spmv_k)(cl::EnqueueArgs(*queue, cl::NDRange(total_work_items), cl::NDRange(work_group_size)),
                        vals, cols, rows, Nb, x, b, cl::Local(lmem_per_work_group));
        } else {
            event = (*spmv_noreset_k)(cl::EnqueueArgs(*queue, cl::NDRange(total_work_items), cl::NDRange(work_group_size)),
                        vals, cols, rows, Nb, x, b, cl::Local(lmem_per_work_group));
        }
    }

    if (verbosity >= 4) {
        event.wait();
        std::ostringstream oss;
        oss << std::scientific << "OpenclKernels spmv_blocked() time: " << t_spmv.stop() << " s";
        OpmLog::info(oss.str());
    }
}

void OpenclKernels::residual(cl::Buffer& vals, cl::Buffer& cols, cl::Buffer& rows,
                            cl::Buffer& x, const cl::Buffer& rhs,
                            cl::Buffer& out, int Nb, unsigned int block_size)
{
    const unsigned int work_group_size = 32;
    const unsigned int num_work_groups = ceilDivision(Nb, work_group_size);
    const unsigned int total_work_items = num_work_groups * work_group_size;
    const unsigned int lmem_per_work_group = sizeof(double) * work_group_size;
    Timer t_residual;
    cl::Event event;

    if (block_size > 1) {
        event = (*residual_blocked_k)(cl::EnqueueArgs(*queue, cl::NDRange(total_work_items), cl::NDRange(work_group_size)),
                    vals, cols, rows, Nb, x, rhs, out, block_size, cl::Local(lmem_per_work_group));
    } else {
        event = (*residual_k)(cl::EnqueueArgs(*queue, cl::NDRange(total_work_items), cl::NDRange(work_group_size)),
                    vals, cols, rows, Nb, x, rhs, out, cl::Local(lmem_per_work_group));
    }

    if (verbosity >= 4) {
        event.wait();
        std::ostringstream oss;
        oss << std::scientific << "OpenclKernels residual_blocked() time: " << t_residual.stop() << " s";
        OpmLog::info(oss.str());
    }
}

void OpenclKernels::ILU_apply1(cl::Buffer& rowIndices, cl::Buffer& vals, cl::Buffer& cols,
                               cl::Buffer& rows, cl::Buffer& diagIndex,
                               const cl::Buffer& y, cl::Buffer& x,
                               cl::Buffer& rowsPerColor, int color,
                               int rowsThisColor, unsigned int block_size)
{
    const unsigned int work_group_size = preferred_workgroup_size_multiple;
    const unsigned int num_work_groups = rowsThisColor;
    const unsigned int total_work_items = num_work_groups * work_group_size;
    const unsigned int lmem_per_work_group = sizeof(double) * work_group_size;
    Timer t_ilu_apply1;

    cl::Event event = (*ILU_apply1_k)(cl::EnqueueArgs(*queue, cl::NDRange(total_work_items), cl::NDRange(work_group_size)),
                          rowIndices, vals, cols, rows, diagIndex,
                          y, x, rowsPerColor, color, block_size,
                          cl::Local(lmem_per_work_group));

    if (verbosity >= 5) {
        event.wait();
        std::ostringstream oss;
        oss << std::scientific << "OpenclKernels ILU_apply1() time: " << t_ilu_apply1.stop() << " s";
        OpmLog::info(oss.str());
    }
}

void OpenclKernels::ILU_apply2(cl::Buffer& rowIndices, cl::Buffer& vals, cl::Buffer& cols,
                               cl::Buffer& rows, cl::Buffer& diagIndex,
                               cl::Buffer& invDiagVals, cl::Buffer& x,
                               cl::Buffer& rowsPerColor, int color,
                               int rowsThisColor, unsigned int block_size)
{
    const unsigned int work_group_size = preferred_workgroup_size_multiple;
    const unsigned int num_work_groups = rowsThisColor;
    const unsigned int total_work_items = num_work_groups * work_group_size;
    const unsigned int lmem_per_work_group = sizeof(double) * work_group_size;
    Timer t_ilu_apply2;

    cl::Event event = (*ILU_apply2_k)(cl::EnqueueArgs(*queue, cl::NDRange(total_work_items), cl::NDRange(work_group_size)),
                          rowIndices, vals, cols, rows, diagIndex,
                          invDiagVals, x, rowsPerColor, color, block_size,
                          cl::Local(lmem_per_work_group));

    if (verbosity >= 5) {
        event.wait();
        std::ostringstream oss;
        oss << std::scientific << "OpenclKernels ILU_apply2() time: " << t_ilu_apply2.stop() << " s";
        OpmLog::info(oss.str());
    }
}

void OpenclKernels::ILU_decomp(int firstRow, int lastRow, cl::Buffer& rowIndices,
                               cl::Buffer& vals, cl::Buffer& cols, cl::Buffer& rows,
                               cl::Buffer& diagIndex, cl::Buffer& invDiagVals,
                               int rowsThisColor, unsigned int block_size)
{
    const unsigned int work_group_size = 128;
    const unsigned int num_work_groups = rowsThisColor;
    const unsigned int total_work_items = num_work_groups * work_group_size;
    const unsigned int num_hwarps_per_group = work_group_size / 16;
    const unsigned int lmem_per_work_group = num_hwarps_per_group * block_size * block_size * sizeof(double);           // each block needs a pivot
    Timer t_ilu_decomp;

    cl::Event event = (*ilu_decomp_k)(cl::EnqueueArgs(*queue, cl::NDRange(total_work_items), cl::NDRange(work_group_size)),
                          firstRow, lastRow, rowIndices,
                          vals, cols, rows,
                          invDiagVals, diagIndex, rowsThisColor,
                          cl::Local(lmem_per_work_group));

    if (verbosity >= 4) {
        event.wait();
        std::ostringstream oss;
        oss << std::scientific << "OpenclKernels ILU_decomp() time: " << t_ilu_decomp.stop() << " s";
        OpmLog::info(oss.str());
    }
}

void OpenclKernels::apply_stdwells(cl::Buffer& d_Cnnzs_ocl, cl::Buffer &d_Dnnzs_ocl, cl::Buffer &d_Bnnzs_ocl,
    cl::Buffer &d_Ccols_ocl, cl::Buffer &d_Bcols_ocl, cl::Buffer &d_x, cl::Buffer &d_y,
    int dim, int dim_wells, cl::Buffer &d_val_pointers_ocl, int num_std_wells)
{
    const unsigned int work_group_size = 32;
    const unsigned int total_work_items = num_std_wells * work_group_size;
    const unsigned int lmem1 = sizeof(double) * work_group_size;
    const unsigned int lmem2 = sizeof(double) * dim_wells;
    Timer t_apply_stdwells;

    cl::Event event = (*stdwell_apply_k)(cl::EnqueueArgs(*queue, cl::NDRange(total_work_items), cl::NDRange(work_group_size)),
                          d_Cnnzs_ocl, d_Dnnzs_ocl, d_Bnnzs_ocl, d_Ccols_ocl, d_Bcols_ocl, d_x, d_y, dim, dim_wells, d_val_pointers_ocl,
                          cl::Local(lmem1), cl::Local(lmem2), cl::Local(lmem2));

    if (verbosity >= 4) {
        event.wait();
        std::ostringstream oss;
        oss << std::scientific << "OpenclKernels apply_stdwells() time: " << t_apply_stdwells.stop() << " s";
        OpmLog::info(oss.str());
    }
}

void OpenclKernels::isaiL(cl::Buffer& diagIndex, cl::Buffer& colPointers, cl::Buffer& mapping, cl::Buffer& nvc,
    cl::Buffer& luIdxs, cl::Buffer& xxIdxs, cl::Buffer& dxIdxs, cl::Buffer& LUvals, cl::Buffer& invLvals, unsigned int Nb)
{
    const unsigned int work_group_size = 256;
    const unsigned int num_work_groups = ceilDivision(Nb, work_group_size);
    const unsigned int total_work_items = num_work_groups * work_group_size;

    Timer t_isaiL;
    cl::Event event = (*isaiL_k)(cl::EnqueueArgs(*queue, cl::NDRange(total_work_items), cl::NDRange(work_group_size)),
            diagIndex, colPointers, mapping, nvc, luIdxs, xxIdxs, dxIdxs, LUvals, invLvals, Nb);

    if (verbosity >= 4) {
        event.wait();
        std::ostringstream oss;
        oss << std::scientific << "OpenclKernels isaiL() time: " << t_isaiL.stop() << " s";
        OpmLog::info(oss.str());
    }
}

void OpenclKernels::isaiU(cl::Buffer& diagIndex, cl::Buffer& colPointers, cl::Buffer& rowIndices, cl::Buffer& mapping,
        cl::Buffer& nvc, cl::Buffer& luIdxs, cl::Buffer& xxIdxs, cl::Buffer& dxIdxs, cl::Buffer& LUvals,
        cl::Buffer& invDiagVals, cl::Buffer& invUvals, unsigned int Nb)
{
    const unsigned int work_group_size = 256;
    const unsigned int num_work_groups = ceilDivision(Nb, work_group_size);
    const unsigned int total_work_items = num_work_groups * work_group_size;

    Timer t_isaiU;
    cl::Event event = (*isaiU_k)(cl::EnqueueArgs(*queue, cl::NDRange(total_work_items), cl::NDRange(work_group_size)),
            diagIndex, colPointers, rowIndices, mapping, nvc, luIdxs, xxIdxs, dxIdxs, LUvals, invDiagVals, invUvals, Nb);

    if (verbosity >= 4) {
        event.wait();
        std::ostringstream oss;
        oss << std::scientific << "OpenclKernels isaiU() time: " << t_isaiU.stop() << " s";
        OpmLog::info(oss.str());
    }
}

} // namespace Accelerator
} // namespace Opm
