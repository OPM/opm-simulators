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
#include <iostream>

#include <opm/common/OpmLog/OpmLog.hpp>
#include <opm/common/ErrorMacros.hpp>
#include <dune/common/timer.hh>

#include <opm/simulators/linalg/bda/openclSolverBackend.hpp>
#include <opm/simulators/linalg/bda/openclKernels.hpp>

#include <opm/simulators/linalg/bda/BdaResult.hpp>
#include <opm/simulators/linalg/bda/Reorder.hpp>


// iff true, the nonzeroes of the matrix are copied row-by-row into a contiguous, pinned memory array, then a single GPU memcpy is done
// otherwise, the nonzeroes of the matrix are assumed to be in a contiguous array, and a single GPU memcpy is enough
#define COPY_ROW_BY_ROW 0

// Level Scheduling respects the depencies in the original matrix
// Graph Coloring is more aggresive and is likely to change the number of linearizations and linear iterations to converge, but can still be faster on GPU because it results in more parallelism
#define LEVEL_SCHEDULING 0
#define GRAPH_COLORING   1

namespace bda
{

using Opm::OpmLog;
using Dune::Timer;

template <unsigned int block_size>
openclSolverBackend<block_size>::openclSolverBackend(int verbosity_, int maxit_, double tolerance_, unsigned int platformID_, unsigned int deviceID_) : BdaSolver<block_size>(verbosity_, maxit_, tolerance_, platformID_, deviceID_) {
    prec = new Preconditioner(LEVEL_SCHEDULING, GRAPH_COLORING, verbosity_);
}


template <unsigned int block_size>
openclSolverBackend<block_size>::~openclSolverBackend() {
    finalize();
}


// divide A by B, and round up: return (int)ceil(A/B)
template <unsigned int block_size>
unsigned int openclSolverBackend<block_size>::ceilDivision(const unsigned int A, const unsigned int B)
{
    return A / B + (A % B > 0);
}


template <unsigned int block_size>
double openclSolverBackend<block_size>::dot_w(cl::Buffer in1, cl::Buffer in2, cl::Buffer out)
{
    const unsigned int work_group_size = 256;
    const unsigned int num_work_groups = ceilDivision(N, work_group_size);
    const unsigned int total_work_items = num_work_groups * work_group_size;
    const unsigned int lmem_per_work_group = sizeof(double) * work_group_size;
    Timer t_dot;

    cl::Event event = (*dot_k)(cl::EnqueueArgs(*queue, cl::NDRange(total_work_items), cl::NDRange(work_group_size)), in1, in2, out, N, cl::Local(lmem_per_work_group));

    queue->enqueueReadBuffer(out, CL_TRUE, 0, sizeof(double) * num_work_groups, tmp);

    double gpu_sum = 0.0;
    for (unsigned int i = 0; i < num_work_groups; ++i) {
        gpu_sum += tmp[i];
    }

    if (verbosity >= 4) {
        event.wait();
        std::ostringstream oss;
        oss << std::scientific << "openclSolver dot_w time: " << t_dot.stop() << " s";
        OpmLog::info(oss.str());
    }

    return gpu_sum;
}

template <unsigned int block_size>
double openclSolverBackend<block_size>::norm_w(cl::Buffer in, cl::Buffer out)
{
    const unsigned int work_group_size = 256;
    const unsigned int num_work_groups = ceilDivision(N, work_group_size);
    const unsigned int total_work_items = num_work_groups * work_group_size;
    const unsigned int lmem_per_work_group = sizeof(double) * work_group_size;
    Timer t_norm;

    cl::Event event = (*norm_k)(cl::EnqueueArgs(*queue, cl::NDRange(total_work_items), cl::NDRange(work_group_size)), in, out, N, cl::Local(lmem_per_work_group));

    queue->enqueueReadBuffer(out, CL_TRUE, 0, sizeof(double) * num_work_groups, tmp);

    double gpu_norm = 0.0;
    for (unsigned int i = 0; i < num_work_groups; ++i) {
        gpu_norm += tmp[i];
    }
    gpu_norm = sqrt(gpu_norm);

    if (verbosity >= 4) {
        event.wait();
        std::ostringstream oss;
        oss << std::scientific << "openclSolver norm_w time: " << t_norm.stop() << " s";
        OpmLog::info(oss.str());
    }

    return gpu_norm;
}

template <unsigned int block_size>
void openclSolverBackend<block_size>::axpy_w(cl::Buffer in, const double a, cl::Buffer out)
{
    const unsigned int work_group_size = 32;
    const unsigned int num_work_groups = ceilDivision(N, work_group_size);
    const unsigned int total_work_items = num_work_groups * work_group_size;
    Timer t_axpy;

    cl::Event event = (*axpy_k)(cl::EnqueueArgs(*queue, cl::NDRange(total_work_items), cl::NDRange(work_group_size)), in, a, out, N);

    if (verbosity >= 4) {
        event.wait();
        std::ostringstream oss;
        oss << std::scientific << "openclSolver axpy_w time: " << t_axpy.stop() << " s";
        OpmLog::info(oss.str());
    }
}

template <unsigned int block_size>
void openclSolverBackend<block_size>::custom_w(cl::Buffer p, cl::Buffer v, cl::Buffer r, const double omega, const double beta)
{
    const unsigned int work_group_size = 32;
    const unsigned int num_work_groups = ceilDivision(N, work_group_size);
    const unsigned int total_work_items = num_work_groups * work_group_size;
    Timer t_custom;

    cl::Event event = (*custom_k)(cl::EnqueueArgs(*queue, cl::NDRange(total_work_items), cl::NDRange(work_group_size)), p, v, r, omega, beta, N);

    if (verbosity >= 4) {
        event.wait();
        std::ostringstream oss;
        oss << std::scientific << "openclSolver custom_w time: " << t_custom.stop() << " s";
        OpmLog::info(oss.str());
    }
}

template <unsigned int block_size>
void openclSolverBackend<block_size>::spmv_blocked_w(cl::Buffer vals, cl::Buffer cols, cl::Buffer rows, cl::Buffer x, cl::Buffer b)
{
    const unsigned int work_group_size = 32;
    const unsigned int num_work_groups = ceilDivision(N, work_group_size);
    const unsigned int total_work_items = num_work_groups * work_group_size;
    const unsigned int lmem_per_work_group = sizeof(double) * work_group_size;
    Timer t_spmv;

    cl::Event event = (*spmv_blocked_k)(cl::EnqueueArgs(*queue, cl::NDRange(total_work_items), cl::NDRange(work_group_size)), vals, cols, rows, Nb, x, b, block_size, cl::Local(lmem_per_work_group));

    if (verbosity >= 4) {
        event.wait();
        std::ostringstream oss;
        oss << std::scientific << "openclSolver spmv_blocked_w time: " << t_spmv.stop() << " s";
        OpmLog::info(oss.str());
    }
}

template <unsigned int block_size>
void openclSolverBackend<block_size>::gpu_pbicgstab(WellContributions& wellContribs, BdaResult& res) {
    float it;
    double rho, rhop, beta, alpha, omega, tmp1, tmp2;
    double norm, norm_0;

    Timer t_total, t_prec(false), t_spmv(false), t_well(false), t_rest(false);

    // set r to the initial residual
    // if initial x guess is not 0, must call applyblockedscaleadd(), not implemented
    //applyblockedscaleadd(-1.0, mat, x, r);

    // set initial values
    cl::Event event;
    queue->enqueueFillBuffer(d_p, 0, 0, sizeof(double) * N);
    queue->enqueueFillBuffer(d_v, 0, 0, sizeof(double) * N, nullptr, &event);
    event.wait();
    rho = 1.0;
    alpha = 1.0;
    omega = 1.0;

    queue->enqueueCopyBuffer(d_b, d_r, 0, 0, sizeof(double) * N, nullptr, &event);
    event.wait();
    queue->enqueueCopyBuffer(d_r, d_rw, 0, 0, sizeof(double) * N, nullptr, &event);
    event.wait();
    queue->enqueueCopyBuffer(d_r, d_p, 0, 0, sizeof(double) * N, nullptr, &event);
    event.wait();

    wellContribs.setReordering(toOrder, true);

    norm = norm_w(d_r, d_tmp);
    norm_0 = norm;

    if (verbosity > 1) {
        std::ostringstream out;
        out << std::scientific << "openclSolver initial norm: " << norm_0;
        OpmLog::info(out.str());
    }

    t_rest.start();
    for (it = 0.5; it < maxit; it += 0.5) {
        rhop = rho;
        rho = dot_w(d_rw, d_r, d_tmp);

        if (it > 1) {
            beta = (rho / rhop) * (alpha / omega);
            custom_w(d_p, d_v, d_r, omega, beta);
        }
        t_rest.stop();

        // pw = prec(p)
        t_prec.start();
        prec->apply(d_p, d_pw);
        t_prec.stop();

        // v = A * pw
        t_spmv.start();
        spmv_blocked_w(d_Avals, d_Acols, d_Arows, d_pw, d_v);
        t_spmv.stop();

        // apply wellContributions
        t_well.start();
        wellContribs.apply(d_pw, d_v);
        t_well.stop();

        t_rest.start();
        tmp1 = dot_w(d_rw, d_v, d_tmp);
        alpha = rho / tmp1;
        axpy_w(d_v, -alpha, d_r);      // r = r - alpha * v
        axpy_w(d_pw, alpha, d_x);      // x = x + alpha * pw
        norm = norm_w(d_r, d_tmp);
        t_rest.stop();

        if (norm < tolerance * norm_0) {
            break;
        }

        it += 0.5;

        // s = prec(r)
        t_prec.start();
        prec->apply(d_r, d_s);
        t_prec.stop();

        // t = A * s
        t_spmv.start();
        spmv_blocked_w(d_Avals, d_Acols, d_Arows, d_s, d_t);
        t_spmv.stop();

        // apply wellContributions
        t_well.start();
        wellContribs.apply(d_s, d_t);
        t_well.stop();

        t_rest.start();
        tmp1 = dot_w(d_t, d_r, d_tmp);
        tmp2 = dot_w(d_t, d_t, d_tmp);
        omega = tmp1 / tmp2;
        axpy_w(d_s, omega, d_x);     // x = x + omega * s
        axpy_w(d_t, -omega, d_r);    // r = r - omega * t
        norm = norm_w(d_r, d_tmp);
        t_rest.stop();

        if (norm < tolerance * norm_0) {
            break;
        }

        if (verbosity > 1) {
            std::ostringstream out;
            out << "it: " << it << std::scientific << ", norm: " << norm;
            OpmLog::info(out.str());
        }
    }

    res.iterations = std::min(it, (float)maxit);
    res.reduction = norm / norm_0;
    res.conv_rate  = static_cast<double>(pow(res.reduction, 1.0 / it));
    res.elapsed = t_total.stop();
    res.converged = (it != (maxit + 0.5));

    if (verbosity > 0) {
        std::ostringstream out;
        out << "=== converged: " << res.converged << ", conv_rate: " << res.conv_rate << ", time: " << res.elapsed << \
            ", time per iteration: " << res.elapsed / it << ", iterations: " << it;
        OpmLog::info(out.str());
    }
    if (verbosity >= 4) {
        std::ostringstream out;
        out << "openclSolver::ilu_apply:   " << t_prec.elapsed() << " s\n";
        out << "wellContributions::apply:  " << t_well.elapsed() << " s\n";
        out << "openclSolver::spmv:        " << t_spmv.elapsed() << " s\n";
        out << "openclSolver::rest:        " << t_rest.elapsed() << " s\n";
        out << "openclSolver::total_solve: " << res.elapsed << " s\n";
        OpmLog::info(out.str());
    }
}


template <unsigned int block_size>
void openclSolverBackend<block_size>::initialize(int N_, int nnz_, int dim, double *vals, int *rows, int *cols) {
    this->N = N_;
    this->nnz = nnz_;
    this->nnzb = nnz_ / block_size / block_size;

    Nb = (N + dim - 1) / dim;
    std::ostringstream out;
    out << "Initializing GPU, matrix size: " << N << " blocks, nnzb: " << nnzb << "\n";
    out << "Maxit: " << maxit << std::scientific << ", tolerance: " << tolerance << "\n";
    out << "PlatformID: " << platformID << ", deviceID: " << deviceID << "\n";
    OpmLog::info(out.str());
    out.str("");
    out.clear();

    cl_int err = CL_SUCCESS;
    try {
        std::vector<cl::Platform> platforms;
        cl::Platform::get(&platforms);
        if (platforms.size() == 0) {
            OPM_THROW(std::logic_error, "Error openclSolver is selected but no OpenCL platforms are found");
        }
        out << "Found " << platforms.size() << " OpenCL platforms" << "\n";

        if (verbosity >= 1) {
            std::string platform_info;
            for (unsigned int i = 0; i < platforms.size(); ++i) {
                platforms[i].getInfo(CL_PLATFORM_NAME, &platform_info);
                out << "Platform name      : " << platform_info << "\n";
                platforms[i].getInfo(CL_PLATFORM_VENDOR, &platform_info);
                out << "Platform vendor    : " << platform_info << "\n";
                platforms[i].getInfo(CL_PLATFORM_VERSION, &platform_info);
                out << "Platform version   : " << platform_info << "\n";
                platforms[i].getInfo(CL_PLATFORM_PROFILE, &platform_info);
                out << "Platform profile   : " << platform_info << "\n";
                platforms[i].getInfo(CL_PLATFORM_EXTENSIONS, &platform_info);
                out << "Platform extensions: " << platform_info << "\n\n";
            }
        }
        OpmLog::info(out.str());
        out.str("");
        out.clear();

        if (platforms.size() <= platformID) {
            OPM_THROW(std::logic_error, "Error chosen too high OpenCL platform ID");
        } else {
            std::string platform_info;
            out << "Chosen:\n";
            platforms[platformID].getInfo(CL_PLATFORM_NAME, &platform_info);
            out << "Platform name      : " << platform_info << "\n";
            platforms[platformID].getInfo(CL_PLATFORM_VERSION, &platform_info);
            out << "Platform version   : " << platform_info << "\n";
            OpmLog::info(out.str());
            out.str("");
            out.clear();
        }

        cl_context_properties properties[] = {CL_CONTEXT_PLATFORM, (cl_context_properties)(platforms[platformID])(), 0};
        context.reset(new cl::Context(CL_DEVICE_TYPE_GPU, properties));

        std::vector<cl::Device> devices = context->getInfo<CL_CONTEXT_DEVICES>();
        if (devices.size() == 0){
            OPM_THROW(std::logic_error, "Error openclSolver is selected but no OpenCL devices are found");
        }
        out << "Found " << devices.size() << " OpenCL devices" << "\n";

        if (verbosity >= 1) {
            for (unsigned int i = 0; i < devices.size(); ++i) {
                std::string device_info;
                std::vector<size_t> work_sizes;
                std::vector<cl_device_partition_property> partitions;

                devices[i].getInfo(CL_DEVICE_NAME, &device_info);
                out << "CL_DEVICE_NAME            : " << device_info << "\n";
                devices[i].getInfo(CL_DEVICE_VENDOR, &device_info);
                out << "CL_DEVICE_VENDOR          : " << device_info << "\n";
                devices[i].getInfo(CL_DRIVER_VERSION, &device_info);
                out << "CL_DRIVER_VERSION         : " << device_info << "\n";
                devices[i].getInfo(CL_DEVICE_BUILT_IN_KERNELS, &device_info);
                out << "CL_DEVICE_BUILT_IN_KERNELS: " << device_info << "\n";
                devices[i].getInfo(CL_DEVICE_PROFILE, &device_info);
                out << "CL_DEVICE_PROFILE         : " << device_info << "\n";
                devices[i].getInfo(CL_DEVICE_OPENCL_C_VERSION, &device_info);
                out << "CL_DEVICE_OPENCL_C_VERSION: " << device_info << "\n";
                devices[i].getInfo(CL_DEVICE_EXTENSIONS, &device_info);
                out << "CL_DEVICE_EXTENSIONS      : " << device_info << "\n";

                devices[i].getInfo(CL_DEVICE_MAX_WORK_ITEM_SIZES, &work_sizes);
                for (unsigned int j = 0; j < work_sizes.size(); ++j) {
                    out << "CL_DEVICE_MAX_WORK_ITEM_SIZES[" << j << "]: " << work_sizes[j] << "\n";
                }
                devices[i].getInfo(CL_DEVICE_PARTITION_PROPERTIES, &partitions);
                for (unsigned int j = 0; j < partitions.size(); ++j) {
                    out << "CL_DEVICE_PARTITION_PROPERTIES[" << j << "]: " << partitions[j] << "\n";
                }
                partitions.clear();
                devices[i].getInfo(CL_DEVICE_PARTITION_TYPE, &partitions);
                for (unsigned int j = 0; j < partitions.size(); ++j) {
                    out << "CL_DEVICE_PARTITION_PROPERTIES[" << j << "]: " << partitions[j] << "\n";
                }

                // C-style properties
                cl_device_id tmp_id = devices[i]();
                cl_ulong size;
                clGetDeviceInfo(tmp_id, CL_DEVICE_LOCAL_MEM_SIZE, sizeof(cl_ulong), &size, 0);
                out << "CL_DEVICE_LOCAL_MEM_SIZE       : " << size / 1024 << " KB\n";
                clGetDeviceInfo(tmp_id, CL_DEVICE_GLOBAL_MEM_SIZE, sizeof(cl_ulong), &size, 0);
                out << "CL_DEVICE_GLOBAL_MEM_SIZE      : " << size / 1024 / 1024 / 1024 << " GB\n";
                clGetDeviceInfo(tmp_id, CL_DEVICE_MAX_COMPUTE_UNITS, sizeof(cl_ulong), &size, 0);
                out << "CL_DEVICE_MAX_COMPUTE_UNITS    : " << size << "\n";
                clGetDeviceInfo(tmp_id, CL_DEVICE_MAX_MEM_ALLOC_SIZE, sizeof(cl_ulong), &size, 0);
                out << "CL_DEVICE_MAX_MEM_ALLOC_SIZE   : " << size / 1024 / 1024 << " MB\n";
                clGetDeviceInfo(tmp_id, CL_DEVICE_MAX_WORK_GROUP_SIZE, sizeof(cl_ulong), &size, 0);
                out << "CL_DEVICE_MAX_WORK_GROUP_SIZE  : " << size << "\n";
                clGetDeviceInfo(tmp_id, CL_DEVICE_GLOBAL_MEM_SIZE, sizeof(cl_ulong), &size, 0);
                out << "CL_DEVICE_GLOBAL_MEM_SIZE      : " << size / 1024 / 1024 / 1024 << " GB\n\n";
            }
        }
        OpmLog::info(out.str());
        out.str("");
        out.clear();

        if (devices.size() <= deviceID){
            OPM_THROW(std::logic_error, "Error chosen too high OpenCL device ID");
        } else {
            std::string device_info;
            out << "Chosen:\n";
            devices[deviceID].getInfo(CL_DEVICE_NAME, &device_info);
            out << "CL_DEVICE_NAME            : " << device_info << "\n";
            devices[deviceID].getInfo(CL_DEVICE_VERSION, &device_info);
            out << "CL_DEVICE_VERSION         : " << device_info << "\n";
            OpmLog::info(out.str());
            out.str("");
            out.clear();
        }

        cl::Program::Sources source(1, std::make_pair(axpy_s, strlen(axpy_s)));  // what does this '1' mean? cl::Program::Sources is of type 'std::vector<std::pair<const char*, long unsigned int> >'
        source.emplace_back(std::make_pair(dot_1_s, strlen(dot_1_s)));
        source.emplace_back(std::make_pair(norm_s, strlen(norm_s)));
        source.emplace_back(std::make_pair(custom_s, strlen(custom_s)));
        source.emplace_back(std::make_pair(spmv_blocked_s, strlen(spmv_blocked_s)));
        source.emplace_back(std::make_pair(ILU_apply1_s, strlen(ILU_apply1_s)));
        source.emplace_back(std::make_pair(ILU_apply2_s, strlen(ILU_apply2_s)));
        source.emplace_back(std::make_pair(add_well_contributions_s, strlen(add_well_contributions_s)));
        program = cl::Program(*context, source);

        program.build(devices);

        cl::Event event;
        queue.reset(new cl::CommandQueue(*context, devices[deviceID], 0, &err));

        prec->setOpenCLContext(context.get());
        prec->setOpenCLQueue(queue.get());

        rb = new double[N];
        tmp = new double[N];
#if COPY_ROW_BY_ROW
        vals_contiguous = new double[N];
#endif
        mat.reset(new BlockedMatrix<block_size>(Nb, nnzb, vals, cols, rows));

        d_x = cl::Buffer(*context, CL_MEM_READ_WRITE, sizeof(double) * N);
        d_b = cl::Buffer(*context, CL_MEM_READ_WRITE, sizeof(double) * N);
        d_rb = cl::Buffer(*context, CL_MEM_READ_WRITE, sizeof(double) * N);
        d_r = cl::Buffer(*context, CL_MEM_READ_WRITE, sizeof(double) * N);
        d_rw = cl::Buffer(*context, CL_MEM_READ_WRITE, sizeof(double) * N);
        d_p = cl::Buffer(*context, CL_MEM_READ_WRITE, sizeof(double) * N);
        d_pw = cl::Buffer(*context, CL_MEM_READ_WRITE, sizeof(double) * N);
        d_s = cl::Buffer(*context, CL_MEM_READ_WRITE, sizeof(double) * N);
        d_t = cl::Buffer(*context, CL_MEM_READ_WRITE, sizeof(double) * N);
        d_v = cl::Buffer(*context, CL_MEM_READ_WRITE, sizeof(double) * N);
        d_tmp = cl::Buffer(*context, CL_MEM_READ_WRITE, sizeof(double) * N);

        d_Avals = cl::Buffer(*context, CL_MEM_READ_WRITE, sizeof(double) * nnz);
        d_Acols = cl::Buffer(*context, CL_MEM_READ_WRITE, sizeof(int) * nnzb);
        d_Arows = cl::Buffer(*context, CL_MEM_READ_WRITE, sizeof(int) * (Nb + 1));

        // queue.enqueueNDRangeKernel() is a blocking/synchronous call, at least for NVIDIA
        // cl::make_kernel<> myKernel(); myKernel(args, arg1, arg2); is also blocking

        // actually creating the kernels
        dot_k.reset(new cl::make_kernel<cl::Buffer&, cl::Buffer&, cl::Buffer&, const unsigned int, cl::LocalSpaceArg>(cl::Kernel(program, "dot_1")));
        norm_k.reset(new cl::make_kernel<cl::Buffer&, cl::Buffer&, const unsigned int, cl::LocalSpaceArg>(cl::Kernel(program, "norm")));
        axpy_k.reset(new cl::make_kernel<cl::Buffer&, const double, cl::Buffer&, const unsigned int>(cl::Kernel(program, "axpy")));
        custom_k.reset(new cl::make_kernel<cl::Buffer&, cl::Buffer&, cl::Buffer&, const double, const double, const unsigned int>(cl::Kernel(program, "custom")));
        spmv_blocked_k.reset(new cl::make_kernel<cl::Buffer&, cl::Buffer&, cl::Buffer&, const unsigned int, cl::Buffer&, cl::Buffer&, const unsigned int, cl::LocalSpaceArg>(cl::Kernel(program, "spmv_blocked")));
        ILU_apply1_k.reset(new cl::make_kernel<cl::Buffer&, cl::Buffer&, cl::Buffer&, const unsigned int, cl::Buffer&, cl::Buffer&, cl::Buffer&, const unsigned int, const unsigned int, cl::LocalSpaceArg>(cl::Kernel(program, "ILU_apply1")));
        ILU_apply2_k.reset(new cl::make_kernel<cl::Buffer&, cl::Buffer&, cl::Buffer&, const unsigned int, cl::Buffer&, cl::Buffer&, cl::Buffer&, const unsigned int, const unsigned int, cl::LocalSpaceArg>(cl::Kernel(program, "ILU_apply2")));

        prec->setKernels(ILU_apply1_k.get(), ILU_apply2_k.get());

    } catch (const cl::Error& error) {
        std::ostringstream oss;
        oss << "OpenCL Error: " << error.what() << "(" << error.err() << ")\n";
        oss << getErrorString(error.err());
        // rethrow exception
        OPM_THROW(std::logic_error, oss.str());
    } catch (const std::logic_error& error) {
        // rethrow exception by OPM_THROW in the try{}, without this, a segfault occurs
        throw error;
    }

    initialized = true;
} // end initialize()

template <unsigned int block_size>
void openclSolverBackend<block_size>::initialize_wellContribs(WellContributions& wellContribs){
    add_well_contributions_k.reset(new cl::make_kernel<cl::Buffer&, cl::Buffer&, cl::Buffer&, cl::Buffer&, cl::Buffer&, cl::Buffer&, cl::Buffer&, const unsigned int, const unsigned int, cl::Buffer&, cl::LocalSpaceArg, cl::LocalSpaceArg, cl::LocalSpaceArg>(cl::Kernel(program, "add_well_contributions")));

    wellContribs.setOpenCLContext(context.get());
    wellContribs.setOpenCLQueue(queue.get());
    wellContribs.init();
    wellContribs.setKernel(add_well_contributions_k.get());
}

template <unsigned int block_size>
void openclSolverBackend<block_size>::finalize() {
    delete[] rb;
    delete[] tmp;
#if COPY_ROW_BY_ROW
    delete[] vals_contiguous;
#endif
    delete prec;
} // end finalize()


template <unsigned int block_size>
void openclSolverBackend<block_size>::copy_system_to_gpu() {
    Timer t;
    cl::Event event;

#if COPY_ROW_BY_ROW
    int sum = 0;
    for (int i = 0; i < Nb; ++i) {
        int size_row = rmat->rowPointers[i + 1] - rmat->rowPointers[i];
        memcpy(vals_contiguous + sum, reinterpret_cast<double*>(rmat->nnzValues) + sum, size_row * sizeof(double) * block_size * block_size);
        sum += size_row * block_size * block_size;
    }
    queue->enqueueWriteBuffer(d_Avals, CL_TRUE, 0, sizeof(double) * nnz, vals_contiguous);
#else
    queue->enqueueWriteBuffer(d_Avals, CL_TRUE, 0, sizeof(double) * nnz, rmat->nnzValues);
#endif

    queue->enqueueWriteBuffer(d_Acols, CL_TRUE, 0, sizeof(int) * nnzb, rmat->colIndices);
    queue->enqueueWriteBuffer(d_Arows, CL_TRUE, 0, sizeof(int) * (Nb + 1), rmat->rowPointers);
    queue->enqueueWriteBuffer(d_b, CL_TRUE, 0, sizeof(double) * N, rb);
    queue->enqueueFillBuffer(d_x, 0, 0, sizeof(double) * N, nullptr, &event);
    event.wait();

    if (verbosity > 2) {
        std::ostringstream out;
        out << "openclSolver::copy_system_to_gpu(): " << t.stop() << " s";
        OpmLog::info(out.str());
    }
} // end copy_system_to_gpu()


// don't copy rowpointers and colindices, they stay the same
template <unsigned int block_size>
void openclSolverBackend<block_size>::update_system_on_gpu() {
    Timer t;
    cl::Event event;

#if COPY_ROW_BY_ROW
    int sum = 0;
    for (int i = 0; i < Nb; ++i) {
        int size_row = rmat->rowPointers[i + 1] - rmat->rowPointers[i];
        memcpy(vals_contiguous + sum, reinterpret_cast<double*>(rmat->nnzValues) + sum, size_row * sizeof(double) * block_size * block_size);
        sum += size_row * block_size * block_size;
    }
    queue->enqueueWriteBuffer(d_Avals, CL_TRUE, 0, sizeof(double) * nnz, vals_contiguous);
#else
    queue->enqueueWriteBuffer(d_Avals, CL_TRUE, 0, sizeof(double) * nnz, rmat->nnzValues);
#endif

    queue->enqueueWriteBuffer(d_b, CL_TRUE, 0, sizeof(double) * N, rb);
    queue->enqueueFillBuffer(d_x, 0, 0, sizeof(double) * N, nullptr, &event);
    event.wait();

    if (verbosity > 2) {
        std::ostringstream out;
        out << "openclSolver::update_system_on_gpu(): " << t.stop() << " s";
        OpmLog::info(out.str());
    }
} // end update_system_on_gpu()


template <unsigned int block_size>
bool openclSolverBackend<block_size>::analyse_matrix() {
    Timer t;

    bool success = prec->init(mat.get());
    int work_group_size = 32;
    int num_work_groups = ceilDivision(N, work_group_size);
    int total_work_items = num_work_groups * work_group_size;
    int lmem_per_work_group = work_group_size * sizeof(double);
    prec->setKernelParameters(work_group_size, total_work_items, lmem_per_work_group);

    toOrder = prec->getToOrder();
    fromOrder = prec->getFromOrder();
    rmat = prec->getRMat();

    if (verbosity > 2) {
        std::ostringstream out;
        out << "openclSolver::analyse_matrix(): " << t.stop() << " s";
        OpmLog::info(out.str());
    }

    analysis_done = true;

    return success;
} // end analyse_matrix()


template <unsigned int block_size>
void openclSolverBackend<block_size>::update_system(double *vals, double *b) {
    Timer t;

    mat->nnzValues = vals;
    reorderBlockedVectorByPattern<block_size>(mat->Nb, b, fromOrder, rb);

    if (verbosity > 2) {
        std::ostringstream out;
        out << "openclSolver::update_system(): " << t.stop() << " s";
        OpmLog::info(out.str());
    }
} // end update_system()


template <unsigned int block_size>
bool openclSolverBackend<block_size>::create_preconditioner() {
    Timer t;

    bool result = prec->create_preconditioner(mat.get());

    if (verbosity > 2) {
        std::ostringstream out;
        out << "openclSolver::create_preconditioner(): " << t.stop() << " s";
        OpmLog::info(out.str());
    }
    return result;
} // end create_preconditioner()


template <unsigned int block_size>
void openclSolverBackend<block_size>::solve_system(WellContributions& wellContribs, BdaResult &res) {
    Timer t;

    // actually solve
    try {
        gpu_pbicgstab(wellContribs, res);
    } catch (const cl::Error& error) {
        std::ostringstream oss;
        oss << "openclSolverBackend::solve_system error: " << error.what() << "(" << error.err() << ")\n";
        oss << getErrorString(error.err());
        // rethrow exception
        OPM_THROW(std::logic_error, oss.str());
    } catch (const std::logic_error& error) {
        // rethrow exception by OPM_THROW in the try{}, without this, a segfault occurs
        throw error;
    }

    if (verbosity > 2) {
        std::ostringstream out;
        out << "openclSolver::solve_system(): " << t.stop() << " s";
        OpmLog::info(out.str());
    }

} // end solve_system()


// copy result to host memory
// caller must be sure that x is a valid array
template <unsigned int block_size>
void openclSolverBackend<block_size>::get_result(double *x) {
    Timer t;

    queue->enqueueReadBuffer(d_x, CL_TRUE, 0, sizeof(double) * N, rb);
    reorderBlockedVectorByPattern<block_size>(mat->Nb, rb, toOrder, x);

    if (verbosity > 2) {
        std::ostringstream out;
        out << "openclSolver::get_result(): " << t.stop() << " s";
        OpmLog::info(out.str());
    }
} // end get_result()


template <unsigned int block_size>
SolverStatus openclSolverBackend<block_size>::solve_system(int N_, int nnz_, int dim, double *vals, int *rows, int *cols, double *b, WellContributions& wellContribs, BdaResult &res) {
    if (initialized == false) {
        initialize(N_, nnz_,  dim, vals, rows, cols);
        if (analysis_done == false) {
            if (!analyse_matrix()) {
                return SolverStatus::BDA_SOLVER_ANALYSIS_FAILED;
            }
        }
        update_system(vals, b);
        if (!create_preconditioner()) {
            return SolverStatus::BDA_SOLVER_CREATE_PRECONDITIONER_FAILED;
        }
        copy_system_to_gpu();
    } else {
        update_system(vals, b);
        if (!create_preconditioner()) {
            return SolverStatus::BDA_SOLVER_CREATE_PRECONDITIONER_FAILED;
        }
        update_system_on_gpu();
    }
    initialize_wellContribs(wellContribs);
    solve_system(wellContribs, res);
    return SolverStatus::BDA_SOLVER_SUCCESS;
}


#define INSTANTIATE_BDA_FUNCTIONS(n)                                                                 \
template openclSolverBackend<n>::openclSolverBackend(int, int, double, unsigned int, unsigned int);  \

INSTANTIATE_BDA_FUNCTIONS(1);
INSTANTIATE_BDA_FUNCTIONS(2);
INSTANTIATE_BDA_FUNCTIONS(3);
INSTANTIATE_BDA_FUNCTIONS(4);

#undef INSTANTIATE_BDA_FUNCTIONS

} // namespace bda

