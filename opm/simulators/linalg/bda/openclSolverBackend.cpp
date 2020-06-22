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

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <iostream>
#include <sstream>

#include <config.h>  // CMake
#include <opm/common/OpmLog/OpmLog.hpp>
#include <opm/common/ErrorMacros.hpp>


#define __CL_ENABLE_EXCEPTIONS
#include <CL/cl.hpp>
#include <opm/simulators/linalg/bda/openclKernels.hpp>


#include <opm/simulators/linalg/bda/openclSolverBackend.hpp>
#include <opm/simulators/linalg/bda/BdaResult.hpp>
#include <opm/simulators/linalg/bda/Reorder.hpp>


// iff true, the nonzeroes of the matrix are copied row-by-row into a contiguous, pinned memory array, then a single GPU memcpy is done
// otherwise, the nonzeroes of the matrix are assumed to be in a contiguous array, and a single GPU memcpy is enough
#define COPY_ROW_BY_ROW 0


#define LEVEL_SCHEDULING 1
#define GRAPH_COLORING   0

namespace bda
{

using Opm::OpmLog;

openclSolverBackend::openclSolverBackend(int verbosity_, int maxit_, double tolerance_) : BdaSolver(verbosity_, maxit_, tolerance_) {
    prec = new Preconditioner(LEVEL_SCHEDULING, GRAPH_COLORING, verbosity_);
}


openclSolverBackend::~openclSolverBackend() {
    finalize();
}


// divide A by B, and round up: return (int)ceil(A/B)
unsigned int openclSolverBackend::ceilDivision(const unsigned int A, const unsigned int B)
{
    return A / B + (A % B > 0);
}

// just for verifying and debugging
bool equal(float a, float b)
{
    const float tol_abs = 1e-2;
    const float tol_rel = 1e-2;
    return std::abs(a - b) <= std::max(tol_rel * std::max(std::abs(a), std::abs(b)), tol_abs);
}

double openclSolverBackend::dot_w(cl::Buffer in1, cl::Buffer in2, cl::Buffer out)
{
    double t1 = 0.0, t2 = 0.0;
    const unsigned int work_group_size = 1024;
    const unsigned int num_work_groups = ceilDivision(N, work_group_size);
    const unsigned int total_work_items = num_work_groups * work_group_size;
    const unsigned int lmem_per_work_group = sizeof(double) * work_group_size;
    if (verbosity >= 4) {
        t1 = second();
    }

    cl::Event event = (*dot_k)(cl::EnqueueArgs(*queue, cl::NDRange(total_work_items), cl::NDRange(work_group_size)), in1, in2, out, N, cl::Local(lmem_per_work_group));

    queue->enqueueReadBuffer(out, CL_TRUE, 0, sizeof(double) * num_work_groups, tmp);

    double gpu_sum = 0.0;
    for (unsigned int i = 0; i < num_work_groups; ++i) {
        gpu_sum += tmp[i];
    }

    if (verbosity >= 4) {
        event.wait();
        t2 = second();
        std::ostringstream oss;
        oss << "openclSolver dot_w time: " << t2 - t1;
        OpmLog::info(oss.str());
    }

    return gpu_sum;
}

double openclSolverBackend::norm_w(cl::Buffer in, cl::Buffer out)
{
    double t1 = 0.0, t2 = 0.0;
    const unsigned int work_group_size = 1024;
    const unsigned int num_work_groups = ceilDivision(N, work_group_size);
    const unsigned int total_work_items = num_work_groups * work_group_size;
    const unsigned int lmem_per_work_group = sizeof(double) * work_group_size;
    if (verbosity >= 4) {
        t1 = second();
    }

    cl::Event event = (*norm_k)(cl::EnqueueArgs(*queue, cl::NDRange(total_work_items), cl::NDRange(work_group_size)), in, out, N, cl::Local(lmem_per_work_group));

    queue->enqueueReadBuffer(out, CL_TRUE, 0, sizeof(double) * num_work_groups, tmp);

    double gpu_norm = 0.0;
    for (unsigned int i = 0; i < num_work_groups; ++i) {
        gpu_norm += tmp[i];
    }
    gpu_norm = sqrt(gpu_norm);

    if (verbosity >= 4) {
        event.wait();
        t2 = second();
        std::ostringstream oss;
        oss << "openclSolver norm_w time: " << t2 - t1;
        OpmLog::info(oss.str());
    }

    return gpu_norm;
}

void openclSolverBackend::axpy_w(cl::Buffer in, const double a, cl::Buffer out)
{
    double t1 = 0.0, t2 = 0.0;
    const unsigned int work_group_size = 32;
    const unsigned int num_work_groups = ceilDivision(N, work_group_size);
    const unsigned int total_work_items = num_work_groups * work_group_size;
    if (verbosity >= 4) {
        t1 = second();
    }

    cl::Event event = (*axpy_k)(cl::EnqueueArgs(*queue, cl::NDRange(total_work_items), cl::NDRange(work_group_size)), in, a, out, N);

    if (verbosity >= 4) {
        event.wait();
        t2 = second();
        std::ostringstream oss;
        oss << "openclSolver axpy_w time: " << t2 - t1;
        OpmLog::info(oss.str());
    }
}

void openclSolverBackend::custom_w(cl::Buffer p, cl::Buffer v, cl::Buffer r, const double omega, const double beta)
{
    double t1 = 0.0, t2 = 0.0;
    const unsigned int work_group_size = 32;
    const unsigned int num_work_groups = ceilDivision(N, work_group_size);
    const unsigned int total_work_items = num_work_groups * work_group_size;
    if (verbosity >= 4) {
        t1 = second();
    }

    cl::Event event = (*custom_k)(cl::EnqueueArgs(*queue, cl::NDRange(total_work_items), cl::NDRange(work_group_size)), p, v, r, omega, beta, N);

    if (verbosity >= 4) {
        event.wait();
        t2 = second();
        std::ostringstream oss;
        oss << "openclSolver custom_w time: " << t2 - t1;
        OpmLog::info(oss.str());
    }
}

void openclSolverBackend::spmv_blocked_w(cl::Buffer vals, cl::Buffer cols, cl::Buffer rows, cl::Buffer x, cl::Buffer b)
{
    double t1 = 0.0, t2 = 0.0;
    const unsigned int work_group_size = 32;
    const unsigned int num_work_groups = ceilDivision(N, work_group_size);
    const unsigned int total_work_items = num_work_groups * work_group_size;
    const unsigned int lmem_per_work_group = sizeof(double) * work_group_size;
    if (verbosity >= 4) {
        t1 = second();
    }

    cl::Event event = (*spmv_blocked_k)(cl::EnqueueArgs(*queue, cl::NDRange(total_work_items), cl::NDRange(work_group_size)), vals, cols, rows, Nb, x, b, block_size, cl::Local(lmem_per_work_group));

    if (verbosity >= 4) {
        event.wait();
        t2 = second();
        std::ostringstream oss;
        oss << "openclSolver spmv_blocked_w time: " << t2 - t1;
        OpmLog::info(oss.str());
    }
}


void openclSolverBackend::gpu_pbicgstab(WellContributions& wellContribs, BdaResult& res) {

    float it;
    double rho, rhop, beta, alpha, omega, tmp1, tmp2;
    double norm, norm_0;

    double t_total1, t_total2, t1 = 0.0, t2 = 0.0;
    double prec_time = 0.0, spmv_time = 0.0, well_time = 0.0, rest_time = 0.0;
    t_total1 = second();

    wellContribs.setOpenCLQueue(queue.get());
    wellContribs.setReordering(toOrder, true);

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

    norm = norm_w(d_r, d_tmp);
    norm_0 = norm;

    if (verbosity > 1) {
        std::ostringstream out;
        out << std::scientific << "openclSolver initial norm: " << norm_0;
        OpmLog::info(out.str());
    }

    t1 = second();
    for (it = 0.5; it < maxit; it += 0.5) {
        rhop = rho;
        rho = dot_w(d_rw, d_r, d_tmp);

        if (it > 1) {
            beta = (rho / rhop) * (alpha / omega);
            custom_w(d_p, d_v, d_r, omega, beta);
        }
        t2 = second();
        rest_time += t2 - t1;

        // pw = prec(p)
        t1 = second();
        prec->apply(d_p, d_pw);
        t2 = second();
        prec_time += t2 - t1;

        // v = A * pw
        t1 = second();
        spmv_blocked_w(d_Avals, d_Acols, d_Arows, d_pw, d_v);
        t2 = second();
        spmv_time += t2 - t1;

        // apply wellContributions
        if (wellContribs.getNumWells() > 0) {
            t1 = second();
            wellContribs.apply(d_pw, d_v);
            t2 = second();
            well_time += t2 - t1;
        }

        t1 = second();
        tmp1 = dot_w(d_rw, d_v, d_tmp);
        alpha = rho / tmp1;
        axpy_w(d_v, -alpha, d_r);      // r = r - alpha * v
        axpy_w(d_pw, alpha, d_x);      // x = x + alpha * pw
        norm = norm_w(d_r, d_tmp);
        t2 = second();
        rest_time += t2 - t1;

        if (norm < tolerance * norm_0) {
            break;
        }

        it += 0.5;

        // s = prec(r)
        t1 = second();
        prec->apply(d_r, d_s);
        t2 = second();
        prec_time += t2 - t1;

        // t = A * s
        t1 = second();
        spmv_blocked_w(d_Avals, d_Acols, d_Arows, d_s, d_t);
        t2 = second();
        spmv_time += t2 - t1;

        // apply wellContributions
        if (wellContribs.getNumWells() > 0) {
            t1 = second();
            wellContribs.apply(d_s, d_t);
            t2 = second();
            well_time += t2 - t1;
        }

        t1 = second();
        tmp1 = dot_w(d_t, d_r, d_tmp);
        tmp2 = dot_w(d_t, d_t, d_tmp);
        omega = tmp1 / tmp2;
        axpy_w(d_s, omega, d_x);     // x = x + omega * s
        axpy_w(d_t, -omega, d_r);    // r = r - omega * t
        norm = norm_w(d_r, d_tmp);
        t2 = second();
        rest_time += t2 - t1;

        if (norm < tolerance * norm_0) {
            break;
        }

        if (verbosity > 1) {
            std::ostringstream out;
            out << "it: " << it << std::scientific << ", norm: " << norm;
            OpmLog::info(out.str());
        }
    }

    t2 = second();
    t_total2 = second();
    rest_time += t2 - t1;

    res.iterations = std::min(it, (float)maxit);
    res.reduction = norm / norm_0;
    res.conv_rate  = static_cast<double>(pow(res.reduction, 1.0 / it));
    res.elapsed = t_total2 - t_total1;
    res.converged = (it != (maxit + 0.5));

    if (verbosity > 0) {
        std::ostringstream out;
        out << "=== converged: " << res.converged << ", conv_rate: " << res.conv_rate << ", time: " << res.elapsed << \
            ", time per iteration: " << res.elapsed / it << ", iterations: " << it;
        OpmLog::info(out.str());
    }
    if (verbosity >= 4) {
        std::ostringstream out;
        out << "openclSolver::ily_apply:   " << prec_time << "s\n";
        out << "openclSolver::spmv:        " << spmv_time << "s\n";
        out << "openclSolver::rest:        " << rest_time << "s\n";
        out << "openclSolver::total_solve: " << res.elapsed << "s\n";
        OpmLog::info(out.str());
    }
}


void openclSolverBackend::initialize(int N_, int nnz_, int dim, double *vals, int *rows, int *cols) {
    this->N = N_;
    this->nnz = nnz_;
    this->block_size = dim;
    this->nnzb = nnz_ / block_size / block_size;

    Nb = (N + dim - 1) / dim;
    std::ostringstream out;
    out << "Initializing GPU, matrix size: " << N << " blocks, nnzb: " << nnzb << "\n";
    out << "Maxit: " << maxit << std::scientific << ", tolerance: " << tolerance << "\n";
    OpmLog::info(out.str());
    out.str("");
    out.clear();

    int deviceID = 0;

    cl_int err = CL_SUCCESS;
    try {
        std::vector<cl::Platform> platforms;
        cl::Platform::get(&platforms);
        if (platforms.size() == 0)
        {
            OPM_THROW(std::logic_error, "Error openclSolver is selected but no OpenCL platforms are found");
        }
        out << "Found " << platforms.size() << " OpenCL platforms" << "\n\n";

        if (verbosity >= 1) {
            std::string platform_info;
            for (unsigned int platformID = 0; platformID < platforms.size(); ++platformID) {
                platforms[platformID].getInfo(CL_PLATFORM_NAME, &platform_info);
                out << "Platform name      : " << platform_info << "\n";
                platforms[platformID].getInfo(CL_PLATFORM_VENDOR, &platform_info);
                out << "Platform vendor    : " << platform_info << "\n";
                platforms[platformID].getInfo(CL_PLATFORM_VERSION, &platform_info);
                out << "Platform version   : " << platform_info << "\n";
                platforms[platformID].getInfo(CL_PLATFORM_PROFILE, &platform_info);
                out << "Platform profile   : " << platform_info << "\n";
                platforms[platformID].getInfo(CL_PLATFORM_EXTENSIONS, &platform_info);
                out << "Platform extensions: " << platform_info << "\n\n";
            }
        }
        OpmLog::info(out.str());
        out.str("");
        out.clear();

        cl_context_properties properties[] = {CL_CONTEXT_PLATFORM, (cl_context_properties)(platforms[deviceID])(), 0};
        context.reset(new cl::Context(CL_DEVICE_TYPE_GPU, properties));

        std::vector<cl::Device> devices = context->getInfo<CL_CONTEXT_DEVICES>();
        if (devices.size() == 0){
            OPM_THROW(std::logic_error, "Error openclSolver is selected but no OpenCL devices are found");
        }
        out << "Found " << devices.size() << " OpenCL devices" << "\n\n";

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
                devices[i].getInfo(CL_DEVICE_VERSION, &device_info);
                out << "CL_DEVICE_VERSION         : " << device_info << "\n";
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

        cl::Program::Sources source(1, std::make_pair(kernel_1, strlen(kernel_1)));  // what does this '1' mean? cl::Program::Sources is of type 'std::vector<std::pair<const char*, long unsigned int> >'
        source.emplace_back(std::make_pair(kernel_2, strlen(kernel_2)));
        source.emplace_back(std::make_pair(axpy_s, strlen(axpy_s)));
        source.emplace_back(std::make_pair(dot_1_s, strlen(dot_1_s)));
        source.emplace_back(std::make_pair(norm_s, strlen(norm_s)));
        source.emplace_back(std::make_pair(custom_s, strlen(custom_s)));
        source.emplace_back(std::make_pair(spmv_blocked_s, strlen(spmv_blocked_s)));
        source.emplace_back(std::make_pair(ILU_apply1_s, strlen(ILU_apply1_s)));
        source.emplace_back(std::make_pair(ILU_apply2_s, strlen(ILU_apply2_s)));
        cl::Program program_ = cl::Program(*context, source);

        program_.build(devices);

        cl::Event event;
        queue.reset(new cl::CommandQueue(*context, devices[0], 0, &err));

        prec->setOpenCLContext(context.get());
        prec->setOpenCLQueue(queue.get());

        rb = new double[N];
        tmp = new double[N];
#if COPY_ROW_BY_ROW
        vals_contiguous = new double[N];
#endif

        mat = (BlockedMatrix *)malloc(sizeof(BlockedMatrix));
        mat->Nb = Nb;
        mat->nnzbs = nnzb;
        mat->nnzValues = (Block*)vals;
        mat->colIndices = cols;
        mat->rowPointers = rows;

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
        dot_k.reset(new cl::make_kernel<cl::Buffer&, cl::Buffer&, cl::Buffer&, const unsigned int, cl::LocalSpaceArg>(cl::Kernel(program_, "dot_1")));
        norm_k.reset(new cl::make_kernel<cl::Buffer&, cl::Buffer&, const unsigned int, cl::LocalSpaceArg>(cl::Kernel(program_, "norm")));
        axpy_k.reset(new cl::make_kernel<cl::Buffer&, const double, cl::Buffer&, const unsigned int>(cl::Kernel(program_, "axpy")));
        custom_k.reset(new cl::make_kernel<cl::Buffer&, cl::Buffer&, cl::Buffer&, const double, const double, const unsigned int>(cl::Kernel(program_, "custom")));
        spmv_blocked_k.reset(new cl::make_kernel<cl::Buffer&, cl::Buffer&, cl::Buffer&, const unsigned int, cl::Buffer&, cl::Buffer&, const unsigned int, cl::LocalSpaceArg>(cl::Kernel(program_, "spmv_blocked")));
        ILU_apply1_k.reset(new cl::make_kernel<cl::Buffer&, cl::Buffer&, cl::Buffer&, const unsigned int, cl::Buffer&, cl::Buffer&, cl::Buffer&, const unsigned int, const unsigned int, cl::LocalSpaceArg>(cl::Kernel(program_, "ILU_apply1")));
        ILU_apply2_k.reset(new cl::make_kernel<cl::Buffer&, cl::Buffer&, cl::Buffer&, const unsigned int, cl::Buffer&, cl::Buffer&, cl::Buffer&, const unsigned int, const unsigned int, cl::LocalSpaceArg>(cl::Kernel(program_, "ILU_apply2")));

        prec->setKernels(ILU_apply1_k.get(), ILU_apply2_k.get());

    } catch (cl::Error error) {
        std::ostringstream oss;
        oss << "OpenCL Error: " << error.what() << "(" << error.err() << ")";
        OpmLog::error(oss.str());
    }


    initialized = true;
} // end initialize()

void openclSolverBackend::finalize() {
    delete[] rb;
    delete[] tmp;
#if COPY_ROW_BY_ROW
    delete[] vals_contiguous;
#endif
} // end finalize()


void openclSolverBackend::copy_system_to_gpu() {

    double t1 = 0.0, t2 = 0.0;
    if (verbosity > 2) {
        t1 = second();
    }

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
        t2 = second();
        std::ostringstream out;
        out << "openclSolver::copy_system_to_gpu(): " << t2 - t1 << " s";
        OpmLog::info(out.str());
    }
} // end copy_system_to_gpu()


// don't copy rowpointers and colindices, they stay the same
void openclSolverBackend::update_system_on_gpu() {

    double t1 = 0.0, t2 = 0.0;
    if (verbosity > 2) {
        t1 = second();
    }

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
        t2 = second();
        std::ostringstream out;
        out << "openclSolver::update_system_on_gpu(): " << t2 - t1 << " s";
        OpmLog::info(out.str());
    }
} // end update_system_on_gpu()


bool openclSolverBackend::analyse_matrix() {

    double t1 = 0.0, t2 = 0.0;

    if (verbosity > 2) {
        t1 = second();
    }

    bool success = prec->init(mat, block_size);
    int work_group_size = 32;
    int num_work_groups = ceilDivision(N, work_group_size);
    int total_work_items = num_work_groups * work_group_size;
    int lmem_per_work_group = work_group_size * sizeof(double);
    prec->setKernelParameters(work_group_size, total_work_items, lmem_per_work_group);

    toOrder = prec->getToOrder();
    fromOrder = prec->getFromOrder();
    rmat = prec->getRMat();

    if (verbosity > 2) {
        t2 = second();
        std::ostringstream out;
        out << "openclSolver::analyse_matrix(): " << t2 - t1 << " s";
        OpmLog::info(out.str());
    }

    analysis_done = true;

    return success;
} // end analyse_matrix()


void openclSolverBackend::update_system(double *vals, double *b) {
    double t1 = 0.0, t2 = 0.0;
    if (verbosity > 2) {
        t1 = second();
    }

    mat->nnzValues = (Block*)vals;
    //mat->nnzValues = static_cast<Block*>(vals);
    blocked_reorder_vector_by_pattern(mat->Nb, b, fromOrder, rb);

    if (verbosity > 2) {
        t2 = second();
        std::ostringstream out;
        out << "openclSolver::update_system(): " << t2 - t1 << " s";
        OpmLog::info(out.str());
    }
} // end update_system()


bool openclSolverBackend::create_preconditioner() {

    double t1 = 0.0, t2 = 0.0;
    if (verbosity > 2) {
        t1 = second();
    }

    bool result = prec->create_preconditioner(mat);

    if (verbosity > 2) {
        t2 = second();
        std::ostringstream out;
        out << "openclSolver::create_preconditioner(): " << t2 - t1 << " s";
        OpmLog::info(out.str());
    }
    return result;
} // end create_preconditioner()


void openclSolverBackend::solve_system(WellContributions& wellContribs, BdaResult &res) {
    // actually solve
    double t1 = 0.0, t2 = 0.0;
    if (verbosity > 2) {
        t1 = second();
    }

    gpu_pbicgstab(wellContribs, res);

    if (verbosity > 2) {
        t2 = second();
        std::ostringstream out;
        out << "openclSolver::solve_system(): " << t2 - t1 << " s";
        OpmLog::info(out.str());
    }

} // end solve_system()


// copy result to host memory
// caller must be sure that x is a valid array
void openclSolverBackend::get_result(double *x) {

    double t1 = 0.0, t2 = 0.0;
    if (verbosity > 2) {
        t1 = second();
    }

    queue->enqueueReadBuffer(d_x, CL_TRUE, 0, sizeof(double) * N, rb);
    blocked_reorder_vector_by_pattern(mat->Nb, rb, toOrder, x);

    if (verbosity > 2) {
        t2 = second();
        std::ostringstream out;
        out << "openclSolver::get_result(): " << t2 - t1 << " s";
        OpmLog::info(out.str());
    }
} // end get_result()



typedef BdaSolver::BdaSolverStatus BdaSolverStatus;

BdaSolverStatus openclSolverBackend::solve_system(int N_, int nnz_, int dim, double *vals, int *rows, int *cols, double *b, WellContributions& wellContribs, BdaResult &res) {
    if (initialized == false) {
        initialize(N_, nnz_,  dim, vals, rows, cols);
        if (analysis_done == false) {
            if (!analyse_matrix()) {
                return BdaSolverStatus::BDA_SOLVER_ANALYSIS_FAILED;
            }
        }
        update_system(vals, b);
        if (!create_preconditioner()) {
            return BdaSolverStatus::BDA_SOLVER_CREATE_PRECONDITIONER_FAILED;
        }
        copy_system_to_gpu();
    } else {
        update_system(vals, b);
        if (!create_preconditioner()) {
            return BdaSolverStatus::BDA_SOLVER_CREATE_PRECONDITIONER_FAILED;
        }
        update_system_on_gpu();
    }
    solve_system(wellContribs, res);
    return BdaSolverStatus::BDA_SOLVER_SUCCESS;
}

}

