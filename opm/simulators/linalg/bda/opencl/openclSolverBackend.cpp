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

#include <opm/simulators/linalg/bda/opencl/opencl.hpp>
#include <opm/simulators/linalg/bda/opencl/openclKernels.hpp>
#include <opm/simulators/linalg/bda/opencl/openclSolverBackend.hpp>
#include <opm/simulators/linalg/bda/opencl/openclWellContributions.hpp>

#include <opm/simulators/linalg/bda/BdaResult.hpp>


// iff true, the nonzeroes of the matrix are copied row-by-row into a contiguous, pinned memory array, then a single GPU memcpy is done
// otherwise, the nonzeroes of the matrix are assumed to be in a contiguous array, and a single GPU memcpy is enough
#define COPY_ROW_BY_ROW 0

namespace Opm
{
namespace Accelerator
{

using Opm::OpmLog;
using Dune::Timer;

template <unsigned int block_size>
openclSolverBackend<block_size>::openclSolverBackend(int verbosity_, int maxit_, double tolerance_, unsigned int platformID_, unsigned int deviceID_, bool opencl_ilu_parallel_, std::string linsolver) : BdaSolver<block_size>(verbosity_, maxit_, tolerance_, platformID_, deviceID_), opencl_ilu_parallel(opencl_ilu_parallel_) {

    bool use_cpr, use_isai;

    if (linsolver.compare("ilu0") == 0) {
        use_cpr = false;
        use_isai = false;
    } else if (linsolver.compare("cpr_quasiimpes") == 0) {
        use_cpr = true;
        use_isai = false;
    } else if (linsolver.compare("isai") == 0) {
        use_cpr = false;
        use_isai = true;
    } else if (linsolver.compare("cpr_trueimpes") == 0) {
        OPM_THROW(std::logic_error, "Error openclSolver does not support --linerar-solver=cpr_trueimpes");
    } else {
        OPM_THROW(std::logic_error, "Error unknown value for argument --linear-solver, " + linsolver);
    }

    using PreconditionerType = typename Preconditioner<block_size>::PreconditionerType;
    if (use_cpr) {
        prec = Preconditioner<block_size>::create(PreconditionerType::CPR, verbosity, opencl_ilu_parallel);
    } else if (use_isai) {
        prec = Preconditioner<block_size>::create(PreconditionerType::BISAI, verbosity, opencl_ilu_parallel);
    } else {
        prec = Preconditioner<block_size>::create(PreconditionerType::BILU0, verbosity, opencl_ilu_parallel);
    }

    std::ostringstream out;
    try {
        std::vector<cl::Platform> platforms;
        cl::Platform::get(&platforms);
        if (platforms.empty()) {
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

        platforms[platformID].getDevices(CL_DEVICE_TYPE_ALL, &devices);

        if (devices.empty()) {
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

        // removed all unused devices
        if (deviceID != 0)
        {
            devices[0] = devices[deviceID];
        }
        devices.resize(1);

        context = std::make_shared<cl::Context>(devices[0]);
        queue.reset(new cl::CommandQueue(*context, devices[0], 0, &err));

        OpenclKernels::init(context.get(), queue.get(), devices, verbosity);

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
}

template <unsigned int block_size>
openclSolverBackend<block_size>::openclSolverBackend(int verbosity_, int maxit_, double tolerance_, bool opencl_ilu_parallel_) :
    BdaSolver<block_size>(verbosity_, maxit_, tolerance_), opencl_ilu_parallel(opencl_ilu_parallel_)
{
    // prec = std::make_unique<BILU0<block_size> >(opencl_ilu_parallel, verbosity_);
    // cpr = std::make_unique<CPR<block_size> >(verbosity_, opencl_ilu_parallel, /*use_amg=*/false);
}

template <unsigned int block_size>
void openclSolverBackend<block_size>::setOpencl(std::shared_ptr<cl::Context>& context_, std::shared_ptr<cl::CommandQueue>& queue_) {
    context = context_;
    queue = queue_;
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
    events.resize(5);
    queue->enqueueFillBuffer(d_p, 0, 0, sizeof(double) * N, nullptr, &events[0]);
    queue->enqueueFillBuffer(d_v, 0, 0, sizeof(double) * N, nullptr, &events[1]);
    rho = 1.0;
    alpha = 1.0;
    omega = 1.0;

    queue->enqueueCopyBuffer(d_b, d_r, 0, 0, sizeof(double) * N, nullptr, &events[2]);
    queue->enqueueCopyBuffer(d_r, d_rw, 0, 0, sizeof(double) * N, nullptr, &events[3]);
    queue->enqueueCopyBuffer(d_r, d_p, 0, 0, sizeof(double) * N, nullptr, &events[4]);

    cl::WaitForEvents(events);
    events.clear();
    if (err != CL_SUCCESS) {
        // enqueueWriteBuffer is C and does not throw exceptions like C++ OpenCL
        OPM_THROW(std::logic_error, "openclSolverBackend OpenCL enqueue[Fill|Copy]Buffer error");
    }

    norm = OpenclKernels::norm(d_r, d_tmp, N);
    norm_0 = norm;

    if (verbosity > 1) {
        std::ostringstream out;
        out << std::scientific << "openclSolver initial norm: " << norm_0;
        OpmLog::info(out.str());
    }

    if (verbosity >= 3) {
        t_rest.start();
    }
    for (it = 0.5; it < maxit; it += 0.5) {
        rhop = rho;
        rho = OpenclKernels::dot(d_rw, d_r, d_tmp, N);

        if (it > 1) {
            beta = (rho / rhop) * (alpha / omega);
            OpenclKernels::custom(d_p, d_v, d_r, omega, beta, N);
        }
        if (verbosity >= 3) {
            queue->finish();
            t_rest.stop();
            t_prec.start();
        }

        // pw = prec(p)
        prec->apply(d_p, d_pw);
        if (verbosity >= 3) {
            queue->finish();
            t_prec.stop();
            t_spmv.start();
        }

        // v = A * pw
        OpenclKernels::spmv(d_Avals, d_Acols, d_Arows, d_pw, d_v, Nb, block_size);
        if (verbosity >= 3) {
            queue->finish();
            t_spmv.stop();
            t_well.start();
        }

        // apply wellContributions
        if(wellContribs.getNumWells() > 0){
            static_cast<WellContributionsOCL&>(wellContribs).apply(d_pw, d_v);
        }
        if(verbosity >= 3) {
            queue->finish();
            t_well.stop();
            t_rest.start();
        }

        tmp1 = OpenclKernels::dot(d_rw, d_v, d_tmp, N);
        alpha = rho / tmp1;
        OpenclKernels::axpy(d_v, -alpha, d_r, N);      // r = r - alpha * v
        OpenclKernels::axpy(d_pw, alpha, d_x, N);      // x = x + alpha * pw
        norm = OpenclKernels::norm(d_r, d_tmp, N);
        if (verbosity >= 3) {
            queue->finish();
            t_rest.stop();
        }

        if (norm < tolerance * norm_0) {
            break;
        }

        it += 0.5;

        // s = prec(r)
        if (verbosity >= 3) {
            t_prec.start();
        }
        prec->apply(d_r, d_s);
        if (verbosity >= 3) {
            queue->finish();
            t_prec.stop();
            t_spmv.start();
        }

        // t = A * s
        OpenclKernels::spmv(d_Avals, d_Acols, d_Arows, d_s, d_t, Nb, block_size);
        if(verbosity >= 3){
            queue->finish();
            t_spmv.stop();
            t_well.start();
        }

        // apply wellContributions
        if(wellContribs.getNumWells() > 0){
            static_cast<WellContributionsOCL&>(wellContribs).apply(d_s, d_t);
        }
        if (verbosity >= 3) {
            queue->finish();
            t_well.stop();
            t_rest.start();
        }

        tmp1 = OpenclKernels::dot(d_t, d_r, d_tmp, N);
        tmp2 = OpenclKernels::dot(d_t, d_t, d_tmp, N);
        omega = tmp1 / tmp2;
        OpenclKernels::axpy(d_s, omega, d_x, N);     // x = x + omega * s
        OpenclKernels::axpy(d_t, -omega, d_r, N);    // r = r - omega * t
        norm = OpenclKernels::norm(d_r, d_tmp, N);
        if (verbosity >= 3) {
            queue->finish();
            t_rest.stop();
        }

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
    if (verbosity >= 3) {
        std::ostringstream out;
        out << "openclSolver::prec_apply:  " << t_prec.elapsed() << " s\n";
        out << "wellContributions::apply:  " << t_well.elapsed() << " s\n";
        out << "openclSolver::spmv:        " << t_spmv.elapsed() << " s\n";
        out << "openclSolver::rest:        " << t_rest.elapsed() << " s\n";
        out << "openclSolver::total_solve: " << res.elapsed << " s\n";
        OpmLog::info(out.str());
    }
}


template <unsigned int block_size>
void openclSolverBackend<block_size>::initialize(std::shared_ptr<BlockedMatrix> matrix, std::shared_ptr<BlockedMatrix> jacMatrix) {
    this->Nb = matrix->Nb;
    this->N = Nb * block_size;
    this->nnzb = matrix->nnzbs;
    this->nnz = nnzb * block_size * block_size;

    if (jacMatrix) {
        useJacMatrix = true;
    }

    std::ostringstream out;
    out << "Initializing GPU, matrix size: " << Nb << " blockrows, nnzb: " << nnzb << "\n";
    if (useJacMatrix) {
        out << "Blocks in ILU matrix: " << jacMatrix->nnzbs << "\n";
    }
    out << "Maxit: " << maxit << std::scientific << ", tolerance: " << tolerance << "\n";
    out << "PlatformID: " << platformID << ", deviceID: " << deviceID << "\n";
    OpmLog::info(out.str());
    out.str("");
    out.clear();

    try {
        prec->setOpencl(context, queue);

#if COPY_ROW_BY_ROW
        vals_contiguous.resize(nnz);
#endif
        mat = matrix;
        jacMat = jacMatrix;

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
void openclSolverBackend<block_size>::copy_system_to_gpu() {
    Timer t;
    events.resize(5);

#if COPY_ROW_BY_ROW
    int sum = 0;
    for (int i = 0; i < Nb; ++i) {
        int size_row = mat->rowPointers[i + 1] - mat->rowPointers[i];
        memcpy(vals_contiguous.data() + sum, mat->nnzValues + sum, size_row * sizeof(double) * block_size * block_size);
        sum += size_row * block_size * block_size;
    }
    err = queue->enqueueWriteBuffer(d_Avals, CL_TRUE, 0, sizeof(double) * nnz, vals_contiguous.data(), nullptr, &events[0]);
#else
    err = queue->enqueueWriteBuffer(d_Avals, CL_TRUE, 0, sizeof(double) * nnz, mat->nnzValues, nullptr, &events[0]);
#endif

    err |= queue->enqueueWriteBuffer(d_Acols, CL_TRUE, 0, sizeof(int) * nnzb, mat->colIndices, nullptr, &events[1]);
    err |= queue->enqueueWriteBuffer(d_Arows, CL_TRUE, 0, sizeof(int) * (Nb + 1), mat->rowPointers, nullptr, &events[2]);
    err |= queue->enqueueWriteBuffer(d_b, CL_TRUE, 0, sizeof(double) * N, h_b, nullptr, &events[3]);
    err |= queue->enqueueFillBuffer(d_x, 0, 0, sizeof(double) * N, nullptr, &events[4]);

    cl::WaitForEvents(events);
    events.clear();
    if (err != CL_SUCCESS) {
        // enqueueWriteBuffer is C and does not throw exceptions like C++ OpenCL
        OPM_THROW(std::logic_error, "openclSolverBackend OpenCL enqueueWriteBuffer error");
    }

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
    events.resize(3);

#if COPY_ROW_BY_ROW
    int sum = 0;
    for (int i = 0; i < Nb; ++i) {
        int size_row = mat->rowPointers[i + 1] - mat->rowPointers[i];
        memcpy(vals_contiguous.data() + sum, mat->nnzValues + sum, size_row * sizeof(double) * block_size * block_size);
        sum += size_row * block_size * block_size;
    }
    err = queue->enqueueWriteBuffer(d_Avals, CL_TRUE, 0, sizeof(double) * nnz, vals_contiguous.data(), nullptr, &events[0]);
#else
    err = queue->enqueueWriteBuffer(d_Avals, CL_TRUE, 0, sizeof(double) * nnz, mat->nnzValues, nullptr, &events[0]);
#endif

    err |= queue->enqueueWriteBuffer(d_b, CL_TRUE, 0, sizeof(double) * N, h_b, nullptr, &events[1]);
    err |= queue->enqueueFillBuffer(d_x, 0, 0, sizeof(double) * N, nullptr, &events[2]);

    cl::WaitForEvents(events);
    events.clear();
    if (err != CL_SUCCESS) {
        // enqueueWriteBuffer is C and does not throw exceptions like C++ OpenCL
        OPM_THROW(std::logic_error, "openclSolverBackend OpenCL enqueueWriteBuffer error");
    }

    if (verbosity > 2) {
        std::ostringstream out;
        out << "openclSolver::update_system_on_gpu(): " << t.stop() << " s";
        OpmLog::info(out.str());
    }
} // end update_system_on_gpu()


template <unsigned int block_size>
bool openclSolverBackend<block_size>::analyze_matrix() {
    Timer t;

    bool success;
    if (useJacMatrix)
        success = prec->analyze_matrix(mat.get(), jacMat.get());
    else
        success = prec->analyze_matrix(mat.get());

    if (verbosity > 2) {
        std::ostringstream out;
        out << "openclSolver::analyze_matrix(): " << t.stop() << " s";
        OpmLog::info(out.str());
    }

    analysis_done = true;

    return success;
} // end analyze_matrix()


template <unsigned int block_size>
void openclSolverBackend<block_size>::update_system(double *vals, double *b) {
    Timer t;

    mat->nnzValues = vals;
    h_b = b;

    if (verbosity > 2) {
        std::ostringstream out;
        out << "openclSolver::update_system(): " << t.stop() << " s";
        OpmLog::info(out.str());
    }
} // end update_system()


template <unsigned int block_size>
bool openclSolverBackend<block_size>::create_preconditioner() {
    Timer t;

    bool result;
    if (useJacMatrix)
        result = prec->create_preconditioner(mat.get(), jacMat.get());
    else
        result = prec->create_preconditioner(mat.get());

    if (verbosity > 2) {
        std::ostringstream out;
        out << "openclSolver::create_preconditioner(): " << t.stop() << " s";
        OpmLog::info(out.str());
    }
    return result;
} // end create_preconditioner()


template <unsigned int block_size>
void openclSolverBackend<block_size>::solve_system(WellContributions &wellContribs, BdaResult &res) {
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

    queue->enqueueReadBuffer(d_x, CL_TRUE, 0, sizeof(double) * N, x);

    if (verbosity > 2) {
        std::ostringstream out;
        out << "openclSolver::get_result(): " << t.stop() << " s";
        OpmLog::info(out.str());
    }
} // end get_result()


template <unsigned int block_size>
SolverStatus openclSolverBackend<block_size>::solve_system(std::shared_ptr<BlockedMatrix> matrix,
                                                           double *b,
                                                           std::shared_ptr<BlockedMatrix> jacMatrix,
                                                           WellContributions& wellContribs,
                                                           BdaResult &res)
{
    if (initialized == false) {
        initialize(matrix, jacMatrix);
        if (analysis_done == false) {
            if (!analyze_matrix()) {
                return SolverStatus::BDA_SOLVER_ANALYSIS_FAILED;
            }
        }
        update_system(matrix->nnzValues, b);
        if (!create_preconditioner()) {
            return SolverStatus::BDA_SOLVER_CREATE_PRECONDITIONER_FAILED;
        }
        copy_system_to_gpu();
    } else {
        update_system(matrix->nnzValues, b);
        if (!create_preconditioner()) {
            return SolverStatus::BDA_SOLVER_CREATE_PRECONDITIONER_FAILED;
        }
        update_system_on_gpu();
    }
    solve_system(wellContribs, res);
    return SolverStatus::BDA_SOLVER_SUCCESS;
}


#define INSTANTIATE_BDA_FUNCTIONS(n)                                        \
template openclSolverBackend<n>::openclSolverBackend(                       \
    int, int, double, unsigned int, unsigned int, bool, std::string); \
template openclSolverBackend<n>::openclSolverBackend(int, int, double, bool); \
template void openclSolverBackend<n>::setOpencl(std::shared_ptr<cl::Context>&, std::shared_ptr<cl::CommandQueue>&);

INSTANTIATE_BDA_FUNCTIONS(1);
INSTANTIATE_BDA_FUNCTIONS(2);
INSTANTIATE_BDA_FUNCTIONS(3);
INSTANTIATE_BDA_FUNCTIONS(4);
INSTANTIATE_BDA_FUNCTIONS(5);
INSTANTIATE_BDA_FUNCTIONS(6);

#undef INSTANTIATE_BDA_FUNCTIONS

} // namespace Accelerator
} // namespace Opm
