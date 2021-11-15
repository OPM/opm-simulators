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

#include <config.h> // CMake
#include <cstdlib>
#include <cstring>
#include <opm/common/OpmLog/OpmLog.hpp>
#include <opm/common/ErrorMacros.hpp>

#include <opm/simulators/linalg/bda/WellContributions.hpp>

namespace Opm
{

#if HAVE_OPENCL
using Opm::Accelerator::OpenclKernels;
#endif

WellContributions::WellContributions(std::string accelerator_mode, bool useWellConn){
    if(accelerator_mode.compare("cusparse") == 0){
        cuda_gpu = true;
    }
    else if(accelerator_mode.compare("opencl") == 0){
        opencl_gpu = true;
    }
    else if(accelerator_mode.compare("fpga") == 0){
        // unused for FPGA, but must be defined to avoid error
    }
    else if(accelerator_mode.compare("amgcl") == 0){
        if (!useWellConn) {
            OPM_THROW(std::logic_error, "Error amgcl requires --matrix-add-well-contributions=true");
        }
    }
    else{
        OPM_THROW(std::logic_error, "Invalid accelerator mode");
    }
}

WellContributions::~WellContributions()
{
    multisegments.clear();

#if HAVE_CUDA
    if(cuda_gpu){
        freeCudaMemory(); // should come before 'delete[] h_x'
    }
#endif

#if HAVE_OPENCL
    if(opencl_gpu){
        if(num_ms_wells > 0){
            delete[] h_x;
            delete[] h_y;
        }
    }
#endif
}

#if HAVE_OPENCL
void WellContributions::setOpenCLEnv(cl::Context *context_, cl::CommandQueue *queue_){
    this->context = context_;
    this->queue = queue_;
}

void WellContributions::setKernel(Opm::Accelerator::stdwell_apply_kernel_type *kernel_,
                                  Opm::Accelerator::stdwell_apply_no_reorder_kernel_type *kernel_no_reorder_){
    this->kernel = kernel_;
    this->kernel_no_reorder = kernel_no_reorder_;
}

void WellContributions::setReordering(int *h_toOrder_, bool reorder_)
{
    this->h_toOrder = h_toOrder_;
    this->reorder = reorder_;
}

void WellContributions::apply_stdwells(cl::Buffer d_x, cl::Buffer d_y, cl::Buffer d_toOrder){
    if (reorder) {
        OpenclKernels::apply_stdwells_reorder(*d_Cnnzs_ocl, *d_Dnnzs_ocl, *d_Bnnzs_ocl, *d_Ccols_ocl, *d_Bcols_ocl,
            d_x, d_y, d_toOrder, dim, dim_wells, *d_val_pointers_ocl, num_std_wells);
    } else {
        OpenclKernels::apply_stdwells_no_reorder(*d_Cnnzs_ocl, *d_Dnnzs_ocl, *d_Bnnzs_ocl, *d_Ccols_ocl, *d_Bcols_ocl,
            d_x, d_y, dim, dim_wells, *d_val_pointers_ocl, num_std_wells);
    }
}

void WellContributions::apply_mswells(cl::Buffer d_x, cl::Buffer d_y){
    if(h_x == nullptr){
        h_x = new double[N];
        h_y = new double[N];
    }

    events.resize(2);
    queue->enqueueReadBuffer(d_x, CL_FALSE, 0, sizeof(double) * N, h_x, nullptr, &events[0]);
    queue->enqueueReadBuffer(d_y, CL_FALSE, 0, sizeof(double) * N, h_y, nullptr, &events[1]);
    cl::WaitForEvents(events);
    events.clear();

    // actually apply MultisegmentWells
    for (auto& well : multisegments) {
        well->setReordering(h_toOrder, reorder);
        well->apply(h_x, h_y);
    }

    // copy vector y from CPU to GPU
    events.resize(1);
    queue->enqueueWriteBuffer(d_y, CL_FALSE, 0, sizeof(double) * N, h_y, nullptr, &events[0]);
    events[0].wait();
    events.clear();
}

void WellContributions::apply(cl::Buffer d_x, cl::Buffer d_y, cl::Buffer d_toOrder){
    if(num_std_wells > 0){
        apply_stdwells(d_x, d_y, d_toOrder);
    }

    if(num_ms_wells > 0){
        apply_mswells(d_x, d_y);
    }
}
#endif

void WellContributions::addMatrix([[maybe_unused]] MatrixType type, [[maybe_unused]] int *colIndices, [[maybe_unused]] double *values, [[maybe_unused]] unsigned int val_size)
{
    if (!allocated) {
        OPM_THROW(std::logic_error, "Error cannot add wellcontribution before allocating memory in WellContributions");
    }

#if HAVE_CUDA
    if(cuda_gpu){
        addMatrixGpu(type, colIndices, values, val_size);
    }
#endif

#if HAVE_OPENCL
    if(opencl_gpu){
        switch (type) {
        case MatrixType::C:
            events.resize(2);
            queue->enqueueWriteBuffer(*d_Cnnzs_ocl, CL_FALSE, sizeof(double) * num_blocks_so_far * dim * dim_wells, sizeof(double) * val_size * dim * dim_wells, values, nullptr, &events[0]);
            queue->enqueueWriteBuffer(*d_Ccols_ocl, CL_FALSE, sizeof(int) * num_blocks_so_far, sizeof(int) * val_size, colIndices, nullptr, &events[1]);
            cl::WaitForEvents(events);
            events.clear();
            break;

        case MatrixType::D:
            events.resize(1);
            queue->enqueueWriteBuffer(*d_Dnnzs_ocl, CL_FALSE, sizeof(double) * num_std_wells_so_far * dim_wells * dim_wells, sizeof(double) * dim_wells * dim_wells, values, nullptr, &events[0]);
            events[0].wait();
            events.clear();
            break;

        case MatrixType::B:
            events.resize(2);
            queue->enqueueWriteBuffer(*d_Bnnzs_ocl, CL_FALSE, sizeof(double) * num_blocks_so_far * dim * dim_wells, sizeof(double) * val_size * dim * dim_wells, values, nullptr, &events[0]);
            queue->enqueueWriteBuffer(*d_Bcols_ocl, CL_FALSE, sizeof(int) * num_blocks_so_far, sizeof(int) * val_size, colIndices, nullptr, &events[1]);
            cl::WaitForEvents(events);
            events.clear();

            val_pointers[num_std_wells_so_far] = num_blocks_so_far;
            if (num_std_wells_so_far == num_std_wells - 1) {
                val_pointers[num_std_wells] = num_blocks;
                events.resize(1);
                queue->enqueueWriteBuffer(*d_val_pointers_ocl, CL_FALSE, 0, sizeof(unsigned int) * (num_std_wells + 1), val_pointers.data(), nullptr, &events[0]);
                events[0].wait();
                events.clear();
            }
            break;

        default:
            OPM_THROW(std::logic_error, "Error unsupported matrix ID for WellContributions::addMatrix()");
        }
    }
#endif

    if(MatrixType::B == type) {
        num_blocks_so_far += val_size;
        num_std_wells_so_far++;
    }

#if !HAVE_CUDA && !HAVE_OPENCL
    OPM_THROW(std::logic_error, "Error cannot add StandardWell matrix on GPU because neither CUDA nor OpenCL were found by cmake");
#endif
}

void WellContributions::setBlockSize(unsigned int dim_, unsigned int dim_wells_)
{
    dim = dim_;
    dim_wells = dim_wells_;

    if(dim != 3 || dim_wells != 4){
        std::ostringstream oss;
        oss << "WellContributions::setBlockSize error: dim and dim_wells must be equal to 3 and 4, repectivelly, otherwise the add well contributions kernel won't work.\n";
        OPM_THROW(std::logic_error, oss.str());
    }
}

void WellContributions::addNumBlocks(unsigned int numBlocks)
{
    if (allocated) {
        OPM_THROW(std::logic_error, "Error cannot add more sizes after allocated in WellContributions");
    }
    num_blocks += numBlocks;
    num_std_wells++;
}

void WellContributions::alloc()
{
    if (num_std_wells > 0) {
        val_pointers.resize(num_std_wells+1);

#if HAVE_CUDA
        if(cuda_gpu){
            allocStandardWells();
        }
#endif

#if HAVE_OPENCL
        if(opencl_gpu){
            d_Cnnzs_ocl = std::make_unique<cl::Buffer>(*context, CL_MEM_READ_WRITE, sizeof(double) * num_blocks * dim * dim_wells);
            d_Dnnzs_ocl = std::make_unique<cl::Buffer>(*context, CL_MEM_READ_WRITE, sizeof(double) * num_std_wells * dim_wells * dim_wells);
            d_Bnnzs_ocl = std::make_unique<cl::Buffer>(*context, CL_MEM_READ_WRITE, sizeof(double) * num_blocks * dim * dim_wells);
            d_Ccols_ocl = std::make_unique<cl::Buffer>(*context, CL_MEM_READ_WRITE, sizeof(int) * num_blocks);
            d_Bcols_ocl = std::make_unique<cl::Buffer>(*context, CL_MEM_READ_WRITE, sizeof(int) * num_blocks);
            d_val_pointers_ocl = std::make_unique<cl::Buffer>(*context, CL_MEM_READ_WRITE, sizeof(unsigned int) * (num_std_wells + 1));
        }
#endif
        allocated = true;
    }
}

void WellContributions::addMultisegmentWellContribution(unsigned int dim_,
                                                        unsigned int dim_wells_,
                                                        unsigned int Mb,
                                                        std::vector<double>& Bvalues,
                                                        std::vector<unsigned int>& BcolIndices,
                                                        std::vector<unsigned int>& BrowPointers,
                                                        unsigned int DnumBlocks,
                                                        double* Dvalues,
                                                        UMFPackIndex* DcolPointers,
                                                        UMFPackIndex* DrowIndices,
                                                        std::vector<double>& Cvalues)
{
    assert(dim==dim_);
    multisegments.push_back(std::make_unique<MultisegmentWellContribution>(dim_,
                                                                           dim_wells_,
                                                                           Mb,
                                                                           Bvalues,
                                                                           BcolIndices,
                                                                           BrowPointers,
                                                                           DnumBlocks,
                                                                           Dvalues,
                                                                           DcolPointers,
                                                                           DrowIndices,
                                                                           Cvalues));
    ++num_ms_wells;
}

} //namespace Opm
