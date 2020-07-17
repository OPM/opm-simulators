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
#include <opm/simulators/linalg/bda/openclKernels.hpp>
#include "opm/simulators/linalg/bda/WellContributions.hpp"

namespace Opm
{

WellContributions::WellContributions(std::string gpu_mode){
    if(gpu_mode.compare("cusparse") == 0){
        cuda_gpu = true;
    }

    if(gpu_mode.compare("opencl") == 0){
        opencl_gpu = true;
    }
}

void WellContributions::alloc()
{
    if (num_std_wells > 0) {
#if HAVE_CUDA
        if(cuda_gpu){
            allocStandardWells();
        }
#endif

#if HAVE_OPENCL
        if(opencl_gpu){
            h_Cnnzs_ocl = new double[num_blocks * dim * dim_wells];
            h_Dnnzs_ocl = new double[num_std_wells * dim_wells * dim_wells];
            h_Bnnzs_ocl = new double[num_blocks * dim * dim_wells];
            h_Ccols_ocl = new int[num_blocks];
            h_Bcols_ocl = new int[num_blocks];
            val_pointers = new unsigned int[num_std_wells + 1];

            allocated = true;
        }
#endif

#if !HAVE_CUDA && !HAVE_OPENCL
        OPM_THROW(std::logic_error, "Error cannot allocate on GPU because neither CUDA nor OpenCL were found by cmake");
#endif
    }
}

WellContributions::~WellContributions()
{
    // delete MultisegmentWellContributions
    for (auto ms : multisegments) {
        delete ms;
    }
    multisegments.clear();

#if HAVE_CUDA
    if(cuda_gpu){
        freeCudaMemory(); // should come before 'delete[] h_x'
    }
#endif

#if HAVE_OPENCL
    if (h_x_ocl) {
        delete[] h_x_ocl;
        delete[] h_y_ocl;
    }

    if(opencl_gpu){
        if(num_std_wells > 0){
            delete[] h_Cnnzs_ocl;
            delete[] h_Dnnzs_ocl;
            delete[] h_Bnnzs_ocl;
            delete[] h_Ccols_ocl;
            delete[] h_Bcols_ocl;
            delete[] val_pointers;
        }
    }
#endif
}

#if HAVE_OPENCL

void WellContributions::init(cl::Context *context){
    d_Cnnzs_ocl = cl::Buffer(*context, CL_MEM_READ_WRITE, sizeof(double) * num_blocks * dim * dim_wells);
    d_Dnnzs_ocl = cl::Buffer(*context, CL_MEM_READ_WRITE, sizeof(double) * num_std_wells * dim_wells * dim_wells);
    d_Bnnzs_ocl = cl::Buffer(*context, CL_MEM_READ_WRITE, sizeof(double) * num_blocks * dim * dim_wells);
    d_Ccols_ocl = cl::Buffer(*context, CL_MEM_READ_WRITE, sizeof(int) * num_blocks);
    d_Bcols_ocl = cl::Buffer(*context, CL_MEM_READ_WRITE, sizeof(int) * num_blocks);
    d_val_pointers_ocl = cl::Buffer(*context, CL_MEM_READ_WRITE, sizeof(unsigned int) * (num_std_wells + 1));
}

void WellContributions::copyDataToGPU(cl::CommandQueue *queue){
    cl::Event event;

    queue->enqueueWriteBuffer(d_Cnnzs_ocl, CL_TRUE, 0, sizeof(double) * num_blocks * dim * dim_wells, h_Cnnzs_ocl);
    queue->enqueueWriteBuffer(d_Dnnzs_ocl, CL_TRUE, 0, sizeof(double) * num_std_wells * dim_wells * dim_wells, h_Dnnzs_ocl);
    queue->enqueueWriteBuffer(d_Bnnzs_ocl, CL_TRUE, 0, sizeof(double) * num_blocks * dim * dim_wells, h_Bnnzs_ocl);
    queue->enqueueWriteBuffer(d_Ccols_ocl, CL_TRUE, 0, sizeof(int) * num_blocks, h_Ccols_ocl);
    queue->enqueueWriteBuffer(d_Bcols_ocl, CL_TRUE, 0, sizeof(int) * num_blocks, h_Bcols_ocl);
    queue->enqueueWriteBuffer(d_val_pointers_ocl, CL_TRUE, 0, sizeof(unsigned int) * (num_std_wells + 1), val_pointers, nullptr, &event);
    event.wait();
}

void WellContributions::applyMSWell(cl::CommandQueue *queue, cl::Buffer& d_x, cl::Buffer& d_y) {
    // apply MultisegmentWells
    if (num_ms_wells > 0) {
        // allocate pinned memory on host if not yet done
        if (h_x_ocl == nullptr) {
            h_x_ocl = new double[N];
            h_y_ocl = new double[N];
        }

        // copy vectors x and y from GPU to CPU
        queue->enqueueReadBuffer(d_x, CL_TRUE, 0, sizeof(double) * N, h_x_ocl);
        queue->enqueueReadBuffer(d_y, CL_TRUE, 0, sizeof(double) * N, h_y_ocl);

        // actually apply MultisegmentWells
        for (MultisegmentWellContribution *well : multisegments) {
            well->apply(h_x_ocl, h_y_ocl);
        }

        // copy vector y from CPU to GPU
        queue->enqueueWriteBuffer(d_y, CL_TRUE, 0, sizeof(double) * N, h_y_ocl);
    }
}

void WellContributions::applyStdWell(cl::CommandQueue *queue, cl::Buffer& d_x, cl::Buffer& d_y, kernel_type *kernel){
    const unsigned int work_group_size = 32;
    const unsigned int total_work_items = num_std_wells * work_group_size;
    const unsigned int lmem1 = sizeof(double) * work_group_size;
    const unsigned int lmem2 = sizeof(double) * dim_wells;

    cl::Event event;
    event = (*kernel)(cl::EnqueueArgs(*queue, cl::NDRange(total_work_items), cl::NDRange(work_group_size)),
                      d_Cnnzs_ocl, d_Dnnzs_ocl, d_Bnnzs_ocl, d_Ccols_ocl, d_Bcols_ocl, d_x, d_y, dim, dim_wells,
                      d_val_pointers_ocl, cl::Local(lmem1), cl::Local(lmem2), cl::Local(lmem2));
    event.wait();
}

void WellContributions::apply(cl::CommandQueue *queue, cl::Buffer& d_x, cl::Buffer& d_y, kernel_type *kernel){
    if(num_std_wells > 0){
        applyStdWell(queue, d_x, d_y, kernel);
    }

    if(num_ms_wells > 0){
        applyMSWell(queue, d_x, d_y);
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
                std::copy(colIndices, colIndices + val_size, h_Ccols_ocl + num_blocks_so_far);
                std::copy(values, values + val_size*dim*dim_wells, h_Cnnzs_ocl + num_blocks_so_far*dim*dim_wells);
                break;

            case MatrixType::D:
                std::copy(values, values + dim_wells*dim_wells, h_Dnnzs_ocl + num_std_wells_so_far*dim_wells*dim_wells);
                break;

            case MatrixType::B:
                std::copy(colIndices, colIndices + val_size, h_Bcols_ocl + num_blocks_so_far);
                std::copy(values, values + val_size*dim*dim_wells, h_Bnnzs_ocl + num_blocks_so_far*dim*dim_wells);
                val_pointers[num_std_wells_so_far] = num_blocks_so_far;

                if(num_std_wells_so_far == num_std_wells - 1){
                    val_pointers[num_std_wells] = num_blocks;
                }
                break;

            default:
                OPM_THROW(std::logic_error, "Error unsupported matrix ID for WellContributions::addMatrix()");
        }

        if (MatrixType::B == type) {
            num_blocks_so_far += val_size;
            num_std_wells_so_far++;
        }
    }

#endif

#if !HAVE_CUDA && !HAVE_OPENCL
    OPM_THROW(std::logic_error, "Error cannot add StandardWell matrix on GPU because neither CUDA nor OpenCL were found by cmake");
#endif
}


void WellContributions::setBlockSize(unsigned int dim_, unsigned int dim_wells_)
{
    dim = dim_;
    dim_wells = dim_wells_;
}

void WellContributions::addNumBlocks(unsigned int numBlocks)
{
    if (allocated) {
        OPM_THROW(std::logic_error, "Error cannot add more sizes after allocated in WellContributions");
    }
    num_blocks += numBlocks;
    num_std_wells++;
}

void WellContributions::addMultisegmentWellContribution(unsigned int dim, unsigned int dim_wells,
        unsigned int Nb, unsigned int Mb,
        unsigned int BnumBlocks, std::vector<double> &Bvalues, std::vector<unsigned int> &BcolIndices, std::vector<unsigned int> &BrowPointers,
        unsigned int DnumBlocks, double *Dvalues, int *DcolPointers, int *DrowIndices,
        std::vector<double> &Cvalues)
{
    this->N = Nb * dim;
    MultisegmentWellContribution *well = new MultisegmentWellContribution(dim, dim_wells, Nb, Mb, BnumBlocks, Bvalues, BcolIndices, BrowPointers, DnumBlocks, Dvalues, DcolPointers, DrowIndices, Cvalues);
    multisegments.emplace_back(well);
    ++num_ms_wells;
}


void WellContributions::setReordering(int *toOrder_, bool reorder_)
{
    this->toOrder = toOrder_;
    this->reorder = reorder_;
    for (auto& ms : multisegments) {
        ms->setReordering(toOrder_, reorder_);
    }
}

} //namespace Opm

