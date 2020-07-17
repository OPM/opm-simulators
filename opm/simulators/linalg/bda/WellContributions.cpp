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


void WellContributions::alloc()
{
    if (num_std_wells > 0) {
#if HAVE_CUDA
        allocStandardWells();
#elif HAVE_OPENCL
        d_Cnnzs = cl::Buffer(*context, CL_MEM_READ_ONLY, sizeof(double) * num_blocks * dim * dim_wells);
        d_Dnnzs = cl::Buffer(*context, CL_MEM_READ_ONLY, sizeof(double) * num_std_wells * dim_wells * dim_wells);
        d_Bnnzs = cl::Buffer(*context, CL_MEM_READ_ONLY, sizeof(double) * num_blocks * dim * dim_wells);
        d_Ccols = cl::Buffer(*context, CL_MEM_READ_ONLY, sizeof(int) * num_blocks);
        d_Bcols = cl::Buffer(*context, CL_MEM_READ_ONLY, sizeof(int) * num_blocks);
        d_val_pointers = cl::Buffer(*context, CL_MEM_READ_ONLY, sizeof(double) * (num_std_wells + 1));
#else
        OPM_THROW(std::logic_error, "Error cannot allocate on GPU because neither CUDA nor OpenCL were found by cmake");
#endif
    }
}

WellContributions::~WellContributions()
{
    if (h_x) {
        delete[] h_x;
        delete[] h_y;
    }

    // delete MultisegmentWellContributions
    for (auto ms : multisegments) {
        delete ms;
    }
    multisegments.clear();

#if HAVE_CUDA
    freeStandardWells();
#elif HAVE_OPENCL
    cl::ReleaseMemObject(d_Cnnzs);
    cl::ReleaseMemObject(d_Dnnzs);
    cl::ReleaseMemObject(d_Bnnzs);
    cl::ReleaseMemObject(d_Ccols);
    cl::ReleaseMemObject(d_Bcols);
    cl::ReleaseMemObject(d_val_pointers);
#endif
}


#if HAVE_OPENCL
void WellContributions::apply(cl::Buffer& d_x, cl::Buffer& d_y) {

    // apply MultisegmentWells
    if (num_ms_wells > 0) {
        // allocate pinned memory on host if not yet done
        if (h_x == nullptr) {
            h_x = new double[N];
            h_y = new double[N];
        }

        // copy vectors x and y from GPU to CPU
        queue->enqueueReadBuffer(d_x, CL_TRUE, 0, sizeof(double) * N, h_x);
        queue->enqueueReadBuffer(d_y, CL_TRUE, 0, sizeof(double) * N, h_y);

        // actually apply MultisegmentWells
        for (MultisegmentWellContribution *well : multisegments) {
            well->apply(h_x, h_y);
        }

        // copy vector y from CPU to GPU
        queue->enqueueWriteBuffer(d_y, CL_TRUE, 0, sizeof(double) * N, h_y);
    }

    queue->enqueueWriteBuffer(d_Cnnzs, CL_TRUE, 0, sizeof(double) * h_Cnnzs.size(), h_Cnnzs.data());
    queue->enqueueWriteBuffer(d_Dnnzs, CL_TRUE, 0, sizeof(double) * h_Dnnzs.size(), h_Dnnzs.data());
    queue->enqueueWriteBuffer(d_Bnnzs, CL_TRUE, 0, sizeof(double) * h_Bnnzs.size(), h_Bnnzs.data());
    queue->enqueueWriteBuffer(d_Ccols, CL_TRUE, 0, sizeof(int) * h_Ccols.size(), h_Ccols.data());
    queue->enqueueWriteBuffer(d_Bcols, CL_TRUE, 0, sizeof(int) * h_Bcols.size(), h_Bcols.data());
    queue->enqueueWriteBuffer(d_val_pointers, CL_TRUE, 0, sizeof(int) * (num_std_wells + 1), val_pointers);

    // apply StandardWells
    if (num_std_wells > 0) {
        cl::Event event;
        event = (*add_well_contributions)(cl::EnqueueArgs(*queue, cl::NDRange(total_work_items), cl::NDRange(work_group_size)), d_Cnnzs, d_Dnnzs, d_Bnnzs, d_Ccols, d_Bcols, d_x, d_y, dim, dim_wells, d_val_pointers, cl::Local(work_group_size*sizeof(double)), cl::Local(2*sizeof(double)*dim_wells), cl::Local(2*sizeof(double)*dim_wells));
    }
}
#endif

void WellContributions::addMatrix(MatrixType type, int *colIndices, double *values, unsigned int val_size)
{
    if (!allocated) {
        OPM_THROW(std::logic_error, "Error cannot add wellcontribution before allocating memory in WellContributions");
    }
#if HAVE_CUDA
    addMatrixGpu(type, colIndices, values, val_size);
#elif HAVE_OPENCL
    switch (type) {
    case MatrixType::C:
        h_Ccols.insert(h_Ccols.end(), colIndices, colIndices + val_size);
        h_Cnnzs.insert(h_Cnnzs.end(), values, values + val_size*dim*dim_wells);
        break;

    case MatrixType::D:
        h_Dnnzs.insert(h_Dnnzs.end(), values, values + val_size*dim_wells*dim_wells);
        break;

    case MatrixType::B:
        h_Bcols.insert(h_Bcols.end(), colIndices, colIndices + val_size);
        h_Bnnzs.insert(h_Bnnzs.end(), values, values + val_size*dim*dim_wells);
        val_pointers[num_std_wells_so_far] = num_blocks_so_far;

        if(num_std_wells_so_far == num_std_wells - 1){
            val_pointers[num_std_wells] = num_blocks;
        }
        break;
    default:
        OPM_THROW(std::logic_error, "Error unsupported matrix ID for WellContributions::addMatrix()");
    }   
#else
    OPM_THROW(std::logic_error, "Error cannot add StandardWell matrix on GPU because neither CUDA nor OpenCL were found by cmake");
#endif
    if (MatrixType::B == type) {
        num_blocks_so_far += val_size;
        num_std_wells_so_far++;
    }
}


void WellContributions::setBlockSize(unsigned int dim_, unsigned int dim_wells_)
{
    dim = dim_;
    dim_wells = dim_wells_;
}

void WellContributions::addNumBlocks(unsigned int nnz)
{
    if (allocated) {
        OPM_THROW(std::logic_error, "Error cannot add more sizes after allocated in WellContributions");
    }
    num_blocks += nnz;
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

#if HAVE_OPENCL
void WellContributions::setOpenCLContext(cl::Context *context_)
{
    this->context = context_;
}

void WellContributions::setOpenCLQueue(cl::CommandQueue *queue_)
{
    this->queue = queue_;
}

void WellContributions::setKernelParameters(const unsigned int work_group_size_, const unsigned int total_work_items_)
{
    this->work_group_size = work_group_size_;
    this->total_work_items = total_work_items_;
}

void WellContributions::setKernel(cl::make_kernel<cl::Buffer&, cl::Buffer&, cl::Buffer&, cl::Buffer&, cl::Buffer&, cl::Buffer&, cl::Buffer&, const unsigned int, const unsigned int, cl::Buffer&, cl::LocalSpaceArg, cl::LocalSpaceArg, cl::LocalSpaceArg> *add_well_contributions_)
{
    this->add_well_contributions = add_well_contributions_;
}
#endif

} //namespace Opm

