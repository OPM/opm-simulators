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
}

/*
#if HAVE_OPENCL
void WellContributions::applyMSWell(cl::Buffer& d_x, cl::Buffer& d_y) {
    // apply MultisegmentWells
    if (num_ms_wells > 0) {
        h_x_ocl.reserve(N);
        h_y_ocl.reserve(N);

        // copy vectors x and y from GPU to CPU
        queue->enqueueReadBuffer(d_x, CL_TRUE, 0, sizeof(double) * N, h_x_ocl.data());
        queue->enqueueReadBuffer(d_y, CL_TRUE, 0, sizeof(double) * N, h_y_ocl.data());

        // actually apply MultisegmentWells
        for (MultisegmentWellContribution *well : multisegments) {
            well->apply(h_x_ocl.data(), h_y_ocl.data());
        }

        // copy vector y from CPU to GPU
        queue->enqueueWriteBuffer(d_y, CL_TRUE, 0, sizeof(double) * N, h_y_ocl.data());
    }
}
#endif
*/

void WellContributions::addMatrix([[maybe_unused]] MatrixType type, [[maybe_unused]] int *colIndices, [[maybe_unused]] double *values, [[maybe_unused]] unsigned int val_size)
{

#if HAVE_CUDA
    if(cuda_gpu){
        if (!allocated) {
            OPM_THROW(std::logic_error, "Error cannot add wellcontribution before allocating memory in WellContributions");
        }
        addMatrixGpu(type, colIndices, values, val_size);
    }
#endif

#if HAVE_OPENCL
    if(opencl_gpu){
        if(h_val_pointers_ocl.empty()){
            h_val_pointers_ocl.push_back(0);
        }

        switch (type) {
        case MatrixType::C:
            h_Ccols_ocl.insert(h_Ccols_ocl.end(), colIndices, colIndices + val_size);
            h_Cnnzs_ocl.insert(h_Cnnzs_ocl.end(), values, values + val_size * dim * dim_wells);
            break;

        case MatrixType::D:
            h_Dnnzs_ocl.insert(h_Dnnzs_ocl.end(), values, values + dim_wells * dim_wells);
            break;

        case MatrixType::B:
            h_Bcols_ocl.insert(h_Bcols_ocl.end(), colIndices, colIndices + val_size);
            h_Bnnzs_ocl.insert(h_Bnnzs_ocl.end(), values, values + val_size * dim * dim_wells);
            h_val_pointers_ocl.push_back(h_val_pointers_ocl.back() + val_size);
            break;

        default:
            OPM_THROW(std::logic_error, "Error unsupported matrix ID for WellContributions::addMatrix()");
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

#if HAVE_CUDA
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
        allocStandardWells();
        allocated = true;
    }
}
#endif

void WellContributions::addMultisegmentWellContribution(unsigned int dim_, unsigned int dim_wells_,
        unsigned int Nb, unsigned int Mb,
        unsigned int BnumBlocks, std::vector<double> &Bvalues, std::vector<unsigned int> &BcolIndices, std::vector<unsigned int> &BrowPointers,
        unsigned int DnumBlocks, double *Dvalues, int *DcolPointers, int *DrowIndices,
        std::vector<double> &Cvalues)
{
    assert(dim==dim_);
    this->N = Nb * dim_;
    MultisegmentWellContribution *well = new MultisegmentWellContribution(dim_, dim_wells_, Nb, Mb, BnumBlocks, Bvalues, BcolIndices, BrowPointers, DnumBlocks, Dvalues, DcolPointers, DrowIndices, Cvalues);
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

