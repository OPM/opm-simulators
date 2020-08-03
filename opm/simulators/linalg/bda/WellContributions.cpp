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
#include <iostream>
#include <algorithm>

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
#endif

#if HAVE_OPENCL
        d_Cnnzs_ocl = new double[num_blocks * dim * dim_wells];
        d_Dnnzs_ocl = new double[num_std_wells * dim_wells * dim_wells];
        d_Bnnzs_ocl = new double[num_blocks * dim * dim_wells];
        d_Ccols_ocl = new int[num_blocks];
        d_Bcols_ocl = new int[num_blocks];
        val_pointers = new unsigned int[num_std_wells + 1];

        //d_Cnnzs_ocl = new double[num_blocks * dim * dim_wells];
        //d_Dnnzs_ocl = new double[num_std_wells * dim_wells * dim_wells];
        //d_Bnnzs_ocl = new double[num_blocks * dim * dim_wells];
        //d_Ccols_ocl = new int[num_blocks];
        //d_Bcols_ocl = new int[num_blocks];
        //val_pointers = new unsigned int[num_std_wells + 1];
        
        allocated = true;
#endif

#if !HAVE_CUDA && !HAVE_OPENCL
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

    if(num_std_wells > 0){
#if HAVE_CUDA
        freeStandardWells();
#endif

#if HAVE_OPENCL
        delete[] d_Cnnzs_ocl;
        delete[] d_Dnnzs_ocl;
        delete[] d_Bnnzs_ocl;
        delete[] d_Ccols_ocl;
        delete[] d_Bcols_ocl; 
        delete[] val_pointers;
        
        //delete[] d_Cnnzs_ocl;
        //delete[] d_Dnnzs_ocl;
        //delete[] d_Bnnzs_ocl;
        //delete[] d_Ccols_ocl;
        //delete[] d_Bcols_ocl; 
        //delete[] val_pointers;
#endif
    }
}


#if HAVE_OPENCL
/*
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
}
*/

void WellContributions::getParams(unsigned int *num_blocks_, unsigned int *num_std_wells_, unsigned int *dim_, unsigned int *dim_wells_){
    *num_blocks_ = num_blocks;
    *num_std_wells_ = num_std_wells;
    *dim_ = dim;
    *dim_wells_ = dim_wells;
}

void WellContributions::getData(double **valsC, double **valsD, double **valsB, int **colsC, int **colsB, unsigned int **val_pointers_){
    std::copy(d_Cnnzs_ocl, d_Cnnzs_ocl + num_blocks*dim*dim_wells, *valsC);
    std::copy(d_Dnnzs_ocl, d_Dnnzs_ocl + num_std_wells*dim_wells*dim_wells, *valsD);
    std::copy(d_Bnnzs_ocl, d_Bnnzs_ocl + num_blocks*dim*dim_wells, *valsB);
    std::copy(d_Ccols_ocl, d_Ccols_ocl + num_blocks, *colsC);
    std::copy(d_Bcols_ocl, d_Bcols_ocl + num_blocks, *colsB);
    std::copy(val_pointers, val_pointers + num_std_wells + 1, *val_pointers_);
}

#endif

void WellContributions::addMatrix([[maybe_unused]] MatrixType type, [[maybe_unused]]int *colIndices, [[maybe_unused]] double *values, [[maybe_unused]] unsigned int val_size)
{
    if (!allocated) {
        OPM_THROW(std::logic_error, "Error cannot add wellcontribution before allocating memory in WellContributions");
    }
#if HAVE_CUDA
    addMatrixGpu(type, colIndices, values, val_size);
#elif HAVE_OPENCL
    switch (type) {
        case MatrixType::C:
            std::copy(colIndices, colIndices + val_size, d_Ccols_ocl + num_blocks_so_far);
            std::copy(values, values + val_size*dim*dim_wells, d_Cnnzs_ocl + num_blocks_so_far*dim*dim_wells);
            break;

        case MatrixType::D:
            std::copy(values, values + dim_wells*dim_wells, d_Dnnzs_ocl + num_std_wells_so_far*dim_wells*dim_wells);
            break;

        case MatrixType::B:
            std::copy(colIndices, colIndices + val_size, d_Bcols_ocl + num_blocks_so_far);
            std::copy(values, values + val_size*dim*dim_wells, d_Bnnzs_ocl + num_blocks_so_far*dim*dim_wells);
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

void WellContributions::addNumBlocks(unsigned int numBlocks)
{
    if (allocated) {
        OPM_THROW(std::logic_error, "Error cannot add more sizes after allocated in WellContributions");
    }
    num_blocks += numBlocks;
    num_std_wells++;
}

void WellContributions::addMultisegmentWellContribution(unsigned int dim_, unsigned int dim_wells_,
        unsigned int Nb, unsigned int Mb,
        unsigned int BnumBlocks, std::vector<double> &Bvalues, std::vector<unsigned int> &BcolIndices, std::vector<unsigned int> &BrowPointers,
        unsigned int DnumBlocks, double *Dvalues, int *DcolPointers, int *DrowIndices,
        std::vector<double> &Cvalues)
{
    this->N = Nb * dim_;
    MultisegmentWellContribution *well = new MultisegmentWellContribution(dim, dim_wells_, Nb, Mb, BnumBlocks, Bvalues, BcolIndices, BrowPointers, DnumBlocks, Dvalues, DcolPointers, DrowIndices, Cvalues);
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

