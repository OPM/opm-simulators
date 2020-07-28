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
#elif HAVE_OPENCL
        h_valsC = new double[num_blocks * dim * dim_wells];
        h_valsD = new double[num_std_wells * dim_wells * dim_wells];
        h_valsB = new double[num_blocks * dim * dim_wells];
        h_colsC = new int[num_blocks];
        h_colsB = new int[num_blocks];
        val_pointers = new unsigned int[num_std_wells + 1];
        
        allocated = true;
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

    if(num_std_wells > 0){
#if HAVE_CUDA
        freeStandardWells();
#elif HAVE_OPENCL
        delete[] h_valsC;
        delete[] h_valsD;
        delete[] h_valsB;
        delete[] h_colsC;
        delete[] h_colsB; 
        delete[] val_pointers;
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
    std::copy(h_valsC, h_valsC + num_blocks*dim*dim_wells, *valsC);
    std::copy(h_valsD, h_valsD + num_std_wells*dim_wells*dim_wells, *valsD);
    std::copy(h_valsB, h_valsB + num_blocks*dim*dim_wells, *valsB);
    std::copy(h_colsC, h_colsC + num_blocks, *colsC);
    std::copy(h_colsB, h_colsB + num_blocks, *colsB);
    std::copy(val_pointers, val_pointers + num_std_wells + 1, *val_pointers_);
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
            std::copy(colIndices, colIndices + val_size, h_colsC + num_blocks_so_far);
            std::copy(values, values + val_size*dim*dim_wells, h_valsC + num_blocks_so_far*dim*dim_wells);
            break;

        case MatrixType::D:
            std::copy(values, values + dim_wells*dim_wells, h_valsD + num_std_wells_so_far*dim_wells*dim_wells);
            break;

        case MatrixType::B:
            std::copy(colIndices, colIndices + val_size, h_colsB + num_blocks_so_far);
            std::copy(values, values + val_size*dim*dim_wells, h_valsB + num_blocks_so_far*dim*dim_wells);
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

