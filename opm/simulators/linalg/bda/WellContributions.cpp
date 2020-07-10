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
#include "opm/simulators/linalg/bda/WellContributions.hpp"

namespace Opm
{


void WellContributions::alloc()
{
    if (num_std_wells > 0) {
#if HAVE_CUDA
        allocStandardWells();
#else
        OPM_THROW(std::logic_error, "Error cannot allocate on GPU for StandardWells because CUDA was not found by cmake");
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


    // apply StandardWells
    if (num_std_wells > 0) {
        OPM_THROW(std::logic_error, "Error StandardWells are not supported by openclSolver");
    }
}
#endif

void WellContributions::addMatrix(MatrixType type, int *colIndices, double *values, unsigned int val_size)
{
#if HAVE_CUDA
        addMatrixGpu(type, colIndices, values, val_size);
#else
        OPM_THROW(std::logic_error, "Error cannot add StandardWell matrix on GPU because CUDA was not found by cmake");
#endif
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
void WellContributions::setOpenCLQueue(cl::CommandQueue *queue_)
{
    this->queue = queue_;
}
#endif

} //namespace Opm

