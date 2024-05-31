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
#include <opm/simulators/linalg/bda/opencl/openclWellContributions.hpp>

#include <opm/common/OpmLog/OpmLog.hpp>
#include <opm/common/ErrorMacros.hpp>

#include <opm/simulators/linalg/bda/MultisegmentWellContribution.hpp>

namespace Opm {

using Accelerator::OpenclKernels;

template<class Scalar>
void WellContributionsOCL<Scalar>::
setOpenCLEnv(cl::Context* context_, cl::CommandQueue* queue_)
{
    this->context = context_;
    this->queue = queue_;
}

template<class Scalar>
void WellContributionsOCL<Scalar>::apply_stdwells(cl::Buffer d_x, cl::Buffer d_y)
{
    OpenclKernels<Scalar>::apply_stdwells(*d_Cnnzs_ocl, *d_Dnnzs_ocl, *d_Bnnzs_ocl,
                                          *d_Ccols_ocl, *d_Bcols_ocl,
                                          d_x, d_y, this->dim, this->dim_wells,
                                          *d_val_pointers_ocl, this->num_std_wells);
}

template<class Scalar>
void WellContributionsOCL<Scalar>::apply_mswells(cl::Buffer d_x, cl::Buffer d_y)
{
    if (h_x.empty()) {
        h_x.resize(this->N);
        h_y.resize(this->N);
    }

    events.resize(2);
    queue->enqueueReadBuffer(d_x, CL_FALSE, 0, sizeof(Scalar) * this->N,
                             h_x.data(), nullptr, &events[0]);
    queue->enqueueReadBuffer(d_y, CL_FALSE, 0, sizeof(Scalar) * this->N,
                             h_y.data(), nullptr, &events[1]);
    cl::WaitForEvents(events);
    events.clear();

    // actually apply MultisegmentWells
    for (auto& well : this->multisegments) {
        well->apply(h_x.data(), h_y.data());
    }

    // copy vector y from CPU to GPU
    events.resize(1);
    queue->enqueueWriteBuffer(d_y, CL_FALSE, 0, sizeof(Scalar) * this->N,
                              h_y.data(), nullptr, &events[0]);
    events[0].wait();
    events.clear();
}

template<class Scalar>
void WellContributionsOCL<Scalar>::apply(cl::Buffer d_x, cl::Buffer d_y)
{
    if (this->num_std_wells > 0){
        apply_stdwells(d_x, d_y);
    }

    if (this->num_ms_wells > 0) {
        apply_mswells(d_x, d_y);
    }
}

template<class Scalar>
void WellContributionsOCL<Scalar>::
APIaddMatrix(MatrixType type,
             int* colIndices,
             Scalar* values,
             unsigned int val_size)
{
    if (!this->allocated) {
        OPM_THROW(std::logic_error, "Error cannot add wellcontribution before allocating memory in WellContributions");
    }

    switch (type) {
    case MatrixType::C:
        events.resize(2);
        queue->enqueueWriteBuffer(*d_Cnnzs_ocl, CL_FALSE,
                                  sizeof(Scalar) * this->num_blocks_so_far * this->dim * this->dim_wells,
                                  sizeof(Scalar) * val_size * this->dim * this->dim_wells,
                                  values, nullptr, &events[0]);
        queue->enqueueWriteBuffer(*d_Ccols_ocl, CL_FALSE,
                                  sizeof(int) * this->num_blocks_so_far,
                                  sizeof(int) * val_size, colIndices, nullptr, &events[1]);
        cl::WaitForEvents(events);
        events.clear();
        break;

    case MatrixType::D:
        events.resize(1);
        queue->enqueueWriteBuffer(*d_Dnnzs_ocl, CL_FALSE,
                                  sizeof(Scalar) * this->num_std_wells_so_far * this->dim_wells * this->dim_wells,
                                  sizeof(Scalar) * this->dim_wells * this->dim_wells,
                                  values, nullptr, &events[0]);
        events[0].wait();
        events.clear();
        break;

    case MatrixType::B:
        events.resize(2);
        queue->enqueueWriteBuffer(*d_Bnnzs_ocl, CL_FALSE,
                                  sizeof(Scalar) * this->num_blocks_so_far * this->dim * this->dim_wells,
                                  sizeof(Scalar) * val_size * this->dim * this->dim_wells,
                                  values, nullptr, &events[0]);
        queue->enqueueWriteBuffer(*d_Bcols_ocl, CL_FALSE,
                                  sizeof(int) * this->num_blocks_so_far, sizeof(int) * val_size,
                                  colIndices, nullptr, &events[1]);
        cl::WaitForEvents(events);
        events.clear();

        this->val_pointers[this->num_std_wells_so_far] = this->num_blocks_so_far;
        if (this->num_std_wells_so_far == this->num_std_wells - 1) {
            this->val_pointers[this->num_std_wells] = this->num_blocks;
            events.resize(1);
            queue->enqueueWriteBuffer(*d_val_pointers_ocl, CL_FALSE, 0,
                                      sizeof(unsigned int) * (this->num_std_wells + 1),
                                      this->val_pointers.data(), nullptr, &events[0]);
            events[0].wait();
            events.clear();
        }
        break;

    default:
        OPM_THROW(std::logic_error, "Error unsupported matrix ID for WellContributionsOCL::addMatrix()");
    }
}

template<class Scalar>
void WellContributionsOCL<Scalar>::APIalloc()
{
    d_Cnnzs_ocl = std::make_unique<cl::Buffer>(*context, CL_MEM_READ_WRITE,
                                               sizeof(Scalar) * this->num_blocks * this->dim * this->dim_wells);
    d_Dnnzs_ocl = std::make_unique<cl::Buffer>(*context, CL_MEM_READ_WRITE,
                                               sizeof(Scalar) * this->num_std_wells * this->dim_wells * this->dim_wells);
    d_Bnnzs_ocl = std::make_unique<cl::Buffer>(*context, CL_MEM_READ_WRITE,
                                               sizeof(Scalar) * this->num_blocks * this->dim * this->dim_wells);
    d_Ccols_ocl = std::make_unique<cl::Buffer>(*context, CL_MEM_READ_WRITE, sizeof(int) * this->num_blocks);
    d_Bcols_ocl = std::make_unique<cl::Buffer>(*context, CL_MEM_READ_WRITE, sizeof(int) * this->num_blocks);
    d_val_pointers_ocl = std::make_unique<cl::Buffer>(*context, CL_MEM_READ_WRITE,
                                                      sizeof(unsigned int) * (this->num_std_wells + 1));
}

template class WellContributionsOCL<double>;

} // namespace Opm
