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


namespace Opm
{

using Accelerator::OpenclKernels;

void WellContributionsOCL::setOpenCLEnv(cl::Context* context_, cl::CommandQueue* queue_) {
    this->context = context_;
    this->queue = queue_;
}


void WellContributionsOCL::apply_stdwells(cl::Buffer d_x, cl::Buffer d_y){
    OpenclKernels::apply_stdwells(*d_Cnnzs_ocl, *d_Dnnzs_ocl, *d_Bnnzs_ocl, *d_Ccols_ocl, *d_Bcols_ocl,
        d_x, d_y, dim, dim_wells, *d_val_pointers_ocl, num_std_wells);
}

void WellContributionsOCL::apply_mswells(cl::Buffer d_x, cl::Buffer d_y){
    if (h_x.empty()) {
        h_x.resize(N);
        h_y.resize(N);
    }

    events.resize(2);
    queue->enqueueReadBuffer(d_x, CL_FALSE, 0, sizeof(double) * N, h_x.data(), nullptr, &events[0]);
    queue->enqueueReadBuffer(d_y, CL_FALSE, 0, sizeof(double) * N, h_y.data(), nullptr, &events[1]);
    cl::WaitForEvents(events);
    events.clear();

    // actually apply MultisegmentWells
    for (auto& well : multisegments) {
        well->apply(h_x.data(), h_y.data());
    }

    // copy vector y from CPU to GPU
    events.resize(1);
    queue->enqueueWriteBuffer(d_y, CL_FALSE, 0, sizeof(double) * N, h_y.data(), nullptr, &events[0]);
    events[0].wait();
    events.clear();
}

void WellContributionsOCL::apply(cl::Buffer d_x, cl::Buffer d_y){
    if(num_std_wells > 0){
        apply_stdwells(d_x, d_y);
    }

    if(num_ms_wells > 0){
        apply_mswells(d_x, d_y);
    }
}

void WellContributionsOCL::APIaddMatrix(MatrixType type,
                                        int* colIndices,
                                        double* values,
                                        unsigned int val_size)
{
    if (!allocated) {
        OPM_THROW(std::logic_error, "Error cannot add wellcontribution before allocating memory in WellContributions");
    }

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
        OPM_THROW(std::logic_error, "Error unsupported matrix ID for WellContributionsOCL::addMatrix()");
    }
}

void WellContributionsOCL::APIalloc()
{
    d_Cnnzs_ocl = std::make_unique<cl::Buffer>(*context, CL_MEM_READ_WRITE, sizeof(double) * num_blocks * dim * dim_wells);
    d_Dnnzs_ocl = std::make_unique<cl::Buffer>(*context, CL_MEM_READ_WRITE, sizeof(double) * num_std_wells * dim_wells * dim_wells);
    d_Bnnzs_ocl = std::make_unique<cl::Buffer>(*context, CL_MEM_READ_WRITE, sizeof(double) * num_blocks * dim * dim_wells);
    d_Ccols_ocl = std::make_unique<cl::Buffer>(*context, CL_MEM_READ_WRITE, sizeof(int) * num_blocks);
    d_Bcols_ocl = std::make_unique<cl::Buffer>(*context, CL_MEM_READ_WRITE, sizeof(int) * num_blocks);
    d_val_pointers_ocl = std::make_unique<cl::Buffer>(*context, CL_MEM_READ_WRITE, sizeof(unsigned int) * (num_std_wells + 1));
}

} //namespace Opm
