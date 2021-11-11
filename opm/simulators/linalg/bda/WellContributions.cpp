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
#include <opm/simulators/linalg/bda/WellContributions.hpp>

#include <opm/common/OpmLog/OpmLog.hpp>
#include <opm/common/ErrorMacros.hpp>

#if HAVE_OPENCL
#include <opm/simulators/linalg/bda/opencl.hpp>
#include <opm/simulators/linalg/bda/openclKernels.hpp>
#include <opm/simulators/linalg/bda/openclSolverBackend.hpp>
#endif

#include <cstdlib>
#include <cstring>

namespace Opm
{

#if HAVE_OPENCL
using Opm::Accelerator::OpenclKernels;

struct OpenCLData {
    OpenCLData(std::vector<MultisegmentWellContribution*>& msegs) : multisegments(msegs) {}

    std::vector<MultisegmentWellContribution*>& multisegments;
    cl::Context* context = nullptr;
    cl::CommandQueue* queue = nullptr;
    Opm::Accelerator::stdwell_apply_kernel_type* kernel = nullptr;
    Opm::Accelerator::stdwell_apply_no_reorder_kernel_type* kernel_no_reorder = nullptr;
    std::vector<cl::Event> events;

    std::unique_ptr<cl::Buffer> d_Cnnzs_ocl, d_Dnnzs_ocl, d_Bnnzs_ocl;
    std::unique_ptr<cl::Buffer> d_Ccols_ocl, d_Bcols_ocl;
    std::unique_ptr<cl::Buffer> d_val_pointers_ocl;

    std::vector<double> h_x;
    std::vector<double> h_y;

    bool reorder = false;
    int* h_toOrder = nullptr;

    void setKernel(Opm::Accelerator::stdwell_apply_kernel_type *kernel_,
                   Opm::Accelerator::stdwell_apply_no_reorder_kernel_type *kernel_no_reorder_)
    {
        this->kernel = kernel_;
        this->kernel_no_reorder = kernel_no_reorder_;
    }

    void setOpenCLEnv(cl::Context *context_, cl::CommandQueue *queue_)
    {
        this->context = context_;
        this->queue = queue_;
    }

    void apply_stdwells(cl::Buffer d_x, cl::Buffer d_y, cl::Buffer d_toOrder, const WellContributions::Dimensions& dims)
    {
        if (reorder) {
            OpenclKernels::apply_stdwells_reorder(*d_Cnnzs_ocl, *d_Dnnzs_ocl, *d_Bnnzs_ocl, *d_Ccols_ocl, *d_Bcols_ocl,
                                                  d_x, d_y, d_toOrder, dims.dim, dims.dim_wells, *d_val_pointers_ocl, dims.num_std_wells);
        } else {
            OpenclKernels::apply_stdwells_no_reorder(*d_Cnnzs_ocl, *d_Dnnzs_ocl, *d_Bnnzs_ocl, *d_Ccols_ocl, *d_Bcols_ocl,
                                                     d_x, d_y, dims.dim, dims.dim_wells, *d_val_pointers_ocl, dims.num_std_wells);
        }
    }

    void apply_mswells(cl::Buffer d_x, cl::Buffer d_y, const WellContributions::Dimensions& dims)
    {
        if(h_x.empty()) {
            h_x.resize(dims.N);
            h_y.resize(dims.N);
        }

        events.resize(2);
        queue->enqueueReadBuffer(d_x, CL_FALSE, 0, sizeof(double) * dims.N, h_x.data(), nullptr, &events[0]);
        queue->enqueueReadBuffer(d_y, CL_FALSE, 0, sizeof(double) * dims.N, h_y.data(), nullptr, &events[1]);
        cl::WaitForEvents(events);
        events.clear();

        // actually apply MultisegmentWells
        for(Opm::MultisegmentWellContribution *well: multisegments){
            well->setReordering(h_toOrder, reorder);
            well->apply(h_x.data(), h_y.data());
        }

        // copy vector y from CPU to GPU
        events.resize(1);
        queue->enqueueWriteBuffer(d_y, CL_FALSE, 0, sizeof(double) * dims.N, h_y.data(), nullptr, &events[0]);
        events[0].wait();
        events.clear();
    }

    void apply(cl::Buffer& d_x, cl::Buffer& d_y, cl::Buffer& d_toOrder, const WellContributions::Dimensions& dims)
    {
        if (dims.num_std_wells > 0) {
            apply_stdwells(d_x, d_y, d_toOrder, dims);
        }

        if (dims.num_ms_wells > 0) {
            apply_mswells(d_x, d_y, dims);
        }
    }
};
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
    // delete MultisegmentWellContributions
    for (auto ms: multisegments) {
        delete ms;
    }
    multisegments.clear();

#if HAVE_CUDA
    if(cuda_gpu){
        freeCudaMemory(); // should come before 'delete[] h_x'
    }
#endif

    if (dimensions.num_std_wells > 0) {
        delete[] val_pointers;
    }
}

#if HAVE_OPENCL
template<unsigned int block_size>
void WellContributions::setOpenCLEnv(Accelerator::BdaSolver<block_size>& backend)
{
    const auto& openclBackend = static_cast<const Opm::Accelerator::openclSolverBackend<block_size>&>(backend);
    this->ocldata->setOpenCLEnv(openclBackend.context.get(), openclBackend.queue.get());
}

template void WellContributions::setOpenCLEnv<1>(Accelerator::BdaSolver<1>&);
template void WellContributions::setOpenCLEnv<2>(Accelerator::BdaSolver<2>&);
template void WellContributions::setOpenCLEnv<3>(Accelerator::BdaSolver<3>&);
template void WellContributions::setOpenCLEnv<4>(Accelerator::BdaSolver<4>&);
template void WellContributions::setOpenCLEnv<5>(Accelerator::BdaSolver<5>&);
template void WellContributions::setOpenCLEnv<6>(Accelerator::BdaSolver<6>&);

void WellContributions::apply(cl::Buffer& d_x, cl::Buffer& d_y, cl::Buffer& d_toOrder)
{
    this->ocldata->apply(d_x, d_y, d_toOrder, dimensions);
}

void WellContributions::setReordering(int* toOrder, bool reorder_)
{
    this->ocldata->h_toOrder = toOrder;
    this->ocldata->reorder = reorder_;
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
            ocldata->events.resize(2);
            ocldata->queue->enqueueWriteBuffer(*ocldata->d_Cnnzs_ocl, CL_FALSE, sizeof(double) * num_blocks_so_far * dimensions.dim * dimensions.dim_wells, sizeof(double) * val_size * dimensions.dim * dimensions.dim_wells, values, nullptr, &ocldata->events[0]);
            ocldata->queue->enqueueWriteBuffer(*ocldata->d_Ccols_ocl, CL_FALSE, sizeof(int) * num_blocks_so_far, sizeof(int) * val_size, colIndices, nullptr, &ocldata->events[1]);
            cl::WaitForEvents(ocldata->events);
            ocldata->events.clear();
            break;

        case MatrixType::D:
            ocldata->events.resize(1);
            ocldata->queue->enqueueWriteBuffer(*ocldata->d_Dnnzs_ocl, CL_FALSE, sizeof(double) * num_std_wells_so_far * dimensions.dim_wells * dimensions.dim_wells, sizeof(double) * dimensions.dim_wells * dimensions.dim_wells, values, nullptr, &ocldata->events[0]);
            ocldata->events[0].wait();
            ocldata->events.clear();
            break;

        case MatrixType::B:
            ocldata->events.resize(2);
            ocldata->queue->enqueueWriteBuffer(*ocldata->d_Bnnzs_ocl, CL_FALSE, sizeof(double) * num_blocks_so_far * dimensions.dim * dimensions.dim_wells, sizeof(double) * val_size * dimensions.dim * dimensions.dim_wells, values, nullptr, &ocldata->events[0]);
            ocldata->queue->enqueueWriteBuffer(*ocldata->d_Bcols_ocl, CL_FALSE, sizeof(int) * num_blocks_so_far, sizeof(int) * val_size, colIndices, nullptr, &ocldata->events[1]);
            cl::WaitForEvents(ocldata->events);
            ocldata->events.clear();

            val_pointers[num_std_wells_so_far] = num_blocks_so_far;
            if (num_std_wells_so_far == dimensions.num_std_wells - 1) {
                val_pointers[dimensions.num_std_wells] = dimensions.num_blocks;
                ocldata->events.resize(1);
                ocldata->queue->enqueueWriteBuffer(*ocldata->d_val_pointers_ocl, CL_FALSE, 0, sizeof(unsigned int) * (dimensions.num_std_wells + 1), val_pointers, nullptr, &ocldata->events[0]);
                ocldata->events[0].wait();
                ocldata->events.clear();
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
    dimensions.dim = dim_;
    dimensions.dim_wells = dim_wells_;

    if (dimensions.dim != 3 || dimensions.dim_wells != 4) {
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
    dimensions.num_blocks += numBlocks;
    dimensions.num_std_wells++;
}

void WellContributions::alloc()
{
    if (dimensions.num_std_wells > 0) {
        val_pointers = new unsigned int[dimensions.num_std_wells + 1];

#if HAVE_CUDA
        if(cuda_gpu){
            allocStandardWells();
        }
#endif

#if HAVE_OPENCL
        if(opencl_gpu){
            ocldata = std::make_unique<OpenCLData>(multisegments);
            ocldata->d_Cnnzs_ocl = std::make_unique<cl::Buffer>(*ocldata->context, CL_MEM_READ_WRITE, sizeof(double) * dimensions.num_blocks * dimensions.dim * dimensions.dim_wells);
            ocldata->d_Dnnzs_ocl = std::make_unique<cl::Buffer>(*ocldata->context, CL_MEM_READ_WRITE, sizeof(double) * dimensions.num_std_wells * dimensions.dim_wells * dimensions.dim_wells);
            ocldata->d_Bnnzs_ocl = std::make_unique<cl::Buffer>(*ocldata->context, CL_MEM_READ_WRITE, sizeof(double) * dimensions.num_blocks * dimensions.dim * dimensions.dim_wells);
            ocldata->d_Ccols_ocl = std::make_unique<cl::Buffer>(*ocldata->context, CL_MEM_READ_WRITE, sizeof(int) * dimensions.num_blocks);
            ocldata->d_Bcols_ocl = std::make_unique<cl::Buffer>(*ocldata->context, CL_MEM_READ_WRITE, sizeof(int) * dimensions.num_blocks);
            ocldata->d_val_pointers_ocl = std::make_unique<cl::Buffer>(*ocldata->context, CL_MEM_READ_WRITE, sizeof(unsigned int) * (dimensions.num_std_wells + 1));
        }
#endif
        allocated = true;
    }
}

void WellContributions::addMultisegmentWellContribution(unsigned int dim_, unsigned int dim_wells_,
        unsigned int Mb,
        std::vector<double> &Bvalues, std::vector<unsigned int> &BcolIndices, std::vector<unsigned int> &BrowPointers,
        unsigned int DnumBlocks, double *Dvalues, UMFPackIndex *DcolPointers, UMFPackIndex *DrowIndices,
        std::vector<double> &Cvalues)
{
    assert(dimensions.dim==dim_);
    MultisegmentWellContribution *well = new MultisegmentWellContribution(dim_, dim_wells_, Mb, Bvalues, BcolIndices, BrowPointers, DnumBlocks, Dvalues, DcolPointers, DrowIndices, Cvalues);
    multisegments.emplace_back(well);
    ++dimensions.num_ms_wells;
}

} //namespace Opm
