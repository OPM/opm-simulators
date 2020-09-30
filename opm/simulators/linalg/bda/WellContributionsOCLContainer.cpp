/*
  Copyright 2019 Equinor ASA

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

#include <config.h>

#include <opm/common/OpmLog/OpmLog.hpp>
#include <opm/common/ErrorMacros.hpp>
#include <dune/common/timer.hh>

#include <opm/simulators/linalg/bda/WellContributionsOCLContainer.hpp>

namespace bda
{
    using Opm::OpmLog;
    using Dune::Timer;

    WellContributionsOCLContainer::WellContributionsOCLContainer(){
        multisegments.reset(new mswVecT());
    }

    WellContributionsOCLContainer::~WellContributionsOCLContainer(){
        if(num_ms_wells > 0){
            for (auto ms : *multisegments) {
                delete ms;
            }
        }
    }

    void WellContributionsOCLContainer::init(Opm::WellContributions &wellContribs, int N_, int Nb_){
        N = N_;
        Nb = Nb_;
        dim = wellContribs.dim;
        dim_wells = wellContribs.dim_wells;

        if(!wellContribs.h_val_pointers_ocl.empty()){
            num_blocks = wellContribs.h_Ccols_ocl.size();
            num_std_wells = wellContribs.h_val_pointers_ocl.size() - 1;

            s.Cnnzs = cl::Buffer(*context, CL_MEM_READ_WRITE, sizeof(double) * wellContribs.h_Cnnzs_ocl.size());
            s.Dnnzs = cl::Buffer(*context, CL_MEM_READ_WRITE, sizeof(double) * wellContribs.h_Dnnzs_ocl.size());
            s.Bnnzs = cl::Buffer(*context, CL_MEM_READ_WRITE, sizeof(double) * wellContribs.h_Bnnzs_ocl.size());
            s.Ccols = cl::Buffer(*context, CL_MEM_READ_WRITE, sizeof(int) * wellContribs.h_Ccols_ocl.size());
            s.Bcols = cl::Buffer(*context, CL_MEM_READ_WRITE, sizeof(int) * wellContribs.h_Bcols_ocl.size());
            s.val_pointers = cl::Buffer(*context, CL_MEM_READ_WRITE, sizeof(unsigned int) * wellContribs.h_val_pointers_ocl.size());
            s.toOrder = cl::Buffer(*context, CL_MEM_READ_WRITE, sizeof(int) * Nb);
        }
    }

    void WellContributionsOCLContainer::reinit(Opm::WellContributions &wellContribs){
        num_blocks = wellContribs.h_Ccols_ocl.size();
        num_std_wells = wellContribs.h_val_pointers_ocl.size() - 1;

        s.Cnnzs = cl::Buffer(*context, CL_MEM_READ_WRITE, sizeof(double) * wellContribs.h_Cnnzs_ocl.size());
        s.Dnnzs = cl::Buffer(*context, CL_MEM_READ_WRITE, sizeof(double) * wellContribs.h_Dnnzs_ocl.size());
        s.Bnnzs = cl::Buffer(*context, CL_MEM_READ_WRITE, sizeof(double) * wellContribs.h_Bnnzs_ocl.size());
        s.Ccols = cl::Buffer(*context, CL_MEM_READ_WRITE, sizeof(int) * wellContribs.h_Ccols_ocl.size());
        s.Bcols = cl::Buffer(*context, CL_MEM_READ_WRITE, sizeof(int) * wellContribs.h_Bcols_ocl.size());
        s.val_pointers = cl::Buffer(*context, CL_MEM_READ_WRITE, sizeof(unsigned int) * wellContribs.h_val_pointers_ocl.size());
    }

    void WellContributionsOCLContainer::copy_to_gpu(Opm::WellContributions &wellContribs, int *toOrder_){
        toOrder.insert(toOrder.end(), toOrder_, toOrder_ + Nb);

        if(num_std_wells > 0){
            cl::Event event;
            std::vector<cl::Event> events(7);
            queue->enqueueWriteBuffer(s.Cnnzs, CL_FALSE, 0, sizeof(double) * wellContribs.h_Cnnzs_ocl.size(), wellContribs.h_Cnnzs_ocl.data(), nullptr, &events[0]);
            queue->enqueueWriteBuffer(s.Dnnzs, CL_FALSE, 0, sizeof(double) * wellContribs.h_Dnnzs_ocl.size(), wellContribs.h_Dnnzs_ocl.data(), nullptr, &events[1]);
            queue->enqueueWriteBuffer(s.Bnnzs, CL_FALSE, 0, sizeof(double) * wellContribs.h_Bnnzs_ocl.size(), wellContribs.h_Bnnzs_ocl.data(), nullptr, &events[2]);
            queue->enqueueWriteBuffer(s.Ccols, CL_FALSE, 0, sizeof(int) * wellContribs.h_Ccols_ocl.size(), wellContribs.h_Ccols_ocl.data(), nullptr, &events[3]);
            queue->enqueueWriteBuffer(s.Bcols, CL_FALSE, 0, sizeof(int) * wellContribs.h_Bcols_ocl.size(), wellContribs.h_Bcols_ocl.data(), nullptr, &events[4]);
            queue->enqueueWriteBuffer(s.val_pointers, CL_FALSE, 0, sizeof(unsigned int) * wellContribs.h_val_pointers_ocl.size(), wellContribs.h_val_pointers_ocl.data(), nullptr, &events[5]);
            queue->enqueueWriteBuffer(s.toOrder, CL_FALSE, 0, sizeof(int) * toOrder.size(), toOrder.data(), nullptr, &events[6]);
            event.waitForEvents(events);
        }

        if(!wellContribs.multisegments.empty()){
            multisegments = std::move(std::make_unique<mswVecT>(wellContribs.multisegments));
            num_ms_wells = multisegments->size();
            x_msw.reserve(N);
            y_msw.reserve(N);
        }
    }

    void WellContributionsOCLContainer::update_on_gpu(Opm::WellContributions &wellContribs){
        if(num_std_wells > 0){
            if(num_std_wells != wellContribs.h_val_pointers_ocl.size() || num_blocks != wellContribs.h_Ccols_ocl.size()){
                reinit(wellContribs);
            }

            cl::Event event;
            std::vector<cl::Event> events(6);
            queue->enqueueWriteBuffer(s.Cnnzs, CL_FALSE, 0, sizeof(double) * wellContribs.h_Cnnzs_ocl.size(), wellContribs.h_Cnnzs_ocl.data(), nullptr, &events[0]);
            queue->enqueueWriteBuffer(s.Dnnzs, CL_FALSE, 0, sizeof(double) * wellContribs.h_Dnnzs_ocl.size(), wellContribs.h_Dnnzs_ocl.data(), nullptr, &events[1]);
            queue->enqueueWriteBuffer(s.Bnnzs, CL_FALSE, 0, sizeof(double) * wellContribs.h_Bnnzs_ocl.size(), wellContribs.h_Bnnzs_ocl.data(), nullptr, &events[2]);
            queue->enqueueWriteBuffer(s.Ccols, CL_FALSE, 0, sizeof(int) * wellContribs.h_Ccols_ocl.size(), wellContribs.h_Ccols_ocl.data(), nullptr, &events[3]);
            queue->enqueueWriteBuffer(s.Bcols, CL_FALSE, 0, sizeof(int) * wellContribs.h_Bcols_ocl.size(), wellContribs.h_Bcols_ocl.data(), nullptr, &events[4]);
            queue->enqueueWriteBuffer(s.val_pointers, CL_FALSE, 0, sizeof(unsigned int) * wellContribs.h_val_pointers_ocl.size(), wellContribs.h_val_pointers_ocl.data(), nullptr, &events[5]);
            event.waitForEvents(events);
        }

        if(!wellContribs.multisegments.empty()){
            multisegments = std::move(std::make_unique<mswVecT>(wellContribs.multisegments));
        }
    }

    void WellContributionsOCLContainer::setOpenCLContext(cl::Context *context_){
        this->context = context_;
    }

    void WellContributionsOCLContainer::setOpenCLQueue(cl::CommandQueue *queue_){
        this->queue = queue_;
    }

    void WellContributionsOCLContainer::setKernel(cl::make_kernel<cl::Buffer&, cl::Buffer&, cl::Buffer&, cl::Buffer&,
                                                                  cl::Buffer&, cl::Buffer&, cl::Buffer&, cl::Buffer&,
                                                                  const unsigned int, const unsigned int, cl::Buffer&,
                                                                  cl::LocalSpaceArg, cl::LocalSpaceArg, cl::LocalSpaceArg> *stdwell_apply_){
        this->stdwell_apply = stdwell_apply_;
    }

    void WellContributionsOCLContainer::applyStdWells(cl::Buffer& x, cl::Buffer& y){
        const unsigned int work_group_size = 32;
        const unsigned int total_work_items = num_std_wells * work_group_size;
        const unsigned int lmem1 = sizeof(double) * work_group_size;
        const unsigned int lmem2 = sizeof(double) * dim_wells;

        cl::Event event;
        event = (*stdwell_apply)(cl::EnqueueArgs(*queue, cl::NDRange(total_work_items), cl::NDRange(work_group_size)),
                                 s.Cnnzs, s.Dnnzs, s.Bnnzs, s.Ccols, s.Bcols, x, y, s.toOrder, dim, dim_wells, s.val_pointers,
                                 cl::Local(lmem1), cl::Local(lmem2), cl::Local(lmem2));
    }


    void WellContributionsOCLContainer::applyMSWells(cl::Buffer& x, cl::Buffer& y) {
        cl::Event event;
        std::vector<cl::Event> events(2);

        // copy vectors x and y from GPU to CPU
        queue->enqueueReadBuffer(x, CL_FALSE, 0, sizeof(double) * N, x_msw.data(), nullptr, &events[0]);
        queue->enqueueReadBuffer(y, CL_FALSE, 0, sizeof(double) * N, y_msw.data(), nullptr, &events[1]);
        event.waitForEvents(events);

        // actually apply MultisegmentWells
        for(auto well: *multisegments){
            well->setReordering(toOrder.data(), reorder);
            well->apply(x_msw.data(), y_msw.data());
        }


        // copy vector y from CPU to GPU
        queue->enqueueWriteBuffer(y, CL_FALSE, 0, sizeof(double) * N, y_msw.data(), nullptr, &event);
        event.wait();
    }

    void WellContributionsOCLContainer::apply(cl::Buffer& x, cl::Buffer& y){
        if(num_std_wells > 0){
            applyStdWells(x, y);
        }

        if(num_ms_wells > 0){
            applyMSWells(x, y);
        }
    }
} // end namespace bda
