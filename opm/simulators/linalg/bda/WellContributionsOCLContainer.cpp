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
#include<iostream>

namespace bda
{
    using Opm::OpmLog;
    using Dune::Timer;

    void WellContributionsOCLContainer::init(Opm::WellContributions &wellContribs, int Nb_){
        Nb = Nb_;
        dim = wellContribs.dim;
        dim_wells = wellContribs.dim_wells;

        if(!wellContribs.h_val_pointers_ocl.empty()){
            num_std_wells = wellContribs.h_val_pointers_ocl.size() - 1;

            s.Cnnzs = cl::Buffer(*context, CL_MEM_READ_WRITE, sizeof(double) * wellContribs.h_Cnnzs_ocl.size());
            s.Dnnzs = cl::Buffer(*context, CL_MEM_READ_WRITE, sizeof(double) * wellContribs.h_Dnnzs_ocl.size());
            s.Bnnzs = cl::Buffer(*context, CL_MEM_READ_WRITE, sizeof(double) * wellContribs.h_Bnnzs_ocl.size());
            s.Ccols = cl::Buffer(*context, CL_MEM_READ_WRITE, sizeof(int) * wellContribs.h_Ccols_ocl.size());
            s.Bcols = cl::Buffer(*context, CL_MEM_READ_WRITE, sizeof(int) * wellContribs.h_Bcols_ocl.size());
            s.val_pointers = cl::Buffer(*context, CL_MEM_READ_WRITE, sizeof(unsigned int) * wellContribs.h_val_pointers_ocl.size());
            s.toOrder = cl::Buffer(*context, CL_MEM_READ_WRITE, sizeof(int) * Nb);
        }
        else{
            num_std_wells = 0;
        }
    }

    void WellContributionsOCLContainer::copy_to_gpu(Opm::WellContributions &wellContribs){
        if(num_std_wells > 0){
            toOrder.insert(toOrder.end(), wellContribs.toOrder, wellContribs.toOrder + Nb);

            cl::Event event;
            queue->enqueueWriteBuffer(s.Cnnzs, CL_TRUE, 0, sizeof(double) * wellContribs.h_Cnnzs_ocl.size(), wellContribs.h_Cnnzs_ocl.data());
            queue->enqueueWriteBuffer(s.Dnnzs, CL_TRUE, 0, sizeof(double) * wellContribs.h_Dnnzs_ocl.size(), wellContribs.h_Dnnzs_ocl.data());
            queue->enqueueWriteBuffer(s.Bnnzs, CL_TRUE, 0, sizeof(double) * wellContribs.h_Bnnzs_ocl.size(), wellContribs.h_Bnnzs_ocl.data());
            queue->enqueueWriteBuffer(s.Ccols, CL_TRUE, 0, sizeof(int) * wellContribs.h_Ccols_ocl.size(), wellContribs.h_Ccols_ocl.data());
            queue->enqueueWriteBuffer(s.Bcols, CL_TRUE, 0, sizeof(int) * wellContribs.h_Bcols_ocl.size(), wellContribs.h_Bcols_ocl.data());
            queue->enqueueWriteBuffer(s.val_pointers, CL_TRUE, 0, sizeof(unsigned int) * wellContribs.h_val_pointers_ocl.size(), wellContribs.h_val_pointers_ocl.data());
            queue->enqueueWriteBuffer(s.toOrder, CL_TRUE, 0, sizeof(int) * toOrder.size(), toOrder.data(), nullptr, &event);
            event.wait();
        }
    }

    void WellContributionsOCLContainer::update_on_gpu(Opm::WellContributions &wellContribs){
        if(num_std_wells > 0){
            cl::Event event;
            queue->enqueueWriteBuffer(s.Cnnzs, CL_TRUE, 0, sizeof(double) * wellContribs.h_Cnnzs_ocl.size(), wellContribs.h_Cnnzs_ocl.data());
            queue->enqueueWriteBuffer(s.Dnnzs, CL_TRUE, 0, sizeof(double) * wellContribs.h_Dnnzs_ocl.size(), wellContribs.h_Dnnzs_ocl.data());
            queue->enqueueWriteBuffer(s.Bnnzs, CL_TRUE, 0, sizeof(double) * wellContribs.h_Bnnzs_ocl.size(), wellContribs.h_Bnnzs_ocl.data(), nullptr, &event);
            event.wait();
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

    void WellContributionsOCLContainer::apply(cl::Buffer& x, cl::Buffer& y){
        if(num_std_wells > 0){
            applyStdWells(x, y);
        }
    }
} // end namespace bda
