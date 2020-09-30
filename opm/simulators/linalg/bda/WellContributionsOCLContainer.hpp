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

#ifndef WELLCONTRIBUTIONSOCLCONTAINER_HPP
#define WELLCONTRIBUTIONSOCLCONTAINER_HPP

#include <opm/simulators/linalg/bda/opencl.hpp>
#include <opm/simulators/linalg/bda/WellContributions.hpp>
#include <opm/simulators/linalg/bda/MultisegmentWellContribution.hpp>

namespace bda
{
    class WellContributionsOCLContainer
    {
    private:
        typedef std::vector<Opm::MultisegmentWellContribution*> mswVecT;

        unsigned int dim, dim_wells;
        unsigned int num_blocks = 0;
        unsigned int num_std_wells = 0;
        unsigned int num_ms_wells = 0;           // number of MultisegmentWells in this object, must equal multisegments.size()
        int N, Nb;
        std::vector<int> toOrder;
        std::vector<double> x_msw, y_msw;
        std::unique_ptr<mswVecT> multisegments;

        typedef struct {
            cl::Buffer Cnnzs, Dnnzs, Bnnzs;
            cl::Buffer Ccols, Bcols;
            cl::Buffer val_pointers, toOrder;
        } GPU_storage;

        GPU_storage s;
        cl::Context *context;
        cl::CommandQueue *queue;
        cl::make_kernel<cl::Buffer&, cl::Buffer&, cl::Buffer&, cl::Buffer&,
                        cl::Buffer&, cl::Buffer&, cl::Buffer&, cl::Buffer&,
                        const unsigned int, const unsigned int, cl::Buffer&,
                        cl::LocalSpaceArg, cl::LocalSpaceArg, cl::LocalSpaceArg> *stdwell_apply;

        void reinit(Opm::WellContributions &wellContribs);
        void applyStdWells(cl::Buffer& x, cl::Buffer& y);
        void applyMSWells(cl::Buffer& x, cl::Buffer& y);

    public:
        WellContributionsOCLContainer();
        ~WellContributionsOCLContainer();

        void apply(cl::Buffer& x, cl::Buffer& y);
        void init(Opm::WellContributions &wellContribs, int N, int Nb);
        void copy_to_gpu(Opm::WellContributions &wellContribs, int *toOrder_);
        void update_on_gpu(Opm::WellContributions &wellContribs);
        void setOpenCLContext(cl::Context *context);
        void setOpenCLQueue(cl::CommandQueue *queue);
        void setKernel(cl::make_kernel<cl::Buffer&, cl::Buffer&, cl::Buffer&, cl::Buffer&,
                                       cl::Buffer&, cl::Buffer&, cl::Buffer&, cl::Buffer&,
                                       const unsigned int, const unsigned int, cl::Buffer&,
                                       cl::LocalSpaceArg, cl::LocalSpaceArg, cl::LocalSpaceArg> *stdwell_apply_);
    };
} // end namespace bda

#endif
