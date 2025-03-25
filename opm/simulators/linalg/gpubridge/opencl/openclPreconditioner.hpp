/*
  Copyright 2021 Equinor ASA

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

#ifndef OPM_OPENCLPRECONDITIONER_HEADER_INCLUDED
#define OPM_OPENCLPRECONDITIONER_HEADER_INCLUDED

#include <opm/simulators/linalg/gpubridge/opencl/opencl.hpp>
#include <opm/simulators/linalg/gpubridge/Preconditioner.hpp>

namespace Opm::Accelerator {

template<class Scalar> class BlockedMatrix;

template <class Scalar, unsigned int block_size>
class openclPreconditioner : public Preconditioner<Scalar, block_size, cl::Buffer>
{

protected:
    std::shared_ptr<cl::Context> context;
    std::shared_ptr<cl::CommandQueue> queue;
    std::vector<cl::Event> events;
    cl_int err;

    openclPreconditioner(int verbosity_)
        : Preconditioner<Scalar, block_size, cl::Buffer>(verbosity_)
    {}

public:
    static std::unique_ptr<openclPreconditioner<Scalar, block_size>> create(PreconditionerType type, int verbosity, bool opencl_ilu_parallel);

    // nested Preconditioners might need to override this
    virtual void setOpencl(std::shared_ptr<cl::Context>& context, std::shared_ptr<cl::CommandQueue>& queue);
};
} //namespace Opm

#endif
