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

#include <opm/simulators/linalg/bda/opencl/opencl.hpp>
#include <opm/simulators/linalg/bda/Preconditioner.hpp>

namespace Opm::Accelerator {

template<class Scalar> class BlockedMatrix;

template <class Scalar, unsigned int block_size>
class openclPreconditioner : public Preconditioner<Scalar, block_size>
{

protected:
    std::shared_ptr<cl::Context> context;
    std::shared_ptr<cl::CommandQueue> queue;
    std::vector<cl::Event> events;
    cl_int err;

    openclPreconditioner(int verbosity_) :
    Preconditioner<Scalar, block_size>(verbosity_)
    {};

public:
    virtual ~openclPreconditioner() = default;

    static std::unique_ptr<openclPreconditioner<Scalar, block_size>> create(PreconditionerType type, int verbosity, bool opencl_ilu_parallel);

    // nested Preconditioners might need to override this
    virtual void setOpencl(std::shared_ptr<cl::Context>& context, std::shared_ptr<cl::CommandQueue>& queue);

    // apply preconditioner, x = prec(y)
    virtual void apply(const cl::Buffer& y, cl::Buffer& x) = 0;
 
    // create/update preconditioner, probably used every linear solve
    // the version with two params can be overloaded, if not, it will default to using the one param version
    virtual bool create_preconditioner(BlockedMatrix<Scalar> *mat) = 0;
    virtual bool create_preconditioner(BlockedMatrix<Scalar> *mat, BlockedMatrix<Scalar> *jacMat) = 0;
};
} //namespace Opm

#endif
