// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*
  This file is part of the Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 2 of the License, or
  (at your option) any later version.

  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.

  Consult the COPYING file in the top-level source directory of this
  module for the precise wording of the license and the list of
  copyright holders.
*/
/*!
 * \file
 *
 * \brief Classes required for dynamic convective mixing.
 */
#ifndef OPM_CONVECTIVEMIXING_MODULE_PARAM_HPP
#define OPM_CONVECTIVEMIXING_MODULE_PARAM_HPP

#include <opm/common/utility/VectorWithDefaultAllocator.hpp>

#if HAVE_CUDA
#if USE_HIP
#include <opm/simulators/linalg/gpuistl_hip/GpuBuffer.hpp>
#include <opm/simulators/linalg/gpuistl_hip/GpuView.hpp>
#else // !USE_HIP
#include <opm/simulators/linalg/gpuistl/GpuBuffer.hpp>
#include <opm/simulators/linalg/gpuistl/GpuView.hpp>
#endif // USE_HIP
#endif // HAVE_CUDA

namespace Opm {

template<class Scalar, template<class> class Storage = VectorWithDefaultAllocator>
struct ConvectiveMixingModuleParam
{
    Storage<bool> active_;
    Storage<Scalar> Xhi_;
    Storage<Scalar> Psi_;
};

#ifdef HAVE_CUDA
namespace gpuistl
{

template <class Scalar>
ConvectiveMixingModuleParam<Scalar, GpuView>
make_view(ConvectiveMixingModuleParam<Scalar, GpuBuffer>& params)
{
    ConvectiveMixingModuleParam<Scalar, GpuView> view;
    view.active_ = gpuistl::make_view(params.active_);
    view.Xhi_ = gpuistl::make_view(params.Xhi_);
    view.Psi_ = gpuistl::make_view(params.Psi_);
    return view;
}

template <class Scalar>
ConvectiveMixingModuleParam<Scalar, GpuBuffer>
copy_to_gpu(const ConvectiveMixingModuleParam<Scalar, VectorWithDefaultAllocator>& params)
{
    return ConvectiveMixingModuleParam<Scalar, GpuBuffer>{
        gpuistl::GpuBuffer(params.active_),
        gpuistl::GpuBuffer(params.Xhi_),
        gpuistl::GpuBuffer(params.Psi_)
    };
}

}
#endif

}

#endif
