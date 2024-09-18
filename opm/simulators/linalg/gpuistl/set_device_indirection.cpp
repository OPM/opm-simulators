/*
  Copyright 2024 SINTEF AS

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

#include <opm/simulators/flow/FlowGenericVanguard.hpp>
#if HAVE_CUDA
#include <opm/simulators/linalg/gpuistl/set_device.hpp>
#endif

namespace Opm::gpuistl {

    void printDevice()
    {
#if HAVE_CUDA
#if HAVE_MPI
        Opm::gpuistl::printDevice(FlowGenericVanguard::comm().rank(), FlowGenericVanguard::comm().size());
#else
        Opm::gpuistl::printDevice(0, 1);
#endif
#endif
    }

    void setDevice()
    {
#if HAVE_CUDA
#if HAVE_MPI
        Opm::gpuistl::setDevice(FlowGenericVanguard::comm().rank(), FlowGenericVanguard::comm().size());
#else
        Opm::gpuistl::setDevice(0, 1);
#endif
#endif
    }

} // namespace Opm::gpuistl
