/*
  Copyright 2026 Equinor AS

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
#ifndef OPM_GPUISTL_ASSEMBLY_IF_AVAILABLE_HPP
#define OPM_GPUISTL_ASSEMBLY_IF_AVAILABLE_HPP

#include <opm/common/utility/gpuDecorators.hpp>

// Safe to include from CPU code. Backend selection is centralized here;
// GpuAssembly.hpp and its dependencies are generated for HIP by the build.
#if HAVE_CUDA
#if USE_HIP
#include <opm/simulators/linalg/gpuistl_hip/GpuAssembly.hpp>
#else
#include <opm/simulators/linalg/gpuistl/GpuAssembly.hpp>
#endif
#endif

#endif // OPM_GPUISTL_ASSEMBLY_IF_AVAILABLE_HPP
