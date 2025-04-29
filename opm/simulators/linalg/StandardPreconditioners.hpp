/*
  Copyright 2025 Equinor ASA

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

#ifndef OPM_STANDARDPRECONDITIONERS_HEADER
#define OPM_STANDARDPRECONDITIONERS_HEADER

// This is a convenience header to include all standard preconditioners.
// It includes the serial and MPI versions of the standard preconditioners.

// Note that you probably should not include this header directly, but rather
// include the PreconditionerFactory.hpp header, which will handle this.

#include <opm/simulators/linalg/is_gpu_operator.hpp>
#include <opm/simulators/linalg/StandardPreconditioners_mpi.hpp>
#include <opm/simulators/linalg/StandardPreconditioners_serial.hpp>

#if HAVE_CUDA
#include <opm/simulators/linalg/StandardPreconditioners_gpu_mpi.hpp>
#include <opm/simulators/linalg/StandardPreconditioners_gpu_serial.hpp>
#endif // HAVE_CUDA

#endif // OPM_STANDARDPRECONDITIONERS_HEADER
