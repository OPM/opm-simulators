/*
  Copyright 2022-2023 SINTEF AS

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

#ifndef OPM_CUISTL_SET_DEVICE_HEADER
#define OPM_CUISTL_SET_DEVICE_HEADER

namespace Opm::cuistl
{
//! @brief Sets the correct CUDA device in the setting of MPI
//!
//! @note This assumes that every node has equally many GPUs, all of the same caliber
//!
//! @note This probably needs to be called *before* MPI_Init if one uses GPUDirect transfers (see eg.
//! https://devtalk.nvidia.com/default/topic/752046/teaching-and-curriculum-support/multi-gpu-system-running-mpi-cuda-/
//! )
//!
//! @note If no CUDA device is present, this does nothing.
void setDevice(int mpiRank, int numberOfMpiRanks);
} // namespace Opm::cuistl
#endif
