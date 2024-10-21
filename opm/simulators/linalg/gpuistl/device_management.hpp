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

#ifndef OPM_GPUISTL_DEVICE_MANAGEMENT
#define OPM_GPUISTL_DEVICE_MANAGEMENT

/*
  This file should not be hipified, and serves as a layer between main and gpuistl/set_device
  that does not depend on the library such that the simulatorobjects to not depend
  on the library and can be built in parallel.
*/

namespace Opm::gpuistl {
    void printDevice();
    void setDevice();
}

#endif // namespace Opm::gpuistl
