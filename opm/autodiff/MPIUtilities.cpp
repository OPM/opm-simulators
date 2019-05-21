/*
  Copyright 2014 SINTEF ICT, Applied Mathematics.
  Copyright 2014 IRIS AS

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
#include <boost/any.hpp>
#include <opm/simulators/linalg/ParallelIstlInformation.hpp>

namespace Opm
{
    /// Return true if this is a serial run, or rank zero on an MPI run.
    bool isIORank(const boost::any& parallel_info)
    {
#if HAVE_MPI
        if (parallel_info.type() == typeid(ParallelISTLInformation)) {
            const ParallelISTLInformation& info =
                boost::any_cast<const ParallelISTLInformation&>(parallel_info);
            return info.communicator().rank() == 0;
        } else {
            return true;
        }
#else
        static_cast<void>(parallel_info); // Suppress unused argument warning.
        return true;
#endif
    }
} // namespace Opm

