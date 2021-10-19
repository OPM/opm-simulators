/*
  Copyright 2019 SINTEF Digital, Mathematics and Cybernetics.
  Copyright 2019 Equinor.

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

#ifndef OPM_GATHERDEFERREDLOGGER_HEADER_INCLUDED
#define OPM_GATHERDEFERREDLOGGER_HEADER_INCLUDED

#include <opm/simulators/utils/DeferredLogger.hpp>

namespace Opm
{

    /// Create a global log combining local logs
    Opm::DeferredLogger gatherDeferredLogger(const Opm::DeferredLogger& local_deferredlogger, Parallel::Communication communicator);

} // namespace Opm


#endif // OPM_GATHERDEFERREDLOGGER_HEADER_INCLUDED
