/*
  Copyright 2013, 2015, 2020 SINTEF Digital, Mathematics and Cybernetics.
  Copyright 2015 Andreas Lauser
  Copyright 2017 IRIS

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

#ifndef OPM_SIMULATOR_REPORT_BANNERS_HEADER_INCLUDED
#define OPM_SIMULATOR_REPORT_BANNERS_HEADER_INCLUDED

namespace Opm { class SimulatorTimer; }

namespace Opm::details {

void outputReportStep(const SimulatorTimer& timer);

} // namespace Opm::details

#endif // OPM_SIMULATOR_REPORT_BANNERS_HEADER_INCLUDED
