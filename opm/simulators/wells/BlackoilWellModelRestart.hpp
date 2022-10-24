/*
  Copyright 2016 SINTEF ICT, Applied Mathematics.
  Copyright 2016 - 2017 Statoil ASA.
  Copyright 2017 Dr. Blatt - HPC-Simulation-Software & Services
  Copyright 2016 - 2018 IRIS AS

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

#ifndef OPM_BLACKOILWELLMODEL_RESTART_HEADER_INCLUDED
#define OPM_BLACKOILWELLMODEL_RESTART_HEADER_INCLUDED

#include <opm/output/data/Wells.hpp>

#include <vector>

namespace Opm {

class BlackoilWellModelGeneric;
struct PerforationData;
class SingleWellState;

/// Class for restarting the blackoil well model.
class BlackoilWellModelRestart
{
public:
    //! \brief Constructor initializes reference to the well model.
    BlackoilWellModelRestart(const BlackoilWellModelGeneric& wellModel)
        : wellModel_(wellModel)
    {}

    //! \brief Loads per-connection data from restart structures.
    void loadRestartConnectionData(const std::vector<data::Rates::opt>& phs,
                                   const data::Well&                    rst_well,
                                   const std::vector<PerforationData>&  old_perf_data,
                                   SingleWellState&                     ws) const;

    //! \brief Loads per-segment data from restart structures.
    void loadRestartSegmentData(const std::string&                   well_name,
                                const std::vector<data::Rates::opt>& phs,
                                const data::Well&                    rst_well,
                                SingleWellState&                     ws) const;

private:
    const BlackoilWellModelGeneric& wellModel_; //!< Reference to well model
};


} // namespace Opm

#endif
