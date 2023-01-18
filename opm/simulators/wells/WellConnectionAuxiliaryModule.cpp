/*
  Copyright 2017 Dr. Blatt - HPC-Simulation-Software & Services
  Copyright 2017 Statoil ASA.

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
#include <opm/simulators/wells/WellConnectionAuxiliaryModule.hpp>

#include <opm/grid/CpGrid.hpp>

#include <opm/input/eclipse/Schedule/Schedule.hpp>
#include <opm/input/eclipse/Schedule/Well/Well.hpp>
#include <opm/input/eclipse/Schedule/Well/WellConnections.hpp>

#include <algorithm>

namespace Opm
{

WellConnectionAuxiliaryModuleGeneric::
WellConnectionAuxiliaryModuleGeneric(const Schedule& schedule,
                                     const Dune::CpGrid& grid)
{
    // Create cartesian to compressed mapping
    const auto& globalCell = grid.globalCell();
    const auto& cartesianSize = grid.logicalCartesianSize();

    auto size = cartesianSize[0] * cartesianSize[1] * cartesianSize[2];

    std::vector<int> cartesianToCompressed(size, -1);
    auto begin = globalCell.begin();

    for (auto cell = begin, end = globalCell.end(); cell != end; ++cell)
    {
      cartesianToCompressed[ *cell ] = cell - begin;
    }

    const auto& schedule_wells = schedule.getWellsatEnd();
    wells_.reserve(schedule_wells.size());

    // initialize the additional cell connections introduced by wells.
    for (const auto& well : schedule_wells)
    {
        std::vector<int> compressed_well_perforations;
        // All possible completions of the well
        const auto& completionSet = well.getConnections();
        compressed_well_perforations.reserve(completionSet.size());

        for (const auto& completion : completionSet)
        {
            int compressed_idx = cartesianToCompressed[completion.global_index()];
            if (compressed_idx >= 0) // Ignore completions in inactive/remote cells.
            {
                compressed_well_perforations.push_back(compressed_idx);
            }
        }

        if (!compressed_well_perforations.empty())
        {
            std::sort(compressed_well_perforations.begin(),
                      compressed_well_perforations.end());

            wells_.push_back(compressed_well_perforations);
        }
    }
}

} // end namespace Opm
