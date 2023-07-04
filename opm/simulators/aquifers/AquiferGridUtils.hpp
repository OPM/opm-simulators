/*
  Copyright 2013, 2015 SINTEF ICT, Applied Mathematics.
  Copyright 2014, 2015 Dr. Blatt - HPC-Simulation-Software & Services
  Copyright 2014, 2015 Statoil ASA.
  Copyright 2015 NTNU
  Copyright 2015, 2016, 2017 IRIS AS

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

#ifndef AQUIFER_GRID_UTILS_HEADER_INCLUDED
#define AQUIFER_GRID_UTILS_HEADER_INCLUDED

#include <opm/grid/CpGrid.hpp>

#include <algorithm>
#include <vector>

namespace Opm {

template<class Grid>
struct IsNumericalAquiferCell {
    IsNumericalAquiferCell(const Grid&)
    {}

    template<class T>
    bool operator()(const T&) const { return false; }
};

template<>
struct IsNumericalAquiferCell<Dune::CpGrid> {
    IsNumericalAquiferCell(const Dune::CpGrid& grid)
        : grid_(grid)
    {}

    template<class T>
    bool operator()(const T& elem) const
    {
        const auto& aquiferCells = grid_.sortedNumAquiferCells();
        if (aquiferCells.empty())
        {
          return false;
        }
        auto candidate = std::lower_bound(aquiferCells.begin(),
                                          aquiferCells.end(), elem.index());
        return candidate != aquiferCells.end() && *candidate == elem.index();
    }

private:
    const Dune::CpGrid& grid_;
};

} // namespace Opm

#endif
