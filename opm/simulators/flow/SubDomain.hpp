/*
  Copyright 2021 Total SE

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

#ifndef OPM_SUBDOMAIN_HEADER_INCLUDED
#define OPM_SUBDOMAIN_HEADER_INCLUDED

#include <opm/grid/common/SubGridPart.hpp>

#include <vector>

namespace Opm
{

    template <class Grid>
    struct SubDomain
    {
        int index;
        std::vector<int> cells;
        std::vector<bool> interior;
        Dune::SubGridPart<Grid> view;
        SubDomain(const int i, std::vector<int>&& c, std::vector<bool>&& in, Dune::SubGridPart<Grid>&& v)
            : index(i), cells(std::move(c)), interior(std::move(in)), view(std::move(v))
        {}
    };

} // namespace Opm


#endif // OPM_SUBDOMAIN_HEADER_INCLUDED
