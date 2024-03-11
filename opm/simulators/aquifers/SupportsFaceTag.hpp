/*
  File adapted from BlackoilWellModel.hpp

  Copyright 2017 TNO - Heat Transfer & Fluid Dynamics, Modelling & Optimization of the Subsurface
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


#ifndef OPM_SUPPORTS_FACETAG_HEADER_INCLUDED
#define OPM_SUPPORTS_FACETAG_HEADER_INCLUDED

namespace Dune { class CpGrid; }

namespace Opm {

template<class Grid>
class SupportsFaceTag
    : public std::bool_constant<false>
{};


template<>
class SupportsFaceTag<Dune::CpGrid>
    : public std::bool_constant<true>
{};


} // namespace Opm

#endif // OPM_SUPPORT_FACETAG_HEADER_INCLUDED
