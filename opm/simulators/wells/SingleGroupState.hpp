/*
  Copyright 2012, 2014 SINTEF Digital, Mathematics and Cybernetics.
  Copyright 2017, 2019 NORCE.
  Copyright 2020 Equinor ASA.

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

#ifndef OPM_SINGLEGROUPSTATE_HEADER_INCLUDED
#define OPM_SINGLEGROUPSTATE_HEADER_INCLUDED

#include <opm/parser/eclipse/EclipseState/Schedule/Group/Group.hpp>

#include <array>

namespace Opm
{

/// Encapsulate all data needed to represent the state of a well group
/// for persistence purposes, i.e. from one timestep to the next or
/// for restarting, and for accumulating rates to get correct group
/// rate limits.
template <int NumActivePhases>
struct SingleGroupState
{
    // ------ Types ------

    using PhaseRates = std::array<double, NumActivePhases>;

    // ------ Data members ------

    // Flags and statuses.
    Group::InjectionCMode current_injection_control = Group::InjectionCMode::NONE;
    Group::ProductionCMode current_production_control = Group::ProductionCMode::NONE;

    // Quantities.
    PhaseRates production_reduction_surface_rates = {0.0};
    PhaseRates injection_reduction_surface_rates = {0.0};
    PhaseRates injection_potentials = {0.0};
    PhaseRates injection_vrep_rates = {0.0};
    PhaseRates injection_rein_rates = {0.0};
};


} // namespace Opm


#endif // OPM_SINGLEGROUPSTATE_HEADER_INCLUDED
