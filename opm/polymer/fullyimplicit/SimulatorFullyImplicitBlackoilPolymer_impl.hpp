/*
  Copyright 2013 SINTEF ICT, Applied Mathematics.
  Copyright 2014 IRIS AS
  Copyright 2014 STATOIL ASA.

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

namespace Opm
{
    SimulatorFullyImplicitBlackoilPolymer::
    SimulatorFullyImplicitBlackoilPolymer(const parameter::ParameterGroup& param,
                                          const Grid& grid,
                                          const DerivedGeology& geo,
                                          BlackoilPropsAdInterface& props,
                                          const PolymerPropsAd& polymer_props,
                                          const RockCompressibility* rock_comp_props,
                                          NewtonIterationBlackoilInterface& linsolver,
                                          const double* gravity,
                                          const bool has_disgas,
                                          const bool has_vapoil,
                                          const bool has_polymer,
                                          std::shared_ptr<EclipseState> eclipse_state,
                                          BlackoilOutputWriter& output_writer,
                                          Opm::DeckConstPtr& deck,
                                          const std::vector<double>& threshold_pressures_by_face)
    : BaseType(param,
               grid,
               geo,
               props,
               rock_comp_props,
               linsolver,
               gravity,
               has_disgas,
               has_vapoil,
               eclipse_state,
               output_writer,
               threshold_pressures_by_face)
        , polymer_props_(polymer_props)
        , has_polymer_(has_polymer)
        , deck_(deck)
    {
    }
} // namespace Opm
