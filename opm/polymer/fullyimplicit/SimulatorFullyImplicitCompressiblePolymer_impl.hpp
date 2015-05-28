/*
  Copyright 2014 SINTEF ICT, Applied Mathematics.
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
#ifndef OPM_SIMULATORFULLYIMPLICITCOMPRESSIBLEPOLYMER_IMPL_HEADER_INCLUDED
#define OPM_SIMULATORFULLYIMPLICITCOMPRESSIBLEPOLYMER_IMPL_HEADER_INCLUDED

namespace Opm
{

/// Class collecting all necessary components for a two-phase simulation.
SimulatorFullyImplicitCompressiblePolymer::
SimulatorFullyImplicitCompressiblePolymer(const parameter::ParameterGroup& param,
                                          const UnstructuredGrid& grid,
                                          const DerivedGeology& geo,
                                          BlackoilPropsAdInterface& props,
                                          const PolymerPropsAd&    polymer_props,
                                          const RockCompressibility* rock_comp_props,
                                          std::shared_ptr<EclipseState> eclipse_state,
                                          BlackoilOutputWriter& output_writer,
                                          Opm::DeckConstPtr& deck,
                                          NewtonIterationBlackoilInterface& linsolver,
                                          const double* gravity)
: BaseType(param,
           grid,
           geo,
           props,
           rock_comp_props,
           linsolver,
           gravity,
           /*disgas=*/false,
           /*vapoil=*/false,
           eclipse_state,
           output_writer,
           /*threshold_pressures_by_face=*/std::vector<double>())
    , deck_(deck)
    , polymer_props_(polymer_props)

{
}

} // namespace Opm

#endif // OPM_SIMULATORFULLYIMPLICITCOMPRESSIBLEPOLYMER_HEADER_INCLUDED
