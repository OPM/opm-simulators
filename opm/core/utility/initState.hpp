/*
  Copyright 2012 SINTEF ICT, Applied Mathematics.

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

#ifndef OPM_INITSTATE_HEADER_INCLUDED
#define OPM_INITSTATE_HEADER_INCLUDED

struct UnstructuredGrid;

namespace Opm
{

    namespace parameter { class ParameterGroup; }
    class EclipseGridParser;
    class IncompPropertiesInterface;

    /// Initialize a state from parameters.
    /// The following parameters are accepted (defaults):
    ///    num_phases         (2)         Must be 1 or 2.
    ///    relperm_func       ("Linear")  Must be "Constant", "Linear" or "Quadratic".
    ///    rho1 [rho2, rho3]  (1.0e3)     Density in kg/m^3
    ///    mu1 [mu2, mu3]     (1.0)       Viscosity in cP
    ///    porosity           (1.0)       Porosity
    ///    permeability       (100.0)     Permeability in mD
    template <class State>
    void initStateTwophaseBasic(const UnstructuredGrid& grid,
                                const IncompPropertiesInterface& props,
                                const parameter::ParameterGroup& param,
                                const double gravity,
                                State& state);

    /// Initialize a state from input deck.
    /// If EQUIL is present:
    ///   - saturation is set according to the water-oil contact,
    ///   - pressure is set to hydrostatic equilibrium.
    /// Otherwise:
    ///   - saturation is set according to SWAT,
    ///   - pressure is set according to PRESSURE.
    template <class State>
    void initStateTwophaseFromDeck(const UnstructuredGrid& grid,
                                   const IncompPropertiesInterface& props,
                                   const EclipseGridParser& deck,
                                   const double gravity,
                                   State& state);

} // namespace Opm

#include <opm/core/utility/initState_impl.hpp>

#endif // OPM_INITSTATE_HEADER_INCLUDED
