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
    class BlackoilPropertiesInterface;

    /// \file
    ///
    /// Functions for initializing a reservoir state.

    /// Initialize a two-phase state from parameters.
    /// The following parameters are accepted (defaults):
    ///    - convection_testcase   (false) -- Water in the 'left' part of the grid.
    ///    - ref_pressure          (100)   -- Initial pressure in bar for all cells
    ///                                     (if convection_testcase is true),
    ///                                     or pressure at woc depth.
    ///    - segregation_testcase  (false) -- Water above the woc instead of below.
    ///    - water_oil_contact     (none)  -- Depth of water-oil contact (woc).
    ///    - init_saturation       (none)  -- Initial water saturation for all cells.
    ///
    /// If convection_testcase is true, the saturation is initialised
    /// as indicated, and pressure is initialised to a constant value
    /// ('ref_pressure').
    /// If segregation_testcase is true, the saturation is initialised
    /// as indicated, and pressure is initialised hydrostatically.
    /// Otherwise we have 3 cases:
    ///   1. If 'water_oil_contact' is given, saturation is initialised
    ///      accordingly.
    ///   2. If 'water_oil_contact' is not given, but 'init_saturation'
    ///      is given, water saturation is set to that value everywhere.
    ///   3. If neither are given, water saturation is set to minimum.
    ///
    /// In all three cases, pressure is initialised hydrostatically.
    /// In case 2) and 3), the depth of the first cell is used as reference depth.
    template <class State>
    void initStateBasic(const UnstructuredGrid& grid,
                        const IncompPropertiesInterface& props,
                        const parameter::ParameterGroup& param,
                        const double gravity,
                        State& state);

    /// Initialize a blackoil state from parameters.
    /// The following parameters are accepted (defaults):
    ///    - convection_testcase   (false) -- Water in the 'left' part of the grid.
    ///    - ref_pressure          (100)   -- Initial pressure in bar for all cells
    ///                                     (if convection_testcase is true),
    ///                                     or pressure at woc depth.
    ///    - water_oil_contact     (none)  -- Depth of water-oil contact (woc).
    /// If convection_testcase is true, the saturation is initialised
    /// as indicated, and pressure is initialised to a constant value
    /// ('ref_pressure').
    /// Otherwise we have 2 cases:
    ///   1. If 'water_oil_contact' is given, saturation is initialised
    ///      accordingly.
    ///   2. Water saturation is set to minimum.
    /// In both cases, pressure is initialised hydrostatically.
    /// In case 2., the depth of the first cell is used as reference depth.
    template <class State>
    void initStateBasic(const UnstructuredGrid& grid,
                        const BlackoilPropertiesInterface& props,
                        const parameter::ParameterGroup& param,
                        const double gravity,
                        State& state);

    /// Initialize a two-phase state from input deck.
    /// If EQUIL is present:
    ///   - saturation is set according to the water-oil contact,
    ///   - pressure is set to hydrostatic equilibrium.
    /// Otherwise:
    ///   - saturation is set according to SWAT,
    ///   - pressure is set according to PRESSURE.
    template <class Props, class State>
    void initStateFromDeck(const UnstructuredGrid& grid,
                           const Props& props,
                           const EclipseGridParser& deck,
                           const double gravity,
                           State& state);

    /// Initialize a two-phase water-oil blackoil state from input deck.
    /// If EQUIL is present:
    ///   - saturation is set according to the water-oil contact,
    ///   - pressure is set to hydrostatic equilibrium.
    /// Otherwise:
    ///   - saturation is set according to SWAT,
    ///   - pressure is set according to PRESSURE.
    /// In addition, this function sets surfacevol.
    template <class Props, class State>
    void initBlackoilStateFromDeck(const UnstructuredGrid& grid,
                                   const Props& props,
                                   const EclipseGridParser& deck,
                                   const double gravity,
                                   State& state);


} // namespace Opm

#include <opm/core/simulator/initState_impl.hpp>

#endif // OPM_INITSTATE_HEADER_INCLUDED
