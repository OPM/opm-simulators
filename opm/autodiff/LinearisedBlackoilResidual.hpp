/*
  Copyright 2014 SINTEF ICT, Applied Mathematics.

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

#ifndef OPM_LINEARISEDBLACKOILRESIDUAL_HEADER_INCLUDED
#define OPM_LINEARISEDBLACKOILRESIDUAL_HEADER_INCLUDED

#include <opm/autodiff/AutoDiffBlock.hpp>

namespace Opm
{

    /// Residual structure of the fully implicit solver.
    /// All equations are given as AD types, with multiple
    /// jacobian blocks corresponding to the primary unknowns. The
    /// primary unknowns are for a three-phase simulation, in order:
    ///    p    (pressure)
    ///    sw   (water saturation)
    ///    xvar (gas saturation, gas-oil ratio or oil-gas ratio)
    ///    qs   (well outflows by well and phase)
    ///    bhp  (bottom hole pressures)
    /// In the above, the xvar variable will have a different
    /// meaning from cell to cell, corresponding to the state in
    /// that cell (saturated, undersaturated oil or undersaturated
    /// gas). In a two-phase simulation, either sw or xvar is not
    /// used, depending on which phase is missing.
    ///
    /// Note: this class is strongly coupled to the class
    /// FullyImplicitBlackoilSolver, and is separated from that
    /// class to facilitate the development of linear solver
    /// strategies outside that class.
    struct LinearisedBlackoilResidual {
        /// A type alias for the automatic differentiation type.
        typedef AutoDiffBlock<double> ADB;
        /// The material_balance_eq vector has one element for each
        /// active phase, each of which has size equal to the number
        /// of cells. Each material balance equation is given in terms
        /// of surface volumes (in SI units, that is standard m^3).
        std::vector<ADB> material_balance_eq;
        /// The well_flux_eq has size equal to the number of wells
        /// times the number of phases. It contains the well flow
        /// equations, relating the total well flows to
        /// bottom-hole pressures and reservoir conditions.
        ADB well_flux_eq;
        /// The well_eq has size equal to the number of wells. It
        /// contains the well control equations, that is for each
        /// well either a rate specification or bottom hole
        /// pressure specification.
        ADB well_eq;

        std::vector<double> matbalscale;

        /// The size of the non-linear system.
        int sizeNonLinear() const;
    };

} // namespace Opm


#endif // OPM_LINEARISEDBLACKOILRESIDUAL_HEADER_INCLUDED
