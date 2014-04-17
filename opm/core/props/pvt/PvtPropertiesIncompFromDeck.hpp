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

#ifndef OPM_PVTPROPERTIESINCOMPFROMDECK_HEADER_INCLUDED
#define OPM_PVTPROPERTIESINCOMPFROMDECK_HEADER_INCLUDED

#include <opm/core/io/eclipse/EclipseGridParser.hpp>
#include <opm/parser/eclipse/Deck/Deck.hpp>
#include <array>

namespace Opm
{

    /// Class collecting pvt properties for 2 phases, reading from
    /// eclipse input (keywords DENSITY, PVTW, PVCDO).
    ///
    /// All phases are incompressible and have constant viscosities.
    /// NOTE: This class is intentionally similar to BlackoilPvtProperties.
    class PvtPropertiesIncompFromDeck
    {
    public:
        /// Default constructor.
        PvtPropertiesIncompFromDeck();

        /// Initialize from deck.
        void init(const EclipseGridParser& deck);

        /// Initialize from deck.
        void init(Opm::DeckConstPtr newParserDeck);

        /// Number of active phases.
        int numPhases() const;

        /// Densities of stock components at surface conditions.
        /// \return  Array of size numPhases().
        const double* surfaceDensities() const;

        /// Densities of stock components at reservoir conditions.
        /// Note: a reasonable question to ask is why there can be
        /// different densities at surface and reservoir conditions,
        /// when the phases are assumed incompressible. The answer is
        /// that even if we approximate the phases as being
        /// incompressible during simulation, the density difference
        /// between surface and reservoir may be larger. For accurate
        /// reporting and using data given in terms of surface values,
        /// we need to handle this difference.
        /// \return  Array of size numPhases().
        const double* reservoirDensities() const;

        /// Viscosities.
        const double* viscosity() const;

    private:
        std::array<double, 2> surface_density_;
        std::array<double, 2> reservoir_density_;
        std::array<double, 2> viscosity_;
    };

}



#endif // OPM_PVTPROPERTIESINCOMPFROMDECK_HEADER_INCLUDED
