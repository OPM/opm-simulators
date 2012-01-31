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


#include <opm/core/eclipse/EclipseGridParser.hpp>
#include <tr1/array>

namespace Opm
{

    /// Class collecting pvt properties for 2 phases, reading from
    /// eclipse input (keywords DENSITY, PVTW, PVCDO).
    ///
    /// All phases are incompressible and have constant viscosities.
    /// For all the methods, the following apply: p and z are unused.
    /// Output arrays shall be of size n*numPhases(), and must be valid
    /// before calling the method.
    /// NOTE: This class is intentionally similar to BlackoilPvtProperties.
    class PvtPropertiesIncompFromDeck
    {
    public:
        /// Default constructor.
        PvtPropertiesIncompFromDeck();

        /// Initialize from deck.
	void init(const EclipseGridParser& deck);

        /// Number of active phases.
        int numPhases() const;

        /// Densities of stock components at surface conditions.
        /// \return  Array of size numPhases().
	const double* surfaceDensities() const;

        /// Viscosities.
        const double* viscosity() const;

    private:
	std::tr1::array<double, 2> density_;
	std::tr1::array<double, 2> viscosity_;
    };

}



#endif // OPM_PVTPROPERTIESINCOMPFROMDECK_HEADER_INCLUDED
