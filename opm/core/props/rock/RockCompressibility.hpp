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

#ifndef OPM_ROCKCOMPRESSIBILITY_HEADER_INCLUDED
#define OPM_ROCKCOMPRESSIBILITY_HEADER_INCLUDED

#include <opm/parser/eclipse/Deck/Deck.hpp>

#include <vector>

namespace Opm
{

    class EclipseGridParser;
    namespace parameter { class ParameterGroup; }

    class RockCompressibility
    {
    public:
        /// Construct from input deck.
        /// Looks for the keywords ROCK and ROCKTAB.
        RockCompressibility(const EclipseGridParser& deck);

        /// Construct from input deck.
        /// Looks for the keywords ROCK and ROCKTAB.
        RockCompressibility(Opm::DeckConstPtr newParserDeck);

        /// Construct from parameters.
        /// Accepts the following parameters (with defaults).
        ///    rock_compressibility_pref (100.0)   [given in bar]
        ///    rock_compressibility      (0.0)     [given in bar^{-1}]
        RockCompressibility(const parameter::ParameterGroup& param);

        /// Returns true if there are compressibility effects.
        bool isActive() const;

        /// Porosity multiplier.
        double poroMult(double pressure) const;

        /// Derivative of porosity multiplier with respect to pressure.
        double poroMultDeriv(double pressure) const;

        /// Transmissibility multiplier.
        double transMult(double pressure) const;

        /// Derivative of transmissibility multiplier with respect to pressure.
        double transMultDeriv(double pressure) const;

        /// Rock compressibility = (d poro / d p)*(1 / poro).
        double rockComp(double pressure) const;

    private:
        std::vector<double> p_;
        std::vector<double> poromult_;
        std::vector<double> transmult_;
        double pref_;
        double rock_comp_;
    };

} // namespace Opm


#endif // OPM_ROCKCOMPRESSIBILITY_HEADER_INCLUDED
