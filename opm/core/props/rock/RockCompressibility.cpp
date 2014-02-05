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

#include "config.h"
#include <opm/core/props/rock/RockCompressibility.hpp>
#include <opm/core/io/eclipse/EclipseGridParser.hpp>
#include <opm/core/utility/parameters/ParameterGroup.hpp>
#include <opm/core/utility/Units.hpp>
#include <opm/core/utility/ErrorMacros.hpp>
#include <opm/core/utility/linearInterpolation.hpp>

#include <opm/parser/eclipse/Utility/RocktabTable.hpp>

#include <iostream>

namespace Opm
{

    RockCompressibility::RockCompressibility(const parameter::ParameterGroup& param)
        : pref_(0.0),
          rock_comp_(0.0)
    {
        pref_ = param.getDefault("rock_compressibility_pref", 100.0)*unit::barsa;
        rock_comp_ = param.getDefault("rock_compressibility", 0.0)/unit::barsa;
    }

    RockCompressibility::RockCompressibility(const EclipseGridParser& deck)
        : pref_(0.0),
          rock_comp_(0.0)
    {
        if (deck.hasField("ROCKTAB")) {
            const table_t& rt = deck.getROCKTAB().rocktab_;
            if (rt.size() != 1) {
                OPM_THROW(std::runtime_error, "Can only handle a single region in ROCKTAB.");
            }
            const int n = rt[0][0].size();
            p_.resize(n);
            poromult_.resize(n);
            transmult_.resize(n);
            for (int i = 0; i < n; ++i) {
                p_[i] = rt[0][0][i];
                poromult_[i] = rt[0][1][i];
                transmult_[i] = rt[0][2][i];
            }
        } else if (deck.hasField("ROCK")) {
            const ROCK& r = deck.getROCK();
            pref_ = r.rock_compressibilities_[0][0];
            rock_comp_ = r.rock_compressibilities_[0][1];
        } else {
            std::cout << "**** warning: no rock compressibility data found in deck (ROCK or ROCKTAB)." << std::endl;
        }
    }

    RockCompressibility::RockCompressibility(Opm::DeckConstPtr newParserDeck)
        : pref_(0.0),
          rock_comp_(0.0)
    {
        if (newParserDeck->hasKeyword("ROCKTAB")) {
            Opm::DeckKeywordConstPtr rtKeyword = newParserDeck->getKeyword("ROCKTAB");
            if (rtKeyword->size() != 1)
                OPM_THROW(std::runtime_error, "Can only handle a single region in ROCKTAB.");

            // the number of colums of the "ROCKTAB" keyword
            // depends on the presence of the "RKTRMDIR"
            // keyword. Messy stuff...
            bool isDirectional = newParserDeck->hasKeyword("RKTRMDIR");
            if (isDirectional)
            {
                // well, okay. we don't support non-isotropic
                // transmisibility multipliers yet
                OPM_THROW(std::runtime_error, "Support for non-isotropic "
                          "transmisibility multipliers is not implemented yet.");
            };

            Opm::RocktabTable rocktabTable(rtKeyword, isDirectional);

            p_ = rocktabTable.getPressureColumn();
            poromult_ = rocktabTable.getPoreVolumeMultiplierColumn();
            transmult_ =  rocktabTable.getTransmisibilityMultiplierColumn();
        } else if (newParserDeck->hasKeyword("ROCK")) {
            Opm::DeckKeywordConstPtr rockKeyword = newParserDeck->getKeyword("ROCK");
            if (rockKeyword->size() != 1)
                OPM_THROW(std::runtime_error, "Can only handle a single region in ROCK.");

            Opm::DeckRecordConstPtr rockRecord = rockKeyword->getRecord(0);
            pref_ = rockRecord->getItem(0)->getSIDouble(0);
            rock_comp_ = rockRecord->getItem(1)->getSIDouble(0);
        } else {
            std::cout << "**** warning: no rock compressibility data found in deck (ROCK or ROCKTAB)." << std::endl;
        }
    }

    bool RockCompressibility::isActive() const
    {
        return !p_.empty() || (rock_comp_ != 0.0);
    }

    double RockCompressibility::poroMult(double pressure) const
    {
        if (p_.empty()) {
            // Approximating with a quadratic curve.
            const double cpnorm = rock_comp_*(pressure - pref_);
            return (1.0 + cpnorm + 0.5*cpnorm*cpnorm);
        } else {
            return Opm::linearInterpolation(p_, poromult_, pressure);
        }
    }

    double RockCompressibility::poroMultDeriv(double pressure) const
    {
        if (p_.empty()) {
            // Approximating poro multiplier with a quadratic curve,
            // we must use its derivative.
            return rock_comp_ + 2 * rock_comp_ * rock_comp_ * (pressure - pref_);
        } else {
            return Opm::linearInterpolationDerivative(p_, poromult_, pressure);
        }
    }

    double RockCompressibility::transMult(double pressure) const
    {
        if (p_.empty()) {
            return 1.0;
        } else {
            return Opm::linearInterpolation(p_, transmult_, pressure);
        }
    }

    double RockCompressibility::transMultDeriv(double pressure) const
    {
        if (p_.empty()) {
            return 0.0;
        } else {
            return Opm::linearInterpolationDerivative(p_, transmult_, pressure);
        }
    }

    double RockCompressibility::rockComp(double pressure) const
    {
        if (p_.empty()) {
            return rock_comp_;
        } else {
            //const double poromult = Opm::linearInterpolation(p_, poromult_, pressure);
            //const double dporomultdp = Opm::linearInterpolationDerivative(p_, poromult_, pressure);
            const double poromult = Opm::linearInterpolation(p_, poromult_, pressure);
            const double dporomultdp = Opm::linearInterpolationDerivative(p_, poromult_, pressure);

            return dporomultdp/poromult;
        }
    }

} // namespace Opm

