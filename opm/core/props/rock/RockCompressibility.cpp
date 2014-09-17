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
#include <opm/core/utility/parameters/ParameterGroup.hpp>
#include <opm/core/utility/Units.hpp>
#include <opm/core/utility/ErrorMacros.hpp>
#include <opm/core/utility/linearInterpolation.hpp>


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

    RockCompressibility::RockCompressibility(Opm::DeckConstPtr deck,
                                             Opm::EclipseStateConstPtr eclipseState)
        : pref_(0.0),
          rock_comp_(0.0)
    {
        const auto& rocktabTables = eclipseState->getRocktabTables();
        if (rocktabTables.size() > 0) {
            if (rocktabTables.size() != 1)
                OPM_THROW(std::runtime_error, "Can only handle a single region in ROCKTAB.");

            p_ = rocktabTables[0].getPressureColumn();
            poromult_ = rocktabTables[0].getPoreVolumeMultiplierColumn();
            transmult_ =  rocktabTables[0].getTransmissibilityMultiplierColumn();
        } else if (deck->hasKeyword("ROCK")) {
            Opm::DeckKeywordConstPtr rockKeyword = deck->getKeyword("ROCK");
            if (rockKeyword->size() != 1) {
                // here it would be better not to use std::cout directly but to add the
                // warning to some "warning list"...
                std::cout << "Can only handle a single region in ROCK ("<<rockKeyword->size()<<" regions specified)."
                          << " Ignoring all except for the first.\n";
            }

            pref_ = rockKeyword->getRecord(0)->getItem("PREF")->getSIDouble(0);
            rock_comp_ = rockKeyword->getRecord(0)->getItem("COMPRESSIBILITY")->getSIDouble(0);
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

