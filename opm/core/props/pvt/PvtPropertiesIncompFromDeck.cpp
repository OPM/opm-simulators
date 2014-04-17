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
#include <opm/core/props/pvt/PvtPropertiesIncompFromDeck.hpp>
#include <opm/core/props/phaseUsageFromDeck.hpp>
#include <opm/core/io/eclipse/EclipseGridParser.hpp>
#include <opm/core/utility/Units.hpp>
#include <opm/core/utility/ErrorMacros.hpp>
#include <opm/core/props/BlackoilPhases.hpp>


namespace Opm
{

    PvtPropertiesIncompFromDeck::PvtPropertiesIncompFromDeck()
    {
    }


    void PvtPropertiesIncompFromDeck::init(const EclipseGridParser& deck)
    {
        // If we need multiple regions, this class and the SinglePvt* classes must change.
        int region_number = 0;

        PhaseUsage phase_usage = phaseUsageFromDeck(deck);
        if (phase_usage.phase_used[PhaseUsage::Vapour] ||
            !phase_usage.phase_used[PhaseUsage::Aqua] ||
            !phase_usage.phase_used[PhaseUsage::Liquid]) {
            OPM_THROW(std::runtime_error, "PvtPropertiesIncompFromDeck::init() -- must have gas and oil phases (only) in deck input.\n");
        }

        // Surface densities. Accounting for different orders in eclipse and our code.
        if (deck.hasField("DENSITY")) {
            const std::vector<double>& d = deck.getDENSITY().densities_[region_number];
            enum { ECL_oil = 0, ECL_water = 1, ECL_gas = 2 };
            surface_density_[phase_usage.phase_pos[PhaseUsage::Aqua]]   = d[ECL_water];
            surface_density_[phase_usage.phase_pos[PhaseUsage::Liquid]] = d[ECL_oil];
        } else {
            OPM_THROW(std::runtime_error, "Input is missing DENSITY\n");
        }

        // Make reservoir densities the same as surface densities initially.
        // We will modify them with formation volume factors if found.
        reservoir_density_ = surface_density_;

        // Water viscosity.
        if (deck.hasField("PVTW")) {
            const std::vector<double>& pvtw = deck.getPVTW().pvtw_[region_number];
            if (pvtw[2] != 0.0 || pvtw[4] != 0.0) {
                OPM_MESSAGE("Compressibility effects in PVTW are ignored.");
            }
            reservoir_density_[phase_usage.phase_pos[PhaseUsage::Aqua]] /= pvtw[1];
            viscosity_[phase_usage.phase_pos[PhaseUsage::Aqua]] = pvtw[3];
        } else {
            // Eclipse 100 default.
            // viscosity_[phase_usage.phase_pos[PhaseUsage::Aqua]] = 0.5*Opm::prefix::centi*Opm::unit::Poise;
            OPM_THROW(std::runtime_error, "Input is missing PVTW\n");
        }

        // Oil viscosity.
        if (deck.hasField("PVCDO")) {
            const std::vector<double>& pvcdo = deck.getPVCDO().pvcdo_[region_number];
            if (pvcdo[2] != 0.0 || pvcdo[4] != 0.0) {
                OPM_MESSAGE("Compressibility effects in PVCDO are ignored.");
            }
            reservoir_density_[phase_usage.phase_pos[PhaseUsage::Liquid]] /= pvcdo[1];
            viscosity_[phase_usage.phase_pos[PhaseUsage::Liquid]] = pvcdo[3];
        } else {
            OPM_THROW(std::runtime_error, "Input is missing PVCDO\n");
        }
    }

    void PvtPropertiesIncompFromDeck::init(Opm::DeckConstPtr newParserDeck )
    {
        // If we need multiple regions, this class and the SinglePvt* classes must change.
        int region_number = 0;

        PhaseUsage phase_usage = phaseUsageFromDeck(newParserDeck);
        if (phase_usage.phase_used[PhaseUsage::Vapour] ||
            !phase_usage.phase_used[PhaseUsage::Aqua] ||
            !phase_usage.phase_used[PhaseUsage::Liquid]) {
            OPM_THROW(std::runtime_error, "PvtPropertiesIncompFromDeck::init() -- must have gas and oil phases (only) in deck input.\n");
        }

        // Surface densities. Accounting for different orders in eclipse and our code.
        if (newParserDeck->hasKeyword("DENSITY")) {
            Opm::DeckRecordConstPtr densityRecord = newParserDeck->getKeyword("DENSITY")->getRecord(region_number);
            surface_density_[phase_usage.phase_pos[PhaseUsage::Aqua]]   = densityRecord->getItem("OIL")->getSIDouble(0);
            surface_density_[phase_usage.phase_pos[PhaseUsage::Liquid]] = densityRecord->getItem("WATER")->getSIDouble(0);
        } else {
            OPM_THROW(std::runtime_error, "Input is missing DENSITY\n");
        }

        // Make reservoir densities the same as surface densities initially.
        // We will modify them with formation volume factors if found.
        reservoir_density_ = surface_density_;

        // Water viscosity.
        if (newParserDeck->hasKeyword("PVTW")) {
            Opm::DeckRecordConstPtr pvtwRecord = newParserDeck->getKeyword("PVTW")->getRecord(region_number);
            if (pvtwRecord->getItem("WATER_COMPRESSIBILITY")->getSIDouble(0) != 0.0 ||
                pvtwRecord->getItem("WATER_VISCOSIBILITY")->getSIDouble(0) != 0.0) {
                OPM_MESSAGE("Compressibility effects in PVTW are ignored.");
            }
            reservoir_density_[phase_usage.phase_pos[PhaseUsage::Aqua]] /= pvtwRecord->getItem("WATER_VOL_FACTOR")->getSIDouble(0);
            viscosity_[phase_usage.phase_pos[PhaseUsage::Aqua]] = pvtwRecord->getItem("WATER_VISCOSITY")->getSIDouble(0);
        } else {
            // Eclipse 100 default.
            // viscosity_[phase_usage.phase_pos[PhaseUsage::Aqua]] = 0.5*Opm::prefix::centi*Opm::unit::Poise;
            OPM_THROW(std::runtime_error, "Input is missing PVTW\n");
        }

        // Oil viscosity.
        if (newParserDeck->hasKeyword("PVCDO")) {
            Opm::DeckRecordConstPtr pvcdoRecord = newParserDeck->getKeyword("PVCDO")->getRecord(region_number);

            if (pvcdoRecord->getItem("OIL_COMPRESSIBILITY")->getSIDouble(0) != 0.0 ||
                pvcdoRecord->getItem("OIL_VISCOSIBILITY")->getSIDouble(0) != 0.0) {
                OPM_MESSAGE("Compressibility effects in PVCDO are ignored.");
            }
            reservoir_density_[phase_usage.phase_pos[PhaseUsage::Liquid]] /= pvcdoRecord->getItem("OIL_VOL_FACTOR")->getSIDouble(0);
            viscosity_[phase_usage.phase_pos[PhaseUsage::Liquid]] = pvcdoRecord->getItem("OIL_VISCOSITY")->getSIDouble(0);
        } else {
            OPM_THROW(std::runtime_error, "Input is missing PVCDO\n");
        }
    }

    const double* PvtPropertiesIncompFromDeck::surfaceDensities() const
    {
        return surface_density_.data();
    }


    const double* PvtPropertiesIncompFromDeck::reservoirDensities() const
    {
        return reservoir_density_.data();
    }


    const double* PvtPropertiesIncompFromDeck::viscosity() const
    {
        return viscosity_.data();
    }


    int PvtPropertiesIncompFromDeck::numPhases() const
    {
        return 2;
    }


} // namespace Opm
