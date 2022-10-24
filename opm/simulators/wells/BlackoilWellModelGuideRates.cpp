/*
  Copyright 2016 SINTEF ICT, Applied Mathematics.
  Copyright 2016 - 2017 Statoil ASA.
  Copyright 2017 Dr. Blatt - HPC-Simulation-Software & Services
  Copyright 2016 - 2018 IRIS AS

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

#include <config.h>
#include <opm/simulators/wells/BlackoilWellModelGuideRates.hpp>

#include <opm/input/eclipse/Schedule/Group/Group.hpp>

#include <opm/output/data/GuideRateValue.hpp>

#include <opm/simulators/wells/BlackoilWellModelGeneric.hpp>
#include <opm/simulators/wells/WellGroupHelpers.hpp>

namespace Opm {

void
BlackoilWellModelGuideRates::
getGuideRateValues(const GuideRate::RateVector& qs,
                   const bool                   is_inj,
                   const std::string&           wgname,
                   data::GuideRateValue&        grval) const
{
    auto getGR = [this, &wgname, &qs](const GuideRateModel::Target t)
    {
        return wellModel_.guideRate().getSI(wgname, t, qs);
    };

    // Note: GuideRate does currently (2020-07-20) not support Target::RES.
    grval.set(data::GuideRateValue::Item::Gas,
              getGR(GuideRateModel::Target::GAS));

    grval.set(data::GuideRateValue::Item::Water,
              getGR(GuideRateModel::Target::WAT));

    if (!is_inj) {
        // Producer.  Extract "all" guiderate values.
        grval.set(data::GuideRateValue::Item::Oil,
                  getGR(GuideRateModel::Target::OIL));
    }
}

data::GuideRateValue
BlackoilWellModelGuideRates::
getGuideRateValues(const Well& well) const
{
    auto grval = data::GuideRateValue{};

    const auto& wname = well.name();
    if (!wellModel_.wellState().has(wname)) {
        // No flow rates for 'wname' -- might be before well comes
        // online (e.g., for the initial condition before simulation
        // starts).
        return grval;
    }

    if (!wellModel_.guideRate().has(wname)) {
        // No guiderates exist for 'wname'.
        return grval;
    }

    const auto qs = WellGroupHelpers::
        getWellRateVector(wellModel_.wellState(), wellModel_.phaseUsage(), wname);

    this->getGuideRateValues(qs, well.isInjector(), wname, grval);

    return grval;
}

data::GuideRateValue
BlackoilWellModelGuideRates::
getGuideRateValues(const Group& group) const
{
    auto grval = data::GuideRateValue{};
    const auto& gname = group.name();

    if (!wellModel_.groupState().has_production_rates(gname)) {
        // No flow rates for production group 'gname' -- might be before
        // group comes online (e.g., for the initial condition before
        // simulation starts).
        return grval;
    }

    if (!wellModel_.guideRate().has(gname)) {
        // No guiderates exist for 'gname'.
        return grval;
    }

    const auto qs = WellGroupHelpers::
        getProductionGroupRateVector(wellModel_.groupState(), wellModel_.phaseUsage(), gname);

    const auto is_inj = false; // This procedure only applies to G*PGR.
    this->getGuideRateValues(qs, is_inj, gname, grval);

    return grval;
}

data::GuideRateValue
BlackoilWellModelGuideRates::
getGuideRateInjectionGroupValues(const Group& group) const
{
    auto grval = data::GuideRateValue{};

    const auto& gname = group.name();
    if (wellModel_.guideRate().has(gname, Phase::GAS)) {
        grval.set(data::GuideRateValue::Item::Gas,
                  wellModel_.guideRate().getSI(gname, Phase::GAS));
    }

    if (wellModel_.guideRate().has(gname, Phase::WATER)) {
        grval.set(data::GuideRateValue::Item::Water,
                  wellModel_.guideRate().getSI(gname, Phase::WATER));
    }

    return grval;
}

}
