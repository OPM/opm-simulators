/*
  Copyright 2017 SINTEF Digital, Mathematics and Cybernetics.
  Copyright 2017 Statoil ASA.
  Copyright 2018 IRIS

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
#include <opm/simulators/wells/WellBhpThpCalculator.hpp>

#include <opm/core/props/BlackoilPhases.hpp>

#include <opm/input/eclipse/Schedule/VFPInjTable.hpp>
#include <opm/input/eclipse/Schedule/Well/Well.hpp>

#include <opm/simulators/utils/DeferredLoggingErrorHelpers.hpp>

#include <opm/simulators/wells/VFPProperties.hpp>
#include <opm/simulators/wells/WellHelpers.hpp>
#include <opm/simulators/wells/WellInterfaceGeneric.hpp>

#include <cassert>

namespace Opm
{

bool
WellBhpThpCalculator::wellHasTHPConstraints(const SummaryState& summaryState) const
{
    const auto& well_ecl = well_.wellEcl();
    if (well_ecl.isInjector()) {
        const auto controls = well_ecl.injectionControls(summaryState);
        if (controls.hasControl(Well::InjectorCMode::THP))
            return true;
    }

    if (well_ecl.isProducer()) {
        const auto controls = well_ecl.productionControls(summaryState);
        if (controls.hasControl(Well::ProducerCMode::THP))
            return true;
    }

    return false;
}

double WellBhpThpCalculator::getTHPConstraint(const SummaryState& summaryState) const
{
    const auto& well_ecl = well_.wellEcl();
    if (well_ecl.isInjector()) {
        const auto& controls = well_ecl.injectionControls(summaryState);
        return controls.thp_limit;
    }

    if (well_ecl.isProducer( )) {
        const auto& controls = well_ecl.productionControls(summaryState);
        return controls.thp_limit;
    }

    return 0.0;
}

double WellBhpThpCalculator::mostStrictBhpFromBhpLimits(const SummaryState& summaryState) const
{
    const auto& well_ecl = well_.wellEcl();
    if (well_ecl.isInjector()) {
        const auto& controls = well_ecl.injectionControls(summaryState);
        return controls.bhp_limit;
    }

    if (well_ecl.isProducer( )) {
        const auto& controls = well_ecl.productionControls(summaryState);
        return controls.bhp_limit;
    }

    return 0.0;
}

double WellBhpThpCalculator::calculateThpFromBhp(const std::vector<double>& rates,
                                                 const double bhp,
                                                 const double rho,
                                                 const double alq,
                                                 DeferredLogger& deferred_logger) const
{
    assert(int(rates.size()) == 3); // the vfp related only supports three phases now.

    static constexpr int Water = BlackoilPhases::Aqua;
    static constexpr int Oil = BlackoilPhases::Liquid;
    static constexpr int Gas = BlackoilPhases::Vapour;

    const double aqua = rates[Water];
    const double liquid = rates[Oil];
    const double vapour = rates[Gas];

    // pick the density in the top layer
    double thp = 0.0;
    if (well_.isInjector()) {
        const int table_id = well_.wellEcl().vfp_table_number();
        const double vfp_ref_depth = well_.vfpProperties()->getInj()->getTable(table_id).getDatumDepth();
        const double dp = wellhelpers::computeHydrostaticCorrection(well_.refDepth(), vfp_ref_depth, rho, well_.gravity());

        thp = well_.vfpProperties()->getInj()->thp(table_id, aqua, liquid, vapour, bhp + dp);
    }
    else if (well_.isProducer()) {
        const int table_id = well_.wellEcl().vfp_table_number();
        const double vfp_ref_depth = well_.vfpProperties()->getProd()->getTable(table_id).getDatumDepth();
        const double dp = wellhelpers::computeHydrostaticCorrection(well_.refDepth(), vfp_ref_depth, rho, well_.gravity());

        thp = well_.vfpProperties()->getProd()->thp(table_id, aqua, liquid, vapour, bhp + dp, alq);
    }
    else {
        OPM_DEFLOG_THROW(std::logic_error, "Expected INJECTOR or PRODUCER well", deferred_logger);
    }

    return thp;
}

} // namespace Opm
