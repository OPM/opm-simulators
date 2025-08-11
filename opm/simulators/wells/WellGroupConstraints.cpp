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

#include <opm/simulators/wells/WellGroupConstraints.hpp>

#include <opm/input/eclipse/Schedule/Schedule.hpp>

#include <opm/material/fluidsystems/BlackOilFluidSystem.hpp>

#include <opm/models/blackoil/blackoilvariableandequationindices.hh>
#include <opm/models/blackoil/blackoilonephaseindices.hh>
#include <opm/models/blackoil/blackoiltwophaseindices.hh>

#include <opm/simulators/wells/WellGroupHelpers.hpp>
#include <opm/simulators/wells/WellInterfaceGeneric.hpp>
#include <opm/simulators/wells/WellState.hpp>

#include <opm/simulators/utils/BlackoilPhases.hpp>

namespace Opm
{

template<typename FluidSystem, typename Indices>
std::pair<bool, typename FluidSystem::Scalar>
WellGroupConstraints<FluidSystem, Indices>::
checkGroupConstraintsInj(const Group& group,
                         const WellState<FluidSystem, Indices>& well_state,
                         const GroupState<Scalar>& group_state,
                         const Scalar efficiencyFactor,
                         const Schedule& schedule,
                         const SummaryState& summaryState,
                         const RateConvFunc& rateConverter,
                         const bool check_guide_rate,
                         DeferredLogger& deferred_logger) const
{
    // Translate injector type from control to Phase.
    const auto& well_controls = well_.wellEcl().injectionControls(summaryState);
    auto injectorType = well_controls.injector_type;
    Phase injectionPhase;
    switch (injectorType) {
    case InjectorType::WATER:
    {
        injectionPhase = Phase::WATER;
        break;
    }
    case InjectorType::OIL:
    {
        injectionPhase = Phase::OIL;
        break;
    }
    case InjectorType::GAS:
    {
        injectionPhase = Phase::GAS;
        break;
    }
    default:
        throw("Expected WATER, OIL or GAS as type for injector " + well_.name());
    }

    // Make conversion factors for RESV <-> surface rates.
    std::vector<Scalar> resv_coeff(Indices::numPhases, 1.0);
    rateConverter(0, well_.pvtRegionIdx(), group.name(), resv_coeff); // FIPNUM region 0 here, should use FIPNUM from WELSPECS.

    const auto& ws = well_state.well(well_.indexOfWell());
    // Call check for the well's injection phase.
    return WellGroupHelpers<FluidSystem, Indices>::checkGroupConstraintsInj(well_.name(),
                                                              well_.wellEcl().groupName(),
                                                              group,
                                                              well_state,
                                                              group_state,
                                                              well_.currentStep(),
                                                              well_.guideRate(),
                                                              ws.surface_rates.data(),
                                                              injectionPhase,
                                                              efficiencyFactor,
                                                              schedule,
                                                              summaryState,
                                                              resv_coeff,
                                                              check_guide_rate,
                                                              deferred_logger);
}

template<typename FluidSystem, typename Indices>
std::pair<bool, typename FluidSystem::Scalar>
WellGroupConstraints<FluidSystem, Indices>::
checkGroupConstraintsProd(const Group& group,
                          const WellState<FluidSystem, Indices>& well_state,
                          const GroupState<Scalar>& group_state,
                          const Scalar efficiencyFactor,
                          const Schedule& schedule,
                          const SummaryState& summaryState,
                          const RateConvFunc& rateConverter,
                          const bool check_guide_rate,
                          DeferredLogger& deferred_logger) const
{
    // Make conversion factors for RESV <-> surface rates.
    std::vector<Scalar> resv_coeff(Indices::numPhases, 1.0);
    rateConverter(0, well_.pvtRegionIdx(), group.name(), resv_coeff); // FIPNUM region 0 here, should use FIPNUM from WELSPECS.

    const auto& ws = well_state.well(well_.indexOfWell());
    return WellGroupHelpers<FluidSystem, Indices>::checkGroupConstraintsProd(well_.name(),
                                                               well_.wellEcl().groupName(),
                                                               group,
                                                               well_state,
                                                               group_state,
                                                               well_.currentStep(),
                                                               well_.guideRate(),
                                                               ws.surface_rates.data(),
                                                               efficiencyFactor,
                                                               schedule,
                                                               summaryState,
                                                               resv_coeff,
                                                               check_guide_rate,
                                                               deferred_logger);
}

template<typename FluidSystem, typename Indices>
bool WellGroupConstraints<FluidSystem, Indices>::
checkGroupConstraints(WellState<FluidSystem, Indices>& well_state,
                      const GroupState<Scalar>& group_state,
                      const Schedule& schedule,
                      const SummaryState& summaryState,
                      const RateConvFunc& rateConverter,
                      const bool check_guide_rate,
                      DeferredLogger& deferred_logger) const
{
    const auto& well = well_.wellEcl();
    const int well_index = well_.indexOfWell();
    auto& ws = well_state.well(well_index);

    if (well.isInjector()) {
        const auto currentControl = ws.injection_cmode;

        if (currentControl != Well::InjectorCMode::GRUP) {
            // This checks only the first encountered group limit,
            // in theory there could be several, and then we should
            // test all but the one currently applied. At that point,
            // this if-statement should be removed and we should always
            // check, skipping over only the single group parent whose
            // control is the active one for the well (if any).
            const auto& group = schedule.getGroup(well.groupName(), well_.currentStep());
            const Scalar efficiencyFactor = well.getEfficiencyFactor() *
                                            well_state[well.name()].efficiency_scaling_factor;
            const std::pair<bool, Scalar> group_constraint =
                this->checkGroupConstraintsInj(group, well_state,
                                               group_state, efficiencyFactor,
                                               schedule, summaryState,
                                               rateConverter,
                                               check_guide_rate,
                                               deferred_logger);
            // If a group constraint was broken, we set the current well control to
            // be GRUP.
            if (group_constraint.first) {
                ws.injection_cmode = Well::InjectorCMode::GRUP;
                const int np = well_state.numPhases();
                for (int p = 0; p<np; ++p) {
                    ws.surface_rates[p] *= group_constraint.second;
                }
            }
            return group_constraint.first;
        }
    }

    if (well.isProducer( )) {
        const auto currentControl = ws.production_cmode;

        if (currentControl != Well::ProducerCMode::GRUP) {
            // This checks only the first encountered group limit,
            // in theory there could be several, and then we should
            // test all but the one currently applied. At that point,
            // this if-statement should be removed and we should always
            // check, skipping over only the single group parent whose
            // control is the active one for the well (if any).
            const auto& group = schedule.getGroup(well.groupName(), well_.currentStep());
            const Scalar efficiencyFactor = well.getEfficiencyFactor() *
                                            well_state[well.name()].efficiency_scaling_factor;
            const std::pair<bool, Scalar> group_constraint =
                this->checkGroupConstraintsProd(group, well_state,
                                                group_state, efficiencyFactor,
                                                schedule, summaryState,
                                                rateConverter, check_guide_rate, deferred_logger);
            // If a group constraint was broken, we set the current well control to
            // be GRUP.
            if (group_constraint.first) {
                ws.production_cmode = Well::ProducerCMode::GRUP;
                const int np = well_state.numPhases();
                for (int p = 0; p<np; ++p) {
                    ws.surface_rates[p] *= group_constraint.second;
                }
            }
            return group_constraint.first;
        }
    }

    return false;
}


#include <opm/simulators/utils/InstantiationIndicesMacros.hpp>

INSTANTIATE_TYPE_INDICES(WellGroupConstraints, double)

#if FLOW_INSTANTIATE_FLOAT
INSTANTIATE_TYPE_INDICES(WellGroupConstraints, float)
#endif

} // namespace Opm
