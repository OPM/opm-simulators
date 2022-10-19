/*
  Copyright 2017 SINTEF Digital, Mathematics and Cybernetics.
  Copyright 2017 Statoil ASA.
  Copyright 2017 IRIS
  Copyright 2019 Norce

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


#ifndef OPM_WELLINTERFACE_FLUID_SYSTEM_HEADER_INCLUDED
#define OPM_WELLINTERFACE_FLUID_SYSTEM_HEADER_INCLUDED

#include <opm/simulators/wells/WellInterfaceGeneric.hpp>
#include <opm/core/props/BlackoilPhases.hpp>

#include <limits>

namespace Opm
{
namespace RateConverter
{
  template <class FluidSystem, class Region> class SurfaceToReservoirVoidage;
}

class Group;
class GroupState;
class Schedule;
class WellState;
class SingleWellState;

template<class FluidSystem>
class WellInterfaceFluidSystem : public WellInterfaceGeneric {
protected:
    using RateConverterType = RateConverter::
    SurfaceToReservoirVoidage<FluidSystem, std::vector<int>>;
    // to indicate a invalid completion
    static constexpr int INVALIDCOMPLETION = std::numeric_limits<int>::max();

public:
    void updateWellTestState(const SingleWellState& ws,
                             const double& simulationTime,
                             const bool& writeMessageToOPMLog,
                             WellTestState& wellTestState,
                             DeferredLogger& deferred_logger) const;

    int flowPhaseToEbosPhaseIdx(const int phaseIdx) const;

    static constexpr int Water = BlackoilPhases::Aqua;
    static constexpr int Oil = BlackoilPhases::Liquid;
    static constexpr int Gas = BlackoilPhases::Vapour;

    const RateConverterType& rateConverter() const
    {
        return rateConverter_;
    }

protected:
    WellInterfaceFluidSystem(const Well& well,
                             const ParallelWellInfo& parallel_well_info,
                             const int time_step,
                             const RateConverterType& rate_converter,
                             const int pvtRegionIdx,
                             const int num_components,
                             const int num_phases,
                             const int index_of_well,
                             const std::vector<PerforationData>& perf_data);

    // updating the voidage rates in well_state when requested
    void calculateReservoirRates(SingleWellState& ws) const;

    bool checkIndividualConstraints(SingleWellState& ws,
                                    const SummaryState& summaryState,
                                    DeferredLogger& deferred_logger) const;

    Well::InjectorCMode activeInjectionConstraint(const SingleWellState& ws,
                                                  const SummaryState& summaryState,
                                                  DeferredLogger& deferred_logger) const;

    Well::ProducerCMode activeProductionConstraint(const SingleWellState& ws,
                                                   const SummaryState& summaryState,
                                                   DeferredLogger& deferred_logger) const;

    std::pair<bool, double> checkGroupConstraintsInj(const Group& group,
                                                     const WellState& well_state,
                                                     const GroupState& group_state,
                                                     const double efficiencyFactor,
                                                     const Schedule& schedule,
                                                     const SummaryState& summaryState,
                                                     DeferredLogger& deferred_logger) const;

    std::pair<bool, double> checkGroupConstraintsProd(const Group& group,
                                                      const WellState& well_state,
                                                      const GroupState& group_state,
                                                      const double efficiencyFactor,
                                                      const Schedule& schedule,
                                                      const SummaryState& summaryState,
                                                      DeferredLogger& deferred_logger) const;

    bool checkGroupConstraints(WellState& well_state,
                               const GroupState& group_state,
                               const Schedule& schedule,
                               const SummaryState& summaryState,
                               DeferredLogger& deferred_logger) const;

    bool checkConstraints(WellState& well_state,
                          const GroupState& group_state,
                          const Schedule& schedule,
                          const SummaryState& summaryState,
                          DeferredLogger& deferred_logger) const;

    bool checkRateEconLimits(const WellEconProductionLimits& econ_production_limits,
                             const double* rates_or_potentials,
                             Opm::DeferredLogger& deferred_logger) const;

    struct RatioLimitCheckReport{
        bool ratio_limit_violated = false;
        int worst_offending_completion = INVALIDCOMPLETION;
        double violation_extent = 0.0;
    };

    void checkMaxWaterCutLimit(const WellEconProductionLimits& econ_production_limits,
                               const SingleWellState& ws,
                               RatioLimitCheckReport& report) const;

    void checkMaxGORLimit(const WellEconProductionLimits& econ_production_limits,
                          const SingleWellState& ws,
                          RatioLimitCheckReport& report) const;

    void checkMaxWGRLimit(const WellEconProductionLimits& econ_production_limits,
                          const SingleWellState& ws,
                          RatioLimitCheckReport& report) const;

    void checkRatioEconLimits(const WellEconProductionLimits& econ_production_limits,
                              const SingleWellState& ws,
                              RatioLimitCheckReport& report,
                              DeferredLogger& deferred_logger) const;

    void updateWellTestStateEconomic(const SingleWellState& ws,
                                     const double simulation_time,
                                     const bool write_message_to_opmlog,
                                     WellTestState& well_test_state,
                                     DeferredLogger& deferred_logger) const;

    std::optional<double>
    getGroupInjectionTargetRate(const Group& group,
                                const WellState& well_state,
                                const GroupState& group_state,
                                const Schedule& schedule,
                                const SummaryState& summaryState,
                                const InjectorType& injectorType,
                                double efficiencyFactor,
                                DeferredLogger& deferred_logger) const;

    double
    getGroupProductionTargetRate(const Group& group,
                                 const WellState& well_state,
                                 const GroupState& group_state,
                                 const Schedule& schedule,
                                 const SummaryState& summaryState,
                                 double efficiencyFactor) const;

    // For the conversion between the surface volume rate and reservoir voidage rate
    const RateConverterType& rateConverter_;

private:
    template <typename RatioFunc>
    void checkMaxRatioLimitCompletions(const SingleWellState& ws,
                                       const double max_ratio_limit,
                                       const RatioFunc& ratioFunc,
                                       RatioLimitCheckReport& report) const;
};

}

#endif // OPM_WELLINTERFACE_FLUID_SYSTEM_HEADER_INCLUDED
