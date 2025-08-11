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

#include <opm/simulators/utils/BlackoilPhases.hpp>

#include <limits>
#include <optional>
#include <vector>

namespace Opm
{
namespace RateConverter
{
  template <class FluidSystem, class Region> class SurfaceToReservoirVoidage;
}

class Group;
template<class Scalar> class GroupState;
class Schedule;
struct RatioLimitCheckReport;
template<class Scalar> class SingleWellState;
template<typename FluidSystem, typename Indices> class WellState;

template<typename FluidSystem, typename Indices>
class WellInterfaceFluidSystem : public WellInterfaceGeneric<FluidSystem, Indices>
{
protected:
    using RateConverterType = RateConverter::
    SurfaceToReservoirVoidage<FluidSystem, std::vector<int>>;
    // to indicate a invalid completion
    static constexpr int INVALIDCOMPLETION = std::numeric_limits<int>::max();

public:
    using Scalar = typename FluidSystem::Scalar;
    using ModelParameters = typename WellInterfaceGeneric<FluidSystem, Indices>::ModelParameters;

    static constexpr int Water = BlackoilPhases::Aqua;
    static constexpr int Oil = BlackoilPhases::Liquid;
    static constexpr int Gas = BlackoilPhases::Vapour;

    const RateConverterType& rateConverter() const
    {
        return rateConverter_;
    }

protected:
    WellInterfaceFluidSystem(const Well& well,
                             const ParallelWellInfo<Scalar>& parallel_well_info,
                             const int time_step,
                             const ModelParameters& param,
                             const RateConverterType& rate_converter,
                             const int pvtRegionIdx,
                             const int num_components,
                             const int num_phases,
                             const int index_of_well,
                             const std::vector<PerforationData<Scalar>>& perf_data);

    // updating the voidage rates in well_state when requested
    void calculateReservoirRates(const bool co2store, SingleWellState<Scalar>& ws) const;

    bool checkIndividualConstraints(SingleWellState<Scalar>& ws,
                                    const SummaryState& summaryState,
                                    DeferredLogger& deferred_logger,
                                    const std::optional<Well::InjectionControls>& inj_controls = std::nullopt,
                                    const std::optional<Well::ProductionControls>& prod_controls = std::nullopt) const;

    bool checkGroupConstraints(WellState<FluidSystem, Indices>& well_state,
                               const GroupState<Scalar>& group_state,
                               const Schedule& schedule,
                               const SummaryState& summaryState,
                               const bool check_guide_rate,
                               DeferredLogger& deferred_logger) const;

    bool checkConstraints(WellState<FluidSystem, Indices>& well_state,
                          const GroupState<Scalar>& group_state,
                          const Schedule& schedule,
                          const SummaryState& summaryState,
                          DeferredLogger& deferred_logger) const;

    std::optional<Scalar>
    getGroupInjectionTargetRate(const Group& group,
                                const WellState<FluidSystem, Indices>& well_state,
                                const GroupState<Scalar>& group_state,
                                const Schedule& schedule,
                                const SummaryState& summaryState,
                                const InjectorType& injectorType,
                                Scalar efficiencyFactor,
                                DeferredLogger& deferred_logger) const;

    Scalar
    getGroupProductionTargetRate(const Group& group,
                                 const WellState<FluidSystem, Indices>& well_state,
                                 const GroupState<Scalar>& group_state,
                                 const Schedule& schedule,
                                 const SummaryState& summaryState,
                                 Scalar efficiencyFactor,
                                 DeferredLogger& deferred_logger) const;

    bool zeroGroupRateTarget(const SummaryState& summary_state,
                             const Schedule& schedule,
                             const WellState<FluidSystem, Indices>& well_state,
                             const GroupState<Scalar>& group_state,
                             DeferredLogger& deferredLogger) const;

    // For the conversion between the surface volume rate and reservoir voidage rate
    const RateConverterType& rateConverter_;
};

}

#endif // OPM_WELLINTERFACE_FLUID_SYSTEM_HEADER_INCLUDED
