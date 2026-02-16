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
template<typename Scalar, typename IndexTraits> class SingleWellState;
template<typename Scalar, typename IndexTraits> class WellState;
template<typename Scalar, typename IndexTraits> class GroupStateHelper;

template<class FluidSystem>
class WellInterfaceFluidSystem : public WellInterfaceGeneric<typename FluidSystem::Scalar, typename FluidSystem::IndexTraitsType>
{
protected:
    using RateConverterType = RateConverter::
    SurfaceToReservoirVoidage<FluidSystem, std::vector<int>>;
    // to indicate a invalid completion
    static constexpr int INVALIDCOMPLETION = std::numeric_limits<int>::max();

public:
    using Scalar = typename FluidSystem::Scalar;
    using IndexTraits = typename FluidSystem::IndexTraitsType;
    using ModelParameters = typename WellInterfaceGeneric<Scalar, IndexTraits>::ModelParameters;
    using GroupStateHelperType = GroupStateHelper<Scalar, IndexTraits>;
    using WellStateType = WellState<Scalar, IndexTraits>;

    static constexpr int Water = IndexTraits::waterPhaseIdx;
    static constexpr int Oil = IndexTraits::oilPhaseIdx;
    static constexpr int Gas = IndexTraits::gasPhaseIdx;

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
                             const int num_conservation_quantities,
                             const int num_phases,
                             const int index_of_well,
                             const std::vector<PerforationData<Scalar>>& perf_data);

    // updating the voidage rates in well_state when requested
    void calculateReservoirRates(const bool use_well_bhp_temperature, SingleWellState<Scalar, IndexTraits>& ws) const;

    bool checkIndividualConstraints(SingleWellState<Scalar, IndexTraits>& ws,
                                    const SummaryState& summaryState,
                                    DeferredLogger& deferred_logger,
                                    const std::optional<Well::InjectionControls>& inj_controls = std::nullopt,
                                    const std::optional<Well::ProductionControls>& prod_controls = std::nullopt) const;

    bool checkGroupConstraints(const GroupStateHelperType& groupStateHelper,
                               const Schedule& schedule,
                               const SummaryState& summaryState,
                               const bool check_guide_rate,
                               WellStateType& well_state) const;

    bool checkConstraints(const GroupStateHelperType& groupStateHelper,
                          const Schedule& schedule,
                          const SummaryState& summaryState,
                          WellStateType& well_state) const;

    std::optional<Scalar>
    getGroupInjectionTargetRate(const Group& group,
                                const GroupStateHelperType& groupStateHelper,
                                const InjectorType& injectorType,
                                Scalar efficiencyFactor) const;

    Scalar
    getGroupProductionTargetRate(const Group& group,
                                 const GroupStateHelperType& groupStateHelper,
                                 Scalar efficiencyFactor) const;

    bool zeroGroupRateTarget(const GroupStateHelperType& groupStateHelper) const;

    // For the conversion between the surface volume rate and reservoir voidage rate
    const RateConverterType& rateConverter_;
};

}

#endif // OPM_WELLINTERFACE_FLUID_SYSTEM_HEADER_INCLUDED
