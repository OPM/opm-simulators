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


#ifndef OPM_WELL_GROUP_CONTROLS_HEADER_INCLUDED
#define OPM_WELL_GROUP_CONTROLS_HEADER_INCLUDED

#include <opm/input/eclipse/Schedule/Group/GuideRate.hpp>
#include <string>
#include <functional>
#include <optional>
#include <vector>

namespace Opm
{

class DeferredLogger;
class Group;
template<class Scalar> class GroupState;
enum class InjectorType;
using RegionId = int;
class Schedule;
class SummaryState;
template<typename Scalar, typename IndexTraits> class WellInterfaceGeneric;
template<typename Scalar, typename IndexTraits> class WellState;

//! \brief Class for computing well group controls.
template<typename Scalar, typename IndexTraits>
class WellGroupControls {
public:
    //! \brief Constructor sets reference to well.
    explicit WellGroupControls(const WellInterfaceGeneric<Scalar, IndexTraits>& well) : well_(well) {}

    using RateConvFunc = std::function<void(const RegionId, const int,
                                            const std::optional<std::string>&, std::vector<Scalar>&)>;

    template<class EvalWell>
    void getGroupInjectionControl(const Group& group,
                                  const WellState<Scalar, IndexTraits>& well_state,
                                  const GroupState<Scalar>& group_state,
                                  const Schedule& schedule,
                                  const SummaryState& summaryState,
                                  const InjectorType& injectorType,
                                  const EvalWell& bhp,
                                  const EvalWell& injection_rate,
                                  const RateConvFunc& rateConverter,
                                  Scalar efficiencyFactor,
                                  EvalWell& control_eq,
                                  DeferredLogger& deferred_logger) const;

    std::optional<Scalar>
    getGroupInjectionTargetRate(const Group& group,
                                const WellState<Scalar, IndexTraits>& well_state,
                                const GroupState<Scalar>& group_state,
                                const Schedule& schedule,
                                const SummaryState& summaryState,
                                const InjectorType& injectorType,
                                const RateConvFunc& rateConverter,
                                Scalar efficiencyFactor,
                                DeferredLogger& deferred_logger) const;

    template<class EvalWell>
    void getGroupProductionControl(const Group& group,
                                   const WellState<Scalar, IndexTraits>& well_state,
                                   const GroupState<Scalar>& group_state,
                                   const Schedule& schedule,
                                   const SummaryState& summaryState,
                                   const EvalWell& bhp,
                                   const std::vector<EvalWell>& rates,
                                   const RateConvFunc& rateConverter,
                                   Scalar efficiencyFactor,
                                   EvalWell& control_eq,
                                   DeferredLogger& deferred_logger) const;

    Scalar getGroupProductionTargetRate(const Group& group,
                                        const WellState<Scalar, IndexTraits>& well_state,
                                        const GroupState<Scalar>& group_state,
                                        const Schedule& schedule,
                                        const SummaryState& summaryState,
                                        const RateConvFunc& rateConverter,
                                        Scalar efficiencyFactor,
                                        DeferredLogger& deferred_logger) const;

    static std::pair<Scalar, Group::ProductionCMode> getAutoChokeGroupProductionTargetRate(const std::string& name,
                                                        const Group& parent,
                                                        const WellState<Scalar, IndexTraits>& well_state,
                                                        const GroupState<Scalar>& group_state,
                                                        const Schedule& schedule,
                                                        const SummaryState& summaryState,
                                                        const std::vector<Scalar>& resv_coeff,
                                                        Scalar efficiencyFactor,
                                                        const int reportStepIdx,
                                                        const GuideRate* guideRate,
                                                        DeferredLogger& deferred_logger);

private:
    const WellInterfaceGeneric<Scalar, IndexTraits>& well_; //!< Reference to well interface
};

}

#endif // OPM_WELL_GROUP_CONTROLS_HEADER_INCLUDED
