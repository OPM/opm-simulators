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

#include <string>
#include <functional>
#include <optional>
#include <vector>

namespace Opm
{

class DeferredLogger;
class Group;
class GroupState;
enum class InjectorType;
using RegionId = int;
class Schedule;
class SummaryState;
class WellInterfaceGeneric;
class WellState;

//! \brief Class for computing well group controls.
class WellGroupControls {
public:
    //! \brief Constructor sets reference to well.
    WellGroupControls(const WellInterfaceGeneric& well) : well_(well) {}

    using RateConvFunc = std::function<void(const RegionId, const int, const std::optional<std::string>&, std::vector<double>&)>;

    template<class EvalWell>
    void getGroupInjectionControl(const Group& group,
                                  const WellState& well_state,
                                  const GroupState& group_state,
                                  const Schedule& schedule,
                                  const SummaryState& summaryState,
                                  const InjectorType& injectorType,
                                  const EvalWell& bhp,
                                  const EvalWell& injection_rate,
                                  const RateConvFunc& rateConverter,
                                  double efficiencyFactor,
                                  EvalWell& control_eq,
                                  DeferredLogger& deferred_logger) const;

    std::optional<double>
    getGroupInjectionTargetRate(const Group& group,
                                const WellState& well_state,
                                const GroupState& group_state,
                                const Schedule& schedule,
                                const SummaryState& summaryState,
                                const InjectorType& injectorType,
                                const RateConvFunc& rateConverter,
                                double efficiencyFactor,
                                DeferredLogger& deferred_logger) const;

    template<class EvalWell>
    void getGroupProductionControl(const Group& group,
                                   const WellState& well_state,
                                   const GroupState& group_state,
                                   const Schedule& schedule,
                                   const SummaryState& summaryState,
                                   const EvalWell& bhp,
                                   const std::vector<EvalWell>& rates,
                                   const RateConvFunc& rateConverter,
                                   double efficiencyFactor,
                                   EvalWell& control_eq,
                                   DeferredLogger& deferred_logger) const;

    double getGroupProductionTargetRate(const Group& group,
                                        const WellState& well_state,
                                        const GroupState& group_state,
                                        const Schedule& schedule,
                                        const SummaryState& summaryState,
                                        const RateConvFunc& rateConverter,
                                        double efficiencyFactor,
                                        DeferredLogger& deferred_logger) const;

private:
    const WellInterfaceGeneric& well_; //!< Reference to well interface
};

}

#endif // OPM_WELL_GROUP_CONTROLS_HEADER_INCLUDED
