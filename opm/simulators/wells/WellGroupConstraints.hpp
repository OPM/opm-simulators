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


#ifndef OPM_WELL_GROUP_CONSTRAINTS_HEADER_INCLUDED
#define OPM_WELL_GROUP_CONSTRAINTS_HEADER_INCLUDED

#include <functional>
#include <utility>
#include <vector>
#include <string>
#include <optional>

namespace Opm
{

class DeferredLogger;
class Group;
template<class Scalar> class GroupState;
enum class InjectorType;
using RegionId = int;
class Schedule;
class SummaryState;
template<class Scalar> class WellInterfaceGeneric;
template<class Scalar> class WellState;

//! \brief Class for computing well group constraints.
template<class Scalar>
class WellGroupConstraints {
public:
    //! \brief Constructor sets reference to well.
    explicit WellGroupConstraints(const WellInterfaceGeneric<Scalar>& well) : well_(well) {}

    using RateConvFunc = std::function<void(const RegionId,
                                            const int,
                                            const std::optional<std::string>&,
                                            std::vector<Scalar>&)>;

    bool checkGroupConstraints(WellState<Scalar>& well_state,
                               const GroupState<Scalar>& group_state,
                               const Schedule& schedule,
                               const SummaryState& summaryState,
                               const RateConvFunc& rateConverter,
                               const bool check_guide_rate,
                               DeferredLogger& deferred_logger) const;

private:
    std::pair<bool, Scalar>
    checkGroupConstraintsInj(const Group& group,
                             const WellState<Scalar>& well_state,
                             const GroupState<Scalar>& group_state,
                             const Scalar efficiencyFactor,
                             const Schedule& schedule,
                             const SummaryState& summaryState,
                             const RateConvFunc& rateConverter,
                             const bool check_guide_rate,
                             DeferredLogger& deferred_logger) const;

    std::pair<bool, Scalar>
    checkGroupConstraintsProd(const Group& group,
                              const WellState<Scalar>& well_state,
                              const GroupState<Scalar>& group_state,
                              const Scalar efficiencyFactor,
                              const Schedule& schedule,
                              const SummaryState& summaryState,
                              const RateConvFunc& rateConverter,
                              const bool check_guide_rate,
                              DeferredLogger& deferred_logger) const;

    const WellInterfaceGeneric<Scalar>& well_; //!< Reference to well interface
};

}

#endif // OPM_WELL_GROUP_CONSTRAINTS_HEADER_INCLUDED
