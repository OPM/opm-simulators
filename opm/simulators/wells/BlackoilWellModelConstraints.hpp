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

#ifndef OPM_BLACKOILWELLMODEL_CONSTRAINTS_HEADER_INCLUDED
#define OPM_BLACKOILWELLMODEL_CONSTRAINTS_HEADER_INCLUDED

#include <opm/input/eclipse/Schedule/Group/Group.hpp>
#include <optional>
#include <utility>

namespace Opm {

template<typename Scalar, typename IndexTraits> class BlackoilWellModelGeneric;
class DeferredLogger;
template<class Scalar> class GroupState;
class SummaryState;
template<typename Scalar, typename IndexTraits> class WellState;
template<typename Scalar, typename IndexTraits> class WellGroupHelpers;

/// Class for handling constraints for the blackoil well model.
template<typename Scalar, typename IndexTraits>
class BlackoilWellModelConstraints
{
public:

    constexpr static int waterPhaseIdx = IndexTraits::waterPhaseIdx;
    constexpr static int oilPhaseIdx = IndexTraits::oilPhaseIdx;
    constexpr static int gasPhaseIdx = IndexTraits::gasPhaseIdx;

    using WellGroupHelpersType = WellGroupHelpers<Scalar, IndexTraits>;

    //! \brief Constructor initializes reference to the well model.
    explicit BlackoilWellModelConstraints(const BlackoilWellModelGeneric<Scalar, IndexTraits>& wellModel)
        : wellModel_(wellModel)
    {}

    //! \brief Check the constraints of a well group.
    bool checkGroupConstraints(const Group& group,
                               const int reportStepIdx,
                               DeferredLogger& deferred_logger) const;

    //! \brief Execute action for broken constraint for an injection well group.
    void actionOnBrokenConstraints(const Group& group,
                                   const Group::InjectionCMode& newControl,
                                   const Phase& controlPhase,
                                   GroupState<Scalar>& group_state,
                                   DeferredLogger& deferred_logger) const;

    //! \brief Execute action on broken constraint for a production well group. Return true if a group control is changed
    bool actionOnBrokenConstraints(const Group& group,
                                   const int reportStepIdx,
                                   const Group::GroupLimitAction group_limit_action,
                                   const Group::ProductionCMode& newControl,
                                   const WellState<Scalar, IndexTraits>& well_state,
                                   std::optional<std::string>& worst_offending_well,
                                   GroupState<Scalar>& group_state,
                                   DeferredLogger& deferred_logger) const;

    //! \brief Update the individual controls for wells in a group. Return true if a group control is changed
    bool updateGroupIndividualControl(const Group& group,
                                      const int reportStepIdx,
                                      const int max_number_of_group_switch,
                                      const bool update_group_switching_log,
                                      std::map<std::string, std::array<std::vector<Group::InjectionCMode>, 3>>& switched_inj,
                                      std::map<std::string, std::vector<Group::ProductionCMode>>& switched_prod,
                                      std::map<std::string, std::pair<std::string, std::string>>& closed_offending_wells,
                                      GroupState<Scalar>& group_state,
                                      WellState<Scalar, IndexTraits>& well_state,
                                      DeferredLogger& deferred_logger) const;

private:
    //! \brief Check and return value and type of constraints for an injection well group.
    std::pair<Group::InjectionCMode, Scalar>
    checkGroupInjectionConstraints(const Group& group,
                                   const int reportStepIdx,
                                   const Phase& phase) const;

    //! \brief Check and return value and type of constraints for a production well group.
    std::pair<Group::ProductionCMode, Scalar>
    checkGroupProductionConstraints(const Group& group,
                                    const int reportStepIdx,
                                    DeferredLogger& deferred_logger) const;

    const BlackoilWellModelGeneric<Scalar, IndexTraits>& wellModel_; //!< Reference to well model
};

} // namespace Opm

#endif
