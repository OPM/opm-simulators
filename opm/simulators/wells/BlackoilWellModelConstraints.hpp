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

#include <utility>

namespace Opm {

class BlackoilWellModelGeneric;
class DeferredLogger;
class GroupState;
class SummaryState;
class WellState;

/// Class for handling constraints for the blackoil well model.
class BlackoilWellModelConstraints
{
public:
    //! \brief Constructor initializes reference to the well model.
    BlackoilWellModelConstraints(const BlackoilWellModelGeneric& wellModel)
        : wellModel_(wellModel)
    {}

    /// Return true if any well has a THP constraint.
    bool hasTHPConstraints() const;

    //! \brief Check the constraints of a well group.
    bool checkGroupConstraints(const Group& group,
                               const int reportStepIdx,
                               DeferredLogger& deferred_logger) const;

    //! \brief Execute action for broken constraint for an injection well group.
    void actionOnBrokenConstraints(const Group& group,
                                   const Group::InjectionCMode& newControl,
                                   const Phase& controlPhase,
                                   GroupState& group_state,
                                   DeferredLogger& deferred_logger) const;

    //! \brief Execute action on broken constraint for a production well group.
    void actionOnBrokenConstraints(const Group& group,
                                   const Group::GroupLimitAction group_limit_action,
                                   const Group::ProductionCMode& newControl,
                                   GroupState& group_state,
                                   DeferredLogger& deferred_logger) const;

    //! \brief Update the individual controls for wells in a group.
    bool updateGroupIndividualControl(const Group& group,
                                      const int reportStepIdx,
                                      std::map<std::pair<std::string,Opm::Phase>,std::string>& switched_inj,
                                      std::map<std::string, std::string>& switched_prod,
                                      GroupState& group_state,
                                      WellState& well_state,
                                      DeferredLogger& deferred_logger) const;

private:
    //! \brief Check and return value and type of constraints for an injection well group.
    std::pair<Group::InjectionCMode, double>
    checkGroupInjectionConstraints(const Group& group,
                                   const int reportStepIdx,
                                   const Phase& phase) const;

    //! \brief Check and return value and type of constraints for a production well group.
    std::pair<Group::ProductionCMode, double>
    checkGroupProductionConstraints(const Group& group,
                                    const int reportStepIdx,
                                    DeferredLogger& deferred_logger) const;

    const BlackoilWellModelGeneric& wellModel_; //!< Reference to well model
};

} // namespace Opm

#endif
