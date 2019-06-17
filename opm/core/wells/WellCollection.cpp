/*
  Copyright 2011 SINTEF ICT, Applied Mathematics.

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

#include "config.h"
#include <opm/core/wells/WellCollection.hpp>

#include <opm/parser/eclipse/EclipseState/Schedule/Well/Well2.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/Group.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/SummaryState.hpp>

#include <boost/lexical_cast.hpp>

#include <memory>

namespace Opm
{
    void WellCollection::addField(const Group& fieldGroup, size_t timeStep, const PhaseUsage& phaseUsage) {
        WellsGroupInterface* fieldNode = findNode(fieldGroup.name());
        if (fieldNode) {
            OPM_THROW(std::runtime_error, "Trying to add FIELD node, but this already exists. Can only have one FIELD node.");
        }

        roots_.push_back(createGroupWellsGroup(fieldGroup, timeStep, phaseUsage));
    }

    void WellCollection::addGroup(const Group& groupChild, std::string parent_name,
                                  size_t timeStep, const PhaseUsage& phaseUsage) {
        WellsGroupInterface* parent = findNode(parent_name);
        if (!parent) {
            OPM_THROW(std::runtime_error, "Trying to add child group to group named " << parent_name << ", but this does not exist in the WellCollection.");
        }

        if (findNode(groupChild.name())) {
            OPM_THROW(std::runtime_error, "Trying to add child group named " << groupChild.name() << ", but this group is already in the WellCollection.");

        }

        if (groupChild.isProductionGroup(timeStep) || groupChild.isInjectionGroup(timeStep)) {
            group_control_active_ = true;
        }

        std::shared_ptr<WellsGroupInterface> child = createGroupWellsGroup(groupChild, timeStep, phaseUsage);

        if (child->injSpec().control_mode_ == InjectionSpecification::VREP) {
            having_vrep_groups_ = true;
        }

        WellsGroup* parent_as_group = static_cast<WellsGroup*> (parent);
        if (!parent_as_group) {
            OPM_THROW(std::runtime_error, "Trying to add child group to group named " << parent->name() << ", but it's not a group.");
        }
        parent_as_group->addChild(child);
        child->setParent(parent);
    }

    void WellCollection::addWell(const Well2& wellChild, const SummaryState& summaryState, size_t timeStep, const PhaseUsage& phaseUsage) {
        if (wellChild.getStatus() == WellCommon::SHUT) {
            //SHUT wells are not added to the well collection
            return;
        }

        WellsGroupInterface* parent = findNode(wellChild.groupName());
        if (!parent) {
            OPM_THROW(std::runtime_error, "Trying to add well " << wellChild.name() << " Step: " << boost::lexical_cast<std::string>(timeStep) << " to group named " << wellChild.groupName() << ", but this group does not exist in the WellCollection.");
        }

        std::shared_ptr<WellsGroupInterface> child = createWellWellsGroup(wellChild, summaryState, timeStep, phaseUsage);

        WellsGroup* parent_as_group = static_cast<WellsGroup*> (parent);
        if (!parent_as_group) {
            OPM_THROW(std::runtime_error, "Trying to add well to group named " << wellChild.groupName() << ", but it's not a group.");
        }
        parent_as_group->addChild(child);

        leaf_nodes_.push_back(static_cast<WellNode*>(child.get()));

        child->setParent(parent);
    }

    const std::vector<WellNode*>& WellCollection::getLeafNodes() const {
        return leaf_nodes_;
    }

    WellsGroupInterface* WellCollection::findNode(const std::string& name)
    {

        for (size_t i = 0; i < roots_.size(); i++) {
            WellsGroupInterface* result = roots_[i]->findGroup(name);
            if (result) {
                return result;
            }
        }
        return NULL;
    }

    const WellsGroupInterface* WellCollection::findNode(const std::string& name) const
    {

        for (size_t i = 0; i < roots_.size(); i++) {
            WellsGroupInterface* result = roots_[i]->findGroup(name);
            if (result) {
                return result;
            }
        }
        return NULL;
    }


    WellNode& WellCollection::findWellNode(const std::string& name) const
    {
        auto well_node = std::find_if(leaf_nodes_.begin(), leaf_nodes_.end(),
                    [&] ( WellNode* w) {
                    return w->name() == name;
                    });

        // Does not find the well
        if (well_node == leaf_nodes_.end()) {
            OPM_THROW(std::runtime_error, "Could not find well " << name << " in the well collection!\n");
        }

        return *(*well_node);
    }

    /// Adds the child to the collection
    /// and appends it to parent's children.
    /// \param[in] child   the child node
    /// \param[in] parent  name of parent node

    void WellCollection::addChild(std::shared_ptr<WellsGroupInterface>& child_node,
                                  const std::string& parent_name)
    {
        WellsGroupInterface* parent = findNode(parent_name);
        if (parent == NULL) {
            OPM_THROW(std::runtime_error, "Parent with name = " << parent_name << " not found.");
        }
        assert(!parent->isLeafNode());
        static_cast<WellsGroup*>(parent)->addChild(child_node);
        if (child_node->isLeafNode()) {
            leaf_nodes_.push_back(static_cast<WellNode*>(child_node.get()));
        }

    }

    /// Adds the node to the collection (as a root node)

    void WellCollection::addChild(std::shared_ptr<WellsGroupInterface>& child_node)
    {
        roots_.push_back(child_node);
        if (child_node->isLeafNode()) {
            leaf_nodes_.push_back(static_cast<WellNode*> (child_node.get()));
        }
    }

    bool WellCollection::conditionsMet(const std::vector<double>& well_bhp,
                                       const std::vector<double>& well_reservoirrates_phase,
                                       const std::vector<double>& well_surfacerates_phase)
    {
        for (size_t i = 0; i < roots_.size(); i++) {
            WellPhasesSummed phases;
            if (!roots_[i]->conditionsMet(well_bhp,
                                          well_reservoirrates_phase,
                                          well_surfacerates_phase,
                                          phases)) {
                return false;
            }
        }
        return true;
    }

    void WellCollection::setWellsPointer(Wells* wells) {
        for(size_t i = 0; i < leaf_nodes_.size(); i++) {
            leaf_nodes_[i]->setWellsPointer(wells, i);
        }
    }

    void WellCollection::applyGroupControls()
    {
        for (size_t i = 0; i < roots_.size(); ++i) {
            roots_[i]->applyProdGroupControls();
            roots_[i]->applyInjGroupControls();
        }

        group_control_applied_ = true;
    }

    /// Applies explicit reinjection controls. This must be called at each timestep to be correct.
    /// \param[in]    well_reservoirrates_phase
    ///                         A vector containing reservoir rates by phase for each well.
    ///                         Is assumed to be ordered the same way as the related Wells-struct,
    ///                         with all phase rates of a single well adjacent in the array.
    /// \param[in]    well_surfacerates_phase
    ///                         A vector containing surface rates by phase for each well.
    ///                         Is assumed to be ordered the same way as the related Wells-struct,
    ///                         with all phase rates of a single well adjacent in the array.

    void WellCollection::applyExplicitReinjectionControls(const std::vector<double>& well_reservoirrates_phase,
                                                          const std::vector<double>& well_surfacerates_phase)
    {
        for (size_t i = 0; i < roots_.size(); ++i) {
            roots_[i]->applyExplicitReinjectionControls(well_reservoirrates_phase, well_surfacerates_phase);
        }
    }


    void WellCollection::applyVREPGroupControls(const std::vector<double>& well_voidage_rates,
                                                const std::vector<double>& conversion_coeffs)
    {
        for (size_t i = 0; i < roots_.size(); ++i) {
            roots_[i]->applyVREPGroupControls(well_voidage_rates, conversion_coeffs);
        }
    }


    // TODO: later, it should be extended to update group targets
    bool WellCollection::needUpdateWellTargets() const
    {
        return needUpdateInjectionTargets() || needUpdateProductionTargets();
    }


    bool WellCollection::needUpdateInjectionTargets() const
    {
        // TODO: it should based on individual group
        // With current approach, it will potentially result in more update,
        // thus more iterations, while it will not cause result wrong.
        // If the group control and individual control is mixed, then it need to
        // update the well targets
        bool any_group_control_node = false;
        bool any_individual_control_node = false;

        for (size_t i = 0; i < leaf_nodes_.size(); ++i) {
            if (leaf_nodes_[i]->isInjector()) {
                if (leaf_nodes_[i]->individualControl()) {
                    any_individual_control_node = true;
                } else {
                    any_group_control_node = true;
                }
            }
        }

        return (any_group_control_node && any_individual_control_node);
    }


    // These two functions should be made one
    bool WellCollection::needUpdateProductionTargets() const
    {
        // TODO: it should based on individual group
        // With current approach, it will potentially result in more update,
        // thus more iterations, while it will not cause result wrong.
        // If the group control and individual control is mixed, then it need to
        // update the well targets
        bool any_group_control_node = false;
        bool any_individual_control_node = false;

        for (size_t i = 0; i < leaf_nodes_.size(); ++i) {
            if (leaf_nodes_[i]->isProducer()) {
                if (leaf_nodes_[i]->individualControl()) {
                    any_individual_control_node = true;
                } else {
                    any_group_control_node = true;
                }
            }
        }

        return (any_group_control_node && any_individual_control_node);
    }


    void WellCollection::updateWellTargets(const std::vector<double>& well_rates)
    {

        // TODO: if it gets converged, should we still update targets?

        // set the target_updated to be false
        for (WellNode* well_node : leaf_nodes_) {
            well_node->setTargetUpdated(false);
        }

        // TODO: currently, we only handle the level of the well groups for the moment, i.e. the level just above wells
        // We believe the relations between groups are similar to the relations between different wells inside the same group.
        // While there will be somre more complication invloved for sure.
        for (size_t i = 0; i < leaf_nodes_.size(); ++i) {
            // find a node needs to update targets, then update targets for all the wellls inside the group.
            if (!leaf_nodes_[i]->targetUpdated()) {
                WellsGroupInterface* parent_node = leaf_nodes_[i]->getParent();
                // update the target within this group.
                if (leaf_nodes_[i]->isProducer()) {
                    if (parent_node->prodSpec().control_mode_ == ProductionSpecification::NONE) {
                        continue;
                    }
                    parent_node->updateWellProductionTargets(well_rates);
                }

                if (leaf_nodes_[i]->isInjector()) {
                    if (parent_node->injSpec().control_mode_ == InjectionSpecification::NONE) {
                        continue;
                    }
                    parent_node->updateWellInjectionTargets(well_rates);
                }
            }
        }
    }


    bool WellCollection::havingVREPGroups() const
    {
        return having_vrep_groups_;
    }


    bool WellCollection::groupControlActive() const
    {
        return group_control_active_;
    }



    bool WellCollection::groupControlApplied() const
    {
        return group_control_applied_;
    }


    bool WellCollection::groupTargetConverged(const std::vector<double>& well_rates) const
    {
        // TODO: eventually, there should be only one root node
        // TODO: we also need to check the injection target, while we have not done that.
        for (const std::shared_ptr<WellsGroupInterface>& root_node : roots_) {
            if ( !root_node->groupProdTargetConverged(well_rates) ) {
                return false;
            }
        }
        return true;
    }



    void WellCollection::
    setGuideRatesWithPotentials(const Wells* wells,
                                const PhaseUsage& phase_usage,
                                const std::vector<double>& well_potentials) const
    {
        // TODO: assuming the order of well_potentials is the same with the order in wells struct
        // TODO: it will overwrite the well potentials from other means. It should be changed after
        // fixing the other part of the code. It makes the current flow only support guide rates based on
        // well potentials.
        const int np = wells->number_of_phases;
        const int nw = wells->number_of_wells;

        for (int w = 0; w < nw; ++w) {
            const std::string well_name = wells->name[w];

            WellNode& well_node = findWellNode(well_name);

            const WellType well_type = wells->type[w];
            // TODO: eventually the following standard will be wrong, it will belong to FIELD group
            if (well_node.getParent() != nullptr) { // If it does not belong a group, will it belong to FIELD?
                const WellsGroupInterface* group = well_node.getParent();
                if (well_type == PRODUCER) {
                    // The guide rates is calculated based on the group control
                    // Currently only supporting WRAT, ORAT and GRAT.
                    ProductionSpecification::ControlMode control_mode = group->prodSpec().control_mode_;
                    if (control_mode == ProductionSpecification::FLD) {
                        if (group->getParent() !=  nullptr) {
                            // TODO: only handle one level FLD control
                            const WellsGroupInterface* higher_group = group->getParent();
                            control_mode = higher_group->prodSpec().control_mode_;
                        } else {
                            OPM_THROW(std::runtime_error, "Group " << group->name() << " is under FLD control while no higher level of group is specified.");
                        }
                    }

                    switch (control_mode) {
                    case ProductionSpecification::WRAT: {
                       if (!phase_usage.phase_used[BlackoilPhases::Aqua]) {
                           OPM_THROW(std::runtime_error, "Water phase not used, yet found water rate controlled well.");
                       }
                       const int water_index = phase_usage.phase_pos[BlackoilPhases::Aqua];
                       well_node.prodSpec().guide_rate_ = well_potentials[np * w + water_index];
                       well_node.prodSpec().guide_rate_type_ = ProductionSpecification::WATER;
                       break;
                    }
                    case ProductionSpecification::ORAT: {
                        if (!phase_usage.phase_used[BlackoilPhases::Liquid]) {
                            OPM_THROW(std::runtime_error, "Oil phase not used, yet found oil rate controlled well.");
                        }
                        const int oil_index = phase_usage.phase_pos[BlackoilPhases::Liquid];
                        well_node.prodSpec().guide_rate_ = well_potentials[np * w + oil_index];
                        well_node.prodSpec().guide_rate_type_ = ProductionSpecification::OIL;
                        break;
                    }
                    case ProductionSpecification::GRAT: {
                        if (!phase_usage.phase_used[BlackoilPhases::Vapour]) {
                            OPM_THROW(std::runtime_error, "Gas phase not used, yet found gas rate controlled well.");
                        }
                        const int gas_index = phase_usage.phase_pos[BlackoilPhases::Vapour];
                        well_node.prodSpec().guide_rate_ = well_potentials[np * w + gas_index];
                        well_node.prodSpec().guide_rate_type_ = ProductionSpecification::GAS;
                        break;
                    }
                    case ProductionSpecification::FLD: {
                        OPM_THROW(std::logic_error, "Not support more than one continous level of FLD control");
                    }
                    case ProductionSpecification::LRAT: {
                        double guide_rate = 0;
                        if (phase_usage.phase_used[BlackoilPhases::Liquid]) {
                            const int oil_index = phase_usage.phase_pos[BlackoilPhases::Liquid];
                            const double potential_oil = well_potentials[np * w + oil_index];
                            guide_rate += potential_oil;
                        }
                        if (phase_usage.phase_used[BlackoilPhases::Aqua]) {
                            const int water_index = phase_usage.phase_pos[BlackoilPhases::Aqua];
                            const double potential_water = well_potentials[np * w + water_index];
                            guide_rate += potential_water;
                        }
                        // not sure if no water and no oil, what will happen here, zero guide_rate?
                        well_node.prodSpec().guide_rate_ = guide_rate;
                        well_node.prodSpec().guide_rate_type_ = ProductionSpecification::LIQ;
                        break;
                    }
                    case ProductionSpecification::NONE: {
                        // Group control is not in use for this group.
                        break;
                    }
                    default:
                        OPM_THROW(std::logic_error, "Not supported control_mode for guide rate computed" <<
                                  " from well potentials: " << ProductionSpecification::toString(group->prodSpec().control_mode_) );
                    }

                } else if (well_type == INJECTOR) {
                    // The guide rates is calculated based on the group injector type
                    switch (group->injSpec().injector_type_) {
                    case InjectionSpecification::WATER: {
                        if (!phase_usage.phase_used[BlackoilPhases::Aqua]) {
                            OPM_THROW(std::runtime_error, "Water phase not used, yet found water injecting well.");
                        }
                        const int water_index = phase_usage.phase_pos[BlackoilPhases::Aqua];
                        well_node.injSpec().guide_rate_ = well_potentials[np * w + water_index];
                        // Guide rates applies to the phase that the well is injecting i.e water
                        well_node.injSpec().guide_rate_type_ = InjectionSpecification::RAT;
                        break;
                    }
                    case InjectionSpecification::OIL: {
                        if (!phase_usage.phase_used[BlackoilPhases::Liquid]) {
                            OPM_THROW(std::runtime_error, "Oil phase not used, yet found oil injecting well.");
                        }
                        const int oil_index = phase_usage.phase_pos[BlackoilPhases::Liquid];
                        well_node.injSpec().guide_rate_ = well_potentials[np * w + oil_index];
                        // Guide rates applies to the phase that the well is injecting i.e. oil
                        well_node.injSpec().guide_rate_type_ = InjectionSpecification::RAT;
                        break;
                    }
                    case InjectionSpecification::GAS: {
                        if (!phase_usage.phase_used[BlackoilPhases::Vapour]) {
                            OPM_THROW(std::runtime_error, "Gas phase not used, yet found gas injecting well.");
                        }
                        const int gas_index = phase_usage.phase_pos[BlackoilPhases::Vapour];
                        well_node.injSpec().guide_rate_ = well_potentials[np * w + gas_index];
                        // Guide rates applies to the phase that the well is injecting i.e gas
                        well_node.injSpec().guide_rate_type_ = InjectionSpecification::RAT;
                        break;
                    }
                    default:
                        OPM_THROW(std::logic_error, "Not supported injector type for guide rate computed" <<
                                  " from well potentials: " << InjectionSpecification::toString(group->injSpec().injector_type_) );
                    }
                } else { // neither injector nor producer
                    OPM_THROW(std::logic_error, "Expected well type to be either INJECTOR or PRODUCER for well " << well_node.name() );

                }
            } // end of if (well_node.getParent() != nullptr)
        } // end of for (int w = 0; w < nw; ++w)
    }


    bool WellCollection::requireWellPotentials() const
    {
        for (const auto& well_node : leaf_nodes_) {
            if (well_node->isGuideRateWellPotential()) {
                return true;
            }
        }
        return false;
    }

}
