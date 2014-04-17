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

#include <opm/parser/eclipse/EclipseState/Schedule/Well.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/Group.hpp>

#include <boost/lexical_cast.hpp>

#include <memory>

namespace Opm
{
    void WellCollection::addField(GroupConstPtr fieldGroup, size_t timeStep, const PhaseUsage& phaseUsage) {
        WellsGroupInterface* fieldNode = findNode(fieldGroup->name());
        if (fieldNode) {
            OPM_THROW(std::runtime_error, "Trying to add FIELD node, but this already exists. Can only have one FIELD node.");
        }

        roots_.push_back(createGroupWellsGroup(fieldGroup, timeStep, phaseUsage));
    }

    void WellCollection::addGroup(GroupConstPtr groupChild, std::string parent_name,
                                  size_t timeStep, const PhaseUsage& phaseUsage) {
        WellsGroupInterface* parent = findNode(parent_name);
        if (!parent) {
            OPM_THROW(std::runtime_error, "Trying to add child group to group named " << parent_name << ", but this does not exist in the WellCollection.");
        }

        if (findNode(groupChild->name())) {
            OPM_THROW(std::runtime_error, "Trying to add child group named " << groupChild->name() << ", but this group is already in the WellCollection.");

        }

        std::shared_ptr<WellsGroupInterface> child = createGroupWellsGroup(groupChild, timeStep, phaseUsage);

        WellsGroup* parent_as_group = static_cast<WellsGroup*> (parent);
        if (!parent_as_group) {
            OPM_THROW(std::runtime_error, "Trying to add child group to group named " << parent->name() << ", but it's not a group.");
        }
        parent_as_group->addChild(child);
        child->setParent(parent);
    }

    void WellCollection::addWell(WellConstPtr wellChild, size_t timeStep, const PhaseUsage& phaseUsage) {
        WellsGroupInterface* parent = findNode(wellChild->getGroupName(timeStep));
        if (!parent) {
            OPM_THROW(std::runtime_error, "Trying to add well " << wellChild->name() << " Step: " << boost::lexical_cast<std::string>(timeStep) << " to group named " << wellChild->getGroupName(timeStep) << ", but this group does not exist in the WellCollection.");
        }

        std::shared_ptr<WellsGroupInterface> child = createWellWellsGroup(wellChild, timeStep, phaseUsage);

        WellsGroup* parent_as_group = static_cast<WellsGroup*> (parent);
        if (!parent_as_group) {
            OPM_THROW(std::runtime_error, "Trying to add well to group named " << wellChild->getGroupName(timeStep) << ", but it's not a group.");
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
}
