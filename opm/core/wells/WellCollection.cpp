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
#include <boost/shared_ptr.hpp>

namespace Opm
{

    void WellCollection::addChild(const std::string& child_name,
                                  const std::string& parent_name,
                                  const EclipseGridParser& deck)
    {   
        WellsGroupInterface* parent = findNode(parent_name);
        if (!parent) {
            roots_.push_back(createWellsGroup(parent_name, deck));
            parent = roots_[roots_.size() - 1].get();
        }
        boost::shared_ptr<WellsGroupInterface> child;

        for (size_t i = 0; i < roots_.size(); ++i) {
            if (roots_[i]->name() == child_name) {
                child = roots_[i];
                // We've found a new parent to the previously thought root, need to remove it
                for(size_t j = i; j < roots_.size() - 1; ++j) {
                    roots_[j] = roots_[j+1];
                }
                
                roots_.resize(roots_.size()-1);
                break;
            }
        }
        if (!child.get()) {
            child = createWellsGroup(child_name, deck);
        }
        
        WellsGroup* parent_as_group = static_cast<WellsGroup*> (parent);
        if (!parent_as_group) {
            THROW("Trying to add child to group named " << parent_name << ", but it's not a group.");
        }
        parent_as_group->addChild(child);

        if(child->isLeafNode()) {
            leaf_nodes_.push_back(static_cast<WellNode*>(child.get()));
        }
        
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

    void WellCollection::addChild(boost::shared_ptr<WellsGroupInterface>& child_node,
                                  const std::string& parent_name)
    {
        WellsGroupInterface* parent = findNode(parent_name);
        if (parent == NULL) {
            THROW("Parent with name = " << parent_name << " not found.");
        }
        ASSERT(!parent->isLeafNode());
        static_cast<WellsGroup*>(parent)->addChild(child_node);
        if (child_node->isLeafNode()) {
            leaf_nodes_.push_back(static_cast<WellNode*>(child_node.get()));
        }

    }

    /// Adds the node to the collection (as a root node)

    void WellCollection::addChild(boost::shared_ptr<WellsGroupInterface>& child_node)
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
