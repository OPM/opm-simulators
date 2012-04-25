/*
Copyright 2011 SINTEF ICT, Applied Mathematics.


This file is part of The Open Reservoir Simulator Project (OpenRS).

OpenRS is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

OpenRS is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with OpenRS.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "WellCollection.hpp"

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
        std::tr1::shared_ptr<WellsGroupInterface> child;

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

    WellsGroupInterface* WellCollection::findNode(std::string name)
    {

        for (size_t i = 0; i < roots_.size(); i++) {
            WellsGroupInterface* result = roots_[i]->findGroup(name);
            if (result) {
                return result;
            }
        }
        return NULL;
    }
    
    const WellsGroupInterface* WellCollection::findNode(std::string name) const
    {

        for (size_t i = 0; i < roots_.size(); i++) {
            WellsGroupInterface* result = roots_[i]->findGroup(name);
            if (result) {
                return result;
            }
        }
        return NULL;
    }
    
    void WellCollection::conditionsMet(const std::vector<double>& well_bhp,
                                       const std::vector<double>& well_rate, 
                                       const UnstructuredGrid& grid,
                                       WellControlResult& result,
                                       double epsilon) const
    {
        for (size_t i = 0; i < leaf_nodes_.size(); i++) {
            leaf_nodes_[i]->conditionsMet(well_bhp, well_rate, grid, result, epsilon);
        }        
    }

    void WellCollection::calculateGuideRates()
    {
        for(size_t i = 0; i < roots_.size(); i++) {
            roots_[i]->calculateGuideRates();
        }
    }
    
    void WellCollection::setWellsPointer(const Wells* wells) {
        for(size_t i = 0; i < leaf_nodes_.size(); i++) {
            leaf_nodes_[i]->setWellsPointer(wells, i);
        }
    }
}
