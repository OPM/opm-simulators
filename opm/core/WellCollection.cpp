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

    WellCollection::WellCollection()
    {
    }

    WellCollection::~WellCollection()
    {
    }

    void WellCollection::addChild(std::string child_name, std::string parent_name,
            const EclipseGridParser& deck)
    {   
        WellsGroupInterface* parent = findNode(parent_name);
        if (!parent) {
            roots_.push_back(createWellsGroup(parent_name, deck));
            parent = roots_[roots_.size() - 1].get();
        }
        std::tr1::shared_ptr<WellsGroupInterface> child;

        for (int i = 0; i < roots_.size(); ++i) {
            if (roots_[i]->name() == child_name) {
                child = roots_[i];
                // We've found a new parent to the previously thought root, need to remove it
                for(int j = i; j < roots_.size() - 1; ++j) {
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
            leaf_nodes_.push_back(child);
        }
    }
    
    const std::vector<std::tr1::shared_ptr<WellsGroupInterface> >& WellCollection::getLeafNodes() const {
        return leaf_nodes_;
    }

    WellsGroupInterface* WellCollection::findNode(std::string name)
    {

        for (int i = 0; i < roots_.size(); i++) {
            WellsGroupInterface* result = roots_[i]->findGroup(name);
            if (result) {
                return result;
            }
        }
        return NULL;
    }
}