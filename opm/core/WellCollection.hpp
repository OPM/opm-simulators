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


#ifndef OPM_WELLCOLLECTION_HPP
#define	OPM_WELLCOLLECTION_HPP

#include <vector>
#include <opm/core/WellsGroup.hpp>
#include <opm/core/grid.h>
#include <opm/core/eclipse/EclipseGridParser.hpp>

namespace Opm
{

    class WellCollection
    {
    public:
        WellCollection();
        virtual ~WellCollection();

        void addChild(std::string child, std::string parent,
                const EclipseGridParser& deck);

        
        bool conditionsMet(const std::vector<double>& well_bhp, const std::vector<double>& well_rate, 
                           const UnstructuredGrid& grid, const std::vector<double>& saturations, double epsilon=1e-8) const;
        
        const std::vector<std::tr1::shared_ptr<WellsGroupInterface> >& getLeafNodes() const;
        
        void calculateGuideRates();
    private:
        // To account for the possibility of a forest
        std::vector<std::tr1::shared_ptr<WellsGroupInterface> > roots_;
        
        // This will be used to traverse the bottom nodes.
        std::vector<std::tr1::shared_ptr<WellsGroupInterface> > leaf_nodes_;
        
        WellsGroupInterface* findNode(std::string name);
    };

} // namespace Opm
#endif	/* OPM_WELLCOLLECTION_HPP */

