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
        /// Adds and creates if necessary the child to the collection
        /// and appends it to parent's children. Also adds and creates the parent
        /// if necessary.
        /// \param[in] child   name of child node
        /// \param[in] parent  name of parent node
        /// \param[in] deck    deck from which we will extract group control data
        void addChild(const std::string& child,
                      const std::string& parent,
                      const EclipseGridParser& deck);

        /// Checks if each condition is met, applies well controls where needed
        /// (that is, it either changes the active control of violating wells, or shuts
        /// down wells). Only one change is applied per invocation. Typical use will be
        /// \code
        /// solve_pressure();
        /// while(!collection.conditionsMet(well_bhp, well_rate, summed_phases)) {
        ///     solve_pressure();
        /// }
        /// \endcode
        ///
        /// \note It's highly recommended to use the conditionsMet found in WellsManager.
        /// \param[in]    well_bhp  A vector containing the bhp for each well. Is assumed 
        ///                         to be ordered the same way as the related Wells-struct.
        /// \param[in]    well_rate A vector containing the rate for each well. Is assumed 
        ///                         to be ordered the same way as the related Wells-struct.
        /// \param[in]    epsilon   The error tolerated for each inequality. Formally, it will accept
        ///                         (a - b <= epsilon) as (a <= b).
        /// \return true if no violations were found, false otherwise (false also implies a change).
        bool conditionsMet(const std::vector<double>& well_bhp,
                           const std::vector<double>& well_rate, 
                           const double epsilon=1e-8);
        
        /// Adds the well pointer to each leaf node (does not take ownership).
        void setWellsPointer(Wells* wells);
        
        /// \return A set of pointers to every well in the collection
        const std::vector<WellNode*>& getLeafNodes() const;
        
        /// This will, for each group \f$G\f$ for each well \f$w\in G\f$ calculate 
        /// the new guide rate \f$\tilde{g}_w\f$ for each well as
        /// \f[ \tilde{g}_w := \frac{g_w}{\sum_{w'\in G} g_{w'}}\f]
        void calculateGuideRates();
        
        /// Finds the group with the given name.
        /// \param[in] the name of the group
        /// \return the pointer to the group if found, NULL otherwise
        WellsGroupInterface* findNode(std::string name);
        
        /// Finds the group with the given name.
        /// \param[in] the name of the group
        /// \return the pointer to the group if found, NULL otherwise
        const WellsGroupInterface* findNode(std::string name) const;
        
    private:
        // To account for the possibility of a forest
        std::vector<std::tr1::shared_ptr<WellsGroupInterface> > roots_;
        
        // This will be used to traverse the bottom nodes.
        std::vector<WellNode*> leaf_nodes_;
        
        
    };

} // namespace Opm
#endif	/* OPM_WELLCOLLECTION_HPP */

