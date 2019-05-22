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
#ifndef OPM_WELLCOLLECTION_HPP
#define	OPM_WELLCOLLECTION_HPP

#include <vector>
#include <memory>

#include <opm/core/wells/WellsGroup.hpp>
#include <opm/grid/UnstructuredGrid.h>
#include <opm/core/props/phaseUsageFromDeck.hpp>

#include <opm/parser/eclipse/EclipseState/Schedule/Well/Well2.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/Group.hpp>

namespace Opm
{

    class WellCollection
    {
    public:

        void addField(const Group& fieldGroup, size_t timeStep, const PhaseUsage& phaseUsage);

        void addWell(const Well2& wellChild, size_t timeStep, const PhaseUsage& phaseUsage);

        void addGroup(const Group& groupChild, std::string parent_name,
                      size_t timeStep, const PhaseUsage& phaseUsage);

        /// Adds the child to the collection
        /// and appends it to parent's children.
        /// \param[in] child   the child node
        /// \param[in] parent  name of parent node
        void addChild(std::shared_ptr<WellsGroupInterface>& child_node,
                      const std::string& parent);

        /// Adds the node to the collection (as a root node)
        void addChild(std::shared_ptr<WellsGroupInterface>& child_node);

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
        /// \param[in]    well_reservoirrates_phase
        ///                         A vector containing reservoir rates by phase for each well.
        ///                         Is assumed to be ordered the same way as the related Wells-struct,
        ///                         with all phase rates of a single well adjacent in the array.
        /// \param[in]    well_surfacerates_phase
        ///                         A vector containing surface rates by phase for each well.
        ///                         Is assumed to be ordered the same way as the related Wells-struct,
        ///                         with all phase rates of a single well adjacent in the array.
        /// \return true if no violations were found, false otherwise (false also implies a change).
        bool conditionsMet(const std::vector<double>& well_bhp,
                           const std::vector<double>& well_reservoirrates_phase,
                           const std::vector<double>& well_surfacerates_phase);

        /// Adds the well pointer to each leaf node (does not take ownership).
        void setWellsPointer(Wells* wells);

        /// \return A set of pointers to every well in the collection
        const std::vector<WellNode*>& getLeafNodes() const;

        /// Finds the group with the given name.
        /// \param[in] the name of the group
        /// \return the pointer to the group if found, NULL otherwise
        WellsGroupInterface* findNode(const std::string& name);

        /// Finds the group with the given name.
        /// \param[in] the name of the group
        /// \return the pointer to the group if found, NULL otherwise
        const WellsGroupInterface* findNode(const std::string& name) const;


        WellNode& findWellNode(const std::string& name) const;


        /// Applies all group controls (injection and production)
        void applyGroupControls();

        /// Applies explicit reinjection controls. This must be called at each timestep to be correct.
        /// \param[in]    well_reservoirrates_phase
        ///                         A vector containing reservoir rates by phase for each well.
        ///                         Is assumed to be ordered the same way as the related Wells-struct,
        ///                         with all phase rates of a single well adjacent in the array.
        /// \param[in]    well_surfacerates_phase
        ///                         A vector containing surface rates by phase for each well.
        ///                         Is assumed to be ordered the same way as the related Wells-struct,
        ///                         with all phase rates of a single well adjacent in the array.
        void applyExplicitReinjectionControls(const std::vector<double>& well_reservoirrates_phase,
                                              const std::vector<double>& well_surfacerates_phase);


        /// applying VREP group control based on calculated voidage rates
        void applyVREPGroupControls(const std::vector<double>& well_voidage_rates,
                                    const std::vector<double>& conversion_coeffs);

        /// Checking whether need to update the targets of the wells / or the groups later
        /// True  need to update well targets within this iteration, no switching control within this iteration.
        /// False no need to update well targets within this iteration, continuing as usual.
        bool needUpdateWellTargets() const;

        /// Checking whether need to update the targets for the injection wells.
        bool needUpdateInjectionTargets() const;

        /// Checking whehter need to update the targets for the production wells.
        bool needUpdateProductionTargets() const;

        /// Updating the well targets based on the well rates.
        void updateWellTargets(const std::vector<double>& well_rates);

        /// When we have VREP group, we need to update the targets based on the updated production voidage rates for each iteration.
        bool havingVREPGroups() const;

        /// Whether we have active group control
        bool groupControlActive() const;

        /// Whether we have applied the group control
        bool groupControlApplied() const;

        /// Whether the group target is converged
        // It is considered converged if eitehr the group targets are matched or the group targets are not matched while the wells are
        // running under their own limits so that they can not produce more
        // It is considered not converged if the group targets are not matched while some of the wells are still running under group control
        // The strategy may need to be adjusted when more complicated multi-layered group control situation applied, not sure about thatyet.
        bool groupTargetConverged(const std::vector<double>& well_rates) const;


        /// Setting the guide rates with well potentials
        void setGuideRatesWithPotentials(const Wells* wells,
                                         const PhaseUsage& phase_usage,
                                         const std::vector<double>& well_potentials) const;


        bool requireWellPotentials() const;

    private:
        // To account for the possibility of a forest
        std::vector<std::shared_ptr<WellsGroupInterface> > roots_;

        // This will be used to traverse the bottom nodes.
        std::vector<WellNode*> leaf_nodes_;

        bool having_vrep_groups_ = false;

        bool group_control_active_ = false;

        // This is used to mark whether apply or update the group control
        bool group_control_applied_ = false;
    };

} // namespace Opm

#endif	/* OPM_WELLCOLLECTION_HPP */

