/*
  Copyright 2012 SINTEF ICT, Applied Mathematics.

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

#ifndef OPM_WELLSMANAGER_HEADER_INCLUDED
#define OPM_WELLSMANAGER_HEADER_INCLUDED


#include <opm/parser/eclipse/EclipseState/EclipseState.hpp>

#include <opm/core/wells/WellCollection.hpp>
#include <opm/core/wells/WellsGroup.hpp>

struct Wells;
struct UnstructuredGrid;


namespace Opm
{

    class EclipseGridParser;

    /// This class manages a Wells struct in the sense that it
    /// encapsulates creation and destruction of the wells
    /// data structure.
    /// The resulting Wells is available through the c_wells() method.
    class WellsManager
    {
    public:
	/// Default constructor -- no wells.
    WellsManager();

    /// Construct from existing wells object.
    /// WellsManager is not properly initialised in the sense that the logic to
    /// manage control switching does not exist.
    ///
    /// @param[in] W Existing wells object.
    WellsManager(struct Wells* W);

	/// Construct from input deck and grid.
	/// The permeability argument may be zero if the input contain
	/// well productivity indices, otherwise it must be given in
	/// order to approximate these by the Peaceman formula.
	WellsManager(const Opm::EclipseGridParser& deck,
		     const UnstructuredGrid& grid,
             const double* permeability);


    WellsManager(const Opm::EclipseStateConstPtr eclipseState,
                 const size_t timeStep,
                 const Opm::EclipseGridParser& deck,
                 const UnstructuredGrid& grid,
                 const double* permeability);

	/// Destructor.
	~WellsManager();

        /// Does the "deck" define any wells?
        bool empty() const;

	/// Access the managed Wells.
	/// The method is named similarly to c_str() in std::string,
	/// to make it clear that we are returning a C-compatible struct.
	const Wells* c_wells() const;

        /// Access the well group hierarchy.
        const WellCollection& wellCollection() const;

        /// Checks if each condition is met, applies well controls where needed
        /// (that is, it either changes the active control of violating wells, or shuts
        /// down wells). Only one change is applied per invocation. Typical use will be
        /// \code
        /// solve_pressure();
        /// while(!wells.conditionsMet(...)) {
        ///     solve_pressure();
        /// }
        /// \endcode
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

    private:
	// Disable copying and assignment.
	WellsManager(const WellsManager& other);
	WellsManager& operator=(const WellsManager& other);

	// Data
	Wells* w_;
        WellCollection well_collection_;




    };

} // namespace Opm


#endif // OPM_WELLSMANAGER_HEADER_INCLUDED
