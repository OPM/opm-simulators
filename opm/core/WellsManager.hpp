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
#include <opm/core/WellCollection.hpp>
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

	/// Construct from input deck and grid.
	/// The permeability argument may be zero if the input contain
	/// well productivity indices, otherwise it must be given in
	/// order to approximate these by the Peaceman formula.
	WellsManager(const Opm::EclipseGridParser& deck,
		     const UnstructuredGrid& grid,
		     const double* permeability);

	/// Destructor.
	~WellsManager();

	/// Access the managed Wells.
	/// The method is named similarly to c_str() in std::string,
	/// to make it clear that we are returning a C-compatible struct.
	const Wells* c_wells() const;

    private:
	// Disable copying and assignment.
	WellsManager(const WellsManager& other);
	WellsManager& operator=(const WellsManager& other);
	// The managed Wells.
	Wells* w_;
        
        WellCollection well_collection_;

    };

} // namespace Opm


#endif // OPM_WELLSMANAGER_HEADER_INCLUDED
