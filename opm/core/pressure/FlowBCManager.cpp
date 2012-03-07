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


#include <opm/core/pressure/FlowBCManager.hpp>
#include <opm/core/utility/ErrorMacros.hpp>

namespace Opm
{



    /// Default constructor sets up empty boundary conditions.
    /// By convention, this is equivalent to all-noflow conditions.
    FlowBCManager::FlowBCManager()
	: bc_(0)
    {
	bc_ = flow_conditions_construct(0);
	if (!bc_) {
	    THROW("Failed to construct FlowBoundaryConditions struct.");
	}
    }


    /// Destructor.
    FlowBCManager::~FlowBCManager()
    {
	flow_conditions_destroy(bc_);
    }


    /// Remove all appended BCs.
    /// By convention, BCs are now equivalent to all-noflow conditions.
    void FlowBCManager::clear()
    {
	flow_conditions_clear(bc_);
    }


    /// Append a single boundary condition.
    /// If the type is BC_NOFLOW the value argument is not used.
    /// If the type is BC_PRESSURE the value argument is a pressure value.
    /// If the type is BC_FLUX_TOTVOL the value argument is a total flux value (m^3/s).
    /// Note: unset boundary conditions are noflow by convention,
    /// so it is normally not necessary to explicitly append
    /// BC_NOFLOW conditions. However, it may make sense to do so
    /// if the bc will change during a simulation run.
    /// Note: if normal velocity bcs are desired, convert to
    /// fluxes by multiplying with face area.
    void FlowBCManager::append(const FlowBCType type,
			       const int face,
			       const double value)
    {
	int ok = flow_conditions_append(type, face, value, bc_);
	if (!ok) {
	    THROW("Failed to append boundary condition for face " << face);
	}
    }


    /// Access the managed boundary conditions.
    /// The method is named similarly to c_str() in std::string,
    /// to make it clear that we are returning a C-compatible struct.
    const FlowBoundaryConditions* FlowBCManager::c_bcs() const
    {
	return bc_;
    }



} // namespace Opm
