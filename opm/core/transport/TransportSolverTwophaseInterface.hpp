/*
  Copyright 2012, 2013 SINTEF ICT, Applied Mathematics.

  This file is part of the Open Porous Media Project (OPM).

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

#ifndef OPM_TRANSPORTSOLVERTWOPHASEINTERFACE_HEADER_INCLUDED
#define OPM_TRANSPORTSOLVERTWOPHASEINTERFACE_HEADER_INCLUDED

#include <opm/core/simulator/TwophaseState.hpp>
#include <opm/core/simulator/WellState.hpp>

namespace Opm
{

    /// Base class for two-phase incompressible transport solvers.
    class TransportSolverTwophaseInterface
    {
    public:
        /// Virtual destructor to enable inheritance.
        virtual ~TransportSolverTwophaseInterface();

        /// Solve for saturation at next timestep.
        /// \param[in]      porevolume   Array of pore volumes.
        /// \param[in]      source       Transport source term. For interpretation see Opm::computeTransportSource().
        /// \param[in]      dt           Time step.
        /// \param[in, out] state        Reservoir state. Calling solve() will read state.faceflux() and
        ///                              read and write state.saturation().
        virtual void solve(const double* porevolume,
                           const double* source,
                           const double dt,
                           TwophaseState& state) = 0;
    };

}

#endif // OPM_TRANSPORTSOLVERTWOPHASEINTERFACE_HEADER_INCLUDED
