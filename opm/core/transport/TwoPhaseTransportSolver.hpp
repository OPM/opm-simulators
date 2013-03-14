/*===========================================================================
//
// File: TwoPhaseTransportSolver.hpp
//
// Author: hnil <hnil@sintef.no>
//
// Created: 9 Nov 2012
//==========================================================================*/
/*
  Copyright 2011 SINTEF ICT, Applied Mathematics.
  Copyright 2011 Statoil ASA.
  
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


#ifndef TWOPHASETRANSPORTSOLVER_HPP
#define TWOPHASETRANSPORTSOLVER_HPP

#include <opm/core/simulator/TwophaseState.hpp>
#include <opm/core/simulator/WellState.hpp>
namespace Opm
{

    /// Base class for tranport solvers
    class TwoPhaseTransportSolver
    {
    public:
        /// Virtual destructor to enable inheritance.
        virtual ~TwoPhaseTransportSolver() {}

        /// Solve for saturation at next timestep.
        /// \param[in] darcyflux         Array of signed face fluxes.
        /// \param[in] porevolume        Array of pore volumes.
        /// \param[in] source            Transport source term.
        /// \param[in] dt                Time step.
        /// \param[in, out] saturation   Phase saturations.
        virtual void solve(const double* porevolume,
                           const double* source,
                           const double dt,
                           TwophaseState& state,
                           WellState& wstate) = 0;
    };

}

#endif // TWOPHASETRANSPORTSOLVER_HPP
