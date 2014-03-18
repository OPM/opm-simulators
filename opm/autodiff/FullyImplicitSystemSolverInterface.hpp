/*
  Copyright 2014 SINTEF ICT, Applied Mathematics.

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

#ifndef OPM_FULLYIMPLICITSYSTEMSOLVERINTERFACE_HEADER_INCLUDED
#define OPM_FULLYIMPLICITSYSTEMSOLVERINTERFACE_HEADER_INCLUDED

#include <opm/autodiff/FullyImplicitBlackoilResidual.hpp>

namespace Opm
{

    /// Interface class for (linear) solvers for the fully implicit black-oil system.
    class FullyImplicitSystemSolverInterface
    {
    public:
        /// Return type for linearSolve(). A simple, non-ad vector type.
        typedef FullyImplicitBlackoilResidual::ADB::V SolutionVector;

        /// Solve the linear system Ax = b, with A being the
        /// combined derivative matrix of the residual and b
        /// being the residual itself.
        /// \param[in] residual   residual object containing A and b.
        /// \return               the solution x
        virtual SolutionVector linearSolve(const FullyImplicitBlackoilResidual& residual) const = 0;
    };

} // namespace Opm


#endif // OPM_FULLYIMPLICITSYSTEMSOLVERINTERFACE_HEADER_INCLUDED
