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

#ifndef OPM_NEWTONITERATIONBLACKOILSIMPLE_HEADER_INCLUDED
#define OPM_NEWTONITERATIONBLACKOILSIMPLE_HEADER_INCLUDED


#include <opm/autodiff/NewtonIterationBlackoilInterface.hpp>
#include <opm/core/linalg/LinearSolverInterface.hpp>

namespace Opm
{

    /// This class solves the fully implicit black-oil system by
    /// simply concatenating the Jacobian matrices and passing the
    /// resulting system to a linear solver. The linear solver used
    /// can be passed in as a constructor argument.
    class NewtonIterationBlackoilSimple : public NewtonIterationBlackoilInterface
    {
    public:
        /// Construct a system solver.
        /// \param[in] linsolver   linear solver to use
        NewtonIterationBlackoilSimple(const LinearSolverInterface& linsolver);

        /// Solve the system of linear equations Ax = b, with A being the
        /// combined derivative matrix of the residual and b
        /// being the residual itself.
        /// \param[in] residual   residual object containing A and b.
        /// \return               the solution x
        virtual SolutionVector linearSolve(const LinearisedBlackoilResidual& residual) const;

    private:
        const LinearSolverInterface& linsolver_;
    };

} // namespace Opm


#endif // OPM_NEWTONITERATIONBLACKOILSIMPLE_HEADER_INCLUDED
