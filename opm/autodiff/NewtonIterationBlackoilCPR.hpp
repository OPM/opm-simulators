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

#ifndef OPM_NEWTONITERATIONBLACKOILCPR_HEADER_INCLUDED
#define OPM_NEWTONITERATIONBLACKOILCPR_HEADER_INCLUDED


#include <opm/autodiff/NewtonIterationBlackoilInterface.hpp>
#include <opm/core/utility/parameters/ParameterGroup.hpp>
#include <opm/core/linalg/LinearSolverInterface.hpp>
#include <memory>

namespace Opm
{

    /// This class solves the fully implicit black-oil system by
    /// applying a Constrained Pressure Residual preconditioning
    /// strategy.
    class NewtonIterationBlackoilCPR : public NewtonIterationBlackoilInterface
    {
    public:
        /// Construct a system solver.
        /// \param[in] param   parameters controlling the behaviour of
        ///                    the preconditioning and choice of
        ///                    linear solvers.
        NewtonIterationBlackoilCPR(const parameter::ParameterGroup& param);

        /// Solve the system of linear equations Ax = b, with A being the
        /// combined derivative matrix of the residual and b
        /// being the residual itself.
        /// \param[in] residual   residual object containing A and b.
        /// \return               the solution x
        virtual SolutionVector computeNewtonIncrement(const LinearisedBlackoilResidual& residual) const;

    private:
        std::unique_ptr<LinearSolverInterface> linsolver_elliptic_;
        std::unique_ptr<LinearSolverInterface> linsolver_full_;
    };

} // namespace Opm

#endif // OPM_NEWTONITERATIONBLACKOILCPR_HEADER_INCLUDED
