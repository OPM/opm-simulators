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
    /// The approach is similar to the one described in
    /// "Preconditioning for Efficiently Applying Algebraic Multigrid
    /// in Fully Implicit Reservoir Simulations" by Gries et al (SPE 163608).
    class NewtonIterationBlackoilCPR : public NewtonIterationBlackoilInterface
    {
    public:
        /// Construct a system solver.
        /// \param[in] param   parameters controlling the behaviour of
        ///                    the preconditioning and choice of
        ///                    linear solvers.
        ///                    Parameters:
        ///                        cpr_relax        (default 1.0) relaxation for the preconditioner
        ///                        cpr_ilu_n        (default 0) use ILU(n) for preconditioning of the linear system
        ///                        cpr_use_amg      (default false) if true, use AMG preconditioner for elliptic part
        ///                        cpr_use_bicgstab (default true)  if true, use BiCGStab (else use CG) for elliptic part
        NewtonIterationBlackoilCPR(const parameter::ParameterGroup& param);

        /// Solve the system of linear equations Ax = b, with A being the
        /// combined derivative matrix of the residual and b
        /// being the residual itself.
        /// \param[in] residual   residual object containing A and b.
        /// \return               the solution x
        virtual SolutionVector computeNewtonIncrement(const LinearisedBlackoilResidual& residual) const;

        /// \copydoc NewtonIterationBlackoilInterface::iterations
        virtual int iterations () const { return iterations_; }
    private:
        mutable int iterations_;
        double cpr_relax_;
        unsigned int cpr_ilu_n_;
        bool cpr_use_amg_;
        bool cpr_use_bicgstab_;
    };

} // namespace Opm

#endif // OPM_NEWTONITERATIONBLACKOILCPR_HEADER_INCLUDED
