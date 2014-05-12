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

#include <config.h>

#include <opm/autodiff/NewtonIterationBlackoilSimple.hpp>
#include <opm/autodiff/AutoDiffHelpers.hpp>
#include <opm/core/utility/ErrorMacros.hpp>
#include <opm/core/linalg/LinearSolverFactory.hpp>

namespace Opm
{

    /// Construct a system solver.
    /// \param[in] linsolver   linear solver to use
    NewtonIterationBlackoilSimple::NewtonIterationBlackoilSimple(const parameter::ParameterGroup& param)
    {
        linsolver_.reset(new LinearSolverFactory(param));
    }

    /// Solve the linear system Ax = b, with A being the
    /// combined derivative matrix of the residual and b
    /// being the residual itself.
    /// \param[in] residual   residual object containing A and b.
    /// \return               the solution x
    NewtonIterationBlackoilSimple::SolutionVector
    NewtonIterationBlackoilSimple::computeNewtonIncrement(const LinearisedBlackoilResidual& residual) const
    {
        typedef LinearisedBlackoilResidual::ADB ADB;
        const int np = residual.material_balance_eq.size();
        ADB mass_res = residual.material_balance_eq[0];
        for (int phase = 1; phase < np; ++phase) {
            mass_res = vertcat(mass_res, residual.material_balance_eq[phase]);
        }
        const ADB well_res = vertcat(residual.well_flux_eq, residual.well_eq);
        const ADB total_residual = collapseJacs(vertcat(mass_res, well_res));

        const Eigen::SparseMatrix<double, Eigen::RowMajor> matr = total_residual.derivative()[0];

        SolutionVector dx(SolutionVector::Zero(total_residual.size()));
        Opm::LinearSolverInterface::LinearSolverReport rep
            = linsolver_->solve(matr.rows(), matr.nonZeros(),
                                matr.outerIndexPtr(), matr.innerIndexPtr(), matr.valuePtr(),
                                total_residual.value().data(), dx.data());
        if (!rep.converged) {
            OPM_THROW(std::runtime_error,
                      "FullyImplicitBlackoilSolver::solveJacobianSystem(): "
                      "Linear solver convergence failure.");
        }
        return dx;
    }

} // namespace Opm

