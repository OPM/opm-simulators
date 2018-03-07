/*
  Copyright 2018 SINTEF Digital, Mathematics and Cybernetics.
  Copyright 2018 Statoil Petroleum AS.

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

namespace Opm
{

    static void LinearSolverAmgcl::solve(const int sz,
                                         const std::vector<int>& ptr,
                                         const std::vector<int>& col,
                                         const std::vector<double>& val,
                                         int& iters,
                                         double& error)
    {
        // Call amgcl to solve system
        typedef amgcl::backend::builtin<double> Backend;
        typedef amgcl::make_solver<
            // Use AMG as preconditioner:
            amgcl::amg<
                Backend,
                amgcl::coarsening::smoothed_aggregation,
                amgcl::relaxation::spai0
                >,
            // And BiCGStab as iterative solver:
            amgcl::solver::bicgstab<Backend>
            > Solver;
        Solver solve( boost::tie(sz, ptr, col, val) );
        std::vector<double> sol(sz, 0.0);
        int    iters;
        double error;
        boost::tie(iters, error) = solve(rhs, sol);
    }

} // namespace Opm
