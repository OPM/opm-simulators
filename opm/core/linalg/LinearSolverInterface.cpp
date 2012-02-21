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


#include <opm/core/linalg/LinearSolverInterface.hpp>
#include <opm/core/linalg/sparse_sys.h>
#include <opm/core/linalg/call_umfpack.h>

namespace Opm
{

    LinearSolverInterface::~LinearSolverInterface()
    {
    }




    LinearSolverInterface::LinearSolverReport
    LinearSolverInterface::solve(const CSRMatrix* A,
                                 const double* rhs,
                                 double* solution) const
    {
        return solve(A->m, A->nnz, A->ia, A->ja, A->sa, rhs, solution);
    }

} // namespace Opm

