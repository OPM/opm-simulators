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

#include "config.h"
#include <opm/core/linalg/LinearSolverUmfpack.hpp>
#include <opm/core/linalg/sparse_sys.h>
#include <opm/core/linalg/call_umfpack.h>

namespace Opm
{

    LinearSolverUmfpack::LinearSolverUmfpack()
    {
    }




    LinearSolverUmfpack::~LinearSolverUmfpack()
    {
    }




    LinearSolverInterface::LinearSolverReport
    LinearSolverUmfpack::solve(const int size,
                               const int nonzeros,
                               const int* ia,
                               const int* ja,
                               const double* sa,
                               const double* rhs,
                               double* solution) const
    {
        CSRMatrix A  = {
            (size_t)size,
            (size_t)nonzeros,
            const_cast<int*>(ia),
            const_cast<int*>(ja),
            const_cast<double*>(sa)
        };
        call_UMFPACK(&A, rhs, solution);
        LinearSolverReport rep;
        rep.converged = true;
        return rep;
    }

    void LinearSolverUmfpack::setTolerance(const double /*tol*/)
    {
    }

    double LinearSolverUmfpack::getTolerance() const
    {
        return -1.;
    }


} // namespace Opm

