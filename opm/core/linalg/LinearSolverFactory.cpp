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

#if HAVE_CONFIG_H
#include "config.h"
#endif

#include <opm/core/linalg/LinearSolverFactory.hpp>

#if HAVE_SUITESPARSE_UMFPACK_H
#include <opm/core/linalg/LinearSolverUmfpack.hpp>
#endif

#if HAVE_DUNE_ISTL
#include <opm/core/linalg/LinearSolverIstl.hpp>
#endif

#include <opm/core/utility/parameters/ParameterGroup.hpp>
#include <opm/core/utility/ErrorMacros.hpp>
#include <string>


namespace Opm
{

    LinearSolverFactory::LinearSolverFactory()
    {
#if HAVE_SUITESPARSE_UMFPACK_H
        solver_.reset(new LinearSolverUmfpack);
#elif HAVE_DUNE_ISTL
        solver_.reset(new LinearSolverIstl);
#else
        THROW("No linear solver available, you must have UMFPACK or dune-istl installed to use LinearSolverFactory.");
#endif
    }




    LinearSolverFactory::LinearSolverFactory(const parameter::ParameterGroup& param)
    {
        const std::string ls =
            param.getDefault<std::string>("linsolver", "umfpack");

        if (ls == "umfpack") {
#if HAVE_SUITESPARSE_UMFPACK_H
            solver_.reset(new LinearSolverUmfpack);
#endif
        }

        else if (ls == "istl") {
#if HAVE_DUNE_ISTL
            solver_.reset(new LinearSolverIstl(param));
#endif
        }

        else {
            THROW("Linear solver " << ls << " is unknown.");
        }

        if (! solver_) {
            THROW("Linear solver " << ls << " is not enabled in "
                  "this configuration.");
        }
    }




    LinearSolverFactory::~LinearSolverFactory()
    {
    }




    LinearSolverInterface::LinearSolverReport
    LinearSolverFactory::solve(const int size,
                               const int nonzeros,
                               const int* ia,
                               const int* ja,
                               const double* sa,
                               const double* rhs,
                               double* solution) const
    {
        return solver_->solve(size, nonzeros, ia, ja, sa, rhs, solution);
    }

    void LinearSolverFactory::setTolerance(const double tol)
    {
        solver_->setTolerance(tol);
    }

    double LinearSolverFactory::getTolerance() const
    {
        return solver_->getTolerance();
    }



} // namespace Opm

