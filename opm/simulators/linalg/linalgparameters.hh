// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*
  This file is part of the Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 2 of the License, or
  (at your option) any later version.

  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.

  Consult the COPYING file in the top-level source directory of this
  module for the precise wording of the license and the list of
  copyright holders.
*/
/*!
 * \file
 * \ingroup BlackOilModel
 *
 * \brief Declares the parameters for the black oil model.
 */
#ifndef EWOMS_LINALG_PARAMETERS_HH
#define EWOMS_LINALG_PARAMETERS_HH

namespace Opm::Parameters {

//! number of iterations between solver restarts for the GMRES solver
struct GMResRestart { static constexpr int value = 10; };

/*!
 * \brief Maximum accepted error of the norm of the residual.
 */
template<class Scalar>
struct LinearSolverAbsTolerance { static constexpr Scalar value = -1.0; };

template<class Scalar>
struct LinearSolverMaxError { static constexpr Scalar value = 1e7; };

//! Maximum number of iterations eyecuted by the linear solver
struct LinearSolverMaxIterations { static constexpr int value = 1000; };

/*!
 * \brief The size of the algebraic overlap of the linear solver.
 *
 * Algebraic overlaps can be thought as being the same as the overlap
 * of a grid, but it is only existant for the linear system of
 * equations.
 */
struct LinearSolverOverlapSize { static constexpr unsigned value = 2; };

/*!
 * \brief Maximum accepted error of the solution of the linear solver.
 */
template<class Scalar>
struct LinearSolverTolerance { static constexpr Scalar value = 1e-3; };

/*!
 * \brief Specifies the verbosity of the linear solver
 *
 * By default it is 0, i.e. it doesn't print anything. Setting this
 * property to 1 prints aggregated convergence rates, 2 prints the
 * convergence rate of every iteration of the scheme.
 */
struct LinearSolverVerbosity { static constexpr int value = 0; };

//! The order of the sequential preconditioner
struct PreconditionerOrder { static constexpr int value = 0; };

//! The relaxation factor of the preconditioner
template<class Scalar>
struct PreconditionerRelaxation { static constexpr Scalar value = 1.0; };

} // namespace Opm::Parameters

#endif
