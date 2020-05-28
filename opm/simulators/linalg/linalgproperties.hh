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
 * \brief Declares the properties required by the black oil model.
 */
#ifndef EWOMS_LINALG_PROPERTIES_HH
#define EWOMS_LINALG_PROPERTIES_HH

#include <opm/models/utils/basicproperties.hh>

BEGIN_PROPERTIES

//! The type of the linear solver to be used
NEW_PROP_TAG(LinearSolverBackend);

//! the preconditioner used by the linear solver
NEW_PROP_TAG(PreconditionerWrapper);


//! The floating point type used internally by the linear solver
NEW_PROP_TAG(LinearSolverScalar);

/*!
 * \brief The size of the algebraic overlap of the linear solver.
 *
 * Algebraic overlaps can be thought as being the same as the overlap
 * of a grid, but it is only existant for the linear system of
 * equations.
 */
NEW_PROP_TAG(LinearSolverOverlapSize);

/*!
 * \brief Maximum accepted error of the solution of the linear solver.
 */
NEW_PROP_TAG(LinearSolverTolerance);

/*!
 * \brief Maximum accepted error of the norm of the residual.
 */
NEW_PROP_TAG(LinearSolverAbsTolerance);

/*!
 * \brief Specifies the verbosity of the linear solver
 *
 * By default it is 0, i.e. it doesn't print anything. Setting this
 * property to 1 prints aggregated convergence rates, 2 prints the
 * convergence rate of every iteration of the scheme.
 */
NEW_PROP_TAG(LinearSolverVerbosity);

//! Maximum number of iterations eyecuted by the linear solver
NEW_PROP_TAG(LinearSolverMaxIterations);

//! The order of the sequential preconditioner
NEW_PROP_TAG(PreconditionerOrder);

//! The relaxation factor of the preconditioner
NEW_PROP_TAG(PreconditionerRelaxation);

//! number of iterations between solver restarts for the GMRES solver
NEW_PROP_TAG(GMResRestart);

//! The class that allows to manipulate sparse matrices
NEW_PROP_TAG(SparseMatrixAdapter);

//! Vector containing a quantity of for equation for each DOF of the whole grid
NEW_PROP_TAG(GlobalEqVector);

NEW_PROP_TAG(AmgCoarsenTarget);
NEW_PROP_TAG(LinearSolverMaxError);
NEW_PROP_TAG(LinearSolverWrapper);
NEW_PROP_TAG(Overlap);
NEW_PROP_TAG(OverlappingLinearOperator);
NEW_PROP_TAG(OverlappingMatrix);
NEW_PROP_TAG(OverlappingScalarProduct);
NEW_PROP_TAG(OverlappingVector);

END_PROPERTIES

#endif
