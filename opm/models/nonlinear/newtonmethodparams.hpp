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
#ifndef OPM_NEWTON_METHOD_PARAMS_HPP
#define OPM_NEWTON_METHOD_PARAMS_HPP

namespace Opm::Parameters {

//! The maximum error which may occur in a simulation before the
//! Newton method for the time step is aborted
template<class Scalar>
struct NewtonMaxError { static constexpr Scalar value = 1e100; };

//! Number of maximum iterations for the Newton method.
struct NewtonMaxIterations { static constexpr int value = 20; };

/*!
 * \brief The number of iterations at which the Newton method
 *        should aim at.
 *
 * This is used to control the time-step size. The heuristic used
 * is to scale the last time-step size by the deviation of the
 * number of iterations used from the target steps.
 */
struct NewtonTargetIterations { static constexpr int value = 10; };

/*!
 * \brief The value for the error below which convergence is declared
 *
 * This value can (and for the porous media models will) be changed to account for grid
 * scaling and other effects.
 */
template<class Scalar>
struct NewtonTolerance { static constexpr Scalar value = 1e-8; };

//! Specifies whether the Newton method should print messages or not
struct NewtonVerbose { static constexpr bool value = true; };

//! Specifies whether the convergence rate and the global residual
//! gets written out to disk for every Newton iteration
struct NewtonWriteConvergence { static constexpr bool value = false; };

} // end namespace Opm::Parameters

#endif // OPM_NEWTON_METHOD_PARAMS_HPP
