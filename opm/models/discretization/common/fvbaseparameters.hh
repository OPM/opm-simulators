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
 * \ingroup FiniteVolumeDiscretizations
 *
 * \brief Declare the properties used by the infrastructure code of
 *        the finite volume discretizations.
 */
#ifndef EWOMS_FV_BASE_PARAMETERS_HH
#define EWOMS_FV_BASE_PARAMETERS_HH

#include <limits>

namespace Opm::Parameters {

/*!
 * \brief Continue with a non-converged solution instead of giving up
 *        if we encounter a time step size smaller than the minimum time
 *        step size.
 */
struct ContinueOnConvergenceError { static constexpr bool value = false; };

/*!
 * \brief Determines if the VTK output is written to disk asynchronously
 *
 * I.e. written to disk using a separate thread. This has only an effect if
 * EnableVtkOutput is true and if the simulation is run sequentially. The reasons for
 * this not being used for MPI-parallel simulations are that Dune's VTK output code does
 * not support multi-threaded multi-process VTK output and even if it would, the result
 * would be slower than when using synchronous output.
 */
struct EnableAsyncVtkOutput { static constexpr bool value = true; };

/*!
 * \brief Switch to enable or disable grid adaptation
 *
 * Currently grid adaptation requires the presence of the dune-FEM module. If it is not
 * available and grid adaptation is enabled, an exception is thrown.
 */
struct EnableGridAdaptation { static constexpr bool value = false; };

/*!
 * \brief Specify whether all intensive quantities for the grid should be
 *        cached in the discretization.
 *
 * This potentially reduces the CPU time, but comes at the cost of
 * higher memory consumption. In turn, the higher memory requirements
 * may cause the simulation to exhibit worse cache coherence behavior
 * which eats some of the computational benefits again.
 */
struct EnableIntensiveQuantityCache { static constexpr bool value = false; };

/*!
 * \brief Specify whether the storage terms for previous solutions should be cached.
 *
 * This potentially reduces the CPU time, but comes at the cost of higher memory
 * consumption.
 */
struct EnableStorageCache { static constexpr bool value = false; };

/*!
 * \brief Specify whether to use the already calculated solutions as
 *        starting values of the intensive quantities.
 *
 * This only makes sense if the calculation of the intensive quantities is
 * very expensive (e.g. for non-linear fugacity functions where the
 * solver converges faster).
 */
struct EnableThermodynamicHints { static constexpr bool value = false; };

/*!
 * \brief Global switch to enable or disable the writing of VTK output files
 *
 * If writing VTK files is disabled, then the WriteVtk$FOO options do
 * not have any effect...
 */
struct EnableVtkOutput { static constexpr bool value = true; };

/*!
 * \brief Specify the maximum size of a time integration [s].
 *
 * The default is to not limit the step size.
 */
template<class Scalar>
struct MaxTimeStepSize { static constexpr Scalar value = std::numeric_limits<Scalar>::infinity(); };

/*!
 * \brief The maximum allowed number of timestep divisions for the
 *        Newton solver.
 */
struct MaxTimeStepDivisions { static constexpr unsigned value = 10; };

/*!
 * \brief Specify the minimal size of a time integration [s].
 *
 * The default is to not limit the step size.
 */
template<class Scalar>
struct MinTimeStepSize { static constexpr Scalar value = 0.0; };

/*!
 * \brief The directory to which simulation output ought to be written to.
 */
struct OutputDir { static constexpr auto value = ""; };

//! \brief Number of threads per process.
struct ThreadsPerProcess { static constexpr int value = 1; };

} // namespace Opm::Parameters

#endif
