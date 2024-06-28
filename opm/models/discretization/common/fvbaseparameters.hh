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

#include <opm/models/utils/propertysystem.hh>

namespace Opm::Parameters {

template<class TypeTag, class MyTypeTag>
struct ThreadsPerProcess { using type = Properties::UndefinedProperty; };

/*!
 * \brief Switch to enable or disable grid adaptation
 *
 * Currently grid adaptation requires the presence of the dune-FEM module. If it is not
 * available and grid adaptation is enabled, an exception is thrown.
 */
template<class TypeTag, class MyTypeTag>
struct EnableGridAdaptation { using type = Properties::UndefinedProperty; };

/*!
 * \brief The directory to which simulation output ought to be written to.
 */
template<class TypeTag, class MyTypeTag>
struct OutputDir { using type = Properties::UndefinedProperty; };

/*!
 * \brief Global switch to enable or disable the writing of VTK output files
 *
 * If writing VTK files is disabled, then the WriteVtk$FOO options do
 * not have any effect...
 */
template<class TypeTag, class MyTypeTag>
struct EnableVtkOutput { using type = Properties::UndefinedProperty; };

/*!
 * \brief Determines if the VTK output is written to disk asynchronously
 *
 * I.e. written to disk using a separate thread. This has only an effect if
 * EnableVtkOutput is true and if the simulation is run sequentially. The reasons for
 * this not being used for MPI-parallel simulations are that Dune's VTK output code does
 * not support multi-threaded multi-process VTK output and even if it would, the result
 * would be slower than when using synchronous output.
 */
template<class TypeTag, class MyTypeTag>
struct EnableAsyncVtkOutput { using type = Properties::UndefinedProperty; };

/*!
 * \brief Specify the maximum size of a time integration [s].
 *
 * The default is to not limit the step size.
 */
template<class TypeTag, class MyTypeTag>
struct MaxTimeStepSize { using type = Properties::UndefinedProperty; };

/*!
 * \brief Specify the minimal size of a time integration [s].
 *
 * The default is to not limit the step size.
 */
template<class TypeTag, class MyTypeTag>
struct MinTimeStepSize { using type = Properties::UndefinedProperty; };

/*!
 * \brief The maximum allowed number of timestep divisions for the
 *        Newton solver.
 */
template<class TypeTag, class MyTypeTag>
struct MaxTimeStepDivisions { using type = Properties::UndefinedProperty; };

/*!
 * \brief Continue with a non-converged solution instead of giving up
 *        if we encounter a time step size smaller than the minimum time
 *        step size.
 */
template<class TypeTag, class MyTypeTag>
struct ContinueOnConvergenceError { using type = Properties::UndefinedProperty; };

} // namespace Opm::Parameters

#endif
