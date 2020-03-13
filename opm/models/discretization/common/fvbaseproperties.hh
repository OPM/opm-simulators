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
#ifndef EWOMS_FV_BASE_PROPERTIES_HH
#define EWOMS_FV_BASE_PROPERTIES_HH

#include "fvbasenewtonmethod.hh"
#include "fvbaseproperties.hh"
#include "fvbasefdlocallinearizer.hh"

#include <opm/models/utils/basicproperties.hh>
#include <opm/models/io/vtkprimaryvarsmodule.hh>
#include <opm/simulators/linalg/parallelbicgstabbackend.hh>

BEGIN_PROPERTIES

//! The type tag for models based on the finite volume schemes
NEW_TYPE_TAG(FvBaseDiscretization,
             INHERITS_FROM(ImplicitModel,
                           FvBaseNewtonMethod,
                           VtkPrimaryVars));


//! set the splices for the finite volume discretizations


//SET_SPLICES(FvBaseDiscretization, LinearSolverSplice, LocalLinearizerSplice);
template<class TypeTag>
struct Splices<TypeTag, TTag::FvBaseDiscretization>
{
//     using type = std::tuple<GetSplicePropType<TypeTag, Properties::LinearSolverSplice>,
//                             GetSplicePropType<TypeTag, Properties::LocalLinearizerSplice>>;
    using type = std::tuple<TTag::ParallelBiCGStabLinearSolver,
                            TTag::FiniteDifferenceLocalLinearizer>;
};

//! use a parallel BiCGStab linear solver by default
SET_TAG_PROP(FvBaseDiscretization, LinearSolverSplice, ParallelBiCGStabLinearSolver);

//! by default, use finite differences to linearize the system of PDEs
SET_TAG_PROP(FvBaseDiscretization, LocalLinearizerSplice, FiniteDifferenceLocalLinearizer);

/*!
 * \brief Representation of a function evaluation and all necessary derivatives with
 *        regard to the intensive quantities of the primary variables.
 *
 * Depending on the chosen linearization method, this property may be the same as the
 * "Scalar" property (if the finite difference linearizer is used), or it may be more
 * complex (for the linearizer which uses automatic differentiation).
 */

//! The type of the DUNE grid
//! The type of the grid view

//! The class describing the stencil of the spatial discretization

//! The class describing the discrete function space when dune-fem is used, otherwise it points to the stencil class

//! The type of the problem
//! The type of the base class for all problems which use this model
//! The type of the model
//! Number of equations in the system of PDEs

//! The type of the spatial discretization used by the model
//! The discretization specific part of the local residual
//! The type of the local residual function
//! The type of the local linearizer
//! Specify if elements that do not belong to the local process' grid partition should be
//! skipped

//! Linearizes the global non-linear system of equations
//! The class that allows to manipulate sparse matrices

//! A vector of holding a quantity for each equation (usually at a given spatial location)
//! A vector of holding a quantity for each equation for each DOF of an element
//! Vector containing a quantity of for equation for each DOF of the whole grid

//! Vector containing volumetric or areal rates of quantities
//! Type of object for specifying boundary conditions
//! The class which represents a constraint degree of freedom

//! Vector containing all primary variables of the grid

//! A vector of primary variables within a sub-control volume
//! The secondary variables within a sub-control volume
//! The discretization specific part of the intensive quantities

//! The secondary variables of all degrees of freedom in an element's stencil
//! The secondary variables of a boundary segment
//! The secondary variables of a constraint degree of freedom
//! Data required to calculate a flux over a face
//! Calculates gradients of arbitrary quantities at flux integration points

//! The part of the intensive quantities which is specific to the spatial discretization

//! The part of the extensive quantities which is specific to the spatial discretization

//! The part of the VTK ouput modules which is specific to the spatial discretization

//! The class to create grid communication handles

/*!
 * \brief The OpenMP threads manager
 */

//! use locking to prevent race conditions when linearizing the global system of
//! equations in multi-threaded mode. (setting this property to true is always save, but
//! it may slightly deter performance in multi-threaded simlations and some
//! discretizations do not need this.)

// high-level simulation control

//! Manages the simulation time

/*!
 * \brief Switch to enable or disable grid adaptation
 *
 * Currently grid adaptation requires the presence of the dune-FEM module. If it is not
 * available and grid adaptation is enabled, an exception is thrown.
 */

/*!
 * \brief The directory to which simulation output ought to be written to.
 */

/*!
 * \brief Global switch to enable or disable the writing of VTK output files
 *
 * If writing VTK files is disabled, then the WriteVtk$FOO options do
 * not have any effect...
 */

/*!
 * \brief Determines if the VTK output is written to disk asynchronously
 *
 * I.e. written to disk using a separate thread. This has only an effect if
 * EnableVtkOutput is true and if the simulation is run sequentially. The reasons for
 * this not being used for MPI-parallel simulations are that Dune's VTK output code does
 * not support multi-threaded multi-process VTK output and even if it would, the result
 * would be slower than when using synchronous output.
 */

/*!
 * \brief Specify the format the VTK output is written to disk
 *
 * Possible values are:
 *   - Dune::VTK::ascii (default)
 *   - Dune::VTK::base64
 *   - Dune::VTK::appendedraw
 *   - Dune::VTK::appendedbase64
 */

//! Specify whether the some degrees of fredom can be constraint

/*!
 * \brief Specify the maximum size of a time integration [s].
 *
 * The default is to not limit the step size.
 */

/*!
 * \brief Specify the minimal size of a time integration [s].
 *
 * The default is to not limit the step size.
 */

/*!
 * \brief The maximum allowed number of timestep divisions for the
 *        Newton solver.
 */

/*!
 * \brief Continue with a non-converged solution instead of giving up
 *        if we encounter a time step size smaller than the minimum time
 *        step size.
 */

/*!
 * \brief Specify whether all intensive quantities for the grid should be
 *        cached in the discretization.
 *
 * This potentially reduces the CPU time, but comes at the cost of
 * higher memory consumption. In turn, the higher memory requirements
 * may cause the simulation to exhibit worse cache coherence behavior
 * which eats some of the computational benefits again.
 */

/*!
 * \brief Specify whether the storage terms for previous solutions should be cached.
 *
 * This potentially reduces the CPU time, but comes at the cost of higher memory
 * consumption.
 */

/*!
 * \brief Specify whether to use the already calculated solutions as
 *        starting values of the intensive quantities.
 *
 * This only makes sense if the calculation of the intensive quantities is
 * very expensive (e.g. for non-linear fugacity functions where the
 * solver converges faster).
 */

// mappers from local to global DOF indices

/*!
 * \brief The mapper to find the global index of a vertex.
 */

/*!
 * \brief The mapper to find the global index of an element.
 */

/*!
 * \brief The mapper to find the global index of a degree of freedom.
 */

/*!
 * \brief The class which marks the border indices associated with the
 *        degrees of freedom on a process boundary.
 *
 * This is required for the algebraic overlap stuff.
 */

/*!
 * \brief The history size required by the time discretization
 */

/*!
 * \brief Specify whether the storage terms use extensive quantities or not.
 *
 * Most models don't need this, but the (Navier-)Stokes ones do...
 */

//! \brief Specify whether to use volumetric residuals or not


//! Specify if experimental features should be enabled or not.

END_PROPERTIES

#endif
