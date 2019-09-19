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
NEW_PROP_TAG(LinearSolverSplice);
NEW_PROP_TAG(ParallelBiCGStabLinearSolver);

NEW_PROP_TAG(LocalLinearizerSplice);
NEW_PROP_TAG(FiniteDifferenceLocalLinearizer);

SET_SPLICES(FvBaseDiscretization, LinearSolverSplice, LocalLinearizerSplice);

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
NEW_PROP_TAG(Evaluation);

//! The type of the DUNE grid
NEW_PROP_TAG(Grid);
//! The type of the grid view
NEW_PROP_TAG(GridView);

//! The class describing the stencil of the spatial discretization
NEW_PROP_TAG(Stencil);

//! The class describing the discrete function space when dune-fem is used, otherwise it points to the stencil class
NEW_PROP_TAG(DiscreteFunctionSpace);

//! The type of the problem
NEW_PROP_TAG(Problem);
//! The type of the base class for all problems which use this model
NEW_PROP_TAG(BaseProblem);
//! The type of the model
NEW_PROP_TAG(Model);
//! Number of equations in the system of PDEs
NEW_PROP_TAG(NumEq);

//! The type of the spatial discretization used by the model
NEW_PROP_TAG(Discretization);
//! The discretization specific part of the local residual
NEW_PROP_TAG(DiscLocalResidual);
//! The type of the local residual function
NEW_PROP_TAG(LocalResidual);
//! The type of the local linearizer
NEW_PROP_TAG(LocalLinearizer);
//! Specify if elements that do not belong to the local process' grid partition should be
//! skipped
NEW_PROP_TAG(LinearizeNonLocalElements);

//! Linearizes the global non-linear system of equations
NEW_PROP_TAG(BaseLinearizer);
//! The class that allows to manipulate sparse matrices
NEW_PROP_TAG(SparseMatrixAdapter);

//! A vector of holding a quantity for each equation (usually at a given spatial location)
NEW_PROP_TAG(EqVector);
//! A vector of holding a quantity for each equation for each DOF of an element
NEW_PROP_TAG(ElementEqVector);
//! Vector containing a quantity of for equation for each DOF of the whole grid
NEW_PROP_TAG(GlobalEqVector);

//! Vector containing volumetric or areal rates of quantities
NEW_PROP_TAG(RateVector);
//! Type of object for specifying boundary conditions
NEW_PROP_TAG(BoundaryRateVector);
//! The class which represents a constraint degree of freedom
NEW_PROP_TAG(Constraints);

//! Vector containing all primary variables of the grid
NEW_PROP_TAG(SolutionVector);

//! A vector of primary variables within a sub-control volume
NEW_PROP_TAG(PrimaryVariables);
//! The secondary variables within a sub-control volume
NEW_PROP_TAG(IntensiveQuantities);
//! The discretization specific part of the intensive quantities
NEW_PROP_TAG(DiscIntensiveQuantities);

//! The secondary variables of all degrees of freedom in an element's stencil
NEW_PROP_TAG(ElementContext);
//! The secondary variables of a boundary segment
NEW_PROP_TAG(BoundaryContext);
//! The secondary variables of a constraint degree of freedom
NEW_PROP_TAG(ConstraintsContext);
//! Data required to calculate a flux over a face
NEW_PROP_TAG(ExtensiveQuantities);
//! Calculates gradients of arbitrary quantities at flux integration points
NEW_PROP_TAG(GradientCalculator);

//! The part of the intensive quantities which is specific to the spatial discretization
NEW_PROP_TAG(DiscBaseIntensiveQuantities);

//! The part of the extensive quantities which is specific to the spatial discretization
NEW_PROP_TAG(DiscExtensiveQuantities);

//! The part of the VTK ouput modules which is specific to the spatial discretization
NEW_PROP_TAG(DiscBaseOutputModule);

//! The class to create grid communication handles
NEW_PROP_TAG(GridCommHandleFactory);

/*!
 * \brief The OpenMP threads manager
 */
NEW_PROP_TAG(ThreadManager);
NEW_PROP_TAG(ThreadsPerProcess);

//! use locking to prevent race conditions when linearizing the global system of
//! equations in multi-threaded mode. (setting this property to true is always save, but
//! it may slightly deter performance in multi-threaded simlations and some
//! discretizations do not need this.)
NEW_PROP_TAG(UseLinearizationLock);

// high-level simulation control

//! Manages the simulation time
NEW_PROP_TAG(Simulator);

/*!
 * \brief Switch to enable or disable grid adaptation
 *
 * Currently grid adaptation requires the presence of the dune-FEM module. If it is not
 * available and grid adaptation is enabled, an exception is thrown.
 */
NEW_PROP_TAG(EnableGridAdaptation);

/*!
 * \brief The directory to which simulation output ought to be written to.
 */
NEW_PROP_TAG(OutputDir);

/*!
 * \brief Global switch to enable or disable the writing of VTK output files
 *
 * If writing VTK files is disabled, then the WriteVtk$FOO options do
 * not have any effect...
 */
NEW_PROP_TAG(EnableVtkOutput);

/*!
 * \brief Determines if the VTK output is written to disk asynchronously
 *
 * I.e. written to disk using a separate thread. This has only an effect if
 * EnableVtkOutput is true and if the simulation is run sequentially. The reasons for
 * this not being used for MPI-parallel simulations are that Dune's VTK output code does
 * not support multi-threaded multi-process VTK output and even if it would, the result
 * would be slower than when using synchronous output.
 */
NEW_PROP_TAG(EnableAsyncVtkOutput);

/*!
 * \brief Specify the format the VTK output is written to disk
 *
 * Possible values are:
 *   - Dune::VTK::ascii (default)
 *   - Dune::VTK::base64
 *   - Dune::VTK::appendedraw
 *   - Dune::VTK::appendedbase64
 */
NEW_PROP_TAG(VtkOutputFormat);

//! Specify whether the some degrees of fredom can be constraint
NEW_PROP_TAG(EnableConstraints);

/*!
 * \brief Specify the maximum size of a time integration [s].
 *
 * The default is to not limit the step size.
 */
NEW_PROP_TAG(MaxTimeStepSize);

/*!
 * \brief Specify the minimal size of a time integration [s].
 *
 * The default is to not limit the step size.
 */
NEW_PROP_TAG(MinTimeStepSize);

/*!
 * \brief The maximum allowed number of timestep divisions for the
 *        Newton solver.
 */
NEW_PROP_TAG(MaxTimeStepDivisions);

/*!
 * \brief Continue with a non-converged solution instead of giving up
 *        if we encounter a time step size smaller than the minimum time
 *        step size.
 */
NEW_PROP_TAG(ContinueOnConvergenceError);

/*!
 * \brief Specify whether all intensive quantities for the grid should be
 *        cached in the discretization.
 *
 * This potentially reduces the CPU time, but comes at the cost of
 * higher memory consumption. In turn, the higher memory requirements
 * may cause the simulation to exhibit worse cache coherence behavior
 * which eats some of the computational benefits again.
 */
NEW_PROP_TAG(EnableIntensiveQuantityCache);

/*!
 * \brief Specify whether the storage terms for previous solutions should be cached.
 *
 * This potentially reduces the CPU time, but comes at the cost of higher memory
 * consumption.
 */
NEW_PROP_TAG(EnableStorageCache);

/*!
 * \brief Specify whether to use the already calculated solutions as
 *        starting values of the intensive quantities.
 *
 * This only makes sense if the calculation of the intensive quantities is
 * very expensive (e.g. for non-linear fugacity functions where the
 * solver converges faster).
 */
NEW_PROP_TAG(EnableThermodynamicHints);

// mappers from local to global DOF indices

/*!
 * \brief The mapper to find the global index of a vertex.
 */
NEW_PROP_TAG(VertexMapper);

/*!
 * \brief The mapper to find the global index of an element.
 */
NEW_PROP_TAG(ElementMapper);

/*!
 * \brief The mapper to find the global index of a degree of freedom.
 */
NEW_PROP_TAG(DofMapper);

/*!
 * \brief The class which marks the border indices associated with the
 *        degrees of freedom on a process boundary.
 *
 * This is required for the algebraic overlap stuff.
 */
NEW_PROP_TAG(BorderListCreator);

/*!
 * \brief The history size required by the time discretization
 */
NEW_PROP_TAG(TimeDiscHistorySize);

/*!
 * \brief Specify whether the storage terms use extensive quantities or not.
 *
 * Most models don't need this, but the (Navier-)Stokes ones do...
 */
NEW_PROP_TAG(ExtensiveStorageTerm);

//! \brief Specify whether to use volumetric residuals or not
NEW_PROP_TAG(UseVolumetricResidual);


//! Specify if experimental features should be enabled or not.
NEW_PROP_TAG(EnableExperiments);

END_PROPERTIES

#endif
