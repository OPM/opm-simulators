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

#include <opm/models/utils/basicproperties.hh>
#include <opm/models/io/dgfvanguard.hh>
#include <opm/simulators/linalg/parallelbicgstabbackend.hh>

namespace Opm::Properties {

namespace TTag {
struct FvBaseNewtonMethod;
struct VtkPrimaryVars;
struct FiniteDifferenceLocalLinearizer;
}

namespace TTag {

//! The type tag for models based on the finite volume schemes
struct FvBaseDiscretization
{ using InheritsFrom = std::tuple<VtkPrimaryVars, FvBaseNewtonMethod, ImplicitModel>; };

} // namespace TTag


//! set the splices for the finite volume discretizations
template<class TypeTag, class MyTypeTag>
struct LinearSolverSplice { using type = UndefinedProperty; };
template<class TypeTag, class MyTypeTag>
struct LocalLinearizerSplice { using type = UndefinedProperty; };
template<class TypeTag>
struct Splices<TypeTag, TTag::FvBaseDiscretization>
{
    using type = std::tuple<GetSplicePropType<TypeTag, TTag::FvBaseDiscretization, Properties::LinearSolverSplice>,
                            GetSplicePropType<TypeTag, TTag::FvBaseDiscretization, Properties::LocalLinearizerSplice>>;
};

//! use a parallel BiCGStab linear solver by default
template<class TypeTag>
struct LinearSolverSplice<TypeTag, TTag::FvBaseDiscretization> { using type = TTag::ParallelBiCGStabLinearSolver; };

//! by default, use finite differences to linearize the system of PDEs
template<class TypeTag>
struct LocalLinearizerSplice<TypeTag, TTag::FvBaseDiscretization> { using type = TTag::FiniteDifferenceLocalLinearizer; };

/*!
 * \brief Representation of a function evaluation and all necessary derivatives with
 *        regard to the intensive quantities of the primary variables.
 *
 * Depending on the chosen linearization method, this property may be the same as the
 * "Scalar" property (if the finite difference linearizer is used), or it may be more
 * complex (for the linearizer which uses automatic differentiation).
 */
template<class TypeTag, class MyTypeTag>
struct Evaluation { using type = UndefinedProperty; };

//! The class describing the stencil of the spatial discretization
template<class TypeTag, class MyTypeTag>
struct Stencil { using type = UndefinedProperty; };

//! The class describing the discrete function space when dune-fem is used, otherwise it points to the stencil class
template<class TypeTag, class MyTypeTag>
struct DiscreteFunctionSpace { using type = UndefinedProperty; };

//! The type of the problem
template<class TypeTag, class MyTypeTag>
struct Problem { using type = UndefinedProperty; };
//! The type of the base class for all problems which use this model
template<class TypeTag, class MyTypeTag>
struct BaseProblem { using type = UndefinedProperty; };

//! The type of the spatial discretization used by the model
template<class TypeTag, class MyTypeTag>
struct Discretization { using type = UndefinedProperty; };
//! The discretization specific part of the local residual
template<class TypeTag, class MyTypeTag>
struct DiscLocalResidual { using type = UndefinedProperty; };
//! The type of the local residual function
template<class TypeTag, class MyTypeTag>
struct LocalResidual { using type = UndefinedProperty; };
//! The type of the local linearizer
template<class TypeTag, class MyTypeTag>
struct LocalLinearizer { using type = UndefinedProperty; };
//! Specify if elements that do not belong to the local process' grid partition should be
//! skipped
template<class TypeTag, class MyTypeTag>
struct LinearizeNonLocalElements { using type = UndefinedProperty; };

//! Linearizes the global non-linear system of equations
template<class TypeTag, class MyTypeTag>
struct BaseLinearizer { using type = UndefinedProperty; };

//! A vector of holding a quantity for each equation (usually at a given spatial location)
template<class TypeTag, class MyTypeTag>
struct EqVector { using type = UndefinedProperty; };
//! A vector of holding a quantity for each equation for each DOF of an element
template<class TypeTag, class MyTypeTag>
struct ElementEqVector { using type = UndefinedProperty; };

//! Vector containing volumetric or areal rates of quantities
template<class TypeTag, class MyTypeTag>
struct RateVector { using type = UndefinedProperty; };
//! Type of object for specifying boundary conditions
template<class TypeTag, class MyTypeTag>
struct BoundaryRateVector { using type = UndefinedProperty; };
//! The class which represents a constraint degree of freedom
template<class TypeTag, class MyTypeTag>
struct Constraints { using type = UndefinedProperty; };

//! Vector containing all primary variables of the grid
template<class TypeTag, class MyTypeTag>
struct SolutionVector { using type = UndefinedProperty; };

//! A vector of primary variables within a sub-control volume
template<class TypeTag, class MyTypeTag>
struct PrimaryVariables { using type = UndefinedProperty; };
//! The secondary variables within a sub-control volume
template<class TypeTag, class MyTypeTag>
struct IntensiveQuantities { using type = UndefinedProperty; };
//! The discretization specific part of the intensive quantities
template<class TypeTag, class MyTypeTag>
struct DiscIntensiveQuantities { using type = UndefinedProperty; };

//! The secondary variables of all degrees of freedom in an element's stencil
template<class TypeTag, class MyTypeTag>
struct ElementContext { using type = UndefinedProperty; };
//! The secondary variables of a boundary segment
template<class TypeTag, class MyTypeTag>
struct BoundaryContext { using type = UndefinedProperty; };
//! The secondary variables of a constraint degree of freedom
template<class TypeTag, class MyTypeTag>
struct ConstraintsContext { using type = UndefinedProperty; };
//! Data required to calculate a flux over a face
template<class TypeTag, class MyTypeTag>
struct ExtensiveQuantities { using type = UndefinedProperty; };
//! Calculates gradients of arbitrary quantities at flux integration points
template<class TypeTag, class MyTypeTag>
struct GradientCalculator { using type = UndefinedProperty; };

//! The part of the intensive quantities which is specific to the spatial discretization
template<class TypeTag, class MyTypeTag>
struct DiscBaseIntensiveQuantities { using type = UndefinedProperty; };

//! The part of the extensive quantities which is specific to the spatial discretization
template<class TypeTag, class MyTypeTag>
struct DiscExtensiveQuantities { using type = UndefinedProperty; };

//! The part of the VTK ouput modules which is specific to the spatial discretization
template<class TypeTag, class MyTypeTag>
struct DiscBaseOutputModule { using type = UndefinedProperty; };

//! The class to create grid communication handles
template<class TypeTag, class MyTypeTag>
struct GridCommHandleFactory { using type = UndefinedProperty; };

/*!
 * \brief The OpenMP threads manager
 */
template<class TypeTag, class MyTypeTag>
struct ThreadManager { using type = UndefinedProperty; };
template<class TypeTag, class MyTypeTag>
struct ThreadsPerProcess { using type = UndefinedProperty; };

//! use locking to prevent race conditions when linearizing the global system of
//! equations in multi-threaded mode. (setting this property to true is always save, but
//! it may slightly deter performance in multi-threaded simlations and some
//! discretizations do not need this.)
template<class TypeTag, class MyTypeTag>
struct UseLinearizationLock { using type = UndefinedProperty; };

// high-level simulation control

/*!
 * \brief Switch to enable or disable grid adaptation
 *
 * Currently grid adaptation requires the presence of the dune-FEM module. If it is not
 * available and grid adaptation is enabled, an exception is thrown.
 */
template<class TypeTag, class MyTypeTag>
struct EnableGridAdaptation { using type = UndefinedProperty; };

/*!
 * \brief The directory to which simulation output ought to be written to.
 */
template<class TypeTag, class MyTypeTag>
struct OutputDir { using type = UndefinedProperty; };

/*!
 * \brief Global switch to enable or disable the writing of VTK output files
 *
 * If writing VTK files is disabled, then the WriteVtk$FOO options do
 * not have any effect...
 */
template<class TypeTag, class MyTypeTag>
struct EnableVtkOutput { using type = UndefinedProperty; };

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
struct EnableAsyncVtkOutput { using type = UndefinedProperty; };

/*!
 * \brief Specify the format the VTK output is written to disk
 *
 * Possible values are:
 *   - Dune::VTK::ascii (default)
 *   - Dune::VTK::base64
 *   - Dune::VTK::appendedraw
 *   - Dune::VTK::appendedbase64
 */
template<class TypeTag, class MyTypeTag>
struct VtkOutputFormat { using type = UndefinedProperty; };

//! Specify whether the some degrees of fredom can be constraint
template<class TypeTag, class MyTypeTag>
struct EnableConstraints { using type = UndefinedProperty; };

/*!
 * \brief Specify the maximum size of a time integration [s].
 *
 * The default is to not limit the step size.
 */
template<class TypeTag, class MyTypeTag>
struct MaxTimeStepSize { using type = UndefinedProperty; };

/*!
 * \brief Specify the minimal size of a time integration [s].
 *
 * The default is to not limit the step size.
 */
template<class TypeTag, class MyTypeTag>
struct MinTimeStepSize { using type = UndefinedProperty; };

/*!
 * \brief The maximum allowed number of timestep divisions for the
 *        Newton solver.
 */
template<class TypeTag, class MyTypeTag>
struct MaxTimeStepDivisions { using type = UndefinedProperty; };

/*!
 * \brief Continue with a non-converged solution instead of giving up
 *        if we encounter a time step size smaller than the minimum time
 *        step size.
 */
template<class TypeTag, class MyTypeTag>
struct ContinueOnConvergenceError { using type = UndefinedProperty; };

/*!
 * \brief Specify whether all intensive quantities for the grid should be
 *        cached in the discretization.
 *
 * This potentially reduces the CPU time, but comes at the cost of
 * higher memory consumption. In turn, the higher memory requirements
 * may cause the simulation to exhibit worse cache coherence behavior
 * which eats some of the computational benefits again.
 */
template<class TypeTag, class MyTypeTag>
struct EnableIntensiveQuantityCache { using type = UndefinedProperty; };

/*!
 * \brief Specify whether the storage terms for previous solutions should be cached.
 *
 * This potentially reduces the CPU time, but comes at the cost of higher memory
 * consumption.
 */
template<class TypeTag, class MyTypeTag>
struct EnableStorageCache { using type = UndefinedProperty; };

/*!
 * \brief Specify whether to use the already calculated solutions as
 *        starting values of the intensive quantities.
 *
 * This only makes sense if the calculation of the intensive quantities is
 * very expensive (e.g. for non-linear fugacity functions where the
 * solver converges faster).
 */
template<class TypeTag, class MyTypeTag>
struct EnableThermodynamicHints { using type = UndefinedProperty; };

// mappers from local to global DOF indices

/*!
 * \brief The mapper to find the global index of a vertex.
 */
template<class TypeTag, class MyTypeTag>
struct VertexMapper { using type = UndefinedProperty; };

/*!
 * \brief The mapper to find the global index of an element.
 */
template<class TypeTag, class MyTypeTag>
struct ElementMapper { using type = UndefinedProperty; };

/*!
 * \brief The mapper to find the global index of a degree of freedom.
 */
template<class TypeTag, class MyTypeTag>
struct DofMapper { using type = UndefinedProperty; };

/*!
 * \brief The history size required by the time discretization
 */
template<class TypeTag, class MyTypeTag>
struct TimeDiscHistorySize { using type = UndefinedProperty; };

/*!
 * \brief Specify whether the storage terms use extensive quantities or not.
 *
 * Most models don't need this, but the (Navier-)Stokes ones do...
 */
template<class TypeTag, class MyTypeTag>
struct ExtensiveStorageTerm { using type = UndefinedProperty; };

//! \brief Specify whether to use volumetric residuals or not
template<class TypeTag, class MyTypeTag>
struct UseVolumetricResidual { using type = UndefinedProperty; };


//! Specify if experimental features should be enabled or not.
template<class TypeTag, class MyTypeTag>
struct EnableExperiments { using type = UndefinedProperty; };

template<class TypeTag>
struct Vanguard<TypeTag, TTag::NumericModel> { using type = Opm::DgfVanguard<TypeTag>; };

} // namespace Opm::Properties

#endif
