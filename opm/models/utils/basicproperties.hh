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
 *
 * \brief Defines a type tags and some fundamental properties all models.
 */
#ifndef EWOMS_BASIC_PROPERTIES_HH
#define EWOMS_BASIC_PROPERTIES_HH

#include <dune/common/parametertree.hh>

// explicitly guard the include so that the property system
// header doesn't need to be opened and checked all the time
#ifndef OPM_PROPERTY_SYSTEM_HH
#include <opm/models/utils/propertysystem.hh>

// remove this after release 3.1 to disable macros per default
#ifndef OPM_ENABLE_OLD_PROPERTY_MACROS
#define OPM_ENABLE_OLD_PROPERTY_MACROS 1
#endif

// remove this after release 3.2 to remove macros completely
#if OPM_ENABLE_OLD_PROPERTY_MACROS
#include <opm/models/utils/propertysystemmacros.hh>
#endif // OPM_ENABLE_OLD_PROPERTY_MACROS
#endif // OPM_PROPERTY_SYSTEM_HH

#include <opm/models/utils/parametersystem.hh>

#if HAVE_DUNE_FEM
#include <dune/fem/gridpart/adaptiveleafgridpart.hh>
#endif

#include <string>

BEGIN_PROPERTIES

///////////////////////////////////
// Type tag definitions:
//
// NumericModel
// |
// +-> ImplicitModel
///////////////////////////////////

NEW_TYPE_TAG(ParameterSystem);
//! Type tag for all models.
NEW_TYPE_TAG(NumericModel, INHERITS_FROM(ParameterSystem));

//! Type tag for all fully coupled models.
NEW_TYPE_TAG(ImplicitModel, INHERITS_FROM(NumericModel));

///////////////////////////////////
// Property names which are always available:
//
// Scalar
///////////////////////////////////

//! Property to specify the type of scalar values.
NEW_PROP_TAG(Scalar);

//! Property which provides a Dune::ParameterTree.
NEW_PROP_TAG(ParameterTree);

//! Property which defines the group that is queried for parameters by default
NEW_PROP_TAG(ModelParameterGroup);

//! Property which provides a Vanguard (manages grids)
NEW_PROP_TAG(Vanguard);

NEW_PROP_TAG(GridView);

NEW_PROP_TAG(Simulator);
NEW_PROP_TAG(Grid);
NEW_PROP_TAG(GridFile);
NEW_PROP_TAG(Model);
NEW_PROP_TAG(Problem);
NEW_PROP_TAG(ThreadManager);
NEW_PROP_TAG(NewtonMethod);
NEW_PROP_TAG(SolutionVector);
NEW_PROP_TAG(GlobalEqVector);
//NEW_PROP_TAG(VtkOutputFormat);

//! Specifies the type of a solution for a single degee of freedom
NEW_PROP_TAG(PrimaryVariables);

//! Specifies whether the problem to be simulated exhibits contraint degrees of freedom
NEW_PROP_TAG(EnableConstraints);

//! Specifies the type of objects which specify constraints for a single degee of freedom
NEW_PROP_TAG(Constraints);

//! Vector containing a quantity of for equation for a single degee of freedom
NEW_PROP_TAG(EqVector);

//! The class which linearizes the non-linear system of equations
NEW_PROP_TAG(Linearizer);

//! Specifies the type of a global Jacobian matrix
NEW_PROP_TAG(SparseMatrixAdapter);

//! Specifies the type of the linear solver to be used
NEW_PROP_TAG(LinearSolverBackend);

//! Specifies whether the Newton method should print messages or not
NEW_PROP_TAG(NewtonVerbose);

//! Specifies the type of the class which writes out the Newton convergence
NEW_PROP_TAG(NewtonConvergenceWriter);

//! Specifies whether the convergence rate and the global residual
//! gets written out to disk for every Newton iteration
NEW_PROP_TAG(NewtonWriteConvergence);

//! Specifies whether the convergence rate and the global residual
//! gets written out to disk for every Newton iteration
NEW_PROP_TAG(ConvergenceWriter);

/*!
 * \brief The value for the error below which convergence is declared
 *
 * This value can (and for the porous media models will) be changed to account for grid
 * scaling and other effects.
 */
NEW_PROP_TAG(NewtonTolerance);

//! The maximum error which may occur in a simulation before the
//! Newton method for the time step is aborted
NEW_PROP_TAG(NewtonMaxError);

/*!
 * \brief The number of iterations at which the Newton method
 *        should aim at.
 *
 * This is used to control the time-step size. The heuristic used
 * is to scale the last time-step size by the deviation of the
 * number of iterations used from the target steps.
 */
NEW_PROP_TAG(NewtonTargetIterations);

//! Number of maximum iterations for the Newton method.
NEW_PROP_TAG(NewtonMaxIterations);



#if HAVE_DUNE_FEM
NEW_PROP_TAG(GridPart);
#endif

NEW_PROP_TAG(LocalLinearizer);
NEW_PROP_TAG(Evaluation);
NEW_PROP_TAG(NumericDifferenceMethod);
NEW_PROP_TAG(BaseEpsilon);
NEW_PROP_TAG(LocalResidual);
NEW_PROP_TAG(ElementContext);

NEW_PROP_TAG(NumPhases);
NEW_PROP_TAG(NumComponents);
NEW_PROP_TAG(NumEq);
NEW_PROP_TAG(FluidSystem);
NEW_PROP_TAG(DiscBaseOutputModule);

// create new type tag for the VTK primary variables output
NEW_PROP_TAG(EnableVtkOutput);

// create the property tags needed for the primary variables module
NEW_PROP_TAG(VtkWritePrimaryVars);
NEW_PROP_TAG(VtkWriteProcessRank);
NEW_PROP_TAG(VtkWriteDofIndex);
NEW_PROP_TAG(VtkWriteExtrusionFactor);
NEW_PROP_TAG(VtkWritePressures);
NEW_PROP_TAG(VtkWriteDensities);
NEW_PROP_TAG(VtkWriteSaturations);
NEW_PROP_TAG(VtkWriteMobilities);
NEW_PROP_TAG(VtkWriteRelativePermeabilities);
NEW_PROP_TAG(VtkWriteViscosities);
NEW_PROP_TAG(VtkWriteAverageMolarMasses);
NEW_PROP_TAG(VtkWritePorosity);
NEW_PROP_TAG(VtkWriteIntrinsicPermeabilities);
NEW_PROP_TAG(VtkWritePotentialGradients);
NEW_PROP_TAG(VtkWriteFilterVelocities);
NEW_PROP_TAG(VtkWriteTemperature);
NEW_PROP_TAG(VtkWriteSolidInternalEnergy);
NEW_PROP_TAG(VtkWriteThermalConductivity);
NEW_PROP_TAG(VtkWriteInternalEnergies);
NEW_PROP_TAG(VtkWriteEnthalpies);
NEW_PROP_TAG(IntensiveQuantities);

NEW_PROP_TAG(BoundaryContext);
NEW_PROP_TAG(BoundaryRateVector);
NEW_PROP_TAG(CellsX);
NEW_PROP_TAG(CellsY);
NEW_PROP_TAG(CellsZ);
NEW_PROP_TAG(ContinueOnConvergenceError);
NEW_PROP_TAG(DiscExtensiveQuantities);
NEW_PROP_TAG(DiscIntensiveQuantities);
NEW_PROP_TAG(DiscLocalResidual);
NEW_PROP_TAG(Discretization);
NEW_PROP_TAG(DofMapper);
NEW_PROP_TAG(DomainSizeX);
NEW_PROP_TAG(DomainSizeY);
NEW_PROP_TAG(DomainSizeZ);
NEW_PROP_TAG(ElementMapper);
NEW_PROP_TAG(EnableAsyncVtkOutput);
NEW_PROP_TAG(EnableEnergy);
NEW_PROP_TAG(EnableGravity);
NEW_PROP_TAG(EnableGridAdaptation);
NEW_PROP_TAG(EnableStorageCache);
NEW_PROP_TAG(ExtensiveQuantities);
NEW_PROP_TAG(ExtensiveStorageTerm);
NEW_PROP_TAG(Fluid);
NEW_PROP_TAG(FluxModule);
NEW_PROP_TAG(GradientCalculator);
NEW_PROP_TAG(GridCommHandleFactory);
NEW_PROP_TAG(Indices);
NEW_PROP_TAG(LinearizeNonLocalElements);
NEW_PROP_TAG(MaterialLaw);
NEW_PROP_TAG(MaterialLawParams);
NEW_PROP_TAG(MaxTimeStepDivisions);
NEW_PROP_TAG(MaxTimeStepSize);
NEW_PROP_TAG(MinTimeStepSize);
NEW_PROP_TAG(OutputDir);
NEW_PROP_TAG(RateVector);
NEW_PROP_TAG(SolidEnergyLaw);
NEW_PROP_TAG(Stencil);
NEW_PROP_TAG(ThermalConductionLaw);
NEW_PROP_TAG(ThreadsPerProcess);
NEW_PROP_TAG(TimeDiscHistorySize);
NEW_PROP_TAG(UseLinearizationLock);
NEW_PROP_TAG(UseP1FiniteElementGradients);
NEW_PROP_TAG(UseVolumetricResidual);
NEW_PROP_TAG(VertexMapper);
NEW_PROP_TAG(SolidEnergyLawParams);
NEW_PROP_TAG(ThermalConductionLawParams);

NEW_PROP_TAG(BaseProblem);
NEW_PROP_TAG(ConstraintsContext);
NEW_PROP_TAG(ElementEqVector);
NEW_PROP_TAG(EnableExperiments);
NEW_PROP_TAG(EnableIntensiveQuantityCache);
NEW_PROP_TAG(EnableThermodynamicHints);
NEW_PROP_TAG(NonwettingPhase);
NEW_PROP_TAG(WettingPhase);

NEW_PROP_TAG(OverlappingMatrix);
NEW_PROP_TAG(OverlappingVector);
NEW_PROP_TAG(PreconditionerOrder);
NEW_PROP_TAG(PreconditionerRelaxation);

NEW_PROP_TAG(AmgCoarsenTarget);

NEW_PROP_TAG(BorderListCreator);
NEW_PROP_TAG(Overlap);
NEW_PROP_TAG(OverlappingScalarProduct);
NEW_PROP_TAG(OverlappingLinearOperator);


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

NEW_PROP_TAG(LinearSolverMaxError);

NEW_PROP_TAG(LinearSolverSplice);
NEW_PROP_TAG(LocalLinearizerSplice);

//! The discretization specific part of he implementing the Newton algorithm
NEW_PROP_TAG(DiscNewtonMethod);

//! Property which tells the Vanguard how often the grid should be refined
//! after creation.
NEW_PROP_TAG(GridGlobalRefinements);

//! Property provides the name of the file from which the additional runtime
//! parameters should to be loaded from
NEW_PROP_TAG(ParameterFile);

/*!
 * \brief Print all properties on startup?
 *
 * 0 means 'no', 1 means 'yes', 2 means 'print only to logfiles'. The
 * default is 2.
 */
NEW_PROP_TAG(PrintProperties);

/*!
 * \brief Print all parameters on startup?
 *
 * 0 means 'no', 1 means 'yes', 2 means 'print only to logfiles'. The
 * default is 2.
 */
NEW_PROP_TAG(PrintParameters);

//! The default value for the simulation's end time
NEW_PROP_TAG(EndTime);

//! The default value for the simulation's initial time step size
NEW_PROP_TAG(InitialTimeStepSize);

//! The default value for the simulation's restart time
NEW_PROP_TAG(RestartTime);

//! The name of the file with a number of forced time step lengths
NEW_PROP_TAG(PredeterminedTimeStepsFile);

NEW_PROP_TAG(ParameterMetaData);

///////////////////////////////////
// Values for the properties
///////////////////////////////////

//! Set the default type of scalar values to double
SET_TYPE_PROP(NumericModel, Scalar, double);

//! Set the ParameterTree property
SET_PROP(NumericModel, ParameterTree)
{
    typedef Dune::ParameterTree type;

    static Dune::ParameterTree& tree()
    {
        static Dune::ParameterTree obj_;
        return obj_;
    }
};

//! use the global group as default for the model's parameter group
SET_STRING_PROP(NumericModel, ModelParameterGroup, "");


//! Set a value for the GridFile property
SET_STRING_PROP(NumericModel, GridFile, "");

#if HAVE_DUNE_FEM
SET_PROP(NumericModel, GridPart)
{
    typedef typename GET_PROP_TYPE(TypeTag, Grid) Grid;
    typedef Dune::Fem::AdaptiveLeafGridPart<Grid> type;
};

SET_TYPE_PROP(NumericModel, GridView, typename GET_PROP_TYPE(TypeTag, GridPart)::GridViewType);
#else
//! Use the leaf grid view by default.
//!
//! Except for spatial refinement, there is rarly a reason to use
//! anything else...
SET_TYPE_PROP(NumericModel, GridView, typename GET_PROP_TYPE(TypeTag, Grid)::LeafGridView);
#endif

//! Set a value for the ParameterFile property
SET_STRING_PROP(NumericModel, ParameterFile, "");

//! Set the number of refinement levels of the grid to 0. This does not belong
//! here, strictly speaking.
SET_INT_PROP(NumericModel, GridGlobalRefinements, 0);

//! By default, print the properties on startup
SET_INT_PROP(NumericModel, PrintProperties, 2);

//! By default, print the values of the run-time parameters on startup
SET_INT_PROP(NumericModel, PrintParameters, 2);

//! The default value for the simulation's end time
SET_SCALAR_PROP(NumericModel, EndTime, -1e35);

//! The default value for the simulation's initial time step size
SET_SCALAR_PROP(NumericModel, InitialTimeStepSize, -1e35);

//! The default value for the simulation's restart time
SET_SCALAR_PROP(NumericModel, RestartTime, -1e35);

//! By default, do not force any time steps
SET_STRING_PROP(NumericModel, PredeterminedTimeStepsFile, "");


END_PROPERTIES

#endif
