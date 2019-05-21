/*
  Copyright 2013, 2015 SINTEF ICT, Applied Mathematics.
  Copyright 2014, 2015 Dr. Blatt - HPC-Simulation-Software & Services
  Copyright 2014, 2015 Statoil ASA.
  Copyright 2015 NTNU
  Copyright 2015, 2016, 2017 IRIS AS

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

#ifndef OPM_BLACKOILMODELEBOS_HEADER_INCLUDED
#define OPM_BLACKOILMODELEBOS_HEADER_INCLUDED

#include <ebos/eclproblem.hh>
#include <ewoms/common/start.hh>

#include <opm/simulators/timestepping/AdaptiveTimeSteppingEbos.hpp>

#include <opm/autodiff/NonlinearSolverEbos.hpp>
#include <opm/autodiff/BlackoilModelParametersEbos.hpp>
#include <opm/simulators/wells/BlackoilWellModel.hpp>
#include <opm/simulators/aquifers/BlackoilAquiferModel.hpp>
#include <opm/simulators/wells/WellConnectionAuxiliaryModule.hpp>
#include <opm/autodiff/BlackoilDetails.hpp>

#include <opm/grid/UnstructuredGrid.h>
#include <opm/core/simulator/SimulatorReport.hpp>
#include <opm/simulators/linalg/ParallelIstlInformation.hpp>
#include <opm/core/props/phaseUsageFromDeck.hpp>
#include <opm/common/ErrorMacros.hpp>
#include <opm/common/Exceptions.hpp>
#include <opm/common/OpmLog/OpmLog.hpp>
#include <opm/parser/eclipse/Units/Units.hpp>
#include <opm/simulators/timestepping/SimulatorTimer.hpp>
#include <opm/common/utility/parameters/ParameterGroup.hpp>
#include <opm/parser/eclipse/EclipseState/EclipseState.hpp>
#include <opm/parser/eclipse/EclipseState/Tables/TableManager.hpp>

#include <opm/simulators/linalg/ISTLSolverEbos.hpp>
#include <opm/common/data/SimulationDataContainer.hpp>

#include <dune/istl/owneroverlapcopy.hh>
#include <dune/common/parallel/collectivecommunication.hh>
#include <dune/common/timer.hh>
#include <dune/common/unused.hh>

#include <cassert>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <limits>
#include <vector>
#include <algorithm>
//#include <fstream>


BEGIN_PROPERTIES

NEW_TYPE_TAG(EclFlowProblem, INHERITS_FROM(BlackOilModel, EclBaseProblem, FlowNonLinearSolver, FlowModelParameters, FlowTimeSteppingParameters));
SET_STRING_PROP(EclFlowProblem, OutputDir, "");
SET_BOOL_PROP(EclFlowProblem, EnableDebuggingChecks, false);
// default in flow is to formulate the equations in surface volumes
SET_BOOL_PROP(EclFlowProblem, BlackoilConserveSurfaceVolume, true);
SET_BOOL_PROP(EclFlowProblem, UseVolumetricResidual, false);

SET_TYPE_PROP(EclFlowProblem, EclAquiferModel, Opm::BlackoilAquiferModel<TypeTag>);

// disable all extensions supported by black oil model. this should not really be
// necessary but it makes things a bit more explicit
SET_BOOL_PROP(EclFlowProblem, EnablePolymer, false);
SET_BOOL_PROP(EclFlowProblem, EnableSolvent, false);
SET_BOOL_PROP(EclFlowProblem, EnableTemperature, true);
SET_BOOL_PROP(EclFlowProblem, EnableEnergy, false);

SET_TYPE_PROP(EclFlowProblem, EclWellModel, Opm::BlackoilWellModel<TypeTag>);
SET_TAG_PROP(EclFlowProblem, LinearSolverSplice, FlowIstlSolver);


END_PROPERTIES

namespace Opm {
    /// A model implementation for three-phase black oil.
    ///
    /// The simulator is capable of handling three-phase problems
    /// where gas can be dissolved in oil and vice versa. It
    /// uses an industry-standard TPFA discretization with per-phase
    /// upwind weighting of mobilities.
    template <class TypeTag>
    class BlackoilModelEbos
    {
    public:
        // ---------  Types and enums  ---------
        typedef WellStateFullyImplicitBlackoil WellState;
        typedef BlackoilModelParametersEbos<TypeTag> ModelParameters;

        typedef typename GET_PROP_TYPE(TypeTag, Simulator)         Simulator;
        typedef typename GET_PROP_TYPE(TypeTag, Grid)              Grid;
        typedef typename GET_PROP_TYPE(TypeTag, ElementContext)    ElementContext;
        typedef typename GET_PROP_TYPE(TypeTag, SparseMatrixAdapter) SparseMatrixAdapter;
        typedef typename GET_PROP_TYPE(TypeTag, SolutionVector)    SolutionVector ;
        typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables)  PrimaryVariables ;
        typedef typename GET_PROP_TYPE(TypeTag, FluidSystem)       FluidSystem;
        typedef typename GET_PROP_TYPE(TypeTag, Indices)           Indices;
        typedef typename GET_PROP_TYPE(TypeTag, MaterialLaw)       MaterialLaw;
        typedef typename GET_PROP_TYPE(TypeTag, MaterialLawParams) MaterialLawParams;

        typedef double Scalar;
        static const int numEq = Indices::numEq;
        static const int contiSolventEqIdx = Indices::contiSolventEqIdx;
        static const int contiPolymerEqIdx = Indices::contiPolymerEqIdx;
        static const int contiEnergyEqIdx = Indices::contiEnergyEqIdx;
        static const int contiPolymerMWEqIdx = Indices::contiPolymerMWEqIdx;
        static const int solventSaturationIdx = Indices::solventSaturationIdx;
        static const int polymerConcentrationIdx = Indices::polymerConcentrationIdx;
        static const int polymerMoleWeightIdx = Indices::polymerMoleWeightIdx;
        static const int temperatureIdx = Indices::temperatureIdx;

        typedef Dune::FieldVector<Scalar, numEq >        VectorBlockType;
        typedef typename SparseMatrixAdapter::MatrixBlock MatrixBlockType;
        typedef typename SparseMatrixAdapter::IstlMatrix Mat;
        typedef Dune::BlockVector<VectorBlockType>      BVector;

        typedef ISTLSolverEbos<TypeTag> ISTLSolverType;
        //typedef typename SolutionVector :: value_type            PrimaryVariables ;

        // ---------  Public methods  ---------

        /// Construct the model. It will retain references to the
        /// arguments of this functions, and they are expected to
        /// remain in scope for the lifetime of the solver.
        /// \param[in] param            parameters
        /// \param[in] grid             grid data structure
        /// \param[in] wells            well structure
        /// \param[in] vfp_properties   Vertical flow performance tables
        /// \param[in] linsolver        linear solver
        /// \param[in] eclState         eclipse state
        /// \param[in] terminal_output  request output to cout/cerr
        BlackoilModelEbos(Simulator& ebosSimulator,
                          const ModelParameters& param,
                          BlackoilWellModel<TypeTag>& well_model,
                          const bool terminal_output)
        : ebosSimulator_(ebosSimulator)
        , grid_(ebosSimulator_.vanguard().grid())
        , phaseUsage_(phaseUsageFromDeck(eclState()))
        , has_disgas_(FluidSystem::enableDissolvedGas())
        , has_vapoil_(FluidSystem::enableVaporizedOil())
        , has_solvent_(GET_PROP_VALUE(TypeTag, EnableSolvent))
        , has_polymer_(GET_PROP_VALUE(TypeTag, EnablePolymer))
        , has_polymermw_(GET_PROP_VALUE(TypeTag, EnablePolymerMW))
        , has_energy_(GET_PROP_VALUE(TypeTag, EnableEnergy))
        , param_( param )
        , well_model_ (well_model)
        , terminal_output_ (terminal_output)
        , current_relaxation_(1.0)
        , dx_old_(UgGridHelpers::numCells(grid_))
        {
            // compute global sum of number of cells
            global_nc_ = detail::countGlobalCells(grid_);
            convergence_reports_.reserve(300); // Often insufficient, but avoids frequent moves.
        }

        bool isParallel() const
        { return  grid_.comm().size() > 1; }

        const EclipseState& eclState() const
        { return ebosSimulator_.vanguard().eclState(); }

        /// Called once before each time step.
        /// \param[in] timer                  simulation timer
        void prepareStep(const SimulatorTimerInterface& timer)
        {
            // update the solution variables in ebos
            if ( timer.lastStepFailed() ) {
                ebosSimulator_.model().updateFailed();
            } else {
                ebosSimulator_.model().advanceTimeLevel();
            }

            // set the timestep size and episode index for ebos explicitly. ebos needs to
            // know the report step/episode index because of timing dependend data
            // despide the fact that flow uses its own time stepper. (The length of the
            // episode does not matter, though.)
            ebosSimulator_.setTime(timer.simulationTimeElapsed());
            ebosSimulator_.setTimeStepSize(timer.currentStepLength());
            ebosSimulator_.problem().beginTimeStep();

            unsigned numDof = ebosSimulator_.model().numGridDof();
            wasSwitched_.resize(numDof);
            std::fill(wasSwitched_.begin(), wasSwitched_.end(), false);

            if (param_.update_equations_scaling_) {
                std::cout << "equation scaling not suported yet" << std::endl;
                //updateEquationsScaling();
            }
        }


        /// Called once per nonlinear iteration.
        /// This model will perform a Newton-Raphson update, changing reservoir_state
        /// and well_state. It will also use the nonlinear_solver to do relaxation of
        /// updates if necessary.
        /// \param[in] iteration              should be 0 for the first call of a new timestep
        /// \param[in] timer                  simulation timer
        /// \param[in] nonlinear_solver       nonlinear solver used (for oscillation/relaxation control)
        /// \param[in, out] reservoir_state   reservoir state variables
        /// \param[in, out] well_state        well state variables
        template <class NonlinearSolverType>
        SimulatorReport nonlinearIteration(const int iteration,
                                           const SimulatorTimerInterface& timer,
                                           NonlinearSolverType& nonlinear_solver)
        {
            SimulatorReport report;
            failureReport_ = SimulatorReport();
            Dune::Timer perfTimer;

            perfTimer.start();
            if (iteration == 0) {
                // For each iteration we store in a vector the norms of the residual of
                // the mass balance for each active phase, the well flux and the well equations.
                residual_norms_history_.clear();
                current_relaxation_ = 1.0;
                dx_old_ = 0.0;
                convergence_reports_.push_back({timer.reportStepNum(), timer.currentStepNum(), {}});
                convergence_reports_.back().report.reserve(11);
            }

            report.total_linearizations = 1;

            try {
                report += assembleReservoir(timer, iteration);
                report.assemble_time += perfTimer.stop();
            }
            catch (...) {
                report.assemble_time += perfTimer.stop();
                failureReport_ += report;
                // todo (?): make the report an attribute of the class
                throw; // continue throwing the stick
            }

            std::vector<double> residual_norms;
            perfTimer.reset();
            perfTimer.start();
            // the step is not considered converged until at least minIter iterations is done
            {
                auto convrep = getConvergence(timer, iteration,residual_norms);
                report.converged = convrep.converged()  && iteration > nonlinear_solver.minIter();;
                ConvergenceReport::Severity severity = convrep.severityOfWorstFailure();
                convergence_reports_.back().report.push_back(std::move(convrep));

                // Throw if any NaN or too large residual found.
                if (severity == ConvergenceReport::Severity::NotANumber) {
                    OPM_THROW(Opm::NumericalIssue, "NaN residual found!");
                } else if (severity == ConvergenceReport::Severity::TooLarge) {
                    OPM_THROW(Opm::NumericalIssue, "Too large residual found!");
                }
            }

             // checking whether the group targets are converged
             if (wellModel().wellCollection().groupControlActive()) {
                  report.converged = report.converged && wellModel().wellCollection().groupTargetConverged(wellModel().wellState().wellRates());
             }

            report.update_time += perfTimer.stop();
            residual_norms_history_.push_back(residual_norms);
            if (!report.converged) {
                perfTimer.reset();
                perfTimer.start();
                report.total_newton_iterations = 1;

                // enable single precision for solvers when dt is smaller then 20 days
                //residual_.singlePrecision = (unit::convert::to(dt, unit::day) < 20.) ;

                // Compute the nonlinear update.
                const int nc = UgGridHelpers::numCells(grid_);
                BVector x(nc);

                // apply the Schur compliment of the well model to the reservoir linearized
                // equations
                wellModel().linearize(ebosSimulator().model().linearizer().jacobian(),
                                      ebosSimulator().model().linearizer().residual());

                // Solve the linear system.
                linear_solve_setup_time_ = 0.0;
                try {
                    solveJacobianSystem(x);
                    report.linear_solve_setup_time += linear_solve_setup_time_;
                    report.linear_solve_time += perfTimer.stop();
                    report.total_linear_iterations += linearIterationsLastSolve();
                }
                catch (...) {
                    report.linear_solve_setup_time += linear_solve_setup_time_;
                    report.linear_solve_time += perfTimer.stop();
                    report.total_linear_iterations += linearIterationsLastSolve();

                    failureReport_ += report;
                    throw; // re-throw up
                }

                perfTimer.reset();
                perfTimer.start();

                // handling well state update before oscillation treatment is a decision based
                // on observation to avoid some big performance degeneration under some circumstances.
                // there is no theorectical explanation which way is better for sure.
                wellModel().postSolve(x);

                if (param_.use_update_stabilization_) {
                    // Stabilize the nonlinear update.
                    bool isOscillate = false;
                    bool isStagnate = false;
                    nonlinear_solver.detectOscillations(residual_norms_history_, iteration, isOscillate, isStagnate);
                    if (isOscillate) {
                        current_relaxation_ -= nonlinear_solver.relaxIncrement();
                        current_relaxation_ = std::max(current_relaxation_, nonlinear_solver.relaxMax());
                        if (terminalOutputEnabled()) {
                            std::string msg = "    Oscillating behavior detected: Relaxation set to "
                                    + std::to_string(current_relaxation_);
                            OpmLog::info(msg);
                        }
                    }
                    nonlinear_solver.stabilizeNonlinearUpdate(x, dx_old_, current_relaxation_);
                }

                // Apply the update, with considering model-dependent limitations and
                // chopping of the update.
                updateSolution(x);

                report.update_time += perfTimer.stop();
            }

            return report;
        }

        void printIf(int c, double x, double y, double eps, std::string type) {
            if (std::abs(x-y) > eps) {
                std::cout << type << " " <<c << ": "<<x << " " << y << std::endl;
            }
        }


        /// Called once after each time step.
        /// In this class, this function does nothing.
        /// \param[in] timer                  simulation timer
        void afterStep(const SimulatorTimerInterface& timer OPM_UNUSED)
        {
            ebosSimulator_.problem().endTimeStep();
        }

        /// Assemble the residual and Jacobian of the nonlinear system.
        /// \param[in]      reservoir_state   reservoir state variables
        /// \param[in, out] well_state        well state variables
        /// \param[in]      initial_assembly  pass true if this is the first call to assemble() in this timestep
        SimulatorReport assembleReservoir(const SimulatorTimerInterface& /* timer */,
                                          const int iterationIdx)
        {
            // -------- Mass balance equations --------
            ebosSimulator_.model().newtonMethod().setIterationIndex(iterationIdx);
            ebosSimulator_.problem().beginIteration();
            ebosSimulator_.model().linearizer().linearizeDomain();
            ebosSimulator_.problem().endIteration();

            return wellModel().lastReport();
        }

        // compute the "relative" change of the solution between time steps
        double relativeChange() const
        {
            Scalar resultDelta = 0.0;
            Scalar resultDenom = 0.0;

            const auto& elemMapper = ebosSimulator_.model().elementMapper();
            const auto& gridView = ebosSimulator_.gridView();
            auto elemIt = gridView.template begin</*codim=*/0>();
            const auto& elemEndIt = gridView.template end</*codim=*/0>();
            for (; elemIt != elemEndIt; ++elemIt) {
                const auto& elem = *elemIt;
                if (elem.partitionType() != Dune::InteriorEntity)
                    continue;

                unsigned globalElemIdx = elemMapper.index(elem);
                const auto& priVarsNew = ebosSimulator_.model().solution(/*timeIdx=*/0)[globalElemIdx];

                Scalar pressureNew;
                pressureNew = priVarsNew[Indices::pressureSwitchIdx];

                Scalar saturationsNew[FluidSystem::numPhases] = { 0.0 };
                Scalar oilSaturationNew = 1.0;
                if (FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx)) {
                    saturationsNew[FluidSystem::waterPhaseIdx] = priVarsNew[Indices::waterSaturationIdx];
                    oilSaturationNew -= saturationsNew[FluidSystem::waterPhaseIdx];
                }

                if (FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx) && priVarsNew.primaryVarsMeaning() == PrimaryVariables::Sw_po_Sg) {
                    saturationsNew[FluidSystem::gasPhaseIdx] = priVarsNew[Indices::compositionSwitchIdx];
                    oilSaturationNew -= saturationsNew[FluidSystem::gasPhaseIdx];
                }

                if (FluidSystem::phaseIsActive(FluidSystem::oilPhaseIdx)) {
                    saturationsNew[FluidSystem::oilPhaseIdx] = oilSaturationNew;
                }

                const auto& priVarsOld = ebosSimulator_.model().solution(/*timeIdx=*/1)[globalElemIdx];

                Scalar pressureOld;
                pressureOld = priVarsOld[Indices::pressureSwitchIdx];

                Scalar saturationsOld[FluidSystem::numPhases] = { 0.0 };
                Scalar oilSaturationOld = 1.0;
                if (FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx)) {
                    saturationsOld[FluidSystem::waterPhaseIdx] = priVarsOld[Indices::waterSaturationIdx];
                    oilSaturationOld -= saturationsOld[FluidSystem::waterPhaseIdx];
                }

                if (FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx) && priVarsOld.primaryVarsMeaning() == PrimaryVariables::Sw_po_Sg) {
                    saturationsOld[FluidSystem::gasPhaseIdx] = priVarsOld[Indices::compositionSwitchIdx];
                    oilSaturationOld -= saturationsOld[FluidSystem::gasPhaseIdx];
                }

                if (FluidSystem::phaseIsActive(FluidSystem::oilPhaseIdx)) {
                    saturationsOld[FluidSystem::oilPhaseIdx] = oilSaturationOld;
                }

                Scalar tmp = pressureNew - pressureOld;
                resultDelta += tmp*tmp;
                resultDenom += pressureNew*pressureNew;

                for (unsigned phaseIdx = 0; phaseIdx < FluidSystem::numPhases; ++ phaseIdx) {
                    Scalar tmp = saturationsNew[phaseIdx] - saturationsOld[phaseIdx];
                    resultDelta += tmp*tmp;
                    resultDenom += saturationsNew[phaseIdx]*saturationsNew[phaseIdx];
                }
            }

            resultDelta = gridView.comm().sum(resultDelta);
            resultDenom = gridView.comm().sum(resultDenom);

            if (resultDenom > 0.0)
                return resultDelta/resultDenom;
            return 0.0;
        }


        /// Number of linear iterations used in last call to solveJacobianSystem().
        int linearIterationsLastSolve() const
        {
            return ebosSimulator_.model().newtonMethod().linearSolver().iterations ();
        }

        /// Solve the Jacobian system Jx = r where J is the Jacobian and
        /// r is the residual.
        void solveJacobianSystem(BVector& x)
        {

            auto& ebosJac = ebosSimulator_.model().linearizer().jacobian();
            auto& ebosResid = ebosSimulator_.model().linearizer().residual();

            // set initial guess
            x = 0.0;

            auto& ebosSolver = ebosSimulator_.model().newtonMethod().linearSolver();
            Dune::Timer perfTimer;
            perfTimer.start();
            ebosSolver.prepare(ebosJac, ebosResid);
            linear_solve_setup_time_ = perfTimer.stop();
            ebosSolver.setResidual(ebosResid);
            // actually, the error needs to be calculated after setResidual in order to
            // account for parallelization properly. since the residual of ECFV
            // discretizations does not need to be synchronized across processes to be
            // consistent, this is not relevant for OPM-flow...
            ebosSolver.setMatrix(ebosJac);
            ebosSolver.solve(x);
       }



        /// Apply an update to the primary variables.
        void updateSolution(const BVector& dx)
        {
            auto& ebosNewtonMethod = ebosSimulator_.model().newtonMethod();
            SolutionVector& solution = ebosSimulator_.model().solution(/*timeIdx=*/0);

            ebosNewtonMethod.update_(/*nextSolution=*/solution,
                                     /*curSolution=*/solution,
                                     /*update=*/dx,
                                     /*resid=*/dx); // the update routines of the black
                                                    // oil model do not care about the
                                                    // residual

            // if the solution is updated, the intensive quantities need to be recalculated
            ebosSimulator_.model().invalidateIntensiveQuantitiesCache(/*timeIdx=*/0);
        }

        /// Return true if output to cout is wanted.
        bool terminalOutputEnabled() const
        {
            return terminal_output_;
        }

        template <class CollectiveCommunication>
        double convergenceReduction(const CollectiveCommunication& comm,
                                    const double pvSumLocal,
                                    std::vector< Scalar >& R_sum,
                                    std::vector< Scalar >& maxCoeff,
                                    std::vector< Scalar >& B_avg)
        {
            // Compute total pore volume (use only owned entries)
            double pvSum = pvSumLocal;

            if( comm.size() > 1 )
            {
                // global reduction
                std::vector< Scalar > sumBuffer;
                std::vector< Scalar > maxBuffer;
                const int numComp = B_avg.size();
                sumBuffer.reserve( 2*numComp + 1 ); // +1 for pvSum
                maxBuffer.reserve( numComp );
                for( int compIdx = 0; compIdx < numComp; ++compIdx )
                {
                    sumBuffer.push_back( B_avg[ compIdx ] );
                    sumBuffer.push_back( R_sum[ compIdx ] );
                    maxBuffer.push_back( maxCoeff[ compIdx ] );
                }

                // Compute total pore volume
                sumBuffer.push_back( pvSum );

                // compute global sum
                comm.sum( sumBuffer.data(), sumBuffer.size() );

                // compute global max
                comm.max( maxBuffer.data(), maxBuffer.size() );

                // restore values to local variables
                for( int compIdx = 0, buffIdx = 0; compIdx < numComp; ++compIdx, ++buffIdx )
                {
                    B_avg[ compIdx ]    = sumBuffer[ buffIdx ];
                    ++buffIdx;

                    R_sum[ compIdx ]       = sumBuffer[ buffIdx ];
                }

                for( int compIdx = 0; compIdx < numComp; ++compIdx )
                {
                    maxCoeff[ compIdx ] = maxBuffer[ compIdx ];
                }

                // restore global pore volume
                pvSum = sumBuffer.back();
            }

            // return global pore volume
            return pvSum;
        }

        // Get reservoir quantities on this process needed for convergence calculations.
        double localConvergenceData(std::vector<Scalar>& R_sum,
                                    std::vector<Scalar>& maxCoeff,
                                    std::vector<Scalar>& B_avg)
        {
            double pvSumLocal = 0.0;
            const auto& ebosModel = ebosSimulator_.model();
            const auto& ebosProblem = ebosSimulator_.problem();

            const auto& ebosResid = ebosSimulator_.model().linearizer().residual();

            ElementContext elemCtx(ebosSimulator_);
            const auto& gridView = ebosSimulator().gridView();
            const auto& elemEndIt = gridView.template end</*codim=*/0, Dune::Interior_Partition>();

            for (auto elemIt = gridView.template begin</*codim=*/0, Dune::Interior_Partition>();
                 elemIt != elemEndIt;
                 ++elemIt)
            {
                const auto& elem = *elemIt;
                elemCtx.updatePrimaryStencil(elem);
                elemCtx.updatePrimaryIntensiveQuantities(/*timeIdx=*/0);
                const unsigned cell_idx = elemCtx.globalSpaceIndex(/*spaceIdx=*/0, /*timeIdx=*/0);
                const auto& intQuants = elemCtx.intensiveQuantities(/*spaceIdx=*/0, /*timeIdx=*/0);
                const auto& fs = intQuants.fluidState();

                const double pvValue = ebosProblem.referencePorosity(cell_idx, /*timeIdx=*/0) * ebosModel.dofTotalVolume( cell_idx );
                pvSumLocal += pvValue;

                for (unsigned phaseIdx = 0; phaseIdx < FluidSystem::numPhases; ++phaseIdx)
                {
                    if (!FluidSystem::phaseIsActive(phaseIdx)) {
                        continue;
                    }

                    const unsigned compIdx = Indices::canonicalToActiveComponentIndex(FluidSystem::solventComponentIndex(phaseIdx));

                    B_avg[ compIdx ] += 1.0 / fs.invB(phaseIdx).value();
                    const auto R2 = ebosResid[cell_idx][compIdx];

                    R_sum[ compIdx ] += R2;
                    maxCoeff[ compIdx ] = std::max( maxCoeff[ compIdx ], std::abs( R2 ) / pvValue );
                }

                if ( has_solvent_ ) {
                    B_avg[ contiSolventEqIdx ] += 1.0 / intQuants.solventInverseFormationVolumeFactor().value();
                    const auto R2 = ebosResid[cell_idx][contiSolventEqIdx];
                    R_sum[ contiSolventEqIdx ] += R2;
                    maxCoeff[ contiSolventEqIdx ] = std::max( maxCoeff[ contiSolventEqIdx ], std::abs( R2 ) / pvValue );
                }
                if (has_polymer_ ) {
                    B_avg[ contiPolymerEqIdx ] += 1.0 / fs.invB(FluidSystem::waterPhaseIdx).value();
                    const auto R2 = ebosResid[cell_idx][contiPolymerEqIdx];
                    R_sum[ contiPolymerEqIdx ] += R2;
                    maxCoeff[ contiPolymerEqIdx ] = std::max( maxCoeff[ contiPolymerEqIdx ], std::abs( R2 ) / pvValue );
                }

                if (has_polymermw_) {
                    assert(has_polymer_);

                    B_avg[contiPolymerMWEqIdx] += 1.0 / fs.invB(FluidSystem::waterPhaseIdx).value();
                    // the residual of the polymer molecular equation is scaled down by a 100, since molecular weight
                    // can be much bigger than 1, and this equation shares the same tolerance with other mass balance equations
                    // TODO: there should be a more general way to determine the scaling-down coefficient
                    const auto R2 = ebosResid[cell_idx][contiPolymerMWEqIdx] / 100.;
                    R_sum[contiPolymerMWEqIdx] += R2;
                    maxCoeff[contiPolymerMWEqIdx] = std::max( maxCoeff[contiPolymerMWEqIdx], std::abs( R2 ) / pvValue );
                }

                if (has_energy_ ) {
                    B_avg[ contiEnergyEqIdx ] += 1.0;
                    const auto R2 = ebosResid[cell_idx][contiEnergyEqIdx];
                    R_sum[ contiEnergyEqIdx ] += R2;
                    maxCoeff[ contiEnergyEqIdx ] = std::max( maxCoeff[ contiEnergyEqIdx ], std::abs( R2 ) / pvValue );
                }

            }

            // compute local average in terms of global number of elements
            const int bSize = B_avg.size();
            for ( int i = 0; i<bSize; ++i )
            {
                B_avg[ i ] /= Scalar( global_nc_ );
            }

            return pvSumLocal;
        }

        ConvergenceReport getReservoirConvergence(const double dt,
                                                  const int iteration,
                                                  std::vector<Scalar>& B_avg,
                                                  std::vector<Scalar>& residual_norms)
        {
            typedef std::vector< Scalar > Vector;

            const double tol_mb  = param_.tolerance_mb_;
            const double tol_cnv = (iteration < param_.max_strict_iter_) ? param_.tolerance_cnv_ : param_.tolerance_cnv_relaxed_;

            const int numComp = numEq;
            Vector R_sum(numComp, 0.0 );
            Vector maxCoeff(numComp, std::numeric_limits< Scalar >::lowest() );
            const double pvSumLocal = localConvergenceData(R_sum, maxCoeff, B_avg);

            // compute global sum and max of quantities
            const double pvSum = convergenceReduction(grid_.comm(), pvSumLocal,
                                                      R_sum, maxCoeff, B_avg);

            // Finish computation
            std::vector<Scalar> CNV(numComp);
            std::vector<Scalar> mass_balance_residual(numComp);
            for ( int compIdx = 0; compIdx < numComp; ++compIdx )
            {
                CNV[compIdx]                    = B_avg[compIdx] * dt * maxCoeff[compIdx];
                mass_balance_residual[compIdx]  = std::abs(B_avg[compIdx]*R_sum[compIdx]) * dt / pvSum;
                residual_norms.push_back(CNV[compIdx]);
            }

            // Setup component names, only the first time the function is run.
            static std::vector<std::string> compNames;
            if (compNames.empty()) {
                compNames.resize(numComp);
                for (unsigned phaseIdx = 0; phaseIdx < FluidSystem::numPhases; ++phaseIdx) {
                    if (!FluidSystem::phaseIsActive(phaseIdx)) {
                        continue;
                    }
                    const unsigned canonicalCompIdx = FluidSystem::solventComponentIndex(phaseIdx);
                    const unsigned compIdx = Indices::canonicalToActiveComponentIndex(canonicalCompIdx);
                    compNames[compIdx] = FluidSystem::componentName(canonicalCompIdx);
                }
                if (has_solvent_) {
                    compNames[solventSaturationIdx] = "Solvent";
                }
                if (has_polymer_) {
                    compNames[polymerConcentrationIdx] = "Polymer";
                }
                if (has_polymermw_) {
                    assert(has_polymer_);
                    compNames[polymerMoleWeightIdx] = "MolecularWeightP";
                }
                if (has_energy_) {
                    compNames[temperatureIdx] = "Energy";
                }
            }

            // Create convergence report.
            ConvergenceReport report;
            using CR = ConvergenceReport;
            for (int compIdx = 0; compIdx < numComp; ++compIdx) {
                double res[2] = { mass_balance_residual[compIdx], CNV[compIdx] };
                CR::ReservoirFailure::Type types[2] = { CR::ReservoirFailure::Type::MassBalance,
                                                        CR::ReservoirFailure::Type::Cnv };
                double tol[2] = { tol_mb, tol_cnv };
                for (int ii : {0, 1}) {
                    if (std::isnan(res[ii])) {
                        report.setReservoirFailed({types[ii], CR::Severity::NotANumber, compIdx});
                        if ( terminal_output_ ) {
                            OpmLog::debug("NaN residual for " + compNames[compIdx] + " equation.");
                        }
                    } else if (res[ii] > maxResidualAllowed()) {
                        report.setReservoirFailed({types[ii], CR::Severity::TooLarge, compIdx});
                        if ( terminal_output_ ) {
                            OpmLog::debug("Too large residual for " + compNames[compIdx] + " equation.");
                        }
                    } else if (res[ii] < 0.0) {
                        report.setReservoirFailed({types[ii], CR::Severity::Normal, compIdx});
                        if ( terminal_output_ ) {
                            OpmLog::debug("Negative residual for " + compNames[compIdx] + " equation.");
                        }
                    } else if (res[ii] > tol[ii]) {
                        report.setReservoirFailed({types[ii], CR::Severity::Normal, compIdx});
                    }
                }
            }

            // Output of residuals.
            if ( terminal_output_ )
            {
                // Only rank 0 does print to std::cout
                if (iteration == 0) {
                    std::string msg = "Iter";
                    for (int compIdx = 0; compIdx < numComp; ++compIdx) {
                        msg += "    MB(";
                        msg += compNames[compIdx][0];
                        msg += ")  ";
                    }
                    for (int compIdx = 0; compIdx < numComp; ++compIdx) {
                        msg += "    CNV(";
                        msg += compNames[compIdx][0];
                        msg += ") ";
                    }
                    OpmLog::debug(msg);
                }
                std::ostringstream ss;
                const std::streamsize oprec = ss.precision(3);
                const std::ios::fmtflags oflags = ss.setf(std::ios::scientific);
                ss << std::setw(4) << iteration;
                for (int compIdx = 0; compIdx < numComp; ++compIdx) {
                    ss << std::setw(11) << mass_balance_residual[compIdx];
                }
                for (int compIdx = 0; compIdx < numComp; ++compIdx) {
                    ss << std::setw(11) << CNV[compIdx];
                }
                ss.precision(oprec);
                ss.flags(oflags);
                OpmLog::debug(ss.str());
            }

            return report;
        }

        /// Compute convergence based on total mass balance (tol_mb) and maximum
        /// residual mass balance (tol_cnv).
        /// \param[in]   timer       simulation timer
        /// \param[in]   iteration   current iteration number
        /// \param[out]  residual_norms   CNV residuals by phase
        ConvergenceReport getConvergence(const SimulatorTimerInterface& timer,
                                         const int iteration,
                                         std::vector<double>& residual_norms)
        {
            // Get convergence reports for reservoir and wells.
            std::vector<Scalar> B_avg(numEq, 0.0);
            auto report = getReservoirConvergence(timer.currentStepLength(), iteration, B_avg, residual_norms);
            report += wellModel().getWellConvergence(B_avg);

            return report;
        }


        /// The number of active fluid phases in the model.
        int numPhases() const
        {
            return phaseUsage_.num_phases;
        }

        /// Wrapper required due to not following generic API
        template<class T>
        std::vector<std::vector<double> >
        computeFluidInPlace(const T&, const std::vector<int>& fipnum) const
        {
            return computeFluidInPlace(fipnum);
        }

        /// Should not be called
        std::vector<std::vector<double> >
        computeFluidInPlace(const std::vector<int>& /*fipnum*/) const
        {
            //assert(true)
            //return an empty vector
            std::vector<std::vector<double> > regionValues(0, std::vector<double>(0,0.0));
            return regionValues;
        }

        const Simulator& ebosSimulator() const
        { return ebosSimulator_; }

        Simulator& ebosSimulator()
        { return ebosSimulator_; }

        /// return the statistics if the nonlinearIteration() method failed
        const SimulatorReport& failureReport() const
        { return failureReport_; }

        struct StepReport
        {
            int report_step;
            int current_step;
            std::vector<ConvergenceReport> report;
        };

        const std::vector<StepReport>& stepReports() const
        {
            return convergence_reports_;
        }

    protected:
        const ISTLSolverType& istlSolver() const
        {
            assert( istlSolver_ );
            return *istlSolver_;
        }

        // ---------  Data members  ---------

        Simulator& ebosSimulator_;
        const Grid&            grid_;
        const ISTLSolverType*  istlSolver_;
        const PhaseUsage phaseUsage_;
        const bool has_disgas_;
        const bool has_vapoil_;
        const bool has_solvent_;
        const bool has_polymer_;
        const bool has_polymermw_;
        const bool has_energy_;

        ModelParameters                 param_;
        SimulatorReport failureReport_;

        // Well Model
        BlackoilWellModel<TypeTag>& well_model_;

        /// \brief Whether we print something to std::cout
        bool terminal_output_;
        /// \brief The number of cells of the global grid.
        long int global_nc_;

        std::vector<std::vector<double>> residual_norms_history_;
        double current_relaxation_;
        BVector dx_old_;

        std::vector<StepReport> convergence_reports_;
    public:
        /// return the StandardWells object
        BlackoilWellModel<TypeTag>&
        wellModel() { return well_model_; }

        const BlackoilWellModel<TypeTag>&
        wellModel() const { return well_model_; }

        void beginReportStep(bool isRestart)
        {
            ebosSimulator_.problem().beginEpisode(isRestart);
        }

        void endReportStep()
        {
            ebosSimulator_.problem().endEpisode();
        }

    private:

        double dpMaxRel() const { return param_.dp_max_rel_; }
        double dsMax() const { return param_.ds_max_; }
        double drMaxRel() const { return param_.dr_max_rel_; }
        double maxResidualAllowed() const { return param_.max_residual_allowed_; }
        double linear_solve_setup_time_;
    public:
        std::vector<bool> wasSwitched_;
    };
} // namespace Opm

#endif // OPM_BLACKOILMODELBASE_IMPL_HEADER_INCLUDED
