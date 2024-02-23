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

#ifndef OPM_BLACKOILMODEL_HEADER_INCLUDED
#define OPM_BLACKOILMODEL_HEADER_INCLUDED

#include <fmt/format.h>

#include <ebos/eclproblem.hh>

#include <opm/common/ErrorMacros.hpp>
#include <opm/common/Exceptions.hpp>
#include <opm/common/OpmLog/OpmLog.hpp>

#include <opm/core/props/phaseUsageFromDeck.hpp>

#include <opm/input/eclipse/EclipseState/EclipseState.hpp>
#include <opm/input/eclipse/EclipseState/Tables/TableManager.hpp>

#include <opm/simulators/aquifers/AquiferGridUtils.hpp>
#include <opm/simulators/aquifers/BlackoilAquiferModel.hpp>
#include <opm/simulators/flow/BlackoilModelNldd.hpp>
#include <opm/simulators/flow/BlackoilModelParameters.hpp>
#include <opm/simulators/flow/countGlobalCells.hpp>
#include <opm/simulators/flow/NonlinearSolver.hpp>
#include <opm/simulators/flow/RSTConv.hpp>
#include <opm/simulators/timestepping/AdaptiveTimeStepping.hpp>
#include <opm/simulators/timestepping/ConvergenceReport.hpp>
#include <opm/simulators/timestepping/SimulatorReport.hpp>
#include <opm/simulators/timestepping/SimulatorTimer.hpp>
#include <opm/simulators/utils/ComponentName.hpp>
#include <opm/simulators/utils/DeferredLoggingErrorHelpers.hpp>
#include <opm/simulators/utils/ParallelCommunication.hpp>
#include <opm/simulators/wells/BlackoilWellModel.hpp>

#include <dune/common/timer.hh>

#include <fmt/format.h>

#include <algorithm>
#include <cassert>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <ios>
#include <limits>
#include <memory>
#include <sstream>
#include <tuple>
#include <utility>
#include <vector>

namespace Opm::Properties {

namespace TTag {
struct FlowProblem {
    using InheritsFrom = std::tuple<FlowTimeSteppingParameters, FlowModelParameters,
                                    FlowNonLinearSolver, EclBaseProblem, BlackOilModel>;
};
}
template<class TypeTag>
struct OutputDir<TypeTag, TTag::FlowProblem> {
    static constexpr auto value = "";
};
template<class TypeTag>
struct EnableDebuggingChecks<TypeTag, TTag::FlowProblem> {
    static constexpr bool value = false;
};
// default in flow is to formulate the equations in surface volumes
template<class TypeTag>
struct BlackoilConserveSurfaceVolume<TypeTag, TTag::FlowProblem> {
    static constexpr bool value = true;
};
template<class TypeTag>
struct UseVolumetricResidual<TypeTag, TTag::FlowProblem> {
    static constexpr bool value = false;
};

template<class TypeTag>
struct AquiferModel<TypeTag, TTag::FlowProblem> {
    using type = BlackoilAquiferModel<TypeTag>;
};

// disable all extensions supported by black oil model. this should not really be
// necessary but it makes things a bit more explicit
template<class TypeTag>
struct EnablePolymer<TypeTag, TTag::FlowProblem> {
    static constexpr bool value = false;
};
template<class TypeTag>
struct EnableSolvent<TypeTag, TTag::FlowProblem> {
    static constexpr bool value = false;
};
template<class TypeTag>
struct EnableTemperature<TypeTag, TTag::FlowProblem> {
    static constexpr bool value = true;
};
template<class TypeTag>
struct EnableEnergy<TypeTag, TTag::FlowProblem> {
    static constexpr bool value = false;
};
template<class TypeTag>
struct EnableFoam<TypeTag, TTag::FlowProblem> {
    static constexpr bool value = false;
};
template<class TypeTag>
struct EnableBrine<TypeTag, TTag::FlowProblem> {
    static constexpr bool value = false;
};
template<class TypeTag>
struct EnableSaltPrecipitation<TypeTag, TTag::FlowProblem> {
    static constexpr bool value = false;
};
template<class TypeTag>
struct EnableMICP<TypeTag, TTag::FlowProblem> {
    static constexpr bool value = false;
};
template<class TypeTag>
struct EnableDispersion<TypeTag, TTag::FlowProblem> {
    static constexpr bool value = false;
};

template<class TypeTag>
struct WellModel<TypeTag, TTag::FlowProblem> {
    using type = BlackoilWellModel<TypeTag>;
};
template<class TypeTag>
struct LinearSolverSplice<TypeTag, TTag::FlowProblem> {
    using type = TTag::FlowIstlSolver;
};

} // namespace Opm::Properties

namespace Opm {

    /// A model implementation for three-phase black oil.
    ///
    /// The simulator is capable of handling three-phase problems
    /// where gas can be dissolved in oil and vice versa. It
    /// uses an industry-standard TPFA discretization with per-phase
    /// upwind weighting of mobilities.
    template <class TypeTag>
    class BlackoilModel
    {
    public:
        // ---------  Types and enums  ---------
        using ModelParameters = BlackoilModelParameters<TypeTag>;

        using Simulator = GetPropType<TypeTag, Properties::Simulator>;
        using Grid = GetPropType<TypeTag, Properties::Grid>;
        using ElementContext = GetPropType<TypeTag, Properties::ElementContext>;
        using IntensiveQuantities = GetPropType<TypeTag, Properties::IntensiveQuantities>;
        using SparseMatrixAdapter = GetPropType<TypeTag, Properties::SparseMatrixAdapter>;
        using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;
        using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
        using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
        using Indices = GetPropType<TypeTag, Properties::Indices>;
        using MaterialLaw = GetPropType<TypeTag, Properties::MaterialLaw>;
        using MaterialLawParams = GetPropType<TypeTag, Properties::MaterialLawParams>;
        using Scalar = GetPropType<TypeTag, Properties::Scalar>;

        static constexpr int numEq = Indices::numEq;
        static constexpr int contiSolventEqIdx = Indices::contiSolventEqIdx;
        static constexpr int contiZfracEqIdx = Indices::contiZfracEqIdx;
        static constexpr int contiPolymerEqIdx = Indices::contiPolymerEqIdx;
        static constexpr int contiEnergyEqIdx = Indices::contiEnergyEqIdx;
        static constexpr int contiPolymerMWEqIdx = Indices::contiPolymerMWEqIdx;
        static constexpr int contiFoamEqIdx = Indices::contiFoamEqIdx;
        static constexpr int contiBrineEqIdx = Indices::contiBrineEqIdx;
        static constexpr int contiMicrobialEqIdx = Indices::contiMicrobialEqIdx;
        static constexpr int contiOxygenEqIdx = Indices::contiOxygenEqIdx;
        static constexpr int contiUreaEqIdx = Indices::contiUreaEqIdx;
        static constexpr int contiBiofilmEqIdx = Indices::contiBiofilmEqIdx;
        static constexpr int contiCalciteEqIdx = Indices::contiCalciteEqIdx;
        static constexpr int solventSaturationIdx = Indices::solventSaturationIdx;
        static constexpr int zFractionIdx = Indices::zFractionIdx;
        static constexpr int polymerConcentrationIdx = Indices::polymerConcentrationIdx;
        static constexpr int polymerMoleWeightIdx = Indices::polymerMoleWeightIdx;
        static constexpr int temperatureIdx = Indices::temperatureIdx;
        static constexpr int foamConcentrationIdx = Indices::foamConcentrationIdx;
        static constexpr int saltConcentrationIdx = Indices::saltConcentrationIdx;
        static constexpr int microbialConcentrationIdx = Indices::microbialConcentrationIdx;
        static constexpr int oxygenConcentrationIdx = Indices::oxygenConcentrationIdx;
        static constexpr int ureaConcentrationIdx = Indices::ureaConcentrationIdx;
        static constexpr int biofilmConcentrationIdx = Indices::biofilmConcentrationIdx;
        static constexpr int calciteConcentrationIdx = Indices::calciteConcentrationIdx;

        using VectorBlockType = Dune::FieldVector<Scalar, numEq>;
        using MatrixBlockType = typename SparseMatrixAdapter::MatrixBlock;
        using Mat = typename SparseMatrixAdapter::IstlMatrix;
        using BVector = Dune::BlockVector<VectorBlockType>;

        using ComponentName = ::Opm::ComponentName<FluidSystem,Indices>;

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
        BlackoilModel(Simulator& simulator,
                      const ModelParameters& param,
                      BlackoilWellModel<TypeTag>& well_model,
                      const bool terminal_output)
        : simulator_(simulator)
        , grid_(simulator_.vanguard().grid())
        , phaseUsage_(phaseUsageFromDeck(eclState()))
        , param_( param )
        , well_model_ (well_model)
        , rst_conv_(simulator_.problem().eclWriter()->collectToIORank().localIdxToGlobalIdxMapping(),
                    grid_.comm())
        , terminal_output_ (terminal_output)
        , current_relaxation_(1.0)
        , dx_old_(simulator_.model().numGridDof())
        {
            // compute global sum of number of cells
            global_nc_ = detail::countGlobalCells(grid_);
            convergence_reports_.reserve(300); // Often insufficient, but avoids frequent moves.
            // TODO: remember to fix!
            if (param_.nonlinear_solver_ == "nldd") {
                if (!param_.matrix_add_well_contributions_) {
                    OPM_THROW(std::runtime_error, "The --nonlinear-solver=nldd option can only be used with --matrix-add-well-contributions=true");
                }
                if (terminal_output) {
                    OpmLog::info("Using Non-Linear Domain Decomposition solver (nldd).");
                }
                nlddSolver_ = std::make_unique<BlackoilModelNldd<TypeTag>>(*this);
            } else if (param_.nonlinear_solver_ == "newton") {
                if (terminal_output) {
                    OpmLog::info("Using Newton nonlinear solver.");
                }
            } else {
                OPM_THROW(std::runtime_error, "Unknown nonlinear solver option: " + param_.nonlinear_solver_);
            }
        }


        bool isParallel() const
        { return  grid_.comm().size() > 1; }


        const EclipseState& eclState() const
        { return simulator_.vanguard().eclState(); }


        /// Called once before each time step.
        /// \param[in] timer                  simulation timer
        SimulatorReportSingle prepareStep(const SimulatorTimerInterface& timer)
        {
            SimulatorReportSingle report;
            Dune::Timer perfTimer;
            perfTimer.start();
            // update the solution variables in the model
            if ( timer.lastStepFailed() ) {
                simulator_.model().updateFailed();
            } else {
                simulator_.model().advanceTimeLevel();
            }

            // Set the timestep size, episode index, and non-linear iteration index
            // for the model explicitly. The model needs to know the report step/episode index
            // because of timing dependent data despite the fact that flow uses its
            // own time stepper. (The length of the episode does not matter, though.)
            simulator_.setTime(timer.simulationTimeElapsed());
            simulator_.setTimeStepSize(timer.currentStepLength());
            simulator_.model().newtonMethod().setIterationIndex(0);

            simulator_.problem().beginTimeStep();

            unsigned numDof = simulator_.model().numGridDof();
            wasSwitched_.resize(numDof);
            std::fill(wasSwitched_.begin(), wasSwitched_.end(), false);

            if (param_.update_equations_scaling_) {
                OpmLog::error("Equation scaling not supported");
                //updateEquationsScaling();
            }

            if (nlddSolver_) {
                nlddSolver_->prepareStep();
            }

            report.pre_post_time += perfTimer.stop();

            auto getIdx = [](unsigned phaseIdx) -> int
            {
                if (FluidSystem::phaseIsActive(phaseIdx)) {
                    const unsigned sIdx = FluidSystem::solventComponentIndex(phaseIdx);
                    return Indices::canonicalToActiveComponentIndex(sIdx);
                }

                return -1;
            };
            const auto& schedule = simulator_.vanguard().schedule();
            rst_conv_.init(simulator_.vanguard().globalNumCells(),
                           schedule[timer.reportStepNum()].rst_config(),
                           {getIdx(FluidSystem::oilPhaseIdx),
                            getIdx(FluidSystem::gasPhaseIdx),
                            getIdx(FluidSystem::waterPhaseIdx),
                            contiPolymerEqIdx,
                            contiBrineEqIdx,
                            contiSolventEqIdx});

            return report;
        }


        void initialLinearization(SimulatorReportSingle& report,
                                  const int iteration,
                                  const int minIter,
                                  const int maxIter,
                                  const SimulatorTimerInterface& timer)
        {
            // -----------   Set up reports and timer   -----------
            failureReport_ = SimulatorReportSingle();
            Dune::Timer perfTimer;

            perfTimer.start();
            report.total_linearizations = 1;

            // -----------   Assemble   -----------
            try {
                report += assembleReservoir(timer, iteration);
                report.assemble_time += perfTimer.stop();
            }
            catch (...) {
                report.assemble_time += perfTimer.stop();
                failureReport_ += report;
                throw; // continue throwing the stick
            }

            // -----------   Check if converged   -----------
            std::vector<double> residual_norms;
            perfTimer.reset();
            perfTimer.start();
            // the step is not considered converged until at least minIter iterations is done
            {
                auto convrep = getConvergence(timer, iteration, maxIter, residual_norms);
                report.converged = convrep.converged() && iteration >= minIter;
                ConvergenceReport::Severity severity = convrep.severityOfWorstFailure();
                convergence_reports_.back().report.push_back(std::move(convrep));

                // Throw if any NaN or too large residual found.
                if (severity == ConvergenceReport::Severity::NotANumber) {
                    OPM_THROW(NumericalProblem, "NaN residual found!");
                } else if (severity == ConvergenceReport::Severity::TooLarge) {
                    OPM_THROW_NOLOG(NumericalProblem, "Too large residual found!");
                }
            }
            report.update_time += perfTimer.stop();
            residual_norms_history_.push_back(residual_norms);
        }


        /// Called once per nonlinear iteration.
        /// This model will perform a Newton-Raphson update, changing reservoir_state
        /// and well_state. It will also use the nonlinear_solver to do relaxation of
        /// updates if necessary.
        /// \param[in] iteration              should be 0 for the first call of a new timestep
        /// \param[in] timer                  simulation timer
        /// \param[in] nonlinear_solver       nonlinear solver used (for oscillation/relaxation control)
        template <class NonlinearSolverType>
        SimulatorReportSingle nonlinearIteration(const int iteration,
                                                 const SimulatorTimerInterface& timer,
                                                 NonlinearSolverType& nonlinear_solver)
        {
            if (iteration == 0) {
                // For each iteration we store in a vector the norms of the residual of
                // the mass balance for each active phase, the well flux and the well equations.
                residual_norms_history_.clear();
                current_relaxation_ = 1.0;
                dx_old_ = 0.0;
                convergence_reports_.push_back({timer.reportStepNum(), timer.currentStepNum(), {}});
                convergence_reports_.back().report.reserve(11);
            }

            SimulatorReportSingle result;
            if ((this->param_.nonlinear_solver_ != "nldd") ||
                (iteration < this->param_.nldd_num_initial_newton_iter_))
            {
                result = this->nonlinearIterationNewton(iteration, timer, nonlinear_solver);
            }
            else {
                result = this->nlddSolver_->nonlinearIterationNldd(iteration, timer, nonlinear_solver);
            }

            rst_conv_.update(simulator_.model().linearizer().residual());

            return result;
        }


        template <class NonlinearSolverType>
        SimulatorReportSingle nonlinearIterationNewton(const int iteration,
                                                       const SimulatorTimerInterface& timer,
                                                       NonlinearSolverType& nonlinear_solver)
        {

            // -----------   Set up reports and timer   -----------
            SimulatorReportSingle report;
            Dune::Timer perfTimer;

            this->initialLinearization(report, iteration, nonlinear_solver.minIter(), nonlinear_solver.maxIter(), timer);

            // -----------   If not converged, solve linear system and do Newton update  -----------
            if (!report.converged) {
                perfTimer.reset();
                perfTimer.start();
                report.total_newton_iterations = 1;

                // Compute the nonlinear update.
                unsigned nc = simulator_.model().numGridDof();
                BVector x(nc);

                // Solve the linear system.
                linear_solve_setup_time_ = 0.0;
                try {
                    // Apply the Schur complement of the well model to
                    // the reservoir linearized equations.
                    // Note that linearize may throw for MSwells.
                    wellModel().linearize(simulator().model().linearizer().jacobian(),
                                          simulator().model().linearizer().residual());

                    // ---- Solve linear system ----
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
                    nonlinear_solver.detectOscillations(residual_norms_history_, residual_norms_history_.size() - 1, isOscillate, isStagnate);
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

                // ---- Newton update ----
                // Apply the update, with considering model-dependent limitations and
                // chopping of the update.
                updateSolution(x);

                report.update_time += perfTimer.stop();
            }

            return report;
        }


        /// Called once after each time step.
        /// In this class, this function does nothing.
        /// \param[in] timer                  simulation timer
        SimulatorReportSingle afterStep(const SimulatorTimerInterface&)
        {
            SimulatorReportSingle report;
            Dune::Timer perfTimer;
            perfTimer.start();
            simulator_.problem().endTimeStep();
            simulator_.problem().setConvData(rst_conv_.getData());
            report.pre_post_time += perfTimer.stop();
            return report;
        }

        /// Assemble the residual and Jacobian of the nonlinear system.
        SimulatorReportSingle assembleReservoir(const SimulatorTimerInterface& /* timer */,
                                                const int iterationIdx)
        {
            // -------- Mass balance equations --------
            simulator_.model().newtonMethod().setIterationIndex(iterationIdx);
            simulator_.problem().beginIteration();
            simulator_.model().linearizer().linearizeDomain();
            simulator_.problem().endIteration();
            return wellModel().lastReport();
        }

        // compute the "relative" change of the solution between time steps
        double relativeChange() const
        {
            Scalar resultDelta = 0.0;
            Scalar resultDenom = 0.0;

            const auto& elemMapper = simulator_.model().elementMapper();
            const auto& gridView = simulator_.gridView();
            for (const auto& elem : elements(gridView, Dune::Partitions::interior)) {
                unsigned globalElemIdx = elemMapper.index(elem);
                const auto& priVarsNew = simulator_.model().solution(/*timeIdx=*/0)[globalElemIdx];

                Scalar pressureNew;
                pressureNew = priVarsNew[Indices::pressureSwitchIdx];

                Scalar saturationsNew[FluidSystem::numPhases] = { 0.0 };
                Scalar oilSaturationNew = 1.0;
                if (FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx) &&
                    FluidSystem::numActivePhases() > 1 &&
                    priVarsNew.primaryVarsMeaningWater() == PrimaryVariables::WaterMeaning::Sw) {
                    saturationsNew[FluidSystem::waterPhaseIdx] = priVarsNew[Indices::waterSwitchIdx];
                    oilSaturationNew -= saturationsNew[FluidSystem::waterPhaseIdx];
                }

                if (FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx) &&
                    FluidSystem::phaseIsActive(FluidSystem::oilPhaseIdx) &&
                    priVarsNew.primaryVarsMeaningGas() == PrimaryVariables::GasMeaning::Sg) {
                    assert(Indices::compositionSwitchIdx >= 0 );
                    saturationsNew[FluidSystem::gasPhaseIdx] = priVarsNew[Indices::compositionSwitchIdx];
                    oilSaturationNew -= saturationsNew[FluidSystem::gasPhaseIdx];
                }

                if (FluidSystem::phaseIsActive(FluidSystem::oilPhaseIdx)) {
                    saturationsNew[FluidSystem::oilPhaseIdx] = oilSaturationNew;
                }

                const auto& priVarsOld = simulator_.model().solution(/*timeIdx=*/1)[globalElemIdx];

                Scalar pressureOld;
                pressureOld = priVarsOld[Indices::pressureSwitchIdx];

                Scalar saturationsOld[FluidSystem::numPhases] = { 0.0 };
                Scalar oilSaturationOld = 1.0;

                // NB fix me! adding pressures changes to satutation changes does not make sense
                Scalar tmp = pressureNew - pressureOld;
                resultDelta += tmp*tmp;
                resultDenom += pressureNew*pressureNew;

                if (FluidSystem::numActivePhases() > 1) {
                    if (priVarsOld.primaryVarsMeaningWater() == PrimaryVariables::WaterMeaning::Sw) {
                        saturationsOld[FluidSystem::waterPhaseIdx] = priVarsOld[Indices::waterSwitchIdx];
                        oilSaturationOld -= saturationsOld[FluidSystem::waterPhaseIdx];
                    }

                    if (priVarsOld.primaryVarsMeaningGas() == PrimaryVariables::GasMeaning::Sg)
                    {
                        assert(Indices::compositionSwitchIdx >= 0 );
                        saturationsOld[FluidSystem::gasPhaseIdx] = priVarsOld[Indices::compositionSwitchIdx];
                        oilSaturationOld -= saturationsOld[FluidSystem::gasPhaseIdx];
                    }

                    if (FluidSystem::phaseIsActive(FluidSystem::oilPhaseIdx)) {
                        saturationsOld[FluidSystem::oilPhaseIdx] = oilSaturationOld;
                    }
                    for (unsigned phaseIdx = 0; phaseIdx < FluidSystem::numPhases; ++ phaseIdx) {
                        Scalar tmpSat = saturationsNew[phaseIdx] - saturationsOld[phaseIdx];
                        resultDelta += tmpSat*tmpSat;
                        resultDenom += saturationsNew[phaseIdx]*saturationsNew[phaseIdx];
                        assert(std::isfinite(resultDelta));
                        assert(std::isfinite(resultDenom));
                    }
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
            return simulator_.model().newtonMethod().linearSolver().iterations ();
        }


        // Obtain reference to linear solver setup time
        double& linearSolveSetupTime()
        {
            return linear_solve_setup_time_;
        }


        /// Solve the Jacobian system Jx = r where J is the Jacobian and
        /// r is the residual.
        void solveJacobianSystem(BVector& x)
        {
            auto& jacobian = simulator_.model().linearizer().jacobian().istlMatrix();
            auto& residual = simulator_.model().linearizer().residual();
            auto& linSolver = simulator_.model().newtonMethod().linearSolver();

            const int numSolvers = linSolver.numAvailableSolvers();
            if ((numSolvers > 1) && (linSolver.getSolveCount() % 100 == 0)) {

                if ( terminal_output_ ) {
                    OpmLog::debug("\nRunning speed test for comparing available linear solvers.");
                }

                Dune::Timer perfTimer;
                std::vector<double> times(numSolvers);
                std::vector<double> setupTimes(numSolvers);

                x = 0.0;
                std::vector<BVector> x_trial(numSolvers, x);
                for (int solver = 0; solver < numSolvers; ++solver) {
                    BVector x0(x);
                    linSolver.setActiveSolver(solver);
                    perfTimer.start();
                    linSolver.prepare(jacobian, residual);
                    setupTimes[solver] = perfTimer.stop();
                    perfTimer.reset();
                    linSolver.setResidual(residual);
                    perfTimer.start();
                    linSolver.solve(x_trial[solver]);
                    times[solver] = perfTimer.stop();
                    perfTimer.reset();
                    if (terminal_output_) {
                        OpmLog::debug(fmt::format("Solver time {}: {}", solver, times[solver]));
                    }
                }

                int fastest_solver = std::min_element(times.begin(), times.end()) - times.begin();
                // Use timing on rank 0 to determine fastest, must be consistent across ranks.
                grid_.comm().broadcast(&fastest_solver, 1, 0);
                linear_solve_setup_time_ = setupTimes[fastest_solver];
                x = x_trial[fastest_solver];
                linSolver.setActiveSolver(fastest_solver);
            } else {
                // set initial guess
                x = 0.0;

                Dune::Timer perfTimer;
                perfTimer.start();
                linSolver.prepare(jacobian, residual);
                linear_solve_setup_time_ = perfTimer.stop();
                linSolver.setResidual(residual);
                // actually, the error needs to be calculated after setResidual in order to
                // account for parallelization properly. since the residual of ECFV
                // discretizations does not need to be synchronized across processes to be
                // consistent, this is not relevant for OPM-flow...
                linSolver.solve(x);
            }
       }


        /// Apply an update to the primary variables.
        void updateSolution(const BVector& dx)
        {
            OPM_TIMEBLOCK(updateSolution);
            auto& newtonMethod = simulator_.model().newtonMethod();
            SolutionVector& solution = simulator_.model().solution(/*timeIdx=*/0);

            newtonMethod.update_(/*nextSolution=*/solution,
                                 /*curSolution=*/solution,
                                 /*update=*/dx,
                                 /*resid=*/dx); // the update routines of the black
                                                // oil model do not care about the
                                                // residual

            // if the solution is updated, the intensive quantities need to be recalculated
            {
                OPM_TIMEBLOCK(invalidateAndUpdateIntensiveQuantities);
                simulator_.model().invalidateAndUpdateIntensiveQuantities(/*timeIdx=*/0);
                simulator_.problem().eclWriter()->mutableEclOutputModule().invalidateLocalData();
            }
        }

        /// Return true if output to cout is wanted.
        bool terminalOutputEnabled() const
        {
            return terminal_output_;
        }

        std::tuple<double,double> convergenceReduction(Parallel::Communication comm,
                                                       const double pvSumLocal,
                                                       const double numAquiferPvSumLocal,
                                                       std::vector< Scalar >& R_sum,
                                                       std::vector< Scalar >& maxCoeff,
                                                       std::vector< Scalar >& B_avg)
        {
            OPM_TIMEBLOCK(convergenceReduction);
            // Compute total pore volume (use only owned entries)
            double pvSum = pvSumLocal;
            double numAquiferPvSum = numAquiferPvSumLocal;

            if( comm.size() > 1 )
            {
                // global reduction
                std::vector< Scalar > sumBuffer;
                std::vector< Scalar > maxBuffer;
                const int numComp = B_avg.size();
                sumBuffer.reserve( 2*numComp + 2 ); // +2 for (numAquifer)pvSum
                maxBuffer.reserve( numComp );
                for( int compIdx = 0; compIdx < numComp; ++compIdx )
                {
                    sumBuffer.push_back( B_avg[ compIdx ] );
                    sumBuffer.push_back( R_sum[ compIdx ] );
                    maxBuffer.push_back( maxCoeff[ compIdx ] );
                }

                // Compute total pore volume
                sumBuffer.push_back( pvSum );
                sumBuffer.push_back( numAquiferPvSum );

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
                pvSum = sumBuffer[sumBuffer.size()-2];
                numAquiferPvSum = sumBuffer.back();
            }

            // return global pore volume
            return {pvSum, numAquiferPvSum};
        }

        /// \brief Get reservoir quantities on this process needed for convergence calculations.
        /// \return A pair of the local pore volume of interior cells and the pore volumes
        ///         of the cells associated with a numerical aquifer.
        std::pair<double,double> localConvergenceData(std::vector<Scalar>& R_sum,
                                                      std::vector<Scalar>& maxCoeff,
                                                      std::vector<Scalar>& B_avg,
                                                      std::vector<int>& maxCoeffCell)
        {
            OPM_TIMEBLOCK(localConvergenceData);
            double pvSumLocal = 0.0;
            double numAquiferPvSumLocal = 0.0;
            const auto& model = simulator_.model();
            const auto& problem = simulator_.problem();

            const auto& residual = simulator_.model().linearizer().residual();

            ElementContext elemCtx(simulator_);
            const auto& gridView = simulator().gridView();
            IsNumericalAquiferCell isNumericalAquiferCell(gridView.grid());
            OPM_BEGIN_PARALLEL_TRY_CATCH();
            for (const auto& elem : elements(gridView, Dune::Partitions::interior)) {
                elemCtx.updatePrimaryStencil(elem);
                elemCtx.updatePrimaryIntensiveQuantities(/*timeIdx=*/0);

                const unsigned cell_idx = elemCtx.globalSpaceIndex(/*spaceIdx=*/0, /*timeIdx=*/0);
                const auto& intQuants = elemCtx.intensiveQuantities(/*spaceIdx=*/0, /*timeIdx=*/0);
                const auto& fs = intQuants.fluidState();

                const auto pvValue = problem.referencePorosity(cell_idx, /*timeIdx=*/0) *
                                     model.dofTotalVolume(cell_idx);
                pvSumLocal += pvValue;

                if (isNumericalAquiferCell(elem))
                {
                    numAquiferPvSumLocal += pvValue;
                }

                this->getMaxCoeff(cell_idx, intQuants, fs, residual, pvValue,
                                  B_avg, R_sum, maxCoeff, maxCoeffCell);
            }

            OPM_END_PARALLEL_TRY_CATCH("BlackoilModel::localConvergenceData() failed: ", grid_.comm());

            // compute local average in terms of global number of elements
            const int bSize = B_avg.size();
            for ( int i = 0; i<bSize; ++i )
            {
                B_avg[ i ] /= Scalar( global_nc_ );
            }

            return {pvSumLocal, numAquiferPvSumLocal};
        }


        /// \brief Compute the total pore volume of cells violating CNV that are not part
        ///        of a numerical aquifer.
        double computeCnvErrorPv(const std::vector<Scalar>& B_avg, double dt)
        {
            OPM_TIMEBLOCK(computeCnvErrorPv);
            double errorPV{};
            const auto& model = simulator_.model();
            const auto& problem = simulator_.problem();
            const auto& residual = simulator_.model().linearizer().residual();
            const auto& gridView = simulator().gridView();
            ElementContext elemCtx(simulator_);
            IsNumericalAquiferCell isNumericalAquiferCell(gridView.grid());

            OPM_BEGIN_PARALLEL_TRY_CATCH();
            for (const auto& elem : elements(gridView, Dune::Partitions::interiorBorder))
            {
                // Skip cells of numerical Aquifer
                if (isNumericalAquiferCell(elem))
                {
                    continue;
                }
                elemCtx.updatePrimaryStencil(elem);
                // elemCtx.updatePrimaryIntensiveQuantities(/*timeIdx=*/0);
                const unsigned cell_idx = elemCtx.globalSpaceIndex(/*spaceIdx=*/0, /*timeIdx=*/0);
                const double pvValue = problem.referencePorosity(cell_idx, /*timeIdx=*/0) * model.dofTotalVolume(cell_idx);
                const auto& cellResidual = residual[cell_idx];
                bool cnvViolated = false;

                for (unsigned eqIdx = 0; eqIdx < cellResidual.size(); ++eqIdx)
                {
                    using std::abs;
                    Scalar CNV = cellResidual[eqIdx] * dt * B_avg[eqIdx] / pvValue;
                    cnvViolated = cnvViolated || (abs(CNV) > param_.tolerance_cnv_);
                }

                if (cnvViolated)
                {
                    errorPV += pvValue;
                }
            }

            OPM_END_PARALLEL_TRY_CATCH("BlackoilModel::ComputeCnvError() failed: ", grid_.comm());

            return grid_.comm().sum(errorPV);
        }


        void updateTUNING(const Tuning& tuning) {
            param_.tolerance_mb_ = tuning.XXXMBE;
            if ( terminal_output_ ) {
                OpmLog::debug(fmt::format("Setting BlackoilModel mass balance limit (XXXMBE) to {:.2e}", tuning.XXXMBE));
            }
        }


        ConvergenceReport getReservoirConvergence(const double reportTime,
                                                  const double dt,
                                                  const int iteration,
                                                  const int maxIter,
                                                  std::vector<Scalar>& B_avg,
                                                  std::vector<Scalar>& residual_norms)
        {
            OPM_TIMEBLOCK(getReservoirConvergence);
            using Vector = std::vector<Scalar>;

            const int numComp = numEq;
            Vector R_sum(numComp, 0.0 );
            Vector maxCoeff(numComp, std::numeric_limits< Scalar >::lowest() );
            std::vector<int> maxCoeffCell(numComp, -1);
            const auto [ pvSumLocal, numAquiferPvSumLocal] = localConvergenceData(R_sum, maxCoeff, B_avg, maxCoeffCell);

            // compute global sum and max of quantities
            const auto [ pvSum, numAquiferPvSum ] =
                convergenceReduction(grid_.comm(), pvSumLocal,
                                     numAquiferPvSumLocal,
                                     R_sum, maxCoeff, B_avg);

            auto cnvErrorPvFraction = computeCnvErrorPv(B_avg, dt);
            cnvErrorPvFraction /= (pvSum - numAquiferPvSum);


            // For each iteration, we need to determine whether to use the relaxed tolerances.
            // To disable the usage of relaxed tolerances, you can set the relaxed tolerances as the strict tolerances.

            // If min_strict_mb_iter = -1 (default) we use a relaxed tolerance for the mass balance for the last iterations
            // For positive values we use the relaxed tolerance after the given number of iterations
            const bool relax_final_iteration_mb = (param_.min_strict_mb_iter_ < 0) && (iteration == maxIter);
            const bool use_relaxed_mb = relax_final_iteration_mb || (param_.min_strict_mb_iter_ >= 0 && iteration >= param_.min_strict_mb_iter_);


            // If min_strict_cnv_iter = -1 we use a relaxed tolerance for the cnv for the last iterations
            // For positive values we use the relaxed tolerance after the given number of iterations
            // We also use relaxed tolerances for cells with total poro volume less than relaxed_max_pv_fraction_
            // Default value of relaxed_max_pv_fraction_ is 0.03
            const bool relax_final_iteration_cnv = (param_.min_strict_cnv_iter_ < 0) && (iteration == maxIter);
            const bool relax_iter_cnv = param_.min_strict_mb_iter_ >= 0 && iteration >= param_.min_strict_mb_iter_;
            const bool relax_pv_fraction_cnv = cnvErrorPvFraction < param_.relaxed_max_pv_fraction_;
            const bool use_relaxed_cnv = relax_final_iteration_cnv || relax_pv_fraction_cnv || relax_iter_cnv;

            if  (relax_final_iteration_mb || relax_final_iteration_cnv)  {
                if ( terminal_output_ ) {
                    std::string message = "Number of newton iterations reached its maximum try to continue with relaxed tolerances:";
                    if (relax_final_iteration_mb)
                        message += fmt::format("  MB: {:.1e}", param_.tolerance_mb_relaxed_);
                    if (relax_final_iteration_cnv)
                        message += fmt::format(" CNV: {:.1e}", param_.tolerance_cnv_relaxed_);

                    OpmLog::debug(message);
                }
            }
            const double tol_cnv = use_relaxed_cnv ? param_.tolerance_cnv_relaxed_ :  param_.tolerance_cnv_;
            const double tol_mb  = use_relaxed_mb ? param_.tolerance_mb_relaxed_ : param_.tolerance_mb_;

            // Finish computation
            std::vector<Scalar> CNV(numComp);
            std::vector<Scalar> mass_balance_residual(numComp);
            for ( int compIdx = 0; compIdx < numComp; ++compIdx )
            {
                CNV[compIdx]                    = B_avg[compIdx] * dt * maxCoeff[compIdx];
                mass_balance_residual[compIdx]  = std::abs(B_avg[compIdx]*R_sum[compIdx]) * dt / pvSum;
                residual_norms.push_back(CNV[compIdx]);
            }

            // Create convergence report.
            ConvergenceReport report{reportTime};
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
                            OpmLog::debug("NaN residual for " + this->compNames_.name(compIdx) + " equation.");
                        }
                    } else if (res[ii] > maxResidualAllowed()) {
                        report.setReservoirFailed({types[ii], CR::Severity::TooLarge, compIdx});
                        if ( terminal_output_ ) {
                            OpmLog::debug("Too large residual for " + this->compNames_.name(compIdx) + " equation.");
                        }
                    } else if (res[ii] < 0.0) {
                        report.setReservoirFailed({types[ii], CR::Severity::Normal, compIdx});
                        if ( terminal_output_ ) {
                            OpmLog::debug("Negative residual for " + this->compNames_.name(compIdx) + " equation.");
                        }
                    } else if (res[ii] > tol[ii]) {
                        report.setReservoirFailed({types[ii], CR::Severity::Normal, compIdx});
                    }
                    report.setReservoirConvergenceMetric(types[ii], compIdx, res[ii]);
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
                        msg += this->compNames_.name(compIdx)[0];
                        msg += ")  ";
                    }
                    for (int compIdx = 0; compIdx < numComp; ++compIdx) {
                        msg += "    CNV(";
                        msg += this->compNames_.name(compIdx)[0];
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
        /// \param[in]   maxIter     maximum number of iterations
        /// \param[out]  residual_norms   CNV residuals by phase
        ConvergenceReport getConvergence(const SimulatorTimerInterface& timer,
                                         const int iteration,
                                         const int maxIter,
                                         std::vector<double>& residual_norms)
        {
            OPM_TIMEBLOCK(getConvergence);
            // Get convergence reports for reservoir and wells.
            std::vector<Scalar> B_avg(numEq, 0.0);
            auto report = getReservoirConvergence(timer.simulationTimeElapsed(),
                                                  timer.currentStepLength(),
                                                  iteration, maxIter, B_avg, residual_norms);
            {
                OPM_TIMEBLOCK(getWellConvergence);
                report += wellModel().getWellConvergence(B_avg, /*checkWellGroupControls*/report.converged());
            }
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
            OPM_TIMEBLOCK(computeFluidInPlace);
            //assert(true)
            //return an empty vector
            std::vector<std::vector<double> > regionValues(0, std::vector<double>(0,0.0));
            return regionValues;
        }

        const Simulator& simulator() const
        { return simulator_; }

        Simulator& simulator()
        { return simulator_; }

        /// return the statistics if the nonlinearIteration() method failed
        const SimulatorReportSingle& failureReport() const
        { return failureReport_; }

        /// return the statistics if the nonlinearIteration() method failed
        SimulatorReportSingle localAccumulatedReports() const
        {
            return nlddSolver_ ? nlddSolver_->localAccumulatedReports()
                               : SimulatorReportSingle{};
        }

        const std::vector<StepReport>& stepReports() const
        {
            return convergence_reports_;
        }

        void writePartitions(const std::filesystem::path& odir) const
        {
            if (this->nlddSolver_ != nullptr) {
                this->nlddSolver_->writePartitions(odir);
                return;
            }

            const auto& elementMapper = this->simulator().model().elementMapper();
            const auto& cartMapper = this->simulator().vanguard().cartesianIndexMapper();

            const auto& grid = this->simulator().vanguard().grid();
            const auto& comm = grid.comm();
            const auto nDigit = 1 + static_cast<int>(std::floor(std::log10(comm.size())));

            std::ofstream pfile { odir / fmt::format("{1:0>{0}}", nDigit, comm.rank()) };

            for (const auto& cell : elements(grid.leafGridView(), Dune::Partitions::interior)) {
                pfile << comm.rank() << ' '
                      << cartMapper.cartesianIndex(elementMapper.index(cell)) << ' '
                      << comm.rank() << '\n';
            }
        }

        const std::vector<std::vector<int>>& getConvCells() const
        { return rst_conv_.getData(); }

    protected:
        // ---------  Data members  ---------

        Simulator& simulator_;
        const Grid& grid_;
        const PhaseUsage phaseUsage_;
        static constexpr bool has_solvent_ = getPropValue<TypeTag, Properties::EnableSolvent>();
        static constexpr bool has_extbo_ = getPropValue<TypeTag, Properties::EnableExtbo>();
        static constexpr bool has_polymer_ = getPropValue<TypeTag, Properties::EnablePolymer>();
        static constexpr bool has_polymermw_ = getPropValue<TypeTag, Properties::EnablePolymerMW>();
        static constexpr bool has_energy_ = getPropValue<TypeTag, Properties::EnableEnergy>();
        static constexpr bool has_foam_ = getPropValue<TypeTag, Properties::EnableFoam>();
        static constexpr bool has_brine_ = getPropValue<TypeTag, Properties::EnableBrine>();
        static constexpr bool has_micp_ = getPropValue<TypeTag, Properties::EnableMICP>();

        ModelParameters                 param_;
        SimulatorReportSingle failureReport_;

        // Well Model
        BlackoilWellModel<TypeTag>& well_model_;

        RSTConv rst_conv_; //!< Helper class for RPTRST CONV

        /// \brief Whether we print something to std::cout
        bool terminal_output_;
        /// \brief The number of cells of the global grid.
        long int global_nc_;

        std::vector<std::vector<double>> residual_norms_history_;
        double current_relaxation_;
        BVector dx_old_;

        std::vector<StepReport> convergence_reports_;
        ComponentName compNames_{};

        std::unique_ptr<BlackoilModelNldd<TypeTag>> nlddSolver_; //!< Non-linear DD solver

    public:
        /// return the StandardWells object
        BlackoilWellModel<TypeTag>&
        wellModel() { return well_model_; }

        const BlackoilWellModel<TypeTag>&
        wellModel() const { return well_model_; }

        void beginReportStep()
        {
            simulator_.problem().beginEpisode();
        }

        void endReportStep()
        {
            simulator_.problem().endEpisode();
        }

        template<class FluidState, class Residual>
        void getMaxCoeff(const unsigned cell_idx,
                         const IntensiveQuantities& intQuants,
                         const FluidState& fs,
                         const Residual& modelResid,
                         const Scalar pvValue,
                         std::vector<Scalar>& B_avg,
                         std::vector<Scalar>& R_sum,
                         std::vector<Scalar>& maxCoeff,
                         std::vector<int>& maxCoeffCell)
        {
            for (unsigned phaseIdx = 0; phaseIdx < FluidSystem::numPhases; ++phaseIdx)
            {
                if (!FluidSystem::phaseIsActive(phaseIdx)) {
                    continue;
                }

                const unsigned compIdx = Indices::canonicalToActiveComponentIndex(FluidSystem::solventComponentIndex(phaseIdx));

                B_avg[compIdx] += 1.0 / fs.invB(phaseIdx).value();
                const auto R2 = modelResid[cell_idx][compIdx];

                R_sum[compIdx] += R2;
                const double Rval = std::abs(R2) / pvValue;
                if (Rval > maxCoeff[compIdx]) {
                    maxCoeff[compIdx] = Rval;
                    maxCoeffCell[compIdx] = cell_idx;
                }
            }

            if constexpr (has_solvent_) {
                B_avg[contiSolventEqIdx] += 1.0 / intQuants.solventInverseFormationVolumeFactor().value();
                const auto R2 = modelResid[cell_idx][contiSolventEqIdx];
                R_sum[contiSolventEqIdx] += R2;
                maxCoeff[contiSolventEqIdx] = std::max(maxCoeff[contiSolventEqIdx],
                                                       std::abs(R2) / pvValue);
            }
            if constexpr (has_extbo_) {
                B_avg[contiZfracEqIdx] += 1.0 / fs.invB(FluidSystem::gasPhaseIdx).value();
                const auto R2 = modelResid[cell_idx][contiZfracEqIdx];
                R_sum[ contiZfracEqIdx ] += R2;
                maxCoeff[contiZfracEqIdx] = std::max(maxCoeff[contiZfracEqIdx],
                                                     std::abs(R2) / pvValue);
            }
            if constexpr (has_polymer_) {
                B_avg[contiPolymerEqIdx] += 1.0 / fs.invB(FluidSystem::waterPhaseIdx).value();
                const auto R2 = modelResid[cell_idx][contiPolymerEqIdx];
                R_sum[contiPolymerEqIdx] += R2;
                maxCoeff[contiPolymerEqIdx] = std::max(maxCoeff[contiPolymerEqIdx],
                                                       std::abs(R2) / pvValue);
            }
            if constexpr (has_foam_) {
                B_avg[ contiFoamEqIdx ] += 1.0 / fs.invB(FluidSystem::gasPhaseIdx).value();
                const auto R2 = modelResid[cell_idx][contiFoamEqIdx];
                R_sum[contiFoamEqIdx] += R2;
                maxCoeff[contiFoamEqIdx] = std::max(maxCoeff[contiFoamEqIdx],
                                                    std::abs(R2) / pvValue);
            }
            if constexpr (has_brine_) {
                B_avg[ contiBrineEqIdx ] += 1.0 / fs.invB(FluidSystem::waterPhaseIdx).value();
                const auto R2 = modelResid[cell_idx][contiBrineEqIdx];
                R_sum[contiBrineEqIdx] += R2;
                maxCoeff[contiBrineEqIdx] = std::max(maxCoeff[contiBrineEqIdx],
                                                     std::abs(R2) / pvValue);
            }

            if constexpr (has_polymermw_) {
                static_assert(has_polymer_);

                B_avg[contiPolymerMWEqIdx] += 1.0 / fs.invB(FluidSystem::waterPhaseIdx).value();
                // the residual of the polymer molecular equation is scaled down by a 100, since molecular weight
                // can be much bigger than 1, and this equation shares the same tolerance with other mass balance equations
                // TODO: there should be a more general way to determine the scaling-down coefficient
                const auto R2 = modelResid[cell_idx][contiPolymerMWEqIdx] / 100.;
                R_sum[contiPolymerMWEqIdx] += R2;
                maxCoeff[contiPolymerMWEqIdx] = std::max(maxCoeff[contiPolymerMWEqIdx],
                                                         std::abs(R2) / pvValue);
            }

            if constexpr (has_energy_) {
                B_avg[contiEnergyEqIdx] += 1.0 / (4.182e1); // converting J -> RM3 (entalpy / (cp * deltaK * rho) assuming change of 1e-5K of water
                const auto R2 = modelResid[cell_idx][contiEnergyEqIdx];
                R_sum[contiEnergyEqIdx] += R2;
                maxCoeff[contiEnergyEqIdx] = std::max(maxCoeff[contiEnergyEqIdx],
                                                      std::abs(R2) / pvValue);
            }

            if constexpr (has_micp_) {
                B_avg[contiMicrobialEqIdx] += 1.0 / fs.invB(FluidSystem::waterPhaseIdx).value();
                const auto R1 = modelResid[cell_idx][contiMicrobialEqIdx];
                R_sum[contiMicrobialEqIdx] += R1;
                maxCoeff[contiMicrobialEqIdx] = std::max(maxCoeff[contiMicrobialEqIdx],
                                                         std::abs(R1) / pvValue);
                B_avg[contiOxygenEqIdx] += 1.0 / fs.invB(FluidSystem::waterPhaseIdx).value();
                const auto R2 = modelResid[cell_idx][contiOxygenEqIdx];
                R_sum[contiOxygenEqIdx] += R2;
                maxCoeff[contiOxygenEqIdx] = std::max(maxCoeff[contiOxygenEqIdx],
                                                      std::abs(R2) / pvValue);
                B_avg[contiUreaEqIdx] += 1.0 / fs.invB(FluidSystem::waterPhaseIdx).value();
                const auto R3 = modelResid[cell_idx][contiUreaEqIdx];
                R_sum[contiUreaEqIdx] += R3;
                maxCoeff[contiUreaEqIdx] = std::max(maxCoeff[contiUreaEqIdx],
                                                    std::abs(R3) / pvValue);
                B_avg[contiBiofilmEqIdx] += 1.0 / fs.invB(FluidSystem::waterPhaseIdx).value();
                const auto R4 = modelResid[cell_idx][contiBiofilmEqIdx];
                R_sum[contiBiofilmEqIdx] += R4;
                maxCoeff[contiBiofilmEqIdx] = std::max(maxCoeff[contiBiofilmEqIdx],
                                                       std::abs(R4) / pvValue);
                B_avg[contiCalciteEqIdx] += 1.0 / fs.invB(FluidSystem::waterPhaseIdx).value();
                const auto R5 = modelResid[cell_idx][contiCalciteEqIdx];
                R_sum[contiCalciteEqIdx] += R5;
                maxCoeff[contiCalciteEqIdx] = std::max(maxCoeff[contiCalciteEqIdx],
                                                       std::abs(R5) / pvValue);
            }
        }

        //! \brief Returns const reference to model parameters.
        const ModelParameters& param() const
        {
            return param_;
        }

        //! \brief Returns const reference to component names.
        const ComponentName& compNames() const
        {
            return compNames_;
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

#endif // OPM_BLACKOILMODEL_HEADER_INCLUDED
