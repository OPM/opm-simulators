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

#ifndef OPM_BLACKOILMODEL_IMPL_HEADER_INCLUDED
#define OPM_BLACKOILMODEL_IMPL_HEADER_INCLUDED

#ifndef OPM_BLACKOILMODEL_HEADER_INCLUDED
#include <config.h>
#include <opm/simulators/flow/BlackoilModel.hpp>
#endif

#include <dune/common/timer.hh>

#include <opm/common/ErrorMacros.hpp>
#include <opm/common/OpmLog/OpmLog.hpp>

#include <opm/simulators/flow/countGlobalCells.hpp>

#include <algorithm>
#include <iomanip>
#include <limits>
#include <stdexcept>
#include <sstream>

#include <fmt/format.h>

namespace Opm {

template <class TypeTag>
BlackoilModel<TypeTag>::
BlackoilModel(Simulator& simulator,
              const ModelParameters& param,
              BlackoilWellModel<TypeTag>& well_model,
              const bool terminal_output)
    : simulator_(simulator)
    , grid_(simulator_.vanguard().grid())
    , param_( param )
    , well_model_ (well_model)
    , terminal_output_ (terminal_output)
    , current_relaxation_(1.0)
    , dx_old_(simulator_.model().numGridDof())
    , conv_monitor_(param_.monitor_params_)
{
    // compute global sum of number of cells
    global_nc_ = detail::countGlobalCells(grid_);
    convergence_reports_.reserve(300); // Often insufficient, but avoids frequent moves.
    // TODO: remember to fix!
    if (param_.nonlinear_solver_ == "nldd") {
        if (terminal_output) {
            OpmLog::info("Using Non-Linear Domain Decomposition solver (nldd).");
        }
        nlddSolver_ = std::make_unique<BlackoilModelNldd<TypeTag>>(*this);
    } else if (param_.nonlinear_solver_ == "newton") {
        if (terminal_output) {
            OpmLog::info("Using Newton nonlinear solver.");
        }
    } else {
        OPM_THROW(std::runtime_error, "Unknown nonlinear solver option: " +
                                      param_.nonlinear_solver_);
    }
}

template <class TypeTag>
SimulatorReportSingle
BlackoilModel<TypeTag>::
prepareStep(const SimulatorTimerInterface& timer)
{
    SimulatorReportSingle report;
    Dune::Timer perfTimer;
    perfTimer.start();
    // update the solution variables in the model
    int lastStepFailed = timer.lastStepFailed();
    if (grid_.comm().size() > 1 && grid_.comm().max(lastStepFailed) != grid_.comm().min(lastStepFailed)) {
        OPM_THROW(std::runtime_error, "Misalignment of the parallel simulation run in prepareStep " +
                                "- the previous step succeeded on some ranks but failed on others.");
    }
    if (lastStepFailed) {
        simulator_.model().updateFailed();
    }
    else {
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

    if (hasNlddSolver()) {
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
    auto& rst_conv = simulator_.problem().eclWriter().mutableOutputModule().getConv();
    rst_conv.init(simulator_.vanguard().globalNumCells(),
                  schedule[timer.reportStepNum()].rst_config(),
                  {getIdx(FluidSystem::oilPhaseIdx),
                   getIdx(FluidSystem::gasPhaseIdx),
                   getIdx(FluidSystem::waterPhaseIdx),
                   contiPolymerEqIdx,
                   contiBrineEqIdx,
                   contiSolventEqIdx});

    return report;
}

template <class TypeTag>
void
BlackoilModel<TypeTag>::
initialLinearization(SimulatorReportSingle& report,
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
    std::vector<Scalar> residual_norms;
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
            failureReport_ += report;
            OPM_THROW_PROBLEM(NumericalProblem, "NaN residual found!");
        } else if (severity == ConvergenceReport::Severity::TooLarge) {
            failureReport_ += report;
            OPM_THROW_NOLOG(NumericalProblem, "Too large residual found!");
        } else if (severity == ConvergenceReport::Severity::ConvergenceMonitorFailure) {
            failureReport_ += report;
            OPM_THROW_PROBLEM(ConvergenceMonitorFailure,
                              fmt::format(
                                  "Total penalty count exceeded cut-off-limit of {}",
                                  param_.monitor_params_.cutoff_
                              ));
        }
    }
    report.update_time += perfTimer.stop();
    residual_norms_history_.push_back(residual_norms);
}

template <class TypeTag>
template <class NonlinearSolverType>
SimulatorReportSingle
BlackoilModel<TypeTag>::
nonlinearIteration(const int iteration,
                   const SimulatorTimerInterface& timer,
                   NonlinearSolverType& nonlinear_solver)
{
    if (iteration == 0) {
        // For each iteration we store in a vector the norms of the residual of
        // the mass balance for each active phase, the well flux and the well equations.
        residual_norms_history_.clear();
        conv_monitor_.reset();
        current_relaxation_ = 1.0;
        dx_old_ = 0.0;
        convergence_reports_.push_back({timer.reportStepNum(), timer.currentStepNum(), {}});
        convergence_reports_.back().report.reserve(11);
    }

    SimulatorReportSingle result;
    if (this->param_.nonlinear_solver_ != "nldd")
    {
        result = this->nonlinearIterationNewton(iteration,
                                                timer,
                                                nonlinear_solver);
    }
    else {
        result = this->nlddSolver_->nonlinearIterationNldd(iteration,
                                                           timer,
                                                           nonlinear_solver);
    }

    auto& rst_conv = simulator_.problem().eclWriter().mutableOutputModule().getConv();
    rst_conv.update(simulator_.model().linearizer().residual());

    return result;
}

template <class TypeTag>
template <class NonlinearSolverType>
SimulatorReportSingle
BlackoilModel<TypeTag>::
nonlinearIterationNewton(const int iteration,
                         const SimulatorTimerInterface& timer,
                         NonlinearSolverType& nonlinear_solver)
{
    // -----------   Set up reports and timer   -----------
    SimulatorReportSingle report;
    Dune::Timer perfTimer;

    this->initialLinearization(report,
                               iteration,
                               nonlinear_solver.minIter(),
                               nonlinear_solver.maxIter(),
                               timer);

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
            nonlinear_solver.detectOscillations(residual_norms_history_,
                                                residual_norms_history_.size() - 1,
                                                isOscillate,
                                                isStagnate);
            if (isOscillate) {
                current_relaxation_ -= nonlinear_solver.relaxIncrement();
                current_relaxation_ = std::max(current_relaxation_,
                                               nonlinear_solver.relaxMax());
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

template <class TypeTag>
SimulatorReportSingle
BlackoilModel<TypeTag>::
afterStep(const SimulatorTimerInterface&)
{
    SimulatorReportSingle report;
    Dune::Timer perfTimer;
    perfTimer.start();
    simulator_.problem().endTimeStep();
    report.pre_post_time += perfTimer.stop();
    return report;
}

template <class TypeTag>
SimulatorReportSingle
BlackoilModel<TypeTag>::
assembleReservoir(const SimulatorTimerInterface& /* timer */,
                  const int iterationIdx)
{
    // -------- Mass balance equations --------
    simulator_.model().newtonMethod().setIterationIndex(iterationIdx);
    simulator_.problem().beginIteration();
    simulator_.model().linearizer().linearizeDomain();
    simulator_.problem().endIteration();
    return wellModel().lastReport();
}

template <class TypeTag>
typename BlackoilModel<TypeTag>::Scalar
BlackoilModel<TypeTag>::
relativeChange() const
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
            priVarsNew.primaryVarsMeaningWater() == PrimaryVariables::WaterMeaning::Sw)
        {
            saturationsNew[FluidSystem::waterPhaseIdx] = priVarsNew[Indices::waterSwitchIdx];
            oilSaturationNew -= saturationsNew[FluidSystem::waterPhaseIdx];
        }

        if (FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx) &&
            FluidSystem::phaseIsActive(FluidSystem::oilPhaseIdx) &&
            priVarsNew.primaryVarsMeaningGas() == PrimaryVariables::GasMeaning::Sg)
        {
            assert(Indices::compositionSwitchIdx >= 0);
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

        // NB fix me! adding pressures changes to saturation changes does not make sense
        Scalar tmp = pressureNew - pressureOld;
        resultDelta += tmp*tmp;
        resultDenom += pressureNew*pressureNew;

        if (FluidSystem::numActivePhases() > 1) {
            if (priVarsOld.primaryVarsMeaningWater() == PrimaryVariables::WaterMeaning::Sw) {
                saturationsOld[FluidSystem::waterPhaseIdx] =
                    priVarsOld[Indices::waterSwitchIdx];
                oilSaturationOld -= saturationsOld[FluidSystem::waterPhaseIdx];
            }

            if (priVarsOld.primaryVarsMeaningGas() == PrimaryVariables::GasMeaning::Sg)
            {
                assert(Indices::compositionSwitchIdx >= 0 );
                saturationsOld[FluidSystem::gasPhaseIdx] =
                    priVarsOld[Indices::compositionSwitchIdx];
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

    return resultDenom > 0.0 ? resultDelta / resultDenom : 0.0;
}

template <class TypeTag>
void
BlackoilModel<TypeTag>::
solveJacobianSystem(BVector& x)
{
    auto& jacobian = simulator_.model().linearizer().jacobian().istlMatrix();
    auto& residual = simulator_.model().linearizer().residual();
    auto& linSolver = simulator_.model().newtonMethod().linearSolver();

    const int numSolvers = linSolver.numAvailableSolvers();
    if (numSolvers > 1 && (linSolver.getSolveCount() % 100 == 0)) {
        if (terminal_output_) {
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
    }
    else {
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

template <class TypeTag>
void
BlackoilModel<TypeTag>::
updateSolution(const BVector& dx)
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
    }
}

template <class TypeTag>
std::tuple<typename BlackoilModel<TypeTag>::Scalar,
           typename BlackoilModel<TypeTag>::Scalar>
BlackoilModel<TypeTag>::
convergenceReduction(Parallel::Communication comm,
                     const Scalar pvSumLocal,
                     const Scalar numAquiferPvSumLocal,
                     std::vector< Scalar >& R_sum,
                     std::vector< Scalar >& maxCoeff,
                     std::vector< Scalar >& B_avg)
{
    OPM_TIMEBLOCK(convergenceReduction);
    // Compute total pore volume (use only owned entries)
    Scalar pvSum = pvSumLocal;
    Scalar numAquiferPvSum = numAquiferPvSumLocal;

    if (comm.size() > 1) {
        // global reduction
        std::vector< Scalar > sumBuffer;
        std::vector< Scalar > maxBuffer;
        const int numComp = B_avg.size();
        sumBuffer.reserve( 2*numComp + 2 ); // +2 for (numAquifer)pvSum
        maxBuffer.reserve( numComp );
        for (int compIdx = 0; compIdx < numComp; ++compIdx) {
            sumBuffer.push_back(B_avg[compIdx]);
            sumBuffer.push_back(R_sum[compIdx]);
            maxBuffer.push_back(maxCoeff[ compIdx]);
        }

        // Compute total pore volume
        sumBuffer.push_back( pvSum );
        sumBuffer.push_back( numAquiferPvSum );

        // compute global sum
        comm.sum( sumBuffer.data(), sumBuffer.size() );

        // compute global max
        comm.max( maxBuffer.data(), maxBuffer.size() );

        // restore values to local variables
        for (int compIdx = 0, buffIdx = 0; compIdx < numComp; ++compIdx, ++buffIdx) {
            B_avg[compIdx] = sumBuffer[buffIdx];
            ++buffIdx;

            R_sum[compIdx] = sumBuffer[buffIdx];
        }

        for (int compIdx = 0; compIdx < numComp; ++compIdx) {
            maxCoeff[compIdx] = maxBuffer[compIdx];
        }

        // restore global pore volume
        pvSum = sumBuffer[sumBuffer.size()-2];
        numAquiferPvSum = sumBuffer.back();
    }

    // return global pore volume
    return {pvSum, numAquiferPvSum};
}

template <class TypeTag>
std::pair<typename BlackoilModel<TypeTag>::Scalar,
          typename BlackoilModel<TypeTag>::Scalar>
BlackoilModel<TypeTag>::
localConvergenceData(std::vector<Scalar>& R_sum,
                     std::vector<Scalar>& maxCoeff,
                     std::vector<Scalar>& B_avg,
                     std::vector<int>& maxCoeffCell)
{
    OPM_TIMEBLOCK(localConvergenceData);
    Scalar pvSumLocal = 0.0;
    Scalar numAquiferPvSumLocal = 0.0;
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

        if (isNumericalAquiferCell(elem)) {
            numAquiferPvSumLocal += pvValue;
        }

        this->getMaxCoeff(cell_idx, intQuants, fs, residual, pvValue,
                          B_avg, R_sum, maxCoeff, maxCoeffCell);
    }

    OPM_END_PARALLEL_TRY_CATCH("BlackoilModel::localConvergenceData() failed: ", grid_.comm());

    // compute local average in terms of global number of elements
    const int bSize = B_avg.size();
    for (int i = 0; i < bSize; ++i) {
        B_avg[i] /= Scalar(global_nc_);
    }

    return {pvSumLocal, numAquiferPvSumLocal};
}

template <class TypeTag>
std::pair<std::vector<double>, std::vector<int>>
BlackoilModel<TypeTag>::
characteriseCnvPvSplit(const std::vector<Scalar>& B_avg, const double dt)
{
    OPM_TIMEBLOCK(computeCnvErrorPv);

    // 0: cnv <= tolerance_cnv
    // 1: tolerance_cnv < cnv <= tolerance_cnv_relaxed
    // 2: tolerance_cnv_relaxed < cnv
    constexpr auto numPvGroups = std::vector<double>::size_type{3};

    auto cnvPvSplit = std::pair<std::vector<double>, std::vector<int>> {
        std::piecewise_construct,
        std::forward_as_tuple(numPvGroups),
        std::forward_as_tuple(numPvGroups)
    };

    auto maxCNV = [&B_avg, dt](const auto& residual, const double pvol)
    {
        return (dt / pvol) *
            std::inner_product(residual.begin(), residual.end(),
                               B_avg.begin(), Scalar{0},
                               [](const Scalar m, const auto& x)
                               {
                                   using std::abs;
                                   return std::max(m, abs(x));
                               }, std::multiplies<>{});
    };

    auto& [splitPV, cellCntPV] = cnvPvSplit;

    const auto& model = this->simulator().model();
    const auto& problem = this->simulator().problem();
    const auto& residual = model.linearizer().residual();
    const auto& gridView = this->simulator().gridView();

    const IsNumericalAquiferCell isNumericalAquiferCell(gridView.grid());

    ElementContext elemCtx(this->simulator());

    OPM_BEGIN_PARALLEL_TRY_CATCH();
    for (const auto& elem : elements(gridView, Dune::Partitions::interior)) {
        // Skip cells of numerical Aquifer
        if (isNumericalAquiferCell(elem)) {
            continue;
        }

        elemCtx.updatePrimaryStencil(elem);

        const unsigned cell_idx = elemCtx.globalSpaceIndex(/*spaceIdx=*/0, /*timeIdx=*/0);
        const auto pvValue = problem.referencePorosity(cell_idx, /*timeIdx=*/0)
            * model.dofTotalVolume(cell_idx);

        const auto maxCnv = maxCNV(residual[cell_idx], pvValue);

        const auto ix = (maxCnv > this->param_.tolerance_cnv_)
            + (maxCnv > this->param_.tolerance_cnv_relaxed_);

        splitPV[ix] += static_cast<double>(pvValue);
        ++cellCntPV[ix];
    }

    OPM_END_PARALLEL_TRY_CATCH("BlackoilModel::characteriseCnvPvSplit() failed: ",
                               this->grid_.comm());

    this->grid_.comm().sum(splitPV  .data(), splitPV  .size());
    this->grid_.comm().sum(cellCntPV.data(), cellCntPV.size());

    return cnvPvSplit;
}

template <class TypeTag>
void
BlackoilModel<TypeTag>::
updateTUNING(const Tuning& tuning)
{
    this->param_.tolerance_cnv_ = tuning.TRGCNV;
    this->param_.tolerance_cnv_relaxed_ = tuning.XXXCNV;
    this->param_.tolerance_mb_ = tuning.TRGMBE;
    this->param_.tolerance_mb_relaxed_ = tuning.XXXMBE;

    if (terminal_output_) {
        const auto msg = fmt::format("Tuning BlackoilModel tolerances: "
                                     "MB: {:.2e}, CNV: {:.2e}",
                                     tuning.TRGMBE, tuning.TRGCNV);
        OpmLog::debug(msg);
    }
}

template <class TypeTag>
ConvergenceReport
BlackoilModel<TypeTag>::
getReservoirConvergence(const double reportTime,
                        const double dt,
                        const int iteration,
                        const int maxIter,
                        std::vector<Scalar>& B_avg,
                        std::vector<Scalar>& residual_norms)
{
    OPM_TIMEBLOCK(getReservoirConvergence);
    using Vector = std::vector<Scalar>;

    ConvergenceReport report{reportTime};

    const int numComp = numEq;

    Vector R_sum(numComp, Scalar{0});
    Vector maxCoeff(numComp, std::numeric_limits<Scalar>::lowest());
    std::vector<int> maxCoeffCell(numComp, -1);

    const auto [pvSumLocal, numAquiferPvSumLocal] =
        this->localConvergenceData(R_sum, maxCoeff, B_avg, maxCoeffCell);

    // compute global sum and max of quantities
    const auto& [pvSum, numAquiferPvSum] =
        this->convergenceReduction(this->grid_.comm(),
                                   pvSumLocal,
                                   numAquiferPvSumLocal,
                                   R_sum, maxCoeff, B_avg);

    report.setCnvPoreVolSplit(this->characteriseCnvPvSplit(B_avg, dt),
                              pvSum - numAquiferPvSum);

    // For each iteration, we need to determine whether to use the
    // relaxed tolerances.  To disable the usage of relaxed
    // tolerances, you can set the relaxed tolerances as the strict
    // tolerances.  If min_strict_mb_iter = -1 (default) we use a
    // relaxed tolerance for the mass balance for the last
    // iterations.  For positive values we use the relaxed tolerance
    // after the given number of iterations
    const bool relax_final_iteration_mb =
        this->param_.min_strict_mb_iter_ < 0 && iteration == maxIter;

    const bool use_relaxed_mb = relax_final_iteration_mb ||
        (this->param_.min_strict_mb_iter_ >= 0 &&
         iteration >= this->param_.min_strict_mb_iter_);

    // If min_strict_cnv_iter = -1 we use a relaxed tolerance for
    // the cnv for the last iterations.  For positive values we use
    // the relaxed tolerance after the given number of iterations.
    // We also use relaxed tolerances for cells with total
    // pore-volume less than relaxed_max_pv_fraction_.  Default
    // value of relaxed_max_pv_fraction_ is 0.03
    const bool relax_final_iteration_cnv =
        this->param_.min_strict_cnv_iter_ < 0 && iteration == maxIter;

    const bool relax_iter_cnv =
        this->param_.min_strict_cnv_iter_ >= 0 &&
        iteration >= this->param_.min_strict_cnv_iter_;

    // Note trailing parentheses here, just before the final
    // semicolon.  This is an immediately invoked function
    // expression which calculates a single boolean value.
    const auto relax_pv_fraction_cnv =
        [&report, this, eligible = pvSum - numAquiferPvSum]()
    {
        const auto& cnvPvSplit = report.cnvPvSplit().first;

        // [1]: tol < cnv <= relaxed
        // [2]: relaxed < cnv
        return static_cast<Scalar>(cnvPvSplit[1] + cnvPvSplit[2]) <
            this->param_.relaxed_max_pv_fraction_ * eligible;
    }();

    const bool use_relaxed_cnv = relax_final_iteration_cnv
        || relax_pv_fraction_cnv
        || relax_iter_cnv;

    if ((relax_final_iteration_mb || relax_final_iteration_cnv) &&
        this->terminal_output_)
    {
        std::string message =
            "Number of newton iterations reached its maximum "
            "try to continue with relaxed tolerances:";

        if (relax_final_iteration_mb) {
            message += fmt::format("  MB: {:.1e}", param_.tolerance_mb_relaxed_);
        }

        if (relax_final_iteration_cnv) {
            message += fmt::format(" CNV: {:.1e}", param_.tolerance_cnv_relaxed_);
        }

        OpmLog::debug(message);
    }

    const auto tol_cnv = use_relaxed_cnv ? param_.tolerance_cnv_relaxed_ : param_.tolerance_cnv_;
    const auto tol_mb  = use_relaxed_mb ? param_.tolerance_mb_relaxed_ : param_.tolerance_mb_;
    const auto tol_cnv_energy = use_relaxed_cnv ? param_.tolerance_cnv_energy_relaxed_ : param_.tolerance_cnv_energy_;
    const auto tol_eb = use_relaxed_mb ? param_.tolerance_energy_balance_relaxed_ : param_.tolerance_energy_balance_;

    // Finish computation
    std::vector<Scalar> CNV(numComp);
    std::vector<Scalar> mass_balance_residual(numComp);
    for (int compIdx = 0; compIdx < numComp; ++compIdx)
    {
        CNV[compIdx]                    = B_avg[compIdx] * dt * maxCoeff[compIdx];
        mass_balance_residual[compIdx]  = std::abs(B_avg[compIdx]*R_sum[compIdx]) * dt / pvSum;
        residual_norms.push_back(CNV[compIdx]);
    }

    using CR = ConvergenceReport;
    for (int compIdx = 0; compIdx < numComp; ++compIdx) {
        const Scalar res[2] = {
            mass_balance_residual[compIdx], CNV[compIdx],
        };

        const CR::ReservoirFailure::Type types[2] = {
            CR::ReservoirFailure::Type::MassBalance,
            CR::ReservoirFailure::Type::Cnv,
        };

        Scalar tol[2] = { tol_mb, tol_cnv, };
        if (has_energy_ && compIdx == contiEnergyEqIdx) {
            tol[0] = tol_eb;
            tol[1] = tol_cnv_energy;
        }

        for (int ii : {0, 1}) {
            if (std::isnan(res[ii])) {
                report.setReservoirFailed({types[ii], CR::Severity::NotANumber, compIdx});
                if (this->terminal_output_) {
                    OpmLog::debug("NaN residual for " + this->compNames_.name(compIdx) + " equation.");
                }
            }
            else if (res[ii] > maxResidualAllowed()) {
                report.setReservoirFailed({types[ii], CR::Severity::TooLarge, compIdx});
                if (this->terminal_output_) {
                    OpmLog::debug("Too large residual for " + this->compNames_.name(compIdx) + " equation.");
                }
            }
            else if (res[ii] < 0.0) {
                report.setReservoirFailed({types[ii], CR::Severity::Normal, compIdx});
                if (this->terminal_output_) {
                    OpmLog::debug("Negative residual for " + this->compNames_.name(compIdx) + " equation.");
                }
            }
            else if (res[ii] > tol[ii]) {
                report.setReservoirFailed({types[ii], CR::Severity::Normal, compIdx});
            }

            report.setReservoirConvergenceMetric(types[ii], compIdx, res[ii], tol[ii]);
        }
    }

    // Compute the Newton convergence per cell.
    this->convergencePerCell(B_avg, dt, tol_cnv, tol_cnv_energy, iteration);

    // Output of residuals.
    if (this->terminal_output_) {
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

template <class TypeTag>
void
BlackoilModel<TypeTag>::
convergencePerCell(const std::vector<Scalar>& B_avg,
                   const double dt,
                   const double tol_cnv,
                   const double tol_cnv_energy,
                   const int iteration)
{
    auto& rst_conv = simulator_.problem().eclWriter().mutableOutputModule().getConv();
    if (!rst_conv.hasConv()) {
        return;
    }

    if (iteration == 0) {
        rst_conv.prepareConv();
    }

    const auto& residual = simulator_.model().linearizer().residual();
    const auto& gridView = this->simulator().gridView();
    const IsNumericalAquiferCell isNumericalAquiferCell(gridView.grid());
    ElementContext elemCtx(this->simulator());
    std::vector<int> convNewt(residual.size(), 0);
    OPM_BEGIN_PARALLEL_TRY_CATCH();
    unsigned idx = 0;
    const int numComp = B_avg.size();
    for (const auto& elem : elements(gridView, Dune::Partitions::interior)) {
        elemCtx.updatePrimaryStencil(elem);

        const unsigned cell_idx = elemCtx.globalSpaceIndex(/*spaceIdx=*/0, /*timeIdx=*/0);
        const auto pvValue = simulator_.problem().referencePorosity(cell_idx, /*timeIdx=*/0) *
                             simulator_.model().dofTotalVolume(cell_idx);
        for (int compIdx = 0; compIdx < numComp; ++compIdx) {
            const auto tol = (has_energy_ && compIdx == contiEnergyEqIdx) ? tol_cnv_energy : tol_cnv;
            const Scalar cnv = std::abs(B_avg[compIdx] * residual[cell_idx][compIdx]) * dt / pvValue;
            if (std::isnan(cnv) || cnv > maxResidualAllowed() || cnv < 0.0 || cnv > tol) {
                convNewt[idx] = 1;
                break;
            }
        }
        ++idx;
    }
    OPM_END_PARALLEL_TRY_CATCH("BlackoilModel::convergencePerCell() failed: ",
                               this->grid_.comm());
    rst_conv.updateNewton(convNewt);
}

template <class TypeTag>
ConvergenceReport
BlackoilModel<TypeTag>::
getConvergence(const SimulatorTimerInterface& timer,
               const int iteration,
               const int maxIter,
               std::vector<Scalar>& residual_norms)
{
    OPM_TIMEBLOCK(getConvergence);
    // Get convergence reports for reservoir and wells.
    std::vector<Scalar> B_avg(numEq, 0.0);
    auto report = getReservoirConvergence(timer.simulationTimeElapsed(),
                                          timer.currentStepLength(),
                                          iteration, maxIter, B_avg, residual_norms);
    {
        OPM_TIMEBLOCK(getWellConvergence);
        report += wellModel().getWellConvergence(B_avg,
                                                 /*checkWellGroupControls*/report.converged());
    }

    conv_monitor_.checkPenaltyCard(report, iteration);

    return report;
}

template <class TypeTag>
std::vector<std::vector<typename BlackoilModel<TypeTag>::Scalar> >
BlackoilModel<TypeTag>::
computeFluidInPlace(const std::vector<int>& /*fipnum*/) const
{
    OPM_TIMEBLOCK(computeFluidInPlace);
    // assert(true)
    // return an empty vector
    std::vector<std::vector<Scalar> > regionValues(0, std::vector<Scalar>(0,0.0));
    return regionValues;
}

template <class TypeTag>
const SimulatorReport&
BlackoilModel<TypeTag>::
localAccumulatedReports() const
{
    if (!hasNlddSolver()) {
        OPM_THROW(std::runtime_error, "Cannot get local reports from a model without NLDD solver");
    }
    return nlddSolver_->localAccumulatedReports();
}

template <class TypeTag>
const std::vector<SimulatorReport>&
BlackoilModel<TypeTag>::
domainAccumulatedReports() const
{
    if (!nlddSolver_)
        OPM_THROW(std::runtime_error, "Cannot get domain reports from a model without NLDD solver");
    return nlddSolver_->domainAccumulatedReports();
}

template <class TypeTag>
void
BlackoilModel<TypeTag>::
writeNonlinearIterationsPerCell(const std::filesystem::path& odir) const
{
    if (hasNlddSolver()) {
        nlddSolver_->writeNonlinearIterationsPerCell(odir);
    }
}

template <class TypeTag>
void
BlackoilModel<TypeTag>::
writePartitions(const std::filesystem::path& odir) const
{
    if (hasNlddSolver()) {
        nlddSolver_->writePartitions(odir);
        return;
    }

    const auto& elementMapper = this->simulator().model().elementMapper();
    const auto& cartMapper = this->simulator().vanguard().cartesianIndexMapper();

    const auto& grid = this->simulator().vanguard().grid();
    const auto& comm = grid.comm();
    const auto nDigit = 1 + static_cast<int>(std::floor(std::log10(comm.size())));

    std::ofstream pfile {odir / fmt::format("{1:0>{0}}", nDigit, comm.rank())};

    for (const auto& cell : elements(grid.leafGridView(), Dune::Partitions::interior)) {
        pfile << comm.rank() << ' '
              << cartMapper.cartesianIndex(elementMapper.index(cell)) << ' '
              << comm.rank() << '\n';
    }
}

template <class TypeTag>
template<class FluidState, class Residual>
void
BlackoilModel<TypeTag>::
getMaxCoeff(const unsigned cell_idx,
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

        const unsigned sIdx = FluidSystem::solventComponentIndex(phaseIdx);
        const unsigned compIdx = Indices::canonicalToActiveComponentIndex(sIdx);

        B_avg[compIdx] += 1.0 / fs.invB(phaseIdx).value();
        const auto R2 = modelResid[cell_idx][compIdx];

        R_sum[compIdx] += R2;
        const Scalar Rval = std::abs(R2) / pvValue;
        if (Rval > maxCoeff[compIdx]) {
            maxCoeff[compIdx] = Rval;
            maxCoeffCell[compIdx] = cell_idx;
        }
    }

    if constexpr (has_solvent_) {
        B_avg[contiSolventEqIdx] +=
            1.0 / intQuants.solventInverseFormationVolumeFactor().value();
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

} // namespace Opm

#endif // OPM_BLACKOILMODEL_IMPL_HEADER_INCLUDED
