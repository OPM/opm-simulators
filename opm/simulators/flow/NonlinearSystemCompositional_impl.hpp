/*
  Copyright 2026, SINTEF Digital

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

#ifndef OPM_NONLINEAR_SYSTEM_COMPOSITIONAL_IMPL_HEADER_INCLUDED
#define OPM_NONLINEAR_SYSTEM_COMPOSITIONAL_IMPL_HEADER_INCLUDED

#ifndef OPM_NONLINEAR_SYSTEM_COMPOSITIONAL_HEADER_INCLUDED
#include <config.h>
#include <opm/simulators/flow/NonlinearSystemCompositional.hpp>
#endif

#include <dune/common/timer.hh>

#include <opm/common/ErrorMacros.hpp>
#include <opm/common/OpmLog/OpmLog.hpp>

#include <opm/material/common/MathToolbox.hpp>

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <iomanip>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <utility>

namespace Opm {

template <class TypeTag>
NonlinearSystemCompositional<TypeTag>::
NonlinearSystemCompositional(Simulator& simulator,
                             const ModelParameters& param,
                             CompWellModel<TypeTag>& wellModel,
                             const bool terminalOutput)
    : ParentType(simulator, param, wellModel, terminalOutput)
{
    this->convergence_reports_.reserve(64);
}

template <class TypeTag>
SimulatorReportSingle
NonlinearSystemCompositional<TypeTag>::
prepareStep(const SimulatorTimerInterface& timer)
{
    SimulatorReportSingle report;
    Dune::Timer perfTimer;
    perfTimer.start();

    const int lastStepFailed = timer.lastStepFailed();
    if (this->grid_.comm().size() > 1
        && this->grid_.comm().max(lastStepFailed) != this->grid_.comm().min(lastStepFailed)) {
        OPM_THROW(std::runtime_error,
                  "Misalignment of the parallel simulation run in prepareStep "
                  "- the previous step succeeded on some ranks but failed on others.");
    }

    if (lastStepFailed) {
        this->wellModel().restoreLastValidState();
        this->simulator_.model().updateFailed();
    }
    else {
        this->simulator_.model().advanceTimeLevel();
    }

    this->simulator_.setTime(timer.simulationTimeElapsed());
    this->simulator_.setTimeStepSize(timer.currentStepLength());

    this->simulator_.problem().resetIterationForNewTimestep();
    this->simulator_.problem().beginTimeStep();

    report.pre_post_time += perfTimer.stop();
    return report;
}

template <class TypeTag>
void
NonlinearSystemCompositional<TypeTag>::
initialLinearization(SimulatorReportSingle& report,
                     const int minIter,
                     const int maxIter,
                     const SimulatorTimerInterface& timer)
{
    ParentType::initialLinearization(report,
                                     minIter,
                                     maxIter,
                                     timer);
    Dune::Timer perfTimer;
    perfTimer.start();

    // Calculate reservoir and well convergence and store in convergence history
    std::vector<Scalar> residual_norms;
    auto convrep = getConvergence(timer, residual_norms);

    // Report converged flag
    report.converged = convrep.converged()
        && this->simulator_.problem().iterationContext().iteration() >= minIter;

    // Throw for severe failures
    const auto severity = convrep.severityOfWorstFailure();
    this->convergence_reports_.back().report.push_back(std::move(convrep));
    if (severity == ConvergenceReport::Severity::NotANumber) {
        this->failureReport_ += report;
        OPM_THROW_PROBLEM(NumericalProblem, "NaN convergence values found!");
    }

    if (severity == ConvergenceReport::Severity::TooLarge) {
        this->failureReport_ += report;
        OPM_THROW_NOLOG(NumericalProblem, "Too large convergence values found!");
    }

    report.update_time += perfTimer.stop();

    // Store residual norms in history container
    this->residual_norms_history_.push_back(residual_norms);
}

template <class TypeTag>
template <class NonlinearSolverType>
SimulatorReportSingle
NonlinearSystemCompositional<TypeTag>::
nonlinearIteration(const SimulatorTimerInterface& timer,
                   NonlinearSolverType& nonlinearSolver)
{
    if (this->simulator_.problem().iterationContext().needsTimestepInit()) {
        this->residual_norms_history_.clear();
        this->current_relaxation_ = 1.0;
        this->dx_old_ = 0.0;
        this->convergence_reports_.push_back({timer.reportStepNum(), timer.currentStepNum(), {}});
        this->convergence_reports_.back().report.reserve(numEq);
    }

    auto result = this->nonlinearIterationNewton(timer, nonlinearSolver);
    this->simulator_.problem().advanceIteration();
    return result;
}

template <class TypeTag>
template <class NonlinearSolverType>
SimulatorReportSingle
NonlinearSystemCompositional<TypeTag>::
nonlinearIterationNewton(const SimulatorTimerInterface& timer,
                         NonlinearSolverType& nonlinearSolver)
{
    OPM_TIMEFUNCTION();

    SimulatorReportSingle report;
    Dune::Timer perfTimer;

    this->initialLinearization(report,
                               this->param_.newton_min_iter_,
                               this->param_.newton_max_iter_,
                               timer);

    if (!report.converged) {
        perfTimer.reset();
        perfTimer.start();
        report.total_newton_iterations = 1;

        BVector x(this->simulator_.model().numGridDof());
        this->linear_solve_setup_time_ = 0.0;

        try {
            auto& linearizer = this->simulator_.model().linearizer();
            linearizer.linearizeAuxiliaryEquations();
            linearizer.finalize();

            this->solveJacobianSystem(x);

            report.linear_solve_setup_time += this->linear_solve_setup_time_;
            report.linear_solve_time += perfTimer.stop();
            report.total_linear_iterations += this->linearIterationsLastSolve();
        }
        catch (...) {
            report.linear_solve_setup_time += this->linear_solve_setup_time_;
            report.linear_solve_time += perfTimer.stop();
            report.total_linear_iterations += this->linearIterationsLastSolve();

            this->failureReport_ += report;
            throw;
        }

        perfTimer.reset();
        perfTimer.start();

        auto& model = this->simulator_.model();
        for (unsigned auxModIdx = 0; auxModIdx < model.numAuxiliaryModules(); ++auxModIdx) {
            model.auxiliaryModule(auxModIdx)->postSolve(x);
        }

        if (this->param_.use_update_stabilization_) {
            bool isOscillate = false;
            bool isStagnate = false;
            nonlinearSolver.detectOscillations(this->residual_norms_history_,
                                               this->residual_norms_history_.size() - 1,
                                               isOscillate,
                                               isStagnate);

            if (isOscillate) {
                this->current_relaxation_ -= nonlinearSolver.relaxIncrement();
                this->current_relaxation_ = std::max(this->current_relaxation_, nonlinearSolver.relaxMax());

                if (this->terminalOutputEnabled()) {
                    OpmLog::info("    Oscillating behavior detected: Relaxation set to "
                                 + std::to_string(this->current_relaxation_));
                }
            }

            nonlinearSolver.stabilizeNonlinearUpdate(x, this->dx_old_, this->current_relaxation_);
        }

        this->updateSolution(x);
        report.update_time += perfTimer.stop();
    }

    return report;
}

template <class TypeTag>
void
NonlinearSystemCompositional<TypeTag>::
updateTUNING(const Tuning& /*tuning*/)
{
    OpmLog::warning("TUNING convergence records 2 and 3 are ignored compositional simulations!");
}

template <class TypeTag>
void
NonlinearSystemCompositional<TypeTag>::
updateTUNINGDP(const TuningDp& /*tuning_dp*/)
{
    OpmLog::warning("TUNINGDP is ignored in compositional simulations!");
}

template <class TypeTag>
typename NonlinearSystemCompositional<TypeTag>::Scalar
NonlinearSystemCompositional<TypeTag>::
relativeChange() const
{
    Scalar resultDelta = 0.0;
    Scalar resultDenom = 0.0;

    const auto& elemMapper = this->simulator_.model().elementMapper();
    const auto& gridView = this->simulator_.gridView();

    for (const auto& elem : elements(gridView, Dune::Partitions::interior)) {
        const unsigned globalElemIdx = elemMapper.index(elem);
        const auto& priVarsNew = this->simulator_.model().solution(/*timeIdx=*/0)[globalElemIdx];
        const auto& priVarsOld = this->simulator_.model().solution(/*timeIdx=*/1)[globalElemIdx];

        for (int pvIdx = 0; pvIdx < static_cast<int>(priVarsNew.size()); ++pvIdx) {
            const auto delta = priVarsNew[pvIdx] - priVarsOld[pvIdx];
            resultDelta += delta * delta;
            resultDenom += priVarsNew[pvIdx] * priVarsNew[pvIdx];
        }
    }

    resultDelta = gridView.comm().sum(resultDelta);
    resultDenom = gridView.comm().sum(resultDenom);

    return resultDenom > 0.0 ? resultDelta / resultDenom : 0.0;
}

template <class TypeTag>
void
NonlinearSystemCompositional<TypeTag>::
solveJacobianSystem(BVector& x)
{
    auto& jacobian = this->simulator_.model().linearizer().jacobian();
    auto& residual = this->simulator_.model().linearizer().residual();
    auto& linSolver = this->simulator_.model().newtonMethod().linearSolver();

    x = 0.0;

    Dune::Timer perfTimer;
    perfTimer.start();
    linSolver.prepare(jacobian, residual);
    this->linear_solve_setup_time_ = perfTimer.stop();
    linSolver.setResidual(residual);
    linSolver.getResidual(residual);
    linSolver.setMatrix(jacobian);
    linSolver.solve(x);
}

template <class TypeTag>
void
NonlinearSystemCompositional<TypeTag>::
prepareSolutionUpdate()
{
    // Init. solution update vector
    unsigned nc = this->simulator_.model().numGridDof();
    dP_.resize(nc);
    dSeff_.resize(nc);
    dP_ = 0.0;
    dSeff_ = 0.0;
    effSatData_.resize(nc);

    const auto& elemMapper = this->simulator_.model().elementMapper();
    const auto& gridView = this->simulator_.gridView();
    for (const auto& elem : elements(gridView, Dune::Partitions::interior)) {
        // Compute effective saturation before Newton iteration
        unsigned globalElemIdx = elemMapper.index(elem);
        // TODO: use element context?
        const auto* intQuants =
            this->simulator_.model().cachedIntensiveQuantities(globalElemIdx, /*timeIdx=*/0);
        assert(intQuants);
        const auto& fs = intQuants->fluidState();
        effSatData_[globalElemIdx] = computeEffectiveSaturationData_(fs);
    }
}

template <class TypeTag>
void
NonlinearSystemCompositional<TypeTag>::
storeSolutionUpdate(const GlobalEqVector& dx)
{
    const auto& elemMapper = this->simulator_.model().elementMapper();
    const auto& gridView = this->simulator_.gridView();
    for (const auto& elem : elements(gridView, Dune::Partitions::interior)) {
         unsigned globalElemIdx = elemMapper.index(elem);

        // Store pressure update
        const auto& dP = dx[globalElemIdx][Indices::pressure0Idx];
        dP_[globalElemIdx] = dP;

        // Calculate effective saturation after Newton iteration to use in diff.
        // TODO: use element context?
        const auto* intQuants =
            this->simulator_.model().cachedIntensiveQuantities(globalElemIdx, /*timeIdx=*/0);
        assert(intQuants);
        const auto& fs = intQuants->fluidState();
        dSeff_[globalElemIdx] = effectiveSaturationChange_(effSatData_[globalElemIdx],
                                                           computeEffectiveSaturationData_(fs));
    }
}

template <class TypeTag>
ConvergenceReport
NonlinearSystemCompositional<TypeTag>::
getConvergence(const SimulatorTimerInterface& timer,
               std::vector<Scalar>& residual_norms)
{
    // Reservoir compositional convergence report
    auto report = getCompositionalConvergence(timer.simulationTimeElapsed(), residual_norms);

    // Well convergence report
    ConvergenceReport wellReport(timer.simulationTimeElapsed());
    const bool wellConverged = this->wellModel().getWellConvergence();
    using CR = ConvergenceReport;
    if (!wellConverged) {
        // Random failure here since CompWellModel does not return ConvergenceReport
        // TODO: change this when wells return ConvergenceReport
        wellReport.setWellFailed(
            {CR::WellFailure::Type::Unsolvable, CR::Severity::Normal, -1, "Unknown"});
    }

    report += wellReport;
    return report;
}

template <class TypeTag>
ConvergenceReport
NonlinearSystemCompositional<TypeTag>::
getCompositionalConvergence(double reportTime,
                            std::vector<Scalar>& residual_norms)
{
    // Init. nonlinear iteration convergence report
    ConvergenceReport report{reportTime};

    using CR = ConvergenceReport;
    using FailureType = CR::ReservoirFailure::Type;
    const std::array types = {FailureType::MaxDP, FailureType::MaxDSeff};

    // No solution update exists yet in the first iteration of a timestep
    const auto& iterCtx = this->simulator_.problem().iterationContext();
    const bool hasSolutionUpdate = !iterCtx.isFirstGlobalIteration();

    // Local and global convergence data
    Scalar dPmax = 0.0;
    Scalar dSmax = 0.0;
    std::vector<Scalar> residualMaxNorm;
    std::vector<Scalar> residualSum;
    std::vector<Scalar> specificVolumeAvg;
    const Scalar poreVolumeSumLocal = localCompositionalConvergenceData(
        dPmax, dSmax, residualMaxNorm, residualSum, specificVolumeAvg);
    const Scalar poreVolumeSum = compositionalConvergenceReduction(
        poreVolumeSumLocal, dPmax, dSmax, residualMaxNorm, residualSum, specificVolumeAvg);

    // Converged if the solution change criteria are met, or, regardless of them, if the residual
    // max-norm or sum is below its tolerance
    const Scalar tolDp = this->param_.tolerance_max_dp_;
    const Scalar tolDs = this->param_.tolerance_max_ds_;
    const bool solutionChangeConverged = hasSolutionUpdate && dPmax <= tolDp && dSmax <= tolDs;

    // CNV and MB residual metrics: max. pore-volume normalized residual and residual sum
    const Scalar dt = this->simulator_.timeStepSize();
    std::vector<Scalar> cnv(numEq);
    std::vector<Scalar> mb(numEq);
    Scalar resMax = 0.0;
    Scalar resSum = 0.0;
    for (int eqIdx = 0; eqIdx < numEq; ++eqIdx) {
        const Scalar specVol = specificVolumeAvg[eqIdx] > 0.0 ? specificVolumeAvg[eqIdx] : 1.0;
        cnv[eqIdx] = specVol * dt * residualMaxNorm[eqIdx];
        mb[eqIdx] = specVol * dt * std::abs(residualSum[eqIdx]) / poreVolumeSum;
        resMax = std::max(resMax, cnv[eqIdx]);
        resSum = std::max(resSum, mb[eqIdx]);
    }
    const bool residualConverged
        = resMax < this->param_.tolerance_cnv_ || resSum < this->param_.tolerance_mb_;
    const bool converged = solutionChangeConverged || residualConverged;

    // Record the same metrics in every iteration
    const auto logFailure = [this](const std::string& message) {
        if (this->terminal_output_) {
            OpmLog::debug(message);
        }
    };

    // Residual max-norm as CNV and sum as MB for each equation
    const Scalar noTolerance = std::numeric_limits<Scalar>::max();
    const std::array residualTypes = {FailureType::Cnv, FailureType::MassBalance};
    const std::array<Scalar, 2> residualTolerances = converged
        ? std::array<Scalar, 2>{noTolerance, noTolerance}
        : std::array<Scalar, 2>{this->param_.tolerance_cnv_, this->param_.tolerance_mb_};
    for (int eqIdx = 0; eqIdx < numEq; ++eqIdx) {
        const std::array<Scalar, 2> residuals = {cnv[eqIdx], mb[eqIdx]};
        this->addReservoirConvergenceMetrics(report,
                                             eqIdx,
                                             this->compNames_.name(eqIdx),
                                             residuals,
                                             residualTypes,
                                             residualTolerances,
                                             this->param_.max_residual_allowed_,
                                             logFailure);
    }

    // Store max-residual vectors for history
    residual_norms = std::move(cnv);

    // Max. pressure and effective saturation change
    const std::array<Scalar, 2> dSolTolerances = converged
        ? std::array<Scalar, 2>{noTolerance, noTolerance}
        : std::array<Scalar, 2>{tolDp, tolDs};
    if (hasSolutionUpdate) {
        const std::array<Scalar, 2> dSolmax = {dPmax, dSmax};
        const std::array<std::string, 2> dSolnames = {"DPMAX", "DSEFFMAX"};
        Scalar maxDSolAllowed = 1.0e20;
        addCompositionalConvergenceMetrics(report,
                                           dSolmax,
                                           dSolnames,
                                           types,
                                           dSolTolerances,
                                           maxDSolAllowed,
                                           logFailure);
    }
    else {
        // No update yet: record a placeholder, and fail unless the residual criteria are met
        for (std::size_t typeIdx = 0; typeIdx < types.size(); ++typeIdx) {
            if (!converged) {
                report.setReservoirFailed({types[typeIdx], CR::Severity::Normal, -1});
            }
            report.setReservoirConvergenceMetric(types[typeIdx],
                                                 -1,
                                                 std::numeric_limits<Scalar>::quiet_NaN(),
                                                 dSolTolerances[typeIdx]);
        }
    }

    // Output convergence
    if (this->terminal_output_) {
        // Header
        if (iterCtx.isFirstGlobalIteration()) {
            std::string msg = "Iter    DPMAX      DSMAX        CNV         MB  ";
            OpmLog::debug(msg);
        }

        // Print values
        std::ostringstream ss;
        const std::streamsize oprec = ss.precision(3);
        const std::ios::fmtflags oflags = ss.setf(std::ios::scientific);

        ss << std::setw(4) << iterCtx.iteration();
        if (hasSolutionUpdate) {
            ss << std::setw(11) << dPmax;
            ss << std::setw(11) << dSmax;
        }
        else {
            ss << std::setw(11) << "-";
            ss << std::setw(11) << "-";
        }
        ss << std::setw(11) << resMax;
        ss << std::setw(11) << resSum;

        ss.precision(oprec);
        ss.flags(oflags);

        OpmLog::debug(ss.str());
    }

    return report;
}

template <class TypeTag>
typename NonlinearSystemCompositional<TypeTag>::Scalar
NonlinearSystemCompositional<TypeTag>::
localCompositionalConvergenceData(Scalar& dPmax,
                                  Scalar& dSmax,
                                  std::vector<Scalar>& residualMaxNorm,
                                  std::vector<Scalar>& residualSum,
                                  std::vector<Scalar>& specificVolumeAvg) const
{
    // Max. absolute pressure and effective saturation change over (local) cells
    dPmax = dP_.infinity_norm();
    dSmax = dSeff_.infinity_norm();

    // Max. pore-volume normalized residual, residual sum and pore volume sum over the interior
    // cells
    const auto& model = this->simulator_.model();
    const auto& problem = this->simulator_.problem();
    const auto& residual = model.linearizer().residual();
    const auto& constraintsMap = model.linearizer().constraintsMap();

    residualMaxNorm.assign(numEq, 0.0);
    residualSum.assign(numEq, 0.0);
    specificVolumeAvg.assign(numEq, 0.0);
    Scalar poreVolumeSumLocal = 0.0;

    const auto& elemMapper = model.elementMapper();
    for (const auto& elem : elements(this->simulator_.gridView(), Dune::Partitions::interior)) {
        const unsigned dofIdx = elemMapper.index(elem);
        if (dofIdx >= model.numGridDof() || model.dofTotalVolume(dofIdx) <= 0.0) {
            continue;
        }

        if (constraintsMap.count(dofIdx) > 0) {
            continue;
        }

        const Scalar pvValue
            = problem.referencePorosity(dofIdx, /*timeIdx=*/0) * model.dofTotalVolume(dofIdx);
        if (pvValue <= 0.0) {
            continue;
        }
        poreVolumeSumLocal += pvValue;

        // Specific volume sums (i.e., sum over the inverse of component densities)
        const auto* intQuants = model.cachedIntensiveQuantities(dofIdx, /*timeIdx=*/0);
        assert(intQuants);
        const auto& fs = intQuants->fluidState();
        const Scalar L = decay<Scalar>(fs.L());
        const Scalar massPerMole = L * decay<Scalar>(fs.averageMolarMass(FluidSystem::oilPhaseIdx))
            + (1.0 - L) * decay<Scalar>(fs.averageMolarMass(FluidSystem::gasPhaseIdx));
        const Scalar volumePerMole = L / decay<Scalar>(fs.molarDensity(FluidSystem::oilPhaseIdx))
            + (1.0 - L) / decay<Scalar>(fs.molarDensity(FluidSystem::gasPhaseIdx));
        const Scalar hcSpecificVolume = volumePerMole / massPerMole;
        for (int compIdx = 0; compIdx < numComponents; ++compIdx) {
            specificVolumeAvg[Indices::conti0EqIdx + compIdx] += hcSpecificVolume;
        }
        if constexpr (waterEnabled) {
            specificVolumeAvg[Indices::conti0EqIdx + numComponents]
                += 1.0 / decay<Scalar>(fs.density(FluidSystem::waterPhaseIdx));
        }

        // Max. and sum residuals, accounting for use of volumetric residuals
        const Scalar residualScale = useVolumetricResidual ? model.dofTotalVolume(dofIdx) : 1.0;
        const auto& localResidual = residual[dofIdx];
        for (int eqIdx = 0; eqIdx < numEq; ++eqIdx) {
            const Scalar r = localResidual[eqIdx] * residualScale;
            residualMaxNorm[eqIdx] = std::max(residualMaxNorm[eqIdx], std::abs(r) / pvValue);
            residualSum[eqIdx] += r;
        }
    }

    // Local share of the average over all cells of the global grid
    for (auto& value : specificVolumeAvg) {
        value /= Scalar(this->global_nc_);
    }

    return poreVolumeSumLocal;
}

template <class TypeTag>
typename NonlinearSystemCompositional<TypeTag>::Scalar
NonlinearSystemCompositional<TypeTag>::
compositionalConvergenceReduction(const Scalar poreVolumeSumLocal,
                                  Scalar& dPmax,
                                  Scalar& dSmax,
                                  std::vector<Scalar>& residualMaxNorm,
                                  std::vector<Scalar>& residualSum,
                                  std::vector<Scalar>& specificVolumeAvg)
{
    // Max. pressure and effective saturation change
    const auto& comm = this->grid_.comm();
    if (comm.size() > 1) {
        std::array<Scalar, 2> maxBuffer = {dPmax, dSmax};
        comm.max(maxBuffer.data(), maxBuffer.size());
        dPmax = maxBuffer[0];
        dSmax = maxBuffer[1];
    }

    // Residual metrics and pore volume. The second volume of the reduction is not used.
    const auto [poreVolumeSum, unused] = this->convergenceReduction(comm,
                                                                    poreVolumeSumLocal,
                                                                    Scalar{0.0},
                                                                    residualSum,
                                                                    residualMaxNorm,
                                                                    specificVolumeAvg);
    return poreVolumeSum;
}

template <class TypeTag>
template <class LogFailure>
void
NonlinearSystemCompositional<TypeTag>::
addCompositionalConvergenceMetrics(
    ConvergenceReport& report,
    const std::span<const Scalar> dSolmax,
    const std::span<const std::string> dSolnames,
    const std::span<const ConvergenceReport::ReservoirFailure::Type> types,
    const std::span<const Scalar> tolerances,
    const Scalar maxdSolMaxAllowed,
    LogFailure&& logFailure) const
{
    if (dSolmax.size() != dSolnames.size() || dSolmax.size() != types.size()
        || dSolmax.size() != tolerances.size()) {
        OPM_THROW(std::logic_error, "Mismatched compositional convergence metric sizes.");
    }

    using CR = ConvergenceReport;
    for (std::size_t metricIdx = 0; metricIdx < dSolmax.size(); ++metricIdx) {
        const auto dsolmax = dSolmax[metricIdx];
        const auto dsolname = dSolnames[metricIdx];
        const auto type = types[metricIdx];
        const auto tolerance = tolerances[metricIdx];

        // Failures
        if (std::isnan(dsolmax)) {
            report.setReservoirFailed({type, CR::Severity::NotANumber, -1});
            logFailure("NaN value for " + dsolname + " .");
        }
        else if (dsolmax > maxdSolMaxAllowed) {
            report.setReservoirFailed({type, CR::Severity::TooLarge, -1});
            logFailure("Too large value for " + dsolname + " .");
        }
        else if (dsolmax < 0.0) {
            report.setReservoirFailed({type, CR::Severity::Normal, -1});
            logFailure("Negative value for " + dsolname + " .");
        }
        else if (dsolmax > tolerance) {
            report.setReservoirFailed({type, CR::Severity::Normal, -1});
        }

        report.setReservoirConvergenceMetric(type, -1, dsolmax, tolerance);
    }
}

template <class TypeTag>
template <class FluidState>
typename NonlinearSystemCompositional<TypeTag>::EffectiveSaturationData
NonlinearSystemCompositional<TypeTag>::
computeEffectiveSaturationData_(const FluidState& fs) const
{
    EffectiveSaturationData data;

    // Component moles per pore volume
    for (const int phaseIdx : {FluidSystem::oilPhaseIdx, FluidSystem::gasPhaseIdx}) {
        const Scalar Sb = decay<Scalar>(fs.saturation(phaseIdx) * fs.molarDensity(phaseIdx));
        for (int compIdx = 0; compIdx < numComponents; ++compIdx) {
            data.molarDens[compIdx] += Sb * decay<Scalar>(fs.moleFraction(phaseIdx, compIdx));
        }
    }

    // Mixture molar volume and its derivatives w.r.t. z at fixed pressure
    const auto& L = fs.L();
    const auto v = L / fs.molarDensity(FluidSystem::oilPhaseIdx)
        + (1.0 - L) / fs.molarDensity(FluidSystem::gasPhaseIdx);
    data.molarVolume = decay<Scalar>(v);
    for (int compIdx = 0; compIdx < numComponents - 1; ++compIdx) {
        data.z[compIdx] = decay<Scalar>(fs.moleFraction(compIdx));
        data.dMolarVolumeDz[compIdx] = v.derivative(Indices::z0Idx + compIdx);
    }

    // Water volume per pore volume is m_w / rho_w
    if constexpr (waterEnabled) {
        data.waterDensity = decay<Scalar>(fs.density(FluidSystem::waterPhaseIdx));
        data.waterMassDens =
            decay<Scalar>(fs.saturation(FluidSystem::waterPhaseIdx)) * data.waterDensity;
    }

    return data;
}

template <class TypeTag>
typename NonlinearSystemCompositional<TypeTag>::Scalar
NonlinearSystemCompositional<TypeTag>::
effectiveSaturationChange_(const EffectiveSaturationData& oldData,
                           const EffectiveSaturationData& newData)
{
    // Linearized fluid volume change per pore volume of each component at fixed pressure,
    // dS'_c = vbar_c * dm_c, with partial molar volume vbar_c = v + dv/dz_c - sum_j z_j * dv/dz_j
    Scalar zdvdz = 0.0;
    for (int compIdx = 0; compIdx < numComponents - 1; ++compIdx) {
        zdvdz += oldData.z[compIdx] * oldData.dMolarVolumeDz[compIdx];
    }

    Scalar dSeffMax = 0.0;
    for (int compIdx = 0; compIdx < numComponents; ++compIdx) {
        const Scalar dvdz = compIdx < numComponents - 1 ? oldData.dMolarVolumeDz[compIdx] : 0.0;
        const Scalar partialMolarVolume = oldData.molarVolume + dvdz - zdvdz;
        const Scalar dSeffComp
            = partialMolarVolume * (newData.molarDens[compIdx] - oldData.molarDens[compIdx]);
        dSeffMax = std::max(dSeffMax, std::abs(dSeffComp));
    }

    if constexpr (waterEnabled) {
        if (oldData.waterDensity > 0.0) {
            const Scalar dSeffWater
                = (newData.waterMassDens - oldData.waterMassDens) / oldData.waterDensity;
            dSeffMax = std::max(dSeffMax, std::abs(dSeffWater));
        }
    }

    return dSeffMax;
}

} // namespace Opm

#endif // OPM_NONLINEAR_SYSTEM_COMPOSITIONAL_IMPL_HEADER_INCLUDED
