#ifndef OPM_FLOW_NONLINEAR_SYSTEM_IMPL_HEADER_INCLUDED
#define OPM_FLOW_NONLINEAR_SYSTEM_IMPL_HEADER_INCLUDED
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
#ifndef OPM_FLOW_NONLINEAR_SYSTEM_HEADER_INCLUDED
#include <config.h>
#include <opm/simulators/flow/NonlinearSystem.hpp>
#endif

#include <dune/common/timer.hh>

#include <opm/common/ErrorMacros.hpp>
#include <opm/common/TimingMacros.hpp>

#include <cmath>
#include <stdexcept>
#include <string>
#include <utility>

namespace Opm {

template <class TypeTag>
SimulatorReportSingle
NonlinearSystem<TypeTag>::
assembleReservoir(const SimulatorTimerInterface&)
{
    return assembleReservoir(well_model_);
}

template <class TypeTag>
void
NonlinearSystem<TypeTag>::
updateTUNING(const Tuning& tuning)
{
    applyTUNING(param_, tuning);
}

template <class TypeTag>
void
NonlinearSystem<TypeTag>::
updateTUNINGDP(const TuningDp& tuning_dp)
{
    applyTUNINGDP(param_, tuning_dp);
}

template <class TypeTag>
void
NonlinearSystem<TypeTag>::
updateSolution(const GlobalEqVector& dx)
{
    OPM_TIMEBLOCK(updateSolution);

    const bool shouldStore = shouldStoreSolutionUpdate();
    if (shouldStore) {
        prepareSolutionUpdate();
    }

    auto& newtonMethod = simulator_.model().newtonMethod();
    auto& solution = simulator_.model().solution(/*timeIdx=*/0);

    newtonMethod.applyUpdate(/*nextSolution=*/solution,
                             /*curSolution=*/solution,
                             /*update=*/dx,
                             /*resid=*/dx);

    {
        OPM_TIMEBLOCK(invalidateAndUpdateIntensiveQuantities);
        simulator_.model().invalidateAndUpdateIntensiveQuantities(/*timeIdx=*/0);
    }

    if (shouldStore) {
        storeSolutionUpdate(dx);
    }
}

template <class TypeTag>
template <class LogFailure>
void
NonlinearSystem<TypeTag>::
addReservoirConvergenceMetrics(ConvergenceReport& report,
                               const int componentIdx,
                               const std::string_view componentName,
                               const std::span<const Scalar> residuals,
                               const std::span<const ConvergenceReport::ReservoirFailure::Type> types,
                               const std::span<const Scalar> tolerances,
                               const Scalar maxResidualAllowed,
                               LogFailure&& logFailure) const
{
    if (residuals.size() != types.size() || residuals.size() != tolerances.size()) {
        OPM_THROW(std::logic_error, "Mismatched reservoir convergence metric sizes.");
    }

    using CR = ConvergenceReport;
    for (std::size_t metricIdx = 0; metricIdx < residuals.size(); ++metricIdx) {
        const auto residual = residuals[metricIdx];
        const auto type = types[metricIdx];
        const auto tolerance = tolerances[metricIdx];

        if (std::isnan(residual)) {
            report.setReservoirFailed({type, CR::Severity::NotANumber, componentIdx});
            logFailure("NaN residual for " + std::string(componentName) + " equation.");
        }
        else if (residual > maxResidualAllowed) {
            report.setReservoirFailed({type, CR::Severity::TooLarge, componentIdx});
            logFailure("Too large residual for " + std::string(componentName) + " equation.");
        }
        else if (residual < 0.0) {
            report.setReservoirFailed({type, CR::Severity::Normal, componentIdx});
            logFailure("Negative residual for " + std::string(componentName) + " equation.");
        }
        else if (residual > tolerance) {
            report.setReservoirFailed({type, CR::Severity::Normal, componentIdx});
        }

        report.setReservoirConvergenceMetric(type, componentIdx, residual, tolerance);
    }
}

template <class TypeTag>
NonlinearSystem<TypeTag>::
NonlinearSystem(Simulator& simulator,
                const ModelParameters& param,
                WellModel& wellModel,
                const bool terminal_output)
    : simulator_(simulator)
    , grid_(simulator_.vanguard().grid())
    , terminal_output_(terminal_output)
    , param_(param)
    , well_model_(wellModel)
    , current_relaxation_(1.0)
    , dx_old_(simulator_.model().numGridDof())
{}

template <class TypeTag>
void
NonlinearSystem<TypeTag>::
initialLinearization(SimulatorReportSingle& report,
                     const int,
                     const int,
                     const SimulatorTimerInterface&)
{
    failureReport_ = SimulatorReportSingle();

    Dune::Timer perfTimer;
    perfTimer.start();
    report.total_linearizations = 1;

    try {
        report += this->assembleReservoir(well_model_);
        report.assemble_time += perfTimer.stop();

        // Mark timestep initialized after assembling, because well models can use
        // needsTimestepInit() to trigger per-step setup during assembly.
        simulator_.problem().markTimestepInitialized();
    }
    catch (...) {
        report.assemble_time += perfTimer.stop();
        failureReport_ += report;
        throw;
    }
}

template <class TypeTag>
SimulatorReportSingle
NonlinearSystem<TypeTag>::
prepareStep(const SimulatorTimerInterface& timer)
{
    SimulatorReportSingle report;
    Dune::Timer perfTimer;
    perfTimer.start();

    const int lastStepFailed = timer.lastStepFailed();
    if (grid_.comm().size() > 1 && grid_.comm().max(lastStepFailed) != grid_.comm().min(lastStepFailed)) {
        OPM_THROW(std::runtime_error,
                  "Misalignment of the parallel simulation run in prepareStep "
                  "- the previous step succeeded on some ranks but failed on others.");
    }

    if (lastStepFailed) {
        simulator_.problem().updateFailed();
        simulator_.model().newtonMethod().eraseMatrix();
    }
    else {
        simulator_.problem().advanceTimeLevel();
    }

    // The model still needs the report-step time context even though flow owns time stepping.
    simulator_.setTime(timer.simulationTimeElapsed());
    simulator_.setTimeStepSize(timer.currentStepLength());

    simulator_.problem().resetIterationForNewTimestep();
    simulator_.problem().beginTimeStep();

    report.pre_post_time += perfTimer.stop();
    return report;
}

template <class TypeTag>
template <class WellModelType>
SimulatorReportSingle
NonlinearSystem<TypeTag>::
assembleReservoir(WellModelType& wellModel)
{
    simulator_.problem().beginIteration();
    simulator_.model().linearizer().linearizeDomain();
    simulator_.problem().endIteration();
    return wellModel.lastReport();
}

template <class TypeTag>
template <class ModelParametersType>
void
NonlinearSystem<TypeTag>::
applyTUNING(ModelParametersType& param,
            const Tuning& tuning)
{
    param.tolerance_cnv_ = tuning.TRGCNV;
    param.tolerance_cnv_relaxed_ = tuning.XXXCNV;
    param.tolerance_mb_ = tuning.TRGMBE;
    param.tolerance_mb_relaxed_ = tuning.XXXMBE;
    param.newton_max_iter_ = tuning.NEWTMX;
    param.newton_min_iter_ = tuning.NEWTMN;
}

template <class TypeTag>
template <class ModelParametersType>
void
NonlinearSystem<TypeTag>::
applyTUNINGDP(ModelParametersType& param,
              const TuningDp& tuning_dp)
{
    param.tolerance_max_dp_ = tuning_dp.TRGDDP;
    param.tolerance_max_ds_ = tuning_dp.TRGDDS;
    param.tolerance_max_drs_ = tuning_dp.TRGDDRS;
    param.tolerance_max_drv_ = tuning_dp.TRGDDRV;
}

template <class TypeTag>
template <class ValueType>
std::tuple<ValueType, ValueType>
NonlinearSystem<TypeTag>::
convergenceReduction(Parallel::Communication comm,
                     const ValueType primaryVolumeLocal,
                     const ValueType secondaryVolumeLocal,
                     std::vector<ValueType>& sumValues,
                     std::vector<ValueType>& maxValues,
                     std::vector<ValueType>& averagedValues)
{
    OPM_TIMEBLOCK(convergenceReduction);

    ValueType primaryVolume = primaryVolumeLocal;
    ValueType secondaryVolume = secondaryVolumeLocal;

    if (comm.size() > 1) {
        std::vector<ValueType> sumBuffer;
        std::vector<ValueType> maxBuffer;
        const int numComp = averagedValues.size();
        sumBuffer.reserve(2 * numComp + 2);
        maxBuffer.reserve(numComp);

        for (int compIdx = 0; compIdx < numComp; ++compIdx) {
            sumBuffer.push_back(averagedValues[compIdx]);
            sumBuffer.push_back(sumValues[compIdx]);
            maxBuffer.push_back(maxValues[compIdx]);
        }

        sumBuffer.push_back(primaryVolume);
        sumBuffer.push_back(secondaryVolume);

        comm.sum(sumBuffer.data(), sumBuffer.size());
        comm.max(maxBuffer.data(), maxBuffer.size());

        for (int compIdx = 0, buffIdx = 0; compIdx < numComp; ++compIdx, ++buffIdx) {
            averagedValues[compIdx] = sumBuffer[buffIdx];
            ++buffIdx;
            sumValues[compIdx] = sumBuffer[buffIdx];
        }

        for (int compIdx = 0; compIdx < numComp; ++compIdx) {
            maxValues[compIdx] = maxBuffer[compIdx];
        }

        primaryVolume = sumBuffer[sumBuffer.size() - 2];
        secondaryVolume = sumBuffer.back();
    }

    return {primaryVolume, secondaryVolume};
}

} // namespace Opm

#endif // OPM_FLOW_NONLINEAR_SYSTEM_IMPL_HEADER_INCLUDED
