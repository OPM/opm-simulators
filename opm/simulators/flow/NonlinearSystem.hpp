#ifndef OPM_FLOW_NONLINEAR_SYSTEM_HEADER_INCLUDED
#define OPM_FLOW_NONLINEAR_SYSTEM_HEADER_INCLUDED

#include <algorithm>

#include <cstddef>

#include <opm/common/TimingMacros.hpp>

#include <dune/common/timer.hh>

#include <opm/common/ErrorMacros.hpp>
#include <opm/common/OpmLog/OpmLog.hpp>

#include <opm/models/utils/propertysystem.hh>

#include <opm/simulators/timestepping/ConvergenceReport.hpp>
#include <opm/simulators/timestepping/SimulatorReport.hpp>
#include <opm/simulators/timestepping/SimulatorTimerInterface.hpp>
#include <opm/simulators/utils/ComponentName.hpp>
#include <opm/simulators/utils/ParallelCommunication.hpp>

#include <stdexcept>
#include <span>
#include <string>
#include <string_view>
#include <tuple>
#include <utility>
#include <vector>

namespace Opm {

template <class TypeTag>
class NonlinearSystem
{
public:
    using Simulator = GetPropType<TypeTag, Properties::Simulator>;
    using Grid = GetPropType<TypeTag, Properties::Grid>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using Indices = GetPropType<TypeTag, Properties::Indices>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using ComponentName = ::Opm::ComponentName<FluidSystem, Indices>;

    bool isParallel() const
    { return grid_.comm().size() > 1; }

    const Simulator& simulator() const
    { return simulator_; }

    Simulator& simulator()
    { return simulator_; }

    bool terminalOutputEnabled() const
    { return terminal_output_; }

    int numPhases() const
    { return Indices::numPhases; }

    const SimulatorReportSingle& failureReport() const
    { return failureReport_; }

    const std::vector<StepReport>& stepReports() const
    { return convergence_reports_; }

    const ComponentName& compNames() const
    { return compNames_; }

    void beginReportStep()
    { simulator_.problem().beginEpisode(); }

    void endReportStep()
    { simulator_.problem().endEpisode(); }

    template <class LogFailure>
    void addReservoirConvergenceMetrics(ConvergenceReport& report,
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
                std::forward<LogFailure>(logFailure)("NaN residual for " + std::string(componentName) + " equation.");
            }
            else if (residual > maxResidualAllowed) {
                report.setReservoirFailed({type, CR::Severity::TooLarge, componentIdx});
                std::forward<LogFailure>(logFailure)("Too large residual for " + std::string(componentName) + " equation.");
            }
            else if (residual < 0.0) {
                report.setReservoirFailed({type, CR::Severity::Normal, componentIdx});
                std::forward<LogFailure>(logFailure)("Negative residual for " + std::string(componentName) + " equation.");
            }
            else if (residual > tolerance) {
                report.setReservoirFailed({type, CR::Severity::Normal, componentIdx});
            }

            report.setReservoirConvergenceMetric(type, componentIdx, residual, tolerance);
        }
    }

protected:
    explicit NonlinearSystem(Simulator& simulator, const bool terminal_output)
        : simulator_(simulator)
        , grid_(simulator_.vanguard().grid())
        , terminal_output_(terminal_output)
    {}

    SimulatorReportSingle prepareStep(const SimulatorTimerInterface& timer)
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
            simulator_.model().updateFailed();
        }
        else {
            simulator_.model().advanceTimeLevel();
        }

        // The model still needs the report-step time context even though flow owns time stepping.
        simulator_.setTime(timer.simulationTimeElapsed());
        simulator_.setTimeStepSize(timer.currentStepLength());

        simulator_.problem().resetIterationForNewTimestep();
        simulator_.problem().beginTimeStep();

        report.pre_post_time += perfTimer.stop();
        return report;
    }

    template <class AssembleReservoir>
    void initialLinearization(SimulatorReportSingle& report,
                              const SimulatorTimerInterface& timer,
                              AssembleReservoir&& assembleReservoir)
    {
        failureReport_ = SimulatorReportSingle();

        Dune::Timer perfTimer;
        perfTimer.start();
        report.total_linearizations = 1;

        try {
            report += std::forward<AssembleReservoir>(assembleReservoir)(timer);
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

    template <class WellModel>
    SimulatorReportSingle assembleReservoir(WellModel& wellModel)
    {
        simulator_.problem().beginIteration();
        simulator_.model().linearizer().linearizeDomain();
        simulator_.problem().endIteration();
        return wellModel.lastReport();
    }

    template <class BVector, class PrepareSolutionUpdate, class StoreSolutionUpdate>
    void updateSolution(const BVector& dx,
                        const bool shouldStoreSolutionUpdate,
                        PrepareSolutionUpdate&& prepareSolutionUpdate,
                        StoreSolutionUpdate&& storeSolutionUpdate)
    {
        OPM_TIMEBLOCK(updateSolution);

        if (shouldStoreSolutionUpdate) {
            std::forward<PrepareSolutionUpdate>(prepareSolutionUpdate)();
        }

        auto& newtonMethod = simulator_.model().newtonMethod();
        auto& solution = simulator_.model().solution(/*timeIdx=*/0);

        newtonMethod.update_(/*nextSolution=*/solution,
                             /*curSolution=*/solution,
                             /*update=*/dx,
                             /*resid=*/dx);

        {
            OPM_TIMEBLOCK(invalidateAndUpdateIntensiveQuantities);
            simulator_.model().invalidateAndUpdateIntensiveQuantities(/*timeIdx=*/0);
        }

        if (shouldStoreSolutionUpdate) {
            std::forward<StoreSolutionUpdate>(storeSolutionUpdate)(dx);
        }
    }

    template <class StepInit, class IterationBody, class PostIteration>
    SimulatorReportSingle nonlinearIteration(const SimulatorTimerInterface& timer,
                                             const std::size_t reportReserve,
                                             StepInit&& stepInit,
                                             IterationBody&& iterationBody,
                                             PostIteration&& postIteration)
    {
        if (simulator_.problem().iterationContext().needsTimestepInit()) {
            std::forward<StepInit>(stepInit)();
            convergence_reports_.push_back({timer.reportStepNum(), timer.currentStepNum(), {}});
            convergence_reports_.back().report.reserve(reportReserve);
        }

        auto result = std::forward<IterationBody>(iterationBody)();
        std::forward<PostIteration>(postIteration)();

        simulator_.problem().advanceIteration();
        return result;
    }

    template <class BVector,
              class NonlinearSolverType,
              class ModelParameters,
              class WellModel,
              class ResidualNormsHistory,
              class ScalarValue,
              class InitialLinearization,
              class SolveJacobianSystem,
              class LinearIterationsLastSolve,
              class UpdateSolution>
    SimulatorReportSingle
    nonlinearIterationNewton(const SimulatorTimerInterface& timer,
                             NonlinearSolverType& nonlinearSolver,
                             const ModelParameters& param,
                             WellModel& wellModel,
                             ResidualNormsHistory& residualNormsHistory,
                             BVector& dxOld,
                             ScalarValue& currentRelaxation,
                             double& linearSolveSetupTime,
                             InitialLinearization&& initialLinearization,
                             SolveJacobianSystem&& solveJacobianSystem,
                             LinearIterationsLastSolve&& linearIterationsLastSolve,
                             UpdateSolution&& updateSolution)
    {
        OPM_TIMEFUNCTION();

        SimulatorReportSingle report;
        Dune::Timer perfTimer;

        std::forward<InitialLinearization>(initialLinearization)(report,
                                                                 param.newton_min_iter_,
                                                                 param.newton_max_iter_,
                                                                 timer);

        if (!report.converged) {
            perfTimer.reset();
            perfTimer.start();
            report.total_newton_iterations = 1;

            BVector x(simulator_.model().numGridDof());
            linearSolveSetupTime = 0.0;

            try {
                wellModel.linearize(simulator_.model().linearizer().jacobian(),
                                    simulator_.model().linearizer().residual());

                std::forward<SolveJacobianSystem>(solveJacobianSystem)(x);

                report.linear_solve_setup_time += linearSolveSetupTime;
                report.linear_solve_time += perfTimer.stop();
                report.total_linear_iterations += std::forward<LinearIterationsLastSolve>(linearIterationsLastSolve)();
            }
            catch (...) {
                report.linear_solve_setup_time += linearSolveSetupTime;
                report.linear_solve_time += perfTimer.stop();
                report.total_linear_iterations += std::forward<LinearIterationsLastSolve>(linearIterationsLastSolve)();

                failureReport_ += report;
                throw;
            }

            perfTimer.reset();
            perfTimer.start();

            wellModel.postSolve(x);

            if (param.use_update_stabilization_) {
                bool isOscillate = false;
                bool isStagnate = false;
                nonlinearSolver.detectOscillations(residualNormsHistory,
                                                   residualNormsHistory.size() - 1,
                                                   isOscillate,
                                                   isStagnate);

                if (isOscillate) {
                    currentRelaxation -= nonlinearSolver.relaxIncrement();
                    currentRelaxation = std::max(currentRelaxation, nonlinearSolver.relaxMax());

                    if (terminalOutputEnabled()) {
                        OpmLog::info("    Oscillating behavior detected: Relaxation set to "
                                     + std::to_string(currentRelaxation));
                    }
                }

                nonlinearSolver.stabilizeNonlinearUpdate(x, dxOld, currentRelaxation);
            }

            std::forward<UpdateSolution>(updateSolution)(x);
            report.update_time += perfTimer.stop();
        }

        return report;
    }

    template <class ValueType>
    std::tuple<ValueType, ValueType>
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

    void popLastStepReport()
    { convergence_reports_.back().report.pop_back(); }

    Simulator& simulator_;
    const Grid& grid_;
    bool terminal_output_;
    SimulatorReportSingle failureReport_;
    std::vector<StepReport> convergence_reports_;
    ComponentName compNames_{};
};

} // namespace Opm

#endif // OPM_FLOW_NONLINEAR_SYSTEM_HEADER_INCLUDED
