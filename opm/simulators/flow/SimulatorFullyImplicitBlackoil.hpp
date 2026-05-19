/*
  Copyright 2013, 2015, 2020 SINTEF Digital, Mathematics and Cybernetics.
  Copyright 2015 Andreas Lauser
  Copyright 2017 IRIS

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

#ifndef OPM_SIMULATOR_FULLY_IMPLICIT_BLACKOIL_HEADER_INCLUDED
#define OPM_SIMULATOR_FULLY_IMPLICIT_BLACKOIL_HEADER_INCLUDED

#include <opm/common/ErrorMacros.hpp>
#include <opm/simulators/flow/rescoup/ReservoirCouplingEnabled.hpp>

#ifdef RESERVOIR_COUPLING_ENABLED
#include <opm/input/eclipse/Schedule/ResCoup/ReservoirCouplingInfo.hpp>
#include <opm/input/eclipse/Schedule/ResCoup/MasterGroup.hpp>
#include <opm/input/eclipse/Schedule/ResCoup/Slaves.hpp>
#include <opm/simulators/flow/rescoup/ReservoirCouplingMaster.hpp>
#include <opm/simulators/flow/rescoup/ReservoirCouplingSlave.hpp>
#include <opm/common/Exceptions.hpp>
#endif

#include <opm/grid/utility/StopWatch.hpp>

#include <opm/simulators/aquifers/BlackoilAquiferModel.hpp>
#include <opm/simulators/flow/BlackoilModel.hpp>
#include <opm/simulators/flow/BlackoilModelParameters.hpp>
#include <opm/simulators/flow/ConvergenceOutputConfiguration.hpp>
#include <opm/simulators/flow/ExtraConvergenceOutputThread.hpp>
#include <opm/simulators/flow/NonlinearSolver.hpp>
#include <opm/simulators/flow/SimulatorConvergenceOutput.hpp>
#include <opm/simulators/flow/SimulatorReportBanners.hpp>
#include <opm/simulators/flow/SimulatorSerializer.hpp>
#include <opm/simulators/timestepping/AdaptiveTimeStepping.hpp>
#include <opm/simulators/timestepping/ConvergenceReport.hpp>
#include <opm/simulators/wells/WellState.hpp>

#if HAVE_HDF5
#include <opm/simulators/utils/HDF5Serializer.hpp>
#endif

#include <array>
#include <memory>
#include <string>
#include <vector>

namespace Opm::Parameters {

struct EnableAdaptiveTimeStepping { static constexpr bool value = true; };
struct OutputExtraConvergenceInfo { static constexpr auto* value = "none"; };
struct SaveStep { static constexpr auto* value = ""; };
struct SaveFile { static constexpr auto* value = ""; };
struct LoadFile { static constexpr auto* value = ""; };
struct LoadStep { static constexpr int value = -1; };
struct Slave { static constexpr bool value = false; };

} // namespace Opm::Parameters

namespace Opm::detail {

void registerSimulatorParameters();

/// \brief Log tuning parameters.
/// \param tuning Tuning values to log
/// \details Logs warnings if unsupported values are provided.
void logTuning(const Tuning& tuning);

}

namespace Opm {

/** \brief Top-level driver for a fully implicit black-oil simulation.
 *
 * Owns the per-report-step loop: \ref run repeatedly invokes
 * \ref runStep until `timer.done()` is reached.  Each \ref runStep
 * covers one report step (the interval between dates in the deck
 * SCHEDULE), either as a single \ref Solver::step call or, when
 * adaptive time stepping is enabled (the default), by delegating the
 * substep loop to \ref AdaptiveTimeStepping::step.
 *
 * Beyond the report-step loop, this class owns:
 *   - the \ref NonlinearSolver, constructed lazily on the first
 *     \ref runStep call;
 *   - per-report-step TUNING / TUNINGDP application via
 *     \ref updateTUNING and \ref updateTUNINGDP;
 *   - the WCYCLE-aware tuning-update callback handed to
 *     \ref AdaptiveTimeStepping;
 *   - simulation report aggregation and wall-clock timing;
 *   - the convergence-output thread lifecycle (INFOSTEP / INFOITER);
 *   - OPMRST save/load via \ref SimulatorSerializer;
 *   - reservoir-coupling master/slave bring-up and shutdown when
 *     RESERVOIR_COUPLING_ENABLED is set.
 *
 * The class is instantiated once per simulation process.  Static
 * configuration is registered via \ref registerParameters during
 * program startup.
 */
template<class TypeTag>
class SimulatorFullyImplicitBlackoil : private SerializableSim
{
protected:
    struct MPI_Comm_Deleter;
public:
    using Simulator = GetPropType<TypeTag, Properties::Simulator>;
    using Grid = GetPropType<TypeTag, Properties::Grid>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using ElementContext = GetPropType<TypeTag, Properties::ElementContext>;
    using BlackoilIndices = GetPropType<TypeTag, Properties::Indices>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using MaterialLaw = GetPropType<TypeTag, Properties::MaterialLaw>;
    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;
    using MaterialLawParams = GetPropType<TypeTag, Properties::MaterialLawParams>;
    using AquiferModel = GetPropType<TypeTag, Properties::AquiferModel>;
    using Model = GetPropType<TypeTag, Properties::NonlinearSystem>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;

    using TimeStepper = AdaptiveTimeStepping<TypeTag>;
    using PolymerModule = BlackOilPolymerModule<TypeTag>;
    using BioeffectsModule = BlackOilBioeffectsModule<TypeTag>;

    using Solver = NonlinearSolver<TypeTag, Model>;
    using ModelParameters = typename Model::ModelParameters;
    using SolverParameters = typename Solver::SolverParameters;
    using WellModel = BlackoilWellModel<TypeTag>;

    /** \brief Construct from the surrounding eWoms `Simulator`.
     *
     * Initialises the OPMRST serializer from the SaveStep / LoadStep /
     * SaveFile / LoadFile parameters and, on rank 0 with terminal
     * output enabled, starts the background convergence-output thread
     * that writes INFOSTEP / INFOITER files.
     *
     * \param simulator The surrounding eWoms simulator; observed, not owned.
     */
    explicit SimulatorFullyImplicitBlackoil(Simulator& simulator);

    /// Ends the convergence-output thread cleanly on all ranks.
    ~SimulatorFullyImplicitBlackoil() override;

    /** \brief Register all parameters consumed by this class and its
     *         major collaborators.
     *
     * Forwards to \ref ModelParameters::registerParameters,
     * \ref SolverParameters::registerParameters,
     * \ref AdaptiveTimeStepping::registerParameters,
     * \ref Opm::detail::registerSimulatorParameters, and the TPSA
     * Newton-method / linear-solver parameter sets.  Called once
     * during program startup.
     */
    static void registerParameters();

#ifdef RESERVOIR_COUPLING_ENABLED
    /** \brief Run the entire simulation to completion.
     *
     * Loops over report steps, calling \ref runStep on each one; on
     * exit, sends the reservoir-coupling shutdown signal (master sends
     * terminate; slave acknowledges) before MPI_Finalize.  Stops early
     * if the schedule triggers an EXIT keyword.
     *
     * \param timer Outer report-step timer; advanced once per report step.
     * \param argc  Process `argc`, forwarded to spawned slave processes.
     * \param argv  Process `argv`, forwarded to spawned slave processes.
     * \return      Aggregated simulation report (timing + per-step data).
     */
    SimulatorReport run(SimulatorTimer& timer, int argc, char** argv);

    /** \brief Detect whether this process should run as a
     *         reservoir-coupling master.
     *
     * Master mode is enabled when the schedule contains a SLAVES
     * keyword (master allocates rates if GRUPMAST is also present;
     * otherwise master only synchronises time-stepping).  GRUPMAST
     * without SLAVES is rejected as an inconsistent schedule.
     *
     * Should only be called when this process is *not* a slave (i.e.
     * `Parameters::Get<Parameters::Slave>()` is false).
     */
    bool checkRunningAsReservoirCouplingMaster();

    /** \brief One-shot setup performed before the first \ref runStep.
     *
     * Constructs the wall-clock timers, the \ref AdaptiveTimeStepping
     * instance (if adaptive stepping is enabled), and — for reservoir
     * coupling — either the master or slave coordination object.  The
     * slave-mode prologue exchanges initial sync data with the master;
     * the master-mode prologue stores `argc`/`argv` for later
     * `maybeSpawnSlaveProcesses` calls in \ref runStep.
     *
     * \param timer Report-step timer; only its current step is read here.
     * \param argc  Process `argc`, used when later spawning slave processes.
     * \param argv  Process `argv`, used when later spawning slave processes.
     */
    void init(const SimulatorTimer& timer, int argc, char** argv);
#else
    /** \brief Run the entire simulation to completion.
     *
     * Loops over report steps, calling \ref runStep on each one.
     * Stops early if the schedule triggers an EXIT keyword.
     *
     * \param timer Outer report-step timer; advanced once per report step.
     * \return      Aggregated simulation report (timing + per-step data).
     */
    SimulatorReport run(SimulatorTimer& timer);

    /** \brief One-shot setup performed before the first \ref runStep.
     *
     * Constructs the wall-clock timers and the \ref AdaptiveTimeStepping
     * instance (if adaptive stepping is enabled).  For restart runs,
     * the suggested next step is seeded from `Simulator::timeStepSize()`.
     *
     * \param timer Report-step timer; only its current step is read here.
     */
    void init(const SimulatorTimer& timer);
#endif

    /** \brief Apply a TUNING keyword to the cached model parameters.
     *
     * Overwrites convergence tolerances (TRGCNV / XXXCNV / TRGMBE /
     * XXXMBE) and Newton iteration limits (NEWTMX / NEWTMN) on the
     * \ref ModelParameters copy held by this class.  Items that are
     * recognised by the parser but not honoured by this simulator are
     * logged at warning level on rank 0.
     *
     * \note This method only updates the copy held here.  The same
     *       TUNING data is forwarded separately to the model held by
     *       the solver and to \ref AdaptiveTimeStepping by the
     *       tuning-update callback in \ref runStep.
     */
    void updateTUNING(const Tuning& tuning);

    /** \brief Apply a TUNINGDP keyword to the cached model parameters.
     *
     * Overwrites the maximum allowed pressure (TRGDDP), saturation
     * (TRGDDS), and dissolved-gas / vaporised-oil ratio (TRGDDRS,
     * TRGDDRV) changes per Newton iteration.  TRGLCV / XXXLCV are
     * accepted by the parser but logged as unsupported on rank 0.
     */
    void updateTUNINGDP(const TuningDp& tuning_dp);

    /** \brief Advance the simulation by one report step.
     *
     * Called by \ref run once per report step.  Performs:
     *   - early exit if the schedule requested EXIT;
     *   - OPMRST state load on the chosen restart step;
     *   - first-step init-state output;
     *   - lazy construction of the \ref Solver;
     *   - construction and dispatch of the WCYCLE / TUNING / TUNINGDP
     *     tuning-update callback consumed by \ref AdaptiveTimeStepping;
     *   - the actual solve, either via \ref AdaptiveTimeStepping::step
     *     (substep loop) or a single \ref Solver::step call;
     *   - end-of-step report output and OPMRST save.
     *
     * \param timer Report-step timer; advanced by one report step on success.
     * \return false to terminate the outer \ref run loop (e.g. EXIT
     *         keyword); true to continue with the next report step.
     */
    bool runStep(SimulatorTimer& timer);

    /** \brief Stop the timers and emit the final OPMRST output.
     *
     * Called by \ref run after the report-step loop finishes.  Stops
     * the total wall-clock timer and marks the report as converged.
     *
     * \return The final aggregated simulation report.
     */
    SimulatorReport finalize();

    const Grid& grid() const
    { return simulator_.vanguard().grid(); }

    /// Serialize the parts of this class needed for OPMRST round-tripping
    /// (the surrounding simulator state, the report, and the adaptive
    /// time stepper).
    template<class Serializer>
    void serializeOp(Serializer& serializer);

    const Model& model() const { return solver_->model(); }

protected:
    /// Load this simulator's data block from an OPMRST file via HDF5.
    void loadState(HDF5Serializer& serializer, const std::string& groupName) override;

    /// Save this simulator's data block to an OPMRST file via HDF5.
    void saveState(HDF5Serializer& serializer, const std::string& groupName) const override;

    /// Return the OPMRST header tuple: product name, module version,
    /// compile timestamp, deck case name, and parameter dump.
    std::array<std::string,5> getHeader() const override;

    /// Local-to-global cell index mapping.
    const std::vector<int>& getCellMapping() const override {
        return simulator_.vanguard().globalCell();
    }

    /** \brief Build the \ref Solver used during the current report step.
     *
     * Wraps a freshly constructed \ref Model in a \ref Solver.  Side
     * effect: if `write_partitions` is set in the model parameters,
     * the partition layout is written under
     * `<output_dir>/partition/<case>` on the first call only — the
     * flag is cleared after, so subsequent calls skip the dump.
     */
    std::unique_ptr<Solver> createSolver(WellModel& wellModel);

    const EclipseState& eclState() const { return simulator_.vanguard().eclState(); }

    const Schedule& schedule() const { return simulator_.vanguard().schedule(); }

    bool isRestart() const { return eclState().getInitConfig().restartRequested(); }

    WellModel& wellModel_() { return simulator_.problem().wellModel(); }

    const WellModel& wellModel_() const { return simulator_.problem().wellModel(); }

    /// Surrounding eWoms simulator; observed, not owned.
    Simulator& simulator_;

    /// Cached model parameters; mutated by TUNING / TUNINGDP application.
    ModelParameters modelParam_;

    /// Cached nonlinear-solver parameters.
    SolverParameters solverParam_;

    /// Built lazily on the first \ref runStep call; reused thereafter.
    std::unique_ptr<Solver> solver_;

    /// Emit high-level progress to std::cout (rank 0 only).
    bool terminalOutput_;

    /// Aggregated report across the entire simulation.
    SimulatorReport report_;

    /// Wall-clock for the current report step's solve.
    std::unique_ptr<time::StopWatch> solverTimer_;

    /// Wall-clock for the entire simulation.
    std::unique_ptr<time::StopWatch> totalTimer_;

    /// Set iff adaptive time stepping is enabled.
    std::unique_ptr<TimeStepper> adaptiveTimeStepping_;

    /// Background thread for INFOSTEP / INFOITER files.
    SimulatorConvergenceOutput convergence_output_{};

#ifdef RESERVOIR_COUPLING_ENABLED
    /// True iff this process runs as a reservoir-coupling slave.
    bool slaveMode_{false};

    /// Non-null iff this process is a reservoir-coupling master.
    std::unique_ptr<ReservoirCouplingMaster<Scalar>> reservoirCouplingMaster_{nullptr};

    /// Non-null iff this process is a reservoir-coupling slave.
    std::unique_ptr<ReservoirCouplingSlave<Scalar>> reservoirCouplingSlave_{nullptr};
#endif

    /// OPMRST save / load.
    SimulatorSerializer serializer_;
};

} // namespace Opm

#include <opm/simulators/flow/SimulatorFullyImplicitBlackoil_impl.hpp>

#endif // OPM_SIMULATOR_FULLY_IMPLICIT_BLACKOIL_HEADER_INCLUDED
