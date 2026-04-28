/*
*/
#ifndef OPM_ADAPTIVE_TIME_STEPPING_HPP
#define OPM_ADAPTIVE_TIME_STEPPING_HPP

#include <opm/models/utils/basicproperties.hh>
#include <opm/models/utils/propertysystem.hh>

#include <opm/simulators/flow/rescoup/ReservoirCouplingEnabled.hpp>
#include <opm/simulators/timestepping/AdaptiveSimulatorTimer.hpp>
#include <opm/simulators/timestepping/SimulatorReport.hpp>
#include <opm/simulators/timestepping/SimulatorTimer.hpp>
#include <opm/simulators/timestepping/TimeStepControl.hpp>
#include <opm/simulators/timestepping/TimeStepControlInterface.hpp>

#ifdef RESERVOIR_COUPLING_ENABLED
#include <opm/simulators/flow/rescoup/ReservoirCoupling.hpp>
#include <opm/simulators/flow/rescoup/ReservoirCouplingMaster.hpp>
#include <opm/simulators/flow/rescoup/ReservoirCouplingSlave.hpp>
#endif

#include <functional>
#include <memory>
#include <set>
#include <string>
#include <tuple>
#include <vector>

namespace Opm::Parameters {

struct SolverContinueOnConvergenceFailure { static constexpr bool value = false; };
struct SolverMaxRestarts { static constexpr int value = 10; };
struct SolverVerbosity { static constexpr int value = 1; };
struct TimeStepVerbosity { static constexpr int value = 1; };
struct InitialTimeStepInDays { static constexpr double value = 1.0;  };
struct FullTimeStepInitially { static constexpr bool value = false; };
struct TimeStepControl { static constexpr auto value = "pid+newtoniteration"; };
struct TimeStepControlTolerance { static constexpr double value = 1e-1; };
struct TimeStepControlTargetIterations { static constexpr int value = 30; };
struct TimeStepControlTargetNewtonIterations { static constexpr int value = 8; };
struct TimeStepControlDecayRate { static constexpr double value = 0.75; };
struct TimeStepControlGrowthRate { static constexpr double  value = 1.25; };
struct TimeStepControlDecayDampingFactor { static constexpr double value = 1.0;  };
struct TimeStepControlGrowthDampingFactor { static constexpr double value = 3.2; };
struct TimeStepControlFileName { static constexpr auto value = "timesteps"; };
struct MinTimeStepBeforeShuttingProblematicWellsInDays { static constexpr double value = 0.01; };
struct MinTimeStepBasedOnNewtonIterations { static constexpr double value = 0.0; };
struct TimeStepControlSafetyFactor { static constexpr double value = 0.8; };
struct TimeStepControlRejectCompletedStep { static constexpr bool value = false; };
struct TimeStepControlToleranceTestVersion { static constexpr auto value = "standard"; };
struct TimeStepControlMaxReductionTimeStep { static constexpr double value = 0.1; };
struct TimeStepControlParameters { static constexpr auto value = "0.125;0.25;0.125;0.75;0.25"; };

} // namespace Opm::Parameters

namespace Opm {

struct Tuning;
class UnitSystem;
struct StepReport;

namespace detail {
    void logTimer(const AdaptiveSimulatorTimer& substep_timer);

    std::set<std::string>
    consistentlyFailingWells(const std::vector<StepReport>& sr,
                             bool requireRepeatedFailures);
    void registerAdaptiveParameters();

    std::tuple<TimeStepControlType, std::unique_ptr<TimeStepControlInterface>, bool>
    createController(const UnitSystem& unitSystem);
}

/** \brief Adaptive time-stepping coordinator for the black-oil simulator.
 *
 * Drives the substep loop inside each report step.  A report step is the
 * time interval between dates defined in the deck SCHEDULE.  Within a
 * report step, this class chooses a sequence of smaller substeps whose
 * lengths are adjusted to honour:
 *   - the solver's convergence behaviour (chop on failure, grow on
 *     success);
 *   - TUNING keywords and ACTIONX-driven NEXTSTEP updates;
 *   - WCYCLE well-cycling schedules;
 *   - (when reservoir coupling is enabled) master/slave sync dates.
 *
 * The user-facing entry point is \ref step().  The substep loop itself
 * is split across the nested helpers \ref SubStepper (which selects
 * between the original loop and the two reservoir-coupling variants) and
 * \ref SubStepIteration (which runs the inner per-substep logic).
 */
template<class TypeTag>
class AdaptiveTimeStepping
{
public:
    /** \brief Callback invoked at the start of each substep to apply
     *         TUNING, NEXTSTEP (via ACTIONX), and WCYCLE updates.
     *
     * Called once by `SimulatorFullyImplicitBlackoil::runStep` at the
     * start of a report step, and once per substep inside
     * \ref SubStepIteration::maybeUpdateTuningAndTimeStep_ (plus equivalent
     * sites in the reservoir-coupling loops).  The callback may adjust
     * the simulator's suggested next substep via
     * \ref updateNEXTSTEP() or \ref updateTUNING().
     *
     * \param elapsed          Simulation time elapsed at the call point [s].
     * \param substep_length   Candidate length of the substep about to run
     *                         [s].  This is the look-ahead horizon the
     *                         callback uses internally, e.g. WCYCLE chops
     *                         the substep if a well switching event falls
     *                         inside `[elapsed, elapsed + substep_length]`.
     * \param sub_step_number  Index of the substep within its report step;
     *                         0 means "first substep".  WCYCLE's
     *                         `REQUEST_OPEN_WELL` handling is only consulted
     *                         when this is 0.
     * \return true if the callback modified `suggested_next_timestep_`.
     */
    using TuningUpdateCallback = std::function<bool(double elapsed,
                                                    double substep_length,
                                                    int sub_step_number)>;

private:
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;

    template <class Solver>
    class SolutionTimeErrorSolverWrapper : public RelativeChangeInterface
    {
    public:
        explicit SolutionTimeErrorSolverWrapper(const Solver& solver);
        double relativeChange() const;

    private:
        const Solver& solver_;
    };

    // Forward declaration of SubStepIteration
    template <class Solver> class SubStepIteration;

    /** \brief Per-report-step substepping driver.
     *
     * Owns the connection between the outer \ref SimulatorTimer (whose
     * step granularity is the report step) and the per-substep state.
     * `run()` selects the substepping loop to use based on whether
     * reservoir coupling is active for this run:
     *   - \ref runStepOriginal_ for standalone (non-rescoup) runs,
     *   - \ref runStepReservoirCouplingMaster_ on the master process of
     *     a coupled run,
     *   - \ref runStepReservoirCouplingSlave_ on a slave process.
     *
     * All three variants delegate the inner per-substep work to
     * \ref SubStepIteration::run().
     */
    template <class Solver>
    class SubStepper {
    public:
        SubStepper(AdaptiveTimeStepping<TypeTag>& adaptive_time_stepping,
                   const SimulatorTimer& simulator_timer,
                   Solver& solver,
                   const bool is_event,
                   const TuningUpdateCallback& tuning_updater);

        AdaptiveTimeStepping<TypeTag>& getAdaptiveTimerStepper();

        /// Dispatch to the appropriate substep loop and run it to
        /// completion (end of report step).  Returns the aggregated
        /// per-substep report.
        SimulatorReport run();
        friend class SubStepIteration<Solver>;

    private:
        bool isReservoirCouplingMaster_() const;
        bool isReservoirCouplingSlave_() const;

        /** \brief Apply start-of-report-step overrides to
         *         `suggested_next_timestep_`.
         *
         * Three overrides, checked in order:
         *   - if the suggestion is negative (uninitialised), set it to
         *     `restart_factor_ * original_time_step`;
         *   - if `full_timestep_initially_` is set, use the full report-step
         *     length;
         *   - if this report step is an event and `timestep_after_event_` is
         *     positive, use `timestep_after_event_`.
         * No-op otherwise.
         *
         * \param originalTimeStep Full report-step length [s].
         */
        void maybeModifySuggestedTimeStepAtBeginningOfReportStep_(const double originalTimeStep);

        /// Invoke the stored \ref TuningUpdateCallback.  Thin
        /// passthrough; see the callback alias for the argument contract.
        bool maybeUpdateTuning_(double elapsed, double substep_length, int sub_step_number) const;

        double maxTimeStep_() const;

        /// Standalone (non-rescoup) substep loop.  One pass through
        /// \ref SubStepIteration covers the whole report step.
        SimulatorReport runStepOriginal_();
#ifdef RESERVOIR_COUPLING_ENABLED
        // Reservoir coupling master: pick the sync-step length for the next outer-loop
        // iteration of `runStepReservoirCouplingMaster_()`, including the chop
        // against slave-report boundaries.  See the helper's own comment for
        // details.
        double getRcMasterSyncStepLength_(double prev_step, double current_time, double step_end_time);
        ReservoirCouplingMaster<Scalar>& reservoirCouplingMaster_();
        ReservoirCouplingSlave<Scalar>& reservoirCouplingSlave_();

        /// Reservoir-coupling master substep loop.  Each outer iteration
        /// is one "sync step" (one master-slave exchange); inside, a
        /// \ref SubStepIteration handles the substeps within that sync
        /// step.
        SimulatorReport runStepReservoirCouplingMaster_();

        /// Reservoir-coupling slave substep loop.  Receives the next
        /// sync-step length from the master, runs it to completion, and
        /// iterates until the report-step end is reached.
        SimulatorReport runStepReservoirCouplingSlave_();
#endif
        double suggestedNextTimestep_() const;

        AdaptiveTimeStepping<TypeTag>& adaptive_time_stepping_;
        const SimulatorTimer& simulator_timer_;
        Solver& solver_;
        const bool is_event_;
        const TuningUpdateCallback& tuning_updater_;
    };

    /** \brief Inner substep loop.
     *
     * Runs the sequence of substeps required to cover a given interval
     * (either a full report step, for \ref SubStepper::runStepOriginal_,
     * or a single sync step, for the reservoir-coupling variants).
     * Handles:
     *   - solver restarts and time-step chopping on convergence failure;
     *   - time-step growth control on convergence success;
     *   - per-substep tuning/WCYCLE updates via
     *     \ref maybeUpdateTuningAndTimeStep_;
     *   - per-substep SUMMARY and restart output.
     */
    template <class Solver>
    class SubStepIteration {
    public:
        SubStepIteration(SubStepper<Solver>& substepper,
                         AdaptiveSimulatorTimer& substep_timer,
                         const double original_time_step,
                         const bool final_step);

        /// Run substeps until `substep_timer_` is done.  Returns the
        /// aggregated per-substep report.
        SimulatorReport run();

    private:
        bool checkContinueOnUnconvergedSolution_(double dt) const;
        void checkTimeStepMaxRestartLimit_(const int restarts) const;
        void checkTimeStepMinLimit_(const double new_time_step) const;
        void chopTimeStep_(const double new_time_step);
        bool chopTimeStepOrCloseFailingWells_(const double new_time_step);
        boost::posix_time::ptime currentDateTime_() const;
        int getNumIterations_(const SimulatorReportSingle& substep_report) const;
        double growthFactor_() const;
        bool ignoreConvergenceFailure_() const;
        bool isReservoirCouplingMaster_() const;
        bool isReservoirCouplingSlave_() const;
        void markFirstSubStepAsFinished_() const;
        void maybeReportSubStep_(SimulatorReportSingle substep_report) const;
        double maybeRestrictTimeStepGrowth_(const double dt,
                                            double dt_estimate,
                                            const int restarts) const;
        void maybeUpdateLastSubstepOfSyncTimestep_(double dt);

        /** \brief Per-substep tuning/WCYCLE update with save-and-restore
         *         of `suggested_next_timestep_`.
         *
         * Invokes \ref SubStepper::maybeUpdateTuning_ with the current
         * substep's elapsed time, length, and index.  If the callback
         * reports a change (e.g. WCYCLE chopped the step via
         * `updateNEXTSTEP`), the new suggestion is applied to the
         * current substep via `setTimeStep_`, then the original
         * suggestion is restored so it can drive the next substep.
         *
         * \note The name is a historical artefact: by this point in the
         *       report step, TUNING_CHANGE/TUNINGDP_CHANGE events have
         *       already been cleared by the initial callback invocation
         *       in `SimulatorFullyImplicitBlackoil::runStep`, so what
         *       actually runs here is the WCYCLE/NEXTSTEP branch of the
         *       callback.
         */
        void maybeUpdateTuningAndTimeStep_();
        double maxGrowth_() const;
        double minTimeStepBeforeClosingWells_() const;
        double minTimeStep_() const;
        double restartFactor_() const;
        SimulatorReportSingle runSubStep_();
        int solverRestartMax_() const;
        double suggestedNextTimestep_() const;
        void setSuggestedNextStep_(double step);
        void setTimeStep_(double dt_estimate);
        Solver& solver_() const;
        bool solverVerbose_() const;
        const SimulatorTimer& simulatorTimer_() const;
        boost::posix_time::ptime startDateTime_() const;
#ifdef RESERVOIR_COUPLING_ENABLED
        ReservoirCouplingMaster<Scalar>& reservoirCouplingMaster_() const;
        ReservoirCouplingSlave<Scalar>& reservoirCouplingSlave_() const;
#endif
        double timeStepControlComputeEstimate_(const double dt,
                                               const int iterations,
                                               const AdaptiveSimulatorTimer& substepTimer) const;
        bool timeStepVerbose_() const;
        void updateSuggestedNextStep_();
        bool useNewtonIteration_() const;
        double writeOutput_() const;

        SubStepper<Solver>& substepper_;
        AdaptiveSimulatorTimer& substep_timer_;
        const double original_time_step_;
        const bool final_step_;
        std::string cause_of_failure_;
        AdaptiveTimeStepping<TypeTag>& adaptive_time_stepping_;
    };

public:
    AdaptiveTimeStepping() = default;

    AdaptiveTimeStepping(const UnitSystem& unitSystem,
                         const SimulatorReport& full_report,
                         const double max_next_tstep = -1.0,
                         const bool terminalOutput = true);

    AdaptiveTimeStepping(double max_next_tstep,
                         const Tuning& tuning,
                         const UnitSystem& unitSystem,
                         const SimulatorReport& full_report,
                         const bool terminalOutput = true);

    bool operator==(const AdaptiveTimeStepping<TypeTag>& rhs) const;

    static void registerParameters();

    /// Set the suggested length for the next substep [s].
    void setSuggestedNextStep(const double x);

    /// Current suggested length for the next substep [s].  Updated by
    /// the time-step controller, by TUNING/NEXTSTEP application, and by
    /// WCYCLE chopping.
    double suggestedNextStep() const;

    const TimeStepControlInterface& timeStepControl() const;

    /** \brief Run one report step by orchestrating adaptive substepping.
     *
     * Splits the report step into a sequence of smaller substeps whose
     * lengths are chosen adaptively, invoking `solver` once per substep.
     * Called from `SimulatorFullyImplicitBlackoil::runStep` in place of a
     * direct `Solver::step()` call on the non-adaptive path.  The actual
     * substep loop lives one layer down in \ref SubStepper::run(); this
     * method just constructs a \ref SubStepper and delegates.
     *
     * \param simulator_timer  Outer report-step timer.
     * \param solver           Newton solver invoked once per substep.
     * \param is_event         True if this report step carries an event
     *                         (NEW_WELL, INJECTION_UPDATE, etc.); affects
     *                         \ref SubStepper::maybeModifySuggestedTimeStepAtBeginningOfReportStep_.
     * \param tuning_updater   Callback invoked to apply TUNING / NEXTSTEP
     *                         / WCYCLE updates; see
     *                         \ref TuningUpdateCallback.
     */
    template <class Solver>
    SimulatorReport step(const SimulatorTimer& simulator_timer,
                         Solver& solver,
                         const bool is_event,
                         const TuningUpdateCallback& tuning_updater);

    /** \brief Apply TUNING keyword parameters.
     *
     * Overwrites `restart_factor_`, `growth_factor_`, `max_growth_`,
     * `max_time_step_`, and `timestep_after_event_` from `tuning`, then
     * forwards `max_next_tstep` to \ref updateNEXTSTEP().
     *
     * \param max_next_tstep   Next-step suggestion [s]; ignored if <= 0.
     * \param tuning           Source of the new TUNING parameters.
     */
    void updateTUNING(double max_next_tstep, const Tuning& tuning);

    /** \brief Set `suggested_next_timestep_` to `max_next_tstep` iff
     *         `max_next_tstep > 0`.
     *
     * Used by both the TUNING path and the WCYCLE path in the
     * tuning-update callback.  A non-positive value is interpreted as
     * "no suggestion" and is silently ignored.
     */
    void updateNEXTSTEP(double max_next_tstep);

    template<class Serializer>
    void serializeOp(Serializer& serializer);

    SimulatorReport& report();

    static AdaptiveTimeStepping<TypeTag> serializationTestObjectHardcoded();
    static AdaptiveTimeStepping<TypeTag> serializationTestObjectPID();
    static AdaptiveTimeStepping<TypeTag> serializationTestObjectPIDIt();
    static AdaptiveTimeStepping<TypeTag> serializationTestObjectSimple();
    static AdaptiveTimeStepping<TypeTag> serializationTestObject3rdOrder();

private:
    void maybeModifySuggestedTimeStepAtBeginningOfReportStep_(const double original_time_step,
                                                              const bool is_event);

    template<class Controller>
    static AdaptiveTimeStepping<TypeTag> serializationTestObject_();

    template<class T, class Serializer>
    void allocAndSerialize(Serializer& serializer);

    template<class T>
    bool castAndComp(const AdaptiveTimeStepping<TypeTag>& Rhs) const;

protected:
    void init_(const UnitSystem& unitSystem);

    using TimeStepController = std::unique_ptr<TimeStepControlInterface>;

    /// type of time step control object
    TimeStepControlType time_step_control_type_{TimeStepControlType::PIDAndIterationCount};
    TimeStepController time_step_control_{}; //!< time step control object
    double restart_factor_{};                //!< factor to multiply time step with when solver fails to converge
    double growth_factor_{};                 //!< factor to multiply time step when solver recovered from failed convergence
    double max_growth_{};                    //!< factor that limits the maximum growth of a time step
    double max_time_step_{};                 //!< maximal allowed time step size in days
    double min_time_step_{};                 //!< minimal allowed time step size before throwing
    bool ignore_convergence_failure_{false}; //!< continue instead of stop when minimum time step is reached
    int solver_restart_max_{};               //!< how many restart of solver are allowed
    bool solver_verbose_{false};             //!< solver verbosity
    bool timestep_verbose_{false};           //!< timestep verbosity
    double suggested_next_timestep_{};       //!< suggested size of next timestep
    bool full_timestep_initially_{false};    //!< beginning with the size of the time step from data file
    double timestep_after_event_{};          //!< suggested size of timestep after an event
    bool use_newton_iteration_{false};       //!< use newton iteration count for adaptive time step control

    //! < shut problematic wells when time step size in days are less than this
    double min_time_step_before_shutting_problematic_wells_{};
    // We store a copy of the full simulator run report for output purposes,
    // so it can be updated and passed to the summary writing code every
    // substep (not just every report step).
    SimulatorReport report_{};
};

} // namespace Opm

#include <opm/simulators/timestepping/AdaptiveTimeStepping_impl.hpp>

#endif // OPM_ADAPTIVE_TIME_STEPPING_HPP
