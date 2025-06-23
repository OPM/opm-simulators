/*
*/
#ifndef OPM_ADAPTIVE_TIME_STEPPING_HPP
#define OPM_ADAPTIVE_TIME_STEPPING_HPP

#include <opm/models/utils/basicproperties.hh>
#include <opm/models/utils/propertysystem.hh>

#include <opm/simulators/timestepping/AdaptiveSimulatorTimer.hpp>
#include <opm/simulators/timestepping/SimulatorReport.hpp>
#include <opm/simulators/timestepping/SimulatorTimer.hpp>
#include <opm/simulators/timestepping/TimeStepControl.hpp>
#include <opm/simulators/timestepping/TimeStepControlInterface.hpp>

#if HAVE_MPI
#define RESERVOIR_COUPLING_ENABLED
#endif
#ifdef RESERVOIR_COUPLING_ENABLED
#include <opm/simulators/flow/ReservoirCoupling.hpp>
#include <opm/simulators/flow/ReservoirCouplingMaster.hpp>
#include <opm/simulators/flow/ReservoirCouplingSlave.hpp>
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

template<class TypeTag>
class AdaptiveTimeStepping
{
public:
    using TuningUpdateCallback = std::function<bool(const double, const double, const int)>;

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

    template <class Solver>
    class SubStepper {
    public:
        SubStepper(AdaptiveTimeStepping<TypeTag>& adaptive_time_stepping,
                   const SimulatorTimer& simulator_timer,
                   Solver& solver,
                   const bool is_event,
                   const TuningUpdateCallback& tuning_updater);

        AdaptiveTimeStepping<TypeTag>& getAdaptiveTimerStepper();
        SimulatorReport run();
        friend class SubStepIteration<Solver>;

    private:
        bool isReservoirCouplingMaster_() const;
        bool isReservoirCouplingSlave_() const;
        void maybeModifySuggestedTimeStepAtBeginningOfReportStep_(const double originalTimeStep);
        bool maybeUpdateTuning_(double elapsed, double dt, int sub_step_number) const;
        double maxTimeStep_() const;
        SimulatorReport runStepOriginal_();
#ifdef RESERVOIR_COUPLING_ENABLED
        ReservoirCouplingMaster& reservoirCouplingMaster_();
        ReservoirCouplingSlave& reservoirCouplingSlave_();
        SimulatorReport runStepReservoirCouplingMaster_();
        SimulatorReport runStepReservoirCouplingSlave_();
#endif
        double suggestedNextTimestep_() const;

        AdaptiveTimeStepping<TypeTag>& adaptive_time_stepping_;
        const SimulatorTimer& simulator_timer_;
        Solver& solver_;
        const bool is_event_;
        const TuningUpdateCallback& tuning_updater_;
    };

    template <class Solver>
    class SubStepIteration {
    public:
        SubStepIteration(SubStepper<Solver>& substepper,
                         AdaptiveSimulatorTimer& substep_timer,
                         const double original_time_step,
                         const bool final_step);

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
        void maybeReportSubStep_(SimulatorReportSingle substep_report) const;
        double maybeRestrictTimeStepGrowth_(const double dt,
                                            double dt_estimate,
                                            const int restarts) const;
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

    bool operator==(const AdaptiveTimeStepping<TypeTag>& rhs);

    static void registerParameters();
    void setSuggestedNextStep(const double x);
    double suggestedNextStep() const;
    const TimeStepControlInterface& timeStepControl() const;

    template <class Solver>
    SimulatorReport step(const SimulatorTimer& simulator_timer,
                         Solver& solver,
                         const bool is_event,
                         const TuningUpdateCallback& tuning_updater);

    void updateTUNING(double max_next_tstep, const Tuning& tuning);
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
