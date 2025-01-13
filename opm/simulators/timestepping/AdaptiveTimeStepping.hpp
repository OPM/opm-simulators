/*
*/
#ifndef OPM_ADAPTIVE_TIME_STEPPING_HPP
#define OPM_ADAPTIVE_TIME_STEPPING_HPP

#include <dune/common/version.hh>
#include <dune/istl/istlexception.hh>

#include <opm/common/OpmLog/OpmLog.hpp>

#include <opm/input/eclipse/Schedule/Tuning.hpp>

#include <opm/models/utils/basicproperties.hh>
#include <opm/models/utils/propertysystem.hh>

#include <opm/simulators/timestepping/AdaptiveSimulatorTimer.hpp>
#include <opm/simulators/timestepping/EclTimeSteppingParams.hpp>
#include <opm/simulators/timestepping/SimulatorReport.hpp>
#include <opm/simulators/timestepping/SimulatorTimer.hpp>
#include <opm/simulators/timestepping/TimeStepControl.hpp>
#include <opm/simulators/timestepping/TimeStepControlInterface.hpp>

#include <cmath>
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

} // namespace Opm::Parameters

namespace Opm {

class UnitSystem;
struct StepReport;

namespace detail {

void logTimer(const AdaptiveSimulatorTimer& substepTimer);

std::set<std::string> consistentlyFailingWells(const std::vector<StepReport>& sr);

void registerAdaptiveParameters();

std::tuple<TimeStepControlType,
           std::unique_ptr<TimeStepControlInterface>,
           bool>
createController(const UnitSystem& unitSystem);

}

    // AdaptiveTimeStepping
    //---------------------
    template<class TypeTag>
    class AdaptiveTimeStepping
    {
        using Scalar = GetPropType<TypeTag, Properties::Scalar>;
        template <class Solver>
        class SolutionTimeErrorSolverWrapper : public RelativeChangeInterface
        {
            const Solver& solver_;
        public:
            SolutionTimeErrorSolverWrapper(const Solver& solver)
              : solver_(solver)
            {}

            /// return || u^n+1 - u^n || / || u^n+1 ||
            double relativeChange() const
            { return solver_.model().relativeChange(); }
        };

        template<class E>
        void logException_(const E& exception, bool verbose)
        {
            if (verbose) {
                std::string message;
                message = "Caught Exception: ";
                message += exception.what();
                OpmLog::debug(message);
            }
        }

    public:
        AdaptiveTimeStepping() = default;

        //! \brief contructor taking parameter object
        explicit AdaptiveTimeStepping(const UnitSystem& unitSystem,
                                      const double max_next_tstep = -1.0,
                                      const bool terminalOutput = true);

        //! \brief contructor taking parameter object
        //! \param tuning Pointer to ecl TUNING keyword
        //! \param timeStep current report step
        AdaptiveTimeStepping(double max_next_tstep,
                             const Tuning& tuning,
                             const UnitSystem& unitSystem,
                             const bool terminalOutput = true);

        static void registerParameters();

        /** \brief  step method that acts like the solver::step method
                    in a sub cycle of time steps
            \param tuningUpdater Function used to update TUNING parameters before each
                                 time step. ACTIONX might change tuning.
        */
        template <class Solver>
        SimulatorReport step(const SimulatorTimer& simulatorTimer,
                             Solver& solver,
                             const bool isEvent,
                             const std::function<bool(const double, const double, const int)> tuningUpdater);

        /** \brief Returns the simulator report for the failed substeps of the last
         *         report step.
         */
        double suggestedNextStep() const
        { return suggestedNextTimestep_; }

        void setSuggestedNextStep(const double x)
        { suggestedNextTimestep_ = x; }

        void updateTUNING(double max_next_tstep, const Tuning& tuning);

        void updateNEXTSTEP(double max_next_tstep);

        template<class Serializer>
        void serializeOp(Serializer& serializer);

        static AdaptiveTimeStepping<TypeTag> serializationTestObjectHardcoded();
        static AdaptiveTimeStepping<TypeTag> serializationTestObjectPID();
        static AdaptiveTimeStepping<TypeTag> serializationTestObjectPIDIt();
        static AdaptiveTimeStepping<TypeTag> serializationTestObjectSimple();

        bool operator==(const AdaptiveTimeStepping<TypeTag>& rhs) const;

    private:
        template<class Controller>
        static AdaptiveTimeStepping<TypeTag> serializationTestObject_();

        template<class T, class Serializer>
        void allocAndSerialize(Serializer& serializer)
        {
            if (!serializer.isSerializing()) {
                timeStepControl_ = std::make_unique<T>();
            }
            serializer(*static_cast<T*>(timeStepControl_.get()));
        }

        template<class T>
        bool castAndComp(const AdaptiveTimeStepping<TypeTag>& Rhs) const
        {
            const T* lhs = static_cast<const T*>(timeStepControl_.get());
            const T* rhs = static_cast<const T*>(Rhs.timeStepControl_.get());
            return *lhs == *rhs;
        }

    protected:
        void init_(const UnitSystem& unitSystem);

        using TimeStepController = std::unique_ptr<TimeStepControlInterface>;

        TimeStepControlType timeStepControlType_{TimeStepControlType::PIDAndIterationCount}; //!< type of time step control object
        TimeStepController timeStepControl_{}; //!< time step control object
        double restartFactor_{};               //!< factor to multiply time step with when solver fails to converge
        double growthFactor_{};                //!< factor to multiply time step when solver recovered from failed convergence
        double maxGrowth_{};                   //!< factor that limits the maximum growth of a time step
        double maxTimeStep_{};                 //!< maximal allowed time step size in days
        double minTimeStep_{};                 //!< minimal allowed time step size before throwing
        bool ignoreConvergenceFailure_{false}; //!< continue instead of stop when minimum time step is reached
        int solverRestartMax_{};               //!< how many restart of solver are allowed
        bool solverVerbose_{false};            //!< solver verbosity
        bool timestepVerbose_{false};          //!< timestep verbosity
        double suggestedNextTimestep_{};       //!< suggested size of next timestep
        bool fullTimestepInitially_{false};    //!< beginning with the size of the time step from data file
        double timestepAfterEvent_{};          //!< suggested size of timestep after an event
        bool useNewtonIteration_{false};       //!< use newton iteration count for adaptive time step control
        double minTimeStepBeforeShuttingProblematicWells_{}; //! < shut problematic wells when time step size in days are less than this
    };
}

#include <opm/simulators/timestepping/AdaptiveTimeStepping_impl.hpp>

#endif // OPM_ADAPTIVE_TIME_STEPPING_HPP
