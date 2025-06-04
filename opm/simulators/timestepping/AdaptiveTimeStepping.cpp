/*
*/

#include <config.h>
#include <opm/simulators/timestepping/AdaptiveTimeStepping.hpp>

#include <opm/common/OpmLog/OpmLog.hpp>

#include <opm/input/eclipse/Units/Units.hpp>

#include <opm/simulators/timestepping/AdaptiveSimulatorTimer.hpp>
#include <opm/simulators/timestepping/ConvergenceReport.hpp>

#include <boost/date_time.hpp>

#include <set>
#include <sstream>

namespace Opm {
namespace detail {

void logTimer(const AdaptiveSimulatorTimer& substepTimer)
{
    std::ostringstream ss;
    boost::posix_time::time_facet* facet = new boost::posix_time::time_facet("%d-%b-%Y");
    ss.imbue(std::locale(std::locale::classic(), facet));
    ss <<"\nStarting time step " << substepTimer.currentStepNum() << ", stepsize "
                         << unit::convert::to(substepTimer.currentStepLength(), unit::day) << " days,"
                         << " at day " << (double)unit::convert::to(substepTimer.simulationTimeElapsed(), unit::day)
                         << "/" << (double)unit::convert::to(substepTimer.totalTime(), unit::day)
                         << ", date = " << substepTimer.currentDateTime();
    OpmLog::info(ss.str());
}

std::set<std::string>
consistentlyFailingWells(const std::vector<StepReport>& sr, bool requireRepeatedFailures)
{
    // If there are wells that cause repeated failures, we
    // close them, and restart the un-chopped timestep.
    std::ostringstream msg;
    msg << "    Excessive chopping detected in report step "
        << sr.back().report_step << ", substep " << sr.back().current_step << "\n";

    std::set<std::string> failing_wells;

    // return empty set if no report exists
    // well failures in assembly is not yet registred
    if (sr.back().report.empty())
        return failing_wells;

    const auto& wfs = sr.back().report.back().wellFailures();
    for (const auto& wf : wfs) {
        msg << "        Well that failed: " << wf.wellName() << "\n";
    }
    msg.flush();
    OpmLog::debug(msg.str());

    // Check the last few step reports.
    const int num_steps = 3;
    const int rep_step = sr.back().report_step;
    const int sub_step = sr.back().current_step;
    const int sr_size = sr.size();
    for (const auto& wf : wfs) {
        failing_wells.insert(wf.wellName());
    }
    if (requireRepeatedFailures && sr_size >= num_steps) {
        for (int s = 1; s < num_steps; ++s) {
            const auto& srep = sr[sr_size - 1 - s];
            // Report must be from same report step and substep, otherwise we have
            // not chopped/retried enough times on this step.
            if (srep.report_step != rep_step || srep.current_step != sub_step) {
                break;
            }
            // Get the failing wells for this step, that also failed all other steps.
            std::set<std::string> failing_wells_step;
            for (const auto& wf : srep.report.back().wellFailures()) {
                if (failing_wells.count(wf.wellName()) > 0) {
                    failing_wells_step.insert(wf.wellName());
                }
            }
            failing_wells.swap(failing_wells_step);
        }
    }
    return failing_wells;
}

void registerAdaptiveParameters()
{
    // TODO: make sure the help messages are correct (and useful)
    Parameters::Register<Parameters::SolverContinueOnConvergenceFailure>
        ("Continue instead of stop when minimum solver time step is reached");
    Parameters::Register<Parameters::SolverMaxRestarts>
        ("The maximum number of breakdowns before a substep is given up and "
         "the simulator is terminated");
    Parameters::Register<Parameters::SolverVerbosity>
        ("Specify the \"chattiness\" of the non-linear solver itself");
    Parameters::Register<Parameters::TimeStepVerbosity>
        ("Specify the \"chattiness\" during the time integration");
    Parameters::Register<Parameters::InitialTimeStepInDays>
        ("The size of the initial time step in days");
    Parameters::Register<Parameters::FullTimeStepInitially>
        ("Always attempt to finish a report step using a single substep");
    Parameters::Register<Parameters::TimeStepControl>
        ("The algorithm used to determine time-step sizes. "
         "Valid options are: "
         "'pid' (default), "
         "'pid+iteration', "
         "'pid+newtoniteration', "
         "'iterationcount', "
        "'newtoniterationcount' "
        "and 'hardcoded'");
    Parameters::Register<Parameters::TimeStepControlTolerance>
        ("The tolerance used by the time step size control algorithm");
    Parameters::Register<Parameters::TimeStepControlTargetIterations>
        ("The number of linear iterations which the time step control scheme "
         "should aim for (if applicable)");
    Parameters::Register<Parameters::TimeStepControlTargetNewtonIterations>
        ("The number of Newton iterations which the time step control scheme "
         "should aim for (if applicable)");
    Parameters::Register<Parameters::TimeStepControlDecayRate>
        ("The decay rate of the time step size of the number of "
         "target iterations is exceeded");
    Parameters::Register<Parameters::TimeStepControlGrowthRate>
        ("The growth rate of the time step size of the number of "
         "target iterations is undercut");
    Parameters::Register<Parameters::TimeStepControlDecayDampingFactor>
        ("The decay rate of the time step decrease when the "
         "target iterations is exceeded");
    Parameters::Register<Parameters::TimeStepControlGrowthDampingFactor>
        ("The growth rate of the time step increase when the "
         "target iterations is undercut");
    Parameters::Register<Parameters::TimeStepControlFileName>
        ("The name of the file which contains the hardcoded time steps sizes");
    Parameters::Register<Parameters::MinTimeStepBeforeShuttingProblematicWellsInDays>
        ("The minimum time step size in days for which problematic wells are not shut");
    Parameters::Register<Parameters::MinTimeStepBasedOnNewtonIterations>
        ("The minimum time step size (in days for field and metric unit and hours for lab unit) "
         "can be reduced to based on newton iteration counts");
    Parameters::Register<Parameters::TimeStepControlSafetyFactor>
        ("Value to be multiplied with the time step control tolerance to ensure that the target "
         "relative change is lower than the tolerance");
    Parameters::Register<Parameters::TimeStepControlRejectCompletedStep>
        ("(Only applicable for the general 3rd order controller.) Include rejection of completed "
         "time steps if the relative change is larger than the time step control tolerance");
    Parameters::Register<Parameters::TimeStepControlToleranceTestVersion>
        ("(Only applicable for the general 3rd order controller.) Ways to decide if the time step "
         "should be rejected. Options: 'standard' and 'control-error-filtering'. The standard "
         "version compares relative change to tolerance directly to decide if the time step should "
         "be rejected, while the control-error-filtering version compares the relative change "
         "in time step size to the max-reduction-time-step parameter.");
    Parameters::Register<Parameters::TimeStepControlMaxReductionTimeStep>
        ("(Only applicable for the general 3rd order controller, using 'control-error-filtering' "
         "as time-step-control-tolerance-test-version) If the (proposed) relative change in time "
         "step size is larger than this parameter, the time step will be rejected.");
    Parameters::Register<Parameters::TimeStepControlParameters>
        ("(Only applicable for the general 3rd order controller.) Parameters for the general "
         "3rd order controller. Should be given as 'beta_1;beta_2;beta_3;alpha_2;alpha_3'.");
}

std::tuple<TimeStepControlType, std::unique_ptr<TimeStepControlInterface>, bool>
createController(const UnitSystem& unitSystem)
{
    const double tol =  Parameters::Get<Parameters::TimeStepControlTolerance>(); // 1e-1
    using RetVal = std::tuple<TimeStepControlType, std::unique_ptr<TimeStepControlInterface>, bool>;
    using Func = std::function<RetVal()>;
    const auto creators = std::unordered_map<std::string, Func> {
        {"pid",
         [tol]() {
             return RetVal{
                 TimeStepControlType::PID,
                 std::make_unique<PIDTimeStepControl>(tol),
                 false
             };
        }},
        {"pid+iteration",
         [tol]() {
             const int iterations =  Parameters::Get<Parameters::TimeStepControlTargetIterations>(); // 30
             const double decayDampingFactor = Parameters::Get<Parameters::TimeStepControlDecayDampingFactor>(); // 1.0
             const double growthDampingFactor = Parameters::Get<Parameters::TimeStepControlGrowthDampingFactor>(); // 3.2
             return RetVal{
                 TimeStepControlType::PIDAndIterationCount,
                 std::make_unique<PIDAndIterationCountTimeStepControl>(iterations,
                                                                       decayDampingFactor,
                                                                       growthDampingFactor,
                                                                       tol),
                 false
             };
         }},
        {"pid+newtoniteration",
         [tol, &unitSystem]() {
             const int iterations =  Parameters::Get<Parameters::TimeStepControlTargetNewtonIterations>(); // 8
             const double decayDampingFactor = Parameters::Get<Parameters::TimeStepControlDecayDampingFactor>(); // 1.0
             const double growthDampingFactor = Parameters::Get<Parameters::TimeStepControlGrowthDampingFactor>(); // 3.2
             const double nonDimensionalMinTimeStepIterations = Parameters::Get<Parameters::MinTimeStepBasedOnNewtonIterations>(); // 0.0 by default
             // the min time step can be reduced by the newton iteration numbers
             double minTimeStepReducedByIterations = unitSystem.to_si(UnitSystem::measure::time,
                                                                      nonDimensionalMinTimeStepIterations);
             return RetVal{
                 TimeStepControlType::PIDAndIterationCount,
                 std::make_unique<PIDAndIterationCountTimeStepControl>(iterations,
                                                                       decayDampingFactor,
                                                                       growthDampingFactor,
                                                                       tol,
                                                                       minTimeStepReducedByIterations),
                 true
             };
         }},
        {"iterationcount",
         []() {
              const int iterations =  Parameters::Get<Parameters::TimeStepControlTargetIterations>(); // 30
              const double decayrate = Parameters::Get<Parameters::TimeStepControlDecayRate>(); // 0.75
              const double growthrate = Parameters::Get<Parameters::TimeStepControlGrowthRate>(); // 1.25
              return RetVal{
                  TimeStepControlType::SimpleIterationCount,
                  std::make_unique<SimpleIterationCountTimeStepControl>(iterations,
                                                                        decayrate,
                                                                        growthrate),
                  false
              };
         }},
        {"newtoniterationcount",
         []() {
             const int iterations =  Parameters::Get<Parameters::TimeStepControlTargetNewtonIterations>(); // 8
             const double decayrate = Parameters::Get<Parameters::TimeStepControlDecayRate>(); // 0.75
             const double growthrate = Parameters::Get<Parameters::TimeStepControlGrowthRate>(); // 1.25
             return RetVal{
                 TimeStepControlType::SimpleIterationCount,
                 std::make_unique<SimpleIterationCountTimeStepControl>(iterations,
                                                                       decayrate,
                                                                       growthrate),
                 true
             };
         }},
        {"general3rdorder",
         [tol]() {
             const double safetyFactor = Parameters::Get<Parameters::TimeStepControlSafetyFactor>();
             const bool rejectCompletedStep = Parameters::Get<Parameters::TimeStepControlRejectCompletedStep>();
             const std::string toleranceTestVersion = Parameters::Get<Parameters::TimeStepControlToleranceTestVersion>();
             const double maxReductionTimeStep = Parameters::Get<Parameters::TimeStepControlMaxReductionTimeStep>();
             const std::string parameters = Parameters::Get<Parameters::TimeStepControlParameters>();
             const bool verbose = Parameters::Get<Parameters::TimeStepVerbosity>();
             return RetVal{
                 TimeStepControlType::General3rdOrder,
                 std::make_unique<General3rdOrderController>(tol,
                                                             safetyFactor,
                                                             rejectCompletedStep,
                                                             toleranceTestVersion,
                                                             maxReductionTimeStep,
                                                             parameters,
                                                             verbose),
                 false
             };
        }},
        {"hardcoded",
         []() {
             const std::string filename = Parameters::Get<Parameters::TimeStepControlFileName>(); // "timesteps"
             return RetVal{
                 TimeStepControlType::HardCodedTimeStep,
                 std::make_unique<HardcodedTimeStepControl>(filename),
                 false
             };
         }},
    };

    const std::string control = Parameters::Get<Parameters::TimeStepControl>(); // "pid"
    const auto it = creators.find(control);
    if (it == creators.end()) {
        OPM_THROW(std::runtime_error,
                  "Unsupported time step control selected " + control);
    }

    // invoke creator
    return it->second();
}

} // namespace detail
} // namespace Opm
