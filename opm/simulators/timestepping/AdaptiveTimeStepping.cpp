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
consistentlyFailingWells(const std::vector<StepReport>& sr)
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
    if (sr_size >= num_steps) {
        for (const auto& wf : wfs) {
            failing_wells.insert(wf.wellName());
        }
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
    Parameters::Register<Parameters::TimeStepSafetyFactor>
        ("Safety factor in the formula for the time step cutting after a "
        "time step has failed to satisfy the tolerance criterion");
}

} // namespace detail
} // namespace Opm
