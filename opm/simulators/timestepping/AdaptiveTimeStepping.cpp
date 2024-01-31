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

} // namespace detail
} // namespace Opm
