/*
*/

#include <config.h>
#include <opm/simulators/timestepping/AdaptiveTimeSteppingEbos.hpp>

#include <opm/common/OpmLog/OpmLog.hpp>
#include <opm/input/eclipse/Units/Units.hpp>
#include <opm/simulators/timestepping/AdaptiveSimulatorTimer.hpp>
#include <boost/date_time.hpp>
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

} // namespace detail
} // namespace Opm
