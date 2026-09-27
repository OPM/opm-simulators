/*
  Copyright 2026 Equinor ASA

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

#ifndef OPM_RESERVOIR_COUPLING_SLAVE_TARGET_CHECK_HPP
#define OPM_RESERVOIR_COUPLING_SLAVE_TARGET_CHECK_HPP

#include <compare>
#include <map>
#include <set>
#include <string>
#include <vector>

namespace Opm {
class Schedule;
class SummaryConfig;
}

namespace Opm::ReservoirCoupling {

/// Startup check of a reservoir coupling slave run's input for group
/// target vectors the slave cannot report for its slave groups.
///
/// The targets of a slave group (a group named in GRUPSLAV) are set by the
/// master run.  The slave reports the master's injection targets through
/// GGIRT and GWIRT.  For the other group target vectors -- GOPRT, GWPRT,
/// GGPRT, GLPRT, GVPRT and GVIRT, and their field versions when FIELD is a
/// slave group -- it only knows its own schedule's value, which is usually
/// zero for a group the master controls.  A UDQ or an ACTIONX condition
/// built on such a vector would silently use the wrong number.
///
/// There is no problem when the slave group's GRUPSLAV flag for the
/// vector's phase is SLAV: the slave's own limit then applies, and its
/// schedule value is the right target.

/// One use of an unsupported target vector for a slave group.
struct UnsupportedTargetUse
{
    /// Where the vector is used, e.g., "UDQ GUOPT (DEFINE in X.DATA, line 12)".
    std::string source;

    /// Summary vector as written, e.g., GOPRT or FOPRT.
    std::string vector;

    /// Slave group the vector is evaluated for.
    std::string group;

    /// Memberwise ordering and equality.  The ordering lets the uses found
    /// at several report steps be collected in a std::set, which removes
    /// the duplicates and sorts the error message.
    auto operator<=>(const UnsupportedTargetUse&) const = default;
};

/// Stop a slave run whose input uses a target vector it cannot report.
///
/// Throws std::runtime_error, listing every use, when a UDQ DEFINE or an
/// ACTIONX condition uses one (see findUnsupportedSlaveGroupTargetUses()).
/// Otherwise logs one warning per unsupported vector requested in the
/// SUMMARY section (see findUnsupportedSlaveGroupTargetsInSummary()),
/// since those only affect the output.
///
/// \param[in] schedule Slave run's schedule.
/// \param[in] summary_config Slave run's summary vector configuration.
void checkSlaveGroupTargetVectors(const Schedule& schedule,
                                  const SummaryConfig& summary_config);

/// Uses of unsupported target vectors for slave groups in UDQ DEFINEs and
/// ACTIONX conditions, over all report steps.
///
/// A vector followed by a group name or a group name root (e.g., 'MANI*')
/// counts for the slave groups it matches.  A vector followed by no group
/// name counts for every slave group, and a field vector counts for FIELD.
///
/// \param[in] schedule Slave run's schedule.
///
/// \return Sorted uses without duplicates.  Empty when the schedule has no
/// GRUPSLAV, as in a history mode slave.
std::vector<UnsupportedTargetUse>
findUnsupportedSlaveGroupTargetUses(const Schedule& schedule);

/// Unsupported target vectors for slave groups requested in the SUMMARY
/// section.
///
/// \param[in] schedule Slave run's schedule.
/// \param[in] summary_config Slave run's summary vector configuration.
///
/// \return For each summary vector, the slave groups it is requested for.
std::map<std::string, std::set<std::string>>
findUnsupportedSlaveGroupTargetsInSummary(const Schedule& schedule,
                                          const SummaryConfig& summary_config);

} // namespace Opm::ReservoirCoupling

#endif // OPM_RESERVOIR_COUPLING_SLAVE_TARGET_CHECK_HPP
