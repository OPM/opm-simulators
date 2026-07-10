/*
  Copyright 2026 Equinor ASA.

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

#ifndef OPM_ECONOMIC_LIMITS_MESSAGE_HEADER_INCLUDED
#define OPM_ECONOMIC_LIMITS_MESSAGE_HEADER_INCLUDED

#include <opm/common/utility/TimeService.hpp>

#include <fmt/chrono.h>
#include <fmt/format.h>

#include <ctime>
#include <string>

namespace Opm {

//! \brief Calendar date (DD-Mon-YYYY, UTC) reached at \p start_time plus
//!        \p sim_time seconds.
//!
//! Shared by the well (WECON) and group (GECON) economic-limit workover
//! messages so both report the same date basis. UTC (gmtime) is used
//! deliberately so the printed date does not depend on the host time zone.
inline std::string economicLimitDateString(const std::time_t start_time, const double sim_time)
{
    const std::time_t cur_time = TimeService::advance(start_time, sim_time);
    return fmt::format("{:%d-%b-%Y}", fmt::gmtime(cur_time));
}

//! \brief Separator line used to frame economic-limit workover messages.
//!        Built once and returned by reference (all callers use the same line).
inline const std::string& economicLimitMessageSeparator()
{
    static const std::string separator(110, '*');
    return separator;
}

} // namespace Opm

#endif // OPM_ECONOMIC_LIMITS_MESSAGE_HEADER_INCLUDED
