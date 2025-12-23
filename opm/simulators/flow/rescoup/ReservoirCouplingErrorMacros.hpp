/*
  Copyright 2025 Equinor ASA

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

#ifndef OPM_RESERVOIR_COUPLING_ERROR_MACROS_HPP
#define OPM_RESERVOIR_COUPLING_ERROR_MACROS_HPP

#include <opm/common/ErrorMacros.hpp>
#include <opm/simulators/utils/DeferredLoggingErrorHelpers.hpp>
#define RCOUP_LOG_THROW(Exception, message)                                       \
    do {                                                                          \
        if (this->logger().haveDeferredLogger()) {                                \
            OPM_DEFLOG_THROW(Exception, message, this->logger().deferredLogger()); \
        }                                                                         \
        else {                                                                    \
            OPM_THROW(Exception, message);                                        \
        }                                                                         \
    } while (false)
#endif // OPM_RESERVOIR_COUPLING_ERROR_MACROS_HPP
