/*
  Copyright 2019 SINTEF Digital, Mathematics and Cybernetics.

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

#ifndef OPM_DEFERREDLOGGINGERRORHELPERS_HPP
#define OPM_DEFERREDLOGGINGERRORHELPERS_HPP

#include <opm/simulators/utils/DeferredLogger.hpp>

#include <dune/common/version.hh>
#include <dune/common/parallel/mpihelper.hh>

#include <string>
#include <sstream>
#include <exception>
#include <stdexcept>

// Macro to throw an exception.
// Inspired by ErrorMacros.hpp in opm-common.
// NOTE: For this macro to work, the
// exception class must exhibit a constructor with the signature
// (const std::string &message). Since this condition is not fulfilled
// for the std::exception, you should use this macro with some
// exception class derived from either std::logic_error or
// std::runtime_error.
//
// Usage: OPM_DEFLOG_THROW(ExceptionClass, "Error message " << value, DeferredLogger);
#define OPM_DEFLOG_THROW(Exception, message, deferred_logger)                             \
    do {                                                                \
        std::ostringstream oss__;                                       \
        oss__ << "[" << __FILE__ << ":" << __LINE__ << "] " << message; \
        deferred_logger.error(oss__.str());                               \
        throw Exception(oss__.str());                                   \
    } while (false)

inline void checkForExceptionsAndThrow(int exception_thrown, const std::string& message)
{
    const auto& cc = Dune::MPIHelper::getCollectiveCommunication();
    if (cc.max(exception_thrown) == 1) {
        throw std::runtime_error(message);
    }
}

inline void logAndCheckForExceptionsAndThrow(Opm::DeferredLogger& deferred_logger, int exception_thrown, const std::string& message, const bool terminal_output)
{
    Opm::DeferredLogger global_deferredLogger = gatherDeferredLogger(deferred_logger);
    if (terminal_output) {
        global_deferredLogger.logMessages();
    }
    // Now that all messages have been logged, they are automatically
    // cleared from the global logger, but we must also clear them
    // from the local logger.
    deferred_logger.clearMessages();
    const auto& cc = Dune::MPIHelper::getCollectiveCommunication();
    if (cc.max(exception_thrown) == 1) {
        throw std::runtime_error(message);
    }
}

#endif // OPM_DEFERREDLOGGINGERRORHELPERS_HPP
