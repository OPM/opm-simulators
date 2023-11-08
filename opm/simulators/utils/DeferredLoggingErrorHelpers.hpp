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

#include <opm/common/Exceptions.hpp>

#include <opm/simulators/utils/DeferredLogger.hpp>
#include <opm/simulators/utils/gatherDeferredLogger.hpp>

#include <opm/simulators/utils/ParallelCommunication.hpp>

#include <string>
#include <exception>
#include <stdexcept>

// Macro to log an error and throw an exception.
// Inspired by ErrorMacros.hpp in opm-common.
// NOTE: For this macro to work, the
// exception class must exhibit a constructor with the signature
// (const std::string &message). Since this condition is not fulfilled
// for the std::exception, you should use this macro with some
// exception class derived from either std::logic_error or
// std::runtime_error.
//
// Usage: OPM_DEFLOG_THROW(ExceptionClass, "Error message", DeferredLogger);
#define OPM_DEFLOG_THROW(Exception, message, deferred_logger)  \
    do {                                                       \
        std::string oss_ = std::string{"["} + __FILE__ + ":" + \
                           std::to_string(__LINE__) + "] " +   \
                           message;                            \
        deferred_logger.error(oss_);                           \
        throw Exception(oss_);                                 \
    } while (false)

// Macro to log a problem and throw an exception.
// Idenitical to OPM_DEFLOG_THROW() except for using
// the "problem" category instead of "error" for the
// log message. The Exception argument will typically
// be NumericalProblem.
//
// Usage: OPM_DEFLOG_PROBLEM(ExceptionClass, "Error message", DeferredLogger);
#define OPM_DEFLOG_PROBLEM(Exception, message, deferred_logger)  \
    do {                                                         \
        std::string oss_ = std::string{"["} + __FILE__ + ":" +   \
                           std::to_string(__LINE__) + "] " +     \
                           message;                              \
        deferred_logger.problem(oss_);                           \
        throw Exception(oss_);                                   \
    } while (false)


namespace {

void _throw(Opm::ExceptionType::ExcEnum exc_type,
            const std::string& message,
            Opm::Parallel::Communication comm)
{
    auto global_exc = comm.max(exc_type);

    switch (global_exc) {
    case Opm::ExceptionType::NONE:
        break;
    case Opm::ExceptionType::RUNTIME_ERROR:
        throw std::runtime_error(message);
        break;
    case Opm::ExceptionType::INVALID_ARGUMENT:
        throw std::invalid_argument(message);
        break;
    case Opm::ExceptionType::NUMERICAL_ISSUE:
        throw Opm::NumericalProblem(message);
        break;
    case Opm::ExceptionType::DEFAULT:
    case Opm::ExceptionType::LOGIC_ERROR:
        throw std::logic_error(message);
        break;
    default:
        throw std::logic_error(message);
    }
}

} // anonymous namespace



inline void checkForExceptionsAndThrow(Opm::ExceptionType::ExcEnum exc_type,
                                       const std::string& message,
                                       Opm::Parallel::Communication comm)
{
    _throw(exc_type, message, comm);
}

inline void logAndCheckForExceptionsAndThrow(Opm::DeferredLogger& deferred_logger,
                                             Opm::ExceptionType::ExcEnum exc_type,
                                             const std::string& message,
                                             const bool terminal_output,
                                             Opm::Parallel::Communication comm)
{
    Opm::DeferredLogger global_deferredLogger = gatherDeferredLogger(deferred_logger, comm);

    if (terminal_output) {
        global_deferredLogger.logMessages();
    }
    // Now that all messages have been logged, they are automatically
    // cleared from the global logger, but we must also clear them
    // from the local logger.
    deferred_logger.clearMessages();
    _throw(exc_type, message, comm);
}


/// \brief Macro to setup the try of a parallel try-catch
///
/// Use OPM_END_PARALLEL_TRY_CATCH or OPM_END_PARALLEL_TRY_CATCH_LOG
/// fot the catch part.
#define OPM_BEGIN_PARALLEL_TRY_CATCH()    \
std::string obptc_exc_msg;                \
auto obptc_exc_type = Opm::ExceptionType::NONE; \
try {

/// \brief Inserts catch classes for the parallel try-catch
///
/// There is a clause that will catch anything
#define OPM_PARALLEL_CATCH_CLAUSE(obptc_exc_type,          \
                                  obptc_exc_msg)           \
catch (const Opm::NumericalProblem& e){                    \
    obptc_exc_type = Opm::ExceptionType::NUMERICAL_ISSUE;  \
    obptc_exc_msg = e.what();                              \
} catch (const std::runtime_error& e) {                    \
    obptc_exc_type = Opm::ExceptionType::RUNTIME_ERROR;    \
    obptc_exc_msg = e.what();                              \
} catch (const std::invalid_argument& e) {                 \
    obptc_exc_type = Opm::ExceptionType::INVALID_ARGUMENT; \
    obptc_exc_msg = e.what();                              \
} catch (const std::logic_error& e) {                      \
    obptc_exc_type = Opm::ExceptionType::LOGIC_ERROR;      \
    obptc_exc_msg = e.what();                              \
} catch (const std::exception& e) {                        \
    obptc_exc_type = Opm::ExceptionType::DEFAULT;          \
    obptc_exc_msg = e.what();                              \
} catch (...) {                                            \
    obptc_exc_type = Opm::ExceptionType::DEFAULT;          \
    obptc_exc_msg = "Unknown exception was thrown";        \
}

/// \brief Catch exception and throw in a parallel try-catch clause
///
/// Assumes that OPM_BEGIN_PARALLEL_TRY_CATCH() was called to initiate
/// the try-catch clause
#define OPM_END_PARALLEL_TRY_CATCH(prefix, comm)         \
}                                                        \
OPM_PARALLEL_CATCH_CLAUSE(obptc_exc_type, obptc_exc_msg);\
checkForExceptionsAndThrow(obptc_exc_type,               \
                           prefix + obptc_exc_msg, comm);                        

/// \brief Catch exception, log, and throw in a parallel try-catch clause
///
/// Assumes that OPM_BEGIN_PARALLEL_TRY_CATCH() was called to initiate
/// the try-catch clause
#define OPM_END_PARALLEL_TRY_CATCH_LOG(obptc_logger, obptc_prefix, obptc_output, comm)\
}                                                              \
OPM_PARALLEL_CATCH_CLAUSE(obptc_exc_type, obptc_exc_msg);      \
logAndCheckForExceptionsAndThrow(obptc_logger, obptc_exc_type, \
    obptc_prefix + obptc_exc_msg, obptc_output, comm);
#endif // OPM_DEFERREDLOGGINGERRORHELPERS_HPP
