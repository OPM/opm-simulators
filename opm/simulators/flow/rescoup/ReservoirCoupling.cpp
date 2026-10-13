/*
  Copyright 2024 Equinor ASA

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

#include <config.h>

#include <opm/simulators/flow/rescoup/ReservoirCoupling.hpp>
#include <opm/simulators/utils/ParallelCommunication.hpp>
#include <opm/common/OpmLog/OpmLog.hpp>
#include <opm/input/eclipse/Parser/ParserKeywords/G.hpp>

#include <fmt/format.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <vector>

namespace Opm {
namespace ReservoirCoupling {

// free functions alphabetically
// -----------------------------

Phase convertPhaseToReservoirCouplingPhase(::Opm::Phase phase)
{
    switch (phase) {
    case ::Opm::Phase::OIL:
        return Phase::Oil;
    case ::Opm::Phase::GAS:
        return Phase::Gas;
    case ::Opm::Phase::WATER:
        return Phase::Water;
    default:
        throw std::invalid_argument{"Unsupported Opm::Phase value for reservoir coupling."};
    }
}

::Opm::Phase convertToOpmPhase(const Phase phase)
{
    switch (phase) {
    case Phase::Oil:
        return ::Opm::Phase::OIL;
    case Phase::Gas:
        return ::Opm::Phase::GAS;
    case Phase::Water:
        return ::Opm::Phase::WATER;
    default:
        throw std::invalid_argument{"Unsupported ReservoirCoupling::Phase value for conversion to Opm::Phase."};
    }
}

void customErrorHandler_(MPI_Comm* comm, int* err, const std::string &msg)
{
    // It can be useful to have a custom error handler for debugging purposes.
    // If not, MPI will use the default error handler which aborts the program, and
    // it can be difficult to determine exactly where the error occurred. With a custom
    // error handler you can at least set a breakpoint here to get a backtrace.
    int rank;
    MPI_Comm_rank(*comm, &rank);
    char err_string[MPI_MAX_ERROR_STRING];
    int len;
    MPI_Error_string(*err, err_string, &len);
    const std::string error_msg = fmt::format(
        "Reservoir coupling MPI error handler {} rank {}: {}", msg, rank, err_string);
    // NOTE: The output to terminal vie stderr or stdout has been redirected to files for
    //   the slaves, see Main.cpp. So the following output will not be visible in the terminal
    //   if we are called from a slave process.
    // std::cerr << error_msg << std::endl;
    OpmLog::error(error_msg);  // Output to log file, also shows the message on stdout for master.
    MPI_Abort(*comm, *err);
}

void customErrorHandlerSlave_(MPI_Comm* comm, int* err, ...)
{
    customErrorHandler_(comm, err, "slave");
}

void customErrorHandlerMaster_(MPI_Comm* comm, int* err, ...)
{
    customErrorHandler_(comm, err, "master");
}

std::pair<std::vector<char>, std::size_t> serializeStrings(const std::vector<std::string>& data)
{
    std::size_t total_size = 0;
    std::vector<char> serialized_data;
    for (const auto& str: data) {
        std::ranges::copy(str, std::back_inserter(serialized_data));
        serialized_data.push_back('\0');
        total_size += str.size() + 1;
    }
    return {serialized_data, total_size};
}

void setErrhandler(MPI_Comm comm, bool is_master)
{
    MPI_Errhandler errhandler;
    // NOTE: Lambdas with captures cannot be used as C function pointers, also
    //   converting lambdas that use ellipsis "..." as last argument to a C function pointer
    //   is not currently possible, so we need to use static functions instead.
    if (is_master) {
        MPI_Comm_create_errhandler(customErrorHandlerMaster_, &errhandler);
    }
    else {
        MPI_Comm_create_errhandler(customErrorHandlerSlave_, &errhandler);
    }
    // NOTE: The errhandler is a handle (an integer) that is associated with the communicator
    //   that is why we pass this by value below. And it is safe to free the errhandler after it has
    //   been associated with the communicator.
    MPI_Comm_set_errhandler(comm, errhandler);
    // Mark the error handler for deletion. According to the documentation: "The error handler will
    //  be deallocated after all the objects associated with it (communicator, window, or file) have
    // been deallocated." So the error handler will still be in effect until the communicator is
    // deallocated.
    MPI_Errhandler_free(&errhandler);
}

// Logger class alphabetically
// ---------------------------

void Logger::debug(const std::string &msg) const {
    forward_(msg,
        [](DeferredLogger &dl, const std::string &m) { dl.debug(m); },
        [](const std::string &m) { OpmLog::debug(m); });
}

void Logger::error(const std::string &msg) const {
    forward_(msg,
        [](DeferredLogger &dl, const std::string &m) { dl.error(m); },
        [](const std::string &m) { OpmLog::error(m); });
}

void Logger::info(const std::string &msg) const {
    forward_(msg,
        [](DeferredLogger &dl, const std::string &m) { dl.info(m); },
        [](const std::string &m) { OpmLog::info(m); });
}

void Logger::warning(const std::string &msg) const {
    forward_(msg,
        [](DeferredLogger &dl, const std::string &m) { dl.warning(m); },
        [](const std::string &m) { OpmLog::warning(m); });
}

// Private methods of Logger
// -------------------------

// Forward a log message to either DeferredLogger (all ranks) or OpmLog (rank 0 only).
// DeferredLogger being null is the expected normal state - it's only temporarily set
// during beginTimeStep() via ScopedLoggerGuard. Messages logged here are identical on
// all ranks, so rank-0-only OpmLog fallback preserves all information.
template<typename DeferredFn, typename OpmLogFn>
void Logger::forward_(const std::string &msg, DeferredFn deferred_fn, OpmLogFn opmlog_fn) const {
    if (haveDeferredLogger()) {
        deferred_fn(*this->deferred_logger_, msg);
    } else {
        if (comm_.rank() == 0) {
            opmlog_fn(msg);
        }
    }
}


// Seconds class alphabetically
// ----------------------------

bool Seconds::compare_eq(double a, double b)
{
    // Are a and b equal?
    return std::abs(a - b) < std::max(abstol, reltol * std::max(std::abs(a), std::abs(b)));
}

bool Seconds::compare_gt_or_eq(double a, double b)
{
    // Is a greater than or equal to b?
    if (compare_eq(a, b)) {
        return true;
    }
    return a > b;
}

bool Seconds::compare_lt(double a, double b)
{
    // Is a less than b?
    return !compare_eq(a, b) && a < b;
}

bool Seconds::compare_gt(double a, double b)
{
    // Is a greater than b?
    return !compare_eq(a, b) && a > b;
}

bool Seconds::compare_lt_or_eq(double a, double b)
{
    // Is a less than or equal to b?
    if (compare_eq(a, b)) {
        return true;
    }
    return a < b;
}

template <class Scalar>
FlowChangeLimit<Scalar>
flowChangeTimeStepLimit(const Scalar limiting_fraction,
                        const Scalar step_length,
                        const Scalar min_step,
                        const std::vector<Scalar>& start_flows,
                        const std::vector<Scalar>& end_flows)
{
    constexpr Scalar infinity = std::numeric_limits<Scalar>::infinity();
    FlowChangeLimit<Scalar> result{infinity, 0};
    constexpr Scalar default_fraction = ParserKeywords::GRUPMAST::LIMITING_FRACTION::defaultValue;
    if (limiting_fraction >= default_fraction || !(step_length > 0)) {
        return result;
    }
    for (std::size_t i = 0; i < start_flows.size(); ++i) {
        const Scalar before = start_flows[i];
        const Scalar after = end_flows[i];
        if (before == after) {
            continue;
        }
        // A flow starting from zero is an infinite fractional change
        const Scalar change = (before == 0) ? infinity : std::abs(after - before) / std::abs(before);
        result.max_change = std::max(result.max_change, change);
    }
    if (result.max_change > 0) {
        // The master and the slaves add up the sync step lengths independently and compare the
        // sums with a tolerance of a few ulps (see Seconds). Whole seconds add up exactly, so
        // the limit is rounded to those, as the report step dates are.
        result.time_step_limit = std::max(std::round(step_length * limiting_fraction / result.max_change),
                                          std::max<Scalar>(1, std::ceil(min_step)));
    }
    return result;
}

template struct InjectionGroupTarget<double>;
template struct ProductionGroupConstraints<double>;
template FlowChangeLimit<double> flowChangeTimeStepLimit(
    double, double, double, const std::vector<double>&, const std::vector<double>&);

#if FLOW_INSTANTIATE_FLOAT
template struct InjectionGroupTarget<float>;
template struct ProductionGroupConstraints<float>;
template FlowChangeLimit<float> flowChangeTimeStepLimit(
    float, float, float, const std::vector<float>&, const std::vector<float>&);
#endif

} // namespace ReservoirCoupling
} // namespace Opm
