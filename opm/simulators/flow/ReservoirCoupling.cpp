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

#include <opm/simulators/flow/ReservoirCoupling.hpp>
#include <opm/simulators/utils/ParallelCommunication.hpp>
#include <opm/common/OpmLog/OpmLog.hpp>

#include <fmt/format.h>

namespace Opm {
namespace ReservoirCoupling {

// free functions alphabetically
// -----------------------------

void custom_error_handler_(MPI_Comm* comm, int* err, const std::string &msg)
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

void custom_error_handler_slave_(MPI_Comm* comm, int* err, ...)
{
    custom_error_handler_(comm, err, "slave");
}

void custom_error_handler_master_(MPI_Comm* comm, int* err, ...)
{
    custom_error_handler_(comm, err, "master");
}

std::pair<std::vector<char>, std::size_t> serializeStrings(const std::vector<std::string>& data)
{
    std::size_t total_size = 0;
    std::vector<char> serialized_data;
    for (const auto& str: data) {
        std::copy(str.begin(), str.end(), std::back_inserter(serialized_data));
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
        MPI_Comm_create_errhandler(custom_error_handler_master_, &errhandler);
    }
    else {
        MPI_Comm_create_errhandler(custom_error_handler_slave_, &errhandler);
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

void Logger::info(const std::string &msg) const {
    if (haveDeferredLogger()) {
        this->deferred_logger_->info(msg);
    } else {
        OpmLog::info(msg);
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


} // namespace ReservoirCoupling
} // namespace Opm
