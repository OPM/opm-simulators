/*
  Copyright 2024 Equinor AS

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

#ifndef OPM_RESERVOIR_COUPLING_HPP
#define OPM_RESERVOIR_COUPLING_HPP
#include <mpi.h>
#include <cmath>
#include <iostream>
#include <memory>

namespace Opm {
namespace ReservoirCoupling {

enum class MessageTag : int {
    SlaveSimulationStartDate,
    SlaveActivationDate,
    SlaveProcessTermination,
    SlaveNextReportDate,
    SlaveNextTimeStep,
    MasterGroupNames,
    MasterGroupNamesSize,
};

// Helper functions
void custom_error_handler_(MPI_Comm* comm, int* err, const std::string &msg);
void setErrhandler(MPI_Comm comm, bool is_master);

// Utility class for comparing double values representing epoch dates (seconds since
// unix epoch) or elapsed time (seconds since the start of the simulation).
// NOTE: It is important that when comparing against start of a report step or similar, that
//   that we do not miss these due to numerical issues. This is because communication between
//   master and slave processes are based on these points in time.
// NOTE: Epoch values in this century (2000-2100) lies in the range of [1e9,4e9], and a double variable cannot
//  represent such large values with high precision. For example, the date 01-01-2020 is equal
//  to 1.5778368e9 seconds and adding 1e-7 seconds to this value will not change the value.
//  So microseconds (1e-6) is approximately the smallest time unit we can represent for such a number.
// NOTE: Report steps seems to have a maximum resolution of whole seconds, see stepLength() in
//   Schedule.cpp in opm-common, which returns the step length in seconds.
struct Seconds {
    static constexpr double abstol = 1e-15;
    static constexpr double reltol = 1e-15;
    // We will will use the following expression to determine if two values a and b are equal:
    //   |a - b| <= tol = abstol + reltol * max(|a|, |b|)
    // For example, assume abstol = reltol = 1e-15, then the following holds:
    // - If |a| and |b| are below 1, then the absolute tolerance applies.
    // - If a and b are above 1, then the relative tolerance applies.
    // For example, for dates in the range 01-01-2000 to 01-01-2100, epoch values will be in the range
    //  [1e9, 4e9]. And we have 1e-15 * 1e9 = 1e-6, so numbers differing below one microsecond will
    // be considered equal.
    // NOTE: The above is not true for numbers close to zero, but we do not expect to compare such numbers.
    static bool compare_eq(double a, double b);
    static bool compare_gt(double a, double b);
    static bool compare_gt_or_eq(double a, double b);
    static bool compare_lt_or_eq(double a, double b);
};

} // namespace ReservoirCoupling
} // namespace Opm

#endif // OPM_RESERVOIR_COUPLING_HPP
