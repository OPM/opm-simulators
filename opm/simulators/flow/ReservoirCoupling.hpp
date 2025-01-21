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

/// \brief Utility class for comparing double values representing epoch dates or elapsed time.
///
/// This class is used to compare double values that represent:
/// - Epoch dates (seconds since the Unix epoch)
/// - Elapsed time (seconds since the start of the simulation)
///
/// \note When comparing against the start of a report step or similar, it is important not to miss
/// these points due to numerical issues. Communication between master and slave processes is based
/// on these specific points in time.
///
/// \note Epoch values in this century (2000-2100) fall within the range [1e9, 4e9]. A double variable
/// cannot represent such large values with high precision. For example, the date 01-01-2020 corresponds
/// to 1.5778368e9 seconds, and adding 1e-7 seconds to this value does not change it. Microseconds (1e-6)
/// are approximately the smallest time units that can be represented accurately for such numbers.
///
/// \note Report steps appear to have a maximum resolution of whole seconds. See the `stepLength()`
/// function in `Schedule.cpp` in the `opm-common` module, which returns the step length in seconds.
struct Seconds {
    /// \brief Absolute tolerance used for comparisons.
    static constexpr double abstol = 1e-15;

    /// \brief Relative tolerance used for comparisons.
    static constexpr double reltol = 1e-15;

    /// \brief Determines if two double values are equal within a specified tolerance.
    ///
    /// Two values \a a and \a b are considered equal if:
    /// \f[ |a - b| \leq \text{tol} = \text{abstol} + \text{reltol} \times \max(|a|, |b|) \f]
    ///
    /// For example, if \a abstol = \a reltol = 1e-15:
    /// - If \a |a| and \a |b| are below 1, the absolute tolerance applies.
    /// - If \a a and \a b are above 1, the relative tolerance applies.
    ///
    /// \note For epoch dates between 01-01-2000 and 01-01-2100, epoch values are in the range [1e9, 4e9].
    /// Therefore, \f$ 10^{-15} \times 10^9 = 10^{-6} \f$, meaning differences smaller than one microsecond
    /// will be considered equal.
    ///
    /// \note This approach is not accurate for numbers close to zero, but such comparisons are not expected.
    ///
    /// \param a First double value.
    /// \param b Second double value.
    /// \return True if \a a and \a b are considered equal within the tolerance.
    static bool compare_eq(double a, double b);

   /// \brief Determines if \a a is greater than \a b within the specified tolerance.
    ///
    /// \param a First double value.
    /// \param b Second double value.
    /// \return True if \a a is greater than \a b.
    static bool compare_gt(double a, double b);

   /// \brief Determines if \a a is greater than \a b within the specified tolerance.
    ///
    /// \param a First double value.
    /// \param b Second double value.
    /// \return True if \a a is greater than \a b.
    static bool compare_gt_or_eq(double a, double b);

    /// \brief Determines if \a a is less than or equal to \a b within the specified tolerance.
    ///
    /// \param a First double value.
    /// \param b Second double value.
    /// \return True if \a a is less than or equal to \a b.
    static bool compare_lt_or_eq(double a, double b);
};

} // namespace ReservoirCoupling
} // namespace Opm

#endif // OPM_RESERVOIR_COUPLING_HPP
