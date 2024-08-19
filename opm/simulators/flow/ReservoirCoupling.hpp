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
    SlaveProcessTermination,
    SlaveNextReportDate,
    SlaveNextTimeStep,
};

// Custom deleter for MPI_Comm
struct MPI_Comm_Deleter {
    void operator()(MPI_Comm* comm) const {
        if (*comm != MPI_COMM_NULL) {
            MPI_Comm_free(comm);
        }
        delete comm;
    }
};

using MPI_Comm_Ptr = std::unique_ptr<MPI_Comm, MPI_Comm_Deleter>;

// This class represents a time point.
// It is currently used to represent an epoch time (a double value in seconds since the epoch),
// or an elapsed time (a double value in seconds since the start of the simulation).
// To avoid numerical issues when adding or subtracting time points and then later comparing
// for equality with for example a given report date, we use a tolerance value.
class TimePoint {
private:
    double time;
    // TODO: Epoch values often lies in the range of [1e9,1e11], so a tolerance value of 1e-10
    //     might be a little too small. However, for elapsed time values, the range is often
    //     in the range of [0, 1e8], so a tolerance value of 1e-10 should be sufficient.
    // NOTE: 1 nano-second = 1e-9 seconds
    static constexpr double tol = 1e-10; // Tolerance value

public:
    TimePoint() : time(0.0) {}
    explicit TimePoint(double t) : time(t) {}
    TimePoint(const TimePoint& other) : time(other.time) {}

    // Assignment operator for double
    TimePoint& operator=(double t) {
        time = t;
        return *this;
    }

    // Copy assignment operator
    TimePoint& operator=(const TimePoint& other) {
        if (this != &other) {
            time = other.time;
        }
        return *this;
    }

    double getTime() const { return time; }

    // Equality operator
    bool operator==(const TimePoint& other) const {
        return std::abs(time - other.time) < tol;
    }

    // Inequality operator
    bool operator!=(const TimePoint& other) const {
        return !(*this == other);
    }

    // Less than operator
    bool operator<(const TimePoint& other) const {
        return (time < other.time) && !(*this == other);
    }

    // Comparison operator: double < TimePoint
    friend bool operator<(double lhs, const TimePoint& rhs) {
        return lhs < rhs.time;
    }

   // Comparison operator: TimePoint < double
    bool operator<(double rhs) const {
        return time < rhs;
    }

    // Less than or equal to operator
    bool operator<=(const TimePoint& other) const {
        return (time < other.time) || (*this == other);
    }

    // Comparison operator: double <= TimePoint
    friend bool operator<=(double lhs, const TimePoint& rhs) {
        return lhs <= rhs.time;
    }

    // Comparison operator: TimePoint <= double
    bool operator<=(double rhs) const {
        return time <= rhs;
    }

    // Greater than operator
    bool operator>(const TimePoint& other) const {
        return (time > other.time) && !(*this == other);
    }

    // Comparison operator: double > TimePoint
    friend bool operator>(double lhs, const TimePoint& rhs) {
        return lhs > rhs.time;
    }

    // Comparison operator: TimePoint > double
    bool operator>(double rhs) const {
        return time > rhs;
    }

    // Greater than or equal to operator
    bool operator>=(const TimePoint& other) const {
        return (time > other.time) || (*this == other);
    }

    // Comparison operator: TimePoint >= double
    bool operator>=(double rhs) const {
        return time >= rhs;
    }

    // Comparison operator: double >= TimePoint
    friend bool operator>=(double lhs, const TimePoint& rhs) {
        return lhs >= rhs.time;
    }

    // Addition operator: TimePoint + TimePoint (summing their times)
    TimePoint operator+(const TimePoint& other) const {
        return TimePoint(time + other.time);
    }

    // Addition operator: TimePoint + double
    TimePoint operator+(double delta) const {
        return TimePoint(time + delta);
    }

    // Friend addition operator: double + TimePoint
    friend TimePoint operator+(double lhs, const TimePoint& rhs) {
        return TimePoint(lhs + rhs.time);
    }

    // Overload += operator for adding a double
    TimePoint& operator+=(double delta) {
        time += delta;
        return *this;
    }

    // Subtraction operator: TimePoint - TimePoint (resulting in a new TimePoint)
    TimePoint operator-(const TimePoint& other) const {
        return TimePoint(time - other.time);
    }

    // Subtraction operator: TimePoint - double
    TimePoint operator-(double delta) const {
        return TimePoint(time - delta);
    }

    // Friend subtraction operator: double - TimePoint
    friend TimePoint operator-(double lhs, const TimePoint& rhs) {
        return TimePoint(lhs - rhs.time);
    }

    // Stream insertion operator for easy printing
    friend std::ostream& operator<<(std::ostream& os, const TimePoint& tp) {
        os << tp.time;
        return os;
    }
};

} // namespace ReservoirCoupling
} // namespace Opm

#endif // OPM_RESERVOIR_COUPLING_HPP
