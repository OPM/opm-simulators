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
#include <opm/simulators/utils/DeferredLogger.hpp>
#include <opm/simulators/utils/ParallelCommunication.hpp>
#include <opm/input/eclipse/Schedule/Group/Group.hpp>
#include <opm/input/eclipse/Schedule/Group/GuideRate.hpp>
#include <opm/input/eclipse/Units/Units.hpp>

#include <dune/common/parallel/mpitraits.hh>

#include <fmt/format.h>

#include <mpi.h>
#include <cmath>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <vector>

namespace Opm {
namespace ReservoirCoupling {

class Logger {
public:
    Logger() = delete;  // No default constructor - must have comm
    explicit Logger(const Parallel::Communication& comm) : comm_(comm) {}

    void clearDeferredLogger() { deferred_logger_ = nullptr; }
    void debug(const std::string &msg) const;
    DeferredLogger& deferredLogger() { return *deferred_logger_; }
    DeferredLogger& deferredLogger() const { return *deferred_logger_; }
    bool haveDeferredLogger() const { return deferred_logger_ != nullptr; }
    void info(const std::string &msg) const;
    void warning(const std::string &msg) const;
    void setDeferredLogger(DeferredLogger *deferred_logger) { deferred_logger_ = deferred_logger; }

private:
    const Parallel::Communication& comm_;
    DeferredLogger *deferred_logger_ = nullptr;
};

/// @brief Guard for managing DeferredLogger lifecycle in ReservoirCoupling
///
/// @details This class ensures that a DeferredLogger pointer is properly set and cleared
/// in the ReservoirCoupling logger. It follows the RAII (Resource Acquisition Is Initialization)
/// pattern to prevent dangling pointer issues when local DeferredLogger objects go out of scope.
///
/// - Uses a pointer instead of a reference to enable move semantics
/// - Implements move constructor/assignment to allow returning from functions and use in std::optional
/// - Deletes copy constructor/assignment to ensure single ownership
/// - The moved-from object is "nullified" so only the final owner clears the logger
///
/// Lifecycle guarantees:
/// - Constructor: Sets the DeferredLogger pointer in the ReservoirCoupling logger
/// - Destructor: Clears the pointer (only if not moved-from)
/// - Move: Transfers ownership; moved-from object becomes inactive (logger_ = nullptr)
///
/// @param logger The ReservoirCoupling logger (either master or slave)
/// @param deferred_logger The DeferredLogger to bind to the ReservoirCoupling logger
///
/// Example usage:
/// @code
///   DeferredLogger local_logger;
///   {
///       ReservoirCoupling::ScopedLoggerGuard guard{rc_logger, &local_logger};
///       // logger is active here
///   } // guard destructor automatically clears the logger
/// @endcode
class ScopedLoggerGuard {
public:
    explicit ScopedLoggerGuard(Logger& logger, DeferredLogger* deferred_logger)
        : logger_(&logger)
    {
        logger_->setDeferredLogger(deferred_logger);
    }

    ~ScopedLoggerGuard() {
        // Only clear if we still own the logger (not moved-from)
        if (logger_) {
            logger_->clearDeferredLogger();
        }
    }

    // Prevent copying to ensure single ownership
    ScopedLoggerGuard(const ScopedLoggerGuard&) = delete;
    ScopedLoggerGuard& operator=(const ScopedLoggerGuard&) = delete;

    // Enable moving - required for returning from functions and std::optional
    ScopedLoggerGuard(ScopedLoggerGuard&& other) noexcept
        : logger_(other.logger_)
    {
        // Transfer ownership: moved-from object becomes inactive and won't clear the logger
        other.logger_ = nullptr;
    }

    ScopedLoggerGuard& operator=(ScopedLoggerGuard&& other) noexcept {
        if (this != &other) {
            // Clean up current logger before taking ownership of new one
            if (logger_) {
                logger_->clearDeferredLogger();
            }
            // Transfer ownership from other
            logger_ = other.logger_;
            other.logger_ = nullptr;
        }
        return *this;
    }

private:
    // Use pointer instead of reference to enable move semantics
    // (references cannot be reassigned, which is required for move operations)
    Logger* logger_{nullptr};
};

enum class MessageTag : int {
    InjectionGroupTargets,
    MasterGroupNames,
    MasterGroupNamesSize,
    MasterStartOfReportStep,
    NumSlaveGroupTargets,
    ProductionGroupTargets,
    SlaveActivationDate,
    SlaveActivationHandshake,
    SlaveInjectionData,
    SlaveProcessTermination,
    SlaveName,
    SlaveNameSize,
    SlaveNextReportDate,
    SlaveNextTimeStep,
    SlaveProductionData,
    SlaveSimulationStartDate,
    SlaveStartOfReportStep,
};

/// @brief Phase indices for reservoir coupling, we currently only support black-oil phases
/// (oil, gas, and water).
enum class Phase : std::size_t {
    Oil = 0,  // Matches Opm::Phase::OIL
    Gas,      // Matches Opm::Phase::GAS
    Water,    // Matches Opm::Phase::WATER
    Count
};

/// @brief Selects which kind of rate to retrieve from slave group data.
///
/// Replaces multiple bool parameters (res_rates, is_injection, network) with a single
/// enum that names the 5 valid rate categories directly. Only these 5 combinations
/// of the bools are meaningful; the enum prevents callers from constructing invalid states.
enum class RateKind {
    InjectionSurface,
    InjectionReservoir,
    ProductionSurface,
    ProductionNetworkSurface,
    ProductionReservoir
};

template <class Scalar>
struct InjectionRates {
    InjectionRates() = default;

    std::array<Scalar, static_cast<std::size_t>(Phase::Count)> rate{};
    [[nodiscard]] Scalar& operator[](Phase p)       noexcept { return rate[static_cast<std::size_t>(p)]; }
    [[nodiscard]] Scalar  operator[](Phase p) const noexcept { return rate[static_cast<std::size_t>(p)]; }
};


// Used to communicate potentials for oil, gas, and water rates between slave and master processes
template <class Scalar>
struct Potentials {
    std::array<Scalar, static_cast<std::size_t>(Phase::Count)> rate{};

    [[nodiscard]] Scalar& operator[](Phase p)       noexcept { return rate[static_cast<std::size_t>(p)]; }
    [[nodiscard]] Scalar  operator[](Phase p) const noexcept { return rate[static_cast<std::size_t>(p)]; }
};

template <class Scalar>
struct ProductionRates {
    ProductionRates() = default;

    explicit ProductionRates(const GuideRate::RateVector& rate_vector)
        : rate{static_cast<Scalar>(rate_vector.oil_rat),
               static_cast<Scalar>(rate_vector.gas_rat),
               static_cast<Scalar>(rate_vector.wat_rat)}
    {}

    std::array<Scalar, static_cast<std::size_t>(Phase::Count)> rate{};
    [[nodiscard]] Scalar& operator[](Phase p)       noexcept { return rate[static_cast<std::size_t>(p)]; }
    [[nodiscard]] Scalar  operator[](Phase p) const noexcept { return rate[static_cast<std::size_t>(p)]; }
};

// Slave group production data sent to the corresponding master group for target calculation.
template <class Scalar>
struct SlaveGroupProductionData {
    // Group production potentials are used by the master group for guiderate calculations
    Potentials<Scalar> potentials;
    // Production rates are used by the master group in guiderate calculations
    // when converting the guide rate target to the phase of the master group.
    ProductionRates<Scalar> surface_rates;  // Surface production rates by phase (network=false)
    // Network surface rates - computed with network=true, meaning efficiency factors
    // are 1.0 for groups/wells with GEFAC/WEFAC item 3 = "NO"
    ProductionRates<Scalar> network_surface_rates;  // Surface rates for network calculations
    // Individual phase reservoir production rates - needed when master's parent group
    // has RESV control mode, so the conversion uses slave's PVT properties
    ProductionRates<Scalar> reservoir_rates;  // Reservoir production rates by phase
    Scalar voidage_rate{0.0};               // Reservoir voidage replacement rate
    Scalar gas_reinjection_rate{0.0};       // Reinjection (surface) rate for the gas phase
};

// Slave group injection data sent to the corresponding master group for target calculation.
template <class Scalar>
struct SlaveGroupInjectionData {
    InjectionRates<Scalar> surface_rates;    // Surface injection rates by phase
    InjectionRates<Scalar> reservoir_rates;  // Reservoir injection rates by phase
};

template <class Scalar>
struct InjectionGroupTarget {
    // To save memory and avoid varying size of the struct when serializing
    // and deserializing the group name, we use an index instead of the full name.
    std::size_t group_name_idx;   // Index of group name in the master group names vector
    Scalar target;                // Target rate for the group
    Group::InjectionCMode cmode;  // Control mode for the group
    Phase phase;                  // Phase the target applies to
};

template <class Scalar>
struct ProductionGroupTarget {
    // To save memory and avoid varying size of the struct when serializing
    // and deserializing the group name, we use an index instead of the full name.
    std::size_t group_name_idx;   // Index of group name in the master group names vector
    Scalar target;                // Target rate for the group
    Group::ProductionCMode cmode;  // Control mode for the group
};

// Helper functions
Phase convertPhaseToReservoirCouplingPhase(::Opm::Phase phase);
void customErrorHandler_(MPI_Comm* comm, int* err, const std::string &msg);
void customErrorHandlerSlave_(MPI_Comm* comm, int* err, ...);
void customErrorHandlerMaster_(MPI_Comm* comm, int* err, ...);
void setErrhandler(MPI_Comm comm, bool is_master);
std::pair<std::vector<char>, std::size_t> serializeStrings(const std::vector<std::string>& data);

/// \brief Format seconds as a human-readable string showing both seconds and days.
/// \param seconds The time value in seconds.
/// \return A string like "864000s (10.00 days)".
inline std::string formatDays(double seconds) {
    double days = seconds / unit::day;
    return fmt::format(fmt::runtime("{:.0f}s ({:.2f} days)"), seconds, days);
}

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

    /// \brief Determines if \a a is less than \a b within the specified tolerance.
    ///
    /// \param a First double value.
    /// \param b Second double value.
    /// \return True if \a a is less than \a b.
    static bool compare_lt(double a, double b);

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
