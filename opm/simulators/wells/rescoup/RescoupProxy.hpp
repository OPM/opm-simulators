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

#ifndef OPM_RESCOUP_PROXY_HPP
#define OPM_RESCOUP_PROXY_HPP

#include <opm/simulators/flow/rescoup/ReservoirCouplingEnabled.hpp>
#ifdef RESERVOIR_COUPLING_ENABLED
#include <opm/simulators/flow/rescoup/ReservoirCoupling.hpp>
#include <opm/simulators/flow/rescoup/ReservoirCouplingMaster.hpp>
#include <opm/simulators/flow/rescoup/ReservoirCouplingSlave.hpp>
#else
#include <stdexcept>
#include <string>
#endif

namespace Opm {

#ifndef RESERVOIR_COUPLING_ENABLED
// Forward declarations for non-MPI builds
template <class Scalar> class ReservoirCouplingMaster;
template <class Scalar> class ReservoirCouplingSlave;
#endif

namespace ReservoirCoupling {

/// @brief Thin proxy for reservoir coupling master/slave pointers.
///
/// This class encapsulates pointers to ReservoirCouplingMaster and ReservoirCouplingSlave,
/// providing mode queries that work regardless of MPI availability. This eliminates the need
/// for #ifdef RESERVOIR_COUPLING_ENABLED guards in consumer code.
///
/// The proxy does not own the master/slave objects - they are owned by SimulatorFullyImplicitBlackoil.
/// Pointers are set after construction via setMaster()/setSlave().
///
/// Master and slave modes are mutually exclusive - setting one clears the other.
template <class Scalar>
class Proxy {
public:
    Proxy() = default;

    // Copyable and movable (just holds non-owning pointers)
    Proxy(const Proxy&) = default;
    Proxy& operator=(const Proxy&) = default;
    Proxy(Proxy&&) noexcept = default;
    Proxy& operator=(Proxy&&) noexcept = default;

#ifdef RESERVOIR_COUPLING_ENABLED
    // === Mode Queries ===

    /// @brief Check if reservoir coupling is enabled (master or slave mode)
    bool isEnabled() const noexcept { return master_ || slave_; }

    /// @brief Check if this process is the reservoir coupling master
    bool isMaster() const noexcept { return master_ != nullptr; }

    /// @brief Check if this process is a reservoir coupling slave
    bool isSlave() const noexcept { return slave_ != nullptr; }

    // === Setters (called during init, after construction) ===

    /// @brief Set the master pointer (clears slave pointer)
    void setMaster(ReservoirCouplingMaster<Scalar>* master) {
        master_ = master;
        slave_ = nullptr;
    }

    /// @brief Set the slave pointer (clears master pointer)
    void setSlave(ReservoirCouplingSlave<Scalar>* slave) {
        slave_ = slave;
        master_ = nullptr;
    }

    // === Pointer Access ===

    /// @brief Get raw pointer to master (nullptr if not in master mode)
    ReservoirCouplingMaster<Scalar>* masterPtr() noexcept { return master_; }
    const ReservoirCouplingMaster<Scalar>* masterPtr() const noexcept { return master_; }

    /// @brief Get raw pointer to slave (nullptr if not in slave mode)
    ReservoirCouplingSlave<Scalar>* slavePtr() noexcept { return slave_; }
    const ReservoirCouplingSlave<Scalar>* slavePtr() const noexcept { return slave_; }

    // === Reference Access (caller must ensure correct mode) ===

    /// @brief Get reference to master (undefined behavior if not in master mode)
    ReservoirCouplingMaster<Scalar>& master() { return *master_; }
    const ReservoirCouplingMaster<Scalar>& master() const { return *master_; }

    /// @brief Get reference to slave (undefined behavior if not in slave mode)
    ReservoirCouplingSlave<Scalar>& slave() { return *slave_; }
    const ReservoirCouplingSlave<Scalar>& slave() const { return *slave_; }

    // === Facade Methods (safe to call regardless of mode) ===

    /// @brief Check if a group is a master group (controls a slave reservoir)
    /// @param group_name The name of the group to check
    /// @return true if in master mode and the group is a master group, false otherwise
    bool isMasterGroup(const std::string& group_name) const {
        return master_ && master_->isMasterGroup(group_name);
    }

    /// @brief Check if a group is a slave group
    /// @param group_name The name of the group to check
    /// @return true if in slave mode and the group is a slave group, false otherwise
    bool isSlaveGroup(const std::string& group_name) const {
        return slave_ && slave_->isSlaveGroup(group_name);
    }

private:
    ReservoirCouplingMaster<Scalar>* master_ = nullptr;
    ReservoirCouplingSlave<Scalar>* slave_ = nullptr;

#else // !RESERVOIR_COUPLING_ENABLED

    // === Mode Queries (always false in non-MPI builds) ===

    bool isEnabled() const noexcept { return false; }
    bool isMaster() const noexcept { return false; }
    bool isSlave() const noexcept { return false; }

    // Stub implementations for non-MPI builds.
    // These should never be called since isMaster()/isSlave() always return false.

    void setMaster(ReservoirCouplingMaster<Scalar>*) {
        throw std::logic_error("ReservoirCoupling::Proxy::setMaster() called in non-MPI build");
    }

    void setSlave(ReservoirCouplingSlave<Scalar>*) {
        throw std::logic_error("ReservoirCoupling::Proxy::setSlave() called in non-MPI build");
    }

    ReservoirCouplingMaster<Scalar>* masterPtr() noexcept { return nullptr; }
    const ReservoirCouplingMaster<Scalar>* masterPtr() const noexcept { return nullptr; }

    ReservoirCouplingSlave<Scalar>* slavePtr() noexcept { return nullptr; }
    const ReservoirCouplingSlave<Scalar>* slavePtr() const noexcept { return nullptr; }

    ReservoirCouplingMaster<Scalar>& master() {
        throw std::logic_error("ReservoirCoupling::Proxy::master() called in non-MPI build");
    }
    const ReservoirCouplingMaster<Scalar>& master() const {
        throw std::logic_error("ReservoirCoupling::Proxy::master() called in non-MPI build");
    }

    ReservoirCouplingSlave<Scalar>& slave() {
        throw std::logic_error("ReservoirCoupling::Proxy::slave() called in non-MPI build");
    }
    const ReservoirCouplingSlave<Scalar>& slave() const {
        throw std::logic_error("ReservoirCoupling::Proxy::slave() called in non-MPI build");
    }

    // === Facade Methods (always return false/no-op in non-MPI builds) ===

    bool isMasterGroup(const std::string& /*group_name*/) const noexcept {
        return false;
    }

    bool isSlaveGroup(const std::string& /*group_name*/) const noexcept {
        return false;
    }
#endif // !RESERVOIR_COUPLING_ENABLED
};

} // namespace ReservoirCoupling
} // namespace Opm

#endif // OPM_RESCOUP_PROXY_HPP
