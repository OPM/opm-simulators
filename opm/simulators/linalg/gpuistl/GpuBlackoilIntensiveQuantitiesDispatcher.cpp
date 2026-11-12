/*
  Copyright 2026 Equinor ASA

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
/*!
 * \file
 *
 * The host side of \c GpuBlackoilIntensiveQuantitiesDispatcher: it chooses the
 * configuration from the PVT approaches of the CPU fluid system and forwards to
 * that configuration's \c Impl. The configurations themselves are in the
 * \c GpuBlackoilIntensiveQuantitiesDispatcher<Configuration>.cu files.
 */
#include <config.h>

#include <opm/common/ErrorMacros.hpp>
#include <opm/common/OpmLog/OpmLog.hpp>

#include <opm/material/fluidsystems/blackoilpvt/BlackOilPvtApproaches.hpp>

#include <opm/simulators/flow/FlowGasWaterEnergyTypeTag.hpp>
#include <opm/simulators/flow/FlowProblemBlackoilProperties.hpp>
#include <opm/simulators/flow/FlowProblemNohystTypeTag.hpp>
// After the TypeTags: the header uses the properties without declaring them.
#include <opm/simulators/linalg/gpuistl/GpuBlackoilIntensiveQuantitiesDispatcher.hpp>

#include <algorithm>
#include <cstddef>
#include <format>
#include <memory>
#include <stdexcept>
#include <string>
#include <string_view>
#include <tuple>
#include <vector>

namespace Opm::gpuistl {

namespace {

/// The PVT approaches that have a configuration for a CPU TypeTag, in the
/// order they are tried. Each entry needs a \c GpuDispatcherTraits
/// specialization and an explicit instantiation of \c makeImpl.
template <class CpuTypeTag>
struct GpuDispatcherConfigurations;

template <>
struct GpuDispatcherConfigurations<Properties::TTag::FlowProblemNohyst>
{
    using type = std::tuple<LiveOilDryGasPvtApproaches, LiveOilWetGasPvtApproaches>;
};

template <>
struct GpuDispatcherConfigurations<Properties::TTag::FlowGasWaterEnergyProblem>
{
    using type = std::tuple<Co2StorePvtApproaches>;
};

template <>
struct GpuDispatcherConfigurations<Properties::TTag::FlowGasWaterEnergyProblemGPU>
{
    using type = std::tuple<Co2StorePvtApproaches>;
};

std::string_view pvtName(GasPvtApproach approach)
{
    switch (approach) {
    case GasPvtApproach::NoGas: return "none (no gas phase)";
    case GasPvtApproach::DryGas: return "PVDG";
    case GasPvtApproach::DryHumidGas: return "PVTGW";
    case GasPvtApproach::WetHumidGas: return "PVTG with PVTGW";
    case GasPvtApproach::WetGas: return "PVTG";
    case GasPvtApproach::ThermalGas: return "thermal";
    case GasPvtApproach::Co2Gas: return "CO2STORE";
    case GasPvtApproach::H2Gas: return "H2STORE";
    }
    return "unknown";
}

std::string_view pvtName(OilPvtApproach approach)
{
    switch (approach) {
    case OilPvtApproach::NoOil: return "none (no oil phase)";
    case OilPvtApproach::LiveOil: return "PVTO";
    case OilPvtApproach::DeadOil: return "PVDO";
    case OilPvtApproach::ConstantCompressibilityOil: return "PVCDO";
    case OilPvtApproach::ThermalOil: return "thermal";
    case OilPvtApproach::BrineCo2: return "CO2STORE";
    case OilPvtApproach::BrineH2: return "H2STORE";
    case OilPvtApproach::ConstantRsDeadOil: return "RSCONST";
    }
    return "unknown";
}

std::string_view pvtName(WaterPvtApproach approach)
{
    switch (approach) {
    case WaterPvtApproach::NoWater: return "none (no water phase)";
    case WaterPvtApproach::ConstantCompressibilityBrine: return "PVTWSALT";
    case WaterPvtApproach::ConstantCompressibilityWater: return "PVTW";
    case WaterPvtApproach::ThermalWater: return "thermal";
    case WaterPvtApproach::BrineCo2: return "CO2STORE";
    case WaterPvtApproach::BrineH2: return "H2STORE";
    }
    return "unknown";
}

/// The error message for a deck whose PVT approaches match no configuration.
/// It names each phase whose PVT model is in no configuration, or, if each one
/// is in some configuration, the combination.
template <class FluidSystem, class... Configurations>
std::string unsupportedPvtMessage()
{
    const auto oil = FluidSystem::oilPvt().approach();
    const auto gas = FluidSystem::gasPvt().approach();
    const auto water = FluidSystem::waterPvt().approach();

    std::string message;
    const auto checkPhase = [&message](std::string_view phase, auto deckApproach,
                                       auto... supportedApproaches) {
        if (((deckApproach == supportedApproaches) || ...)) {
            return;
        }
        std::vector<std::string_view> names;
        std::string supported;
        for (const auto name : {pvtName(supportedApproaches)...}) {
            if (std::find(names.begin(), names.end(), name) == names.end()) {
                supported += (names.empty() ? "" : ", ") + std::string(name);
                names.push_back(name);
            }
        }
        message += std::format("{}{} PVT model {} (supported: {})",
                               message.empty() ? "" : "; ", phase, pvtName(deckApproach),
                               supported);
    };
    checkPhase("oil", oil, Configurations::oil...);
    checkPhase("gas", gas, Configurations::gas...);
    checkPhase("water", water, Configurations::water...);

    if (message.empty()) {
        message = std::format("combination of oil PVT model {}, gas PVT model {} and water "
                              "PVT model {}", pvtName(oil), pvtName(gas), pvtName(water));
    }
    return "GPU intensive quantities do not support the " + message + " for this simulator";
}

template <class FluidSystem, class PvtApproaches>
bool deckHasPvtApproaches()
{
    return FluidSystem::gasPvt().approach() == PvtApproaches::gas
        && FluidSystem::oilPvt().approach() == PvtApproaches::oil
        && FluidSystem::waterPvt().approach() == PvtApproaches::water;
}

/// The input checks that do not need the problem. The number of DoFs is checked
/// against the grid by the configuration.
template <class PrimaryVariables, class IntensiveQuantities>
void validateDispatcherInputs(const PrimaryVariables* const* primaryVariables,
                              IntensiveQuantities* const* intensiveQuantities,
                              std::size_t numDof,
                              unsigned timeIdx)
{
    if (timeIdx != 0) {
        OPM_THROW(std::invalid_argument,
                  "GPU intensive-quantities dispatcher supports only timeIdx 0");
    }

    if (primaryVariables == nullptr || intensiveQuantities == nullptr) {
        OPM_THROW(std::invalid_argument,
                  "GPU intensive-quantities dispatcher received a null pointer array");
    }

    for (std::size_t i = 0; i < numDof; ++i) {
        if (primaryVariables[i] == nullptr || intensiveQuantities[i] == nullptr) {
            OPM_THROW(std::invalid_argument,
                      "GPU intensive-quantities dispatcher received a null DoF entry");
        }
    }
}

} // namespace

template <class CpuTypeTag>
GpuBlackoilIntensiveQuantitiesDispatcher<CpuTypeTag>::GpuBlackoilIntensiveQuantitiesDispatcher()
    = default;

template <class CpuTypeTag>
GpuBlackoilIntensiveQuantitiesDispatcher<CpuTypeTag>::~GpuBlackoilIntensiveQuantitiesDispatcher()
    = default;

template <class CpuTypeTag>
const char* GpuBlackoilIntensiveQuantitiesDispatcher<CpuTypeTag>::configurationName() const
{
    return impl_ ? impl_->name() : nullptr;
}

template <class CpuTypeTag>
void GpuBlackoilIntensiveQuantitiesDispatcher<CpuTypeTag>::update(
    const Problem& cpuProblem,
    const PrimaryVariables* const* cpuPriVars,
    IntensiveQuantities* const* outIQ,
    std::size_t numDof,
    unsigned timeIdx)
{
    if (numDof == 0u) {
        return;
    }
    validateDispatcherInputs(cpuPriVars, outIQ, numDof, timeIdx);

    if (!impl_) {
        using FluidSystem = GetPropType<CpuTypeTag, Properties::FluidSystem>;
        // The first configuration with the deck's PVT approaches validates the
        // rest of the deck in makeImpl.
        [&]<class... Configurations>(std::tuple<Configurations...>*) {
            const auto tryConfiguration = [&]<class Approaches>(Approaches*) {
                if (!impl_ && deckHasPvtApproaches<FluidSystem, Approaches>()) {
                    impl_ = makeImpl<Approaches>(cpuProblem);
                }
            };
            (tryConfiguration(static_cast<Configurations*>(nullptr)), ...);
            if (!impl_) {
                const auto message = unsupportedPvtMessage<FluidSystem, Configurations...>();
                OPM_THROW(std::logic_error, message);
            }
        }(static_cast<typename GpuDispatcherConfigurations<CpuTypeTag>::type*>(nullptr));

        OpmLog::info(std::format(
            "[GpuBlackoilIntensiveQuantitiesDispatcher] initialized for {} cells ({})",
            numDof, impl_->name()));
    }

    impl_->update(cpuProblem, cpuPriVars, outIQ, numDof, timeIdx);
}

template class GpuBlackoilIntensiveQuantitiesDispatcher<Properties::TTag::FlowProblemNohyst>;
template class GpuBlackoilIntensiveQuantitiesDispatcher<Properties::TTag::FlowGasWaterEnergyProblem>;
template class GpuBlackoilIntensiveQuantitiesDispatcher<Properties::TTag::FlowGasWaterEnergyProblemGPU>;

} // namespace Opm::gpuistl
