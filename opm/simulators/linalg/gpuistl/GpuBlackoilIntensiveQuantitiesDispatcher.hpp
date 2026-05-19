// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
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
#ifndef OPM_GPU_BLACKOIL_INTENSIVE_QUANTITIES_DISPATCHER_HPP
#define OPM_GPU_BLACKOIL_INTENSIVE_QUANTITIES_DISPATCHER_HPP

#if HAVE_CUDA

#include <opm/models/utils/propertysystem.hh>

#include <cstddef>
#include <memory>

namespace Opm::Properties::TTag {
    struct FlowGasWaterEnergyProblem;
    struct FlowGasWaterEnergyProblemGPU;
}

namespace Opm::gpuistl {

/// Compile-time predicate: does this CPU \c TypeTag describe the specific
/// CO2STORE configuration (FlowGasWaterEnergyProblem) that the GPU
/// intensive-quantities dispatcher currently supports? Restricted to a
/// single TypeTag so that explicit instantiations of the dispatcher
/// remain manageable for the experimental prototype.
template <class CpuTypeTag>
struct GpuBlackoilIntensiveQuantitiesDispatcherSupport {
    static constexpr bool value = false;
};

template <>
struct GpuBlackoilIntensiveQuantitiesDispatcherSupport<
    Opm::Properties::TTag::FlowGasWaterEnergyProblem> {
    static constexpr bool value = true;
};

/// Enable the dispatcher for the GPU-assembly simulation TypeTag
/// (\c FlowGasWaterEnergyProblemGPU, declared in FlowGasWaterEnergyTypeTag.hpp).
/// This tag inherits all physics from \c FlowGasWaterEnergyProblem and adds
/// GPU-specific assembly properties in \c flow_gpu.cu.
template <>
struct GpuBlackoilIntensiveQuantitiesDispatcherSupport<
    Opm::Properties::TTag::FlowGasWaterEnergyProblemGPU> {
    static constexpr bool value = true;
};

/// Runs the BlackOil intensive-quantities update on the GPU for a given set
/// of degrees of freedom. The dispatcher is process-wide (one singleton per
/// CPU \c TypeTag). It lazily constructs the GPU-side \c GpuFlowProblem
/// from the supplied CPU problem on the first call.
///
/// On every call, primary variables for the requested DoFs are uploaded to
/// the device, the per-cell update kernel from
/// \c test_blackoilintensivequantities_gpu.cu is launched (one thread per
/// DoF), and the resulting intensive quantities are read back to host
/// memory, as the design plan in \c PLAN.md requires.
///
/// The class is a template on the CPU \c TypeTag and is explicitly
/// instantiated in the \c .cu translation unit.
template <class CpuTypeTag>
class GpuBlackoilIntensiveQuantitiesDispatcher
{
public:
    using Problem            = Opm::GetPropType<CpuTypeTag, Opm::Properties::Problem>;
    using PrimaryVariables   = Opm::GetPropType<CpuTypeTag, Opm::Properties::PrimaryVariables>;
    using IntensiveQuantities = Opm::GetPropType<CpuTypeTag, Opm::Properties::IntensiveQuantities>;

    GpuBlackoilIntensiveQuantitiesDispatcher();
    ~GpuBlackoilIntensiveQuantitiesDispatcher();

    GpuBlackoilIntensiveQuantitiesDispatcher(const GpuBlackoilIntensiveQuantitiesDispatcher&) = delete;
    GpuBlackoilIntensiveQuantitiesDispatcher&
    operator=(const GpuBlackoilIntensiveQuantitiesDispatcher&) = delete;

    /// Run the per-cell intensive-quantities update kernel on \p numDof
    /// DoFs. \p cpuPriVars[i] points at the CPU primary variables for DoF
    /// \c i. The GPU-computed BlackOil intensive quantities are written
    /// onto \p outIQ[i] field-by-field via
    /// \c BlackOilIntensiveQuantities::overlayBlackOilFieldsFrom (so the
    /// caller is expected to have run the CPU update first to fill in any
    /// fields the dispatcher does not overwrite, e.g. \c mobility_).
    ///
    /// Per-call host/device timings are accumulated and printed every call
    /// to \c std::cout, prefixed with
    /// \c "[GpuBlackoilIntensiveQuantitiesDispatcher]".
    void update(const Problem& cpuProblem,
                const PrimaryVariables* const* cpuPriVars,
                IntensiveQuantities* const* outIQ,
                std::size_t numDof);

private:
    struct Impl;
    std::unique_ptr<Impl> impl_;
};

} // namespace Opm::gpuistl

#endif // HAVE_CUDA

#endif // OPM_GPU_BLACKOIL_INTENSIVE_QUANTITIES_DISPATCHER_HPP
