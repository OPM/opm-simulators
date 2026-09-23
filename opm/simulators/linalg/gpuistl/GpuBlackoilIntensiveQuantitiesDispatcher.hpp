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
#include <opm/models/blackoil/blackoilnewtonmethodparams.hpp>

#include <cstddef>
#include <memory>
#include <vector>

namespace Opm::Parameters {

// Experimental: route the per-element BlackOilIntensiveQuantities update through the
// GPU dispatcher instead of computing it on the CPU. This currently supports
// only the gas-water energy configuration.
struct ExperimentalComputePropertiesOnGpu { static constexpr bool value = false; };

} // namespace Opm::Parameters

namespace Opm::Properties::TTag {
    struct FlowGasWaterEnergyProblem;
    struct FlowGasWaterEnergyProblemGPU;
    template <template <class> class Storage>
    struct FlowGasWaterEnergyDeviceTypeTag;
}

namespace Opm::gpuistl {

template <class T>
class GpuView;
template <class T>
class GpuVector;

template <class CpuTypeTag,
          class DeviceTypeTag =
              Properties::TTag::FlowGasWaterEnergyDeviceTypeTag<GpuView>>
class GpuFlowGasWaterEnergyBridge;

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

/// Runs the supported BlackOil intensive-quantities update on the GPU for all
/// grid degrees of freedom. Each dispatcher instance owns a persistent,
/// typed property/assembly bridge which is lazily constructed from the CPU
/// problem on the first call.
///
/// On every call, primary variables for the requested DoFs are uploaded to
/// reusable device storage and the per-cell update kernel is launched (one
/// thread per DoF). The resulting intensive quantities remain device-resident
/// until an explicit CPU materialization request. Device allocations are
/// rebuilt only if the number of DoFs changes.
/// The supported gas-water thermal configuration computes the complete
/// intensive-quantity state needed by this dispatcher, including mobility.
///
/// The class is a template on the CPU \c TypeTag and is explicitly
/// instantiated in the \c .cu translation unit.
template <class CpuTypeTag>
class GpuBlackoilIntensiveQuantitiesDispatcher
{
public:
    using Scalar = Opm::GetPropType<CpuTypeTag, Opm::Properties::Scalar>;
    using Problem            = Opm::GetPropType<CpuTypeTag, Opm::Properties::Problem>;
    using PrimaryVariables = Opm::GetPropType<CpuTypeTag, Opm::Properties::PrimaryVariables>;
    using SolutionVector = Opm::GetPropType<CpuTypeTag, Opm::Properties::SolutionVector>;
    using IntensiveQuantities = Opm::GetPropType<CpuTypeTag, Opm::Properties::IntensiveQuantities>;
    using Bridge = GpuFlowGasWaterEnergyBridge<CpuTypeTag>;

    GpuBlackoilIntensiveQuantitiesDispatcher();
    ~GpuBlackoilIntensiveQuantitiesDispatcher();

    GpuBlackoilIntensiveQuantitiesDispatcher(const GpuBlackoilIntensiveQuantitiesDispatcher&) = delete;
    GpuBlackoilIntensiveQuantitiesDispatcher&
    operator=(const GpuBlackoilIntensiveQuantitiesDispatcher&) = delete;

    /// Run the per-cell intensive-quantities update kernel on the complete
    /// CPU solution. The primary-variable transfer and kernel are ordered on
    /// the bridge stream, and the typed device-IQ result is recorded there.
    void update(const Problem& cpuProblem,
                const SolutionVector& solution,
                unsigned timeIdx);

    void evaluateResident(unsigned timeIdx);
    unsigned applyNewtonUpdate(const Problem& problem,
                               const BlackoilNewtonParams<Scalar>& params,
                               Scalar relaxation, bool useSOR, bool stabilize,
                               bool validate);
    bool hasBridge() const;
    std::vector<Scalar> compactConvergenceFactors();
    std::vector<Scalar> compactRelativeChange();
    std::vector<Scalar> compactRockCompactionState();
    void reportTransferCounters() const;

    /// Explicit CPU-boundary materialization for legacy CPU consumers.
    void materializeHostIntensiveQuantities(unsigned timeIdx,
                                            IntensiveQuantities* const* destination,
                                            std::size_t numDof);

    /// Materialize one source cell without downloading the complete IQ slot.
    void materializeHostIntensiveQuantity(unsigned timeIdx, unsigned globalIdx,
                                         IntensiveQuantities& destination);

    /// Form true-IMPES CPR weights directly from resident storage derivatives.
    /// Return false for a singular block so the CPU path retains its diagnostics.
    bool computeTrueImpesWeights(GpuVector<Scalar>& weights, Scalar timeStepSize);

    bool hasDeviceModelView() const;
    const Bridge& bridge() const;
    Bridge& bridge();

private:
    struct Impl;
    std::unique_ptr<Impl> impl_;
};

} // namespace Opm::gpuistl

#endif // HAVE_CUDA

#endif // OPM_GPU_BLACKOIL_INTENSIVE_QUANTITIES_DISPATCHER_HPP
