// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*
  Copyright 2026 Equinor ASA

  This file is part of the Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.
*/
/*!
 * \file
 *
 * \brief Typed owner for the GPU property-evaluation to TPFA-assembly handoff.
 */
#ifndef OPM_GPU_FLOW_GASWATER_ENERGY_BRIDGE_HPP
#define OPM_GPU_FLOW_GASWATER_ENERGY_BRIDGE_HPP

#if HAVE_CUDA

#include <opm/models/blackoil/blackoilintensivequantities.hh>
#include <opm/models/utils/propertysystem.hh>

#include <opm/material/fluidsystems/BlackOilFluidSystem.hpp>
#include <opm/simulators/flow/SimpleFIBlackOilModel.hpp>
#include <opm/simulators/linalg/gpuistl/GpuFlowGasWaterEnergyContract.hpp>

#include <array>
#include <memory>
#include <optional>
#include <stdexcept>
#include <type_traits>
#include <vector>

#if USE_HIP
#include <hip/hip_runtime.h>
#else
#include <cuda_runtime.h>
#endif

namespace Opm::gpuistl {

/*!
 * \brief Owns the persistent, typed device state shared by GPU IQ evaluation
 *        and TPFA assembly.
 *
 * The bridge retains all allocations backing its views.  It is deliberately
 * model-owned (through the dispatcher) so assembly cannot outlive the IQ
 * cache, fluid-system object, or problem data it consumes.
 */
template <class CpuTypeTag, class DeviceTypeTag>
class GpuFlowGasWaterEnergyBridge
{
    static constexpr unsigned numTimeSlots = 2;

    using CpuPrototypeTypeTag = Properties::TTag::FlowGasWaterEnergyCpuDeviceContract<CpuTypeTag>;
    using Scalar = GetPropType<CpuTypeTag, Properties::Scalar>;
    using Problem = GetPropType<CpuTypeTag, Properties::Problem>;
    using SolutionVector = GetPropType<CpuTypeTag, Properties::SolutionVector>;
    using HostIntensiveQuantities = GetPropType<CpuTypeTag, Properties::IntensiveQuantities>;
    using CpuFluidSystem = GetPropType<CpuTypeTag, Properties::FluidSystem>;

    using DevicePrimaryVariables = GetPropType<DeviceTypeTag, Properties::PrimaryVariables>;
    using DeviceIntensiveQuantities = BlackOilIntensiveQuantities<DeviceTypeTag>;
    using DeviceFluidSystem = GetPropType<DeviceTypeTag, Properties::FluidSystem>;
    using DeviceProblem = GetPropType<DeviceTypeTag, Properties::Problem>;
    using DeviceModelView = SimpleFIBlackOilModel<DeviceTypeTag, gpuistl::GpuView>;

    using DynamicCpuFluidSystem =
        std::remove_reference_t<decltype(CpuFluidSystem::getNonStaticInstance())>;
    using FluidSystemBuffer =
        decltype(::Opm::gpuistl::copy_to_gpu(std::declval<DynamicCpuFluidSystem&>()));
    using MaterialLawManager =
        typename GetProp<DeviceTypeTag, Properties::MaterialLaw>::EclMaterialLawManager;
    using GpuMaterialLawManager =
        EclMaterialLaw::GpuManager<typename MaterialLawManager::Traits,
                                   typename MaterialLawManager::GasOilLaw,
                                   typename MaterialLawManager::OilWaterLaw,
                                   gpuistl::GpuBuffer,
                                   typename MaterialLawManager::MaterialLaw>;
    using GpuThermalLawManager =
        EclThermalLaw::GpuManager<Scalar, DeviceFluidSystem, gpuistl::GpuBuffer, gpuistl::GpuView>;
    using GpuProblemBuffer =
        GpuFlowProblem<Scalar, GpuMaterialLawManager, gpuistl::GpuBuffer, GpuThermalLawManager>;
    using GpuProblemView = decltype(::Opm::gpuistl::make_view(std::declval<GpuProblemBuffer&>()));

    static_assert(std::is_same_v<DeviceIntensiveQuantities,
                                 GetPropType<DeviceTypeTag, Properties::IntensiveQuantities>>);
    static_assert(std::is_same_v<DeviceProblem, GpuProblemView>);

public:
    using DeviceTypeTagPublic = DeviceTypeTag;
    using DeviceIntensiveQuantitiesView = gpuistl::GpuView<DeviceIntensiveQuantities>;
    using DevicePrimaryVariablesView = gpuistl::GpuView<const DevicePrimaryVariables>;
    using DeviceModelViewPublic = DeviceModelView;

    GpuFlowGasWaterEnergyBridge()
    {
#if USE_HIP
        OPM_GPU_SAFE_CALL(hipStreamCreate(&stream_));
        OPM_GPU_SAFE_CALL(hipEventCreate(&propertyReady_));
#else
        OPM_GPU_SAFE_CALL(cudaStreamCreate(&stream_));
        OPM_GPU_SAFE_CALL(cudaEventCreate(&propertyReady_));
#endif
    }

    GpuFlowGasWaterEnergyBridge(const GpuFlowGasWaterEnergyBridge&) = delete;
    GpuFlowGasWaterEnergyBridge& operator=(const GpuFlowGasWaterEnergyBridge&) = delete;

    ~GpuFlowGasWaterEnergyBridge()
    {
        waitForLastWriter_();
#if USE_HIP
        OPM_GPU_WARN_IF_ERROR(hipEventDestroy(propertyReady_));
        OPM_GPU_WARN_IF_ERROR(hipStreamDestroy(stream_));
#else
        OPM_GPU_WARN_IF_ERROR(cudaEventDestroy(propertyReady_));
        OPM_GPU_WARN_IF_ERROR(cudaStreamDestroy(stream_));
#endif
    }

    void updatePrimaryVariables(const Problem& cpuProblem,
                                const SolutionVector& solution,
                                unsigned timeIdx)
    {
        validateTimeIdx_(timeIdx);
        ensureInitialized_(cpuProblem, timeIdx, solution.size());

        if (solution.size() != numDof_) {
            OPM_THROW(std::invalid_argument,
                      "GPU property bridge received a solution with a mismatched number of DoFs");
        }

        auto& hostPrimaryVariables = hostPrimaryVariables_[timeIdx];
        for (std::size_t i = 0; i < numDof_; ++i) {
            hostPrimaryVariables[i] = DevicePrimaryVariables(solution[i]);
        }

#if USE_HIP
        OPM_GPU_SAFE_CALL(hipMemcpyAsync(primaryVariablesBuffer_[timeIdx]->data(),
                                         hostPrimaryVariables.data(),
                                         numDof_ * sizeof(DevicePrimaryVariables),
                                         hipMemcpyHostToDevice,
                                         stream_));
#else
        OPM_GPU_SAFE_CALL(cudaMemcpyAsync(primaryVariablesBuffer_[timeIdx]->data(),
                                          hostPrimaryVariables.data(),
                                          numDof_ * sizeof(DevicePrimaryVariables),
                                          cudaMemcpyHostToDevice,
                                          stream_));
#endif
        deviceIqValid_[timeIdx] = false;
        lastWriterTimeIdx_ = timeIdx;
    }

    auto stream() const
    {
        return stream_;
    }

    DevicePrimaryVariablesView primaryVariablesView(unsigned timeIdx) const
    {
        validateReadySlot_(timeIdx);
        return DevicePrimaryVariablesView(primaryVariablesBuffer_[timeIdx]->data(), numDof_);
    }

    DeviceIntensiveQuantitiesView intensiveQuantitiesView(unsigned timeIdx) const
    {
        validateReadySlot_(timeIdx);
        return DeviceIntensiveQuantitiesView(intensiveQuantitiesBuffer_[timeIdx]->data(), numDof_);
    }

    void recordPropertyReady(unsigned timeIdx)
    {
        validateReadySlot_(timeIdx);
#if USE_HIP
        OPM_GPU_SAFE_CALL(hipEventRecord(propertyReady_, stream_));
#else
        OPM_GPU_SAFE_CALL(cudaEventRecord(propertyReady_, stream_));
#endif
        deviceIqValid_[timeIdx] = true;
        lastWriterTimeIdx_ = timeIdx;
    }

    /*!
     * \brief Order default-stream TPFA kernels after the requested IQ writer.
     *
     * No host wait or device-wide synchronize is performed here.
     */
    void waitForAssembly(unsigned timeIdx) const
    {
        validateValidSlot_(timeIdx);
        if (lastWriterTimeIdx_ != invalidTimeIdx) {
#if USE_HIP
            OPM_GPU_SAFE_CALL(hipStreamWaitEvent(nullptr, propertyReady_, 0));
#else
            OPM_GPU_SAFE_CALL(cudaStreamWaitEvent(nullptr, propertyReady_, 0));
#endif
        }
    }

    DeviceModelView modelView() const
    {
        validateValidSlot_(0);
        validateValidSlot_(1);
        return DeviceModelView(::Opm::gpuistl::make_view(*intensiveQuantitiesBuffer_[0]),
                               ::Opm::gpuistl::make_view(*intensiveQuantitiesBuffer_[1]),
                               ::Opm::gpuistl::make_view(*volumesBuffer_));
    }

    GpuProblemView flowProblemView() const
    {
        if (!problemBuffer_) {
            OPM_THROW(std::logic_error, "GPU property bridge has not been initialized");
        }
        return ::Opm::gpuistl::make_view(*problemBuffer_);
    }

    /*!
     * \brief Returns the device-resident fluid-system object used by all IQs.
     *
     * This reference is only passed to typed conversion helpers.  It is never
     * dereferenced on the host; the resulting IQ stores its device address.
     */
    const DeviceFluidSystem& deviceFluidSystem() const
    {
        if (!deviceFluidSystem_) {
            OPM_THROW(std::logic_error, "GPU property bridge has not been initialized");
        }
        return *deviceFluidSystem_;
    }

    bool hasModelView() const
    {
        return deviceIqValid_[0] && deviceIqValid_[1];
    }

    bool hasIntensiveQuantities(unsigned timeIdx) const
    {
        return timeIdx < numTimeSlots && deviceIqValid_[timeIdx];
    }

    /*!
     * \brief Explicit CPU-boundary materialization of one IQ time slot.
     */
    void materializeHostIntensiveQuantities(unsigned timeIdx,
                                            HostIntensiveQuantities* const* destination,
                                            std::size_t numDof)
    {
        validateValidSlot_(timeIdx);
        if (destination == nullptr || numDof != numDof_) {
            OPM_THROW(std::invalid_argument,
                      "GPU property bridge received an invalid host IQ materialization destination");
        }
        for (std::size_t i = 0; i < numDof_; ++i) {
            if (destination[i] == nullptr) {
                OPM_THROW(std::invalid_argument,
                          "GPU property bridge received a null host IQ materialization entry");
            }
        }

        if (lastWriterTimeIdx_ != invalidTimeIdx) {
#if USE_HIP
            OPM_GPU_SAFE_CALL(hipEventSynchronize(propertyReady_));
#else
            OPM_GPU_SAFE_CALL(cudaEventSynchronize(propertyReady_));
#endif
        }
        hostIntensiveQuantities_.assign(numDof_, *prototype_);
        intensiveQuantitiesBuffer_[timeIdx]->copyToHost(hostIntensiveQuantities_);
        for (std::size_t i = 0; i < numDof_; ++i) {
            destination[i]->overlayBlackOilFieldsFrom(hostIntensiveQuantities_[i]);
        }
    }

    /*!
     * \brief Mirrors FvBaseDiscretization's time-level shift for device IQs.
     */
    void advanceTimeLevel()
    {
        if (!deviceIqValid_[0]) {
            deviceIqValid_[1] = false;
            return;
        }
#if USE_HIP
        OPM_GPU_SAFE_CALL(hipMemcpyAsync(intensiveQuantitiesBuffer_[1]->data(),
                                         intensiveQuantitiesBuffer_[0]->data(),
                                         numDof_ * sizeof(DeviceIntensiveQuantities),
                                         hipMemcpyDeviceToDevice,
                                         stream_));
#else
        OPM_GPU_SAFE_CALL(cudaMemcpyAsync(intensiveQuantitiesBuffer_[1]->data(),
                                          intensiveQuantitiesBuffer_[0]->data(),
                                          numDof_ * sizeof(DeviceIntensiveQuantities),
                                          cudaMemcpyDeviceToDevice,
                                          stream_));
#endif
        recordPropertyReady(/*timeIdx=*/1);
    }

private:
    static constexpr unsigned invalidTimeIdx = numTimeSlots;

    void ensureInitialized_(const Problem& cpuProblem, unsigned currentTimeIdx, std::size_t numDof)
    {
        if (numDof == 0u) {
            OPM_THROW(std::invalid_argument, "GPU property bridge cannot be initialized with zero DoFs");
        }
        if (numDof_ == numDof && problemBuffer_) {
            return;
        }

        waitForLastWriter_();
        reset_();
        numDof_ = numDof;

        auto fluidSystemBuffer =
            ::Opm::gpuistl::copy_to_gpu(CpuFluidSystem::getNonStaticInstance());
        fluidSystemBuffer_ = std::make_unique<FluidSystemBuffer>(std::move(fluidSystemBuffer));
        const auto fluidSystemView = ::Opm::gpuistl::make_view(*fluidSystemBuffer_);
        deviceFluidSystem_ = ::Opm::gpuistl::make_gpu_shared_ptr<DeviceFluidSystem>(fluidSystemView);

        BlackOilIntensiveQuantities<CpuPrototypeTypeTag> cpuPrototype;
        prototype_ =
            cpuPrototype.template withOtherFluidSystem<DeviceTypeTag>(deviceFluidSystem());

        std::vector<DeviceIntensiveQuantities> initialIq(numDof_, *prototype_);
        for (unsigned timeIdx = 0; timeIdx < numTimeSlots; ++timeIdx) {
            hostPrimaryVariables_[timeIdx].resize(numDof_);
            primaryVariablesBuffer_[timeIdx] =
                std::make_unique<gpuistl::GpuBuffer<DevicePrimaryVariables>>(numDof_);
            intensiveQuantitiesBuffer_[timeIdx] =
                std::make_unique<gpuistl::GpuBuffer<DeviceIntensiveQuantities>>(initialIq);
        }

        std::vector<Scalar> volumes(numDof_);
        for (std::size_t i = 0; i < numDof_; ++i) {
            volumes[i] = cpuProblem.model().dofTotalVolume(static_cast<unsigned>(i));
        }
        volumesBuffer_ = std::make_unique<gpuistl::GpuBuffer<Scalar>>(volumes);
        problemBuffer_ = std::make_unique<GpuProblemBuffer>(cpuProblem);

        // The first device update overwrites currentTimeIdx.  The previous slot
        // can be seeded once from an already-valid CPU history cache, avoiding
        // any IQ transfer in steady-state GPU/GPU iterations.
        initializePreviousSlot_(cpuProblem, currentTimeIdx);
    }

    void initializePreviousSlot_(const Problem& cpuProblem, unsigned currentTimeIdx)
    {
        const unsigned previousTimeIdx = currentTimeIdx == 0 ? 1 : 0;
        const auto& cpuModel = cpuProblem.model();
        if (cpuModel.intensiveQuantityCache().size() <= previousTimeIdx) {
            return;
        }

        if constexpr (!getPropValue<CpuTypeTag, Properties::EnableDiffusion>()
                      && !getPropValue<CpuTypeTag, Properties::EnableDispersion>()) {
            std::vector<DeviceIntensiveQuantities> initialIq;
            initialIq.reserve(numDof_);
            bool hostCacheIsValid = true;
            for (std::size_t i = 0; i < numDof_; ++i) {
                hostCacheIsValid =
                    hostCacheIsValid
                    && cpuModel.cachedIntensiveQuantities(static_cast<unsigned>(i), previousTimeIdx) != nullptr;
                initialIq.emplace_back(
                    cpuModel.intensiveQuantityCache()[previousTimeIdx][i]
                        .template withOtherFluidSystem<DeviceTypeTag>(deviceFluidSystem()));
            }
            intensiveQuantitiesBuffer_[previousTimeIdx]->copyFromHost(initialIq);
            deviceIqValid_[previousTimeIdx] = hostCacheIsValid;
        }
    }

    void validateTimeIdx_(unsigned timeIdx) const
    {
        if (timeIdx >= numTimeSlots) {
            OPM_THROW(std::invalid_argument,
                      "GPU property bridge supports only current and previous IQ time slots");
        }
    }

    void validateReadySlot_(unsigned timeIdx) const
    {
        validateTimeIdx_(timeIdx);
        if (!intensiveQuantitiesBuffer_[timeIdx] || !primaryVariablesBuffer_[timeIdx]) {
            OPM_THROW(std::logic_error, "GPU property bridge has not been initialized");
        }
    }

    void validateValidSlot_(unsigned timeIdx) const
    {
        validateReadySlot_(timeIdx);
        if (!deviceIqValid_[timeIdx]) {
            OPM_THROW(std::logic_error,
                      "GPU property bridge attempted to consume an invalid intensive-quantity slot");
        }
    }

    void waitForLastWriter_() const
    {
        if (lastWriterTimeIdx_ != invalidTimeIdx) {
#if USE_HIP
            OPM_GPU_WARN_IF_ERROR(hipEventSynchronize(propertyReady_));
#else
            OPM_GPU_WARN_IF_ERROR(cudaEventSynchronize(propertyReady_));
#endif
        }
    }

    void reset_()
    {
        for (unsigned timeIdx = 0; timeIdx < numTimeSlots; ++timeIdx) {
            primaryVariablesBuffer_[timeIdx].reset();
            intensiveQuantitiesBuffer_[timeIdx].reset();
            hostPrimaryVariables_[timeIdx].clear();
            deviceIqValid_[timeIdx] = false;
        }
        volumesBuffer_.reset();
        problemBuffer_.reset();
        deviceFluidSystem_.reset();
        fluidSystemBuffer_.reset();
        hostIntensiveQuantities_.clear();
        prototype_.reset();
        lastWriterTimeIdx_ = invalidTimeIdx;
        numDof_ = 0;
    }

#if USE_HIP
    hipStream_t stream_{nullptr};
    hipEvent_t propertyReady_{nullptr};
#else
    cudaStream_t stream_{nullptr};
    cudaEvent_t propertyReady_{nullptr};
#endif
    std::size_t numDof_{0};
    std::array<std::unique_ptr<gpuistl::GpuBuffer<DevicePrimaryVariables>>, numTimeSlots>
        primaryVariablesBuffer_{};
    std::array<std::unique_ptr<gpuistl::GpuBuffer<DeviceIntensiveQuantities>>, numTimeSlots>
        intensiveQuantitiesBuffer_{};
    std::array<std::vector<DevicePrimaryVariables>, numTimeSlots> hostPrimaryVariables_{};
    std::vector<DeviceIntensiveQuantities> hostIntensiveQuantities_{};
    std::optional<DeviceIntensiveQuantities> prototype_;
    std::unique_ptr<gpuistl::GpuBuffer<Scalar>> volumesBuffer_;
    std::unique_ptr<FluidSystemBuffer> fluidSystemBuffer_;
    std::shared_ptr<DeviceFluidSystem> deviceFluidSystem_;
    std::unique_ptr<GpuProblemBuffer> problemBuffer_;
    std::array<bool, numTimeSlots> deviceIqValid_{false, false};
    mutable unsigned lastWriterTimeIdx_{invalidTimeIdx};
};

} // namespace Opm::gpuistl

#endif // HAVE_CUDA

#endif // OPM_GPU_FLOW_GASWATER_ENERGY_BRIDGE_HPP
