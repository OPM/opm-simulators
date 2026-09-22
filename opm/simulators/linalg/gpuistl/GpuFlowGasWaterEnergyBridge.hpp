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

#include <opm/simulators/linalg/gpuistl/GpuVector.hpp>

#include <array>
#include <cstdint>
#include <memory>
#include <optional>
#include <stdexcept>
#include <string>
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
public:
    using Scalar = GetPropType<CpuTypeTag, Properties::Scalar>;

private:
    using Problem = GetPropType<CpuTypeTag, Properties::Problem>;
    using SolutionVector = GetPropType<CpuTypeTag, Properties::SolutionVector>;
    using HostPrimaryVariables = GetPropType<CpuTypeTag, Properties::PrimaryVariables>;
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
    using DevicePrimaryVariablesPublic = DevicePrimaryVariables;
    using DeviceIntensiveQuantitiesView = gpuistl::GpuView<DeviceIntensiveQuantities>;
    using DevicePrimaryVariablesView = gpuistl::GpuView<const DevicePrimaryVariables>;
    using DeviceModelViewPublic = DeviceModelView;

    GpuFlowGasWaterEnergyBridge()
    {
#if USE_HIP
        OPM_GPU_SAFE_CALL(hipStreamCreate(&stream_));
        OPM_GPU_SAFE_CALL(hipEventCreate(&propertyReady_));
        OPM_GPU_SAFE_CALL(hipEventCreate(&primaryReady_));
        OPM_GPU_SAFE_CALL(hipEventCreate(&solverReady_));
#else
        OPM_GPU_SAFE_CALL(cudaStreamCreate(&stream_));
        OPM_GPU_SAFE_CALL(cudaEventCreate(&propertyReady_));
        OPM_GPU_SAFE_CALL(cudaEventCreate(&primaryReady_));
        OPM_GPU_SAFE_CALL(cudaEventCreate(&solverReady_));
#endif
    }

    GpuFlowGasWaterEnergyBridge(const GpuFlowGasWaterEnergyBridge&) = delete;
    GpuFlowGasWaterEnergyBridge& operator=(const GpuFlowGasWaterEnergyBridge&) = delete;

    ~GpuFlowGasWaterEnergyBridge()
    {
        waitForLastWriter_();
#if USE_HIP
        OPM_GPU_WARN_IF_ERROR(hipEventDestroy(propertyReady_));
        OPM_GPU_WARN_IF_ERROR(hipEventDestroy(primaryReady_));
        OPM_GPU_WARN_IF_ERROR(hipEventDestroy(solverReady_));
        OPM_GPU_WARN_IF_ERROR(hipStreamDestroy(stream_));
#else
        OPM_GPU_WARN_IF_ERROR(cudaEventDestroy(propertyReady_));
        OPM_GPU_WARN_IF_ERROR(cudaEventDestroy(primaryReady_));
        OPM_GPU_WARN_IF_ERROR(cudaEventDestroy(solverReady_));
        OPM_GPU_WARN_IF_ERROR(cudaStreamDestroy(stream_));
#endif
    }

    struct TransferCounters {
        std::uint64_t primaryVariableUploads{0}, primaryVariableUploadBytes{0};
        std::uint64_t primaryVariableDownloads{0}, primaryVariableDownloadBytes{0};
        std::uint64_t correctionDownloads{0}, correctionDownloadBytes{0};
        std::uint64_t correctionHistoryUploads{0}, correctionHistoryUploadBytes{0};
        std::uint64_t intensiveQuantityDownloads{0}, intensiveQuantityDownloadBytes{0};
        std::uint64_t compactConvergenceDownloads{0}, compactConvergenceDownloadBytes{0};
        std::uint64_t relativeChangeDownloads{0}, relativeChangeDownloadBytes{0};
        std::uint64_t ownedBufferAllocations{0};
        std::uint64_t ownedBufferAllocationBytes{0};
        std::uint64_t staticUploadBatches{0};
        std::uint64_t staticUploadCalls{0}, staticUploadBytes{0};
        std::uint64_t solverCorrectionAllocations{0};
        std::uint64_t successfulNewtonUpdates{0};
        std::string lastCorrectionDownloadCause;
        std::string lastPrimaryVariableDownloadCause;
    };

    const TransferCounters& transferCounters() const { return counters_; }

    void recordCompatibilityCorrectionDownload(std::size_t bytes,
                                               const char* cause = "host solver interface")
    {
        ++counters_.correctionDownloads;
        counters_.correctionDownloadBytes += bytes;
        ++counters_.solverCorrectionAllocations;
        counters_.lastCorrectionDownloadCause = cause;
    }
    bool hasInitialized() const { return problemBuffer_ != nullptr; }
    bool initializedFor(std::size_t numDof) const { return hasInitialized() && numDof_ == numDof; }
    std::size_t numDof() const { return numDof_; }

    bool hasPrimaryVariables(unsigned timeIdx) const
    { return timeIdx < numTimeSlots && primaryGeneration_[timeIdx] != 0; }

    std::uint64_t primaryVariablesGeneration(unsigned timeIdx) const
    { validateTimeIdx_(timeIdx); return primaryGeneration_[timeIdx]; }

    bool hostPrimaryVariablesCurrent(unsigned timeIdx) const
    { return hasPrimaryVariables(timeIdx) && hostPrimaryGeneration_[timeIdx] == primaryGeneration_[timeIdx]; }

    // Only an explicit CPU mutation invalidates resident state. Taking a
    // non-const solution reference or refreshing its mirror does not.
    void invalidateHostEditedPrimaryVariables(unsigned timeIdx)
    {
        validateTimeIdx_(timeIdx);
        primaryGeneration_[timeIdx] = hostPrimaryGeneration_[timeIdx] = 0;
        deviceIqValid_[timeIdx] = false;
        iqGeneration_[timeIdx] = 0;
    }

    void updatePrimaryVariables(const Problem& cpuProblem,
                                const SolutionVector& solution,
                                unsigned timeIdx)
    { importHostPrimaryVariables(cpuProblem, solution, timeIdx); }

    void importHostPrimaryVariables(const Problem& cpuProblem,
                                    const SolutionVector& solution,
                                    unsigned timeIdx)
    {
        validateTimeIdx_(timeIdx);
        ensureInitialized_(cpuProblem, timeIdx, solution.size());
        importSlot_(solution, timeIdx);
    }

    void materializeHostPrimaryVariables(unsigned timeIdx, SolutionVector& destination)
    {
        validatePrimarySlot_(timeIdx);
        if (destination.size() != numDof_) {
            OPM_THROW(std::invalid_argument, "Invalid host primary-variable destination size");
        }
        if (hostPrimaryVariablesCurrent(timeIdx)) {
            return;
        }
        auto& mirror = hostPrimaryVariables_[timeIdx];
        copyDeviceToHost_(mirror.data(), primaryVariablesBuffer_[timeIdx]->data(),
                          numDof_ * sizeof(DevicePrimaryVariables));
        synchronizeStream_();
        for (std::size_t i = 0; i < numDof_; ++i) {
            destination[i] = HostPrimaryVariables(mirror[i]);
        }
        hostPrimaryGeneration_[timeIdx] = primaryGeneration_[timeIdx];
        ++counters_.primaryVariableDownloads;
        counters_.primaryVariableDownloadBytes += numDof_ * sizeof(DevicePrimaryVariables);
        counters_.lastPrimaryVariableDownloadCause = "CPU consumer";
    }

    gpuistl::GpuVector<Scalar>& correction() { return *correction_; }
    gpuistl::GpuVector<Scalar>& previousCorrection() { return *previousCorrection_; }
    gpuistl::GpuView<DevicePrimaryVariables> scratchPrimaryVariablesView()
    { return {scratchPrimaryVariables_->data(), numDof_}; }
    gpuistl::GpuView<std::uint8_t> switchedLastIterationView()
    { return {switchHistory_->data(), numDof_}; }
    gpuistl::GpuView<std::uint8_t> scratchSwitchHistoryView()
    { return {scratchSwitchHistory_->data(), numDof_}; }
    gpuistl::GpuView<std::uint32_t> updateStatusView()
    { return {updateStatus_->data(), 4}; }

    void clearUpdateStatus()
    {
        updateStatusChecked_ = false;
        zeroAsync_(updateStatus_->data(), 4 * sizeof(std::uint32_t));
    }

    std::array<std::uint32_t, 4> readUpdateStatus()
    {
        std::array<std::uint32_t, 4> status{};
        copyDeviceToHost_(status.data(), updateStatus_->data(), sizeof(status));
        synchronizeStream_();
        lastUpdateStatus_ = status;
        updateStatusChecked_ = true;
        return status;
    }

    // The solver uses its existing default stream. Both handoffs use events,
    // including the reset writes queued before the next solve.
    void orderSolverAfterProperties()
    {
#if USE_HIP
        OPM_GPU_SAFE_CALL(hipEventRecord(primaryReady_, stream_));
        OPM_GPU_SAFE_CALL(hipStreamWaitEvent(nullptr, primaryReady_, 0));
#else
        OPM_GPU_SAFE_CALL(cudaEventRecord(primaryReady_, stream_));
        OPM_GPU_SAFE_CALL(cudaStreamWaitEvent(nullptr, primaryReady_, 0));
#endif
    }

    void orderUpdateAfterSolver()
    {
#if USE_HIP
        OPM_GPU_SAFE_CALL(hipEventRecord(solverReady_, nullptr));
        OPM_GPU_SAFE_CALL(hipStreamWaitEvent(stream_, solverReady_, 0));
#else
        OPM_GPU_SAFE_CALL(cudaEventRecord(solverReady_, nullptr));
        OPM_GPU_SAFE_CALL(cudaStreamWaitEvent(stream_, solverReady_, 0));
#endif
    }

    // Call only after checking readUpdateStatus(). Scratch data is never
    // visible to properties or CPU consumers on a failed update.
    void commitDeviceUpdate(unsigned timeIdx = 0)
    {
        validatePrimarySlot_(timeIdx);
        if (!updateStatusChecked_ || lastUpdateStatus_[0] != 0) {
            OPM_THROW(std::logic_error, "Cannot publish an unchecked or failed GPU Newton update");
        }
        primaryVariablesBuffer_[timeIdx].swap(scratchPrimaryVariables_);
        switchHistory_.swap(scratchSwitchHistory_);
        primaryGeneration_[timeIdx] = ++generation_;
        deviceIqValid_[timeIdx] = false;
        iqGeneration_[timeIdx] = 0;
        recordPrimaryReady_();
        updateStatusChecked_ = false;
        ++counters_.successfulNewtonUpdates;
    }

    void resetCorrectionHistory()
    {
        if (previousCorrection_) {
            zeroAsync_(previousCorrection_->data(), previousCorrection_->dim() * sizeof(Scalar));
            recordPrimaryReady_();
        }
    }

    // CPU diagnostics and shadow checks must request and account for their
    // correction copy explicitly; ordinary Newton iterations never call this.
    // Explicit activation/fallback boundary. CPU stabilization may have
    // advanced dx_old_ while the resident path was inactive, including SOR.
    template<class HostVector>
    void importPreviousCorrection(const HostVector& source)
    {
        synchronizeStream_();
        previousCorrection_->copyFromHost(source);
        ++counters_.correctionHistoryUploads;
        counters_.correctionHistoryUploadBytes += previousCorrection_->dim() * sizeof(Scalar);
    }

    template<class HostVector>
    void materializeHostCorrection(HostVector& destination, const char* cause = "diagnostic")
    {
        synchronizeStream_();
        correction_->copyToHost(destination);
        counters_.lastCorrectionDownloadCause = cause;
        ++counters_.correctionDownloads;
        counters_.correctionDownloadBytes += correction_->dim() * sizeof(Scalar);
    }

    template<class HostVector>
    void materializeHostPreviousCorrection(HostVector& destination, const char* cause = "validation")
    {
        synchronizeStream_();
        previousCorrection_->copyToHost(destination);
        counters_.lastCorrectionDownloadCause = cause;
        ++counters_.correctionDownloads;
        counters_.correctionDownloadBytes += previousCorrection_->dim() * sizeof(Scalar);
    }

    template<class HostVector>
    void materializePreviousCorrection(HostVector& destination, const char* cause = "validation")
    { materializeHostPreviousCorrection(destination, cause); }

    void materializeCandidatePrimaryVariables(SolutionVector& destination)
    {
        if (destination.size() != numDof_) {
            OPM_THROW(std::invalid_argument, "Invalid shadow primary-variable destination size");
        }
        std::vector<DevicePrimaryVariables> values(numDof_);
        copyDeviceToHost_(values.data(), scratchPrimaryVariables_->data(),
                          numDof_ * sizeof(DevicePrimaryVariables));
        synchronizeStream_();
        for (std::size_t i = 0; i < numDof_; ++i) {
            destination[i] = HostPrimaryVariables(values[i]);
        }
        ++counters_.primaryVariableDownloads;
        counters_.primaryVariableDownloadBytes += numDof_ * sizeof(DevicePrimaryVariables);
        counters_.lastPrimaryVariableDownloadCause = "validation";
    }

    std::vector<std::uint8_t> materializeCandidateSwitchHistory()
    {
        std::vector<std::uint8_t> result(numDof_);
        copyDeviceToHost_(result.data(), scratchSwitchHistory_->data(), numDof_);
        synchronizeStream_();
        return result;
    }

    std::vector<std::uint8_t> materializeSwitchHistory()
    {
        std::vector<std::uint8_t> result(numDof_);
        copyDeviceToHost_(result.data(), switchHistory_->data(), numDof_);
        synchronizeStream_();
        return result;
    }

    void importSwitchHistory(const std::vector<std::uint8_t>& values)
    {
        if (values.size() != numDof_) {
            OPM_THROW(std::invalid_argument, "Invalid GPU Newton switch-history size");
        }
        copyHostToDevice_(switchHistory_->data(), values.data(), numDof_);
        // The caller owns this staging vector, so consume it before returning.
        synchronizeStream_();
    }

    auto stream() const
    {
        return stream_;
    }

    DevicePrimaryVariablesView primaryVariablesView(unsigned timeIdx) const
    {
        validatePrimarySlot_(timeIdx);
        return DevicePrimaryVariablesView(primaryVariablesBuffer_[timeIdx]->data(), numDof_);
    }

    DeviceIntensiveQuantitiesView intensiveQuantitiesView(unsigned timeIdx) const
    {
        validateReadySlot_(timeIdx);
        return DeviceIntensiveQuantitiesView(intensiveQuantitiesBuffer_[timeIdx]->data(), numDof_);
    }

    gpuistl::GpuView<Scalar> compactConvergenceView()
    {
        validateValidSlot_(0);
        return {compactConvergenceBuffer_->data(), compactConvergenceBuffer_->size()};
    }

    std::vector<Scalar> downloadCompactConvergence()
    {
        std::vector<Scalar> result(compactConvergenceBuffer_->size());
        copyDeviceToHost_(result.data(), compactConvergenceBuffer_->data(),
                          result.size() * sizeof(Scalar));
        synchronizeStream_();
        ++counters_.compactConvergenceDownloads;
        counters_.compactConvergenceDownloadBytes += result.size() * sizeof(Scalar);
        return result;
    }

    gpuistl::GpuView<Scalar> relativeChangeView()
    {
        validatePrimarySlot_(0);
        validatePrimarySlot_(1);
        return {relativeChangeBuffer_->data(), relativeChangeBuffer_->size()};
    }

    std::vector<Scalar> downloadRelativeChange()
    {
        std::vector<Scalar> result(relativeChangeBuffer_->size());
        copyDeviceToHost_(result.data(), relativeChangeBuffer_->data(),
                          result.size() * sizeof(Scalar));
        synchronizeStream_();
        ++counters_.relativeChangeDownloads;
        counters_.relativeChangeDownloadBytes += result.size() * sizeof(Scalar);
        return result;
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
        iqGeneration_[timeIdx] = primaryGeneration_[timeIdx];
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
        return hasIntensiveQuantities(0) && hasIntensiveQuantities(1);
    }

    bool hasIntensiveQuantities(unsigned timeIdx) const
    {
        return timeIdx < numTimeSlots && deviceIqValid_[timeIdx]
               && iqGeneration_[timeIdx] == primaryGeneration_[timeIdx];
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
        ++counters_.intensiveQuantityDownloads;
        counters_.intensiveQuantityDownloadBytes += numDof_ * sizeof(DeviceIntensiveQuantities);
        for (std::size_t i = 0; i < numDof_; ++i) {
            destination[i]->overlayBlackOilFieldsFrom(hostIntensiveQuantities_[i]);
        }
    }

    /*!
     * \brief Mirrors FvBaseDiscretization's time-level shift for device IQs.
     */
    void advanceTimeLevel()
    {
        if (!hasInitialized()) { return; }
        orderUpdateAfterSolver();
        if (hasPrimaryVariables(0)) {
            copyDeviceToDevice_(primaryVariablesBuffer_[1]->data(), primaryVariablesBuffer_[0]->data(),
                                numDof_ * sizeof(DevicePrimaryVariables));
            primaryGeneration_[1] = primaryGeneration_[0];
            hostPrimaryGeneration_[1] = 0;
            recordPrimaryReady_();
        }
        else {
            primaryGeneration_[1] = hostPrimaryGeneration_[1] = 0;
        }
        deviceIqValid_[1] = false;
        iqGeneration_[1] = 0;
        if (deviceIqValid_[0]) {
            copyDeviceToDevice_(intensiveQuantitiesBuffer_[1]->data(), intensiveQuantitiesBuffer_[0]->data(),
                                numDof_ * sizeof(DeviceIntensiveQuantities));
            recordPropertyReady(1);
        }
    }

    void restorePreviousSolution()
    {
        validatePrimarySlot_(1);
        orderUpdateAfterSolver();
        copyDeviceToDevice_(primaryVariablesBuffer_[0]->data(), primaryVariablesBuffer_[1]->data(),
                            numDof_ * sizeof(DevicePrimaryVariables));
        primaryGeneration_[0] = ++generation_;
        hostPrimaryGeneration_[0] = 0;
        deviceIqValid_[0] = false;
        iqGeneration_[0] = 0;
        recordPrimaryReady_();
        // Switch history deliberately survives rollback, like CPU wasSwitched_.
        resetCorrectionHistory();
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

        detail::GpuMemoryCounts allocations, staticData;
        const detail::ScopedGpuMemoryAccounting ownerAccounting(allocations);
        {
            const detail::ScopedGpuMemoryAccounting staticAccounting(staticData);
            auto fluidSystemBuffer =
                ::Opm::gpuistl::copy_to_gpu(CpuFluidSystem::getNonStaticInstance());
            fluidSystemBuffer_ = std::make_unique<FluidSystemBuffer>(std::move(fluidSystemBuffer));
            const auto fluidSystemView = ::Opm::gpuistl::make_view(*fluidSystemBuffer_);
            deviceFluidSystem_ = ::Opm::gpuistl::make_gpu_shared_ptr<DeviceFluidSystem>(fluidSystemView);
        }

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

        scratchPrimaryVariables_ = std::make_unique<gpuistl::GpuBuffer<DevicePrimaryVariables>>(numDof_);
        switchHistory_ = std::make_unique<gpuistl::GpuBuffer<std::uint8_t>>(numDof_);
        scratchSwitchHistory_ = std::make_unique<gpuistl::GpuBuffer<std::uint8_t>>(numDof_);
        updateStatus_ = std::make_unique<gpuistl::GpuBuffer<std::uint32_t>>(4);
        constexpr unsigned numEq = getPropValue<CpuTypeTag, Properties::NumEq>();
        correction_ = std::make_unique<gpuistl::GpuVector<Scalar>>(numDof_ * numEq);
        previousCorrection_ = std::make_unique<gpuistl::GpuVector<Scalar>>(numDof_ * numEq);
        compactConvergenceBuffer_ =
            std::make_unique<gpuistl::GpuBuffer<Scalar>>(numDof_ * numEq);
        relativeChangeBuffer_ = std::make_unique<gpuistl::GpuBuffer<Scalar>>(numDof_ * 2);
        zeroAsync_(switchHistory_->data(), numDof_);
        zeroAsync_(scratchSwitchHistory_->data(), numDof_);
        resetCorrectionHistory();

        std::vector<Scalar> volumes(numDof_);
        for (std::size_t i = 0; i < numDof_; ++i) {
            volumes[i] = cpuProblem.model().dofTotalVolume(static_cast<unsigned>(i));
        }
        {
            const detail::ScopedGpuMemoryAccounting staticAccounting(staticData);
            volumesBuffer_ = std::make_unique<gpuistl::GpuBuffer<Scalar>>(volumes);
            problemBuffer_ = std::make_unique<GpuProblemBuffer>(cpuProblem);
        }

        // The first device update overwrites currentTimeIdx.  The previous slot
        // can be seeded once from an already-valid CPU history cache, avoiding
        // any IQ transfer in steady-state GPU/GPU iterations.
        const unsigned previousTimeIdx = currentTimeIdx == 0 ? 1 : 0;
        importSlot_(cpuProblem.model().solution(previousTimeIdx), previousTimeIdx);
        initializePreviousSlot_(cpuProblem, currentTimeIdx);
        counters_.ownedBufferAllocations += allocations.allocations;
        counters_.ownedBufferAllocationBytes += allocations.allocationBytes;
        ++counters_.staticUploadBatches;
        counters_.staticUploadCalls += staticData.hostToDeviceCalls;
        counters_.staticUploadBytes += staticData.hostToDeviceBytes;
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
            iqGeneration_[previousTimeIdx] = hostCacheIsValid ? primaryGeneration_[previousTimeIdx] : 0;
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
        if (!hasIntensiveQuantities(timeIdx)) {
            OPM_THROW(std::logic_error,
                      "GPU property bridge attempted to consume an invalid intensive-quantity slot");
        }
    }

    void validatePrimarySlot_(unsigned timeIdx) const
    {
        validateReadySlot_(timeIdx);
        if (!hasPrimaryVariables(timeIdx)) {
            OPM_THROW(std::logic_error, "GPU bridge attempted to consume invalid primary variables");
        }
    }

    void importSlot_(const SolutionVector& solution, unsigned timeIdx)
    {
        if (solution.size() != numDof_) {
            OPM_THROW(std::invalid_argument, "GPU bridge received a mismatched solution size");
        }
        orderUpdateAfterSolver();
        // Host import is an explicit compatibility boundary. Wait before
        // rewriting staging memory that may still back an earlier upload.
        synchronizeStream_();
        auto& mirror = hostPrimaryVariables_[timeIdx];
        for (std::size_t i = 0; i < numDof_; ++i) {
            mirror[i] = DevicePrimaryVariables(solution[i]);
        }
        copyHostToDevice_(primaryVariablesBuffer_[timeIdx]->data(), mirror.data(),
                          numDof_ * sizeof(DevicePrimaryVariables));
        primaryGeneration_[timeIdx] = ++generation_;
        hostPrimaryGeneration_[timeIdx] = primaryGeneration_[timeIdx];
        deviceIqValid_[timeIdx] = false;
        iqGeneration_[timeIdx] = 0;
        recordPrimaryReady_();
        ++counters_.primaryVariableUploads;
        counters_.primaryVariableUploadBytes += numDof_ * sizeof(DevicePrimaryVariables);
    }

    void recordPrimaryReady_()
    {
#if USE_HIP
        OPM_GPU_SAFE_CALL(hipEventRecord(primaryReady_, stream_));
#else
        OPM_GPU_SAFE_CALL(cudaEventRecord(primaryReady_, stream_));
#endif
    }

    void synchronizeStream_() const
    {
#if USE_HIP
        OPM_GPU_SAFE_CALL(hipStreamSynchronize(stream_));
#else
        OPM_GPU_SAFE_CALL(cudaStreamSynchronize(stream_));
#endif
    }

    void zeroAsync_(void* destination, std::size_t bytes)
    {
#if USE_HIP
        OPM_GPU_SAFE_CALL(hipMemsetAsync(destination, 0, bytes, stream_));
#else
        OPM_GPU_SAFE_CALL(cudaMemsetAsync(destination, 0, bytes, stream_));
#endif
    }

    void copyHostToDevice_(void* destination, const void* source, std::size_t bytes)
    {
#if USE_HIP
        OPM_GPU_SAFE_CALL(hipMemcpyAsync(destination, source, bytes, hipMemcpyHostToDevice, stream_));
#else
        OPM_GPU_SAFE_CALL(cudaMemcpyAsync(destination, source, bytes, cudaMemcpyHostToDevice, stream_));
#endif
    }

    void copyDeviceToHost_(void* destination, const void* source, std::size_t bytes)
    {
#if USE_HIP
        OPM_GPU_SAFE_CALL(hipMemcpyAsync(destination, source, bytes, hipMemcpyDeviceToHost, stream_));
#else
        OPM_GPU_SAFE_CALL(cudaMemcpyAsync(destination, source, bytes, cudaMemcpyDeviceToHost, stream_));
#endif
    }

    void copyDeviceToDevice_(void* destination, const void* source, std::size_t bytes)
    {
#if USE_HIP
        OPM_GPU_SAFE_CALL(hipMemcpyAsync(destination, source, bytes, hipMemcpyDeviceToDevice, stream_));
#else
        OPM_GPU_SAFE_CALL(cudaMemcpyAsync(destination, source, bytes, cudaMemcpyDeviceToDevice, stream_));
#endif
    }

    void waitForLastWriter_() const
    {
        // Capture default-stream assembly readers as well as bridge writers
        // before replacing owners. This is a lifetime boundary, not a stage
        // transition, and does not synchronize unrelated device streams.
#if USE_HIP
        OPM_GPU_WARN_IF_ERROR(hipEventRecord(solverReady_, nullptr));
        OPM_GPU_WARN_IF_ERROR(hipStreamWaitEvent(stream_, solverReady_, 0));
        OPM_GPU_WARN_IF_ERROR(hipStreamSynchronize(stream_));
#else
        OPM_GPU_WARN_IF_ERROR(cudaEventRecord(solverReady_, nullptr));
        OPM_GPU_WARN_IF_ERROR(cudaStreamWaitEvent(stream_, solverReady_, 0));
        OPM_GPU_WARN_IF_ERROR(cudaStreamSynchronize(stream_));
#endif
    }

    void reset_()
    {
        for (unsigned timeIdx = 0; timeIdx < numTimeSlots; ++timeIdx) {
            primaryVariablesBuffer_[timeIdx].reset();
            intensiveQuantitiesBuffer_[timeIdx].reset();
            hostPrimaryVariables_[timeIdx].clear();
            deviceIqValid_[timeIdx] = false;
            primaryGeneration_[timeIdx] = hostPrimaryGeneration_[timeIdx] = iqGeneration_[timeIdx] = 0;
        }
        scratchPrimaryVariables_.reset();
        correction_.reset();
        previousCorrection_.reset();
        switchHistory_.reset();
        scratchSwitchHistory_.reset();
        updateStatus_.reset();
        compactConvergenceBuffer_.reset();
        relativeChangeBuffer_.reset();
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
    hipEvent_t primaryReady_{nullptr};
    hipEvent_t solverReady_{nullptr};
#else
    cudaStream_t stream_{nullptr};
    cudaEvent_t propertyReady_{nullptr};
    cudaEvent_t primaryReady_{nullptr};
    cudaEvent_t solverReady_{nullptr};
#endif
    std::size_t numDof_{0};
    TransferCounters counters_{};
    std::array<std::uint32_t, 4> lastUpdateStatus_{};
    bool updateStatusChecked_{false};
    std::uint64_t generation_{0};
    std::array<std::uint64_t, numTimeSlots> primaryGeneration_{}, hostPrimaryGeneration_{}, iqGeneration_{};
    std::unique_ptr<gpuistl::GpuBuffer<DevicePrimaryVariables>> scratchPrimaryVariables_;
    std::unique_ptr<gpuistl::GpuVector<Scalar>> correction_, previousCorrection_;
    std::unique_ptr<gpuistl::GpuBuffer<std::uint8_t>> switchHistory_, scratchSwitchHistory_;
    std::unique_ptr<gpuistl::GpuBuffer<std::uint32_t>> updateStatus_;
    std::unique_ptr<gpuistl::GpuBuffer<Scalar>> compactConvergenceBuffer_;
    std::unique_ptr<gpuistl::GpuBuffer<Scalar>> relativeChangeBuffer_;
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
