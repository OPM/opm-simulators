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
/*!
 * \file
 *
 * GPU-portable, simplified version of \c Opm::EclThermalLawManager.
 *
 * Mirrors the design of \c Opm::EclMaterialLaw::GpuManager:
 *
 *   - Templated on the outer storage container (per-cell / per-region
 *     bulk arrays of params): \c VectorWithDefaultAllocator on the CPU,
 *     \c GpuBuffer for owning device memory, \c GpuView for non-owning
 *     device storage usable from a kernel.
 *   - Templated on the inner storage used by the per-region SPECROCK
 *     sample tables (one tiny array per region).
 *   - The CPU \c Opm::EclThermalLawManager is treated as a *builder*: a
 *     dedicated constructor extracts the per-region SPECROCK sample
 *     tables and the per-cell THCONR coefficients and uploads them to
 *     device memory.
 *
 * Currently only the SPECROCK solid-energy approach and the THCONR
 * thermal-conduction approach are supported by the builder, since these
 * are the only ones used by the CO2STORE+THERMAL setup that the GPU
 * dispatcher targets. The builder throws via \c OPM_THROW for any other
 * approach.
 */
#ifndef OPM_GPU_ECL_THERMAL_LAW_MANAGER_HPP
#define OPM_GPU_ECL_THERMAL_LAW_MANAGER_HPP

#include <opm/common/ErrorMacros.hpp>
#include <opm/common/utility/VectorWithDefaultAllocator.hpp>
#include <opm/common/utility/gpuDecorators.hpp>

#include <opm/material/thermal/EclSolidEnergyLawMultiplexerParams.hpp>
#include <opm/material/thermal/EclThconrLawParams.hpp>
#include <opm/material/thermal/EclThermalConductionLawMultiplexerParams.hpp>
#include <opm/material/thermal/EclSpecrockLawParams.hpp>
#include <opm/material/thermal/EclThermalLawManager.hpp>

#include <opm/simulators/linalg/gpuistl/GpuBuffer.hpp>
#include <opm/simulators/linalg/gpuistl/GpuView.hpp>

#include <algorithm>
#include <cstddef>
#include <stdexcept>
#include <type_traits>
#include <utility>
#include <vector>

namespace Opm::EclThermalLaw {

namespace detail {

    /*!
     * \brief Holder for the per-region SPECROCK sample buffers when the
     *        manager owns its device memory (i.e. the outer storage is
     *        \c GpuBuffer). Each per-region
     *        \c EclSpecrockLawParams<Scalar, GpuView> stored in the
     *        bulk array references one of the buffers in this holder, so
     *        the holder must outlive every kernel invocation.
     *
     *        For non-owning outer storage (host or GpuView), the holder
     *        is empty.
     */
    template <class Scalar, bool Owning>
    struct SpecrockSampleHolder {
    };

    template <class Scalar>
    struct SpecrockSampleHolder<Scalar, true> {
        std::vector<::Opm::gpuistl::GpuBuffer<Scalar>> sampleBuffers {};
    };

} // namespace detail

/*!
 * \brief Minimal, GPU-portable thermal-law manager.
 *
 * \tparam ScalarT       Floating point type used for the per-cell and
 *                       per-region sample tables.
 * \tparam FluidSystemT  Fluid system associated with this manager.
 * \tparam OuterStorage  Container used for the per-cell / per-region
 *                       bulk arrays of params and the SATNUM index map.
 * \tparam SampleStorage Container used for the per-region SPECROCK
 *                       internal-energy / temperature sample tables.
 *
 * Typical instantiations:
 *   - CPU-only:                       \c GpuManager<S, FS>
 *   - GPU owning (built from CPU):    \c GpuManager<S, FS, GpuBuffer, GpuView>
 *   - GPU non-owning (kernel arg):    \c GpuManager<S, FS, GpuView,   GpuView>
 */
template <class ScalarT,
          class FluidSystemT,
          template <class> class OuterStorage = ::Opm::VectorWithDefaultAllocator,
          template <class> class SampleStorage = ::Opm::VectorWithDefaultAllocator>
class GpuManager
    : private detail::SpecrockSampleHolder<
          ScalarT,
          std::is_same_v<OuterStorage<int>, ::Opm::gpuistl::GpuBuffer<int>>>
{
public:
    using Scalar = ScalarT;
    using FluidSystem = FluidSystemT;
    using SolidEnergyLawParams = ::Opm::EclSpecrockLawParams<ScalarT, SampleStorage>;
    using ThermalConductionLawParams = ::Opm::EclThconrLawParams<ScalarT>;

    static constexpr bool isOwningGpu
        = std::is_same_v<OuterStorage<int>, ::Opm::gpuistl::GpuBuffer<int>>;

    OPM_HOST_DEVICE GpuManager() = default;

    OPM_HOST_DEVICE
    GpuManager(OuterStorage<SolidEnergyLawParams> solidEnergyParams,
               OuterStorage<int> elementToSolidRegionIdx,
               OuterStorage<ThermalConductionLawParams> thermalConductionParams)
        : solidEnergyParams_(std::move(solidEnergyParams))
        , elementToSolidRegionIdx_(std::move(elementToSolidRegionIdx))
        , thermalConductionParams_(std::move(thermalConductionParams))
    {
    }

    /*!
     * \brief Owning constructor that also takes ownership of the
     *        per-region SPECROCK sample buffers (only available for the
     *        owning GpuBuffer outer-storage instantiation).
     */
    template <bool Enabled = isOwningGpu, std::enable_if_t<Enabled, int> = 0>
    GpuManager(std::vector<::Opm::gpuistl::GpuBuffer<Scalar>> sampleBuffers,
               OuterStorage<SolidEnergyLawParams> solidEnergyParams,
               OuterStorage<int> elementToSolidRegionIdx,
               OuterStorage<ThermalConductionLawParams> thermalConductionParams)
        : detail::SpecrockSampleHolder<Scalar, true>{std::move(sampleBuffers)}
        , solidEnergyParams_(std::move(solidEnergyParams))
        , elementToSolidRegionIdx_(std::move(elementToSolidRegionIdx))
        , thermalConductionParams_(std::move(thermalConductionParams))
    {
    }

    /*! \brief Solid-energy law parameters for an active cell. */
    OPM_HOST_DEVICE SolidEnergyLawParams solidEnergyLawParams(unsigned elemIdx) const
    {
        const int regionIdx = elementToSolidRegionIdx_[elemIdx];
        return solidEnergyParams_[regionIdx];
    }

    /*! \brief Thermal-conduction law parameters for an active cell. */
    OPM_HOST_DEVICE ThermalConductionLawParams
    thermalConductionLawParams(unsigned elemIdx) const
    {
        return thermalConductionParams_[elemIdx];
    }

    /*! \brief Direct storage accessors used by copy_to_gpu / make_view. */
    const OuterStorage<SolidEnergyLawParams>& solidEnergyParamsStorage() const
    {
        return solidEnergyParams_;
    }
    OuterStorage<SolidEnergyLawParams>& solidEnergyParamsStorage()
    {
        return solidEnergyParams_;
    }

    const OuterStorage<int>& elementToSolidRegionIdxStorage() const
    {
        return elementToSolidRegionIdx_;
    }
    OuterStorage<int>& elementToSolidRegionIdxStorage()
    {
        return elementToSolidRegionIdx_;
    }

    const OuterStorage<ThermalConductionLawParams>& thermalConductionParamsStorage() const
    {
        return thermalConductionParams_;
    }
    OuterStorage<ThermalConductionLawParams>& thermalConductionParamsStorage()
    {
        return thermalConductionParams_;
    }

    /*! \brief Mutable accessor for the per-region SPECROCK sample
     *         buffer holder; used by \c copy_to_gpu to populate the
     *         owning device buffers. Only available when this manager
     *         owns its device memory.
     */
    template <bool Enabled = isOwningGpu, std::enable_if_t<Enabled, int> = 0>
    std::vector<::Opm::gpuistl::GpuBuffer<Scalar>>& specrockSampleBuffers()
    {
        return this->sampleBuffers;
    }

private:
    OuterStorage<SolidEnergyLawParams> solidEnergyParams_ {};
    OuterStorage<int> elementToSolidRegionIdx_ {};
    OuterStorage<ThermalConductionLawParams> thermalConductionParams_ {};
};

/*!
 * \brief Build a CPU \c GpuManager from a CPU \c Opm::FlowProblem.
 *
 * Iterates over the active cells, querying the CPU problem for each
 * cell's \c SolidEnergyLawParams (multiplexer) and
 * \c ThermalConductionLawParams (multiplexer). The per-cell SPECROCK
 * multiplexer references are deduplicated by address to recover the
 * underlying region index (so the resulting GPU manager stores one
 * \c SolidEnergyLawParams per region rather than one per cell).
 *
 * Throws via \c OPM_THROW for any solid-energy approach other than
 * SPECROCK or any thermal-conduction approach other than THCONR.
 */
template <class ScalarT, class FluidSystemT, class CpuFlowProblemT>
GpuManager<ScalarT, FluidSystemT>
buildCpuManagerFromFlowProblem(const CpuFlowProblemT& cpu, std::size_t numElements)
{
    using ManagerCpu = GpuManager<ScalarT, FluidSystemT>;
    using SolidEnergyLawParams = typename ManagerCpu::SolidEnergyLawParams;
    using ThermalConductionLawParams = typename ManagerCpu::ThermalConductionLawParams;

    std::vector<const void*> regionMultiplexerAddress;
    std::vector<int> hostElementToSolidRegionIdx(numElements);
    std::vector<SolidEnergyLawParams> hostSolidEnergyParams;
    std::vector<ThermalConductionLawParams> hostThermalConductionParams(numElements);

    for (std::size_t i = 0; i < numElements; ++i) {
        const auto& cpuSolidMultiplexer
            = cpu.solidEnergyLawParams(static_cast<unsigned>(i), 0u);
        if (cpuSolidMultiplexer.solidEnergyApproach()
            != ::Opm::EclSolidEnergyApproach::Specrock) {
            OPM_THROW(std::logic_error,
                      "Opm::EclThermalLaw::GpuManager only supports the SPECROCK "
                      "solid-energy approach.");
        }
        const void* address = static_cast<const void*>(&cpuSolidMultiplexer);
        int regionIdx = -1;
        for (std::size_t r = 0; r < regionMultiplexerAddress.size(); ++r) {
            if (regionMultiplexerAddress[r] == address) {
                regionIdx = static_cast<int>(r);
                break;
            }
        }
        if (regionIdx < 0) {
            const auto& cpuSpecrock = cpuSolidMultiplexer.template getRealParams<
                ::Opm::EclSolidEnergyApproach::Specrock>();
            const auto& xValues = cpuSpecrock.temperatureSamples();
            const auto& yValues = cpuSpecrock.internalEnergySamples();
            std::vector<ScalarT> temperatureSamples(xValues.begin(), xValues.end());
            std::vector<ScalarT> internalEnergySamples(yValues.begin(), yValues.end());
            SolidEnergyLawParams params;
            params.setSamples(temperatureSamples, internalEnergySamples);
            regionIdx = static_cast<int>(hostSolidEnergyParams.size());
            hostSolidEnergyParams.emplace_back(std::move(params));
            regionMultiplexerAddress.push_back(address);
        }
        hostElementToSolidRegionIdx[i] = regionIdx;

        const auto& cpuThermalMultiplexer
            = cpu.thermalConductionLawParams(static_cast<unsigned>(i), 0u);
        if (cpuThermalMultiplexer.thermalConductionApproach()
            != ::Opm::EclThermalConductionApproach::Thconr) {
            OPM_THROW(std::logic_error,
                      "Opm::EclThermalLaw::GpuManager only supports the THCONR "
                      "thermal-conduction approach.");
        }
        const auto& cpuThconr = cpuThermalMultiplexer.template getRealParams<
            ::Opm::EclThermalConductionApproach::Thconr>();
        hostThermalConductionParams[i].setReferenceTotalThermalConductivity(
            cpuThconr.referenceTotalThermalConductivity());
        hostThermalConductionParams[i].setDTotalThermalConductivity_dSg(
            cpuThconr.dTotalThermalConductivity_dSg());
        hostThermalConductionParams[i].finalize();
    }

    return ManagerCpu(std::move(hostSolidEnergyParams),
                      std::move(hostElementToSolidRegionIdx),
                      std::move(hostThermalConductionParams));
}

} // namespace Opm::EclThermalLaw

namespace Opm::gpuistl {

/*!
 * \brief Copy a CPU \c GpuManager (plain host storage) to GPU-resident
 *        \c GpuBuffer storage.
 *
 * Each per-region SPECROCK sample table is uploaded to its own
 * \c GpuBuffer<Scalar>; those buffers are owned by the returned
 * manager's \c SpecrockSampleHolder base. The corresponding per-region
 * \c EclSpecrockLawParams<Scalar, GpuView> objects are then assembled
 * on the host and uploaded as a single bulk
 * \c GpuBuffer<SolidEnergyLawParams>.
 */
template <class ScalarT, class FluidSystemT>
::Opm::EclThermalLaw::GpuManager<ScalarT, FluidSystemT, GpuBuffer, GpuView>
copy_to_gpu(const ::Opm::EclThermalLaw::GpuManager<ScalarT, FluidSystemT>& cpu)
{
    using ManagerBuf
        = ::Opm::EclThermalLaw::GpuManager<ScalarT, FluidSystemT, GpuBuffer, GpuView>;
    using SolidEnergyLawParamsView = typename ManagerBuf::SolidEnergyLawParams;
    using ThermalConductionLawParams = typename ManagerBuf::ThermalConductionLawParams;

    std::vector<GpuBuffer<ScalarT>> sampleBuffers;
    sampleBuffers.reserve(2u * cpu.solidEnergyParamsStorage().size());

    std::vector<SolidEnergyLawParamsView> hostStaging;
    hostStaging.reserve(cpu.solidEnergyParamsStorage().size());

    for (const auto& cpuRegion : cpu.solidEnergyParamsStorage()) {
        sampleBuffers.emplace_back(GpuBuffer<ScalarT>(cpuRegion.temperatureSamples()));
        auto& tBuf = sampleBuffers.back();
        sampleBuffers.emplace_back(GpuBuffer<ScalarT>(cpuRegion.internalEnergySamples()));
        auto& eBuf = sampleBuffers.back();
        SolidEnergyLawParamsView params(GpuView<ScalarT>(tBuf.data(), tBuf.size()),
                                        GpuView<ScalarT>(eBuf.data(), eBuf.size()));
        hostStaging.emplace_back(std::move(params));
    }

    return ManagerBuf(std::move(sampleBuffers),
                      GpuBuffer<SolidEnergyLawParamsView>(hostStaging),
                      GpuBuffer<int>(cpu.elementToSolidRegionIdxStorage()),
                      GpuBuffer<ThermalConductionLawParams>(cpu.thermalConductionParamsStorage()));
}

/*!
 * \brief Make a non-owning \c GpuView based \c GpuManager from an owning
 *        \c GpuBuffer based \c GpuManager. The per-region solid-energy
 *        params element type is unchanged (\c GpuView sample storage in
 *        both cases).
 */
template <class ScalarT, class FluidSystemT>
::Opm::EclThermalLaw::GpuManager<ScalarT, FluidSystemT, GpuView, GpuView>
make_view(::Opm::EclThermalLaw::GpuManager<ScalarT, FluidSystemT, GpuBuffer, GpuView>& buf)
{
    using ManagerView
        = ::Opm::EclThermalLaw::GpuManager<ScalarT, FluidSystemT, GpuView, GpuView>;
    using SolidEnergyLawParamsView = typename ManagerView::SolidEnergyLawParams;
    using ThermalConductionLawParams = typename ManagerView::ThermalConductionLawParams;
    return ManagerView(GpuView<SolidEnergyLawParamsView>(buf.solidEnergyParamsStorage().data(),
                                                          buf.solidEnergyParamsStorage().size()),
                       GpuView<int>(buf.elementToSolidRegionIdxStorage().data(),
                                    buf.elementToSolidRegionIdxStorage().size()),
                       GpuView<ThermalConductionLawParams>(
                           buf.thermalConductionParamsStorage().data(),
                           buf.thermalConductionParamsStorage().size()));
}

} // namespace Opm::gpuistl

#endif // OPM_GPU_ECL_THERMAL_LAW_MANAGER_HPP
