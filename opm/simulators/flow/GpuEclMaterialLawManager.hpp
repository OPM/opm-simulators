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
 * \brief GPU-compatible, simplified version of Opm::EclMaterialLaw::Manager.
 *
 * This class only supports the actions needed by FlowProblem,
 * FlowProblemBlackoil and the BlackOilIntensiveQuantities. It re-uses the
 * generic two-phase material types from opm-common
 * (\c Opm::EclTwoPhaseMaterial and \c Opm::EclTwoPhaseMaterialParams) with
 * a value-storage policy (\c Opm::gpuistl::ValueAsPointer) so the per-cell
 * material law parameters are trivially copyable to the device.
 *
 * The chain of multiplexers used by the CPU manager
 * (\c EclMultiplexerMaterial \f$\to\f$ \c EclEpsTwoPhaseLaw \f$\to\f$
 * \c EclHysteresisTwoPhaseLaw \f$\to\f$ \c SatCurveMultiplexer) is bypassed:
 * the GPU manager unwraps the CPU multiplexer parameters down to the
 * underlying \c PiecewiseLinearTwoPhaseMaterialParams and uploads only the
 * piecewise-linear sample tables.
 */
#ifndef OPM_GPU_ECL_MATERIAL_LAW_MANAGER_HPP
#define OPM_GPU_ECL_MATERIAL_LAW_MANAGER_HPP

#include <opm/common/utility/VectorWithDefaultAllocator.hpp>
#include <opm/common/utility/gpuDecorators.hpp>

#include <opm/material/fluidmatrixinteractions/EclMaterialLawManager.hpp>
#include <opm/material/fluidmatrixinteractions/EclMultiplexerMaterialParams.hpp>
#include <opm/material/fluidmatrixinteractions/EclTwoPhaseMaterial.hpp>
#include <opm/material/fluidmatrixinteractions/EclTwoPhaseMaterialParams.hpp>
#include <opm/material/fluidmatrixinteractions/PiecewiseLinearTwoPhaseMaterialParams.hpp>
#include <opm/material/fluidmatrixinteractions/SatCurveMultiplexerParams.hpp>

#include <opm/simulators/linalg/gpuistl/GpuBuffer.hpp>
#include <opm/simulators/linalg/gpuistl/GpuView.hpp>
#include <opm/simulators/linalg/gpuistl/gpu_smart_pointer.hpp>

#include <cassert>
#include <cstddef>
#include <stdexcept>
#include <type_traits>
#include <utility>
#include <vector>

namespace Opm::EclMaterialLaw
{

namespace detail
{

    /*!
     * \brief Walk down a CPU material-law parameter object until the enclosed
     *        \c PiecewiseLinearTwoPhaseMaterialParams is reached.
     *
     * Supports the standard layering used by the CPU \c EclMaterialLaw::Manager:
     *   \c EclHysteresisTwoPhaseLawParams \f$\to\f$ \c EclEpsTwoPhaseLawParams
     *   \f$\to\f$ \c SatCurveMultiplexerParams
     *   \f$\to\f$ \c PiecewiseLinearTwoPhaseMaterialParams.
     */
    template <class CpuParams>
    const auto& extractCpuPlParams(const CpuParams& params)
    {
        if constexpr (requires { params.drainageParams(); }) {
            return extractCpuPlParams(params.drainageParams());
        } else if constexpr (requires { params.effectiveLawParams(); }) {
            return extractCpuPlParams(params.effectiveLawParams());
        } else if constexpr (requires {
                                 params.template getRealParams<
                                     ::Opm::SatCurveMultiplexerApproach::PiecewiseLinear>();
                             }) {
            return params
                .template getRealParams<::Opm::SatCurveMultiplexerApproach::PiecewiseLinear>();
        } else {
            return params;
        }
    }

    /*!
     * \brief Optional host-side base used by \c GpuManager when its storage is
     *        owning device memory: keeps the per-cell sample buffers alive so
     *        that the \c GpuView pointers stored inside the cells'
     *        \c PiecewiseLinearTwoPhaseMaterialParams remain valid for as long
     *        as the manager exists. Empty for any non-owning storage so that
     *        \c GpuView based managers stay device-trivially-copyable.
     */
    template <class Scalar, bool Owning>
    struct GpuPiecewiseLinearSampleHolder {
    };

    template <class Scalar>
    struct GpuPiecewiseLinearSampleHolder<Scalar, true> {
        std::vector<::Opm::gpuistl::GpuBuffer<Scalar>> sampleBuffers {};
    };

} // namespace detail

/*!
 * \brief A minimal, GPU-compatible material-law manager.
 *
 * Only the EclDefaultMaterial three-phase law and the EclTwoPhaseMaterial
 * two-phase law are supported. The two-phase gas/oil and oil/water sub-law
 * types are template parameters so the caller can choose any GPU compatible
 * implementation (e.g. PiecewiseLinear).
 *
 * \tparam TraitsT       Three-phase material traits.
 * \tparam GasOilLawT    Two-phase gas/oil law type (must implement the
 *                       saturation-only API).
 * \tparam OilWaterLawT  Two-phase oil/water law type (must implement the
 *                       saturation-only API).
 * \tparam Storage       Storage container template; defaults to a CPU vector.
 *                       Use \c Opm::gpuistl::GpuBuffer for owning GPU storage
 *                       and \c Opm::gpuistl::GpuView for non-owning GPU
 *                       storage.
 * \tparam MaterialLawT  The three-phase material-law type to use; defaults to
 *                       \c Opm::EclDefaultMaterial.
 */
template <class TraitsT,
          class GasOilLawT,
          class OilWaterLawT,
          template <class> class Storage = ::Opm::VectorWithDefaultAllocator,
          class MaterialLawT = ::Opm::EclDefaultMaterial<TraitsT, GasOilLawT, OilWaterLawT>>
class GpuManager : private detail::GpuPiecewiseLinearSampleHolder<
                       typename TraitsT::Scalar,
                       std::is_same_v<Storage<int>, ::Opm::gpuistl::GpuBuffer<int>>>
{
public:
    using Traits = TraitsT;
    using Scalar = typename Traits::Scalar;
    using GasOilLaw = GasOilLawT;
    using OilWaterLaw = OilWaterLawT;
    using GasOilParams = typename GasOilLaw::Params;
    using OilWaterParams = typename OilWaterLaw::Params;

    /*!
     * \brief The actual material-law type used by the IQ.
     *
     * Defaults to the three-phase \c EclDefaultMaterial. Pass an
     * \c Opm::EclTwoPhaseMaterial<...> instantiation (with a value-storage
     * \c EclTwoPhaseMaterialParams, e.g. one using
     * \c Opm::gpuistl::ValueAsPointer) to instead use the two-phase law.
     * Only the GasWater sub-approach is currently supported by the
     * from-CPU constructors below; this matches the CO2STORE setup.
     */
    using MaterialLaw = MaterialLawT;
    using MaterialLawParams = typename MaterialLaw::Params;

    static constexpr int waterPhaseIdx = Traits::wettingPhaseIdx;
    static constexpr int oilPhaseIdx = Traits::nonWettingPhaseIdx;
    static constexpr int gasPhaseIdx = Traits::gasPhaseIdx;
    static constexpr int numPhases = Traits::numPhases;

private:
    /*! \brief Trait: true if the material-law parameters are an
     *         \c EclTwoPhaseMaterialParams (i.e. expose a nested
     *         \c GasWaterParams type). */
    template <class T, class = void>
    struct IsTwoPhaseMaterialLawParams : std::false_type {
    };

    template <class T>
    struct IsTwoPhaseMaterialLawParams<T, std::void_t<typename T::GasWaterParams>>
        : std::true_type {
    };

    static constexpr bool isTwoPhase = IsTwoPhaseMaterialLawParams<MaterialLawParams>::value;

public:
    OPM_HOST_DEVICE GpuManager() = default;

    OPM_HOST_DEVICE GpuManager(Storage<MaterialLawParams> materialLawParams,
                               Storage<int> satnumRegionArray)
        : materialLawParams_(std::move(materialLawParams))
        , satnumRegionArray_(std::move(satnumRegionArray))
    {
    }

    /*!
     * \brief Construct from a CPU \c Opm::EclMaterialLaw::Manager directly
     *        into device-resident \c GpuBuffer storage.
     *
     * For each cell, the per-cell piecewise-linear sample arrays are
     * uploaded to the GPU as individual \c GpuBuffer<Scalar> instances kept
     * alive by this manager (via its private
     * \c detail::GpuPiecewiseLinearSampleHolder base). The cell's
     * \c MaterialLawParams is populated with \c GpuView<const Scalar> views
     * referencing those device buffers and is then bulk-copied to the
     * device alongside the satnum array.
     *
     * Only enabled when this manager itself uses \c GpuBuffer storage and
     * the two-phase laws use \c GpuView<const Scalar> sample storage.
     */
    template <class CpuTraits = TraitsT,
              class CpuManager = ::Opm::EclMaterialLaw::Manager<CpuTraits>,
              class GasOilParamsArg = GasOilParams,
              class OilWaterParamsArg = OilWaterParams,
              class IntegerStorage = Storage<int>,
              class MaterialLawParamsStorageT = Storage<MaterialLawParams>,
              std::enable_if_t<std::is_same_v<IntegerStorage, ::Opm::gpuistl::GpuBuffer<int>>
                                   && std::is_same_v<MaterialLawParamsStorageT,
                                                     ::Opm::gpuistl::GpuBuffer<MaterialLawParams>>
                                   && std::is_same_v<typename GasOilParamsArg::ValueVector,
                                                     ::Opm::gpuistl::GpuView<const Scalar>>
                                   && std::is_same_v<typename OilWaterParamsArg::ValueVector,
                                                     ::Opm::gpuistl::GpuView<const Scalar>>,
                               int>
              = 0>
    explicit GpuManager(const CpuManager& cpu, std::size_t numElements)
        : detail::GpuPiecewiseLinearSampleHolder<Scalar, true> {}
        , materialLawParams_(
              buildHostMaterialLawParamsForGpu(cpu, numElements, this->sampleBuffers))
        , satnumRegionArray_(buildHostSatnumRegionArray(cpu, numElements))
    {
    }

    /*! \brief Material-law parameters of an active cell. */
    OPM_HOST_DEVICE MaterialLawParams materialLawParams(unsigned elemIdx) const
    {
        // Return by value: GpuView::operator[] const already returns by
        // value, so binding a reference to materialLawParams_[idx] inside
        // this function would dangle. Returning by value lets the caller
        // (which typically does `const auto& mp = ...`) safely extend the
        // temporary's lifetime.
        return materialLawParams_[elemIdx];
    }

    OPM_HOST_DEVICE MaterialLawParams& materialLawParams(unsigned elemIdx)
    {
        return materialLawParams_[elemIdx];
    }

    /*! \brief Saturation-region index of an active cell. */
    OPM_HOST_DEVICE int satnumRegionIdx(unsigned elemIdx) const
    {
        return satnumRegionArray_[elemIdx];
    }

    /*! \brief Direct access to the underlying storages
     *         (used by copy_to_gpu / make_view). */
    const Storage<MaterialLawParams>& materialLawParamsStorage() const
    {
        return materialLawParams_;
    }
    Storage<MaterialLawParams>& materialLawParamsStorage()
    {
        return materialLawParams_;
    }

    const Storage<int>& satnumRegionArrayStorage() const
    {
        return satnumRegionArray_;
    }
    Storage<int>& satnumRegionArrayStorage()
    {
        return satnumRegionArray_;
    }

private:
    /*!
     * \brief Build the per-cell host-side \c MaterialLawParams vector for
     *        the GpuBuffer-storage from-CPU constructor. Allocates one
     *        \c GpuBuffer<Scalar> per per-cell sample array and wraps each
     *        in a \c GpuView<const Scalar> stored inside the cell's MLP.
     *        The owning device buffers are pushed into \p sampleBuffers,
     *        which is the manager's holder member.
     */
    template <class CpuManager>
    static std::vector<MaterialLawParams>
    buildHostMaterialLawParamsForGpu(const CpuManager& cpu,
                                     std::size_t numElements,
                                     std::vector<::Opm::gpuistl::GpuBuffer<Scalar>>& sampleBuffers)
    {
        // The CPU manager produces one piecewise-linear sample table per
        // sub-law per cell. For a three-phase EclDefaultMaterial that is up
        // to twelve buffers per cell (six per gas/oil and oil/water sub-law);
        // for the two-phase EclTwoPhaseMaterial it is six (a single
        // gas/water sub-law). Reserving the upper bound avoids reallocation.
        sampleBuffers.reserve(numElements * 12);
        std::vector<MaterialLawParams> materialLawParams(numElements);

        auto pushSampleBuffer = [&](const auto& sampleVector) {
            std::vector<Scalar> hostCopy(sampleVector.begin(), sampleVector.end());
            sampleBuffers.emplace_back(hostCopy);
            const auto& deviceBuffer = sampleBuffers.back();
            return ::Opm::gpuistl::GpuView<const Scalar>(deviceBuffer.data(), deviceBuffer.size());
        };

        for (std::size_t i = 0; i < numElements; ++i) {
            const auto& cpuMaterialParams = cpu.materialLawParams(static_cast<unsigned>(i));

            if constexpr (isTwoPhase) {
                buildTwoPhaseCellParams(cpuMaterialParams, materialLawParams[i], pushSampleBuffer);
            } else {
                buildThreePhaseCellParams(
                    cpuMaterialParams, materialLawParams[i], pushSampleBuffer);
            }
        }
        return materialLawParams;
    }

    /*!
     * \brief Build a single cell's two-phase \c EclTwoPhaseMaterialParams
     *        instance with GPU-resident piecewise-linear sample tables.
     *
     * Only the \c GasWater sub-approach is currently supported, since this
     * is the only configuration produced by CO2STORE-style decks.
     */
    template <class CpuMaterialParams, class PushBuffer>
    static void buildTwoPhaseCellParams(const CpuMaterialParams& cpuMaterialParams,
                                        MaterialLawParams& cellParams,
                                        PushBuffer& pushSampleBuffer)
    {
        if (cpuMaterialParams.approach() != ::Opm::EclMultiplexerApproach::TwoPhase) {
            throw std::logic_error("Opm::EclMaterialLaw::GpuManager configured for the two-phase "
                                   "material law requires a CPU material-law parameter set that "
                                   "uses EclMultiplexerApproach::TwoPhase.");
        }
        const auto& cpuTwoPhaseParams
            = cpuMaterialParams.template getRealParams<::Opm::EclMultiplexerApproach::TwoPhase>();
        if (cpuTwoPhaseParams.approach() != ::Opm::EclTwoPhaseApproach::GasWater) {
            throw std::logic_error("Opm::EclMaterialLaw::GpuManager only supports the GasWater "
                                   "two-phase sub-approach.");
        }

        using GasWaterParams = typename MaterialLawParams::GasWaterParams;
        const auto& cpuGasWaterPiecewiseLinear
            = detail::extractCpuPlParams(cpuTwoPhaseParams.gasWaterParams());
        GasWaterParams gasWaterParams(pushSampleBuffer(cpuGasWaterPiecewiseLinear.SwPcwnSamples()),
                                      pushSampleBuffer(cpuGasWaterPiecewiseLinear.pcwnSamples()),
                                      pushSampleBuffer(cpuGasWaterPiecewiseLinear.SwKrwSamples()),
                                      pushSampleBuffer(cpuGasWaterPiecewiseLinear.krwSamples()),
                                      pushSampleBuffer(cpuGasWaterPiecewiseLinear.SwKrnSamples()),
                                      pushSampleBuffer(cpuGasWaterPiecewiseLinear.krnSamples()));

        cellParams.setGasWaterParams(
            typename MaterialLawParams::GasWaterParamsStorage(std::move(gasWaterParams)));
        cellParams.setApproach(::Opm::EclTwoPhaseApproach::GasWater);
        cellParams.finalize();
    }

    /*!
     * \brief Build a single cell's three-phase \c EclDefaultMaterialParams
     *        instance with GPU-resident piecewise-linear sample tables.
     */
    template <class CpuMaterialParams, class PushBuffer>
    static void buildThreePhaseCellParams(const CpuMaterialParams& cpuMaterialParams,
                                          MaterialLawParams& cellParams,
                                          PushBuffer& pushSampleBuffer)
    {
        if (cpuMaterialParams.approach() != ::Opm::EclMultiplexerApproach::Default) {
            throw std::logic_error("Opm::EclMaterialLaw::GpuManager only supports the Default "
                                   "three-phase material law approach.");
        }
        const auto& cpuDefaultParams
            = cpuMaterialParams.template getRealParams<::Opm::EclMultiplexerApproach::Default>();
        const auto& cpuGasOilPiecewiseLinear
            = detail::extractCpuPlParams(cpuDefaultParams.gasOilParams());
        const auto& cpuOilWaterPiecewiseLinear
            = detail::extractCpuPlParams(cpuDefaultParams.oilWaterParams());

        GasOilParams gasOilParams(pushSampleBuffer(cpuGasOilPiecewiseLinear.SwPcwnSamples()),
                                  pushSampleBuffer(cpuGasOilPiecewiseLinear.pcwnSamples()),
                                  pushSampleBuffer(cpuGasOilPiecewiseLinear.SwKrwSamples()),
                                  pushSampleBuffer(cpuGasOilPiecewiseLinear.krwSamples()),
                                  pushSampleBuffer(cpuGasOilPiecewiseLinear.SwKrnSamples()),
                                  pushSampleBuffer(cpuGasOilPiecewiseLinear.krnSamples()));
        OilWaterParams oilWaterParams(pushSampleBuffer(cpuOilWaterPiecewiseLinear.SwPcwnSamples()),
                                      pushSampleBuffer(cpuOilWaterPiecewiseLinear.pcwnSamples()),
                                      pushSampleBuffer(cpuOilWaterPiecewiseLinear.SwKrwSamples()),
                                      pushSampleBuffer(cpuOilWaterPiecewiseLinear.krwSamples()),
                                      pushSampleBuffer(cpuOilWaterPiecewiseLinear.SwKrnSamples()),
                                      pushSampleBuffer(cpuOilWaterPiecewiseLinear.krnSamples()));

        cellParams.setGasOilParams(std::make_shared<GasOilParams>(std::move(gasOilParams)));
        cellParams.setOilWaterParams(std::make_shared<OilWaterParams>(std::move(oilWaterParams)));
        cellParams.setSwl(cpuDefaultParams.Swl());
        cellParams.finalize();
    }

    template <class CpuManager>
    static std::vector<int> buildHostSatnumRegionArray(const CpuManager& cpu,
                                                       std::size_t numElements)
    {
        std::vector<int> satnumRegionArray(numElements);
        for (std::size_t i = 0; i < numElements; ++i) {
            satnumRegionArray[i] = static_cast<int>(cpu.satnumRegionIdx(static_cast<unsigned>(i)));
        }
        return satnumRegionArray;
    }

    Storage<MaterialLawParams> materialLawParams_ {};
    Storage<int> satnumRegionArray_ {};
};

} // namespace Opm::EclMaterialLaw

namespace Opm::gpuistl
{

/*!
 * \brief Copy a CPU GpuManager to GPU-resident GpuBuffer storage.
 *
 * The MaterialLawParams element type is assumed to be the same on the CPU
 * and the GPU, i.e. the caller is responsible for making the GasOilLaw and
 * OilWaterLaw GPU-compatible (typically by templating their parameter type
 * on a GPU storage).
 */
template <class TraitsT, class GasOilLawT, class OilWaterLawT, class MaterialLawT>
::Opm::EclMaterialLaw::GpuManager<TraitsT, GasOilLawT, OilWaterLawT, GpuBuffer, MaterialLawT>
copy_to_gpu(const ::Opm::EclMaterialLaw::GpuManager<TraitsT,
                                                    GasOilLawT,
                                                    OilWaterLawT,
                                                    ::Opm::VectorWithDefaultAllocator,
                                                    MaterialLawT>& cpu)
{
    using GpuManagerBuffer = ::Opm::EclMaterialLaw::
        GpuManager<TraitsT, GasOilLawT, OilWaterLawT, GpuBuffer, MaterialLawT>;
    using MaterialLawParams = typename GpuManagerBuffer::MaterialLawParams;
    return GpuManagerBuffer(GpuBuffer<MaterialLawParams>(cpu.materialLawParamsStorage()),
                            GpuBuffer<int>(cpu.satnumRegionArrayStorage()));
}

/*!
 * \brief Make a non-owning GpuView based GpuManager from an owning GpuBuffer
 *        based GpuManager.
 */
template <class TraitsT, class GasOilLawT, class OilWaterLawT, class MaterialLawT>
::Opm::EclMaterialLaw::GpuManager<TraitsT, GasOilLawT, OilWaterLawT, GpuView, MaterialLawT>
make_view(
    ::Opm::EclMaterialLaw::GpuManager<TraitsT, GasOilLawT, OilWaterLawT, GpuBuffer, MaterialLawT>&
        buf)
{
    using GpuManagerView = ::Opm::EclMaterialLaw::
        GpuManager<TraitsT, GasOilLawT, OilWaterLawT, GpuView, MaterialLawT>;
    using MaterialLawParams = typename GpuManagerView::MaterialLawParams;
    return GpuManagerView(
        GpuView<MaterialLawParams>(buf.materialLawParamsStorage().data(),
                                   buf.materialLawParamsStorage().size()),
        GpuView<int>(buf.satnumRegionArrayStorage().data(), buf.satnumRegionArrayStorage().size()));
}

} // namespace Opm::gpuistl

#endif // OPM_GPU_ECL_MATERIAL_LAW_MANAGER_HPP
