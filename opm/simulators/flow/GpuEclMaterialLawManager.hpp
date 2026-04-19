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
 * FlowProblemBlackoil and the BlackOilIntensiveQuantities. It only supports
 * the EclDefaultMaterial three-phase law (no multiplexer, no hysteresis, no
 * directional/SWATINIT/PPCWMAX features).
 */
#ifndef OPM_GPU_ECL_MATERIAL_LAW_MANAGER_HPP
#define OPM_GPU_ECL_MATERIAL_LAW_MANAGER_HPP

#include <opm/common/utility/VectorWithDefaultAllocator.hpp>
#include <opm/common/utility/gpuDecorators.hpp>

#include <opm/material/fluidmatrixinteractions/EclDefaultMaterial.hpp>
#include <opm/material/fluidmatrixinteractions/EclDefaultMaterialParams.hpp>
#include <opm/material/fluidmatrixinteractions/EclMaterialLawManager.hpp>
#include <opm/material/fluidmatrixinteractions/EclMultiplexerMaterialParams.hpp>
#include <opm/material/fluidmatrixinteractions/EclTwoPhaseMaterial.hpp>
#include <opm/material/fluidmatrixinteractions/EclTwoPhaseMaterialParams.hpp>
#include <opm/material/fluidmatrixinteractions/PiecewiseLinearTwoPhaseMaterialParams.hpp>
#include <opm/material/fluidmatrixinteractions/SatCurveMultiplexerParams.hpp>

#include <opm/simulators/linalg/gpuistl/GpuBuffer.hpp>
#include <opm/simulators/linalg/gpuistl/GpuView.hpp>

#include <cassert>
#include <cstddef>
#include <memory>
#include <stdexcept>
#include <type_traits>
#include <utility>
#include <vector>

namespace Opm::EclMaterialLaw {

namespace detail {

/*!
 * \brief Walk down a CPU material-law parameter object until the enclosed
 *        \c PiecewiseLinearTwoPhaseMaterialParams is reached.
 *
 * Supports the standard layering used by the CPU
 * \c EclMaterialLaw::Manager:
 *   \c EclHysteresisTwoPhaseLawParams \c -> \c EclEpsTwoPhaseLawParams
 *   \c -> \c SatCurveMultiplexerParams \c -> \c PiecewiseLinearTwoPhaseMaterialParams.
 */
template <class CpuParams>
const auto& extractCpuPlParams(const CpuParams& params)
{
    if constexpr (requires { params.drainageParams(); }) {
        return extractCpuPlParams(params.drainageParams());
    }
    else if constexpr (requires { params.effectiveLawParams(); }) {
        return extractCpuPlParams(params.effectiveLawParams());
    }
    else if constexpr (requires {
                           params.template getRealParams<
                               ::Opm::SatCurveMultiplexerApproach::PiecewiseLinear>();
                       }) {
        return params.template getRealParams<
            ::Opm::SatCurveMultiplexerApproach::PiecewiseLinear>();
    }
    else {
        return params;
    }
}

/*!
 * \brief Build a vector-storage \c PiecewiseLinearTwoPhaseMaterialParams from
 *        any CPU two-phase parameter object that wraps a piecewise-linear law.
 */
template <class TargetParams, class CpuParams>
TargetParams makeVectorPlParams(const CpuParams& src)
{
    const auto& pl = extractCpuPlParams(src);
    using Vec = typename TargetParams::ValueVector;
    auto toVec = [](const auto& v) { return Vec(v.begin(), v.end()); };
    TargetParams out(toVec(pl.SwPcwnSamples()),
                     toVec(pl.pcwnSamples()),
                     toVec(pl.SwKrwSamples()),
                     toVec(pl.krwSamples()),
                     toVec(pl.SwKrnSamples()),
                     toVec(pl.krnSamples()));
    out.finalize();
    return out;
}

/*!
 * \brief Optional host-side base used by \c GpuManager when its storage is
 *        owning device memory: keeps the per-cell sample buffers alive so that
 *        the \c GpuView pointers stored inside the cells'
 *        \c PiecewiseLinearTwoPhaseMaterialParams remain valid for as long as
 *        the manager exists.  Empty for any non-owning storage so that
 *        GpuView-based managers stay device-trivially-copyable.
 */
template <class Scalar, bool Owning>
struct GpuPlSampleHolder {};

template <class Scalar>
struct GpuPlSampleHolder<Scalar, true> {
    std::vector<::Opm::gpuistl::GpuBuffer<Scalar>> plSampleBuffers_{};
};

/*!
 * \brief Compile-time trait: true if the given material-law parameter type
 *        looks like an \c EclTwoPhaseMaterialParams (i.e. has a nested
 *        \c GasWaterParams type).
 */
template <class T, class = void>
struct is_two_phase_mlp : std::false_type {};

template <class T>
struct is_two_phase_mlp<T, std::void_t<typename T::GasWaterParams>>
    : std::true_type {};

template <class T>
inline constexpr bool is_two_phase_mlp_v = is_two_phase_mlp<T>::value;

} // namespace detail

/*!
 * \brief Device-friendly stand-in for \c EclTwoPhaseMaterialParams.
 *
 * Stores all three sub-parameter objects by value (no \c std::shared_ptr) so
 * the whole object is trivially copyable to/from the device, exactly as
 * required by \c GpuView::operator[]. Provides the subset of the
 * \c EclTwoPhaseMaterialParams interface used by \c EclTwoPhaseMaterial.
 */
template <class TraitsT,
          class GasOilParamsT,
          class OilWaterParamsT,
          class GasWaterParamsT>
class GpuTwoPhaseMaterialParams
{
public:
    using Traits = TraitsT;
    using Scalar = typename Traits::Scalar;
    using GasOilParams = GasOilParamsT;
    using OilWaterParams = OilWaterParamsT;
    using GasWaterParams = GasWaterParamsT;

    OPM_HOST_DEVICE GpuTwoPhaseMaterialParams() = default;
    OPM_HOST_DEVICE GpuTwoPhaseMaterialParams(const GpuTwoPhaseMaterialParams&) = default;
    OPM_HOST_DEVICE GpuTwoPhaseMaterialParams(GpuTwoPhaseMaterialParams&&) = default;
    OPM_HOST_DEVICE GpuTwoPhaseMaterialParams& operator=(const GpuTwoPhaseMaterialParams&) = default;
    OPM_HOST_DEVICE GpuTwoPhaseMaterialParams& operator=(GpuTwoPhaseMaterialParams&&) = default;
    OPM_HOST_DEVICE ~GpuTwoPhaseMaterialParams() = default;

    OPM_HOST_DEVICE ::Opm::EclTwoPhaseApproach approach() const { return approach_; }
    void setApproach(::Opm::EclTwoPhaseApproach a) { approach_ = a; }

    OPM_HOST_DEVICE const GasOilParams& gasOilParams() const { return gasOilParams_; }
    OPM_HOST_DEVICE GasOilParams& gasOilParams() { return gasOilParams_; }
    void setGasOilParams(GasOilParams p) { gasOilParams_ = std::move(p); }

    OPM_HOST_DEVICE const OilWaterParams& oilWaterParams() const { return oilWaterParams_; }
    OPM_HOST_DEVICE OilWaterParams& oilWaterParams() { return oilWaterParams_; }
    void setOilWaterParams(OilWaterParams p) { oilWaterParams_ = std::move(p); }

    OPM_HOST_DEVICE const GasWaterParams& gasWaterParams() const { return gasWaterParams_; }
    OPM_HOST_DEVICE GasWaterParams& gasWaterParams() { return gasWaterParams_; }
    void setGasWaterParams(GasWaterParams p)
    {
        gasWaterParams_ = std::move(p);
        // Re-mark the inner sub-params as finalized: the assignment may
        // overwrite the EnsureFinalized base from a default-constructed
        // value. Calling EnsureFinalized::finalize directly only sets the
        // flag (no host-side touch of the GpuView samples).
        static_cast<::Opm::EnsureFinalized&>(gasWaterParams_).finalize();
    }

    void setSwl(Scalar) {}
    void finalize() {}

private:
    ::Opm::EclTwoPhaseApproach approach_{::Opm::EclTwoPhaseApproach::GasWater};
    GasOilParams gasOilParams_{};
    OilWaterParams oilWaterParams_{};
    GasWaterParams gasWaterParams_{};
};

/*!
 * \brief Device-friendly minimal two-phase material law.
 *
 * Mirrors the subset of \c EclTwoPhaseMaterial actually invoked by
 * \c BlackOilIntensiveQuantities (only \c capillaryPressures), and only
 * supports the GasWater sub-approach used by CO2STORE-style decks. All
 * methods are \c OPM_HOST_DEVICE so the law can be invoked from a kernel.
 */
template <class TraitsT,
          class GasOilMaterialLawT,
          class OilWaterMaterialLawT,
          class GasWaterMaterialLawT,
          class ParamsT = GpuTwoPhaseMaterialParams<TraitsT,
                                                    typename GasOilMaterialLawT::Params,
                                                    typename OilWaterMaterialLawT::Params,
                                                    typename GasWaterMaterialLawT::Params>>
class GpuTwoPhaseMaterial : public TraitsT
{
public:
    using GasOilMaterialLaw = GasOilMaterialLawT;
    using OilWaterMaterialLaw = OilWaterMaterialLawT;
    using GasWaterMaterialLaw = GasWaterMaterialLawT;

    using Traits = TraitsT;
    using Params = ParamsT;
    using Scalar = typename Traits::Scalar;

    static constexpr int numPhases = 3;
    static constexpr int waterPhaseIdx = Traits::wettingPhaseIdx;
    static constexpr int oilPhaseIdx = Traits::nonWettingPhaseIdx;
    static constexpr int gasPhaseIdx = Traits::gasPhaseIdx;

    static constexpr bool implementsTwoPhaseApi = false;
    static constexpr bool implementsTwoPhaseSatApi = false;
    static constexpr bool isSaturationDependent = true;
    static constexpr bool isPressureDependent = false;
    static constexpr bool isTemperatureDependent = false;
    static constexpr bool isCompositionDependent = false;

    template <class ContainerT, class FluidState, class... Args>
    OPM_HOST_DEVICE static void capillaryPressures(ContainerT& values,
                                                   const Params& params,
                                                   const FluidState& fluidState)
    {
        using Evaluation = typename std::remove_reference<decltype(values[0])>::type;
        const Evaluation& Sw = ::Opm::decay<Evaluation>(fluidState.saturation(waterPhaseIdx));
        values[waterPhaseIdx] = Evaluation(0.0);
        values[oilPhaseIdx] = Evaluation(0.0);
        values[gasPhaseIdx] =
            GasWaterMaterialLaw::template twoPhaseSatPcnw<Evaluation, Args...>(
                params.gasWaterParams(), Sw);
    }
};

/*!
 * \brief A minimal, GPU-compatible material-law manager.
 *
 * Only the EclDefaultMaterial three-phase law is supported. The two-phase
 * gas/oil and oil/water sub-law types are template parameters so the caller
 * can choose any GPU compatible implementation (e.g. PiecewiseLinear).
 *
 * \tparam TraitsT       Three-phase material traits.
 * \tparam GasOilLawT    Two-phase gas/oil law type (must implement the
 *                       saturation-only API).
 * \tparam OilWaterLawT  Two-phase oil/water law type (must implement the
 *                       saturation-only API).
 * \tparam Storage       Storage container template; defaults to a CPU vector.
 *                       Use Opm::gpuistl::GpuBuffer for owning GPU storage and
 *                       Opm::gpuistl::GpuView for non-owning GPU storage.
 */
template <class TraitsT,
          class GasOilLawT,
          class OilWaterLawT,
          template <class> class Storage = ::Opm::VectorWithDefaultAllocator,
          class MaterialLawT = ::Opm::EclDefaultMaterial<TraitsT, GasOilLawT, OilWaterLawT>>
class GpuManager
    : private detail::GpuPlSampleHolder<
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
     * Defaults to the three-phase \c EclDefaultMaterial. Pass
     * \c EclTwoPhaseMaterial<Traits, ..., ..., ...> to instead use the
     * two-phase law (only the GasWater sub-approach is currently supported by
     * the from-CPU constructors below; this matches the CO2STORE setup).
     */
    using MaterialLaw = MaterialLawT;
    using MaterialLawParams = typename MaterialLaw::Params;

    static constexpr int waterPhaseIdx = Traits::wettingPhaseIdx;
    static constexpr int oilPhaseIdx = Traits::nonWettingPhaseIdx;
    static constexpr int gasPhaseIdx = Traits::gasPhaseIdx;
    static constexpr int numPhases = Traits::numPhases;

    OPM_HOST_DEVICE GpuManager() = default;

    OPM_HOST_DEVICE GpuManager(Storage<MaterialLawParams> materialLawParams,
                               Storage<int> satnumRegionArray)
        : materialLawParams_(std::move(materialLawParams))
        , satnumRegionArray_(std::move(satnumRegionArray))
    {
    }

    /*!
     * \brief Construct from a fully initialised CPU
     *        \c Opm::EclMaterialLaw::Manager.
     *
     * For each cell the multiplexer parameter object is unwrapped down to
     * its inner \c PiecewiseLinearTwoPhaseMaterialParams; only the
     * \c EclMultiplexerApproach::Default approach is supported. The result
     * uses CPU vector storage and can subsequently be moved to the GPU via
     * \c Opm::gpuistl::copy_to_gpu followed by \c Opm::gpuistl::make_view.
     *
     * Only enabled when this manager itself uses CPU vector storage and the
     * two-phase laws use \c std::vector<Scalar> sample storage.
     */
    template <class CpuTraits = TraitsT,
              class CpuMgr = ::Opm::EclMaterialLaw::Manager<CpuTraits>,
              class GoP = GasOilParams,
              class OwP = OilWaterParams,
              class StI = Storage<int>,
              std::enable_if_t<
                  std::is_same_v<StI, ::Opm::VectorWithDefaultAllocator<int>>
                      && std::is_same_v<typename GoP::ValueVector, std::vector<Scalar>>
                      && std::is_same_v<typename OwP::ValueVector, std::vector<Scalar>>,
                  int> = 0>
    explicit GpuManager(const CpuMgr& cpu, std::size_t numElements)
    {
        materialLawParams_.resize(numElements);
        satnumRegionArray_.resize(numElements);

        for (std::size_t i = 0; i < numElements; ++i) {
            const auto& cpuMp = cpu.materialLawParams(static_cast<unsigned>(i));
            auto& mlp = materialLawParams_[i];

            if constexpr (detail::is_two_phase_mlp_v<MaterialLawParams>) {
                // Two-phase law (CO2STORE-style GasWater) -------------------
                if (cpuMp.approach() != ::Opm::EclMultiplexerApproach::TwoPhase) {
                    throw std::logic_error(
                        "Opm::EclMaterialLaw::GpuManager configured for the "
                        "two-phase material law requires a CPU material-law "
                        "parameter set that uses EclMultiplexerApproach::TwoPhase.");
                }
                const auto& cpuTp = cpuMp.template getRealParams<
                    ::Opm::EclMultiplexerApproach::TwoPhase>();
                if (cpuTp.approach() != ::Opm::EclTwoPhaseApproach::GasWater) {
                    throw std::logic_error(
                        "Opm::EclMaterialLaw::GpuManager only supports the "
                        "GasWater two-phase sub-approach.");
                }
                using GasWaterParams = typename MaterialLawParams::GasWaterParams;
                auto gwParams = detail::makeVectorPlParams<GasWaterParams>(
                    cpuTp.gasWaterParams());
                mlp.setGasWaterParams(std::move(gwParams));
                mlp.setApproach(::Opm::EclTwoPhaseApproach::GasWater);
                mlp.finalize();
            }
            else {
                // Default three-phase law ----------------------------------
                if (cpuMp.approach() != ::Opm::EclMultiplexerApproach::Default) {
                    throw std::logic_error(
                        "Opm::EclMaterialLaw::GpuManager only supports the Default "
                        "three-phase material law approach.");
                }
                const auto& cpuDef = cpuMp.template getRealParams<
                    ::Opm::EclMultiplexerApproach::Default>();

                auto goParams = detail::makeVectorPlParams<GasOilParams>(cpuDef.gasOilParams());
                auto owParams = detail::makeVectorPlParams<OilWaterParams>(cpuDef.oilWaterParams());

                mlp.gasOilParams() = std::move(goParams);
                mlp.oilWaterParams() = std::move(owParams);
                mlp.setSwl(cpuDef.Swl());
                mlp.finalize();
            }

            satnumRegionArray_[i] = static_cast<int>(
                cpu.satnumRegionIdx(static_cast<unsigned>(i)));
        }
    }

    /*!
     * \brief Construct from a CPU \c Opm::EclMaterialLaw::Manager directly
     *        into device-resident \c GpuBuffer storage.
     *
     * For each cell, the per-cell piecewise-linear sample arrays are uploaded
     * to the GPU as individual \c GpuBuffer<Scalar> instances kept alive by
     * this manager (via its private \c detail::GpuPlSampleHolder base).  The
     * cell's \c MaterialLawParams is populated with \c GpuView<const Scalar>
     * views referencing those device buffers and is then bulk-copied to the
     * device alongside the satnum array.
     *
     * Only enabled when this manager itself uses \c GpuBuffer storage and the
     * two-phase laws use \c GpuView<const Scalar> sample storage.
     */
    template <class CpuTraits = TraitsT,
              class CpuMgr = ::Opm::EclMaterialLaw::Manager<CpuTraits>,
              class GoP = GasOilParams,
              class OwP = OilWaterParams,
              class StI = Storage<int>,
              class StM = Storage<MaterialLawParams>,
              std::enable_if_t<
                  std::is_same_v<StI, ::Opm::gpuistl::GpuBuffer<int>>
                      && std::is_same_v<StM, ::Opm::gpuistl::GpuBuffer<MaterialLawParams>>
                      && std::is_same_v<typename GoP::ValueVector,
                                        ::Opm::gpuistl::GpuView<const Scalar>>
                      && std::is_same_v<typename OwP::ValueVector,
                                        ::Opm::gpuistl::GpuView<const Scalar>>,
                  int> = 0>
    explicit GpuManager(const CpuMgr& cpu, std::size_t numElements)
        : detail::GpuPlSampleHolder<Scalar, true>{}
        , materialLawParams_(buildHostMlpsForGpu(cpu, numElements,
                                                 this->plSampleBuffers_))
        , satnumRegionArray_(buildHostSatnum(cpu, numElements))
    {
    }

    /*! \brief Material-law parameters of an active cell. */
    OPM_HOST_DEVICE MaterialLawParams materialLawParams(unsigned elemIdx) const
    {
        // Return by value: GpuView::operator[] const already returns by value,
        // so binding a reference to materialLawParams_[idx] inside this function
        // would dangle. Returning by value lets the caller (which typically does
        // `const auto& mp = ...`) safely extend the temporary's lifetime.
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

    /*! \brief Direct access to the underlying storages (used by copy_to_gpu / make_view). */
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
     * \brief Build the per-cell host-side \c MaterialLawParams vector for the
     *        GpuBuffer-storage from-CPU constructor.  Allocates one
     *        \c GpuBuffer<Scalar> per per-cell sample array and wraps each in
     *        a \c GpuView<const Scalar> stored inside the cell's MLP. The
     *        owning device buffers are pushed into \p sampleBuffers, which is
     *        the manager's holder member.
     */
    template <class CpuMgr>
    static std::vector<MaterialLawParams>
    buildHostMlpsForGpu(const CpuMgr& cpu,
                        std::size_t numElements,
                        std::vector<::Opm::gpuistl::GpuBuffer<Scalar>>& sampleBuffers)
    {
        sampleBuffers.reserve(numElements * 12);
        std::vector<MaterialLawParams> mlps(numElements);

        auto pushBuf = [&](const auto& v) {
            std::vector<Scalar> tmp(v.begin(), v.end());
            sampleBuffers.emplace_back(tmp);
            const auto& back = sampleBuffers.back();
            return ::Opm::gpuistl::GpuView<const Scalar>(back.data(), back.size());
        };

        for (std::size_t i = 0; i < numElements; ++i) {
            const auto& cpuMp = cpu.materialLawParams(static_cast<unsigned>(i));

            if constexpr (detail::is_two_phase_mlp_v<MaterialLawParams>) {
                // Two-phase law (CO2STORE-style GasWater) -------------------
                if (cpuMp.approach() != ::Opm::EclMultiplexerApproach::TwoPhase) {
                    throw std::logic_error(
                        "Opm::EclMaterialLaw::GpuManager configured for the "
                        "two-phase material law requires a CPU material-law "
                        "parameter set that uses EclMultiplexerApproach::TwoPhase.");
                }
                const auto& cpuTp = cpuMp.template getRealParams<
                    ::Opm::EclMultiplexerApproach::TwoPhase>();
                if (cpuTp.approach() != ::Opm::EclTwoPhaseApproach::GasWater) {
                    throw std::logic_error(
                        "Opm::EclMaterialLaw::GpuManager only supports the "
                        "GasWater two-phase sub-approach.");
                }
                using GasWaterParams = typename MaterialLawParams::GasWaterParams;
                const auto& gwPl = detail::extractCpuPlParams(cpuTp.gasWaterParams());
                GasWaterParams gwParams(pushBuf(gwPl.SwPcwnSamples()),
                                        pushBuf(gwPl.pcwnSamples()),
                                        pushBuf(gwPl.SwKrwSamples()),
                                        pushBuf(gwPl.krwSamples()),
                                        pushBuf(gwPl.SwKrnSamples()),
                                        pushBuf(gwPl.krnSamples()));
                // The GpuView-backed ctor of PiecewiseLinearTwoPhaseMaterialParams
                // already invokes EnsureFinalized::finalize(); calling the full
                // finalize() here would dereference the device samples on host.
                mlps[i].setGasWaterParams(std::move(gwParams));
                mlps[i].setApproach(::Opm::EclTwoPhaseApproach::GasWater);
                mlps[i].finalize();
            }
            else {
                // Default three-phase law ----------------------------------
                if (cpuMp.approach() != ::Opm::EclMultiplexerApproach::Default) {
                    throw std::logic_error(
                        "Opm::EclMaterialLaw::GpuManager only supports the Default "
                        "three-phase material law approach.");
                }
                const auto& cpuDef = cpuMp.template getRealParams<
                    ::Opm::EclMultiplexerApproach::Default>();
                const auto& goPl = detail::extractCpuPlParams(cpuDef.gasOilParams());
                const auto& owPl = detail::extractCpuPlParams(cpuDef.oilWaterParams());

                GasOilParams goParams(pushBuf(goPl.SwPcwnSamples()),
                                      pushBuf(goPl.pcwnSamples()),
                                      pushBuf(goPl.SwKrwSamples()),
                                      pushBuf(goPl.krwSamples()),
                                      pushBuf(goPl.SwKrnSamples()),
                                      pushBuf(goPl.krnSamples()));
                OilWaterParams owParams(pushBuf(owPl.SwPcwnSamples()),
                                        pushBuf(owPl.pcwnSamples()),
                                        pushBuf(owPl.SwKrwSamples()),
                                        pushBuf(owPl.krwSamples()),
                                        pushBuf(owPl.SwKrnSamples()),
                                        pushBuf(owPl.krnSamples()));

                mlps[i].setGasOilParams(std::make_shared<GasOilParams>(goParams));
                mlps[i].setOilWaterParams(std::make_shared<OilWaterParams>(owParams));
                mlps[i].setSwl(cpuDef.Swl());
                mlps[i].finalize();
            }
        }
        return mlps;
    }

    template <class CpuMgr>
    static std::vector<int>
    buildHostSatnum(const CpuMgr& cpu, std::size_t numElements)
    {
        std::vector<int> satnum(numElements);
        for (std::size_t i = 0; i < numElements; ++i) {
            satnum[i] = static_cast<int>(
                cpu.satnumRegionIdx(static_cast<unsigned>(i)));
        }
        return satnum;
    }

    Storage<MaterialLawParams> materialLawParams_{};
    Storage<int> satnumRegionArray_{};
};

} // namespace Opm::EclMaterialLaw

namespace Opm::gpuistl {

/*!
 * \brief Copy a CPU GpuManager to GPU-resident GpuBuffer storage.
 *
 * The MaterialLawParams element type is assumed to be the same on the CPU and
 * the GPU, i.e. the caller is responsible for making the GasOilLaw /
 * OilWaterLaw GPU-compatible (typically by templating their parameter type on
 * a GPU storage).
 */
template <class TraitsT,
          class GasOilLawT,
          class OilWaterLawT,
          class MaterialLawT>
::Opm::EclMaterialLaw::GpuManager<TraitsT, GasOilLawT, OilWaterLawT, GpuBuffer, MaterialLawT>
copy_to_gpu(const ::Opm::EclMaterialLaw::GpuManager<TraitsT, GasOilLawT, OilWaterLawT,
                                                    ::Opm::VectorWithDefaultAllocator,
                                                    MaterialLawT>& cpu)
{
    using ManagerGpu = ::Opm::EclMaterialLaw::GpuManager<TraitsT, GasOilLawT, OilWaterLawT,
                                                         GpuBuffer, MaterialLawT>;
    using MLP = typename ManagerGpu::MaterialLawParams;
    return ManagerGpu(GpuBuffer<MLP>(cpu.materialLawParamsStorage()),
                      GpuBuffer<int>(cpu.satnumRegionArrayStorage()));
}

/*!
 * \brief Make a non-owning GpuView based GpuManager from an owning GpuBuffer
 *        based GpuManager.
 */
template <class TraitsT,
          class GasOilLawT,
          class OilWaterLawT,
          class MaterialLawT>
::Opm::EclMaterialLaw::GpuManager<TraitsT, GasOilLawT, OilWaterLawT, GpuView, MaterialLawT>
make_view(::Opm::EclMaterialLaw::GpuManager<TraitsT, GasOilLawT, OilWaterLawT, GpuBuffer,
                                            MaterialLawT>& buf)
{
    using ManagerView = ::Opm::EclMaterialLaw::GpuManager<TraitsT, GasOilLawT, OilWaterLawT,
                                                          GpuView, MaterialLawT>;
    using MLP = typename ManagerView::MaterialLawParams;
    return ManagerView(GpuView<MLP>(buf.materialLawParamsStorage().data(),
                                    buf.materialLawParamsStorage().size()),
                       GpuView<int>(buf.satnumRegionArrayStorage().data(),
                                    buf.satnumRegionArrayStorage().size()));
}

} // namespace Opm::gpuistl

#endif // OPM_GPU_ECL_MATERIAL_LAW_MANAGER_HPP
