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
#include <opm/material/fluidmatrixinteractions/PiecewiseLinearTwoPhaseMaterialParams.hpp>
#include <opm/material/fluidmatrixinteractions/SatCurveMultiplexerParams.hpp>

#include <cassert>
#include <cstddef>
#include <stdexcept>
#include <type_traits>
#include <utility>
#include <vector>

namespace Opm::gpuistl
{
template <typename T>
class GpuBuffer;
template <typename T>
class GpuView;
} // namespace Opm::gpuistl

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
    return TargetParams(toVec(pl.SwPcwnSamples()),
                        toVec(pl.pcwnSamples()),
                        toVec(pl.SwKrwSamples()),
                        toVec(pl.krwSamples()),
                        toVec(pl.SwKrnSamples()),
                        toVec(pl.krnSamples()));
}

} // namespace detail

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
          template <class> class Storage = ::Opm::VectorWithDefaultAllocator>
class GpuManager
{
public:
    using Traits = TraitsT;
    using Scalar = typename Traits::Scalar;
    using GasOilLaw = GasOilLawT;
    using OilWaterLaw = OilWaterLawT;
    using GasOilParams = typename GasOilLaw::Params;
    using OilWaterParams = typename OilWaterLaw::Params;

    // The three-phase material law and its parameter pack. The default
    // EclDefaultMaterialParams already exposes Storage-compatible storage
    // because it stores its sub-parameter objects by value.
    using MaterialLaw = ::Opm::EclDefaultMaterial<Traits, GasOilLaw, OilWaterLaw>;
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
            if (cpuMp.approach() != ::Opm::EclMultiplexerApproach::Default) {
                throw std::logic_error(
                    "Opm::EclMaterialLaw::GpuManager only supports the Default "
                    "three-phase material law approach.");
            }
            const auto& cpuDef = cpuMp.template getRealParams<
                ::Opm::EclMultiplexerApproach::Default>();

            auto goParams = detail::makeVectorPlParams<GasOilParams>(cpuDef.gasOilParams());
            auto owParams = detail::makeVectorPlParams<OilWaterParams>(cpuDef.oilWaterParams());

            auto& mlp = materialLawParams_[i];
            mlp.gasOilParams() = std::move(goParams);
            mlp.oilWaterParams() = std::move(owParams);
            mlp.setSwl(cpuDef.Swl());
            mlp.finalize();

            satnumRegionArray_[i] = static_cast<int>(
                cpu.satnumRegionIdx(static_cast<unsigned>(i)));
        }
    }

    /*! \brief Material-law parameters of an active cell. */
    OPM_HOST_DEVICE const MaterialLawParams& materialLawParams(unsigned elemIdx) const
    {
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
          class OilWaterLawT>
::Opm::EclMaterialLaw::GpuManager<TraitsT, GasOilLawT, OilWaterLawT, GpuBuffer>
copy_to_gpu(const ::Opm::EclMaterialLaw::GpuManager<TraitsT, GasOilLawT, OilWaterLawT>& cpu)
{
    using ManagerGpu = ::Opm::EclMaterialLaw::GpuManager<TraitsT, GasOilLawT, OilWaterLawT, GpuBuffer>;
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
          class OilWaterLawT>
::Opm::EclMaterialLaw::GpuManager<TraitsT, GasOilLawT, OilWaterLawT, GpuView>
make_view(::Opm::EclMaterialLaw::GpuManager<TraitsT, GasOilLawT, OilWaterLawT, GpuBuffer>& buf)
{
    using ManagerView = ::Opm::EclMaterialLaw::GpuManager<TraitsT, GasOilLawT, OilWaterLawT, GpuView>;
    using MLP = typename ManagerView::MaterialLawParams;
    return ManagerView(GpuView<MLP>(buf.materialLawParamsStorage().data(),
                                    buf.materialLawParamsStorage().size()),
                       GpuView<int>(buf.satnumRegionArrayStorage().data(),
                                    buf.satnumRegionArrayStorage().size()));
}

} // namespace Opm::gpuistl

#endif // OPM_GPU_ECL_MATERIAL_LAW_MANAGER_HPP
