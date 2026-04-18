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
 * \brief GPU-compatible, simplified version of Opm::FlowProblem.
 *
 * Only exposes the subset of the FlowProblem / FlowProblemBlackoil interface
 * that BlackOilIntensiveQuantities needs in order to call ::update(). Anything
 * that is not used by the intensive quantities (well model, aquifers, source
 * terms, boundary conditions, ...) is intentionally absent.
 */
#ifndef OPM_GPU_FLOW_PROBLEM_HPP
#define OPM_GPU_FLOW_PROBLEM_HPP

#include <opm/common/utility/VectorWithDefaultAllocator.hpp>
#include <opm/common/utility/gpuDecorators.hpp>

#include <opm/models/discretization/common/linearizationtype.hh>

#include <opm/simulators/flow/GpuEclMaterialLawManager.hpp>

#include <cstddef>
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

namespace Opm {

/*!
 * \brief Minimal, GPU-compatible problem class.
 *
 * Provides exactly the interface required by BlackOilIntensiveQuantities to
 * compute a primary-variable update on the device. All cell-wise rock /
 * fluid bookkeeping data is kept in templated Storage so the same class can
 * be used as a CPU object (default \c std::vector) or as a GPU object
 * (\c GpuBuffer / \c GpuView).
 *
 * The material law manager is held by value (so it is trivially copyable
 * into a CUDA kernel together with the rest of the problem state).
 *
 * \tparam ScalarT             Floating point type used for rock data.
 * \tparam MaterialLawManagerT GPU material law manager type (e.g.
 *                             \c Opm::EclMaterialLaw::GpuManager).
 * \tparam Storage             Container template used for per-cell data.
 */
template <class ScalarT,
          class MaterialLawManagerT,
          template <class> class Storage = ::Opm::VectorWithDefaultAllocator>
class GpuFlowProblem
{
public:
    using Scalar = ScalarT;
    using EclMaterialLawManager = MaterialLawManagerT;
    using MaterialLawParams = typename EclMaterialLawManager::MaterialLawParams;

    /*! \brief Trivial nested model() helper that satisfies the
     *         \c problem.model().linearizer().getLinearizationType() chain
     *         that BlackOilIntensiveQuantities walks. */
    struct ModelView
    {
        struct LinearizerView
        {
            OPM_HOST_DEVICE LinearizationType getLinearizationType() const
            {
                return LinearizationType{};
            }
        };

        OPM_HOST_DEVICE LinearizerView linearizer() const
        {
            return LinearizerView{};
        }
    };

    OPM_HOST_DEVICE GpuFlowProblem() = default;

    OPM_HOST_DEVICE GpuFlowProblem(EclMaterialLawManager materialLawManager,
                                   Storage<Scalar> porosity,
                                   Storage<Scalar> rockCompressibility,
                                   Storage<Scalar> rockReferencePressure,
                                   Storage<Scalar> maxOilSaturation,
                                   Storage<Scalar> maxOilVaporizationFactor,
                                   Storage<Scalar> maxGasDissolutionFactor)
        : materialLawManager_(std::move(materialLawManager))
        , porosity_(std::move(porosity))
        , rockCompressibility_(std::move(rockCompressibility))
        , rockReferencePressure_(std::move(rockReferencePressure))
        , maxOilSaturation_(std::move(maxOilSaturation))
        , maxOilVaporizationFactor_(std::move(maxOilVaporizationFactor))
        , maxGasDissolutionFactor_(std::move(maxGasDissolutionFactor))
    {
    }

    /*!
     * \brief Construct from any CPU \c FlowProblem (or \c FlowProblemBlackoil).
     *
     * Iterates over all grid degrees of freedom and copies the rock /
     * mixing-control fields used by \c BlackOilIntensiveQuantities into local
     * vectors. The inner material-law manager is constructed from
     * \c *cpu.materialLawManager(), so the GPU manager type held by this
     * problem must in turn be constructible from
     * \c Opm::EclMaterialLaw::Manager (i.e. use CPU vector storage and
     * piecewise-linear vector-storage two-phase laws).
     *
     * Only enabled when this problem itself uses CPU vector storage.
     */
    template <class CpuProblem,
              class StS = Storage<Scalar>,
              std::enable_if_t<
                  std::is_same_v<StS, ::Opm::VectorWithDefaultAllocator<Scalar>>,
                  int> = 0>
    explicit GpuFlowProblem(const CpuProblem& cpu)
        : materialLawManager_(*cpu.materialLawManager(), cpu.model().numGridDof())
    {
        const std::size_t n = cpu.model().numGridDof();
        porosity_.resize(n);
        rockCompressibility_.resize(n);
        rockReferencePressure_.resize(n);
        maxOilSaturation_.resize(n);
        maxOilVaporizationFactor_.resize(n);
        maxGasDissolutionFactor_.resize(n);
        for (std::size_t i = 0; i < n; ++i) {
            const unsigned u = static_cast<unsigned>(i);
            porosity_[i]                = cpu.porosity(u, 0u);
            rockCompressibility_[i]     = cpu.rockCompressibility(u);
            rockReferencePressure_[i]   = cpu.rockReferencePressure(u);
            maxOilSaturation_[i]        = cpu.maxOilSaturation(u);
            maxOilVaporizationFactor_[i] = cpu.maxOilVaporizationFactor(0u, u);
            maxGasDissolutionFactor_[i]  = cpu.maxGasDissolutionFactor(0u, u);
        }
    }

    OPM_HOST_DEVICE ModelView model() const
    {
        return ModelView{};
    }

    OPM_HOST_DEVICE int satnumRegionIndex(std::size_t elemIdx) const
    {
        return materialLawManager_.satnumRegionIdx(static_cast<unsigned>(elemIdx));
    }

    OPM_HOST_DEVICE const MaterialLawParams& materialLawParams(std::size_t elemIdx) const
    {
        return materialLawManager_.materialLawParams(static_cast<unsigned>(elemIdx));
    }

    OPM_HOST_DEVICE Scalar rockCompressibility(std::size_t elemIdx) const
    {
        return rockCompressibility_.size() == 0 ? Scalar(0) : rockCompressibility_[elemIdx];
    }

    OPM_HOST_DEVICE Scalar rockReferencePressure(std::size_t elemIdx) const
    {
        return rockReferencePressure_.size() == 0 ? Scalar(0) : rockReferencePressure_[elemIdx];
    }

    OPM_HOST_DEVICE Scalar porosity(std::size_t elemIdx, unsigned /*timeIdx*/) const
    {
        return porosity_.size() == 0 ? Scalar(0) : porosity_[elemIdx];
    }

    OPM_HOST_DEVICE Scalar maxOilVaporizationFactor(unsigned /*timeIdx*/, std::size_t elemIdx) const
    {
        return maxOilVaporizationFactor_.size() == 0 ? Scalar(0) : maxOilVaporizationFactor_[elemIdx];
    }

    OPM_HOST_DEVICE Scalar maxGasDissolutionFactor(unsigned /*timeIdx*/, std::size_t elemIdx) const
    {
        return maxGasDissolutionFactor_.size() == 0 ? Scalar(0) : maxGasDissolutionFactor_[elemIdx];
    }

    OPM_HOST_DEVICE Scalar maxOilSaturation(std::size_t elemIdx) const
    {
        return maxOilSaturation_.size() == 0 ? Scalar(0) : maxOilSaturation_[elemIdx];
    }

    /*! \brief Default rock-pore-volume multiplier (1 if no compressibility). */
    template <class Evaluation>
    OPM_HOST_DEVICE Evaluation rockCompPoroMultiplier(const auto& /*intQuants*/, std::size_t /*elemIdx*/) const
    {
        return Evaluation(1.0);
    }

    /*! \brief Default rock-trans multiplier (1 if no compressibility). */
    template <class Evaluation>
    OPM_HOST_DEVICE Evaluation rockCompTransMultiplier(const auto& /*intQuants*/, std::size_t /*elemIdx*/) const
    {
        return Evaluation(1.0);
    }

    /*! \brief No-op relperm update; intensive quantities are expected to be
     *         updated through the material law manager directly when needed. */
    template <class FluidState, class... Args>
    OPM_HOST_DEVICE void updateRelperms(auto& /*mobility*/,
                                        auto& /*dirMob*/,
                                        FluidState& /*fluidState*/,
                                        unsigned /*globalSpaceIdx*/) const
    {
    }

    /*! \name Direct accessors for serialization to the GPU. */
    //!\{
    const EclMaterialLawManager& materialLawManager() const { return materialLawManager_; }
    EclMaterialLawManager& materialLawManager() { return materialLawManager_; }
    const Storage<Scalar>& porosityStorage() const { return porosity_; }
    Storage<Scalar>& porosityStorage() { return porosity_; }
    const Storage<Scalar>& rockCompressibilityStorage() const { return rockCompressibility_; }
    Storage<Scalar>& rockCompressibilityStorage() { return rockCompressibility_; }
    const Storage<Scalar>& rockReferencePressureStorage() const { return rockReferencePressure_; }
    Storage<Scalar>& rockReferencePressureStorage() { return rockReferencePressure_; }
    const Storage<Scalar>& maxOilSaturationStorage() const { return maxOilSaturation_; }
    Storage<Scalar>& maxOilSaturationStorage() { return maxOilSaturation_; }
    const Storage<Scalar>& maxOilVaporizationFactorStorage() const { return maxOilVaporizationFactor_; }
    Storage<Scalar>& maxOilVaporizationFactorStorage() { return maxOilVaporizationFactor_; }
    const Storage<Scalar>& maxGasDissolutionFactorStorage() const { return maxGasDissolutionFactor_; }
    Storage<Scalar>& maxGasDissolutionFactorStorage() { return maxGasDissolutionFactor_; }
    //!\}

private:
    EclMaterialLawManager materialLawManager_{};
    Storage<Scalar> porosity_{};
    Storage<Scalar> rockCompressibility_{};
    Storage<Scalar> rockReferencePressure_{};
    Storage<Scalar> maxOilSaturation_{};
    Storage<Scalar> maxOilVaporizationFactor_{};
    Storage<Scalar> maxGasDissolutionFactor_{};
};

} // namespace Opm

namespace Opm::gpuistl {

/*!
 * \brief Copy a CPU GpuFlowProblem to GPU-resident GpuBuffer storage.
 *
 *  The material law manager is also moved to the GPU (via copy_to_gpu).
 */
template <class ScalarT, class CpuMaterialLawManager>
auto copy_to_gpu(const ::Opm::GpuFlowProblem<ScalarT, CpuMaterialLawManager>& cpu)
{
    using GpuMaterialLawManagerBuffer = decltype(::Opm::gpuistl::copy_to_gpu(cpu.materialLawManager()));
    using GpuProblemType = ::Opm::GpuFlowProblem<ScalarT, GpuMaterialLawManagerBuffer, GpuBuffer>;

    return GpuProblemType(::Opm::gpuistl::copy_to_gpu(cpu.materialLawManager()),
                          GpuBuffer<ScalarT>(cpu.porosityStorage()),
                          GpuBuffer<ScalarT>(cpu.rockCompressibilityStorage()),
                          GpuBuffer<ScalarT>(cpu.rockReferencePressureStorage()),
                          GpuBuffer<ScalarT>(cpu.maxOilSaturationStorage()),
                          GpuBuffer<ScalarT>(cpu.maxOilVaporizationFactorStorage()),
                          GpuBuffer<ScalarT>(cpu.maxGasDissolutionFactorStorage()));
}

/*!
 * \brief Make a non-owning GpuView based GpuFlowProblem from an owning
 *        GpuBuffer based GpuFlowProblem.
 */
template <class ScalarT, class GpuBufferMaterialLawManager>
auto make_view(::Opm::GpuFlowProblem<ScalarT, GpuBufferMaterialLawManager, GpuBuffer>& buf)
{
    using GpuMaterialLawManagerView = decltype(::Opm::gpuistl::make_view(buf.materialLawManager()));
    using GpuProblemView = ::Opm::GpuFlowProblem<ScalarT, GpuMaterialLawManagerView, GpuView>;

    auto toView = [](auto& storage) {
        using T = typename std::decay_t<decltype(storage)>::value_type;
        return GpuView<T>(storage.data(), storage.size());
    };

    return GpuProblemView(::Opm::gpuistl::make_view(buf.materialLawManager()),
                          toView(buf.porosityStorage()),
                          toView(buf.rockCompressibilityStorage()),
                          toView(buf.rockReferencePressureStorage()),
                          toView(buf.maxOilSaturationStorage()),
                          toView(buf.maxOilVaporizationFactorStorage()),
                          toView(buf.maxGasDissolutionFactorStorage()));
}

} // namespace Opm::gpuistl

#endif // OPM_GPU_FLOW_PROBLEM_HPP
