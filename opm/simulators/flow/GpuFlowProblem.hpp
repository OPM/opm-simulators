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
#include <opm/simulators/flow/GpuEclThermalLawManager.hpp>

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
 * \brief No-op thermal-law manager used as the default 4th template arg of
 *        \c GpuFlowProblem.
 *
 * Carries no data and exposes the minimum surface needed by the
 * \c GpuFlowProblem builder/copy/view machinery: a \c FluidSystem alias
 * fixed to \c void (so the \c hasThermal compile-time switch evaluates
 * to false), trivial empty params types, and trivial accessors that
 * return default-constructed values.
 *
 * Intentionally not callable on the device when actually queried: it is
 * only used by non-thermal \c GpuFlowProblem instantiations whose dispatch
 * paths never call \c solidEnergyLawParams or \c thermalConductionLawParams.
 */
struct NoThermalLawManager
{
    using FluidSystem = void;
    struct SolidEnergyLawParams {};
    struct ThermalConductionLawParams {};

    OPM_HOST_DEVICE SolidEnergyLawParams solidEnergyLawParams(unsigned /*elemIdx*/) const
    {
        return SolidEnergyLawParams{};
    }
    OPM_HOST_DEVICE ThermalConductionLawParams
    thermalConductionLawParams(unsigned /*elemIdx*/) const
    {
        return ThermalConductionLawParams{};
    }
};

} // namespace Opm

namespace Opm::gpuistl {

/*! \brief copy_to_gpu overload for the no-op thermal manager. */
inline ::Opm::NoThermalLawManager
copy_to_gpu(const ::Opm::NoThermalLawManager& /*cpu*/)
{
    return ::Opm::NoThermalLawManager{};
}

/*! \brief make_view overload for the no-op thermal manager. */
inline ::Opm::NoThermalLawManager
make_view(::Opm::NoThermalLawManager& /*buf*/)
{
    return ::Opm::NoThermalLawManager{};
}

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
 * \tparam ThermalLawManagerT  GPU thermal-law manager type (e.g.
 *                             \c Opm::EclThermalLaw::GpuManager). Defaults
 *                             to \c NoThermalLawManager, a no-op stub used
 *                             by non-thermal callers.
 */
template <class ScalarT,
          class MaterialLawManagerT,
          template <class> class Storage = ::Opm::VectorWithDefaultAllocator,
          class ThermalLawManagerT = ::Opm::NoThermalLawManager>
class GpuFlowProblem
{
public:
    using Scalar = ScalarT;
    using EclMaterialLawManager = MaterialLawManagerT;
    using MaterialLawParams = typename EclMaterialLawManager::MaterialLawParams;
    using EclThermalLawManager = ThermalLawManagerT;
    using SolidEnergyLawParams = typename EclThermalLawManager::SolidEnergyLawParams;
    using ThermalConductionLawParams = typename EclThermalLawManager::ThermalConductionLawParams;

    /*! \brief True iff the thermal-law manager carries a real fluid-system
     *         tag (anything other than \c void) and therefore the
     *         thermal-extraction code paths should be enabled. */
    static constexpr bool hasThermal
        = !std::is_same_v<typename EclThermalLawManager::FluidSystem, void>;

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
                                   Storage<Scalar> maxGasDissolutionFactor,
                                   EclThermalLawManager thermalLawManager,
                                   Storage<Scalar> rockFraction,
                                   Storage<int> pvtRegionIndex)
        : materialLawManager_(std::move(materialLawManager))
        , porosity_(std::move(porosity))
        , rockCompressibility_(std::move(rockCompressibility))
        , rockReferencePressure_(std::move(rockReferencePressure))
        , maxOilSaturation_(std::move(maxOilSaturation))
        , maxOilVaporizationFactor_(std::move(maxOilVaporizationFactor))
        , maxGasDissolutionFactor_(std::move(maxGasDissolutionFactor))
        , thermalLawManager_(std::move(thermalLawManager))
        , rockFraction_(std::move(rockFraction))
        , pvtRegionIndex_(std::move(pvtRegionIndex))
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
        if constexpr (hasThermal) {
            using FluidSystemTag = typename EclThermalLawManager::FluidSystem;
            thermalLawManager_
                = ::Opm::EclThermalLaw::buildCpuManagerFromFlowProblem<Scalar, FluidSystemTag>(
                    cpu, n);
            rockFraction_.resize(n);
            pvtRegionIndex_.resize(n);
            for (std::size_t i = 0; i < n; ++i) {
                const unsigned u = static_cast<unsigned>(i);
                rockFraction_[i] = cpu.rockFraction(u, 0u);
                pvtRegionIndex_[i] = static_cast<int>(cpu.pvtRegionIndex(u));
            }
        }
    }

    /*!
     * \brief Construct directly into device-resident \c GpuBuffer storage
     *        from any CPU \c FlowProblem (or \c FlowProblemBlackoil).
     *
     * Internally builds the per-cell rock/mixing data on the host as
     * \c std::vector<Scalar>, then uploads each as a single
     * \c GpuBuffer<Scalar>. The inner GPU material law manager is constructed
     * from \c *cpu.materialLawManager() via its own GpuBuffer-storage
     * constructor, which uploads the per-cell sample arrays.
     *
     * Only enabled when this problem itself uses \c GpuBuffer storage.
     */
    template <class CpuProblem,
              class StS = Storage<Scalar>,
              std::enable_if_t<
                  std::is_same_v<StS, ::Opm::gpuistl::GpuBuffer<Scalar>>,
                  int> = 0>
    explicit GpuFlowProblem(const CpuProblem& cpu)
        : materialLawManager_(*cpu.materialLawManager(), cpu.model().numGridDof())
        , porosity_(extractRockField(cpu, [](const CpuProblem& p, unsigned u) {
              return Scalar(p.porosity(u, 0u));
          }))
        , rockCompressibility_(extractRockField(cpu, [](const CpuProblem& p, unsigned u) {
              return Scalar(p.rockCompressibility(u));
          }))
        , rockReferencePressure_(extractRockField(cpu, [](const CpuProblem& p, unsigned u) {
              return Scalar(p.rockReferencePressure(u));
          }))
        , maxOilSaturation_(extractRockField(cpu, [](const CpuProblem& p, unsigned u) {
              return Scalar(p.maxOilSaturation(u));
          }))
        , maxOilVaporizationFactor_(extractRockField(cpu, [](const CpuProblem& p, unsigned u) {
              return Scalar(p.maxOilVaporizationFactor(0u, u));
          }))
        , maxGasDissolutionFactor_(extractRockField(cpu, [](const CpuProblem& p, unsigned u) {
              return Scalar(p.maxGasDissolutionFactor(0u, u));
          }))
        , thermalLawManager_(buildGpuThermalManager(cpu))
        , rockFraction_(extractThermalScalarField(cpu, [](const CpuProblem& p, unsigned u) {
              return Scalar(p.rockFraction(u, 0u));
          }))
        , pvtRegionIndex_(extractThermalIntField(cpu, [](const CpuProblem& p, unsigned u) {
              return static_cast<int>(p.pvtRegionIndex(u));
          }))
    {
    }

    OPM_HOST_DEVICE ModelView model() const
    {
        return ModelView{};
    }

    OPM_HOST_DEVICE int satnumRegionIndex(std::size_t elemIdx) const
    {
        return materialLawManager_.satnumRegionIdx(static_cast<unsigned>(elemIdx));
    }

    OPM_HOST_DEVICE MaterialLawParams materialLawParams(std::size_t elemIdx) const
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

    /*! \brief Per-cell PVT region index. Returns 0 when the GpuFlowProblem
     *         instantiation has no thermal support. */
    OPM_HOST_DEVICE unsigned pvtRegionIndex(std::size_t elemIdx) const
    {
        return pvtRegionIndex_.size() == 0
                   ? 0u
                   : static_cast<unsigned>(pvtRegionIndex_[elemIdx]);
    }

    /*! \brief Per-cell rock fraction (1 - effective porosity). Returns 0
     *         when the GpuFlowProblem instantiation has no thermal
     *         support. */
    OPM_HOST_DEVICE Scalar rockFraction(std::size_t elemIdx, unsigned /*timeIdx*/) const
    {
        return rockFraction_.size() == 0 ? Scalar(0) : rockFraction_[elemIdx];
    }

    /*! \brief Solid-energy law parameters for a single cell. Forwards to
     *         the embedded thermal-law manager. */
    OPM_HOST_DEVICE SolidEnergyLawParams
    solidEnergyLawParams(std::size_t elemIdx, unsigned /*timeIdx*/) const
    {
        return thermalLawManager_.solidEnergyLawParams(static_cast<unsigned>(elemIdx));
    }

    /*! \brief Thermal-conduction law parameters for a single cell.
     *         Forwards to the embedded thermal-law manager. */
    OPM_HOST_DEVICE ThermalConductionLawParams
    thermalConductionLawParams(std::size_t elemIdx, unsigned /*timeIdx*/) const
    {
        return thermalLawManager_.thermalConductionLawParams(static_cast<unsigned>(elemIdx));
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

    /*!
     * \brief Update the relative permeabilities of all phases for a single
     *        cell, in the same way as the CPU \c FlowProblem does.
     *
     * Calls \c MaterialLaw::relativePermeabilities on the cell's
     * \c MaterialLawParams. The \c BlackOilIntensiveQuantities later divides
     * the result by the phase viscosity to obtain the mobility.
     *
     * Directional relative permeabilities are not supported by the GPU
     * problem, so the \p dirMob output is left untouched.
     */
    template <class FluidState, class... Args>
    OPM_HOST_DEVICE void updateRelperms(auto& mobility,
                                        auto& /*dirMob*/,
                                        FluidState& fluidState,
                                        unsigned globalSpaceIdx) const
    {
        using ContainerT = std::decay_t<decltype(mobility)>;
        const auto materialParams = materialLawParams(globalSpaceIdx);
        EclMaterialLawManager::MaterialLaw::template relativePermeabilities<
            ContainerT, FluidState, Args...>(mobility, materialParams, fluidState);
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
    const EclThermalLawManager& thermalLawManager() const { return thermalLawManager_; }
    EclThermalLawManager& thermalLawManager() { return thermalLawManager_; }
    const Storage<Scalar>& rockFractionStorage() const { return rockFraction_; }
    Storage<Scalar>& rockFractionStorage() { return rockFraction_; }
    const Storage<int>& pvtRegionIndexStorage() const { return pvtRegionIndex_; }
    Storage<int>& pvtRegionIndexStorage() { return pvtRegionIndex_; }
    //!\}

private:
    template <class CpuProblem, class F>
    static Storage<Scalar> extractRockField(const CpuProblem& cpu, F f)
    {
        const std::size_t n = cpu.model().numGridDof();
        std::vector<Scalar> v(n);
        for (std::size_t i = 0; i < n; ++i) {
            v[i] = f(cpu, static_cast<unsigned>(i));
        }
        if constexpr (std::is_same_v<Storage<Scalar>,
                                     ::Opm::VectorWithDefaultAllocator<Scalar>>) {
            return Storage<Scalar>(v.begin(), v.end());
        } else {
            return Storage<Scalar>(v);
        }
    }

    /*! \brief Build the per-cell rockFraction storage; empty when the
     *         GpuFlowProblem instantiation has no thermal support. */
    template <class CpuProblem, class F>
    static Storage<Scalar> extractThermalScalarField(const CpuProblem& cpu, F f)
    {
        if constexpr (!hasThermal) {
            return Storage<Scalar>{};
        } else {
            return extractRockField(cpu, f);
        }
    }

    /*! \brief Build the per-cell pvtRegionIndex storage; empty when the
     *         GpuFlowProblem instantiation has no thermal support. */
    template <class CpuProblem, class F>
    static Storage<int> extractThermalIntField(const CpuProblem& cpu, F f)
    {
        if constexpr (!hasThermal) {
            return Storage<int>{};
        } else {
            const std::size_t n = cpu.model().numGridDof();
            std::vector<int> v(n);
            for (std::size_t i = 0; i < n; ++i) {
                v[i] = f(cpu, static_cast<unsigned>(i));
            }
            if constexpr (std::is_same_v<Storage<int>,
                                         ::Opm::VectorWithDefaultAllocator<int>>) {
                return Storage<int>(v.begin(), v.end());
            } else {
                return Storage<int>(v);
            }
        }
    }

    /*! \brief Build the GPU thermal-law manager from a CPU FlowProblem.
     *         Returns a default-constructed (empty) manager when this
     *         GpuFlowProblem instantiation has no thermal support. */
    template <class CpuProblem>
    static EclThermalLawManager buildGpuThermalManager(const CpuProblem& cpu)
    {
        if constexpr (!hasThermal) {
            return EclThermalLawManager{};
        } else {
            using FluidSystemTag = typename EclThermalLawManager::FluidSystem;
            auto cpuMgr
                = ::Opm::EclThermalLaw::buildCpuManagerFromFlowProblem<Scalar, FluidSystemTag>(
                    cpu, cpu.model().numGridDof());
            return ::Opm::gpuistl::copy_to_gpu(cpuMgr);
        }
    }

    EclMaterialLawManager materialLawManager_{};
    Storage<Scalar> porosity_{};
    Storage<Scalar> rockCompressibility_{};
    Storage<Scalar> rockReferencePressure_{};
    Storage<Scalar> maxOilSaturation_{};
    Storage<Scalar> maxOilVaporizationFactor_{};
    Storage<Scalar> maxGasDissolutionFactor_{};
    EclThermalLawManager thermalLawManager_{};
    Storage<Scalar> rockFraction_{};
    Storage<int> pvtRegionIndex_{};
};

} // namespace Opm

namespace Opm::gpuistl {

/*!
 * \brief Copy a CPU GpuFlowProblem to GPU-resident GpuBuffer storage.
 *
 *  The material law manager is also moved to the GPU (via copy_to_gpu).
 */
template <class ScalarT, class CpuMaterialLawManager, class CpuThermalLawManager>
auto copy_to_gpu(
    const ::Opm::GpuFlowProblem<ScalarT,
                                CpuMaterialLawManager,
                                ::Opm::VectorWithDefaultAllocator,
                                CpuThermalLawManager>& cpu)
{
    using GpuMaterialLawManagerBuffer
        = decltype(::Opm::gpuistl::copy_to_gpu(cpu.materialLawManager()));
    using GpuThermalLawManagerBuffer
        = decltype(::Opm::gpuistl::copy_to_gpu(cpu.thermalLawManager()));
    using GpuProblemType = ::Opm::GpuFlowProblem<ScalarT,
                                                 GpuMaterialLawManagerBuffer,
                                                 GpuBuffer,
                                                 GpuThermalLawManagerBuffer>;

    return GpuProblemType(::Opm::gpuistl::copy_to_gpu(cpu.materialLawManager()),
                          GpuBuffer<ScalarT>(cpu.porosityStorage()),
                          GpuBuffer<ScalarT>(cpu.rockCompressibilityStorage()),
                          GpuBuffer<ScalarT>(cpu.rockReferencePressureStorage()),
                          GpuBuffer<ScalarT>(cpu.maxOilSaturationStorage()),
                          GpuBuffer<ScalarT>(cpu.maxOilVaporizationFactorStorage()),
                          GpuBuffer<ScalarT>(cpu.maxGasDissolutionFactorStorage()),
                          ::Opm::gpuistl::copy_to_gpu(cpu.thermalLawManager()),
                          GpuBuffer<ScalarT>(cpu.rockFractionStorage()),
                          GpuBuffer<int>(cpu.pvtRegionIndexStorage()));
}

/*!
 * \brief Make a non-owning GpuView based GpuFlowProblem from an owning
 *        GpuBuffer based GpuFlowProblem.
 */
template <class ScalarT, class GpuBufferMaterialLawManager, class GpuBufferThermalLawManager>
auto make_view(::Opm::GpuFlowProblem<ScalarT,
                                     GpuBufferMaterialLawManager,
                                     GpuBuffer,
                                     GpuBufferThermalLawManager>& buf)
{
    using GpuMaterialLawManagerView
        = decltype(::Opm::gpuistl::make_view(buf.materialLawManager()));
    using GpuThermalLawManagerView
        = decltype(::Opm::gpuistl::make_view(buf.thermalLawManager()));
    using GpuProblemView = ::Opm::GpuFlowProblem<ScalarT,
                                                 GpuMaterialLawManagerView,
                                                 GpuView,
                                                 GpuThermalLawManagerView>;

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
                          toView(buf.maxGasDissolutionFactorStorage()),
                          ::Opm::gpuistl::make_view(buf.thermalLawManager()),
                          toView(buf.rockFractionStorage()),
                          toView(buf.pvtRegionIndexStorage()));
}

} // namespace Opm::gpuistl

#endif // OPM_GPU_FLOW_PROBLEM_HPP
