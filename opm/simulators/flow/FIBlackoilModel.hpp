// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*
  This file is part of the Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 2 of the License, or
  (at your option) any later version.

  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.

  Consult the COPYING file in the top-level source directory of this
  module for the precise wording of the license and the list of
  copyright holders.
*/
/*!
 * \file
 *
 * \copydoc Opm::FIBlackOilModel
 */
#ifndef FI_BLACK_OIL_MODEL_HPP
#define FI_BLACK_OIL_MODEL_HPP

#include <opm/models/blackoil/blackoilmodel.hh>
#include <opm/models/utils/propertysystem.hh>

#include <opm/common/ErrorMacros.hpp>
#include <opm/models/utils/propertysystem.hh>

#include <opm/grid/CpGrid.hpp>
#include <opm/grid/utility/ElementChunks.hpp>

#include <opm/models/parallel/threadmanager.hpp>
#include <opm/simulators/utils/DeferredLoggingErrorHelpers.hpp>

#include <opm/material/fluidmatrixinteractions/EclMultiplexerMaterialParams.hpp>

#include <opm/common/OpmLog/OpmLog.hpp>

#include <opm/simulators/flow/FlowProblemParameters.hpp>

#if HAVE_CUDA
#include <opm/simulators/linalg/gpuistl/GpuBlackoilIntensiveQuantitiesDispatcher.hpp>
#include <memory>
#include <variant>
#include <vector>
#endif

#include <chrono>
#include <cstddef>
#include <format>
#include <stdexcept>
#include <type_traits>

namespace Opm {

template <typename TypeTag>
class FIBlackOilModel : public BlackOilModel<TypeTag>
{
    using ParentType = BlackOilModel<TypeTag>;
    using Simulator = GetPropType<TypeTag, Properties::Simulator>;
    using IntensiveQuantities = GetPropType<TypeTag, Properties::IntensiveQuantities>;
    using ElementContext = GetPropType<TypeTag, Properties::ElementContext>;
    using ThreadManager = GetPropType<TypeTag, Properties::ThreadManager>;
    using GridView = GetPropType<TypeTag, Properties::GridView>;
    using Element = typename GridView::template Codim<0>::Entity;
    using ElementIterator = typename GridView::template Codim<0>::Iterator;
    enum {
        numEq = getPropValue<TypeTag, Properties::NumEq>(),
        historySize = getPropValue<TypeTag, Properties::TimeDiscHistorySize>(),
    };

    // The chunked and threaded iteration over elements in this class assumes that the number
    // and order of elements is fixed, and is therefore constrained to only work with CpGrid.
    // For example, ALUGrid supports refinement and can not assume that.
    static constexpr bool gridIsUnchanging =
        std::is_same_v<GetPropType<TypeTag, Properties::Grid>, Dune::CpGrid>;

    static constexpr bool avoidElementContext = getPropValue<TypeTag, Properties::AvoidElementContext>();

public:
    explicit FIBlackOilModel(Simulator& simulator)
        : BlackOilModel<TypeTag>(simulator)
        , element_chunks_(this->gridView_,
                          Dune::Partitions::all,
                          ThreadManager::maxThreads())
    {
        initGpuIntensiveQuantitiesDispatcher_();
    }

    void invalidateAndUpdateIntensiveQuantities(unsigned timeIdx) const
    {
        this->invalidateIntensiveQuantitiesCache(timeIdx);

        // When the experimental GPU dispatcher is active and the CPU has
        // already populated the cache once (including the fields the GPU
        // dispatcher does not overlay, e.g. mobility / energy), skip the
        // expensive CPU loop and rely solely on the GPU dispatcher to
        // refresh the BlackOil intensive-quantities fields. The remaining
        // cached fields keep their previously computed values, which is
        // the intended behaviour for this experimental GPU-only path.
        if (gpuDispatcherActiveAndInitialized_()) {
            markIntensiveQuantitiesCacheValid_(timeIdx);
            maybeRunGpuIntensiveQuantitiesDispatcher_(timeIdx);
            return;
        }

        if constexpr (gridIsUnchanging) {
            if constexpr (avoidElementContext) {
                updateCachedIntQuants(timeIdx);
                cpuIntensiveQuantitiesInitialized_ = true;
                return;
            }
            OPM_BEGIN_PARALLEL_TRY_CATCH();
#ifdef _OPENMP
#pragma omp parallel for
#endif
            for (const auto& chunk : element_chunks_) {
                ElementContext elemCtx(this->simulator_);
                for (const auto& elem : chunk) {
                    elemCtx.updatePrimaryStencil(elem);
                    elemCtx.updatePrimaryIntensiveQuantities(timeIdx);
                }
            }
            OPM_END_PARALLEL_TRY_CATCH("invalidateAndUpdateIntensiveQuantities: state error",
                                       this->simulator_.vanguard().grid().comm());
            // After all cells are CPU-updated and written into the cache,
            // overlay the GPU-computed BlackOil fields in one batched call
            // (no-op when the GPU dispatcher is unavailable or disabled).
            maybeRunGpuIntensiveQuantitiesDispatcher_(timeIdx);
            cpuIntensiveQuantitiesInitialized_ = true;
        } else {
            // Grid is possibly refined or otherwise changed between calls.
            ElementContext elemCtx(this->simulator_);
            for (const auto& elem : elements(this->gridView_)) {
                elemCtx.updatePrimaryStencil(elem);
                elemCtx.updatePrimaryIntensiveQuantities(timeIdx);
            }
        }
    }

    void invalidateAndUpdateIntensiveQuantitiesOverlap(unsigned timeIdx) const
    {
        // loop over all elements
        ThreadedEntityIterator<GridView, /*codim=*/0> threadedElemIt(this->gridView_);
        OPM_BEGIN_PARALLEL_TRY_CATCH()
#ifdef _OPENMP
#pragma omp parallel
#endif
        {
            ElementContext elemCtx(this->simulator_);
            auto elemIt = threadedElemIt.beginParallel();
            for (; !threadedElemIt.isFinished(elemIt); elemIt = threadedElemIt.increment()) {
                if (elemIt->partitionType() != Dune::OverlapEntity) {
                    continue;
                }
                const Element& elem = *elemIt;
                elemCtx.updatePrimaryStencil(elem);
                // Mark cache for this element as invalid.
                const std::size_t numPrimaryDof = elemCtx.numPrimaryDof(timeIdx);
                for (unsigned dofIdx = 0; dofIdx < numPrimaryDof; ++dofIdx) {
                    const unsigned globalIndex = elemCtx.globalSpaceIndex(dofIdx, timeIdx);
                    this->setIntensiveQuantitiesCacheEntryValidity(globalIndex, timeIdx, false);
                }
                // Update for this element.
                elemCtx.updatePrimaryIntensiveQuantities(/*timeIdx=*/0);
            }
        }
        OPM_END_PARALLEL_TRY_CATCH("InvalideAndUpdateIntensiveQuantitiesOverlap: state error",
                                   this->simulator_.vanguard().grid().comm());
    }

    template <class GridSubDomain>
    void invalidateAndUpdateIntensiveQuantities(unsigned timeIdx, const GridSubDomain& gridSubDomain) const
    {
        // loop over all elements in the subdomain
        using GridViewType = decltype(gridSubDomain.view);
        ThreadedEntityIterator<GridViewType, /*codim=*/0> threadedElemIt(gridSubDomain.view);
#ifdef _OPENMP
#pragma omp parallel
#endif
        {
            ElementContext elemCtx(this->simulator_);
            auto elemIt = threadedElemIt.beginParallel();
            for (; !threadedElemIt.isFinished(elemIt); elemIt = threadedElemIt.increment()) {
                if (elemIt->partitionType() != Dune::InteriorEntity) {
                    continue;
                }
                const Element& elem = *elemIt;
                elemCtx.updatePrimaryStencil(elem);
                // Mark cache for this element as invalid.
                const std::size_t numPrimaryDof = elemCtx.numPrimaryDof(timeIdx);
                for (unsigned dofIdx = 0; dofIdx < numPrimaryDof; ++dofIdx) {
                    const unsigned globalIndex = elemCtx.globalSpaceIndex(dofIdx, timeIdx);
                    this->setIntensiveQuantitiesCacheEntryValidity(globalIndex, timeIdx, false);
                }
                // Update for this element.
                elemCtx.updatePrimaryIntensiveQuantities(/*timeIdx=*/0);
            }
        }
    }

    /*!
     * \brief Called by the update() method if it was
     *        unsuccessful. This is primary a hook which the actual
     *        model can overload.
     */
    void updateFailed()
    {
        // Reset the current solution to the one of the
        // previous time step so that we can start the next
        // update at a physically meaningful solution.
        // this->solution(/*timeIdx=*/0) = this->solution(/*timeIdx=*/1);
        ParentType::updateFailed();
        invalidateAndUpdateIntensiveQuantities(/*timeIdx=*/0);
    }

    // standard flow
    const IntensiveQuantities& intensiveQuantities(unsigned globalIdx, unsigned timeIdx) const
    {
        if (!this->enableIntensiveQuantityCache_) {
            OPM_THROW(std::logic_error,
                      "Run without intensive quantites not enabled: "
                      "Use --enable-intensive-quantity=true");
        }

        assert(timeIdx < this->cachedIntensiveQuantityHistorySize());
        const auto* intquant = this->cachedIntensiveQuantities(globalIdx, timeIdx);
        if (!intquant) {
            OPM_THROW(std::logic_error, "Intensive quantites need to be updated in code");
        }
        return *intquant;
    }

protected:

    template <EclMultiplexerApproach ApproachArg>
    using EMD = EclMultiplexerDispatch<ApproachArg>;

    void updateCachedIntQuants(const unsigned timeIdx) const
    {
        // Runtime dispatch on three-phase approach to compile-time fixed approach.
        switch (this->simulator_.problem().materialLawManager()->threePhaseApproach()) {
        case EclMultiplexerApproach::Stone1:
            updateCachedIntQuants1<EclMultiplexerDispatch<EclMultiplexerApproach::Stone1>>(timeIdx);
            break;

        case EclMultiplexerApproach::Stone2:
            updateCachedIntQuants1<EMD<EclMultiplexerApproach::Stone2>>(timeIdx);
            break;

        case EclMultiplexerApproach::Default:
            updateCachedIntQuants1<EMD<EclMultiplexerApproach::Default>>(timeIdx);
            break;

        case EclMultiplexerApproach::TwoPhase:
            updateCachedIntQuants1<EMD<EclMultiplexerApproach::TwoPhase>>(timeIdx);
            break;

        case EclMultiplexerApproach::OnePhase:
            updateCachedIntQuants1<EMD<EclMultiplexerApproach::OnePhase>>(timeIdx);
            break;
        }
    }

    template <class EMDArg>
    void updateCachedIntQuants1(const unsigned timeIdx) const
    {
        // Runtime dispatch on using piecewise linear saturation curves to compile-time fixed approach.
        if (this->simulator_.problem().materialLawManager()->satCurveIsAllPiecewiseLinear()) {
            using PL = SatCurveMultiplexerDispatch<SatCurveMultiplexerApproach::PiecewiseLinear>;
            updateCachedIntQuantsLoop<EMDArg, PL>(timeIdx);
        } else {
            // TODO: Might want to set LET here, but need to check if partial use of LET is possible.
            updateCachedIntQuantsLoop<EMDArg>(timeIdx);
        }
    }


    template <class ...Args>
    void updateCachedIntQuantsLoop(const unsigned timeIdx) const
    {
        const auto& elementMapper = this->simulator_.model().elementMapper();
#ifdef _OPENMP
#pragma omp parallel for
#endif
        for (const auto& chunk : element_chunks_) {
            for (const auto& elem : chunk) {
                this->template updateSingleCachedIntQuantUnchecked<Args...>(elementMapper.index(elem), timeIdx);
            }
        }

        // After the CPU per-cell update has populated all fields, optionally
        // overlay the BlackOil intensive-quantities fields with their GPU
        // counterparts via the experimental dispatcher. The dispatcher
        // overwrites the subset of fields covered by
        // BlackOilIntensiveQuantities::overlayBlackOilFieldsFrom; any field
        // not handled there keeps the CPU-computed value.
        // The call is a no-op when the GPU dispatcher is unavailable or the
        // user has not enabled it.
        maybeRunGpuIntensiveQuantitiesDispatcher_(timeIdx);
    }

    template <class ...Args>
    void updateSingleCachedIntQuantUnchecked(const unsigned globalIdx, const unsigned timeIdx) const
    {
        // Get the cached data.
        auto& intquant = this->intensiveQuantityCache_[timeIdx][globalIdx];
        // Update it.
        intquant.template update<Args...>(this->simulator_.problem(), this->solution(timeIdx)[globalIdx], globalIdx, timeIdx);
        // Set the up-to-date flag.
        this->intensiveQuantityCacheUpToDate_[timeIdx][globalIdx] = 1;
    }

    ElementChunks<GridView, Dune::Partitions::All> element_chunks_;

    // ----------------------------------------------------------------------
    // GPU intensive-quantities dispatcher integration.
    //
    // All knowledge about whether the experimental GPU dispatcher is
    // available, whether it is supported for the current TypeTag, and how
    // to invoke it, is intentionally isolated to the few private members
    // and helper methods below. The rest of the model only sees a single
    // entry point: maybeRunGpuIntensiveQuantitiesDispatcher_(timeIdx).
    // ----------------------------------------------------------------------

    // Compile-time predicates. We use them to avoid sprinkling preprocessor
    // conditionals throughout the rest of the class.
    static constexpr bool gpuDispatcherCompiledIn_ =
#if HAVE_CUDA && OPM_HAVE_GPU_BLACKOIL_INTENSIVE_QUANTITIES_DISPATCHER
        true;
#else
        false;
#endif

    static constexpr bool gpuDispatcherSupportsTypeTag_ =
#if HAVE_CUDA && OPM_HAVE_GPU_BLACKOIL_INTENSIVE_QUANTITIES_DISPATCHER
        Opm::gpuistl::GpuBlackoilIntensiveQuantitiesDispatcherSupport<TypeTag>::value;
#else
        false;
#endif

    bool useGpuIntensiveQuantitiesDispatcher_{false};

    // Tracks whether the CPU intensive-quantities update has been run at
    // least once. Used to ensure non-BlackOil cached fields (mobility,
    // energy, ...) are populated before we switch to the GPU-only path
    // when the experimental dispatcher is enabled.
    mutable bool cpuIntensiveQuantitiesInitialized_{false};

    // True iff the GPU dispatcher is compiled in, supported for the
    // current TypeTag, enabled at runtime, and the CPU has already been
    // run once to populate the cached non-BlackOil fields.
    bool gpuDispatcherActiveAndInitialized_() const noexcept
    {
        if constexpr (gpuDispatcherCompiledIn_ && gpuDispatcherSupportsTypeTag_) {
            return useGpuIntensiveQuantitiesDispatcher_
                && cpuIntensiveQuantitiesInitialized_;
        } else {
            return false;
        }
    }

    // Mark every entry of the intensive-quantity cache for the given time
    // index as valid. Used on the GPU-only path where we skip the CPU
    // update loop (which would normally mark each entry valid as a
    // side-effect) and rely on the GPU dispatcher to refresh the
    // BlackOil fields in-place.
    void markIntensiveQuantitiesCacheValid_(const unsigned timeIdx) const
    {
        const std::size_t numCells = this->intensiveQuantityCache_[timeIdx].size();
        for (std::size_t i = 0; i < numCells; ++i) {
            this->setIntensiveQuantitiesCacheEntryValidity(i, timeIdx, true);
        }
    }

#if HAVE_CUDA && OPM_HAVE_GPU_BLACKOIL_INTENSIVE_QUANTITIES_DISPATCHER
    using GpuDispatcherStorage = std::conditional_t<
        Opm::gpuistl::GpuBlackoilIntensiveQuantitiesDispatcherSupport<TypeTag>::value,
        std::unique_ptr<Opm::gpuistl::GpuBlackoilIntensiveQuantitiesDispatcher<TypeTag>>,
        std::monostate>;
    mutable GpuDispatcherStorage gpuIntensiveQuantitiesDispatcher_{};
#endif

    // Read the user-facing parameter and validate that the requested
    // configuration is actually supported by this build / TypeTag.
    // Throws (with a self-contained reason) when the user requests the GPU
    // dispatcher but it is not available.
    void initGpuIntensiveQuantitiesDispatcher_()
    {
        const bool requested =
            Parameters::Get<Parameters::ExperimentalComputePropertiesOnGpu>();
        if (!requested) {
            useGpuIntensiveQuantitiesDispatcher_ = false;
            return;
        }

        if constexpr (!gpuDispatcherCompiledIn_) {
            OPM_THROW(std::runtime_error,
                      "--experimental-compute-properties-on-gpu=true was "
                      "specified, but this binary was built without a "
                      "compatible GPU back-end. The GPU intensive-quantities "
                      "dispatcher requires either HIP or CUDA >= 13.1; "
                      "rebuild with HIP or with a sufficiently new CUDA "
                      "toolkit (see CMake option "
                      "OPM_HAVE_GPU_BLACKOIL_INTENSIVE_QUANTITIES_DISPATCHER).");
        } else if constexpr (!gpuDispatcherSupportsTypeTag_) {
            OPM_THROW(std::runtime_error,
                      "--experimental-compute-properties-on-gpu=true was "
                      "specified, but the active TypeTag is not supported "
                      "by the GPU BlackOil intensive-quantities dispatcher. "
                      "Only the CO2STORE-compatible TypeTag "
                      "FlowGasWaterEnergyProblem (gas+water+energy, no oil) "
                      "is currently supported by the experimental GPU "
                      "dispatcher.");
        } else {
            useGpuIntensiveQuantitiesDispatcher_ = true;
        }
    }

    // Single entry point that hides all dispatcher-related logic from the
    // rest of the model. No-op when the dispatcher is either not compiled
    // in, not supported for the current TypeTag, or not enabled at runtime.
    void maybeRunGpuIntensiveQuantitiesDispatcher_(const unsigned timeIdx) const
    {
        if constexpr (gpuDispatcherCompiledIn_ && gpuDispatcherSupportsTypeTag_) {
            if (useGpuIntensiveQuantitiesDispatcher_) {
                const auto gpuStartTime = std::chrono::steady_clock::now();
                this->runGpuIntensiveQuantitiesDispatcher_(timeIdx);
                const auto gpuDuration =
                    std::chrono::duration_cast<std::chrono::milliseconds>(
                        std::chrono::steady_clock::now() - gpuStartTime);
                Opm::OpmLog::info(std::format(
                    "GPU intensive-quantities dispatch took {} ms",
                    gpuDuration.count()));
            }
        } else {
            (void)timeIdx;
        }
    }

#if HAVE_CUDA && OPM_HAVE_GPU_BLACKOIL_INTENSIVE_QUANTITIES_DISPATCHER
    void runGpuIntensiveQuantitiesDispatcher_(const unsigned timeIdx) const
    {
        if constexpr (Opm::gpuistl::GpuBlackoilIntensiveQuantitiesDispatcherSupport<TypeTag>::value) {
            if (!gpuIntensiveQuantitiesDispatcher_) {
                gpuIntensiveQuantitiesDispatcher_ =
                    std::make_unique<
                        Opm::gpuistl::GpuBlackoilIntensiveQuantitiesDispatcher<TypeTag>>();
            }
            using PV = GetPropType<TypeTag, Properties::PrimaryVariables>;
            const auto& sol = this->solution(timeIdx);
            const std::size_t numCells = this->intensiveQuantityCache_[timeIdx].size();
            std::vector<const PV*> priVarsPtrs(numCells);
            std::vector<IntensiveQuantities*> outIqPtrs(numCells);
            for (std::size_t i = 0; i < numCells; ++i) {
                priVarsPtrs[i] = &sol[i];
                outIqPtrs[i] = &this->intensiveQuantityCache_[timeIdx][i];
            }
            gpuIntensiveQuantitiesDispatcher_->update(
                this->simulator_.problem(),
                priVarsPtrs.data(),
                outIqPtrs.data(),
                numCells);
        }
    }
#endif
};

} // namespace Opm

#endif // FI_BLACK_OIL_MODEL_HPP
