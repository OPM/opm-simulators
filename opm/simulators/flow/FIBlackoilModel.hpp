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

#if HAVE_CUDA
#include <opm/simulators/flow/FlowProblemParameters.hpp>
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
#include <mutex>
#include <thread>

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
#if HAVE_CUDA
        useGpuIntensiveQuantitiesDispatcher_ =
            Parameters::Get<Parameters::ExperimentalComputePropertiesOnGpu>();
#endif
    }

    void invalidateAndUpdateIntensiveQuantities(unsigned timeIdx) const
    {

        this->invalidateIntensiveQuantitiesCache(timeIdx);
        if constexpr (gridIsUnchanging) {
            if constexpr (avoidElementContext) {
                updateCachedIntQuants(timeIdx);
                return;
            }

            if (updateHasRun < 4) {
                const auto timeBegin = std::chrono::steady_clock::now();
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
                const auto timeEnd = std::chrono::steady_clock::now();
                const auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(timeEnd - timeBegin).count();
                OpmLog::info(std::format("Updated intensive quantities for {} elements in {} ms", this->gridView_.size(0), duration));
                ++updateHasRun;
            } else {
                OpmLog::info("Intensive quantities update skipped since it has already been run once");
            }

#if HAVE_CUDA
            // After all cells are CPU-updated and written into the cache,
            // overlay the GPU-computed BlackOil fields in one batched call.
            if constexpr (Opm::gpuistl::GpuBlackoilIntensiveQuantitiesDispatcherSupport<TypeTag>::value)
            {
                if (useGpuIntensiveQuantitiesDispatcher_) {
                    runGpuIntensiveQuantitiesDispatcher_(timeIdx);

                    const std::size_t numCells = this->intensiveQuantityCache_[timeIdx].size();
                    for (std::size_t i = 0; i < numCells; ++i) {
                        this->setIntensiveQuantitiesCacheEntryValidity(i, timeIdx, true);
                    }
                }
                else {
                    updateHasRun = 0; // reset so that the CPU update runs on the next call, if GPU dispatcher is not used
                }
            } else {
                updateHasRun = 0; // reset so that the CPU update runs on the next call, if GPU dispatcher is not supported
            }
#else
            updateHasRun = 0; // reset so that the CPU update runs on the next call, if GPU dispatcher is not supported
#endif
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
        const auto cpuStartTime = std::chrono::steady_clock::now();
#ifdef _OPENMP
#pragma omp parallel for
#endif
        for (const auto& chunk : element_chunks_) {
            for (const auto& elem : chunk) {
                this->template updateSingleCachedIntQuantUnchecked<Args...>(elementMapper.index(elem), timeIdx);
            }
        }
        const auto cpuDuration =
            std::chrono::duration_cast<std::chrono::milliseconds>(
                std::chrono::steady_clock::now() - cpuStartTime);
        Opm::OpmLog::info(std::format("updateCachedIntQuantsLoop CPU loop took {} ms",
                                      cpuDuration.count()));

#if HAVE_CUDA
        // After the CPU per-cell update has populated all fields, optionally
        // overlay the BlackOil intensive-quantities fields with their GPU
        // counterparts via the experimental dispatcher. The dispatcher only
        // overwrites the subset of fields covered by
        // BlackOilIntensiveQuantities::overlayBlackOilFieldsFrom; everything
        // else (mobility, energy, ...) keeps the CPU-computed value.
        if constexpr (Opm::gpuistl::GpuBlackoilIntensiveQuantitiesDispatcherSupport<TypeTag>::value)
        {
            if (useGpuIntensiveQuantitiesDispatcher_) {
                const auto gpuStartTime = std::chrono::steady_clock::now();
                runGpuIntensiveQuantitiesDispatcher_(timeIdx);
                const auto gpuDuration =
                    std::chrono::duration_cast<std::chrono::milliseconds>(
                        std::chrono::steady_clock::now() - gpuStartTime);
                Opm::OpmLog::info(std::format("updateCachedIntQuantsLoop GPU dispatch took {} ms",
                                              gpuDuration.count()));
            }
        }
#endif
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

#if HAVE_CUDA
    bool useGpuIntensiveQuantitiesDispatcher_{false};
    using GpuDispatcherStorage = std::conditional_t<
        Opm::gpuistl::GpuBlackoilIntensiveQuantitiesDispatcherSupport<TypeTag>::value,
        std::unique_ptr<Opm::gpuistl::GpuBlackoilIntensiveQuantitiesDispatcher<TypeTag>>,
        std::monostate>;
    mutable GpuDispatcherStorage gpuIntensiveQuantitiesDispatcher_{};

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

    mutable int updateHasRun = 0;
};

} // namespace Opm

#endif // FI_BLACK_OIL_MODEL_HPP
