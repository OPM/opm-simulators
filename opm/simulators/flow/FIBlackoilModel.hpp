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

#include <opm/simulators/utils/DeferredLoggingErrorHelpers.hpp>

#include <cstddef>
#include <stdexcept>
#include <type_traits>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace Opm
{
template <typename TypeTag>
class FIBlackOilModel : public BlackOilModel<TypeTag>
{
    using ParentType = BlackOilModel<TypeTag>;
    using Simulator = GetPropType<TypeTag, Properties::Simulator>;
    using IntensiveQuantities = GetPropType<TypeTag, Properties::IntensiveQuantities>;
    using ElementContext = GetPropType<TypeTag, Properties::ElementContext>;
    using ThreadManager = GetPropType<TypeTag, Properties::ThreadManager>;
    using GridView = GetPropType<TypeTag, Properties::GridView>;
    using LocalFluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using Element = typename GridView::template Codim<0>::Entity;
    using ElementIterator = typename GridView::template Codim<0>::Iterator;
    enum {
        numEq = getPropValue<TypeTag, Properties::NumEq>(),
        historySize = getPropValue<TypeTag, Properties::TimeDiscHistorySize>(),
    };

    // The chunked and threaded iteration over elements in this class assumes that the number
    // and order of elements is fixed, and is therefore constrained to only work with CpGrid.
    // For example, ALUGrid supports refinement and can not assume that.
    static constexpr bool gridIsUnchanging = std::is_same_v<GetPropType<TypeTag, Properties::Grid>, Dune::CpGrid>;

public:
    FIBlackOilModel(Simulator& simulator)
        : BlackOilModel<TypeTag>(simulator)
    {
        if constexpr (gridIsUnchanging) {
            const auto& gv = this->gridView_;
#ifdef _OPENMP
            const int nt = omp_get_max_threads();
            if (nt > 1) {
                const auto num_elements = gv.size(0);
                constexpr int max_chunk_size = 1000;
                const int chunk_size = std::clamp(num_elements / nt, 1, max_chunk_size);
                OpmLog::debug("Using chunk size " + std::to_string(chunk_size) +
                              " for property evaluation with " + std::to_string(nt) + " OpenMP threads.");
                grid_chunk_iterators_.reserve(num_elements / chunk_size + 2);
                auto it = gv.template begin<0>();
                const auto end = gv.template end<0>();
                for (int count = 0; it != end; ++it, ++count) {
                    if (count % chunk_size == 0) {
                        grid_chunk_iterators_.push_back(it);
                    }
                }
                grid_chunk_iterators_.push_back(end);
            } else
#endif
            {
                // With one thread, or without OpenMP, we use a single chunk.
                grid_chunk_iterators_.push_back(gv.template begin<0>());
                grid_chunk_iterators_.push_back(gv.template end<0>());
            }
        }
    }

    void invalidateAndUpdateIntensiveQuantities(unsigned timeIdx) const
    {

        using DynamicFluidSystem = std::remove_reference_t<decltype(LocalFluidSystem::getNonStatic())>;
        // static constexpr bool use_dynamic_fluidsystem = is_a_dynamic_blackoil_system<DynamicFluidSystem>;

        DynamicFluidSystem* fluidSystemInstance = &LocalFluidSystem::getNonStatic();
        // if constexpr (use_dynamic_fluidsystem) {
        //     fluidSystemInstance = &LocalFluidSystem::getNonStatic();
        // }

        this->invalidateIntensiveQuantitiesCache(timeIdx);
        OPM_BEGIN_PARALLEL_TRY_CATCH();
        if constexpr (gridIsUnchanging) {
            const int num_chunks = grid_chunk_iterators_.size() - 1;
#ifdef _OPENMP
#pragma omp parallel for
#endif
            for (int chunk = 0; chunk < num_chunks; ++chunk) {
                ElementContext elemCtx(this->simulator_);
                for (auto it = grid_chunk_iterators_[chunk]; it != grid_chunk_iterators_[chunk+1]; ++it) {
                    const Element& elem = *it;
                    elemCtx.updatePrimaryStencil(elem);

                    elemCtx.updatePrimaryIntensiveQuantities(timeIdx, *fluidSystemInstance);
                }
            }
        } else {
            // Grid is possibly refined or otherwise changed between calls.
            const auto& gv = this->gridView_;
            auto it = gv.template begin<0>();
            const auto end = gv.template end<0>();
            ElementContext elemCtx(this->simulator_);
            for (; it != end; ++it) {
                const Element& elem = *it;
                elemCtx.updatePrimaryStencil(elem);
                elemCtx.updatePrimaryIntensiveQuantities(timeIdx, *fluidSystemInstance);
            }
        }
        OPM_END_PARALLEL_TRY_CATCH("InvalideAndUpdateIntensiveQuantities: state error", this->simulator_.vanguard().grid().comm());
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
        OPM_END_PARALLEL_TRY_CATCH("InvalideAndUpdateIntensiveQuantitiesOverlap: state error", this->simulator_.vanguard().grid().comm());
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
                      "Run without intensive quantites not enabled: Use --enable-intensive-quantity=true");
        }
        const auto* intquant = this->cachedIntensiveQuantities(globalIdx, timeIdx);
        if (!intquant) {
            OPM_THROW(std::logic_error, "Intensive quantites need to be updated in code");
        }
        return *intquant;
    }

protected:
    std::vector<ElementIterator> grid_chunk_iterators_;
};
} // namespace Opm
#endif // FI_BLACK_OIL_MODEL_HPP
