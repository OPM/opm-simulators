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

#include <stdexcept>

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
    using Element = typename GridView::template Codim<0>::Entity;
    using ElementIterator = typename GridView::template Codim<0>::Iterator;
    enum {
        numEq = getPropValue<TypeTag, Properties::NumEq>(),
        historySize = getPropValue<TypeTag, Properties::TimeDiscHistorySize>(),
    };

public:
    FIBlackOilModel(Simulator& simulator)
        : BlackOilModel<TypeTag>(simulator)
    {
    }

    void invalidateAndUpdateIntensiveQuantities(unsigned timeIdx) const
    {

        this->invalidateIntensiveQuantitiesCache(timeIdx);
        OPM_BEGIN_PARALLEL_TRY_CATCH()
        // loop over all elements...
        ThreadedEntityIterator<GridView, /*codim=*/0> threadedElemIt(this->gridView_);
#ifdef _OPENMP
#pragma omp parallel
#endif
        {
            ElementContext elemCtx(this->simulator_);
            ElementIterator elemIt = threadedElemIt.beginParallel();
            for (; !threadedElemIt.isFinished(elemIt); elemIt = threadedElemIt.increment()) {
                const Element& elem = *elemIt;
                elemCtx.updatePrimaryStencil(elem);
                elemCtx.updatePrimaryIntensiveQuantities(timeIdx);
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
};
} // namespace Opm
#endif // FI_BLACK_OIL_MODEL_HPP
