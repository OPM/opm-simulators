// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*
  Copyright 2026, SINTEF Digital

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
 * \copydoc Opm::FICompositionalModel
 */
#ifndef OPM_FI_COMPOSITIONAL_MODEL_HPP
#define OPM_FI_COMPOSITIONAL_MODEL_HPP

#include <opm/models/ptflash/flashmodel.hh>
#include <opm/models/utils/propertysystem.hh>

#include <exception>
#include <mutex>

namespace Opm {

/*!
 * \brief PT-flash model that updates the cached intensive quantities cell by
 *        cell, without element contexts.
 */
template <class TypeTag>
class FICompositionalModel : public FlashModel<TypeTag>
{
    using ParentType = FlashModel<TypeTag>;
    using IntensiveQuantities = GetPropType<TypeTag, Properties::IntensiveQuantities>;
    using Simulator = GetPropType<TypeTag, Properties::Simulator>;

public:
    explicit FICompositionalModel(Simulator& simulator)
        : ParentType(simulator)
    {}

    void invalidateAndUpdateIntensiveQuantities(unsigned timeIdx) const
    {
        if (!this->storeIntensiveQuantities() ||
            timeIdx >= this->cachedIntensiveQuantityHistorySize()) {
            ParentType::invalidateAndUpdateIntensiveQuantities(timeIdx);
            return;
        }

        this->invalidateIntensiveQuantitiesCache(timeIdx);

        const auto& problem = this->simulator_.problem();
        const auto& solution = this->solution(timeIdx);
        const unsigned numDof = this->numGridDof();
        const auto flashSettings = IntensiveQuantities::FlashSettings::fromParameters();

        // An exception must not escape the parallel loop, which would call
        // std::terminate(); rethrow it afterwards, so that a failed flash
        // chops the time step.
        std::mutex exceptionLock;
        std::exception_ptr exceptionPtr = nullptr;

        // dynamic scheduling, since the flash costs more in two-phase cells
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 64)
#endif
        for (unsigned globalIdx = 0; globalIdx < numDof; ++globalIdx) {
            try {
                this->intensiveQuantityCache_[timeIdx][globalIdx]
                    .update(problem, solution[globalIdx], globalIdx, timeIdx, flashSettings);
                this->setIntensiveQuantitiesCacheEntryValidity(globalIdx, timeIdx, true);
            }
            catch (...) {
                const std::lock_guard lock(exceptionLock);
                exceptionPtr = std::current_exception();
            }
        }

        if (exceptionPtr) {
            std::rethrow_exception(exceptionPtr);
        }
    }
};

} // namespace Opm

#endif // OPM_FI_COMPOSITIONAL_MODEL_HPP
