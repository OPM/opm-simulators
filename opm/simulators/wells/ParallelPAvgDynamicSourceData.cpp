/*
  Copyright 2023 Equinor ASA.

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

#include <config.h>

#include <opm/simulators/wells/ParallelPAvgDynamicSourceData.hpp>

#include <opm/simulators/utils/ParallelCommunication.hpp>

#include <opm/grid/common/CommunicationUtils.hpp>

#include <algorithm>
#include <cstddef>
#include <functional>
#include <iterator>
#include <numeric>
#include <utility>
#include <vector>

template<class Scalar>
Opm::ParallelPAvgDynamicSourceData<Scalar>::
ParallelPAvgDynamicSourceData(const Parallel::Communication&  comm,
                              const std::vector<std::size_t>& sourceLocations,
                              GlobalToLocal                   localCellIdx)
    : PAvgDynamicSourceData<Scalar> { sourceLocations }
    , comm_                 { comm }
{
    this->finaliseConstruction(sourceLocations, std::move(localCellIdx));
}

template<class Scalar>
void Opm::ParallelPAvgDynamicSourceData<Scalar>::setToZero()
{
    std::fill_n(this->localSrc_.begin(), this->localSrc_.size(), 0.0);
}

template<class Scalar>
void Opm::ParallelPAvgDynamicSourceData<Scalar>::
reconstruct(const std::vector<std::size_t>& sourceLocations,
            GlobalToLocal                   localCellIdx)
{
    PAvgDynamicSourceData<Scalar>::reconstruct(sourceLocations); // Reconstruct base
    this->finaliseConstruction(sourceLocations, std::move(localCellIdx));
}

template<class Scalar>
void Opm::ParallelPAvgDynamicSourceData<Scalar>::collectLocalSources(Evaluator eval)
{
    auto localIx = std::size_t{0};

    for (const auto& location : this->locations_) {
        eval(location.cell, this->localSourceTerm(localIx++));
    }
}

template<class Scalar>
void Opm::ParallelPAvgDynamicSourceData<Scalar>::synchroniseSources()
{
    this->comm_.get()
        .allgatherv(this->localSrc_.data(),       // Input (from)
                    static_cast<int>(this->localSrc_.size()),
                    this->src_.data(),            // Output (to)
                    this->allSizes_.data(),       // #elements per rank
                    this->startPointers_.data()); // Output offsets
}

template<class Scalar>
typename std::vector<Scalar>::size_type
Opm::ParallelPAvgDynamicSourceData<Scalar>::
storageIndex(const typename std::vector<Scalar>::size_type elemIndex) const
{
    return this->storageIndex_[elemIndex];
}

template<class Scalar>
void Opm::ParallelPAvgDynamicSourceData<Scalar>::
finaliseConstruction(const std::vector<std::size_t>& sourceLocations,
                     GlobalToLocal                   localCellIdx)
{
    auto ix = std::size_t{0};

    this->locations_.clear();

    for (const auto& location : sourceLocations) {
        if (const auto cell = localCellIdx(location); cell >= 0) {
            this->locations_.push_back({ ix, cell });
        }

        ix += 1;
    }

    this->localSrc_.assign(this->numSpanItems() * this->locations_.size(), 0.0);

    this->defineCommunication();
}

template<class Scalar>
typename Opm::PAvgDynamicSourceData<Scalar>::template SourceDataSpan<Scalar>
Opm::ParallelPAvgDynamicSourceData<Scalar>::localSourceTerm(const std::size_t localIx)
{
    return this->sourceTerm(localIx, this->localSrc_);
}

template<class Scalar>
void Opm::ParallelPAvgDynamicSourceData<Scalar>::defineCommunication()
{
    // 1) Determine origins/owning ranks for all source terms.
    auto ixVec = std::vector<std::size_t>(this->locations_.size());
    std::transform(this->locations_.begin(), this->locations_.end(),
                   ixVec.begin(),
                   [](const auto& location) { return location.ix; });

    constexpr auto numItems = ParallelPAvgDynamicSourceData<Scalar>::numSpanItems();

    const auto& [allIndices, allIxStart] = allGatherv(ixVec, this->comm_.get());

    // -----------------------------------------------------------------------

    // 2) Determine starting pointers/offsets/displacements for received
    //    basic elements from each rank.  There are 'numItems' basic data
    //    elements for each source term.
    this->startPointers_.resize(allIxStart.size());
    std::transform(allIxStart.begin(), allIxStart.end(),
                   this->startPointers_.begin(),
                   [](const int start)
                   {
                       return numItems * start;
                   });

    // -----------------------------------------------------------------------

    // 3) Determine number of basic data elements to receive from each rank.
    this->allSizes_.resize(allIxStart.size() - 1);
    std::adjacent_difference(this->startPointers_.begin() + 1,
                             this->startPointers_.end(),
                             this->allSizes_.begin());

    // -----------------------------------------------------------------------

    // 4) Build translation mapping from source term element indices to
    //    storage indices.
    //
    //    Note that if the source terms aren't all active--e.g., if a well
    //    is connected in deactivated cells--then allIndices is not a
    //    permutation of 0..allIndices.size()-1 and the maximum source
    //    location may exceed size()-1.  Resize the storageIndex_ according
    //    to the largest source location ID.
    if (auto maxIxPos = std::max_element(allIndices.begin(), allIndices.end());
        maxIxPos != allIndices.end())
    {
        // +1 for zero-based indices.
        this->storageIndex_.resize(*maxIxPos + 1);
    }

    auto storageIx = std::vector<double>::size_type{0};
    for (const auto& elemIndex : allIndices) {
        this->storageIndex_[elemIndex] = storageIx++;
    }
}

template class Opm::ParallelPAvgDynamicSourceData<double>;
