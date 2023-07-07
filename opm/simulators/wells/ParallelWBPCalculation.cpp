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

#include <opm/simulators/wells/ParallelWBPCalculation.hpp>

#include <opm/simulators/wells/PerforationData.hpp>
#include <opm/simulators/wells/ParallelPAvgCalculator.hpp>
#include <opm/simulators/wells/ParallelWellInfo.hpp>

#include <opm/simulators/utils/ParallelCommunication.hpp>

#include <opm/input/eclipse/Schedule/Well/PAvgCalculator.hpp>
#include <opm/input/eclipse/Schedule/Well/PAvgCalculatorCollection.hpp>
#include <opm/input/eclipse/Schedule/Well/Well.hpp>
#include <opm/input/eclipse/Schedule/Well/WellConnections.hpp>

#include <opm/grid/common/CommunicationUtils.hpp>

#include <opm/input/eclipse/EclipseState/Grid/GridDims.hpp>

#include <opm/common/ErrorMacros.hpp>

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <memory>
#include <numeric>
#include <stdexcept>
#include <tuple>
#include <utility>
#include <vector>

Opm::ParallelWBPCalculation::
LocalConnSet::LocalConnSet(const std::vector<int>& localConnIdx)
    : localConnIdx_ { localConnIdx }
{}

int
Opm::ParallelWBPCalculation::
LocalConnSet::localIndex(const std::size_t connIdx) const
{
    return (connIdx >= this->localConnIdx_.size())
        ? -1
        : this->localConnIdx_[connIdx];
}

// ---------------------------------------------------------------------------

Opm::ParallelWBPCalculation::SourceData::SourceData(const Parallel::Communication& comm)
    : comm_ { comm }
{}

Opm::ParallelWBPCalculation::SourceData&
Opm::ParallelWBPCalculation::SourceData::
localIndex(GlobalToLocal localIdx)
{
    this->localIdx_ = std::move(localIdx);
    return *this;
}

Opm::ParallelWBPCalculation::SourceData&
Opm::ParallelWBPCalculation::SourceData::
evaluator(Evaluator eval)
{
    this->eval_ = std::move(eval);
    return *this;
}

Opm::ParallelWBPCalculation::SourceData&
Opm::ParallelWBPCalculation::SourceData::
evaluatorFactory(EvaluatorFactory evalFactory)
{
    this->evalFactory_ = std::move(evalFactory);
    return *this;
}

void
Opm::ParallelWBPCalculation::SourceData::
buildStructure(const std::vector<std::size_t>& sourceLocations)
{
    if (this->srcData_ == nullptr) {
        this->srcData_ = std::make_unique<ParallelPAvgDynamicSourceData>
            (this->comm_, sourceLocations, this->localIdx_);
    }
    else {
        this->srcData_->reconstruct(sourceLocations, this->localIdx_);
    }
}

void Opm::ParallelWBPCalculation::SourceData::collectDynamicValues()
{
    if (this->srcData_ == nullptr) {
        throw std::logic_error {
            "Cannot collect dynamic WBP source values "
            "prior to constructing source object"
        };
    }

    // For safety reasons, especially when this SourceData object pertains
    // to well connections and the user selects the "ALL" connection flag.
    this->srcData_->setToZero();

    if (this->eval_) {
        this->srcData_->collectLocalSources(this->eval_);
    }
    else if (this->evalFactory_) {
        this->srcData_->collectLocalSources(this->evalFactory_());
    }
    else {
        OPM_THROW(std::logic_error,
                  "Collecting WBP inputs requires an evaluation "
                  "function or an evaluation function factory");
    }

    this->srcData_->synchroniseSources();
}

std::vector<int>
Opm::ParallelWBPCalculation::SourceData::
getLocalIndex(const std::vector<std::size_t>& globalIndex) const
{
    auto localIdx = std::vector<int>(globalIndex.size());

    std::transform(globalIndex.begin(), globalIndex.end(), localIdx.begin(),
                   [this](const std::size_t globIx)
                   {
                       return this->localIdx_(globIx);
                   });

    return localIdx;
}

// ===========================================================================

Opm::ParallelWBPCalculation::ParallelWBPCalculation(const GridDims&                cellIndexMap,
                                                    const Parallel::Communication& gridComm)
    : cellIndexMap_{ cellIndexMap }
    , reservoirSrc_{ gridComm }
{}

Opm::ParallelWBPCalculation&
Opm::ParallelWBPCalculation::localCellIndex(GlobalToLocal localCellIdx)
{
    this->reservoirSrc_.localIndex(std::move(localCellIdx));
    return *this;
}

Opm::ParallelWBPCalculation&
Opm::ParallelWBPCalculation::evalCellSource(Evaluator evalCellSrc)
{
    this->reservoirSrc_.evaluator(std::move(evalCellSrc));
    return *this;
}

std::size_t
Opm::ParallelWBPCalculation::
createCalculator(const Well&             well,
                 const ParallelWellInfo& parallelWellInfo,
                 const std::vector<int>& localConnIdx,
                 EvaluatorFactory        makeWellSourceEvaluator)
{
    assert (this->wellConnSrc_.size() == this->localConnSet_.size());

    const auto ix = this->calculators_
        .setCalculator(well.seqIndex(), std::make_unique<ParallelPAvgCalculator>
                       (parallelWellInfo.communication(),
                        this->cellIndexMap_, well.getConnections()));

    assert (ix <= this->wellConnSrc_.size());

    if (ix == this->wellConnSrc_.size()) {
        // New calculator.
        this->wellConnSrc_.emplace_back(parallelWellInfo.communication());
        this->localConnSet_.emplace_back(localConnIdx);
    }
    else {
        // Update existing calculator.  Reset sources and connection sets.
        this->wellConnSrc_[ix] = SourceData { parallelWellInfo.communication() };
        this->localConnSet_[ix] = LocalConnSet { localConnIdx };
    }

    this->wellConnSrc_[ix]
        .evaluatorFactory(std::move(makeWellSourceEvaluator))
        .localIndex([ix = ix, this](const std::size_t connIdx)
        { return this->localConnSet_[ix].localIndex(connIdx); });

    return ix;
}

void Opm::ParallelWBPCalculation::defineCommunication()
{
    assert (this->calculators_.numCalculators() == this->wellConnSrc_.size());

    this->pruneInactiveWBPCells();

    this->defineReservoirCommunication();

    const auto numWells = this->calculators_.numCalculators();
    for (auto well = 0*numWells; well < numWells; ++well) {
        this->defineWellCommunication(well);
    }
}

void Opm::ParallelWBPCalculation::collectDynamicValues()
{
    this->reservoirSrc_.collectDynamicValues();

    for (auto& wellSrc : this->wellConnSrc_) {
        wellSrc.collectDynamicValues();
    }
}

void
Opm::ParallelWBPCalculation::
inferBlockAveragePressures(const std::size_t calcIndex,
                           const PAvg&       controls,
                           const double      gravity,
                           const double      refDepth)
{
    this->calculators_[calcIndex]
        .inferBlockAveragePressures(this->makeEvaluationSources(calcIndex),
                                    controls, gravity, refDepth);
}

const Opm::PAvgCalculator::Result&
Opm::ParallelWBPCalculation::averagePressures(const std::size_t calcIndex) const
{
    return this->calculators_[calcIndex].averagePressures();
}

void Opm::ParallelWBPCalculation::pruneInactiveWBPCells()
{
    if (this->reservoirSrc_.comm().size() == 1) {
        this->pruneInactiveWBPCellsSerial();
    }
    else {
        this->pruneInactiveWBPCellsParallel();
    }
}

void Opm::ParallelWBPCalculation::pruneInactiveWBPCellsSerial()
{
    this->calculators_.pruneInactiveWBPCells
        ([this](const std::vector<std::size_t>& globalWBPCellIdx)
    {
        auto isActive = std::vector<bool>(globalWBPCellIdx.size(), false);

        // Recall: localIndex() is negative if input is inactive or not on rank.
        const auto localWBPCellIdx =
            this->reservoirSrc_.getLocalIndex(globalWBPCellIdx);

        const auto nCells = isActive.size();
        for (auto cellIdx = 0*nCells; cellIdx < nCells; ++cellIdx) {
            isActive[cellIdx] = localWBPCellIdx[cellIdx] >= 0;
        }

        return isActive;
    });
}

void Opm::ParallelWBPCalculation::pruneInactiveWBPCellsParallel()
{
    this->calculators_.pruneInactiveWBPCells(
        [this](const std::vector<std::size_t>& globalWBPCellIdx)
    {
        auto isActive = std::vector<bool>(globalWBPCellIdx.size(), false);

        // AllWBPCells possibly contains repeated indices.  That's okay here.
        const auto& [allWBPCells, rankStart] =
            allGatherv(globalWBPCellIdx, this->reservoirSrc_.comm());

        // Recall: localIndex() is negative if input is inactive or not on rank.
        auto localWBPCellIdx =
            this->reservoirSrc_.getLocalIndex(allWBPCells);

        // The WBP cell is active if it has a non-negative local index on
        // one of the ranks.
        this->reservoirSrc_.comm()
            .max(localWBPCellIdx.data(), localWBPCellIdx.size());

        const auto myRank = this->reservoirSrc_.comm().rank();

        // Recall: This lambda function/callback is responsible only for the
        // range of 'allWBPCells' which applies to 'myRank'.  Consequently,
        // that's the range for which we determine the values in 'isActive'.
        // The rest of 'allWBPCells', or localWBPCellIdx for that matter, is
        // ignored and discarded when the lambda returns.
        auto cellIdx = 0*isActive.size();
        for (auto begin = localWBPCellIdx.begin() + rankStart[myRank + 0],
                 end    = localWBPCellIdx.begin() + rankStart[myRank + 1];
             begin != end; ++begin, ++cellIdx)
        {
            isActive[cellIdx] = *begin >= 0;
        }

        return isActive;
    });
}

void Opm::ParallelWBPCalculation::defineReservoirCommunication()
{
    auto sourceCells = std::vector<std::size_t>{};

    std::tie(sourceCells, std::ignore) =
        allGatherv(this->calculators_.allWBPCells(), this->reservoirSrc_.comm());

    std::sort(sourceCells.begin(), sourceCells.end());
    auto u = std::unique(sourceCells.begin(), sourceCells.end());

    this->reservoirSrc_.buildStructure({sourceCells.begin(), u});
}

void Opm::ParallelWBPCalculation::defineWellCommunication(const std::size_t well)
{
    this->wellConnSrc_[well]
        .buildStructure(this->calculators_[well].allWellConnections());
}

Opm::PAvgCalculator::Sources
Opm::ParallelWBPCalculation::makeEvaluationSources(const WellID well) const
{
    return PAvgCalculator::Sources{}
        .wellBlocks(this->reservoirSrc_)
        .wellConns (this->wellConnSrc_[well]);
}
