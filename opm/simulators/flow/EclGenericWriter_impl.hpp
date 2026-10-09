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
#ifndef OPM_ECL_GENERIC_WRITER_IMPL_HPP
#define OPM_ECL_GENERIC_WRITER_IMPL_HPP

#include <dune/grid/common/mcmgmapper.hh>

#include <opm/grid/cpgrid/LgrOutputHelpers.hpp>
#include <opm/grid/GridHelpers.hpp>
#include <opm/grid/utility/cartesianToCompressed.hpp>

#include <opm/common/OpmLog/OpmLog.hpp>

#include <opm/input/eclipse/EclipseState/EclipseState.hpp>
#include <opm/input/eclipse/EclipseState/Grid/RegionSetMatcher.hpp>
#include <opm/input/eclipse/EclipseState/Grid/NNC.hpp>
#include <opm/input/eclipse/EclipseState/SummaryConfig/SummaryConfig.hpp>

#include <opm/input/eclipse/Schedule/Action/State.hpp>
#include <opm/input/eclipse/Schedule/RPTConfig.hpp>
#include <opm/input/eclipse/Schedule/Schedule.hpp>
#include <opm/input/eclipse/Schedule/SummaryState.hpp>
#include <opm/input/eclipse/Schedule/UDQ/UDQConfig.hpp>
#include <opm/input/eclipse/Schedule/UDQ/UDQState.hpp>
#include <opm/input/eclipse/Schedule/Well/WellConnections.hpp>
#include <opm/input/eclipse/Schedule/Well/WellMatcher.hpp>

#include <opm/input/eclipse/Units/UnitSystem.hpp>

#include <opm/output/data/RegionVariableMapping.hpp>

#include <opm/output/eclipse/EclipseIO.hpp>
#include <opm/output/eclipse/RegionVariableCollection.hpp>
#include <opm/output/eclipse/RestartValue.hpp>
#include <opm/output/eclipse/Summary.hpp>

#include <opm/simulators/flow/EclGenericWriter.hpp>

#include <opm/simulators/utils/DeferredLoggingErrorHelpers.hpp>

#if HAVE_MPI
#include <opm/simulators/utils/MPISerializer.hpp>
#endif

#if HAVE_MPI
#include <mpi.h>
#endif

#include <fmt/format.h>

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <functional>
#include <map>
#include <memory>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace {

/*!
 * \brief Restart values for a declared LGR that the current grid does not refine.
 *
 * The EGRID and INIT describe every LGR the *deck* declares, whether or not the
 * simulation grid refines it. A restart step that omits the corresponding
 * solution section therefore covers fewer cells than the grid has, and
 * post-processors that size results from the grid reject the result outright --
 * ResInsight requires each keyword's total value count to be a whole multiple of
 * the active cell count (RifEclipseOutputFileTools::validKeywordsForPorosityModel).
 *
 * So emit a section for the unrefined LGR too, giving every refined cell the
 * value of its father coarse cell: the same father mapping, and the same
 * father-replication, that the INIT writer already applies to the static
 * properties in this situation (WriteInit.cpp, the !fullProperties branch).
 */
Opm::data::Solution
fatherReplicatedSolution(const Opm::EclipseGrid&     inputGrid,
                         const std::string&          lgrName,
                         const Opm::data::Solution&  coarse)
{
    const auto& lgrGrid = inputGrid.getLGRCell(lgrName);
    // The coarse solution is active-sized, so the father mapping must be too.
    const auto fathers = lgrGrid.getLGRCell_active_father(inputGrid);
    const auto nActive = inputGrid.getNumActive();

    // Copy first: keeps the keys, units, targets and the SI flag of the
    // coarse solution, so only the per-cell arrays need replacing.
    auto refined = coarse;

    for (const auto& [key, cellData] : coarse) {
        cellData.visit([&refined, &key = key, &cellData, &fathers, nActive](const auto& src)
        {
            using Vector = std::decay_t<decltype(src)>;
            if constexpr (! std::is_same_v<Vector, std::monostate>) {
                if (src.size() != nActive) {
                    return;     // not per-cell data; leave it alone
                }
                Vector dst(fathers.size());
                for (auto i = 0*fathers.size(); i < fathers.size(); ++i) {
                    dst[i] = src[fathers[i]];
                }
                if constexpr (std::is_same_v<Vector, std::vector<double>>) {
                    refined[key] = Opm::data::CellData {
                        cellData.dim, std::move(dst), cellData.target
                    };
                }
                else {
                    refined[key] = Opm::data::CellData {
                        std::move(dst), cellData.target
                    };
                }
            }
        });
    }

    return refined;
}

/*!
 * \brief Detect whether two cells are direct vertical neighbours.
 *
 * I.e. have the same i and j index and all cartesian cells between them
 * along the vertical column are inactive.
 *
 * \tparam CM The type of the cartesian index mapper.
 * \param cartMapper The mapper onto cartesian indices.
 * \param cartesianToActive The mapping of cartesian indices to active indices.
 * \param smallGlobalIndex The cartesian cell index of the cell with smaller index
 * \param largeGlobalIndex The cartesian cell index of the cell with larger index
 * \return True if the cells have the same i and j indices and all cartesian cells
 *         between them are inactive.
 */
bool directVerticalNeighbors(const std::array<int, 3>& cartDims,
                             const std::unordered_map<int,int>& cartesianToActive,
                             int smallGlobalIndex, int largeGlobalIndex)
{
    assert(smallGlobalIndex <= largeGlobalIndex);
    std::array<int, 3> ijk1, ijk2;
    auto globalToIjk = [cartDims](int gc) {
                           std::array<int, 3> ijk;
                           ijk[0] = gc % cartDims[0];
                           gc /= cartDims[0];
                           ijk[1] = gc % cartDims[1];
                           ijk[2] = gc / cartDims[1];
                           return ijk;
                       };
    ijk1 = globalToIjk(smallGlobalIndex);
    ijk2 = globalToIjk(largeGlobalIndex);
    assert(ijk2[2]>=ijk1[2]);

    if ( ijk1[0] == ijk2[0] && ijk1[1] == ijk2[1] && (ijk2[2] - ijk1[2]) > 1)
    {
        assert((largeGlobalIndex-smallGlobalIndex)%(cartDims[0]*cartDims[1])==0);
        for ( int gi = smallGlobalIndex + cartDims[0] * cartDims[1]; gi < largeGlobalIndex;
              gi += cartDims[0] * cartDims[1] )
        {
            if ( cartesianToActive.find( gi ) != cartesianToActive.end() )
            {
                return false;
            }
        }
        return true;
    } else
        return false;
}

std::unordered_map<std::string, Opm::data::InterRegFlowMap>
getInterRegFlowsAsMap(const Opm::InterRegFlowMap& map)
{
    auto maps = std::unordered_map<std::string, Opm::data::InterRegFlowMap>{};

    const auto& regionNames = map.names();
    auto flows = map.getInterRegFlows();
    const auto nmap = regionNames.size();

    maps.reserve(nmap);
    for (auto mapID = 0*nmap; mapID < nmap; ++mapID) {
        maps.emplace(regionNames[mapID], std::move(flows[mapID]));
    }

    return maps;
}

struct EclWriteTasklet : public Opm::TaskletInterface
{
    Opm::Action::State actionState_;
    Opm::WellTestState wtestState_;
    Opm::SummaryState summaryState_;
    Opm::UDQState udqState_;
    Opm::EclipseIO& eclIO_;
    int reportStepNum_;
    std::optional<int> timeStepNum_;
    bool isSubStep_;
    double secondsElapsed_;
    std::vector<Opm::RestartValue> restartValue_;
    bool writeDoublePrecision_;
    /// \brief True if there was an EXIT keyword in ACTIONX causing a simulation end
    bool forcedSimulationFinished_;

    explicit EclWriteTasklet(const Opm::Action::State& actionState,
                             const Opm::WellTestState& wtestState,
                             const Opm::SummaryState& summaryState,
                             const Opm::UDQState& udqState,
                             Opm::EclipseIO& eclIO,
                             int reportStepNum,
                             std::optional<int> timeStepNum,
                             bool isSubStep,
                             double secondsElapsed,
                             std::vector<Opm::RestartValue> restartValue,
                             bool writeDoublePrecision,
                             bool forcedSimulationFinished)
        : actionState_(actionState)
        , wtestState_(wtestState)
        , summaryState_(summaryState)
        , udqState_(udqState)
        , eclIO_(eclIO)
        , reportStepNum_(reportStepNum)
        , timeStepNum_(timeStepNum)
        , isSubStep_(isSubStep)
        , secondsElapsed_(secondsElapsed)
        , restartValue_(std::move(restartValue))
        , writeDoublePrecision_(writeDoublePrecision)
        , forcedSimulationFinished_(forcedSimulationFinished)
    {}

    // callback to eclIO serial writeTimeStep method
    void run() override
    {
        if (this->restartValue_.size() == 1) {
            this->eclIO_.writeTimeStep(this->actionState_,
                                       this->wtestState_,
                                       this->summaryState_,
                                       this->udqState_,
                                       this->reportStepNum_,
                                       this->isSubStep_,
                                       this->secondsElapsed_,
                                       std::move(this->restartValue_.back()),
                                       this->writeDoublePrecision_,
                                       this->timeStepNum_,
                                       forcedSimulationFinished_);
        }
        else{
            this->eclIO_.writeTimeStep(this->actionState_,
                                       this->wtestState_,
                                       this->summaryState_,
                                       this->udqState_,
                                       this->reportStepNum_,
                                       this->isSubStep_,
                                       this->secondsElapsed_,
                                       std::move(this->restartValue_),
                                       this->writeDoublePrecision_,
                                       this->timeStepNum_,
                                       forcedSimulationFinished_);
        }
    }
};

}

namespace Opm {

template<class Grid, class EquilGrid, class GridView, class ElementMapper, class Scalar>
EclGenericWriter<Grid,EquilGrid,GridView,ElementMapper,Scalar>::
EclGenericWriter(const Schedule& schedule,
                 const EclipseState& eclState,
                 const SummaryConfig& summaryConfig,
                 const Grid& grid,
                 const EquilGrid* equilGrid,
                 const GridView& gridView,
                 const Dune::CartesianIndexMapper<Grid>& cartMapper,
                 const Dune::CartesianIndexMapper<EquilGrid>* equilCartMapper,
                 bool enableAsyncOutput,
                 bool enableEsmry,
                 const EquilGrid* collectGrid,
                 const Dune::CartesianIndexMapper<EquilGrid>* collectCartMapper )
    : collectOnIORank_(grid,
                       collectGrid ? collectGrid : equilGrid,
                       gridView,
                       cartMapper,
                       collectCartMapper ? collectCartMapper : equilCartMapper,
                       summaryConfig.fip_regions_interreg_flow())
    , grid_           (grid)
    , gridView_       (gridView)
    , schedule_       (schedule)
    , eclState_       (eclState)
    , cartMapper_     (cartMapper)
    , equilCartMapper_(equilCartMapper)
    , collectCartMapper_(collectCartMapper ? collectCartMapper : equilCartMapper)
    , equilGrid_      (equilGrid)
    , collectGrid_    (collectGrid ? collectGrid : equilGrid)
{
    // Make sure outputNnc_ vector has at least 1 entry in all ranks.
    outputNnc_.resize(1);

    if (this->collectOnIORank_.isIORank()) {
        this->eclIO_ = std::make_unique<EclipseIO>
            (this->eclState_,
             UgGridHelpers::createEclipseGrid(*equilGrid, eclState_.getInputGrid()),
             this->schedule_, summaryConfig, "", enableEsmry);
    }

    // create output thread if enabled and rank is I/O rank
    // async output is enabled by default if pthread are enabled
    int numWorkerThreads = 0;
    if (enableAsyncOutput && collectOnIORank_.isIORank()) {
        numWorkerThreads = 1;
    }

    this->taskletRunner_.reset(new TaskletRunner(numWorkerThreads));
}

template<class Grid, class EquilGrid, class GridView, class ElementMapper, class Scalar>
const EclipseIO& EclGenericWriter<Grid,EquilGrid,GridView,ElementMapper,Scalar>::
eclIO() const
{
    assert(eclIO_);
    return *eclIO_;
}

template<class Grid, class EquilGrid, class GridView, class ElementMapper, class Scalar>
void EclGenericWriter<Grid,EquilGrid,GridView,ElementMapper,Scalar>::
writeInit()
{
    if (collectOnIORank_.isIORank()) {
        std::map<std::string, std::vector<int>> integerVectors;
        // globalRanks() is empty when the I/O-rank cell collection is not set up.
        if (collectOnIORank_.isParallel() && !collectOnIORank_.globalRanks().empty()) {
            const auto& leafRanks = collectOnIORank_.globalRanks();
            if (collectGrid_ != nullptr && collectGrid_->maxLevel() > 0) {
                // For a refined grid the gather reference (collectGrid_) is the
                // leaf grid, so globalRanks_ is leaf-sized.  MPI_RANK in the INIT
                // belongs to the level-zero (main) grid, so reduce it: each leaf
                // cell contributes its owning rank to its level-zero origin (all
                // children of a refined cell share one rank - the box is
                // rank-interior).  Without this MPI_RANK is written at leaf size
                // into the main-grid slot, mismatching the main grid.
                std::vector<int> mpiRank(collectGrid_->currentData().front()->size(0), 0);
                const auto leafView = collectGrid_->leafGridView();
                Dune::MultipleCodimMultipleGeomTypeMapper<std::decay_t<decltype(leafView)>>
                    leafMapper(leafView, Dune::mcmgElementLayout());
                for (const auto& elem : elements(leafView)) {
                    const auto leafIdx = leafMapper.index(elem);
                    const auto originIdx = elem.getOrigin().index();
                    if (originIdx >= 0
                        && static_cast<std::size_t>(originIdx) < mpiRank.size()
                        && static_cast<std::size_t>(leafIdx) < leafRanks.size()) {
                        mpiRank[originIdx] = leafRanks[leafIdx];
                    }
                }
                integerVectors.emplace("MPI_RANK", std::move(mpiRank));
            }
            else {
                integerVectors.emplace("MPI_RANK", leafRanks);
            }
        }

        if (const auto& lgrs = this->eclState_.getLgrs(); lgrs.size() > 0) {

            const auto nncCollection = Opm::NNCCollection::fromLGROutputContainers(this->outputNnc_,
                                                                                   this->outputNncGlobalLocal_,
                                                                                   this->outputAmalgamatedNnc_);
            eclIO_->writeInitial(*this->outputTrans_,
                                 integerVectors,
                                 nncCollection);
        } else {
            eclIO_->writeInitial(*this->outputTrans_,
                                 integerVectors,
                                 this->outputNnc_.front());
        }
        this->outputTrans_.reset();
    }
    this->gatheredLgrTrans_.reset();
}

template<class Grid, class EquilGrid, class GridView, class ElementMapper, class Scalar>
void
EclGenericWriter<Grid,EquilGrid,GridView,ElementMapper,Scalar>::
extractOutputTransAndNNC(const std::function<unsigned int(unsigned int)>& map)
{
    // The work below runs on the I/O rank only, but an exception there (e.g. a
    // missing gathered LGR transmissibility) must not leave the other ranks
    // blocked in the broadcast that follows -- so failure is agreed on
    // collectively and every rank throws together.
    OPM_BEGIN_PARALLEL_TRY_CATCH()

    if (collectOnIORank_.isIORank()) {
        constexpr bool equilGridIsCpGrid = std::is_same_v<EquilGrid, Dune::CpGrid>;

        const auto levelCartMapp = this->createLevelCartMapp_<equilGridIsCpGrid>();
        const auto levelCartToLevelCompressed = this->createCartesianToActiveMaps_<equilGridIsCpGrid>(levelCartMapp);
        auto computeLevelIndices = this->computeLevelIndices_<equilGridIsCpGrid>();
        auto computeLevelCartIdx = this->computeLevelCartIdx_<equilGridIsCpGrid>(levelCartMapp, *(this->collectCartMapper_));
        auto computeLevelCartDimensions = this->computeLevelCartDimensions_<equilGridIsCpGrid>(levelCartMapp, *(this->collectCartMapper_));
        auto computeOriginIndices = this->computeOriginIndices_<equilGridIsCpGrid>();

        computeTrans_(levelCartToLevelCompressed, map, computeLevelIndices,
                      computeLevelCartIdx, computeLevelCartDimensions, computeOriginIndices);
        exportNncStructure_(levelCartToLevelCompressed, map, computeLevelIndices, computeLevelCartIdx,
                            computeLevelCartDimensions, computeOriginIndices);
    }

    OPM_END_PARALLEL_TRY_CATCH("EclGenericWriter::extractOutputTransAndNNC() failed: ",
                               grid_.comm());

#if HAVE_MPI
    if (collectOnIORank_.isParallel()) {
        const auto& comm = grid_.comm();
        Parallel::MpiSerializer ser(comm);
        ser.broadcast(Parallel::RootRank{0}, outputNnc_);
    }
#endif
}

template<class Grid, class EquilGrid, class GridView, class ElementMapper, class Scalar>
bool
EclGenericWriter<Grid,EquilGrid,GridView,ElementMapper,Scalar>::
isNumAquCell_(const std::size_t cartIdx) const
{
    const auto& numAquCell = this->eclState_.aquifer().hasNumericalAquifer()
        ? this->eclState_.aquifer().numericalAquifers().allAquiferCellIds()
        : std::vector<std::size_t>{};

    return std::ranges::binary_search(numAquCell.begin(), numAquCell.end(), cartIdx);
}

template<class Grid, class EquilGrid, class GridView, class ElementMapper, class Scalar>
bool
EclGenericWriter<Grid,EquilGrid,GridView,ElementMapper,Scalar>::
isNumAquConn_(const std::size_t cartIdx1,
              const std::size_t cartIdx2) const
{
    return isNumAquCell_(cartIdx1) || isNumAquCell_(cartIdx2);
}

template<class Grid, class EquilGrid, class GridView, class ElementMapper, class Scalar>
template<bool equilGridIsCpGrid>
Opm::LevelCartesianIndexMapper<EquilGrid>
EclGenericWriter<Grid,EquilGrid,GridView,ElementMapper,Scalar>::
createLevelCartMapp_() const
{
    if constexpr (equilGridIsCpGrid) {
        return Opm::LevelCartesianIndexMapper<EquilGrid>(*this->collectGrid_);
    } else {
        return Opm::LevelCartesianIndexMapper<EquilGrid>(*collectCartMapper_); }
}

template<class Grid, class EquilGrid, class GridView, class ElementMapper, class Scalar>
template<bool equilGridIsCpGrid>
std::vector<std::unordered_map<int,int>>
EclGenericWriter<Grid,EquilGrid,GridView,ElementMapper,Scalar>::
createCartesianToActiveMaps_(const Opm::LevelCartesianIndexMapper<EquilGrid>& levelCartMapp) const
{
    if constexpr (equilGridIsCpGrid) {
        if (this->collectGrid_->maxLevel()) {
            return Opm::Lgr::levelCartesianToLevelCompressedMaps(*this->collectGrid_, levelCartMapp); }
        else {
            return std::vector<std::unordered_map<int,int>>{ cartesianToCompressed(collectGrid_->size(0), UgGridHelpers::globalCell(*collectGrid_)) };
        }
    }
    return std::vector<std::unordered_map<int,int>>{ cartesianToCompressed(collectGrid_->size(0), UgGridHelpers::globalCell(*collectGrid_)) };
}

template<class Grid, class EquilGrid, class GridView, class ElementMapper, class Scalar>
template<bool equilGridIsCpGrid>
std::function<std::array<int,3>(int)>
EclGenericWriter<Grid,EquilGrid,GridView,ElementMapper,Scalar>::
computeLevelCartDimensions_(const Opm::LevelCartesianIndexMapper<EquilGrid>& levelCartMapp,
                            const Dune::CartesianIndexMapper<EquilGrid>& equilCartMapp) const
{
    if constexpr (equilGridIsCpGrid) {
        return [&](int level)
        {
            return levelCartMapp.cartesianDimensions(level);
        };
    }
    else {
        return [&]([[maybe_unused]] int level)
        {
            assert(level == 0); // refinement only supported for CpGrid for now
            return equilCartMapp.cartesianDimensions();
        };
    }
}

template<class Grid, class EquilGrid, class GridView, class ElementMapper, class Scalar>
template<bool equilGridIsCpGrid>
std::function<int(int, int)>
EclGenericWriter<Grid,EquilGrid,GridView,ElementMapper,Scalar>::
computeLevelCartIdx_(const Opm::LevelCartesianIndexMapper<EquilGrid>& levelCartMapp,
                     const Dune::CartesianIndexMapper<EquilGrid>& equilCartMapp) const
{
    if constexpr (equilGridIsCpGrid) {
        return [&](int levelCompressedIdx,
                   int level)
        {
            return levelCartMapp.cartesianIndex(levelCompressedIdx, level);
        };
    }
    else {
        return [&](int levelCompressedIdx,
                   [[maybe_unused]] int level)
        {
            assert(level == 0); // refinement only supported for CpGrid for now
            return equilCartMapp.cartesianIndex(levelCompressedIdx);
        };
    }
}

template<class Grid, class EquilGrid, class GridView, class ElementMapper, class Scalar>
template <bool equilGridIsCpGrid>
auto
EclGenericWriter<Grid,EquilGrid,GridView,ElementMapper,Scalar>::
computeLevelIndices_() const
{
    if constexpr (equilGridIsCpGrid) {
        return [](const auto& intersection,
                  const auto&,
                  const auto&)
        {
            return std::pair{intersection.inside().getLevelElem().index(), intersection.outside().getLevelElem().index()};
        };
    }
    else {
        return [](const auto&,
                  const auto& intersectionInsideLeafIdx,
                  const auto& intersectionOutsideLeafIdx)
        {
            return std::pair{intersectionInsideLeafIdx, intersectionOutsideLeafIdx};
        };
    }
}

template<class Grid, class EquilGrid, class GridView, class ElementMapper, class Scalar>
template <bool equilGridIsCpGrid>
auto
EclGenericWriter<Grid,EquilGrid,GridView,ElementMapper,Scalar>::
computeOriginIndices_() const
{
    if constexpr (equilGridIsCpGrid) {
        return [](const auto& intersection,
                  const auto&,
                  const auto&)
        {
            return std::pair{intersection.inside().getOrigin().index(), intersection.outside().getOrigin().index()};
        };
    }
    else {
        return [](const auto&,
                  const auto& intersectionInsideLeafIdx,
                  const auto& intersectionOutsideLeafIdx)
        {
            return std::pair{intersectionInsideLeafIdx, intersectionOutsideLeafIdx};
        };
    }
}

template<class Grid, class EquilGrid, class GridView, class ElementMapper, class Scalar>
void
EclGenericWriter<Grid,EquilGrid,GridView,ElementMapper,Scalar>::
allocateLevelTrans_(const std::array<int,3>& levelCartDims,
                    data::Solution& levelTrans) const
{
    auto createLevelCellData = [&levelCartDims]() {
        return Opm::data::CellData{
            Opm::UnitSystem::measure::transmissibility,
            std::vector<double>(levelCartDims[0] * levelCartDims[1] * levelCartDims[2], 0.0),
            Opm::data::TargetType::INIT
        };
    };

    levelTrans.clear();
    levelTrans.emplace("TRANX", createLevelCellData());
    levelTrans.emplace("TRANY", createLevelCellData());
    levelTrans.emplace("TRANZ", createLevelCellData());
}

template<class Grid, class EquilGrid, class GridView, class ElementMapper, class Scalar>
template<typename LevelIndicesFunction, typename OriginIndicesFunction>
void
EclGenericWriter<Grid,EquilGrid,GridView,ElementMapper,Scalar>::
computeTrans_(const std::vector<std::unordered_map<int,int>>&  levelCartToLevelCompressed,
              const std::function<unsigned int(unsigned int)>& map,
              const LevelIndicesFunction& computeLevelIndices,
              const std::function<int(int, int)>& computeLevelCartIdx,
              const std::function<std::array<int,3>(int)>& computeLevelCartDims,
              const OriginIndicesFunction& computeOriginIndices) const
{
    if (!outputTrans_) {
        outputTrans_ = std::make_unique<std::vector<data::Solution>>(std::vector<data::Solution>{});
    }

    using GlobalGridView = typename EquilGrid::LeafGridView;
    using GlobElementMapper = Dune::MultipleCodimMultipleGeomTypeMapper<GlobalGridView>;
    const GlobalGridView& globalGridView = this->collectGrid_->leafGridView();
    const GlobElementMapper globalElemMapper { globalGridView, Dune::mcmgElementLayout() };

    // Refinement supported only for CpGrid for now.
    int maxLevel = this->collectGrid_->maxLevel();

    outputTrans_->resize(maxLevel+1); // including level zero grid

    for (int level = 0; level <= maxLevel; ++level) {
        allocateLevelTrans_(computeLevelCartDims(level), this->outputTrans_->at(level));
    }

    for (const auto& elem : elements(globalGridView)) {
        for (const auto& is : intersections(globalGridView, elem)) {
            if (!is.neighbor())
                continue; // intersection is on the domain boundary

            if ( is.inside().level() != is.outside().level() ) // Those are treated as NNCs
                continue;

            // Not 'const' because remapped if 'map' is non-null.
            unsigned c1 = globalElemMapper.index(is.inside());
            unsigned c2 = globalElemMapper.index(is.outside());

            if (c1 > c2)
                continue; // we only need to handle each connection once, thank you.

            int level = is.inside().level();

            // For CpGrid with LGRs, level*Idx and c* do not coincide.
            const auto& [levelInIdx, levelOutIdx] = computeLevelIndices(is, c1, c2);

            const int levelCartIdxIn = computeLevelCartIdx(levelInIdx, level);
            const int levelCartIdxOut = computeLevelCartIdx(levelOutIdx, level);

            // For CpGrid with LGRs, the origin cell index refers to the coarsest
            // ancestor cell when the cell is refined. For cells not involved in
            // any refinement, it corresponds to the geometrically equivalent
            // cell in the level-zero grid.
            const auto [originInIdx, originOutIdx] = computeOriginIndices(is, c1, c2);

            const auto originCartIdxIn = computeLevelCartIdx(originInIdx, /* level = */ 0);
            const auto originCartIdxOut = computeLevelCartIdx(originOutIdx, /* level = */ 0);

            // For level-zero grid, level Cartesian indices coincide with the grid Cartesian indices.
            if (isNumAquCell_(originCartIdxIn) || isNumAquCell_(originCartIdxOut)) {
                // Check there are no refined aquifer cells.
                assert(level == 0);
                // Connections involving numerical aquifers are always NNCs
                // for the purpose of file output.  This holds even for
                // connections between cells like (I,J,K) and (I+1,J,K)
                // which are nominally neighbours in the Cartesian grid.
                continue;
            }

            const auto minLevelCartIdx = std::min(levelCartIdxIn, levelCartIdxOut);
            const auto maxLevelCartIdx = std::max(levelCartIdxIn, levelCartIdxOut);

            const auto& levelCartDims = computeLevelCartDims(level);

            // Re-ordering in case of non-empty mapping between equilGrid to grid
            if (map) {
                c1 = map(c1); // equilGridToGrid map
                c2 = map(c2);
            }

            if (maxLevelCartIdx - minLevelCartIdx == 1 && levelCartDims[0] > 1 ) {
                outputTrans_->at(level).at("TRANX").template data<double>()[minLevelCartIdx] =
                    gatheredOrGlobalTrans_(std::array{level, minLevelCartIdx, maxLevelCartIdx}, c1, c2);
                continue; // skip other if clauses as they are false, last one needs some computation
            }

            if (maxLevelCartIdx - minLevelCartIdx == levelCartDims[0] && levelCartDims[1] > 1) {
                outputTrans_->at(level).at("TRANY").template data<double>()[minLevelCartIdx] =
                    gatheredOrGlobalTrans_(std::array{level, minLevelCartIdx, maxLevelCartIdx}, c1, c2);
                continue; // skipt next if clause as it needs some computation
            }

            if ( maxLevelCartIdx - minLevelCartIdx == levelCartDims[0]*levelCartDims[1] ||
                 directVerticalNeighbors(levelCartDims,
                                         levelCartToLevelCompressed[level],
                                         minLevelCartIdx,
                                         maxLevelCartIdx)) {
                outputTrans_->at(level).at("TRANZ").template data<double>()[minLevelCartIdx] =
                    gatheredOrGlobalTrans_(std::array{level, minLevelCartIdx, maxLevelCartIdx}, c1, c2);
            }
        }
    }

    // A refined cell's host keeps its unrefined connections in the global section.
    const auto& dims0 = computeLevelCartDims(0);
    for (const auto& [cells, t] : globalTrans().hostLevelTransmissibilities()) {
        const auto [cell1, cell2] = cells;
        const char* name = nullptr;
        if ((cell2 - cell1 == 1) && (dims0[0] > 1)) {
            name = "TRANX";
        }
        else if ((cell2 - cell1 == dims0[0]) && (dims0[1] > 1)) {
            name = "TRANY";
        }
        else if ((cell2 - cell1 == dims0[0]*dims0[1]) ||
                 directVerticalNeighbors(dims0, levelCartToLevelCompressed[0], cell1, cell2)) {
            name = "TRANZ";
        }
        if (name != nullptr) {
            outputTrans_->at(0).at(name).template data<double>()[cell1] = t;
        }
    }
}

template<class Grid, class EquilGrid, class GridView, class ElementMapper, class Scalar>
bool
EclGenericWriter<Grid,EquilGrid,GridView,ElementMapper,Scalar>::
isCartesianNeighbour_(const std::array<int,3>& levelCartDims,
                      const std::size_t levelCartIdx1,
                      const std::size_t levelCartIdx2) const
{
    const int diff = levelCartIdx2 - levelCartIdx1;

    return (diff == 1)
        || (diff == levelCartDims[0])
        || (diff == (levelCartDims[0] * levelCartDims[1]));
}

template<class Grid, class EquilGrid, class GridView, class ElementMapper, class Scalar>
bool
EclGenericWriter<Grid,EquilGrid,GridView,ElementMapper,Scalar>::
isDirectNeighbours_(const std::unordered_map<int,int>& levelCartesianToActive,
                    const std::array<int,3>& levelCartDims,
                    const std::size_t levelCartIdx1,
                    const std::size_t levelCartIdx2) const
{
    return isCartesianNeighbour_(levelCartDims, levelCartIdx1, levelCartIdx2)
        || directVerticalNeighbors(levelCartDims, levelCartesianToActive, levelCartIdx1, levelCartIdx2);
}

template<class Grid, class EquilGrid, class GridView, class ElementMapper, class Scalar>
auto
EclGenericWriter<Grid,EquilGrid,GridView,ElementMapper,Scalar>::
activeCell_(const std::unordered_map<int,int>& levelCartToLevelCompressed,
            const std::size_t levelCartIdx) const
{
    auto pos = levelCartToLevelCompressed.find(levelCartIdx);
    return (pos == levelCartToLevelCompressed.end()) ? -1 : pos->second;
}

template<class Grid, class EquilGrid, class GridView, class ElementMapper, class Scalar>
void
EclGenericWriter<Grid,EquilGrid,GridView,ElementMapper,Scalar>::
allocateAllNncs_(int maxLevel) const
{
    this->outputNnc_.resize(maxLevel+1); // level 0,1,..., maxLevel

    if (maxLevel) {
        // NNCs between main (level zero) grid and LGRs: level 1, ...., maxLevel.
        // Example: grid with maxLevel == 3, outputNncGlobalLocal_.size() is maxLevel = 3
        //          outputNncGlobalLocal_[0] -> NNCs between level 0 and level 1
        //          outputNncGlobalLocal_[1] -> NNCs between level 0 and level 2
        //          outputAmalgamatedNnc_[2] -> NNCs between level 0 and level 3
        this->outputNncGlobalLocal_.resize(maxLevel);

        // NNCs between different refined level grids: (level1, level2)
        // with 0 < level1 < level2 <= maxLevel
        // Example: grid with maxLevel == 3, outputAmalgamatedNnc_.size() is maxLevel-1 = 2
        //          outputAmalgamatedNnc_[0][0] -> NNCs between level 1 and level 2
        //          outputAmalgamatedNnc_[0][1] -> NNCs between level 1 and level 3
        //          outputAmalgamatedNnc_[1][2] -> NNCs between level 2 and level 3
        this->outputAmalgamatedNnc_.resize(maxLevel-1);
        for (int i = 0; i < maxLevel-1; ++i) {
            this->outputAmalgamatedNnc_[i].resize(maxLevel-1-i);
        }
    }
}

template<class Grid, class EquilGrid, class GridView, class ElementMapper, class Scalar>
template<typename LevelIndicesFunction, typename OriginIndicesFunction>
std::vector<std::vector<NNCdata>>
EclGenericWriter<Grid,EquilGrid,GridView,ElementMapper,Scalar>::
exportNncStructure_(const std::vector<std::unordered_map<int,int>>& levelCartToLevelCompressed,
                    const std::function<unsigned int(unsigned int)>& map,
                    const LevelIndicesFunction& computeLevelIndices,
                    const std::function<int(int, int)>& computeLevelCartIdx,
                    const std::function<std::array<int,3>(int)>& computeLevelCartDims,
                    const OriginIndicesFunction& computeOriginIndices) const
{
    const auto& nncData = this->eclState_.getInputNNC().input();
    const auto& nncEdit = this->eclState_.getInputNNC().edit();
    const auto& nncEditr = this->eclState_.getInputNNC().editr();
    const auto& unitSystem = this->eclState_.getDeckUnitSystem();
    const auto& transMult = this->eclState_.getTransMult();

    // Cartesian index mapper for the serial I/O grid
    const auto& equilCartMapper = *collectCartMapper_;

    const auto& level0CartDims = equilCartMapper.cartesianDimensions();

    int maxLevel = this->collectGrid_->maxLevel();
    allocateAllNncs_(maxLevel);

    using GlobalGridView = typename EquilGrid::LeafGridView;
    using GlobElementMapper = Dune::MultipleCodimMultipleGeomTypeMapper<GlobalGridView>;
    const GlobalGridView& globalGridView = this->collectGrid_->leafGridView();
    const GlobElementMapper globalElemMapper { globalGridView, Dune::mcmgElementLayout() };

    for (const auto& elem : elements(globalGridView)) {
        for (const auto& is : intersections(globalGridView, elem)) {
            if (!is.neighbor())
                continue; // intersection is on the domain boundary

            // Not 'const' because remapped if 'map' is non-null.
            unsigned c1 = globalElemMapper.index(is.inside());
            unsigned c2 = globalElemMapper.index(is.outside());

            if (c1 > c2)
                continue; // we only need to handle each connection once, thank you.

            if ( is.inside().level() != is.outside().level() ) { // TRANGL and TRANLL
                // For CpGrid with LGRs, level*Idx and c* do not coincide.
                const auto& [levelInIdx, levelOutIdx] = computeLevelIndices(is, c1, c2);

                const int levelIn = is.inside().level();
                const int levelOut = is.outside().level();

                auto levelCartIdxIn = computeLevelCartIdx(levelInIdx, levelIn);
                auto levelCartIdxOut = computeLevelCartIdx(levelOutIdx, levelOut);

                // To store correctly and only once the corresponding NNC
                std::pair<int,int> smallerPair = {levelIn, levelCartIdxIn},
                    largerPair = {levelOut, levelCartIdxOut};
                if (smallerPair.first > largerPair.first) {
                    std::swap(smallerPair, largerPair);
                }

                const auto& [smallerLevel, smallerLevelCartIdx] = smallerPair;
                const auto& [largerLevel, largerLevelCartIdx] = largerPair;

                auto t = this->gatheredOrGlobalTrans_(std::array{smallerLevel, smallerLevelCartIdx,
                                                                 largerLevel, largerLevelCartIdx},
                                                      c1, c2);

                // ECLIPSE ignores NNCs with zero transmissibility
                // (different threshold than for NNC with corresponding
                // EDITNNC above).  In addition we do set small
                // transmissibilities to zero when setting up the simulator.
                // These will be ignored here, too.
                const auto tt = unitSystem
                    .from_si(UnitSystem::measure::transmissibility, t);

                if (std::isnormal(tt) && (tt > 1.0e-12)) {
                    // Store always FIRST the level Cartesian index of the cell belonging to the smaller level grid involved.
                    if (smallerLevel == 0) { // NNC between main (level zero) grid and a refined level/local grid
                        this->outputNncGlobalLocal_[largerLevel-1].emplace_back(smallerLevelCartIdx, largerLevelCartIdx, t);
                    }
                    else { // NNC between different refined level/local grids -> amlgamated NNC
                        assert(smallerLevel >= 1);
                        this->outputAmalgamatedNnc_[smallerLevel-1][largerLevel-smallerLevel-1].emplace_back(smallerLevelCartIdx, largerLevelCartIdx, t);
                    }
                }
            }
            else {
                // the cells sharing the intersection belong to the same level
                assert(is.inside().level() == is.outside().level());
                const int level = is.inside().level();

                // For CpGrid with LGRs, the origin cell index refers to the coarsest
                // ancestor cell when the cell is refined. For cells not involved in
                // any refinement, it corresponds to the geometrically equivalent
                // cell in the level-zero grid.
                const auto [originInIdx, originOutIdx] = computeOriginIndices(is, c1, c2);

                const std::size_t originCartIdxIn = computeLevelCartIdx(originInIdx, /* level = */ 0);
                const std::size_t originCartIdxOut = computeLevelCartIdx(originOutIdx, /* level = */ 0);

                // For CpGrid with LGRs, level*Idx and c* do not coincide.
                const auto& [levelInIdx, levelOutIdx] = computeLevelIndices(is, c1, c2);

                auto levelCartIdxIn = computeLevelCartIdx(levelInIdx, level);
                auto levelCartIdxOut = computeLevelCartIdx(levelOutIdx, level);

                if ( levelCartIdxOut < levelCartIdxIn )
                    std::swap(levelCartIdxIn, levelCartIdxOut);

                // Re-ordering in case of non-empty mapping between equilGrid to grid
                if (map) {
                    c1 = map(c1); // equilGridToGrid map
                    c2 = map(c2);
                }

                const auto& levelCartDims = computeLevelCartDims(level);

                // Check there are no refined aquifer connections
                assert(!isNumAquConn_(originCartIdxIn, originCartIdxOut) || level == 0);

                if (isNumAquConn_(originCartIdxIn, originCartIdxOut) ||
                    ! isDirectNeighbours_(levelCartToLevelCompressed[level],
                                          levelCartDims,
                                          levelCartIdxIn, levelCartIdxOut)) {
                    // We need to check whether an NNC for this face was also
                    // specified via the NNC keyword in the deck.
                    // (levelCartIdxIn/Out are already swapped into min/max order above.)
                    auto t = this->gatheredOrGlobalTrans_(std::array{level, levelCartIdxIn, levelCartIdxOut},
                                                          c1, c2);

                    if (level == 0) {
                        auto candidate = std::lower_bound(nncData.begin(), nncData.end(),
                                                          NNCdata { originCartIdxIn, originCartIdxOut, 0.0 });
                        const auto transMlt = transMult.getRegionMultiplierNNC(originCartIdxIn, originCartIdxOut);
                        bool foundNncEditr = false;

                        while ((candidate != nncData.end()) &&
                               (candidate->cell1 == originCartIdxIn) &&
                               (candidate->cell2 == originCartIdxOut))
                        {
                            auto trans = candidate->trans;
                            trans *= transMlt;
                            if (! nncEditr.empty()) {
                                auto it = std::lower_bound(nncEditr.begin(), nncEditr.end(),
                                                           NNCdata { originCartIdxIn, originCartIdxOut, 0.0 });
                                foundNncEditr = it != nncEditr.end() && it->cell1 == originCartIdxIn && it->cell2 == originCartIdxOut;
                            }
                            if (foundNncEditr) {
                                // Only write one value for EDITNNCR, then skip it here and add it on the second loop below
                                break;
                            }
                            if (! nncEdit.empty()) {
                                auto it = std::lower_bound(nncEdit.begin(), nncEdit.end(),
                                                           NNCdata { originCartIdxIn, originCartIdxOut, 0.0 });
                                if (it != nncEdit.end() && it->cell1 == originCartIdxIn && it->cell2 == originCartIdxOut) {
                                    trans *= it->trans;
                                }
                            }
                            t -= trans;
                            ++candidate;
                        }
                        if (foundNncEditr) {
                            // Only write one value for EDITNNCR, then skip it here and add it on the second loop below
                            continue;
                        }
                    }

                    // ECLIPSE ignores NNCs with zero transmissibility
                    // (different threshold than for NNC with corresponding
                    // EDITNNC above).  In addition we do set small
                    // transmissibilities to zero when setting up the simulator.
                    // These will be ignored here, too.
                    const auto tt = unitSystem
                        .from_si(UnitSystem::measure::transmissibility, t);

                    if (std::isnormal(tt) && (tt > 1.0e-12)) {
                        this->outputNnc_[level].emplace_back(levelCartIdxIn, levelCartIdxOut, t);
                    }
                }
            }
        }
    }

    // Do not include the generated NNCs transsmisibilities in the input NNCs
    // The hosts' unrefined connections that are not Cartesian neighbours.
    for (const auto& [cells, t] : this->globalTrans().hostLevelTransmissibilities()) {
        const auto [cell1, cell2] = cells;
        if (isDirectNeighbours_(levelCartToLevelCompressed[0], level0CartDims, cell1, cell2)) {
            continue;
        }
        const auto tt = unitSystem.from_si(UnitSystem::measure::transmissibility, t);
        if (std::isnormal(tt) && (tt > 1.0e-12)) {
            this->outputNnc_[0].emplace_back(cell1, cell2, t);
        }
    }

    std::vector<NNCdata> inputedNnc{};
    const auto generatedNnc = outputNnc_[0];

    // Deck NNCs and numerical-aquifer connections name level-zero cells, but
    // globalTrans() is indexed on the leaf.  Refinement separates the two, so
    // translate; a refined coarse cell has no single leaf cell and drops out.
    std::unordered_map<int,int> level0CartToLeaf{};
    if (maxLevel > 0) {
        for (const auto& elem : elements(globalGridView)) {
            if (elem.level() == 0) {
                const auto leafIdx = globalElemMapper.index(elem);
                level0CartToLeaf.emplace(equilCartMapper.cartesianIndex(leafIdx), leafIdx);
            }
        }
    }
    const auto& cartToLeaf = (maxLevel > 0)
        ? level0CartToLeaf
        : levelCartToLevelCompressed[/* level */ 0];
    int nncInRefinedBlock = 0;

    // The NNC keyword in the deck is defined only for faces in the level-0 grid.
    // The same limitation applies to aquifer data.
    for (const auto& entry : nncData) {
        // Ignore most explicit NNCs between otherwise neighbouring cells.
        // We keep NNCs that involve cells with numerical aquifers even if
        // these might be between neighbouring cells in the Cartesian
        // grid--e.g., between cells (I,J,K) and (I+1,J,K).  All such
        // connections should be written to NNC output arrays provided the
        // transmissibility value is sufficiently large.
        //
        // The condition cell2 >= cell1 holds by construction of nncData.
        assert (entry.cell2 >= entry.cell1);

        if (! isCartesianNeighbour_(level0CartDims, entry.cell1, entry.cell2) ||
            isNumAquConn_(entry.cell1, entry.cell2))
        {
            bool foundNncEdit = false;
            auto trans = entry.trans;
            if (! nncEdit.empty()) {
                auto it = std::lower_bound(nncEdit.begin(), nncEdit.end(),
                                           NNCdata {entry.cell1, entry.cell2, 0.0 });
                if (it != nncEdit.end() && it->cell1 == entry.cell1 && it->cell2 == entry.cell2) {
                    trans *= it->trans;
                    foundNncEdit = true;
                }
            }
            if (! foundNncEdit) {
                // Pick up transmissibility value from 'globalTrans()' since
                // multiplier keywords like MULTREGT might have impacted the
                // values entered in primary sources like NNC/EDITNNC/EDITNNCR.
                const auto c1 = activeCell_(cartToLeaf, entry.cell1);
                const auto c2 = activeCell_(cartToLeaf, entry.cell2);

                if ((c1 < 0) || (c2 < 0)) {
                    // Connection between inactive cells?  Unexpected at this
                    // level.  Might consider 'throw'ing if this happens...
                    // Under refinement it also means a refined cell, which the
                    // deck cannot name -- counted and reported below.
                    if ((maxLevel > 0) &&
                        (activeCell_(levelCartToLevelCompressed[/* level */ 0], entry.cell1) >= 0) &&
                        (activeCell_(levelCartToLevelCompressed[/* level */ 0], entry.cell2) >= 0))
                    {
                        ++nncInRefinedBlock;
                    }
                    continue;
                }

                // A deck NNC names a pair of level-zero cells. When one of them is
                // refined by an LGR the pair is no longer a leaf connection, so the
                // gathered records (correctly) hold no value for it -- keep the
                // deck-specified transmissibility in that case.
                const auto key = std::array{0, static_cast<int>(entry.cell1),
                                            static_cast<int>(entry.cell2)};
                const double* gathered = this->findGatheredTrans_(key);
                if (!this->gatheredLgrTrans_.has_value() || gathered != nullptr) {
                    // The refined output grid of a parallel run does not carry the
                    // numerical aquifer's connections; keep the input value there.
                    const auto found = (gathered != nullptr)
                        ? std::optional<double>{*gathered}
                        : std::optional<double>{this->globalTrans().findTransmissibility(c1, c2)};
                    if (!found.has_value() && (maxLevel == 0)) {
                        trans = this->globalTrans().transmissibility(c1, c2);   // reports the miss
                    }
                    trans = found.value_or(trans);

                    if (found.has_value() && ! generatedNnc.empty()) {
                        for (const auto& generated : generatedNnc) {
                            if (entry.cell1 == generated.cell1 && entry.cell2 == generated.cell2) {
                                trans -= generated.trans;
                                break;
                            }
                        }
                    }
                }
            }
            const auto tt = unitSystem
                .from_si(UnitSystem::measure::transmissibility, trans);

            // ECLIPSE ignores NNCs (with EDITNNC/EDITNNCR applied) with
            // small transmissibility values.  Seems like the threshold is
            // 1.0e-6 in output units.
            if (std::isnormal(tt) && ! (tt < 1.0e-6)) {
                inputedNnc.emplace_back(entry.cell1, entry.cell2, trans);
            }
        }
    }
    if (nncInRefinedBlock > 0) {
        OpmLog::warning(fmt::format("{} explicit NNC/aquifer connection(s) reach a cell inside "
                                    "a refined block and are omitted from the NNC output arrays.",
                                    nncInRefinedBlock));
    }

    // Write first the inputed NNCs and after the internally computed NNCs
    this->outputNnc_[0].insert(this->outputNnc_[0].begin(), inputedNnc.begin(), inputedNnc.end());
    return this->outputNnc_;
}

template<class Grid, class EquilGrid, class GridView, class ElementMapper, class Scalar>
void EclGenericWriter<Grid,EquilGrid,GridView,ElementMapper,Scalar>::
doWriteOutput(const int                          reportStepNum,
              const std::optional<int>           timeStepNum,
              const bool                         isSubStep,
              const bool                         isForcedFinalOutput,
              data::Solution&&                   localCellData,
              data::Wells&&                      localWellData,
              data::GroupAndNetworkValues&&      localGroupAndNetworkData,
              data::Aquifers&&                   localAquiferData,
              WellTestState&&                    localWTestState,
              const Action::State&               actionState,
              const UDQState&                    udqState,
              const SummaryState&                summaryState,
              const std::vector<Scalar>&         thresholdPressure,
              Scalar                             curTime,
              Scalar                             nextStepSize,
              bool                               doublePrecision,
              bool                               isFlowsn,
              std::array<FlowsData<double>, 3>&& flowsn,
              bool                               isFloresn,
              std::array<FlowsData<double>, 3>&& floresn)
{
    const auto isParallel = this->collectOnIORank_.isParallel();
    const bool needsReordering = this->collectOnIORank_.doesNeedReordering();

    RestartValue restartValue {
        (isParallel || needsReordering)
        ? this->collectOnIORank_.globalCellData()
        : std::move(localCellData),

        isParallel ? this->collectOnIORank_.globalWellData()
                   : std::move(localWellData),

        isParallel ? this->collectOnIORank_.globalGroupAndNetworkData()
                   : std::move(localGroupAndNetworkData),

        isParallel ? this->collectOnIORank_.globalAquiferData()
                   : std::move(localAquiferData)
    };

    if (eclState_.getSimulationConfig().useThresholdPressure()) {
        restartValue.addExtra("THRESHPR", UnitSystem::measure::pressure,
                              thresholdPressure);
    }

    // Add suggested next timestep to extra data.
    if (! isSubStep) {
        restartValue.addExtra("OPMEXTRA", std::vector<double>(1, nextStepSize));
    }

    // Add nnc flows and flores.
    if (isFlowsn) {
        const auto flowsn_global = isParallel ? this->collectOnIORank_.globalFlowsn() : std::move(flowsn);
        for (const auto& flows : flowsn_global) {
            if (flows.name.empty())
                continue;
            if (flows.name == "FLOGASN+") {
                restartValue.addExtra(flows.name, UnitSystem::measure::gas_surface_rate, flows.values);
            } else {
                restartValue.addExtra(flows.name, UnitSystem::measure::liquid_surface_rate, flows.values);
            }
        }
    }
    if (isFloresn) {
        const auto floresn_global = isParallel ? this->collectOnIORank_.globalFloresn() : std::move(floresn);
        for (const auto& flores : floresn_global) {
            if (flores.name.empty()) {
                continue;
            }
            restartValue.addExtra(flores.name, UnitSystem::measure::rate, flores.values);
        }
    }

    std::vector<Opm::RestartValue> restartValues{};
    const bool haveLgrCellOutput = !needsReordering
        && (this->eclState_.getLgrs().size() > 0)
        // The leaf solution is split onto the per-level grids using a grid that
        // carries the level hierarchy. That is the simulation grid in serial /
        // the rank-interior path, but with refine-before-redistribute the
        // distributed simulation grid is the flat leaf (maxLevel() == 0); there
        // the I/O-rank reference grid (collectGrid_) holds the hierarchy and the
        // parallel branch below splits on it instead. So accept either.
        && ( (this->grid_.maxLevel() > 0)
             || (this->collectGrid_ != nullptr && this->collectGrid_->maxLevel() > 0) );
    // Split the leaf solution onto the per-level grids.  Level cells that appear
    // on the leaf grid view get the data::Solution values from there; other
    // cells (parent cells that vanished due to refinement) get rubbish values
    // for now.  Only data::Solution is restricted to the level grids; well,
    // group/network and aquifer data are not modified here.
    if ( !isParallel && haveLgrCellOutput ) {
        Opm::Lgr::extractRestartValueLevelGrids<Grid>(this->grid_, restartValue, restartValues);
    }
    else if ( isParallel && haveLgrCellOutput
              && this->collectOnIORank_.isIORank()
              && (this->collectGrid_ != nullptr)
              && (this->collectGrid_->maxLevel() > 0) ) {
        // In parallel the leaf solution has been gathered onto the I/O rank's
        // refined output grid (collectGrid_); split it there.
        Opm::Lgr::extractRestartValueLevelGrids<EquilGrid>(*this->collectGrid_, restartValue, restartValues);
    }
    else {
        const auto& deckLgrs = this->eclState_.getLgrs();
        restartValues.reserve(1 + deckLgrs.size());
        restartValues.push_back(std::move(restartValue)); // no LGRs-> only one restart value

        // The deck declares LGRs but this grid refines none of them (LGROFF on
        // every one). They are still in the EGRID/INIT, so write a solution
        // section for each, father-replicated from the coarse solution, and
        // keep every report step covering the same cells.
        //
        // Restricted to a genuinely unrefined grid: the other way into this
        // branch is the reordering path, where refined data does exist and has
        // simply not been split onto the levels -- replicating there would
        // overwrite real values with coarse ones.
        const bool noRefinementAtAll = (this->grid_.maxLevel() == 0)
            && ((this->collectGrid_ == nullptr) || (this->collectGrid_->maxLevel() == 0));

        for (std::size_t lgr = 0; noRefinementAtAll && lgr < deckLgrs.size(); ++lgr) {
            const auto& name = deckLgrs.getLgr(lgr).NAME();
            restartValues.emplace_back(fatherReplicatedSolution(this->eclState_.getInputGrid(),
                                                                name,
                                                                restartValues.front().solution),
                                       restartValues.front().wells,
                                       restartValues.front().grp_nwrk,
                                       restartValues.front().aquifer,
                                       static_cast<int>(lgr) + 1);
        }
    }

    // make sure that the previous I/O request has been completed
    // and the number of incomplete tasklets does not increase between
    // time steps
    this->taskletRunner_->barrier();

    // check if there might have been a failure in the TaskletRunner
    if (this->taskletRunner_->failure()) {
        throw std::runtime_error("Failure in the TaskletRunner while writing output.");
    }

    // create a tasklet to write the data for the current time step to disk
    auto eclWriteTasklet = std::make_shared<EclWriteTasklet>(
        actionState,
        isParallel ? this->collectOnIORank_.globalWellTestState() : std::move(localWTestState),
        summaryState, udqState, *this->eclIO_,
        reportStepNum, timeStepNum, isSubStep, curTime, std::move(restartValues), doublePrecision,
        isForcedFinalOutput);

    // finally, start a new output writing job
    this->taskletRunner_->dispatch(std::move(eclWriteTasklet));
}

template<class Grid, class EquilGrid, class GridView, class ElementMapper, class Scalar>
void EclGenericWriter<Grid,EquilGrid,GridView,ElementMapper,Scalar>::
evalSummary(const int                                            reportStepNum,
            const Scalar                                         curTime,
            const data::Wells&                                   localWellData,
            const data::WellBlockAveragePressures&               localWBPData,
            const data::GroupAndNetworkValues&                   localGroupAndNetworkData,
            const std::map<int,data::AquiferData>&               localAquiferData,
            const std::map<std::pair<std::string, int>, double>& blockData,
            const std::map<std::tuple<std::string, int, int>, double>& lgrBlockData,
            const std::map<std::string, double>&                 miscSummaryData,
            const std::map<std::string, std::vector<double>>&    regionData,
            const data::RegionVariableMapping&                   regVarMap,
            const RegionVariableCollection&                      regVars,
            const Inplace&                                       inplace,
            const Inplace*                                       initialInPlace,
            const InterRegFlowMap&                               interRegFlows,
            SummaryState&                                        summaryState,
            UDQState&                                            udqState,
            const data::ReservoirCouplingGroupRates*             rcGroupRates)
{
    if (collectOnIORank_.isIORank()) {
        const auto& wellData = this->collectOnIORank_.isParallel()
            ? this->collectOnIORank_.globalWellData()
            : localWellData;

        const auto& wbpData = this->collectOnIORank_.isParallel()
            ? this->collectOnIORank_.globalWBPData()
            : localWBPData;

        const auto& groupAndNetworkData = this->collectOnIORank_.isParallel()
            ? this->collectOnIORank_.globalGroupAndNetworkData()
            : localGroupAndNetworkData;

        const auto& aquiferData = this->collectOnIORank_.isParallel()
            ? this->collectOnIORank_.globalAquiferData()
            : localAquiferData;

        const auto interreg_flows = getInterRegFlowsAsMap(interRegFlows);

        const auto values = out::Summary::DynamicSimulatorState {
            .well_solution           = &wellData,
            .wbp                     = &wbpData,
            .group_and_nwrk_solution = &groupAndNetworkData,
            .single_values           = &miscSummaryData,
            .region_values           = &regionData,
            .reg_var_map             = &regVarMap,
            .reg_var_coll            = &regVars,
            .block_values            = &blockData,
            .aquifer_values          = &aquiferData,
            .interreg_flows          = &interreg_flows,
            .rc_group_rates          = rcGroupRates,
            .inplace                 = {
                .current = &inplace,
                .initial = initialInPlace
            },
            .lgr_block_values        = &lgrBlockData
        };

        this->eclIO_->summary()
            .eval(reportStepNum, curTime, values, summaryState);

        // Off-by-one-fun: The reportStepNum argument corresponds to the
        // report step these results will be written to, whereas the
        // argument to UDQ function evaluation corresponds to the report
        // step we are currently on.
        const auto udq_step = reportStepNum - 1;

        this->schedule_[udq_step].udq()
            .eval(udq_step,
                  this->schedule_.wellMatcher(udq_step),
                  this->schedule_[udq_step].group_order(),
                  this->schedule_.segmentMatcherFactory(udq_step),
                  [es = std::cref(this->eclState_)]() {
                      return std::make_unique<RegionSetMatcher>
                          (es.get().fipRegionStatistics());
                  },
                  summaryState, udqState);
    }

#if HAVE_MPI
    if (collectOnIORank_.isParallel()) {
        Parallel::MpiSerializer ser(grid_.comm());
        ser.append(summaryState);
    }
#endif
}

template<class Grid, class EquilGrid, class GridView, class ElementMapper, class Scalar>
template<std::size_t N>
const double*
EclGenericWriter<Grid,EquilGrid,GridView,ElementMapper,Scalar>::
findGatheredTrans_(const std::array<int,N>& key) const
{
    static_assert(N == 3 || N == 4, "unknown gathered LGR record shape");

    if (!gatheredLgrTrans_.has_value()) {
        return nullptr;
    }

    if constexpr (N == 3) {
        return gatheredLgrTrans_->sameLevel.find(key);
    } else {
        return gatheredLgrTrans_->crossLevel.find(key);
    }
}

template<class Grid, class EquilGrid, class GridView, class ElementMapper, class Scalar>
template<std::size_t N>
double
EclGenericWriter<Grid,EquilGrid,GridView,ElementMapper,Scalar>::
gatheredOrGlobalTrans_(const std::array<int,N>& key,
                       unsigned c1,
                       unsigned c2) const
{
    if (!gatheredLgrTrans_.has_value()) {
        return this->globalTrans().transmissibility(c1, c2);
    }

    if (const double* value = this->findGatheredTrans_(key); value != nullptr) {
        return *value;
    }

    std::string msg = "Gathered LGR transmissibilities: no value for connection key (";
    for (std::size_t j = 0; j < N; ++j) {
        msg += std::to_string(key[j]);
        msg += (j + 1 < N) ? ", " : ")";
    }
    throw std::logic_error { msg };
}

template<class Grid, class EquilGrid, class GridView, class ElementMapper, class Scalar>
const typename EclGenericWriter<Grid,EquilGrid,GridView,ElementMapper,Scalar>::TransmissibilityType&
EclGenericWriter<Grid,EquilGrid,GridView,ElementMapper,Scalar>::
globalTrans() const
{
    assert (globalTrans_);
    return *globalTrans_;
}

} // namespace Opm

#endif // OPM_ECL_GENERIC_WRITER_IMPL_HPP
