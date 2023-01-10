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

#include <config.h>
#include <ebos/eclgenericwriter.hh>

#include <opm/grid/CpGrid.hpp>
#include <opm/grid/cpgrid/GridHelpers.hpp>
#include <opm/grid/polyhedralgrid.hh>
#include <opm/grid/utility/cartesianToCompressed.hpp>
#if HAVE_DUNE_ALUGRID
#include "eclalugridvanguard.hh"
#include <dune/alugrid/grid.hh>
#include <dune/alugrid/3d/gridview.hh>
#endif // HAVE_DUNE_ALUGRID

#include <opm/output/eclipse/EclipseIO.hpp>
#include <opm/output/eclipse/RestartValue.hpp>
#include <opm/output/eclipse/Summary.hpp>

#include <opm/input/eclipse/EclipseState/EclipseState.hpp>
#include <opm/input/eclipse/EclipseState/SummaryConfig/SummaryConfig.hpp>
#include <opm/input/eclipse/Schedule/Action/State.hpp>
#include <opm/input/eclipse/Schedule/Schedule.hpp>
#include <opm/input/eclipse/Schedule/SummaryState.hpp>
#include <opm/input/eclipse/Schedule/UDQ/UDQConfig.hpp>
#include <opm/input/eclipse/Schedule/UDQ/UDQState.hpp>
#include <opm/input/eclipse/Schedule/Well/WellMatcher.hpp>
#include <opm/input/eclipse/Units/UnitSystem.hpp>

#include <dune/grid/common/mcmgmapper.hh>

#if HAVE_MPI
#include <ebos/eclmpiserializer.hh>
#endif

#if HAVE_DUNE_FEM
#include <dune/fem/gridpart/adaptiveleafgridpart.hh>
#include <dune/fem/gridpart/common/gridpart2gridview.hh>
#include <ebos/femcpgridcompat.hh>
#endif // HAVE_DUNE_FEM

#if HAVE_MPI
#include <mpi.h>
#endif

#include <array>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace {

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
getInterRegFlowsAsMap(const Opm::EclInterRegFlowMap& map)
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
    bool isSubStep_;
    double secondsElapsed_;
    Opm::RestartValue restartValue_;
    bool writeDoublePrecision_;

    explicit EclWriteTasklet(const Opm::Action::State& actionState,
                             const Opm::WellTestState& wtestState,
                             const Opm::SummaryState& summaryState,
                             const Opm::UDQState& udqState,
                             Opm::EclipseIO& eclIO,
                             int reportStepNum,
                             bool isSubStep,
                             double secondsElapsed,
                             Opm::RestartValue restartValue,
                             bool writeDoublePrecision)
        : actionState_(actionState)
        , wtestState_(wtestState)
        , summaryState_(summaryState)
        , udqState_(udqState)
        , eclIO_(eclIO)
        , reportStepNum_(reportStepNum)
        , isSubStep_(isSubStep)
        , secondsElapsed_(secondsElapsed)
        , restartValue_(restartValue)
        , writeDoublePrecision_(writeDoublePrecision)
    { }

    // callback to eclIO serial writeTimeStep method
    void run()
    {
        eclIO_.writeTimeStep(actionState_,
                             wtestState_,
                             summaryState_,
                             udqState_,
                             reportStepNum_,
                             isSubStep_,
                             secondsElapsed_,
                             restartValue_,
                             writeDoublePrecision_);
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
                 bool enableEsmry )
    : collectToIORank_(grid,
                       equilGrid,
                       gridView,
                       cartMapper,
                       equilCartMapper,
                       summaryConfig.fip_regions_interreg_flow())
    , grid_(grid)
    , gridView_(gridView)
    , schedule_(schedule)
    , eclState_(eclState)
    , summaryConfig_(summaryConfig)
    , cartMapper_(cartMapper)
    , equilCartMapper_(equilCartMapper)
    , equilGrid_(equilGrid)
{
    if (collectToIORank_.isIORank()) {
        eclIO_.reset(new EclipseIO(eclState_,
                                   UgGridHelpers::createEclipseGrid(*equilGrid, eclState_.getInputGrid()),
                                   schedule_,
                                   summaryConfig_, "", enableEsmry));

        const auto& wbp_calculators = eclIO_->summary().wbp_calculators( schedule.size() - 1 );
        wbp_index_list_ = wbp_calculators.index_list();
    }
    if (collectToIORank_.isParallel()) {
        const auto& comm = grid_.comm();
        unsigned long size = wbp_index_list_.size();
        comm.broadcast(&size, 1, collectToIORank_.ioRank);
        if (!collectToIORank_.isIORank())
            wbp_index_list_.resize( size );
        comm.broadcast(wbp_index_list_.data(), size, collectToIORank_.ioRank);
    }
    // create output thread if enabled and rank is I/O rank
    // async output is enabled by default if pthread are enabled
    int numWorkerThreads = 0;
    if (enableAsyncOutput && collectToIORank_.isIORank())
        numWorkerThreads = 1;
    taskletRunner_.reset(new TaskletRunner(numWorkerThreads));
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
writeInit(const std::function<unsigned int(unsigned int)>& map)
{
    if (collectToIORank_.isIORank()) {
        std::map<std::string, std::vector<int> > integerVectors;
        if (collectToIORank_.isParallel())
            integerVectors.emplace("MPI_RANK", collectToIORank_.globalRanks());
        auto cartMap = cartesianToCompressed(equilGrid_->size(0), UgGridHelpers::globalCell(*equilGrid_));
        eclIO_->writeInitial(computeTrans_(cartMap, map), integerVectors, exportNncStructure_(cartMap, map));
    }
}

template<class Grid, class EquilGrid, class GridView, class ElementMapper, class Scalar>
data::Solution EclGenericWriter<Grid,EquilGrid,GridView,ElementMapper,Scalar>::
computeTrans_(const std::unordered_map<int,int>& cartesianToActive, const std::function<unsigned int(unsigned int)>& map) const
{
    const auto& cartMapper = *equilCartMapper_;
    const auto& cartDims = cartMapper.cartesianDimensions();
    const int globalSize = cartDims[0]*cartDims[1]*cartDims[2];

    data::CellData tranx = {UnitSystem::measure::transmissibility, std::vector<double>(globalSize), data::TargetType::INIT};
    data::CellData trany = {UnitSystem::measure::transmissibility, std::vector<double>(globalSize), data::TargetType::INIT};
    data::CellData tranz = {UnitSystem::measure::transmissibility, std::vector<double>(globalSize), data::TargetType::INIT};

    for (size_t i = 0; i < tranx.data.size(); ++i) {
        tranx.data[0] = 0.0;
        trany.data[0] = 0.0;
        tranz.data[0] = 0.0;
    }

    using GlobalGridView = typename EquilGrid::LeafGridView;
    const GlobalGridView& globalGridView = equilGrid_->leafGridView();
    using GlobElementMapper = Dune::MultipleCodimMultipleGeomTypeMapper<GlobalGridView>;
    GlobElementMapper globalElemMapper(globalGridView, Dune::mcmgElementLayout());

    for (const auto& elem : elements(globalGridView)) {
        for (const auto& is : intersections(globalGridView, elem)) {
            if (!is.neighbor())
                continue; // intersection is on the domain boundary

            unsigned c1 = globalElemMapper.index(is.inside());
            unsigned c2 = globalElemMapper.index(is.outside());

            if (c1 > c2)
                continue; // we only need to handle each connection once, thank you.

            // Ordering of compressed and uncompressed index should be the same
            const int cartIdx1 = cartMapper.cartesianIndex( c1 );
            const int cartIdx2 = cartMapper.cartesianIndex( c2 );
            // Ordering of compressed and uncompressed index should be the same
            assert(cartIdx1 <= cartIdx2);
            int gc1 = std::min(cartIdx1, cartIdx2);
            int gc2 = std::max(cartIdx1, cartIdx2);

            // Re-ordering in case of non-empty mapping between equilGrid to grid
            if (map) {
                c1 = map(c1); // equilGridToGrid map
                c2 = map(c2);
            }

            if (gc2 - gc1 == 1 && cartDims[0] > 1 ) {
                tranx.data[gc1] = globalTrans().transmissibility(c1, c2);
                continue; // skip other if clauses as they are false, last one needs some computation
            }

            if (gc2 - gc1 == cartDims[0] && cartDims[1] > 1) {
                trany.data[gc1] = globalTrans().transmissibility(c1, c2);
                continue; // skipt next if clause as it needs some computation
            }

            if ( gc2 - gc1 == cartDims[0]*cartDims[1] ||
                 directVerticalNeighbors(cartDims, cartesianToActive, gc1, gc2))
                tranz.data[gc1] = globalTrans().transmissibility(c1, c2);
        }
    }

    return {{"TRANX", tranx},
            {"TRANY", trany},
            {"TRANZ", tranz}};
}

template<class Grid, class EquilGrid, class GridView, class ElementMapper, class Scalar>
std::vector<NNCdata> EclGenericWriter<Grid,EquilGrid,GridView,ElementMapper,Scalar>::
exportNncStructure_(const std::unordered_map<int,int>& cartesianToActive, const std::function<unsigned int(unsigned int)>& map) const
{
    std::size_t nx = eclState_.getInputGrid().getNX();
    std::size_t ny = eclState_.getInputGrid().getNY();
    auto nncData = eclState_.getInputNNC().input();
    const auto& unitSystem = eclState_.getDeckUnitSystem();
    std::vector<NNCdata> outputNnc;
    std::size_t index = 0;

    for( const auto& entry : nncData ) {
        // test whether NNC is not a neighboring connection
        // cell2>=cell1 holds due to sortNncAndApplyEditnnc
        assert( entry.cell2 >= entry.cell1 );
        auto cellDiff = entry.cell2 - entry.cell1;

        if (cellDiff != 1 && cellDiff != nx && cellDiff != nx*ny) {
            auto tt = unitSystem.from_si(UnitSystem::measure::transmissibility, entry.trans);
            // Eclipse ignores NNCs (with EDITNNC applied) that are small. Seems like the threshold is 1.0e-6
            if ( tt >= 1.0e-6 )
                outputNnc.emplace_back(entry.cell1, entry.cell2, entry.trans);
        }
        ++index;
    }

    using GlobalGridView = typename EquilGrid::LeafGridView;
    const GlobalGridView& globalGridView = equilGrid_->leafGridView();
    using GlobElementMapper = Dune::MultipleCodimMultipleGeomTypeMapper<GlobalGridView>;
    GlobElementMapper globalElemMapper(globalGridView, Dune::mcmgElementLayout());

    // Cartesian index mapper for the serial I/O grid
    const auto& equilCartMapper =  *equilCartMapper_;
    const auto& cartDims = cartMapper_.cartesianDimensions();
    for (const auto& elem : elements(globalGridView)) {
        for (const auto& is : intersections(globalGridView, elem)) {
            if (!is.neighbor())
                continue; // intersection is on the domain boundary

            unsigned c1 = globalElemMapper.index(is.inside());
            unsigned c2 = globalElemMapper.index(is.outside());

            if (c1 > c2)
                continue; // we only need to handle each connection once, thank you.

            std::size_t cc1 = equilCartMapper.cartesianIndex( c1 );
            std::size_t cc2 = equilCartMapper.cartesianIndex( c2 );

            if ( cc2 < cc1 )
                std::swap(cc1, cc2);

            auto cellDiff = cc2 - cc1;

            // Re-ordering in case of non-empty mapping between equilGrid to grid
            if (map) {
                c1 = map(c1); // equilGridToGrid map
                c2 = map(c2);
            }

            if (cellDiff != 1 &&
                cellDiff != nx &&
                cellDiff != nx*ny &&
                !directVerticalNeighbors(cartDims, cartesianToActive, cc1, cc2)) {
                // We need to check whether an NNC for this face was also specified
                // via the NNC keyword in the deck (i.e. in the first origNncSize entries.
                auto t = globalTrans().transmissibility(c1, c2);
                auto candidate = std::lower_bound(nncData.begin(), nncData.end(), NNCdata(cc1, cc2, 0.0));

                while ( candidate != nncData.end() && candidate->cell1 == cc1
                     && candidate->cell2 == cc2) {
                    t -= candidate->trans;
                    ++candidate;
                }
                // eclipse ignores NNCs with zero transmissibility (different threshold than for NNC
                // with corresponding EDITNNC above). In addition we do set small transmissibilties
                // to zero when setting up the simulator. These will be ignored here, too.
                auto tt = unitSystem.from_si(UnitSystem::measure::transmissibility, std::abs(t));
                if ( tt > 1e-12 )
                    outputNnc.push_back({cc1, cc2, t});
            }
        }
    }
    return outputNnc;
}

template<class Grid, class EquilGrid, class GridView, class ElementMapper, class Scalar>
void EclGenericWriter<Grid,EquilGrid,GridView,ElementMapper,Scalar>::
doWriteOutput(const int                     reportStepNum,
              const bool                    isSubStep,
              data::Solution&&              localCellData,
              data::Wells&&                 localWellData,
              data::GroupAndNetworkValues&& localGroupAndNetworkData,
              data::Aquifers&&              localAquiferData,
              WellTestState&&               localWTestState,
              const Action::State& actionState,
              const UDQState& udqState,
              const SummaryState& summaryState,
              const std::vector<Scalar>& thresholdPressure,
              Scalar curTime,
              Scalar nextStepSize,
              bool doublePrecision)
{
    const auto isParallel = this->collectToIORank_.isParallel();
    const bool needsReordering = this->collectToIORank_.doesNeedReordering();

    RestartValue restartValue {
        (isParallel || needsReordering) ? this->collectToIORank_.globalCellData()
                   : std::move(localCellData),

        isParallel ? this->collectToIORank_.globalWellData()
                   : std::move(localWellData),

        isParallel ? this->collectToIORank_.globalGroupAndNetworkData()
                   : std::move(localGroupAndNetworkData),

        isParallel ? this->collectToIORank_.globalAquiferData()
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

    // first, create a tasklet to write the data for the current time
    // step to disk
    auto eclWriteTasklet = std::make_shared<EclWriteTasklet>(
        actionState,
        isParallel ? this->collectToIORank_.globalWellTestState() : std::move(localWTestState),
        summaryState, udqState, *this->eclIO_,
        reportStepNum, isSubStep, curTime, std::move(restartValue), doublePrecision);

    // then, make sure that the previous I/O request has been completed
    // and the number of incomplete tasklets does not increase between
    // time steps
    this->taskletRunner_->barrier();

    // finally, start a new output writing job
    this->taskletRunner_->dispatch(std::move(eclWriteTasklet));
}

template<class Grid, class EquilGrid, class GridView, class ElementMapper, class Scalar>
void EclGenericWriter<Grid,EquilGrid,GridView,ElementMapper,Scalar>::
evalSummary(const int reportStepNum,
            const Scalar curTime,
            const std::map<std::size_t, double>& wbpData,
            const data::Wells& localWellData,
            const data::GroupAndNetworkValues& localGroupAndNetworkData,
            const std::map<int,data::AquiferData>& localAquiferData,
            const std::map<std::pair<std::string, int>, double>& blockData,
            const std::map<std::string, double>& miscSummaryData,
            const std::map<std::string, std::vector<double>>& regionData,
            const Inplace& inplace,
            const Inplace& initialInPlace,
            const EclInterRegFlowMap& interRegionFlowMap,
            SummaryState& summaryState,
            UDQState& udqState)
{
    if (collectToIORank_.isIORank()) {
        const auto& summary = eclIO_->summary();
        auto wbp_calculators = summary.wbp_calculators(reportStepNum);

        for (const auto& [global_index, pressure] : wbpData)
            wbp_calculators.add_pressure( global_index, pressure );

        const auto& wellData = this->collectToIORank_.isParallel()
            ? this->collectToIORank_.globalWellData()
            : localWellData;

        const auto& groupAndNetworkData = this->collectToIORank_.isParallel()
            ? this->collectToIORank_.globalGroupAndNetworkData()
            : localGroupAndNetworkData;

        const auto& aquiferData = this->collectToIORank_.isParallel()
            ? this->collectToIORank_.globalAquiferData()
            : localAquiferData;

        summary.eval(summaryState,
                     reportStepNum,
                     curTime,
                     wellData,
                     groupAndNetworkData,
                     miscSummaryData,
                     initialInPlace,
                     inplace,
                     wbp_calculators,
                     regionData,
                     blockData,
                     aquiferData,
                     getInterRegFlowsAsMap(interRegionFlowMap));

        // Off-by-one-fun: The reportStepNum argument corresponds to the
        // report step these results will be written to, whereas the
        // argument to UDQ function evaluation corresponds to the report
        // step we are currently on.
        auto udq_step = reportStepNum - 1;
        const auto& udq_config = schedule_.getUDQConfig(udq_step);
        udq_config.eval( udq_step, schedule_.wellMatcher(udq_step), summaryState, udqState);
    }
#if HAVE_MPI
    if (collectToIORank_.isParallel()) {
        EclMpiSerializer ser(grid_.comm());
        ser.append(summaryState);
    }
#endif
}

template<class Grid, class EquilGrid, class GridView, class ElementMapper, class Scalar>
const typename EclGenericWriter<Grid,EquilGrid,GridView,ElementMapper,Scalar>::TransmissibilityType&
EclGenericWriter<Grid,EquilGrid,GridView,ElementMapper,Scalar>::
globalTrans() const
{
    assert (globalTrans_);
    return *globalTrans_;
}

#if HAVE_DUNE_FEM
template class EclGenericWriter<Dune::CpGrid,
                                Dune::CpGrid,
                                Dune::GridView<Dune::Fem::GridPart2GridViewTraits<Dune::Fem::AdaptiveLeafGridPart<Dune::CpGrid, Dune::PartitionIteratorType(4), false>>>, Dune::MultipleCodimMultipleGeomTypeMapper<Dune::GridView<Dune::Fem::GridPart2GridViewTraits<Dune::Fem::AdaptiveLeafGridPart<Dune::CpGrid, Dune::PartitionIteratorType(4), false>>>>,
                                double>;
template class EclGenericWriter<Dune::CpGrid,
                                Dune::CpGrid,
                                Dune::Fem::GridPart2GridViewImpl<
                                    Dune::Fem::AdaptiveLeafGridPart<
                                        Dune::CpGrid,
                                        Dune::PartitionIteratorType(4),
                                        false>>,
                                Dune::MultipleCodimMultipleGeomTypeMapper<
                                    Dune::Fem::GridPart2GridViewImpl<
                                        Dune::Fem::AdaptiveLeafGridPart<
                                            Dune::CpGrid,
                                            Dune::PartitionIteratorType(4),
                                            false>>>,
                                double>;


#ifdef HAVE_DUNE_ALUGRID
#if HAVE_MPI
    using ALUGrid3CN = Dune::ALUGrid<3, 3, Dune::cube, Dune::nonconforming, Dune::ALUGridMPIComm>;
#else
    using ALUGrid3CN = Dune::ALUGrid<3, 3, Dune::cube, Dune::nonconforming, Dune::ALUGridNoComm>;
#endif //HAVE_MPI
                                                               
template class EclGenericWriter<ALUGrid3CN,
                                Dune::CpGrid,
                                Dune::GridView<Dune::Fem::GridPart2GridViewTraits<Dune::Fem::AdaptiveLeafGridPart<ALUGrid3CN, Dune::PartitionIteratorType(4), false>>>, Dune::MultipleCodimMultipleGeomTypeMapper<Dune::GridView<Dune::Fem::GridPart2GridViewTraits<Dune::Fem::AdaptiveLeafGridPart<ALUGrid3CN, Dune::PartitionIteratorType(4), false>>>>,
                                double>;
                                
template class EclGenericWriter<ALUGrid3CN,
                                Dune::CpGrid,
                                Dune::Fem::GridPart2GridViewImpl<
                                    Dune::Fem::AdaptiveLeafGridPart<
                                        ALUGrid3CN,
                                        Dune::PartitionIteratorType(4),
                                        false>>,
                                Dune::MultipleCodimMultipleGeomTypeMapper<
                                    Dune::Fem::GridPart2GridViewImpl<
                                        Dune::Fem::AdaptiveLeafGridPart<
                                            ALUGrid3CN,
                                            Dune::PartitionIteratorType(4),
                                            false>>>,
                                double>;                                
#endif // HAVE_DUNE_ALUGRID

#else // !HAVE_DUNE_FEM
template class EclGenericWriter<Dune::CpGrid,
                                Dune::CpGrid,
                                Dune::GridView<Dune::DefaultLeafGridViewTraits<Dune::CpGrid>>,
                                Dune::MultipleCodimMultipleGeomTypeMapper<Dune::GridView<Dune::DefaultLeafGridViewTraits<Dune::CpGrid>>>,
                                double>;
#if HAVE_DUNE_ALUGRID
#if HAVE_MPI
    using ALUGrid3CN = Dune::ALUGrid<3, 3, Dune::cube, Dune::nonconforming, Dune::ALUGridMPIComm>;
#else
    using ALUGrid3CN = Dune::ALUGrid<3, 3, Dune::cube, Dune::nonconforming, Dune::ALUGridNoComm>;
#endif //HAVE_MPI
template class EclGenericWriter<ALUGrid3CN,
                                Dune::CpGrid,
                                Dune::GridView<Dune::ALU3dLeafGridViewTraits<const ALUGrid3CN,Dune::PartitionIteratorType(4)>>,  Dune::MultipleCodimMultipleGeomTypeMapper<Dune::GridView<Dune::ALU3dLeafGridViewTraits<const ALUGrid3CN, Dune::PartitionIteratorType(4)>>>,
                                double>;

#endif // HAVE_DUNE_ALUGRID
#endif // !HAVE_DUNE_FEM

template class EclGenericWriter<Dune::PolyhedralGrid<3,3,double>,
                                Dune::PolyhedralGrid<3,3,double>,
                                Dune::GridView<Dune::PolyhedralGridViewTraits<3, 3, double, Dune::PartitionIteratorType(4)>>,                              Dune::MultipleCodimMultipleGeomTypeMapper<Dune::GridView<Dune::PolyhedralGridViewTraits<3,3,double,Dune::PartitionIteratorType(4)>>>,
                                double>;
                                                                  
} // namespace Opm
