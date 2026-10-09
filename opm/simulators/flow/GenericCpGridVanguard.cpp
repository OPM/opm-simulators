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
#include <filesystem>
#include <opm/simulators/flow/GenericCpGridVanguard.hpp>

#include <dune/grid/common/mcmgmapper.hh>
#include <dune/grid/common/partitionset.hh>
#include <dune/common/version.hh>

#include <opm/common/TimingMacros.hpp>
#include <opm/common/utility/ActiveGridCells.hpp>

#include <opm/grid/cpgrid/GridHelpers.hpp>
#include <opm/grid/cpgrid/LevelCartesianIndexMapper.hpp>

#include <opm/grid/cpgrid/refinement/GridStateWriter.hpp>

#include <opm/input/eclipse/EclipseState/EclipseState.hpp>
#include <opm/input/eclipse/Schedule/Schedule.hpp>
#include <opm/input/eclipse/Schedule/Well/Well.hpp>
#include <opm/input/eclipse/Schedule/Well/WellConnections.hpp>
#include <opm/input/eclipse/EclipseState/Grid/LgrCollection.hpp>

#include <opm/simulators/utils/ParallelEclipseState.hpp>
#include <opm/simulators/utils/ParallelSerialization.hpp>
#include <opm/simulators/utils/PropsDataHandle.hpp>
#include <opm/simulators/utils/SetupPartitioningParams.hpp>

#if HAVE_MPI
#include <opm/simulators/utils/MPISerializer.hpp>
#endif

#if HAVE_DUNE_FEM
#include <dune/fem/gridpart/adaptiveleafgridpart.hh>
#include <opm/simulators/flow/FemCpGridCompat.hpp>
#endif //HAVE_DUNE_FEM

#include <algorithm>
#include <cassert>
#include <fstream>
#include <iterator>
#include <map>
#include <numeric>
#include <optional>
#include <stdexcept>
#include <string>
#include <vector>

#include <fmt/format.h>

namespace Opm {

#if HAVE_MPI
namespace details {

std::vector<int>
MPIPartitionFromFile::operator()(const Dune::CpGrid& grid) const
{
    // Check if the external partition file exists
    if (!std::filesystem::exists(this->partitionFile_)) {
        OPM_THROW(std::runtime_error,
                  fmt::format("External partition file '{}' does not exist.", this->partitionFile_.string()));
    }

    std::ifstream pfile { this->partitionFile_ };

    auto partition = std::vector<int> {
        std::istream_iterator<int> { pfile },
        std::istream_iterator<int> {}
    };

    const auto nc =
        static_cast<std::vector<int>::size_type>(grid.size(0));

    if (partition.size() == nc) {
        // Input is one process ID for each active cell
        return partition;
    }

    if (partition.size() == 3 * nc) {
        // Partition file is of the form
        //
        //   Process_ID   Cartesian_Idx   NLDD_Domain
        //
        // with one row for each active cell.  Select first column.
        auto g2l = std::unordered_map<int, int>{};
        auto locCell = 0;
        for (const auto& globCell : grid.globalCell()) {
            g2l.insert_or_assign(globCell, locCell++);
        }

        auto filtered_partition = std::vector<int>(nc);
        for (auto c = 0*nc; c < nc; ++c) {
            auto pos = g2l.find(partition[3*c + 1]);
            if (pos != g2l.end()) {
                filtered_partition[pos->second] = partition[3*c + 0];
            }
        }

        return filtered_partition;
    }

    throw std::invalid_argument {
        fmt::format("Partition file '{}' with {} values does "
                    "not match CpGrid instance with {} cells",
                    this->partitionFile_.generic_string(),
                    partition.size(), nc)
    };
}

}
#endif

std::optional<std::function<std::vector<int> (const Dune::CpGrid&)>> externalLoadBalancer;

template<class ElementMapper, class GridView, class Scalar>
GenericCpGridVanguard<ElementMapper,GridView,Scalar>::GenericCpGridVanguard()
{
    this->mpiRank = 0;

#if HAVE_MPI
    this->mpiRank = FlowGenericVanguard::comm().rank();
#endif  // HAVE_MPI
}

template<class ElementMapper, class GridView, class Scalar>
void GenericCpGridVanguard<ElementMapper,GridView,Scalar>::releaseEquilGrid()
{
    this->equilGrid_.reset();
    this->equilCartesianIndexMapper_.reset();
}

#if HAVE_MPI
template<class ElementMapper, class GridView, class Scalar>
void GenericCpGridVanguard<ElementMapper, GridView, Scalar>::
doLoadBalance_(const Dune::EdgeWeightMethod             edgeWeightsMethod,
               const bool                               ownersFirst,
               const bool                               addCorners,
               const int                                numOverlap,
               const Dune::PartitionMethod              partitionMethod,
               const bool                               serialPartitioning,
               const bool                               enableDistributedWells,
               const bool                               allowSplittingInactiveWells,
               const double                             imbalanceTol,
               const GridView&                          gridView,
               const Schedule&                          schedule,
               EclipseState&                            eclState1,
               FlowGenericVanguard::ParallelWellStruct& parallelWells,
               const int                                numJacobiBlocks,
               const bool                               enableEclOutput)
{
    if (((partitionMethod == Dune::PartitionMethod::zoltan) ||
         (partitionMethod == Dune::PartitionMethod::zoltanGoG)) &&
        !this->zoltanParams().empty())
    {
        this->grid_->setPartitioningParams
            (setupZoltanParams(this->zoltanParams(),
                               this->zoltanPhgEdgeSizeThreshold()));
    }
    else if ((partitionMethod == Dune::PartitionMethod::metis) &&
        !this->metisParams().empty())
    {
        this->grid_->setPartitioningParams
            (setupMetisParams(this->metisParams()));
    }

    const auto mpiSize = this->grid_->comm().size();

    const auto partitionJacobiBlocks =
        (numJacobiBlocks > 1) && (mpiSize == 1);

    if ((mpiSize > 1) || (numJacobiBlocks > 1)) {
        // CpGrid's loadBalance() method uses transmissibilities as edge
        // weights.  This is arguably a layering violation and extracting
        // the per-face transmissibilities as a linear array is relatively
        // expensive.  We therefore extract transmissibility values only if
        // the values are actually needed.
        auto loadBalancerSet = static_cast<int>(externalLoadBalancer.has_value());
        this->grid_->comm().broadcast(&loadBalancerSet, 1, 0);

        if ((this->grid_->size(0) > 0) &&
            (enableEclOutput ||
             (loadBalancerSet == 0) ||
             partitionJacobiBlocks))
        {
            this->allocTrans();
        }

        std::vector<double> faceTrans;
        {
            OPM_TIMEBLOCK(extractTrans);
            if (loadBalancerSet == 0 || partitionJacobiBlocks) {
                faceTrans = this->extractFaceTrans(gridView);
            }
        }

        // Skipping inactive wells in partitioning currently does not play nice with restart..
        const bool restart = eclState1.getInitConfig().restartRequested();
        const bool split_inactive = (!restart && allowSplittingInactiveWells);
        const auto wells = ((mpiSize > 1) || partitionJacobiBlocks)
            ? split_inactive
                ? schedule.getActiveWellsAtEnd()
                : schedule.getWellsatEnd()
            : std::vector<Well>{};

        // NOTE: a mutable copy (upstream takes a const&): the LGR well
        // anchoring below inserts additional entries.
        auto possibleFutureConnections = schedule.getPossibleFutureConnections();
        // Mechanics couples cells across zero-transmissibility faces,
        // so the overlap layer must not be pruned by transmissibility
        // in mechanical runs (same reasoning as for thermal runs).
        const bool useTransToFilterOverlap =
                !(eclState1.getSimulationConfig().isThermal() ||
                    eclState1.getSimulationConfig().isTemp() ||
                    eclState1.runspec().mech());
        // Wells completed inside an LGR (COMPDATL) carry LGR-local connection
        // (i,j,k) that WellConnections cannot place in the level-zero load-
        // balance graph, so it skips them and the well is left unanchored. The
        // partitioner may then put the well on a different rank than the one
        // owning its (rank-interior) refinement box; the box's connection cells
        // are then not found on the well's rank (ParallelWellInfo "cells not
        // found", surfacing as a deadlock/failure at higher rank counts - works
        // at np<=4 only by partition coincidence). Anchor each such well to its
        // box's coarse parent cells (which the LGR partition cell groups keep on
        // one rank) so the well lands on the box's rank.
        if (mpiSize > 1) {
            const auto& lgrs = eclState1.getLgrs();
            if (lgrs.size() > 0) {
                const auto cartDims = this->grid_->logicalCartesianSize();
                for (const auto& well : wells) {
                    auto& anchors = possibleFutureConnections[well.name()];
                    for (const auto& conn : well.getConnections()) {
                        // A connection carries its own LGR (deck number); a
                        // coarse connection is already in the graph.
                        const int n = conn.get_lgr_level();
                        if (n <= 0 || static_cast<std::size_t>(n) > lgrs.size()) {
                            continue;
                        }
                        const auto ijk = lgrs.levelZeroIJK(static_cast<std::size_t>(n),
                                                           { conn.getI(), conn.getJ(), conn.getK() });
                        anchors.insert(ijk[0] + cartDims[0] * (ijk[1] + cartDims[1] * ijk[2]));
                    }
                }
            }
        }
        // Distribute the grid and switch to the distributed view.
        if (mpiSize > 1) {
            this->distributeGrid(edgeWeightsMethod, ownersFirst,
                                 addCorners, numOverlap, partitionMethod,
                                 serialPartitioning, enableDistributedWells,
                                                                 imbalanceTol, loadBalancerSet != 0,
                                                                 useTransToFilterOverlap,
                                 faceTrans, wells,
                                 possibleFutureConnections,
                                 eclState1, parallelWells);
        }

        // Add inactive wells to all ranks with connections (not solved, so OK even without distributed wells)
        if (split_inactive) {
            std::unordered_set<unsigned> cellOnRank;
            const auto& global_cells = this->grid_->globalCell();
            for (const auto cell : global_cells) cellOnRank.insert(cell);
            const auto inactive_well_names = schedule.getInactiveWellNamesAtEnd();
            std::size_t num_wells = inactive_well_names.size();
            std::vector<int> well_on_rank(num_wells, 0);
            std::size_t well_idx = 0;
            for (const auto& well_name : inactive_well_names) {
                const auto& well = schedule.getWell(well_name, schedule.size()-1);
                for (const auto& conn: well.getConnections()) {
                    if (cellOnRank.count(conn.global_index()) > 0) {
                        well_on_rank[well_idx] = 1;
                        parallelWells.emplace_back(well_name, true);
                        break;
                    }
                }
                if (well_on_rank[well_idx] == 0) {
                    parallelWells.emplace_back(well_name, false);
                }
                ++well_idx;
            }

            // Provide information message
            const auto& comm = this->grid_->comm();
            const auto nranks = comm.size();
            // // values from rank i will be at indices i*num_wells, i*num_wells + 1, ..., (i+1) * num_wells -1
            std::vector<int> well_on_rank_global(num_wells * nranks, 0);
            comm.allgather(well_on_rank.data(), static_cast<int>(num_wells), well_on_rank_global.data());
            if (comm.rank() == 0) {
                well_idx = 0;
                for (const auto& well_name : inactive_well_names) {
                    std::string msg = fmt::format("Well {} is inactive, with perforations on ranks: ", well_name);
                    for (int i = 0; i < nranks; ++i) {
                        if (well_on_rank_global[i*num_wells + well_idx] != 0) {
                            msg += fmt::format("{} ", i);
                        }
                    }
                    OpmLog::info(msg);
                    ++well_idx;
                }
            }
        }
        std::ranges::sort(parallelWells);

        // Calling Schedule::filterConnections would remove any perforated
        // cells that exist only on other ranks even in the case of
        // distributed wells.  But we need all connections to figure out the
        // first cell of a well (e.g. for pressure).  Hence this is now
        // skipped.  Rank 0 had everything even before.

#if HAVE_OPENCL || HAVE_ROCSPARSE || HAVE_CUDA
        if (partitionJacobiBlocks) {
            this->cell_part_ = this->grid_->
                zoltanPartitionWithoutScatter(&wells, possibleFutureConnections, faceTrans.data(),
                                              numJacobiBlocks,
                                              imbalanceTol);
        }
#endif
    }
}

template<class ElementMapper, class GridView, class Scalar>
void GenericCpGridVanguard<ElementMapper,GridView,Scalar>::
distributeFieldProps_(EclipseState& eclState1)
{
    OPM_TIMEBLOCK(distributeFProps);
    const auto mpiSize = this->grid_->comm().size();

    if (mpiSize == 1) {
        return;
    }

    if (auto* parallelEclState = dynamic_cast<ParallelEclipseState*>(&eclState1);
        parallelEclState != nullptr)
    {
        // Reset Cartesian index mapper for automatic creation of field
        // properties
        // Note: the previous cartesianIndexMapper_ has been replaces by levelCartesianIndexMapper_
        // to support also the case of a distributed level zero grid in a CpGrid with LGRs.
        // This change allows access to the Cartesian indices of the distributed level zero grid,
        // where the field properties are given - for now.
        // In case of supporting LGR field properties, this need to be adapted, to access instead
        // each local/level Cartesian index set.
        parallelEclState->resetCartesianMapper(this->levelCartesianIndexMapper_.get());
        parallelEclState->switchToDistributedProps();
    }
    else {
        const auto message = std::string {
            "Parallel simulator setup is incorrect as "
            "it does not use ParallelEclipseState"
        };

        OpmLog::error(message);

        throw std::invalid_argument { message };
    }
}

template <class ElementMapper, class GridView, class Scalar>
std::vector<double>
GenericCpGridVanguard<ElementMapper, GridView, Scalar>::
extractFaceTrans(const GridView& gridView) const
{
    auto faceTrans = std::vector<double>(this->grid_->numFaces(), 0.0);

    const auto elemMapper = ElementMapper { gridView, Dune::mcmgElementLayout() };

    for (const auto& elem : elements(gridView, Dune::Partitions::interiorBorder)) {
        for (const auto& is : intersections(gridView, elem)) {
            if (!is.neighbor()) {
                continue;
            }

            const auto I = static_cast<unsigned int>(elemMapper.index(is.inside()));
            const auto J = static_cast<unsigned int>(elemMapper.index(is.outside()));

            faceTrans[is.id()] = this->getTransmissibility(I, J);
        }
    }

    return faceTrans;
}

template <class ElementMapper, class GridView, class Scalar>
void
GenericCpGridVanguard<ElementMapper, GridView, Scalar>::
distributeGrid(const Dune::EdgeWeightMethod                          edgeWeightsMethod,
               const bool                                            ownersFirst,
               const bool                                            addCorners,
               const int                                             numOverlap,
               const Dune::PartitionMethod                           partitionMethod,
               const bool                                            serialPartitioning,
               const bool                                            enableDistributedWells,
               const double                                          imbalanceTol,
               const bool                                            loadBalancerSet,
               const bool                                            useTransToFilterOverlap,
               const std::vector<double>&                            faceTrans,
               const std::vector<Well>&                              wells,
               const std::unordered_map<std::string, std::set<int>>& possibleFutureConnections,
               EclipseState&                                         eclState1,
               FlowGenericVanguard::ParallelWellStruct&              parallelWells)
{
    if (auto* eclState = dynamic_cast<ParallelEclipseState*>(&eclState1);
        eclState != nullptr)
    {
        this->distributeGrid(edgeWeightsMethod, ownersFirst, addCorners,
                             numOverlap, partitionMethod,
                             serialPartitioning, enableDistributedWells,
                             imbalanceTol, loadBalancerSet,
                             useTransToFilterOverlap, faceTrans,
                             wells, possibleFutureConnections,
                             eclState, parallelWells);
    }
    else {
        const auto message = std::string {
            "Parallel simulator setup is incorrect as "
            "it does not use ParallelEclipseState"
        };

        OpmLog::error(message);

        throw std::invalid_argument { message };
    }

    this->grid_->switchToDistributedView();
}

template <class ElementMapper, class GridView, class Scalar>
void
GenericCpGridVanguard<ElementMapper, GridView, Scalar>::
distributeGrid(const Dune::EdgeWeightMethod                          edgeWeightsMethod,
               const bool                                            ownersFirst,
               const bool                                            addCorners,
               const int                                             numOverlap,
               const Dune::PartitionMethod                           partitionMethod,
               const bool                                            serialPartitioning,
               const bool                                            enableDistributedWells,
               const double                                          imbalanceTol,
               const bool                                            loadBalancerSet,
               const bool                                            useTransToFilterOverlap,
               const std::vector<double>&                            faceTrans,
               const std::vector<Well>&                              wells,
               const std::unordered_map<std::string, std::set<int>>& possibleFutureConnections,
               ParallelEclipseState*                                 eclState,
               FlowGenericVanguard::ParallelWellStruct&              parallelWells)
{
    OPM_TIMEBLOCK(gridDistribute);
    const auto isIORank = this->grid_->comm().rank() == 0;

    PropsDataHandle<Dune::CpGrid> handle {
        *this->grid_, *eclState
    };

    const auto addCornerCells = addCorners;
    const auto overlapLayers = numOverlap;

    if (loadBalancerSet) {
        auto parts = isIORank
            ? (*externalLoadBalancer)(*this->grid_)
            : std::vector<int>{};

        // For this case, simple partitioning is selected automatically.
        parallelWells = std::get<1>
            (this->grid_->loadBalance(handle, parts, &wells,
                                      possibleFutureConnections, ownersFirst,
                                      addCornerCells, overlapLayers));
    }
    else {
        parallelWells = std::get<1>
            (this->grid_->loadBalance(handle, edgeWeightsMethod,
                                      &wells, possibleFutureConnections,
                                      serialPartitioning,
                                      faceTrans.data(), ownersFirst,
                                      addCornerCells, overlapLayers,
                                      partitionMethod, imbalanceTol,
                                      enableDistributedWells,
                                      useTransToFilterOverlap));
    }
}

#endif  // HAVE_MPI

template<class ElementMapper, class GridView, class Scalar>
void GenericCpGridVanguard<ElementMapper,GridView,Scalar>::
doCreateGrids_(const bool edge_conformal, const bool conformingLgr, EclipseState& eclState)
{
    const auto isRoot = this->mpiRank == 0;

    // At this stage the ParallelEclipseState instance is still in global
    // view; on rank 0 we have undistributed data for the entire grid, on
    // the other ranks the EclipseState is empty.
    const EclipseGrid* input_grid = nullptr;
    std::vector<double> global_porv;
    if (isRoot) {
        input_grid = &eclState.getInputGrid();
        global_porv = eclState.fieldProps().porv(true);
        OpmLog::info("\nProcessing grid");
    }

    // --- Create grid ---
    OPM_TIMEBLOCK(createGrids);

#if HAVE_MPI
    this->grid_ = std::make_unique<Dune::CpGrid>(FlowGenericVanguard::comm());
#else
    this->grid_ = std::make_unique<Dune::CpGrid>();
#endif
    // Chosen before processing, which keeps the corner-point input the Conforming builder needs.
    if (conformingLgr) {
        this->grid_->setLgrBackend(Opm::Refinement::Backend::Conforming);
    }

    // Note: removed_cells is guaranteed to be empty on ranks other than 0.
    auto removed_cells = this->grid_
        ->processEclipseFormat(input_grid,
                               &eclState,
                               /* isPeriodic = */ false,
                               /* flipNormals = */ false,
                               /* clipZ = */ false,
                               edge_conformal);

    if (isRoot) {
        const auto& active_porv = eclState.fieldProps().porv(false);
        const auto& unit_system = eclState.getUnits();
        const auto& volume_unit = unit_system.name(UnitSystem::measure::volume);
        double total_pore_volume = unit_system.from_si(UnitSystem::measure::volume,
                                                       std::accumulate(active_porv.begin(), active_porv.end(), 0.0));
        OpmLog::info(fmt::format("Total number of active cells: {} / total pore volume: {:0.0f} {}",
                                 grid_->numCells(), total_pore_volume , volume_unit));

        double removed_pore_volume =
            std::accumulate(removed_cells.begin(), removed_cells.end(), 0.0,
                            [&active_porv, &eclState](const auto acc, const auto global_index)
                            { return acc + active_porv[eclState.getInputGrid().activeIndex(global_index)]; });

        if (removed_pore_volume > 0) {
            removed_pore_volume = unit_system.from_si(UnitSystem::measure::volume, removed_pore_volume);
            OpmLog::info(fmt::format("Removed {} cells with a pore volume of {:0.0f} "
                                     "{} ({:5.3f} %) due to MINPV/MINPVV",
                                     removed_cells.size(),
                                     removed_pore_volume,
                                     volume_unit,
                                     100 * removed_pore_volume / total_pore_volume));
        }
    }

    cartesianIndexMapper_ = std::make_unique<CartesianIndexMapper>(*grid_);
    levelCartesianIndexMapper_ = std::make_unique<LevelCartesianIndexMapper>(*grid_);

#if HAVE_MPI
    if (this->grid_->comm().size() > 1) {
        // Numerical aquifers generate new NNCs during grid processing.  We
        // need to distribute these NNCs from the root process to all other
        // MPI processes.
        if (eclState.aquifer().hasNumericalAquifer()) {
            auto nnc_input = eclState.getInputNNC();
            Parallel::MpiSerializer ser(this->grid_->comm());
            ser.broadcast(Parallel::RootRank{0}, nnc_input);

            if (! isRoot) {
                eclState.setInputNNC(nnc_input);
            }
        }

        bool hasPinchNnc = eclState.hasPinchNNC();
        grid_->comm().broadcast(&hasPinchNnc, 1, 0);

        if (hasPinchNnc) {
            auto pinch_nnc = eclState.getPinchNNC();
            Parallel::MpiSerializer ser(this->grid_->comm());
            ser.broadcast(Parallel::RootRank{0}, pinch_nnc);

            if (! isRoot) {
                eclState.setPinchNNC(std::move(pinch_nnc));
            }
        }
    }
#endif

    // --- Copy grid to equilGrid_ ---
    // We use separate grid objects: one for the calculation of the initial
    // condition via EQUIL and one for the actual simulation. The reason is
    // that the EQUIL code is allergic to distributed grids and the
    // simulation grid is distributed before the initial condition is
    // calculated.
    //
    // Raised by the I/O reference grid's refinement, which runs on rank 0 alone.
    std::string outputGridError{};

    // After loadbalance, grid_ will contain a global and distribute view.
    // equilGrid_ being a shallow copy only the global view.
    if (isRoot) {
        this->equilGrid_ = std::make_unique<Dune::CpGrid>(*this->grid_);
        this->equilCartesianIndexMapper_ =
            std::make_unique<CartesianIndexMapper>(*this->equilGrid_);

#if HAVE_MPI
        // For a parallel run with LGRs, the simulation grid is refined only
        // after distribution (rank-interior refinement) and equilGrid_ stays
        // coarse.  ECL output gathers cell data onto the I/O rank using a
        // reference grid, which must therefore carry the refined leaf cells.
        // Build a self-communicator copy of the global grid and refine it
        // locally - getCommunicator() in the refinement builder follows the
        // grid's own communicator, so this is collective-free - giving the
        // I/O rank the same full refined grid a serial run would have.
        // equilGrid_ is left coarse for EQUIL; this is a separate grid used
        // only for output.
        if (this->grid_->comm().size() > 1 && input_grid != nullptr && conformingLgr) {
            if (const auto& lgrs = eclState.getLgrs(); lgrs.size() > 0) {
                // This backend refines from a retained corner-point
                // description rather than from the grid itself, so hand the
                // fresh grid the copy that grid_ (still undistributed here)
                // is holding.  Without it there is nothing to refine from and
                // the I/O rank is left without a refined reference grid.
                auto retained = Opm::Refinement::GridStateWriter::retainedCornerPointInput(
                    *this->grid_->currentData().front());
                if (retained) try {
                    auto outGrid = std::make_unique<Dune::CpGrid>(Dune::MPIHelper::getLocalCommunicator());
                    outGrid->setLgrBackend(Opm::Refinement::Backend::Conforming);
                    // Build it from the corner-point description grid_ was
                    // built from, not from the input grid again: MINPV works by
                    // collapsing ZCORN, and re-processing without an
                    // EclipseState skips MINPV and PINCH altogether (their whole
                    // block is gated on it). The reference grid then keeps cells
                    // the simulation grid does not have, and gathering onto the
                    // I/O rank leaves them claimed by no rank -- 399 of them on
                    // Norne, being the 496 MINPV removed less the 97 inside the
                    // refinement box, which the box refines away.
                    const grdecl raw {
                        { retained->dims[0], retained->dims[1], retained->dims[2] },
                        retained->coord.data(),
                        retained->zcorn.data(),
                        retained->actnum.empty() ? nullptr : retained->actnum.data()
                    };
                    outGrid->processEclipseFormat(raw, retained->pinchNnc,
                                                  retained->pinchActive,
                                                  retained->edgeConformal);
                    Opm::Refinement::GridStateWriter::setRetainedCornerPointInput(
                        *outGrid->currentData().front(), retained);
                    this->addLgrsUpdateLeafView(lgrs, lgrs.size(), *outGrid);
                    this->outputGrid_ = std::move(outGrid);
                    this->outputCartesianIndexMapper_ =
                        std::make_unique<CartesianIndexMapper>(*this->outputGrid_);
                }
                catch (const std::exception& e) {
                    outputGridError = e.what();
                }
            }
        }
#endif

        eclState.reset_actnum(UgGridHelpers::createACTNUM(*this->grid_));
        eclState.set_active_indices(this->grid_->globalCell());
    }

#if HAVE_MPI
    // Rank 0's refusal must stop every rank, which would otherwise wait below.
    if (this->grid_->comm().size() > 1) {
        Parallel::MpiSerializer ser(this->grid_->comm());
        ser.broadcast(Parallel::RootRank{0}, outputGridError);
        if (!outputGridError.empty()) {
            OPM_THROW(std::invalid_argument, outputGridError);
        }
    }
#endif

    {
        auto size = removed_cells.size();
        this->grid_->comm().broadcast(&size, 1, 0);

        if (! isRoot) {
            removed_cells.resize(size);
        }

        this->grid_->comm().broadcast(removed_cells.data(), static_cast<int>(size), 0);
    }

    // Inform the aquifer object that we might have removed/deactivated
    // cells as part of minimum pore-volume threshold processing.
    eclState.pruneDeactivatedAquiferConnections(removed_cells);
}

template<class ElementMapper, class GridView, class Scalar>
void GenericCpGridVanguard<ElementMapper,GridView,Scalar>::addLgrsUpdateLeafView(const LgrCollection& lgrCollection,
                                                                                 const int lgrsSize,
                                                                                 Dune::CpGrid& grid)
{
    std::vector<std::array<int,3>> cells_per_dim_vec;
    std::vector<std::array<int,3>> startIJK_vec;
    std::vector<std::array<int,3>> endIJK_vec;
    std::vector<std::string> lgrName_vec;
    std::vector<std::string> lgrParentName_vec;
    cells_per_dim_vec.reserve(lgrsSize);
    startIJK_vec.reserve(lgrsSize);
    endIJK_vec.reserve(lgrsSize);
    lgrName_vec.reserve(lgrsSize);
    lgrParentName_vec.reserve(lgrsSize);
    // A graded box (N*FIN/H*FIN) has no single subdivision factor, so it goes
    // through the request-taking overload with its column tables. Boxes are
    // graded per direction, and a direction the deck did not grade is left
    // empty so it keeps the uniform path.
    std::vector<Opm::Refinement::AxisSubdivision> subdivisions;
    subdivisions.reserve(3*lgrsSize);
    std::vector<std::vector<int>> minpvRemoved;
    minpvRemoved.reserve(lgrsSize);
    bool anyGraded = false;

    for (int lgr = 0; lgr < lgrsSize; ++lgr)
    {
        const auto lgrCarfin = lgrCollection.getLgr(lgr);
        cells_per_dim_vec.push_back({lgrCarfin.NX()/(lgrCarfin.I2() +1 - lgrCarfin.I1()),
                lgrCarfin.NY()/(lgrCarfin.J2() +1 - lgrCarfin.J1()), lgrCarfin.NZ()/(lgrCarfin.K2() +1 - lgrCarfin.K1())});
        startIJK_vec.push_back({lgrCarfin.I1(), lgrCarfin.J1(), lgrCarfin.K1()});
        endIJK_vec.push_back({lgrCarfin.I2()+1, lgrCarfin.J2()+1, lgrCarfin.K2()+1});
        lgrName_vec.emplace_back(lgrCarfin.NAME());
        lgrParentName_vec.emplace_back(lgrCarfin.PARENT_NAME());

        anyGraded = anyGraded || lgrCarfin.isGraded();
        minpvRemoved.push_back(lgrCarfin.minpvRemoved());
        for (std::size_t dim = 0; dim < 3; ++dim) {
            auto columns = lgrCarfin.refinedColumns(dim);
            subdivisions.push_back(Opm::Refinement::AxisSubdivision{
                std::move(columns.parentOffset),
                std::move(columns.fracLo),
                std::move(columns.fracHi) });
        }
    }

    // Common case: no nesting (every CARFIN refines GLOBAL). Take the exact
    // original code path - the 4-arg call, no reordering, no parent vector -
    // so single-level decks are completely unaffected by the nested-LGR
    // machinery below.
    const bool anyNested = std::any_of(lgrParentName_vec.begin(), lgrParentName_vec.end(),
                                       [](const std::string& p) { return p != "GLOBAL"; });
    // Build the requests from the pieces above; only a graded deck carries the
    // column tables, so a uniform deck reaches the builder exactly as before.
    const auto makeRequests = [&](const std::vector<int>& order) {
        std::vector<Opm::Refinement::BlockRefinement> requests;
        requests.reserve(order.size());
        for (const int lgr : order) {
            Opm::Refinement::BlockRefinement request;
            request.name = lgrName_vec[lgr];
            request.parentGridName = lgrParentName_vec[lgr];
            request.cellsPerDim = cells_per_dim_vec[lgr];
            request.startIJK = startIJK_vec[lgr];
            request.endIJK = endIJK_vec[lgr];
            if (anyGraded) {
                for (std::size_t dim = 0; dim < 3; ++dim) {
                    request.subdivision[dim] = subdivisions[3*lgr + dim];
                }
            }
            request.minpvRemoved = minpvRemoved[lgr];
            request.pillarsFromBoxLayer = lgrCollection.pillarsFromBoxLayer();
            requests.push_back(std::move(request));
        }
        return requests;
    };

    if (!anyNested) {
        std::vector<int> deckOrder(lgrsSize);
        std::iota(deckOrder.begin(), deckOrder.end(), 0);
        grid.addLgrsUpdateLeafView(makeRequests(deckOrder));
        return;
    }

    // Order the LGRs so a parent grid is always added before its children: the
    // conforming builder appends level grids in request order and resolves each
    // box's parent by name, so a child must follow its parent. Use a STABLE
    // topological order that preserves the deck order otherwise - a child is
    // only deferred until after its parent. This matters because the level
    // numbering produced here (CpGrid level = position in this order) must match
    // opm-common's LGR numbering, which is the deck order: the per-LGR ECL
    // output pairs each level's data with the LGR at the same position, so a
    // depth-grouping sort (which reorders sibling top-level LGRs around a nested
    // one) would mispair the output. A valid deck (parent before child) is left
    // exactly in deck order.
    std::map<std::string,int> indexOfName;
    for (int i = 0; i < lgrsSize; ++i) {
        indexOfName[lgrName_vec[i]] = i;
    }
    const auto parentIndex = [&](int i) -> int {
        if (lgrParentName_vec[i] == "GLOBAL") {
            return -1;
        }
        const auto it = indexOfName.find(lgrParentName_vec[i]);
        return (it == indexOfName.end()) ? -1 : it->second;  // unknown parent: builder reports it
    };
    std::vector<int> order;
    order.reserve(lgrsSize);
    std::vector<char> placed(lgrsSize, 0);
    bool progress = true;
    while (static_cast<int>(order.size()) < lgrsSize && progress) {
        progress = false;
        for (int i = 0; i < lgrsSize; ++i) {           // deck order within each pass
            if (placed[i]) {
                continue;
            }
            const int p = parentIndex(i);
            if (p < 0 || placed[p]) {                  // parent is GLOBAL/unknown or already placed
                order.push_back(i);
                placed[i] = 1;
                progress = true;
            }
        }
    }
    for (int i = 0; i < lgrsSize; ++i) {               // any leftovers (cycle): builder reports it
        if (!placed[i]) {
            order.push_back(i);
        }
    }

    grid.addLgrsUpdateLeafView(makeRequests(order));
};

template<class ElementMapper, class GridView, class Scalar>
const Dune::CpGrid&
GenericCpGridVanguard<ElementMapper,GridView,Scalar>::equilGrid() const
{
    assert(mpiRank == 0);
    return *equilGrid_;
}

template<class ElementMapper, class GridView, class Scalar>
const Dune::CartesianIndexMapper<Dune::CpGrid>&
GenericCpGridVanguard<ElementMapper,GridView,Scalar>::cartesianIndexMapper() const
{
    return *cartesianIndexMapper_;
}

template<class ElementMapper, class GridView, class Scalar>
const LevelCartesianIndexMapper<Dune::CpGrid>
GenericCpGridVanguard<ElementMapper,GridView,Scalar>::levelCartesianIndexMapper() const
{
    return *levelCartesianIndexMapper_;
}

template<class ElementMapper, class GridView, class Scalar>
const Dune::CartesianIndexMapper<Dune::CpGrid>&
GenericCpGridVanguard<ElementMapper,GridView,Scalar>::equilCartesianIndexMapper() const
{
    assert(mpiRank == 0);
    assert(equilCartesianIndexMapper_);
    return *equilCartesianIndexMapper_;
}

template<class ElementMapper, class GridView, class Scalar>
const Dune::CpGrid&
GenericCpGridVanguard<ElementMapper,GridView,Scalar>::eclOutputGrid() const
{
    assert(mpiRank == 0);
    return outputGrid_ ? *outputGrid_ : *equilGrid_;
}

template<class ElementMapper, class GridView, class Scalar>
const Dune::CartesianIndexMapper<Dune::CpGrid>&
GenericCpGridVanguard<ElementMapper,GridView,Scalar>::eclOutputCartesianIndexMapper() const
{
    assert(mpiRank == 0);
    return outputCartesianIndexMapper_ ? *outputCartesianIndexMapper_
                                       : *equilCartesianIndexMapper_;
}

template<class ElementMapper, class GridView, class Scalar>
Scalar
GenericCpGridVanguard<ElementMapper,GridView,Scalar>::
computeCellThickness(const typename GridView::template Codim<0>::Entity& element) const
{
    typedef typename Element::Geometry Geometry;
    static constexpr int zCoord = Element::dimension - 1;
    Scalar zz1 = 0.0;
    Scalar zz2 = 0.0;

    const Geometry& geometry = element.geometry();
    // This code only works with CP-grid where the
    // number of corners are 8 and
    // also assumes that the first
    // 4 corners are the top surface and
    // the 4 next are the bottomn.
    assert(geometry.corners() == 8);
    for (int i=0; i < 4; ++i){
        zz1 += geometry.corner(i)[zCoord];
        zz2 += geometry.corner(i+4)[zCoord];
    }
    zz1 /=4;
    zz2 /=4;
    return zz2-zz1;
}

#define INSTANTIATE_TYPE(T)                                      \
    template class GenericCpGridVanguard<                        \
        Dune::MultipleCodimMultipleGeomTypeMapper<               \
            Dune::GridView<                                      \
                Dune::DefaultLeafGridViewTraits<Dune::CpGrid>>>, \
        Dune::GridView<                                          \
            Dune::DefaultLeafGridViewTraits<Dune::CpGrid>>,      \
        T>;

INSTANTIATE_TYPE(double)
#if FLOW_INSTANTIATE_FLOAT
INSTANTIATE_TYPE(float)
#endif

#if HAVE_DUNE_FEM
using GV = Dune::Fem::AdaptiveLeafGridPart<Dune::CpGrid,
                                           (Dune::PartitionIteratorType)4,
                                           false>;
template class GenericCpGridVanguard<Dune::MultipleCodimMultipleGeomTypeMapper<GV>,
                                     GV,
                                     double>;
#if FLOW_INSTANTIATE_FLOAT
template class GenericCpGridVanguard<Dune::MultipleCodimMultipleGeomTypeMapper<GV>,
                                     GV,
                                     float>;
#endif

#endif // HAVE_DUNE_FEM

} // namespace Opm
