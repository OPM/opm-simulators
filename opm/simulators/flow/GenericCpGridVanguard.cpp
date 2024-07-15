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
#include <opm/simulators/flow/GenericCpGridVanguard.hpp>

#include <dune/grid/common/mcmgmapper.hh>
#include <dune/grid/common/partitionset.hh>
#include <dune/common/version.hh>

#include <opm/common/TimingMacros.hpp>
#include <opm/common/utility/ActiveGridCells.hpp>

#include <opm/grid/cpgrid/GridHelpers.hpp>

#include <opm/input/eclipse/Schedule/Schedule.hpp>
#include <opm/input/eclipse/Schedule/Well/Well.hpp>
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
#if !DUNE_VERSION_GTE(DUNE_FEM, 2, 9)
#include <dune/fem/gridpart/common/gridpart2gridview.hh>
#endif
#include <opm/simulators/flow/FemCpGridCompat.hpp>
#endif //HAVE_DUNE_FEM

#include <cassert>
#include <fstream>
#include <iterator>
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
               const Dune::PartitionMethod              partitionMethod,
               const bool                               serialPartitioning,
               const bool                               enableDistributedWells,
               const double                             imbalanceTol,
               const GridView&                          gridView,
               const Schedule&                          schedule,
               EclipseState&                            eclState1,
               FlowGenericVanguard::ParallelWellStruct& parallelWells,
               const int                                numJacobiBlocks)
{
    if (partitionMethod == Dune::PartitionMethod::zoltan && !this->zoltanParams().empty())
        this->grid_->setPartitioningParams(setupZoltanParams(this->zoltanParams()));
    if (partitionMethod == Dune::PartitionMethod::metis && !this->metisParams().empty())
        this->grid_->setPartitioningParams(setupMetisParams(this->metisParams()));

    const auto mpiSize = this->grid_->comm().size();

    const auto partitionJacobiBlocks =
        (numJacobiBlocks > 1) && (mpiSize == 1);

    if ((mpiSize > 1) || (numJacobiBlocks > 1)) {
        if (this->grid_->size(0) > 0) {
            // Generally needed in parallel runs both when there is and when
            // there is not an externally defined load-balancing function.
            // In addition to being used in CpGrid::loadBalance(), the
            // transmissibilities are also output to the .INIT file.  Thus,
            // transmissiblity values must exist on the I/O rank for derived
            // classes such as EclCpGridVanguard<>.
            this->allocTrans();
        }

        // CpGrid's loadBalance() method uses transmissibilities as edge
        // weights.  This is arguably a layering violation and extracting
        // the per-face transmissibilities as a linear array is relatively
        // expensive.  We therefore extract transmissibility values only if
        // the values are actually needed.
        auto loadBalancerSet = static_cast<int>(externalLoadBalancer.has_value());
        this->grid_->comm().broadcast(&loadBalancerSet, 1, 0);

        std::vector<double> faceTrans;
        {
            OPM_TIMEBLOCK(extractTrans);
            if (loadBalancerSet == 0 || partitionJacobiBlocks) {
                faceTrans = this->extractFaceTrans(gridView);
            }
        }

        const auto wells = ((mpiSize > 1) || partitionJacobiBlocks)
            ? schedule.getWellsatEnd()
            : std::vector<Well>{};

        // Distribute the grid and switch to the distributed view.
        if (mpiSize > 1) {
            this->distributeGrid(edgeWeightsMethod, ownersFirst, partitionMethod,
                                 serialPartitioning, enableDistributedWells,
                                 imbalanceTol, loadBalancerSet != 0,
                                 faceTrans, wells,
                                 eclState1, parallelWells);
        }

        // Calling Schedule::filterConnections would remove any perforated
        // cells that exist only on other ranks even in the case of
        // distributed wells.  But we need all connections to figure out the
        // first cell of a well (e.g. for pressure).  Hence this is now
        // skipped.  Rank 0 had everything even before.

#if HAVE_OPENCL || HAVE_ROCSPARSE || HAVE_CUDA
        if (partitionJacobiBlocks) {
            this->cell_part_ = this->grid_->
                zoltanPartitionWithoutScatter(&wells, faceTrans.data(),
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
        parallelEclState->resetCartesianMapper(this->cartesianIndexMapper_.get());
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
distributeGrid(const Dune::EdgeWeightMethod             edgeWeightsMethod,
               const bool                               ownersFirst,
               const Dune::PartitionMethod              partitionMethod,
               const bool                               serialPartitioning,
               const bool                               enableDistributedWells,
               const double                             imbalanceTol,
               const bool                               loadBalancerSet,
               const std::vector<double>&               faceTrans,
               const std::vector<Well>&                 wells,
               EclipseState&                            eclState1,
               FlowGenericVanguard::ParallelWellStruct& parallelWells)
{
    if (auto* eclState = dynamic_cast<ParallelEclipseState*>(&eclState1);
        eclState != nullptr)
    {
        this->distributeGrid(edgeWeightsMethod, ownersFirst, partitionMethod,
                             serialPartitioning, enableDistributedWells,
                             imbalanceTol, loadBalancerSet, faceTrans,
                             wells, eclState, parallelWells);
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
distributeGrid(const Dune::EdgeWeightMethod             edgeWeightsMethod,
               const bool                               ownersFirst,
               const Dune::PartitionMethod              partitionMethod,
               const bool                               serialPartitioning,
               const bool                               enableDistributedWells,
               const double                             imbalanceTol,
               const bool                               loadBalancerSet,
               const std::vector<double>&               faceTrans,
               const std::vector<Well>&                 wells,
               ParallelEclipseState*                    eclState,
               FlowGenericVanguard::ParallelWellStruct& parallelWells)
{
    OPM_TIMEBLOCK(gridDistribute);
    const auto isIORank = this->grid_->comm().rank() == 0;

    PropsDataHandle<Dune::CpGrid> handle {
        *this->grid_, *eclState
    };

    const auto addCornerCells = false;
    const auto overlapLayers = 1;

    if (loadBalancerSet) {
        auto parts = isIORank
            ? (*externalLoadBalancer)(*this->grid_)
            : std::vector<int>{};
        //For this case, simple partitioning is selected automatically
        parallelWells =
            std::get<1>(this->grid_->loadBalance(handle, parts, &wells, ownersFirst,
                                                 addCornerCells, overlapLayers));
    }
    else {
        parallelWells =
            std::get<1>(this->grid_->loadBalance(handle, edgeWeightsMethod,
                                                 &wells, serialPartitioning,
                                                 faceTrans.data(), ownersFirst,
                                                 addCornerCells, overlapLayers,
                                                 partitionMethod, imbalanceTol,
                                                 enableDistributedWells));
    }
}

#endif  // HAVE_MPI

template<class ElementMapper, class GridView, class Scalar>
void GenericCpGridVanguard<ElementMapper,GridView,Scalar>::doCreateGrids_(EclipseState& eclState)
{
    const EclipseGrid* input_grid = nullptr;
    std::vector<double> global_porv;
    // At this stage the ParallelEclipseState instance is still in global
    // view; on rank 0 we have undistributed data for the entire grid, on
    // the other ranks the EclipseState is empty.
    if (mpiRank == 0) {
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

    // Note: removed_cells is guaranteed to be empty on ranks other than 0.
    auto removed_cells =
        this->grid_->processEclipseFormat(input_grid,
                                          &eclState,
                                          /*isPeriodic=*/false,
                                          /*flipNormals=*/false,
                                          /*clipZ=*/false);

    if (mpiRank == 0) {
        const auto& active_porv = eclState.fieldProps().porv(false);
        const auto& unit_system = eclState.getUnits();
        const auto& volume_unit = unit_system.name( UnitSystem::measure::volume);
        double total_pore_volume = unit_system.from_si( UnitSystem::measure::volume, std::accumulate(active_porv.begin(), active_porv.end(), 0.0));
        OpmLog::info(fmt::format("Total number of active cells: {} / total pore volume: {:0.0f} {}", grid_->numCells(), total_pore_volume , volume_unit));

        double removed_pore_volume = 0;
        for (const auto& global_index : removed_cells)
            removed_pore_volume += active_porv[ eclState.getInputGrid().activeIndex(global_index) ];

        if (removed_pore_volume > 0) {
            removed_pore_volume = unit_system.from_si( UnitSystem::measure::volume, removed_pore_volume );
            OpmLog::info(fmt::format("Removed {} cells with a pore volume of {:0.0f} {} ({:5.3f} %) due to MINPV/MINPVV",
                                     removed_cells.size(),
                                     removed_pore_volume,
                                     volume_unit,
                                     100 * removed_pore_volume / total_pore_volume));
        }
    }

    cartesianIndexMapper_ = std::make_unique<CartesianIndexMapper>(*grid_);

    // --- Add LGRs and update Leaf Grid View ---
    // Check if input file contains Lgrs.
    const auto& lgrs = eclState.getLgrs();
    const auto lgrsSize = lgrs.size();
    // If there are lgrs, create the grid with them, and update the leaf grid view.
    if (lgrsSize)
    {
        OpmLog::info("\nAdding LGRs to the grid and updating its leaf grid view");
        this->addLgrsUpdateLeafView(lgrs, lgrsSize);
    }

#if HAVE_MPI
    {
        const bool has_numerical_aquifer = eclState.aquifer().hasNumericalAquifer();
        int mpiSize = 1;
        MPI_Comm_size(grid_->comm(), &mpiSize);

        // when there is numerical aquifers, new NNC are generated during
        // grid processing we need to pass the NNC from root process to
        // other processes
        if (has_numerical_aquifer && mpiSize > 1) {
            auto nnc_input = eclState.getInputNNC();
            Parallel::MpiSerializer ser(grid_->comm());
            ser.broadcast(nnc_input);
            if (mpiRank > 0) {
                eclState.setInputNNC(nnc_input);
            }
        }
    }
#endif

    // --- Copy grid with LGRs to equilGrid_ ---
    // We use separate grid objects: one for the calculation of the initial
    // condition via EQUIL and one for the actual simulation. The reason is
    // that the EQUIL code is allergic to distributed grids and the
    // simulation grid is distributed before the initial condition is
    // calculated.
    //
    // After loadbalance, grid_ will contain a global and distribute view.
    // equilGrid_ being a shallow copy only the global view.
    if (mpiRank == 0)
    {
        equilGrid_.reset(new Dune::CpGrid(*grid_));
        equilCartesianIndexMapper_ = std::make_unique<CartesianIndexMapper>(*equilGrid_);

        eclState.reset_actnum(UgGridHelpers::createACTNUM(*grid_));
    }

    {
        auto size = removed_cells.size();

        this->grid_->comm().broadcast(&size, 1, 0);

        if (mpiRank != 0) {
            removed_cells.resize(size);
        }

        this->grid_->comm().broadcast(removed_cells.data(), size, 0);
    }

    // Inform the aquifer object that we might have removed/deactivated
    // cells as part of minimum pore-volume threshold processing.
    eclState.pruneDeactivatedAquiferConnections(removed_cells);
}

template<class ElementMapper, class GridView, class Scalar>
void GenericCpGridVanguard<ElementMapper,GridView,Scalar>::addLgrsUpdateLeafView(const LgrCollection& lgrCollection, const int lgrsSize)
{
    std::vector<std::array<int,3>> cells_per_dim_vec;
    std::vector<std::array<int,3>> startIJK_vec;
    std::vector<std::array<int,3>> endIJK_vec;
    std::vector<std::string> lgrName_vec;
    cells_per_dim_vec.reserve(lgrsSize);
    startIJK_vec.reserve(lgrsSize);
    endIJK_vec.reserve(lgrsSize);
    lgrName_vec.reserve(lgrsSize);
    for (int lgr = 0; lgr < lgrsSize; ++lgr)
    {
        const auto lgrCarfin = lgrCollection.getLgr(lgr);
        cells_per_dim_vec.push_back({lgrCarfin.NX()/(lgrCarfin.I2() +1 - lgrCarfin.I1()),
                lgrCarfin.NY()/(lgrCarfin.J2() +1 - lgrCarfin.J1()), lgrCarfin.NZ()/(lgrCarfin.K2() +1 - lgrCarfin.K1())});
        startIJK_vec.push_back({lgrCarfin.I1(), lgrCarfin.J1(), lgrCarfin.K1()});
        endIJK_vec.push_back({lgrCarfin.I2()+1, lgrCarfin.J2()+1, lgrCarfin.K2()+1});
        lgrName_vec.emplace_back(lgrCarfin.NAME());
    }
    this->grid_->addLgrsUpdateLeafView(cells_per_dim_vec, startIJK_vec, endIJK_vec, lgrName_vec);
};

template<class ElementMapper, class GridView, class Scalar>
void GenericCpGridVanguard<ElementMapper,GridView,Scalar>::
doFilterConnections_(Schedule& schedule)
{
    // We only filter if we hold the global grid. Otherwise the filtering
    // is done after load balancing as in the future the other processes
    // will hold an empty partition for the global grid and hence filtering
    // here would remove all well connections.
    if (this->equilGrid_ != nullptr) {
        ActiveGridCells activeCells(equilGrid().logicalCartesianSize(),
                                    equilGrid().globalCell().data(),
                                    equilGrid().size(0));

        schedule.filterConnections(activeCells);
    }

#if HAVE_MPI
    try {
        // Broadcast another time to remove inactive peforations on
        // slave processors.
        eclBroadcast(FlowGenericVanguard::comm(), schedule);
    }
    catch (const std::exception& broadcast_error) {
        OpmLog::error(fmt::format("Distributing properties to all processes failed\n"
                                  "Internal error message: {}", broadcast_error.what()));
        MPI_Finalize();
        std::exit(EXIT_FAILURE);
    }
#endif  // HAVE_MPI
}

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
const Dune::CartesianIndexMapper<Dune::CpGrid>&
GenericCpGridVanguard<ElementMapper,GridView,Scalar>::equilCartesianIndexMapper() const
{
    assert(mpiRank == 0);
    assert(equilCartesianIndexMapper_);
    return *equilCartesianIndexMapper_;
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
template class GenericCpGridVanguard<
    Dune::MultipleCodimMultipleGeomTypeMapper<
        Dune::GridView<
            Dune::DefaultLeafGridViewTraits<Dune::CpGrid>>>,
    Dune::GridView<
        Dune::DefaultLeafGridViewTraits<Dune::CpGrid>>,
    double>;

#if HAVE_DUNE_FEM
#if DUNE_VERSION_GTE(DUNE_FEM, 2, 9)
using GV = Dune::Fem::AdaptiveLeafGridPart<Dune::CpGrid,
                                           (Dune::PartitionIteratorType)4,
                                           false>;
template class GenericCpGridVanguard<Dune::MultipleCodimMultipleGeomTypeMapper<GV>,
                                     GV,
                                     double>;
#else
template class GenericCpGridVanguard<
    Dune::MultipleCodimMultipleGeomTypeMapper<
        Dune::GridView<
            Dune::Fem::GridPart2GridViewTraits<
                Dune::Fem::AdaptiveLeafGridPart<
                    Dune::CpGrid,
                    Dune::PartitionIteratorType(4),
                    false>>>>,
    Dune::GridView<
        Dune::Fem::GridPart2GridViewTraits<
            Dune::Fem::AdaptiveLeafGridPart<
                Dune::CpGrid,
                Dune::PartitionIteratorType(4),
                false>>>,
    double>;

template class GenericCpGridVanguard<
    Dune::MultipleCodimMultipleGeomTypeMapper<
        Dune::Fem::GridPart2GridViewImpl<
            Dune::Fem::AdaptiveLeafGridPart<
                Dune::CpGrid,
                Dune::PartitionIteratorType(4),
                false>>>,
    Dune::Fem::GridPart2GridViewImpl<
        Dune::Fem::AdaptiveLeafGridPart<
            Dune::CpGrid,
            Dune::PartitionIteratorType(4),
            false> >,
    double>;
#endif
#endif // HAVE_DUNE_FEM

} // namespace Opm
