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
#include <ebos/eclgenericcpgridvanguard.hh>

#if HAVE_MPI
#include <ebos/eclmpiserializer.hh>
#endif

#include <opm/common/utility/ActiveGridCells.hpp>
#include <opm/grid/cpgrid/GridHelpers.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/Schedule.hpp>
#include <opm/simulators/utils/ParallelEclipseState.hpp>
#include <opm/simulators/utils/PropsCentroidsDataHandle.hpp>
#include <opm/simulators/utils/ParallelSerialization.hpp>

#include <dune/common/version.hh>
#include <dune/grid/common/mcmgmapper.hh>

#if HAVE_DUNE_FEM
#include <dune/fem/gridpart/adaptiveleafgridpart.hh>
#include <dune/fem/gridpart/common/gridpart2gridview.hh>
#include <ebos/femcpgridcompat.hh>
#endif

#include <fmt/format.h>

#include <cassert>
#include <numeric>
#include <sstream>

namespace Opm {

std::optional<std::function<std::vector<int> (const Dune::CpGrid&)>> externalLoadBalancer;

template<class ElementMapper, class GridView, class Scalar>
EclGenericCpGridVanguard<ElementMapper,GridView,Scalar>::EclGenericCpGridVanguard()
{
#if HAVE_MPI
    MPI_Comm_rank(EclGenericVanguard::comm(), &mpiRank);
#else
  mpiRank = 0;
#endif
}

template<class ElementMapper, class GridView, class Scalar>
void EclGenericCpGridVanguard<ElementMapper,GridView,Scalar>::releaseEquilGrid()
{
    equilGrid_.reset();
    equilCartesianIndexMapper_.reset();
}

#if HAVE_MPI
template<class ElementMapper, class GridView, class Scalar>
void EclGenericCpGridVanguard<ElementMapper,GridView,Scalar>::doLoadBalance_(Dune::EdgeWeightMethod edgeWeightsMethod,
                                                                             bool ownersFirst,
                                                                             bool serialPartitioning,
                                                                             bool enableDistributedWells,
                                                                             double zoltanImbalanceTol,
                                                                             const GridView& gridv,
                                                                             const Schedule& schedule,
                                                                             std::vector<double>& centroids,
                                                                             EclipseState& eclState1,
                                                                             EclGenericVanguard::ParallelWellStruct& parallelWells)
{
    int mpiSize = 1;
    MPI_Comm_size(grid_->comm(), &mpiSize);

    if (mpiSize > 1) {
        // the CpGrid's loadBalance() method likes to have the transmissibilities as
        // its edge weights. since this is (kind of) a layering violation and
        // transmissibilities are relatively expensive to compute, we only do it if
        // more than a single process is involved in the simulation.
        cartesianIndexMapper_.reset(new CartesianIndexMapper(*grid_));
        if (grid_->size(0))
        {
            this->allocTrans();
        }

        // convert to transmissibility for faces
        // TODO: grid_->numFaces() is not generic. use grid_->size(1) instead? (might
        // not work)
        const auto& gridView = grid_->leafGridView();
        unsigned numFaces = grid_->numFaces();
        std::vector<double> faceTrans;
        int loadBalancerSet = externalLoadBalancer.has_value();
        grid_->comm().broadcast(&loadBalancerSet, 1, 0);
        if (!loadBalancerSet){
            faceTrans.resize(numFaces, 0.0);
            ElementMapper elemMapper(gridv, Dune::mcmgElementLayout());
            auto elemIt = gridView.template begin</*codim=*/0>();
            const auto& elemEndIt = gridView.template end</*codim=*/0>();
            for (; elemIt != elemEndIt; ++ elemIt) {
                const auto& elem = *elemIt;
                auto isIt = gridView.ibegin(elem);
                const auto& isEndIt = gridView.iend(elem);
                for (; isIt != isEndIt; ++ isIt) {
                    const auto& is = *isIt;
                    if (!is.neighbor())
                        continue;

                    unsigned I = elemMapper.index(is.inside());
                    unsigned J = elemMapper.index(is.outside());

                    // FIXME (?): this is not portable!
                    unsigned faceIdx = is.id();

                    faceTrans[faceIdx] = this->getTransmissibility(I,J);
                }
            }
        }

        //distribute the grid and switch to the distributed view.
        {
            const auto wells = schedule.getWellsatEnd();

            try
            {
                auto& eclState = dynamic_cast<ParallelEclipseState&>(eclState1);
                const EclipseGrid* eclGrid = nullptr;

                if (grid_->comm().rank() == 0)
                {
                    eclGrid = &eclState.getInputGrid();
                }

                PropsCentroidsDataHandle<Dune::CpGrid> handle(*grid_, eclState, eclGrid, centroids,
                                                              cartesianIndexMapper());
                if (loadBalancerSet)
                {
                    std::vector<int> parts;
                    if (grid_->comm().rank() == 0)
                    {
                        parts =  (*externalLoadBalancer)(*grid_);
                    }
                    parallelWells = std::get<1>(grid_->loadBalance(handle, parts, &wells, ownersFirst, false, 1));
                }
                else
                {
                    parallelWells =
                        std::get<1>(grid_->loadBalance(handle, edgeWeightsMethod, &wells, serialPartitioning,
                                                       faceTrans.data(), ownersFirst, false, 1, true, zoltanImbalanceTol,
                                                       enableDistributedWells));
                }
            }
            catch(const std::bad_cast& e)
            {
                std::ostringstream message;
                message << "Parallel simulator setup is incorrect as it does not use ParallelEclipseState ("
                        << e.what() <<")"<<std::flush;
                OpmLog::error(message.str());
                std::rethrow_exception(std::current_exception());
            }
        }
        grid_->switchToDistributedView();

        cartesianIndexMapper_.reset();

        // Calling Schedule::filterConnections would remove any perforated
        // cells that exist only on other ranks even in the case of distributed wells
        // But we need all connections to figure out the first cell of a well (e.g. for
        // pressure). Hence this is now skipped. Rank 0 had everything even before.
    }
}

template<class ElementMapper, class GridView, class Scalar>
void EclGenericCpGridVanguard<ElementMapper,GridView,Scalar>::distributeFieldProps_(EclipseState& eclState1)
{
    int mpiSize = 1;
    MPI_Comm_size(grid_->comm(), &mpiSize);

    if (mpiSize > 1) {
        try
        {
            auto& parallelEclState = dynamic_cast<ParallelEclipseState&>(eclState1);
            // reset cartesian index mapper for auto creation of field properties
            parallelEclState.resetCartesianMapper(cartesianIndexMapper_.get());
            parallelEclState.switchToDistributedProps();
        }
        catch(const std::bad_cast& e)
        {
            std::ostringstream message;
            message << "Parallel simulator setup is incorrect as it does not use ParallelEclipseState ("
                            << e.what() <<")"<<std::flush;
            OpmLog::error(message.str());
            std::rethrow_exception(std::current_exception());
        }
    }
}
#endif

template<class ElementMapper, class GridView, class Scalar>
void EclGenericCpGridVanguard<ElementMapper,GridView,Scalar>::allocCartMapper()
{
    this->cartesianIndexMapper_.reset(new CartesianIndexMapper(this->grid()));
}

template<class ElementMapper, class GridView, class Scalar>
void EclGenericCpGridVanguard<ElementMapper,GridView,Scalar>::doCreateGrids_(EclipseState& eclState)
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

#if HAVE_MPI
    grid_.reset(new Dune::CpGrid(EclGenericVanguard::comm()));
#else
    grid_.reset(new Dune::CpGrid());
#endif

    const auto& removed_cells = grid_->processEclipseFormat(input_grid,
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
#if HAVE_MPI
    {
        const bool has_numerical_aquifer = eclState.aquifer().hasNumericalAquifer();
        int mpiSize = 1;
        MPI_Comm_size(grid_->comm(), &mpiSize);
        // when there is numerical aquifers, new NNC are generated during grid processing
        // we need to pass the NNC from root process to other processes
        if (has_numerical_aquifer && mpiSize > 1) {
            auto nnc_input = eclState.getInputNNC();
            EclMpiSerializer ser(grid_->comm());
            ser.broadcast(nnc_input);
            if (mpiRank > 0) {
                eclState.setInputNNC(nnc_input);
            }
        }
    }
#endif

    // we use separate grid objects: one for the calculation of the initial condition
    // via EQUIL and one for the actual simulation. The reason is that the EQUIL code
    // is allergic to distributed grids and the simulation grid is distributed before
    // the initial condition is calculated.
    // After loadbalance grid_ will contain a global and distribute view.
    // equilGrid_being a shallow copy only the global view.
    if (mpiRank == 0)
    {
        equilGrid_.reset(new Dune::CpGrid(*grid_));
        equilCartesianIndexMapper_.reset(new CartesianIndexMapper(*equilGrid_));

        std::vector<int> actnum = UgGridHelpers::createACTNUM(*grid_);
        auto &field_props = eclState.fieldProps();
        const_cast<FieldPropsManager&>(field_props).reset_actnum(actnum);
    }
}

template<class ElementMapper, class GridView, class Scalar>
void EclGenericCpGridVanguard<ElementMapper,GridView,Scalar>::doFilterConnections_(Schedule& schedule)
{
    // We only filter if we hold the global grid. Otherwise the filtering
    // is done after load balancing as in the future the other processes
    // will hold an empty partition for the global grid and hence filtering
    // here would remove all well connections.
    if (equilGrid_)
    {
        ActiveGridCells activeCells(equilGrid().logicalCartesianSize(),
                                    equilGrid().globalCell().data(),
                                    equilGrid().size(0));
        schedule.filterConnections(activeCells);
    }
#if HAVE_MPI
    try
    {
        // Broadcast another time to remove inactive peforations on
        // slave processors.
        eclScheduleBroadcast(EclGenericVanguard::comm(), schedule);
    }
    catch(const std::exception& broadcast_error)
    {
        OpmLog::error(fmt::format("Distributing properties to all processes failed\n"
                                  "Internal error message: {}", broadcast_error.what()));
        MPI_Finalize();
        std::exit(EXIT_FAILURE);
    }
#endif
}

template<class ElementMapper, class GridView, class Scalar>
const Dune::CpGrid& EclGenericCpGridVanguard<ElementMapper,GridView,Scalar>::equilGrid() const
{
    assert(mpiRank == 0);
    return *equilGrid_;
}

template<class ElementMapper, class GridView, class Scalar>
const Dune::CartesianIndexMapper<Dune::CpGrid>& EclGenericCpGridVanguard<ElementMapper,GridView,Scalar>::cartesianIndexMapper() const
{
    return *cartesianIndexMapper_;
}

template<class ElementMapper, class GridView, class Scalar>
const Dune::CartesianIndexMapper<Dune::CpGrid>& EclGenericCpGridVanguard<ElementMapper,GridView,Scalar>::equilCartesianIndexMapper() const
{
    assert(mpiRank == 0);
    assert(equilCartesianIndexMapper_);
    return *equilCartesianIndexMapper_;
}

template<class ElementMapper, class GridView, class Scalar>
Scalar EclGenericCpGridVanguard<ElementMapper,GridView,Scalar>::computeCellThickness(const typename GridView::template Codim<0>::Entity& element) const
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

#if HAVE_DUNE_FEM
template class EclGenericCpGridVanguard<Dune::MultipleCodimMultipleGeomTypeMapper<
                                        Dune::GridView<
                                        Dune::Fem::GridPart2GridViewTraits<
                                        Dune::Fem::AdaptiveLeafGridPart<Dune::CpGrid, Dune::PartitionIteratorType(4), false>>>>,
                                        Dune::GridView<
                                        Dune::Fem::GridPart2GridViewTraits<
                                        Dune::Fem::AdaptiveLeafGridPart<
                                        Dune::CpGrid, Dune::PartitionIteratorType(4), false>>>,
                                        double>;
#else
template class EclGenericCpGridVanguard<Dune::MultipleCodimMultipleGeomTypeMapper<Dune::GridView<Dune::DefaultLeafGridViewTraits<Dune::CpGrid>>>,
                                        Dune::GridView<Dune::DefaultLeafGridViewTraits<Dune::CpGrid>>,
                                        double>;
#endif
} // namespace Opm
