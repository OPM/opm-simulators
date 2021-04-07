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
 * \copydoc Opm::EclCpGridVanguard
 */
#ifndef EWOMS_ECL_CP_GRID_VANGUARD_HH
#define EWOMS_ECL_CP_GRID_VANGUARD_HH

#include "eclbasevanguard.hh"
#include "ecltransmissibility.hh"
#include "femcpgridcompat.hh"

#include <opm/grid/CpGrid.hpp>
#include <opm/grid/cpgrid/GridHelpers.hpp>
#include <opm/simulators/utils/ParallelEclipseState.hpp>
#include <opm/simulators/utils/PropsCentroidsDataHandle.hpp>
#include <opm/simulators/utils/ParallelSerialization.hpp>

#include <ebos/eclmpiserializer.hh>

#include <dune/grid/common/mcmgmapper.hh>

#include <dune/common/version.hh>

#include <functional>
#include <numeric>
#include <sstream>

namespace Opm {
template <class TypeTag>
class EclCpGridVanguard;
}

namespace Opm::Properties {

namespace TTag {
struct EclCpGridVanguard {
    using InheritsFrom = std::tuple<EclBaseVanguard>;
};
}

// declare the properties
template<class TypeTag>
struct Vanguard<TypeTag, TTag::EclCpGridVanguard> {
    using type = Opm::EclCpGridVanguard<TypeTag>;
};
template<class TypeTag>
struct Grid<TypeTag, TTag::EclCpGridVanguard> {
    using type = Dune::CpGrid;
};
template<class TypeTag>
struct EquilGrid<TypeTag, TTag::EclCpGridVanguard> {
    using type = GetPropType<TypeTag, Properties::Grid>;
};

} // namespace Opm::Properties

namespace Opm {

/*!
 * \ingroup EclBlackOilSimulator
 *
 * \brief Helper class for grid instantiation of ECL file-format using problems.
 *
 * This class uses Dune::CpGrid as the simulation grid.
 */
template <class TypeTag>
class EclCpGridVanguard : public EclBaseVanguard<TypeTag>
{
    friend class EclBaseVanguard<TypeTag>;
    typedef EclBaseVanguard<TypeTag> ParentType;

    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Simulator = GetPropType<TypeTag, Properties::Simulator>;
    using ElementMapper = GetPropType<TypeTag, Properties::ElementMapper>;

public:
    using Grid = GetPropType<TypeTag, Properties::Grid>;
    using EquilGrid = GetPropType<TypeTag, Properties::EquilGrid>;
    using GridView = GetPropType<TypeTag, Properties::GridView>;

private:
    typedef Dune::CartesianIndexMapper<Grid> CartesianIndexMapper;
    using Element = typename GridView::template Codim<0>::Entity;

public:
    EclCpGridVanguard(Simulator& simulator)
        : EclBaseVanguard<TypeTag>(simulator), mpiRank()
    {
#if HAVE_MPI
        MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);
#endif
        this->callImplementationInit();
    }

    /*!
     * \brief Return a reference to the simulation grid.
     */
    Grid& grid()
    { return *grid_; }

    /*!
     * \brief Return a reference to the simulation grid.
     */
    const Grid& grid() const
    { return *grid_; }

    /*!
     * \brief Returns a refefence to the grid which should be used by the EQUIL
     *        initialization code.
     *
     * The EQUIL keyword is used to specify the initial condition of the reservoir in
     * hydrostatic equilibrium. Since the code which does this is not accepting arbitrary
     * DUNE grids (the code is part of the opm-core module), this is not necessarily the
     * same as the grid which is used for the actual simulation.
     */
    const EquilGrid& equilGrid() const
    {
        assert(mpiRank == 0);
        return *equilGrid_;
    }

    /*!
     * \brief Indicates that the initial condition has been computed and the memory used
     *        by the EQUIL grid can be released.
     *
     * Depending on the implementation, subsequent accesses to the EQUIL grid lead to
     * crashes.
     */
    void releaseEquilGrid()
    {
        equilGrid_.reset();
        equilCartesianIndexMapper_.reset();
    }

    /*!
     * \brief Distribute the simulation grid over multiple processes
     *
     * (For parallel simulation runs.)
     */
    void loadBalance()
    {
#if HAVE_MPI
        int mpiSize = 1;
        MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);

        if (mpiSize > 1) {
            // the CpGrid's loadBalance() method likes to have the transmissibilities as
            // its edge weights. since this is (kind of) a layering violation and
            // transmissibilities are relatively expensive to compute, we only do it if
            // more than a single process is involved in the simulation.
            cartesianIndexMapper_.reset(new CartesianIndexMapper(*grid_));
            if (grid_->size(0))
            {
                globalTrans_.reset(new EclTransmissibility<TypeTag>(*this));
                globalTrans_->update(false);
            }

            Dune::EdgeWeightMethod edgeWeightsMethod = this->edgeWeightsMethod();
            bool ownersFirst = this->ownersFirst();
            bool serialPartitioning = this->serialPartitioning();
            bool enableDistributedWells = this->enableDistributedWells();
            Scalar zoltanImbalanceTol = this->zoltanImbalanceTol();

            // convert to transmissibility for faces
            // TODO: grid_->numFaces() is not generic. use grid_->size(1) instead? (might
            // not work)
            const auto& gridView = grid_->leafGridView();
            unsigned numFaces = grid_->numFaces();
            std::vector<double> faceTrans;
            int loadBalancerSet = externalLoadBalancer_.has_value();
            grid_->comm().broadcast(&loadBalancerSet, 1, 0);
            if (!loadBalancerSet){
                faceTrans.resize(numFaces, 0.0);
                ElementMapper elemMapper(this->gridView(), Dune::mcmgElementLayout());
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

                        faceTrans[faceIdx] = globalTrans_->transmissibility(I, J);
                    }
                }
            }

            //distribute the grid and switch to the distributed view.
            {
                const auto wells = this->schedule().getWellsatEnd();

                try
                {
                    auto& eclState = dynamic_cast<ParallelEclipseState&>(this->eclState());
                    const EclipseGrid* eclGrid = nullptr;

                    if (grid_->comm().rank() == 0)
                    {
                        eclGrid = &this->eclState().getInputGrid();
                    }

                    PropsCentroidsDataHandle<Dune::CpGrid> handle(*grid_, eclState, eclGrid, this->centroids_,
                                                                  cartesianIndexMapper());
                    if (loadBalancerSet)
                    {
                        std::vector<int> parts;
                        if (grid_->comm().rank() == 0)
                        {
                            parts =  (*externalLoadBalancer_)(*grid_);
                        }
                        this->parallelWells_ = std::get<1>(grid_->loadBalance(handle, parts, &wells, ownersFirst, false, 1));
                    }
                    else
                    {
                        this->parallelWells_ =
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
#endif

        cartesianIndexMapper_.reset(new CartesianIndexMapper(*grid_));
        this->updateGridView_();
        this->updateCartesianToCompressedMapping_();
        this->updateCellDepths_();
        this->updateCellThickness_();

#if HAVE_MPI
        if (mpiSize > 1) {
            try
            {
                auto& parallelEclState = dynamic_cast<ParallelEclipseState&>(this->eclState());
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
#endif
    }

    /*!
     * \brief Free the memory occupied by the global transmissibility object.
     *
     * After writing the initial solution, this array should not be necessary anymore.
     */
    void releaseGlobalTransmissibilities()
    {
        globalTrans_.reset();
    }

    /*!
     * \brief Returns the object which maps a global element index of the simulation grid
     *        to the corresponding element index of the logically Cartesian index.
     */
    const CartesianIndexMapper& cartesianIndexMapper() const
    { return *cartesianIndexMapper_; }

    /*!
     * \brief Returns mapper from compressed to cartesian indices for the EQUIL grid
     */
    const CartesianIndexMapper& equilCartesianIndexMapper() const
    {
        assert(mpiRank == 0);
        assert(equilCartesianIndexMapper_);
        return *equilCartesianIndexMapper_;
    }

    const EclTransmissibility<TypeTag>& globalTransmissibility() const
    {
        assert( globalTrans_ != nullptr );
        return *globalTrans_;
    }

    void releaseGlobalTransmissibility()
    {
        globalTrans_.reset();
    }

    /// \brief Sets a function that returns external load balancing information when passed the grid
    ///
    /// The information is a vector of integers indication the partition index for each cell id.
    static void setExternalLoadBalancer(const std::function<std::vector<int> (const Grid&)>& loadBalancer)
    {
        externalLoadBalancer_ = loadBalancer;
    }
protected:
    void createGrids_()
    {
        const EclipseGrid * input_grid = nullptr;
        std::vector<double> global_porv;
        std::unordered_map<size_t, double> aquifer_cell_volumes;
        // At this stage the ParallelEclipseState instance is still in global
        // view; on rank 0 we have undistributed data for the entire grid, on
        // the other ranks the EclipseState is empty.
        if (mpiRank == 0) {
            input_grid = &this->eclState().getInputGrid();
            global_porv = this->eclState().fieldProps().porv(true);
            aquifer_cell_volumes = this->eclState().aquifer().numericalAquifers().aquiferCellVolumes();
            OpmLog::info("\nProcessing grid");
        }

        grid_.reset(new Dune::CpGrid());
        const auto& removed_cells = grid_->processEclipseFormat(input_grid,
                                                                /*isPeriodic=*/false,
                                                                /*flipNormals=*/false,
                                                                /*clipZ=*/false,
                                                                global_porv,
                                                                this->eclState().getInputNNC(),
                                                                aquifer_cell_volumes);

        if (mpiRank == 0) {
            const auto& active_porv = this->eclState().fieldProps().porv(false);
            const auto& unit_system = this->eclState().getUnits();
            const auto& volume_unit = unit_system.name( UnitSystem::measure::volume);
            double total_pore_volume = unit_system.from_si( UnitSystem::measure::volume, std::accumulate(active_porv.begin(), active_porv.end(), 0.0));
            OpmLog::info(fmt::format("Total number of active cells: {} / total pore volume: {:0.0f} {}", grid_->numCells(), total_pore_volume , volume_unit));

            double removed_pore_volume = 0;
            for (const auto& global_index : removed_cells)
                removed_pore_volume += active_porv[ input_grid->activeIndex(global_index) ];

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
            const bool has_numerical_aquifer = this->eclState().aquifer().hasNumericalAquifer();
            int mpiSize = 1;
            MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);
            // when there is numerical aquifers, new NNC are generated during grid processing
            // we need to pass the NNC from root process to other processes
            if (has_numerical_aquifer && mpiSize > 1) {
                auto nnc_input = this->eclState().getInputNNC();
                Opm::EclMpiSerializer ser(Dune::MPIHelper::getCollectiveCommunication());
                ser.template broadcast(nnc_input);
                if (mpiRank > 0) {
                    this->eclState().setInputNNC(nnc_input);
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

            std::vector<int> actnum = Opm::UgGridHelpers::createACTNUM(*grid_);
            auto &field_props = this->eclState().fieldProps();
            const_cast<FieldPropsManager&>(field_props).reset_actnum(actnum);
        }
    }

    // removing some connection located in inactive grid cells
    void filterConnections_()
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
            this->schedule().filterConnections(activeCells);
        }
#if HAVE_MPI
        try
        {
            // Broadcast another time to remove inactive peforations on
            // slave processors.
            Opm::eclScheduleBroadcast(this->schedule());
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

    Scalar computeCellThickness(const Element& element) const
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

    std::unique_ptr<Grid> grid_;
    std::unique_ptr<EquilGrid> equilGrid_;
    std::unique_ptr<CartesianIndexMapper> cartesianIndexMapper_;
    std::unique_ptr<CartesianIndexMapper> equilCartesianIndexMapper_;

    std::unique_ptr<EclTransmissibility<TypeTag> > globalTrans_;

    /// \brief optional functor returning external load balancing information
    ///
    /// If it is set then this will be used during loadbalance.
    static std::optional<std::function<std::vector<int> (const Grid&)>> externalLoadBalancer_;
    int mpiRank;
};

template<class TypeTag>
std::optional<std::function<std::vector<int>(const typename EclCpGridVanguard<TypeTag>::Grid&)>>
Opm::EclCpGridVanguard<TypeTag>::externalLoadBalancer_;

} // namespace Opm

#endif
