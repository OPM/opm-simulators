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

#include <dune/grid/common/mcmgmapper.hh>

#include <dune/common/version.hh>

#include <sstream>

namespace Opm {
template <class TypeTag>
class EclCpGridVanguard;
}

BEGIN_PROPERTIES

NEW_TYPE_TAG(EclCpGridVanguard, INHERITS_FROM(EclBaseVanguard));

// declare the properties
SET_TYPE_PROP(EclCpGridVanguard, Vanguard, Opm::EclCpGridVanguard<TypeTag>);
SET_TYPE_PROP(EclCpGridVanguard, Grid, Dune::CpGrid);
SET_TYPE_PROP(EclCpGridVanguard, EquilGrid, typename GET_PROP_TYPE(TypeTag, Grid));

END_PROPERTIES

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

    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, Simulator) Simulator;
    typedef typename GET_PROP_TYPE(TypeTag, ElementMapper) ElementMapper;

public:
    typedef typename GET_PROP_TYPE(TypeTag, Grid) Grid;
    typedef typename GET_PROP_TYPE(TypeTag, EquilGrid) EquilGrid;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;

private:
    typedef Dune::CartesianIndexMapper<Grid> CartesianIndexMapper;

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

            // convert to transmissibility for faces
            // TODO: grid_->numFaces() is not generic. use grid_->size(1) instead? (might
            // not work)
            const auto& gridView = grid_->leafGridView();
            unsigned numFaces = grid_->numFaces();
            std::vector<double> faceTrans(numFaces, 0.0);
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
                    defunctWellNames_ = std::get<1>(grid_->loadBalance(handle, edgeWeightsMethod, &wells, faceTrans.data()));
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

            if ( ! equilGrid_ )
            {
                // for processes that do not hold the global grid we filter here using the local grid.
                // If we would filter in filterConnection_ our partition would be empty and the connections of all
                // wells would be removed.
                ActiveGridCells activeCells(grid().logicalCartesianSize(),
                                            grid().globalCell().data(), grid().size(0));
                this->schedule().filterConnections(activeCells);
            }
        }
#endif

        cartesianIndexMapper_.reset(new CartesianIndexMapper(*grid_));
        this->updateGridView_();
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

    std::unordered_set<std::string> defunctWellNames() const
    { return defunctWellNames_; }

    const EclTransmissibility<TypeTag>& globalTransmissibility() const
    {
        assert( globalTrans_ != nullptr );
        return *globalTrans_;
    }

    void releaseGlobalTransmissibility()
    {
        globalTrans_.reset();
    }

protected:
    void createGrids_()
    {
        grid_.reset(new Dune::CpGrid());
        grid_->processEclipseFormat(mpiRank == 0 ? &this->eclState().getInputGrid()
                                                 : nullptr,
                                    /*isPeriodic=*/false,
                                    /*flipNormals=*/false,
                                    /*clipZ=*/false,
                                    mpiRank == 0 ? this->eclState().fieldProps().porv(true)
                                                 : std::vector<double>(),
                                    this->eclState().getInputNNC());

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
    }

    std::unique_ptr<Grid> grid_;
    std::unique_ptr<EquilGrid> equilGrid_;
    std::unique_ptr<CartesianIndexMapper> cartesianIndexMapper_;
    std::unique_ptr<CartesianIndexMapper> equilCartesianIndexMapper_;

    std::unique_ptr<EclTransmissibility<TypeTag> > globalTrans_;
    std::unordered_set<std::string> defunctWellNames_;
    int mpiRank;
};

} // namespace Opm

#endif
