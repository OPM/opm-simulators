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
 * \copydoc Ewoms::EclCpGridVanguard
 */
#ifndef EWOMS_ECL_CP_GRID_VANGUARD_HH
#define EWOMS_ECL_CP_GRID_VANGUARD_HH

#include "eclbasevanguard.hh"
#include "ecltransmissibility.hh"
#include "femcpgridcompat.hh"

#include <opm/grid/CpGrid.hpp>
#include <opm/grid/cpgrid/GridHelpers.hpp>

#include <dune/grid/common/mcmgmapper.hh>

#include <dune/common/version.hh>

namespace Ewoms {
template <class TypeTag>
class EclCpGridVanguard;
}

BEGIN_PROPERTIES

NEW_TYPE_TAG(EclCpGridVanguard, INHERITS_FROM(EclBaseVanguard));

// declare the properties
SET_TYPE_PROP(EclCpGridVanguard, Vanguard, Ewoms::EclCpGridVanguard<TypeTag>);
SET_TYPE_PROP(EclCpGridVanguard, Grid, Dune::CpGrid);
SET_TYPE_PROP(EclCpGridVanguard, EquilGrid, typename GET_PROP_TYPE(TypeTag, Grid));

END_PROPERTIES

namespace Ewoms {

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
    typedef typename GET_PROP_TYPE(TypeTag, ElementMapper) ElementMapper;

public:
    typedef typename GET_PROP_TYPE(TypeTag, Grid) Grid;
    typedef typename GET_PROP_TYPE(TypeTag, EquilGrid) EquilGrid;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;

private:
    typedef Dune::CartesianIndexMapper<Grid> CartesianIndexMapper;

public:
    /*!
     * \brief Inherit the constructors from the base class.
     */
    using EclBaseVanguard<TypeTag>::EclBaseVanguard;

    ~EclCpGridVanguard()
    {
        delete cartesianIndexMapper_;
        delete equilCartesianIndexMapper_;
        delete grid_;
        delete equilGrid_;
        delete globalTrans_;
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
    { return *equilGrid_; }

    /*!
     * \brief Indicates that the initial condition has been computed and the memory used
     *        by the EQUIL grid can be released.
     *
     * Depending on the implementation, subsequent accesses to the EQUIL grid lead to
     * crashes.
     */
    void releaseEquilGrid()
    {
        delete equilGrid_;
        equilGrid_ = 0;

        delete equilCartesianIndexMapper_;
        equilCartesianIndexMapper_ = 0;
    }

    /*!
     * \brief Distribute the simulation grid over multiple processes
     *
     * (For parallel simulation runs.)
     */
    void loadBalance()
    {
#if HAVE_MPI
        int mpiRank = 0;
        int mpiSize = 1;
        MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);
        MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);

        if (mpiSize > 1) {
            // the CpGrid's loadBalance() method likes to have the transmissibilities as
            // its edge weights. since this is (kind of) a layering violation and
            // transmissibilities are relatively expensive to compute, we only do it if
            // more than a single process is involved in the simulation.
            cartesianIndexMapper_ = new CartesianIndexMapper(*grid_);
            globalTrans_ = new EclTransmissibility<TypeTag>(*this);
            globalTrans_->update();

            // convert to transmissibility for faces
            // TODO: grid_->numFaces() is not generic. use grid_->size(1) instead? (might
            // not work)
            const auto& gridView = grid_->leafGridView();
            unsigned numFaces = grid_->numFaces();
            std::vector<double> faceTrans(numFaces, 0.0);
#if DUNE_VERSION_NEWER(DUNE_GRID, 2,6)
            ElementMapper elemMapper(this->gridView(), Dune::mcmgElementLayout());
#else
            ElementMapper elemMapper(this->gridView());
#endif
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
                const auto wells = this->schedule().getWells2atEnd();
                defunctWellNames_ = std::get<1>(grid_->loadBalance(&wells, faceTrans.data()));
            }
            grid_->switchToDistributedView();

            delete cartesianIndexMapper_;
            cartesianIndexMapper_ = nullptr;
        }
#endif

        cartesianIndexMapper_ = new CartesianIndexMapper(*grid_);

        this->updateGridView_();
    }

    /*!
     * \brief Free the memory occupied by the global transmissibility object.
     *
     * After writing the initial solution, this array should not be necessary anymore.
     */
    void releaseGlobalTransmissibilities()
    {
        delete globalTrans_;
        globalTrans_ = nullptr;
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
    { return *equilCartesianIndexMapper_; }

    std::unordered_set<std::string> defunctWellNames() const
    { return defunctWellNames_; }

    const EclTransmissibility<TypeTag>& globalTransmissibility() const
    { return *globalTrans_; }

    void releaseGlobalTransmissibility()
    {
        delete globalTrans_;
        globalTrans_ = nullptr;
    }

protected:
    void createGrids_()
    {
        const auto& gridProps = this->eclState().get3DProperties();
        const std::vector<double>& porv = gridProps.getDoubleGridProperty("PORV").getData();

        grid_ = new Dune::CpGrid();
        grid_->processEclipseFormat(this->eclState().getInputGrid(),
                                    /*isPeriodic=*/false,
                                    /*flipNormals=*/false,
                                    /*clipZ=*/false,
                                    porv,
                                    this->eclState().getInputNNC());

        // we use separate grid objects: one for the calculation of the initial condition
        // via EQUIL and one for the actual simulation. The reason is that the EQUIL code
        // is allergic to distributed grids and the simulation grid is distributed before
        // the initial condition is calculated.
        // After loadbalance grid_ will contain a global and distribute view.
        // equilGrid_being a shallow copy only the global view.
        equilGrid_ = new Dune::CpGrid(*grid_);
        equilCartesianIndexMapper_ = new CartesianIndexMapper(*equilGrid_);

        globalTrans_ = nullptr;
    }

    // removing some connection located in inactive grid cells
    void filterConnections_()
    {
        assert(grid_);
        Grid grid = *grid_;
        grid.switchToGlobalView();
        const auto eclipseGrid = Opm::UgGridHelpers::createEclipseGrid(grid, this->eclState().getInputGrid());
        this->schedule().filterConnections(eclipseGrid);
    }

    Grid* grid_;
    EquilGrid* equilGrid_;
    CartesianIndexMapper* cartesianIndexMapper_;
    CartesianIndexMapper* equilCartesianIndexMapper_;

    EclTransmissibility<TypeTag>* globalTrans_;
    std::unordered_set<std::string> defunctWellNames_;
};

} // namespace Ewoms

#endif
