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
*/
/*!
 * \file
 * \copydoc Ewoms::EclCpGridManager
 */
#ifndef EWOMS_ECL_CP_GRID_MANAGER_HH
#define EWOMS_ECL_CP_GRID_MANAGER_HH

#include "eclbasegridmanager.hh"

#include <dune/grid/CpGrid.hpp>

namespace Ewoms {
template <class TypeTag>
class EclCpGridManager;

namespace Properties {
NEW_TYPE_TAG(EclCpGridManager, INHERITS_FROM(EclBaseGridManager));

// declare the properties
SET_TYPE_PROP(EclCpGridManager, GridManager, Ewoms::EclCpGridManager<TypeTag>);
SET_TYPE_PROP(EclCpGridManager, Grid, Dune::CpGrid);
SET_TYPE_PROP(EclCpGridManager, EquilGrid, typename GET_PROP_TYPE(TypeTag, Grid));
} // namespace Properties

/*!
 * \ingroup EclBlackOilSimulator
 *
 * \brief Helper class for grid instantiation of ECL file-format using problems.
 *
 * This class uses Dune::CpGrid as the simulation grid.
 */
template <class TypeTag>
class EclCpGridManager : public EclBaseGridManager<TypeTag>
{
    friend class EclBaseGridManager<TypeTag>;
    typedef EclBaseGridManager<TypeTag> ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, Simulator) Simulator;

public:
    typedef typename GET_PROP_TYPE(TypeTag, Grid) Grid;
    typedef typename GET_PROP_TYPE(TypeTag, EquilGrid) EquilGrid;

private:
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;

    typedef Grid* GridPointer;
    typedef EquilGrid* EquilGridPointer;
    typedef Dune::CartesianIndexMapper<Grid> CartesianIndexMapper;
    typedef CartesianIndexMapper* CartesianIndexMapperPointer;

public:
    /*!
     * \brief Inherit the constructors from the base class.
     */
    using EclBaseGridManager<TypeTag>::EclBaseGridManager;

    ~EclCpGridManager()
    {
        delete cartesianIndexMapper_;
        delete grid_;
        delete equilGrid_;
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
    }

    /*!
     * \brief Distribute the simulation grid over multiple processes
     *
     * (For parallel simulation runs.)
     */
    void loadBalance()
    {
        // distribute the grid and switch to the distributed view
        grid_->loadBalance();
        grid_->switchToDistributedView();
        cartesianIndexMapper_ = new CartesianIndexMapper(*grid_);
    }

    /*!
     * \brief Returns the object which maps a global element index of the simulation grid
     *        to the corresponding element index of the logically Cartesian index.
     */
    const CartesianIndexMapper& cartesianIndexMapper() const
    { return *cartesianIndexMapper_; }

protected:
    void createGrids_()
    {
        std::vector<double> porv = this->eclState()->getDoubleGridProperty("PORV")->getData();

        grid_ = new Dune::CpGrid();
        grid_->processEclipseFormat(this->eclState()->getEclipseGrid(),
                                    /*isPeriodic=*/false,
                                    /*flipNormals=*/false,
                                    /*clipZ=*/false,
                                    porv);

        // we use separate grid objects: one for the calculation of the initial condition
        // via EQUIL and one for the actual simulation. The reason is that the EQUIL code
        // is allergic to distributed grids and the simulation grid is distributed before
        // the initial condition is calculated.
        equilGrid_ = new Dune::CpGrid();
        equilGrid_->processEclipseFormat(this->eclState()->getEclipseGrid(),
                                         /*isPeriodic=*/false,
                                         /*flipNormals=*/false,
                                         /*clipZ=*/false,
                                         porv);
    }

    GridPointer grid_;
    EquilGridPointer equilGrid_;
    CartesianIndexMapperPointer cartesianIndexMapper_;
};

} // namespace Ewoms

#endif
