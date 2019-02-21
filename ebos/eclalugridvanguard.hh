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
 * \copydoc Ewoms::EclAluGridVanguard
 */
#ifndef EWOMS_ECL_ALU_GRID_VANGUARD_HH
#define EWOMS_ECL_ALU_GRID_VANGUARD_HH

#include "eclbasevanguard.hh"
#include "alucartesianindexmapper.hh"

#include <dune/alugrid/grid.hh>
#include <dune/alugrid/common/fromtogridfactory.hh>
#include <opm/grid/CpGrid.hpp>

namespace Ewoms {
template <class TypeTag>
class EclAluGridVanguard;

} // namespace Ewoms

BEGIN_PROPERTIES

NEW_TYPE_TAG(EclAluGridVanguard, INHERITS_FROM(EclBaseVanguard));

// declare the properties
SET_TYPE_PROP(EclAluGridVanguard, Vanguard, Ewoms::EclAluGridVanguard<TypeTag>);
SET_TYPE_PROP(EclAluGridVanguard, Grid,  Dune::ALUGrid<3, 3, Dune::cube, Dune::nonconforming>);
SET_TYPE_PROP(EclAluGridVanguard, EquilGrid, Dune::CpGrid);

END_PROPERTIES

namespace Ewoms {

/*!
 * \ingroup EclBlackOilSimulator
 *
 * \brief Helper class for grid instantiation of ECL file-format using problems.
 *
 * This class uses Dune::ALUGrid as the simulation grid.
 */
template <class TypeTag>
class EclAluGridVanguard : public EclBaseVanguard<TypeTag>
{
    friend class EclBaseVanguard<TypeTag>;
    typedef EclBaseVanguard<TypeTag> ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, Simulator) Simulator;

public:
    typedef typename GET_PROP_TYPE(TypeTag, Grid) Grid;
    typedef typename GET_PROP_TYPE(TypeTag, EquilGrid) EquilGrid;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;

private:
    typedef Ewoms::AluCartesianIndexMapper<Grid> CartesianIndexMapper;
    typedef Dune::CartesianIndexMapper<EquilGrid> EquilCartesianIndexMapper;

    static const int dimension = Grid::dimension;

public:
    /*!
     * \brief Inherit the constructors from the base class.
     */
    using EclBaseVanguard<TypeTag>::EclBaseVanguard;

    ~EclAluGridVanguard()
    {
        delete cartesianIndexMapper_;
        delete equilCartesianIndexMapper_;
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
        delete equilCartesianIndexMapper_;
        equilCartesianIndexMapper_ = 0;

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
        auto gridView = grid().leafGridView();
        auto dataHandle = cartesianIndexMapper_->dataHandle(gridView);
        grid().loadBalance(*dataHandle);

        // communicate non-interior cells values
        grid().communicate(*dataHandle,
                           Dune::InteriorBorder_All_Interface,
                           Dune::ForwardCommunication );
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
    const EquilCartesianIndexMapper& equilCartesianIndexMapper() const
    { return *equilCartesianIndexMapper_; }

protected:
    void createGrids_()
    {
        const auto& gridProps = this->eclState().get3DProperties();
        const std::vector<double>& porv = gridProps.getDoubleGridProperty("PORV").getData();

        // we use separate grid objects: one for the calculation of the initial condition
        // via EQUIL and one for the actual simulation. The reason is that the EQUIL code
        // cannot cope with arbitrary Dune grids and is also allergic to distributed
        // grids.

        /////
        // create the EQUIL grid
        /////
        equilGrid_ = new EquilGrid();
        equilGrid_->processEclipseFormat(this->eclState().getInputGrid(),
                                         /*isPeriodic=*/false,
                                         /*flipNormals=*/false,
                                         /*clipZ=*/false,
                                         porv);

        cartesianCellId_ = equilGrid_->globalCell();

        for (unsigned i = 0; i < dimension; ++i)
            cartesianDimension_[i] = equilGrid_->logicalCartesianSize()[i];

        equilCartesianIndexMapper_ = new EquilCartesianIndexMapper(*equilGrid_);

        /////
        // create the simulation grid
        /////
        Dune::FromToGridFactory<Grid> factory;
        grid_ = factory.convert(*equilGrid_, cartesianCellId_);

        cartesianIndexMapper_ =
            new CartesianIndexMapper(*grid_, cartesianDimension_, cartesianCellId_);
    }

    void filterConnections_()
    {
        // not handling the removal of completions for this type of grid yet.
    }

    Grid* grid_;
    EquilGrid* equilGrid_;
    std::vector<int> cartesianCellId_;
    std::array<int,dimension> cartesianDimension_;
    CartesianIndexMapper* cartesianIndexMapper_;
    EquilCartesianIndexMapper* equilCartesianIndexMapper_;
};

} // namespace Ewoms

#endif
