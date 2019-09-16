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
 * \copydoc Opm::SimplexGridVanguard
 */
#ifndef EWOMS_SIMPLEX_GRID_VANGUARD_HH
#define EWOMS_SIMPLEX_GRID_VANGUARD_HH

#include <opm/models/utils/basicproperties.hh>
#include <opm/models/utils/propertysystem.hh>
#include <opm/models/utils/parametersystem.hh>

#include <dune/grid/utility/structuredgridfactory.hh>
#include <dune/common/fvector.hh>

#include <memory>

BEGIN_PROPERTIES

NEW_PROP_TAG(Scalar);
NEW_PROP_TAG(Grid);

NEW_PROP_TAG(DomainSizeX);
NEW_PROP_TAG(DomainSizeY);
NEW_PROP_TAG(DomainSizeZ);

NEW_PROP_TAG(CellsX);
NEW_PROP_TAG(CellsY);
NEW_PROP_TAG(CellsZ);

NEW_PROP_TAG(GridGlobalRefinements);

END_PROPERTIES

namespace Opm {
/*!
 * \brief Provides a simulator vanguard which a creates regular grid made of simplices.
 */
template <class TypeTag>
class SimplexGridVanguard
{
    typedef BaseVanguard<TypeTag> ParentType;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, Simulator) Simulator;
    typedef typename GET_PROP_TYPE(TypeTag, Grid) Grid;

    typedef Dune::shared_ptr<Grid> GridPointer;
    typedef typename Grid::ctype CoordScalar;
    enum { dimWorld = Grid::dimensionworld };
    typedef Dune::FieldVector<CoordScalar, dimWorld> GlobalPosition;

public:
    /*!
     * \brief Register all run-time parameters for the grid manager.
     */
    static void registerParameters()
    {
        EWOMS_REGISTER_PARAM(TypeTag, unsigned, GridGlobalRefinements,
                             "The number of global refinements of the grid "
                             "executed after it was loaded");
        EWOMS_REGISTER_PARAM(TypeTag, Scalar, DomainSizeX,
                             "The size of the domain in x direction");
        EWOMS_REGISTER_PARAM(TypeTag, unsigned, CellsX,
                             "The number of intervalls in x direction");
        if (dimWorld > 1) {
            EWOMS_REGISTER_PARAM(TypeTag, Scalar, DomainSizeY,
                                 "The size of the domain in y direction");
            EWOMS_REGISTER_PARAM(TypeTag, unsigned, CellsY,
                                 "The number of intervalls in y direction");
        }
        if (dimWorld > 2) {
            EWOMS_REGISTER_PARAM(TypeTag, Scalar, DomainSizeZ,
                                 "The size of the domain in z direction");
            EWOMS_REGISTER_PARAM(TypeTag, unsigned, CellsZ,
                                 "The number of intervalls in z direction");
        }
    }

    /*!
     * \brief Create the Grid
     */
    SimplexGridVanguard(Simulator& simulator)
        : ParentType(simulator)
    {
        Dune::array<unsigned, dimWorld> cellRes;
        GlobalPosition upperRight;
        GlobalPosition lowerLeft;

        lowerLeft[0] = 0.0;
        upperRight[0] = EWOMS_GET_PARAM(TypeTag, Scalar, DomainSizeX);
        cellRes[0] = EWOMS_GET_PARAM(TypeTag, unsigned, CellsX);
        if (dimWorld > 1) {
            lowerLeft[1] = 0.0;
            upperRight[1] = EWOMS_GET_PARAM(TypeTag, Scalar, DomainSizeY);
            cellRes[1] = EWOMS_GET_PARAM(TypeTag, unsigned, CellsY);
        }
        if (dimWorld > 2) {
            lowerLeft[2] = 0.0;
            upperRight[2] = EWOMS_GET_PARAM(TypeTag, Scalar, DomainSizeZ);
            cellRes[2] = EWOMS_GET_PARAM(TypeTag, unsigned, CellsZ);
        }

        simplexGrid_ = Dune::StructuredGridFactory<Grid>::createSimplexGrid(lowerLeft,
                                                                            upperRight,
                                                                            cellRes);

        unsigned numRefinments = EWOMS_GET_PARAM(TypeTag, unsigned, GridGlobalRefinements);
        simplexGrid_->globalRefine(numRefinments);

        this->finalizeInit_();
    }

    /*!
     * \brief Returns a reference to the grid.
     */
    Grid& grid()
    { return simplexGrid_; }

    /*!
     * \brief Returns a reference to the grid.
     */
    const Grid& grid() const
    { return *simplexGrid_; }

private:
    GridPointer simplexGrid_;
};
} // namespace Opm

#endif
