// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*
  Copyright (C) 2012-2013 by Andreas Lauser
  Copyright (C) 2016      IRIS AS

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
 * \copydoc Ewoms::StructuredGridManager
 */
#ifndef EWOMS_STRUCTURED_GRID_MANAGER_HH
#define EWOMS_STRUCTURED_GRID_MANAGER_HH

#include <ewoms/io/basegridmanager.hh>
#include <ewoms/common/propertysystem.hh>
#include <ewoms/common/parametersystem.hh>

#include <dune/grid/yaspgrid.hh>
#include <dune/grid/io/file/dgfparser/dgfyasp.hh>

//#define TESTS_USE_ALUGRID 1
#if TESTS_USE_ALUGRID
#include <dune/alugrid/grid.hh>
#include <dune/alugrid/dgf.hh>
#endif

#include <dune/common/fvector.hh>
#include <dune/common/version.hh>

#include <vector>
#include <memory>

namespace Ewoms {

template <class TypeTag>
class StructuredGridManager;

namespace Properties {
NEW_TYPE_TAG(StructuredGridManager);

// declare the properties required by the for the lens grid manager
NEW_PROP_TAG(Grid);
NEW_PROP_TAG(Scalar);

NEW_PROP_TAG(DomainSizeX);
NEW_PROP_TAG(DomainSizeY);
NEW_PROP_TAG(DomainSizeZ);

NEW_PROP_TAG(CellsX);
NEW_PROP_TAG(CellsY);
NEW_PROP_TAG(CellsZ);

NEW_PROP_TAG(GridGlobalRefinements);

// GRIDDIM is only set by the finger problem
#ifndef GRIDDIM
static const int dim = 2;
#else
static const int dim = GRIDDIM;
#endif

// set the Grid and GridManager properties
#if TESTS_USE_ALUGRID
SET_TYPE_PROP(StructuredGridManager, Grid, Dune::ALUGrid< dim, dim, Dune::cube, Dune::nonconforming >);
#else
SET_TYPE_PROP(StructuredGridManager, Grid, Dune::YaspGrid< dim >);
#endif

SET_TYPE_PROP(StructuredGridManager, GridManager, Ewoms::StructuredGridManager<TypeTag>);
} // namespace Properties

/*!
 * \ingroup TestProblems
 *
 * \brief Helper class for grid instantiation of the lens problem.
 */
template <class TypeTag>
class StructuredGridManager : public BaseGridManager<TypeTag>
{
    typedef BaseGridManager<TypeTag> ParentType;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, Simulator) Simulator;
    typedef typename GET_PROP_TYPE(TypeTag, Grid) Grid;

    typedef std::unique_ptr<Grid> GridPointer;

    static const int dim = Grid::dimension;

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
        EWOMS_REGISTER_PARAM(TypeTag, int, CellsX,
                             "The number of intervalls in x direction");
        if (dim > 1) {
            EWOMS_REGISTER_PARAM(TypeTag, Scalar, DomainSizeY,
                                 "The size of the domain in y direction");
            EWOMS_REGISTER_PARAM(TypeTag, int, CellsY,
                                 "The number of intervalls in y direction");
        }
        if (dim > 2) {
            EWOMS_REGISTER_PARAM(TypeTag, Scalar, DomainSizeZ,
                                 "The size of the domain in z direction");
            EWOMS_REGISTER_PARAM(TypeTag, int, CellsZ,
                                 "The number of intervalls in z direction");
        }
    }

    /*!
     * \brief Create the grid for the lens problem
     */
    StructuredGridManager(Simulator &simulator)
        : ParentType(simulator)
    {
        Dune::FieldVector<int, dim> cellRes;

        typedef double GridScalar;
        Dune::FieldVector<GridScalar, dim> upperRight;
        Dune::FieldVector<GridScalar, dim> lowerLeft( 0 );

        upperRight[0] = EWOMS_GET_PARAM(TypeTag, Scalar, DomainSizeX);
        upperRight[1] = EWOMS_GET_PARAM(TypeTag, Scalar, DomainSizeY);

        cellRes[0] = EWOMS_GET_PARAM(TypeTag, int, CellsX);
        cellRes[1] = EWOMS_GET_PARAM(TypeTag, int, CellsY);
        if (dim == 3) {
            upperRight[2] = EWOMS_GET_PARAM(TypeTag, Scalar, DomainSizeZ);
            cellRes[2] = EWOMS_GET_PARAM(TypeTag, int, CellsZ);
        }

        std::stringstream dgffile;
        dgffile << "DGF" << std::endl;
        dgffile << "INTERVAL" << std::endl;
        dgffile << lowerLeft  << std::endl;
        dgffile << upperRight << std::endl;
        dgffile << cellRes    << std::endl;
        dgffile << "#" << std::endl;
        dgffile << "GridParameter" << std::endl;
        dgffile << "overlap 1" << std::endl;
        dgffile << "#" << std::endl;
        dgffile << "Simplex" << std::endl;
        dgffile << "#" << std::endl;

        // use DGF parser to create a grid from interval block
        gridPtr_.reset( Dune::GridPtr< Grid >( dgffile ).release() );

        unsigned numRefinements = EWOMS_GET_PARAM(TypeTag, unsigned, GridGlobalRefinements);
        gridPtr_->globalRefine(numRefinements);

        this->finalizeInit_();
    }

    /*!
     * \brief Return a reference to the grid object.
     */
    Grid& grid()
    { return *gridPtr_; }

    /*!
     * \brief Return a constant reference to the grid object.
     */
    const Grid& grid() const
    { return *gridPtr_; }

private:
    GridPointer gridPtr_;
};

} // namespace Ewoms

#endif
