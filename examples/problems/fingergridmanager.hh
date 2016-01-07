// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*
  Copyright (C) 2012-2013 by Andreas Lauser

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
 * \copydoc Ewoms::FingerGridManager
 */
#ifndef EWOMS_FINGER_GRID_MANAGER_HH
#define EWOMS_FINGER_GRID_MANAGER_HH

#include <ewoms/io/basegridmanager.hh>
#include <ewoms/common/propertysystem.hh>
#include <ewoms/common/parametersystem.hh>

#if HAVE_DUNE_ALUGRID
// we cannot use SFC reordering since it messes up the test case logic
#define DISABLE_ALUGRID_SFC_ORDERING
#include <dune/alugrid/grid.hh>
#include <dune/alugrid/dgf.hh>
#else
#include <dune/grid/alugrid.hh>
#endif

#include <dune/grid/yaspgrid.hh>
#include <dune/grid/io/file/dgfparser/dgfyasp.hh>

#include <dune/common/fvector.hh>
#include <dune/common/version.hh>

#include <vector>
#include <memory>

// some hacky defines for the grid manager
#define FINGER_DIM 2
#define FINGER_CUBES 1

namespace Ewoms {
template <class TypeTag>
class FingerGridManager;

template <class TypeTag>
class FingerProblem;

namespace Properties {
// declare the properties required by the for the finger grid manager
NEW_TYPE_TAG(FingerGridManager);

NEW_PROP_TAG(Grid);
NEW_PROP_TAG(Scalar);

NEW_PROP_TAG(DomainSizeX);
NEW_PROP_TAG(DomainSizeY);
NEW_PROP_TAG(DomainSizeZ);

NEW_PROP_TAG(CellsX);
NEW_PROP_TAG(CellsY);
NEW_PROP_TAG(CellsZ);

NEW_PROP_TAG(GridGlobalRefinements);

#if HAVE_DUNE_ALUGRID
SET_TYPE_PROP(FingerGridManager, Grid, Dune::ALUGrid<FINGER_DIM, FINGER_DIM, Dune::cube, Dune::nonconforming>);
#else
SET_TYPE_PROP(FingerGridManager, Grid, Dune::YaspGrid<FINGER_DIM> );
#endif
SET_TYPE_PROP(FingerGridManager, GridManager, Ewoms::FingerGridManager<TypeTag>);

} // namespace Properties

/*!
 * \brief Helper class for grid instantiation of the finger problem.
 */
template <class TypeTag>
class FingerGridManager : public BaseGridManager<TypeTag>
{
    typedef BaseGridManager<TypeTag> ParentType;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, Simulator) Simulator;
    typedef typename GET_PROP_TYPE(TypeTag, Grid) Grid;

    typedef std::unique_ptr<Grid> GridPointer;

    enum { dim = FINGER_DIM };

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
     * \brief Create the grid for the finger problem
     */
    FingerGridManager(Simulator &simulator)
        : ParentType(simulator)
    {
        Dune::FieldVector<int, dim> cellRes;
        Dune::FieldVector<Scalar, dim> upperRight;
        Dune::FieldVector<Scalar, dim> lowerLeft;

        lowerLeft = 0.0;
        upperRight[0] = EWOMS_GET_PARAM(TypeTag, Scalar, DomainSizeX);
        upperRight[1] = EWOMS_GET_PARAM(TypeTag, Scalar, DomainSizeY);

        cellRes[0] = EWOMS_GET_PARAM(TypeTag, int, CellsX);
        cellRes[1] = EWOMS_GET_PARAM(TypeTag, int, CellsY);
        if (dim == 3) {
            upperRight[2] = EWOMS_GET_PARAM(TypeTag, Scalar, DomainSizeZ);
            cellRes[2] = EWOMS_GET_PARAM(TypeTag, int, CellsZ);
        }

        // create a DGF input stream containing an interval block
        std::stringstream dgfFile;
        dgfFile << "DGF"      << std::endl;
        dgfFile << "Interval" << std::endl;
        dgfFile << lowerLeft  << std::endl;
        dgfFile << upperRight << std::endl;
        dgfFile << cellRes    << std::endl;
        dgfFile << "#" << std::endl;
        dgfFile << "GridParameter" << std::endl;
        dgfFile << "overlap 1"     << std::endl;
        dgfFile << "#" << std::endl;
        dgfFile << "Simplex" << std::endl;
        dgfFile << "#" << std::endl;

        gridPtr_.reset( Dune::GridPtr< Grid >(dgfFile).release() );

        unsigned numRefinments = EWOMS_GET_PARAM(TypeTag, unsigned, GridGlobalRefinements);
        gridPtr_->globalRefine(numRefinments);

        this->finalizeInit_();
    }

    /*!
     * \brief Return a reference to the grid.
     */
    Grid& grid()
    { return *gridPtr_; }

    /*!
     * \brief Return a reference to the grid.
     */
    const Grid& grid() const
    { return *gridPtr_; }

private:
    GridPointer gridPtr_;
};

} // namespace Ewoms

#endif
