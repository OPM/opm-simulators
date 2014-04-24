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
 * \copydoc Ewoms::LensGridManager
 */
#ifndef EWOMS_LENS_GRID_MANAGER_HH
#define EWOMS_LENS_GRID_MANAGER_HH

#include <ewoms/parallel/mpihelper.hh>
#include <ewoms/io/basegridmanager.hh>
#include <opm/core/utility/PropertySystem.hpp>
#include <ewoms/common/parametersystem.hh>

#include <dune/grid/yaspgrid.hh>
#include <dune/common/fvector.hh>
#include <dune/common/version.hh>

#include <vector>

namespace Ewoms {
template <class TypeTag>
class LensProblem;

template <class TypeTag>
class LensGridManager;
} // namespace Ewoms

namespace Opm {
namespace Properties {
NEW_TYPE_TAG(LensGridManager);

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

// set the Grid and GridManager properties
SET_TYPE_PROP(LensGridManager, Grid, Dune::YaspGrid<2>);
SET_TYPE_PROP(LensGridManager, GridManager, Ewoms::LensGridManager<TypeTag>);
}} // namespace Opm, Properties

namespace Ewoms {
/*!
 * \ingroup TestProblems
 *
 * \brief Helper class for grid instantiation of the lens problem.
 */
template <class TypeTag>
class LensGridManager : public BaseGridManager<TypeTag>
{
    typedef BaseGridManager<TypeTag> ParentType;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, Simulator) Simulator;
    typedef typename GET_PROP_TYPE(TypeTag, Grid) Grid;

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
    LensGridManager(Simulator &simulator)
        : ParentType(simulator)
    {
#if DUNE_VERSION_NEWER(DUNE_COMMON, 2, 3)
        std::bitset<dim> isPeriodic(false);
        std::array<int, dim> cellRes;
#else
        Dune::FieldVector<bool, dim> isPeriodic(false);
        Dune::FieldVector<int, dim> cellRes;
#endif

        Dune::FieldVector<Scalar, dim> upperRight;
        Dune::FieldVector<Scalar, dim> lowerLeft;

        grid_ = 0;

        lowerLeft[1] = 0.0;
        upperRight[0] = EWOMS_GET_PARAM(TypeTag, Scalar, DomainSizeX);
        upperRight[1] = EWOMS_GET_PARAM(TypeTag, Scalar, DomainSizeY);

        cellRes[0] = EWOMS_GET_PARAM(TypeTag, int, CellsX);
        cellRes[1] = EWOMS_GET_PARAM(TypeTag, int, CellsY);
        if (dim == 3) {
            upperRight[2] = EWOMS_GET_PARAM(TypeTag, Scalar, DomainSizeZ);
            cellRes[2] = EWOMS_GET_PARAM(TypeTag, int, CellsZ);
        }

        unsigned numRefinements
            = EWOMS_GET_PARAM(TypeTag, unsigned, GridGlobalRefinements);

        grid_ = new Dune::YaspGrid<dim>(
#ifdef HAVE_MPI
            /*mpiCommunicator=*/Dune::MPIHelper::getCommunicator(),
#endif
            /*upperRightCorner=*/upperRight,
            /*numCells=*/cellRes, isPeriodic,
            /*overlap=*/1);
        grid_->globalRefine(numRefinements);

        this->finalizeInit_();
    }

    /*!
     * \brief Destroy the grid
     */
    ~LensGridManager()
    { delete grid_; }

    /*!
     * \brief Return a pointer to the grid object.
     */
    Grid* gridPointer()
    { return grid_; }

    /*!
     * \brief Return a constant pointer to the grid object.
     */
    const Grid* gridPointer() const
    { return grid_; }

private:
    Grid *grid_;
};

} // namespace Ewoms

#endif
