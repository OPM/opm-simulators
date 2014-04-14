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
 * \copydoc Ewoms::FingerGridCreator
 */
#ifndef EWOMS_FINGER_GRID_CREATOR_HH
#define EWOMS_FINGER_GRID_CREATOR_HH

#include <ewoms/parallel/mpihelper.hh>
#include <ewoms/io/basegridcreator.hh>
#include <opm/core/utility/PropertySystem.hpp>
#include <ewoms/common/parametersystem.hh>

#include <dune/grid/alugrid.hh>
#include <dune/common/fvector.hh>
#include <dune/common/version.hh>

#include <vector>

// some hacky defines for the grid creator
#define FINGER_DIM 2
#define FINGER_CUBES 1

namespace Ewoms {
template <class TypeTag>
class FingerGridCreator;

template <class TypeTag>
class FingerProblem;
} // namespace Ewoms

namespace Opm {
namespace Properties {
// declare the properties required by the for the finger grid creator
NEW_TYPE_TAG(FingerGridCreator);

NEW_PROP_TAG(Grid);
NEW_PROP_TAG(Scalar);

NEW_PROP_TAG(DomainSizeX);
NEW_PROP_TAG(DomainSizeY);
NEW_PROP_TAG(DomainSizeZ);

NEW_PROP_TAG(CellsX);
NEW_PROP_TAG(CellsY);
NEW_PROP_TAG(CellsZ);

NEW_PROP_TAG(GridGlobalRefinements);

SET_TYPE_PROP(FingerGridCreator, Grid, Dune::ALUGrid<FINGER_DIM, FINGER_DIM, Dune::cube, Dune::nonconforming>);
SET_TYPE_PROP(FingerGridCreator, GridCreator, Ewoms::FingerGridCreator<TypeTag>);

}} // namespace Opm, Properties

namespace Ewoms {
/*!
 * \brief Helper class for grid instantiation of the finger problem.
 */
template <class TypeTag>
class FingerGridCreator : public BaseGridCreator<TypeTag>
{
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;

    enum { dim = FINGER_DIM };

public:
    typedef typename GET_PROP_TYPE(TypeTag, Grid) Grid;

    /*!
     * \brief Register all run-time parameters for the grid creator.
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
    static void makeGrid()
    {
        grid_ = new Grid;

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

        unsigned numRefinments
            = EWOMS_GET_PARAM(TypeTag, unsigned, GridGlobalRefinements);

        Dune::GridFactory<Grid> factory(grid_);

        if (dim == 3) {
            Dune::FieldVector<double, dim> pos;
            for (int k = 0; k <= cellRes[0]; k++) {
                pos[2] = upperRight[2] * double(k) / cellRes[2];

                for (int j = 0; j <= cellRes[1]; j++) {
                    pos[1] = upperRight[1] * double(j) / cellRes[1];

                    for (int i = 0; i <= cellRes[0]; i++) {
                        pos[0] = upperRight[0] * double(i) / cellRes[0];
                        factory.insertVertex(pos);
                    }
                }
            }
        }
        else {
            assert(dim == 2);
            Dune::FieldVector<double, dim> pos;
            for (int j = 0; j <= cellRes[1]; j++) {
                pos[1] = upperRight[1] * double(j) / cellRes[1];

                for (int i = 0; i <= cellRes[0]; i++) {
                    pos[0] = upperRight[0] * double(i) / cellRes[0];
                    factory.insertVertex(pos);
                }
            }
        }

        for (int i = 0; i < cellRes[0]; ++i) {
            for (int j = 0; j < cellRes[1]; ++j) {
#if FINGER_CUBES
                std::vector<unsigned int> v(1 << dim);
#else
                std::vector<unsigned int> v(dim + 1);
#endif
                if (dim == 3) {
                    int m = cellRes[0] + 1;
                    int n = cellRes[1] + 1;
                    for (int k = 0; k < cellRes[2]; ++k) {
                        int i0 = k * m * n + j * m + i;
                        int i1 = k * m * n + j * m + (i + 1);
                        int i2 = k * m * n + (j + 1) * m + i;
                        int i3 = k * m * n + (j + 1) * m + (i + 1);
                        int i4 = (k + 1) * m * n + j * m + i;
                        int i5 = (k + 1) * m * n + j * m + (i + 1);
                        int i6 = (k + 1) * m * n + (j + 1) * m + i;
                        int i7 = (k + 1) * m * n + (j + 1) * m + (i + 1);

#if FINGER_CUBES
                        v[0] = i0;
                        v[1] = i1;
                        v[2] = i2;
                        v[3] = i3;
                        v[4] = i4;
                        v[5] = i5;
                        v[6] = i6;
                        v[7] = i7;
                        factory.insertElement(Dune::GeometryType(Dune::GeometryType::cube, 3), v);

#else
                        v[0] = i0;
                        v[1] = i1;
                        v[2] = i2;
                        v[3] = i4;
                        factory.insertElement(Dune::GeometryType(Dune::GeometryType::simplex, 3), v);

                        v[0] = i4;
                        v[1] = i5;
                        v[2] = i6;
                        v[3] = i2;
                        factory.insertElement(Dune::GeometryType(Dune::GeometryType::simplex, 3), v);

                        v[0] = i2;
                        v[1] = i5;
                        v[2] = i4;
                        v[3] = i1;
                        factory.insertElement(Dune::GeometryType(Dune::GeometryType::simplex, 3), v);

                        v[0] = i2;
                        v[1] = i3;
                        v[2] = i7;
                        v[3] = i5;
                        factory.insertElement(Dune::GeometryType(Dune::GeometryType::simplex, 3), v);

                        v[0] = i5;
                        v[1] = i7;
                        v[2] = i6;
                        v[3] = i2;
                        factory.insertElement(Dune::GeometryType(Dune::GeometryType::simplex, 3), v);

                        v[0] = i1;
                        v[1] = i3;
                        v[2] = i5;
                        v[3] = i2;
                        factory.insertElement(Dune::GeometryType(Dune::GeometryType::simplex, 3), v);
#endif
                    }
                }
                else {
                    assert(dim == 2);

                    int m = cellRes[0] + 1;
                    int i0 = j * m + i;
                    int i1 = j * m + (i + 1);
                    int i2 = (j + 1) * m + i;
                    int i3 = (j + 1) * m + (i + 1);
#if FINGER_CUBES
                    v[0] = i0;
                    v[1] = i1;
                    v[2] = i2;
                    v[3] = i3;
                    factory.insertElement(Dune::GeometryType(Dune::GeometryType::cube, 2), v);
#else
                    v[0] = i0;
                    v[1] = i1;
                    v[2] = i2;
                    factory.insertElement(Dune::GeometryType(Dune::GeometryType::simplex, 2), v);

                    v[0] = i1;
                    v[1] = i3;
                    v[2] = i2;
                    factory.insertElement(Dune::GeometryType(Dune::GeometryType::simplex, 2), v);
#endif
                }
            }
        }

        grid_ = factory.createGrid();
        grid_->globalRefine(numRefinments);
    }

    /*!
     * \brief Return a reference to the grid.
     */
    static Grid* gridPointer()
    { return grid_; }

    /*!
     * \brief Distribute the grid (and attached data) over all
     *        processes.
     */
    static void loadBalance()
    { grid_->loadBalance(); }

    /*!
     * \brief Destroy the grid.
     *
     * This is required to guarantee that the grid is deleted before
     * MPI_Comm_free is called.
     */
    static void deleteGrid()
    { delete grid_; }

private:
    static Grid *grid_;
};

template <class TypeTag>
typename FingerGridCreator<TypeTag>::Grid *FingerGridCreator<TypeTag>::grid_;
} // namespace Ewoms

#endif
