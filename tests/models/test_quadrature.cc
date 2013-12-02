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
 * \brief A test for numerical integration using the vertex-centered finite
 * volume geometries.
 */
#include "config.h"

#include <dune/common/exceptions.hh>
#include <dune/grid/yaspgrid.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>

#if HAVE_ALUGRID
#include <dune/grid/alugrid.hh>
#endif // HAVE_ALUGRID

#include <dune/common/version.hh>
#include <dune/geometry/quadraturerules.hh>

#include <ewoms/disc/vcfv/vcfvelementgeometry.hh>

const unsigned dim = 3;
typedef double Scalar;
typedef Ewoms::QuadrialteralQuadratureGeometry<Scalar, dim> QuadratureGeom;
typedef QuadratureGeom::LocalPosition LocalPosition;
typedef QuadratureGeom::GlobalPosition GlobalPosition;

GlobalPosition::field_type f(const GlobalPosition &pos)
{
    GlobalPosition::field_type result = 1;
    for (int i = 0; i < GlobalPosition::dimension; ++i)
        result *= pos[i];
    return result;
}

void testIdenityMapping()
{
    QuadratureGeom foo;

    Scalar corners[][3] = { { 0, 0, 0 },
                            { 1, 0, 0 },
                            { 0, 1, 0 },
                            { 1, 1, 0 },
                            { 0, 0, 1 },
                            { 1, 0, 1 },
                            { 0, 1, 1 },
                            { 1, 1, 1 } };
    foo.setCorners(corners, 8);

    std::cout << "testing identity mapping...\n";
    int n = 100;
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            for (int k = 0; k < n; ++k) {
                LocalPosition localPos;

                localPos[0] = Scalar(i) / (n - 1);
                localPos[1] = Scalar(j) / (n - 1);
                localPos[2] = Scalar(k) / (n - 1);

                GlobalPosition globalPos = foo.global(localPos);

                GlobalPosition diff(localPos);
                diff -= globalPos;
                assert(diff.two_norm() < 1e-10);
            }
        }
    }
}

template <class Grid>
void writeTetrahedronSubControlVolumes(const Grid &grid)
{
#if HAVE_ALUGRID
    typedef typename Grid::LeafGridView GridView;

    typedef Dune::ALUGrid<dim, dim, Dune::cube, Dune::nonconforming> Grid2;
    typedef typename Grid2::LeafGridView GridView2;
    typedef Dune::GridFactory<Grid2> GridFactory2;

    // instanciate a ElementGeometry
    typedef Ewoms::VcfvElementGeometry<Scalar, GridView> FvElementGeometry;
    FvElementGeometry fvElemGeom;

    GridFactory2 gf2;
    const auto &gridView = grid.leafView();
    auto eIt = gridView.template begin<0>();
    const auto &eEndIt = gridView.template end<0>();
    for (; eIt != eEndIt; ++eIt) {
        fvElemGeom.update(gridView, *eIt);
        for (int scvIdx = 0; scvIdx < fvElemGeom.numVertices; ++scvIdx) {
            const auto &scvLocalGeom
                = *(fvElemGeom.subContVol[scvIdx].localGeometry);

            for (int i = 0; i < scvLocalGeom.numCorners; ++i) {
                GlobalPosition pos(
                    eIt->geometry().global(scvLocalGeom.corner(i)));
                gf2.insertVertex(pos);
            }
        }
    }

    int cornerOffset = 0;
    eIt = gridView.template begin<0>();
    for (; eIt != eEndIt; ++eIt) {
        fvElemGeom.update(gridView, *eIt);
        for (int scvIdx = 0; scvIdx < fvElemGeom.numVertices; ++scvIdx) {
            const auto &scvLocalGeom
                = *fvElemGeom.subContVol[scvIdx].localGeometry;

            std::vector<unsigned int> vertexIndices;
            for (int i = 0; i < scvLocalGeom.numCorners; ++i) {
                vertexIndices.push_back(cornerOffset);
                ++cornerOffset;
            }

            gf2.insertElement(Dune::GeometryType(Dune::GeometryType::cube, dim),
                              vertexIndices);
        }
    }

    const auto &grid2 = *gf2.createGrid();
    typedef Dune::VTKWriter<GridView2> VtkWriter;
    VtkWriter writer(grid2.leafView(), Dune::VTK::conforming);
    writer.write("tetrahedron-scvs", Dune::VTK::ascii);
#endif // HAVE_ALUGRID
}

void testTetrahedron()
{
#if HAVE_ALUGRID
    typedef Dune::ALUGrid<dim, dim, Dune::simplex, Dune::nonconforming> Grid;
    typedef Dune::GridFactory<Grid> GridFactory;
    GridFactory gf;
    Scalar corners[][3] = { { 0, 0, 0 }, { 1, 0, 0 }, { 0, 1, 0 }, { 0, 0, 1 } };

    for (unsigned i = 0; i < sizeof(corners) / sizeof(corners[0]); ++i) {
        GlobalPosition pos;
        for (unsigned j = 0; j < dim; ++j)
            pos[j] = corners[i][j];
        gf.insertVertex(pos);
    }
    std::vector<unsigned int> v = { 0, 1, 2, 3 };
    gf.insertElement(Dune::GeometryType(Dune::GeometryType::simplex, dim), v);
    auto *grid = gf.createGrid();

    // write the sub-control volumes to a VTK file.
    writeTetrahedronSubControlVolumes(*grid);

    delete grid;
#endif // HAVE_ALUGRID
}

template <class Grid>
void writeCubeSubControlVolumes(const Grid &grid)
{
#if HAVE_ALUGRID
    typedef typename Grid::LeafGridView GridView;

    typedef Dune::ALUGrid<dim, dim, Dune::cube, Dune::nonconforming> Grid2;
    typedef typename Grid2::LeafGridView GridView2;
    typedef Dune::GridFactory<Grid2> GridFactory2;

    // instanciate a ElementGeometry
    typedef Ewoms::VcfvElementGeometry<Scalar, GridView> FvElementGeometry;
    FvElementGeometry fvElemGeom;

    GridFactory2 gf2;
    const auto &gridView = grid.leafView();
    auto eIt = gridView.template begin<0>();
    const auto &eEndIt = gridView.template end<0>();
    for (; eIt != eEndIt; ++eIt) {
        fvElemGeom.update(gridView, *eIt);
        for (int scvIdx = 0; scvIdx < fvElemGeom.numVertices; ++scvIdx) {
            const auto &scvLocalGeom
                = *(fvElemGeom.subContVol[scvIdx].localGeometry);

            for (int i = 0; i < scvLocalGeom.numCorners; ++i) {
                GlobalPosition pos(
                    eIt->geometry().global(scvLocalGeom.corner(i)));
                gf2.insertVertex(pos);
            }
        }
    }

    int cornerOffset = 0;
    eIt = gridView.template begin<0>();
    for (; eIt != eEndIt; ++eIt) {
        fvElemGeom.update(gridView, *eIt);
        for (int scvIdx = 0; scvIdx < fvElemGeom.numVertices; ++scvIdx) {
            const auto &scvLocalGeom
                = *fvElemGeom.subContVol[scvIdx].localGeometry;

            std::vector<unsigned int> vertexIndices;
            for (int i = 0; i < scvLocalGeom.numCorners; ++i) {
                vertexIndices.push_back(cornerOffset);
                ++cornerOffset;
            }

            gf2.insertElement(Dune::GeometryType(Dune::GeometryType::cube, dim),
                              vertexIndices);
        }
    }

    const auto &grid2 = *gf2.createGrid();
    typedef Dune::VTKWriter<GridView2> VtkWriter;
    VtkWriter writer(grid2.leafView(), Dune::VTK::conforming);
    writer.write("cube-scvs", Dune::VTK::ascii);
#endif // HAVE_ALUGRID
}

void testCube()
{
#if HAVE_ALUGRID
    typedef Dune::ALUGrid<dim, dim, Dune::cube, Dune::nonconforming> Grid;
    typedef Dune::GridFactory<Grid> GridFactory;
    GridFactory gf;
    Scalar corners[][3] = { { 0, 0, 0 },
                            { 1, 0, 0 },
                            { 0, 2, 0 },
                            { 3, 3, 0 },
                            { 0, 0, 4 },
                            { 5, 0, 5 },
                            { 0, 6, 6 },
                            { 7, 7, 7 }, };

    for (unsigned i = 0; i < sizeof(corners) / sizeof(corners[0]); ++i) {
        GlobalPosition pos;
        for (unsigned j = 0; j < dim; ++j)
            pos[j] = corners[i][j];
        gf.insertVertex(pos);
    }
    std::vector<unsigned int> v = { 0, 1, 2, 3, 4, 5, 6, 7 };
    gf.insertElement(Dune::GeometryType(Dune::GeometryType::cube, dim), v);
    auto *grid = gf.createGrid();

    // write the sub-control volumes to a VTK file.
    writeCubeSubControlVolumes(*grid);

    delete grid;
#endif // HAVE_ALUGRID
}

void testQuadrature()
{
    std::cout << "testing SCV quadrature...\n";

#if DUNE_VERSION_NEWER(DUNE_COMMON, 2, 3)
    std::bitset<dim> isPeriodic(false);
    std::array<int, dim> cellRes;
#else
    Dune::FieldVector<bool, dim> isPeriodic(false);
    Dune::FieldVector<int, dim> cellRes;
#endif
    std::fill(cellRes.begin(), cellRes.end(), 10);

    GlobalPosition upperRight(1.0);

    typedef Dune::YaspGrid<dim> Grid;
    typedef Grid::LeafGridView GridView;
    Grid grid(
#ifdef HAVE_MPI
        Dune::MPIHelper::getCommunicator(),
#endif
        upperRight,     // upper right
        cellRes,        // number of cells
        isPeriodic, 0); // overlap

    // compute approximate integral
    auto gridView = grid.leafView();
    auto eIt = gridView.begin<0>();
    const auto eEndIt = gridView.end<0>();
    Scalar result = 0;
    // instanciate a FvElementGeometry
    typedef Ewoms::VcfvElementGeometry<Scalar, GridView> FvElementGeometry;
    FvElementGeometry fvElemGeom;
    for (; eIt != eEndIt; ++eIt) {
        const auto &elemGeom = eIt->geometry();

        fvElemGeom.update(gridView, *eIt);

        // loop over all sub-control volumes
        for (int scvIdx = 0; scvIdx < fvElemGeom.numVertices; ++scvIdx) {
            const auto &scvLocalGeom
                = *fvElemGeom.subContVol[scvIdx].localGeometry;

            Dune::GeometryType geomType = scvLocalGeom.type();
            static const int quadratureOrder = 2;
            const auto &rule
                = Dune::QuadratureRules<Scalar, dim>::rule(geomType,
                                                           quadratureOrder);

            // integrate f over the sub-control volume
            for (auto it = rule.begin(); it != rule.end(); ++it) {
                auto posScvLocal = it->position();
                auto posElemLocal = scvLocalGeom.global(posScvLocal);
                auto posGlobal = elemGeom.global(posScvLocal);

                Scalar fval = f(posGlobal);
                Scalar weight = it->weight();
                Scalar detjac = scvLocalGeom.integrationElement(posScvLocal)
                                * elemGeom.integrationElement(posElemLocal);

                result += fval * weight * detjac;
            }
        }
    }

    std::cout << "result: " << result
              << " (expected value: " << 1.0 / (1 << dim) << ")\n";
}

int main(int argc, char **argv)
{
    // initialize MPI, finalize is done automatically on exit
    Dune::MPIHelper::instance(argc, argv);

    testIdenityMapping();
// test the quadrature in a tetrahedron. since the CLang compiler
// prior to version 3.2 generates incorrect code here, we do not
// do it if the compiler is clang 3.1 or older.
#if !__clang__ || __clang_major__ > 3                                          \
    || (__clang_major__ == 3 && __clang_minor__ >= 3)
    testTetrahedron();
    testCube();
#endif
    testQuadrature();

    return 0;
}
