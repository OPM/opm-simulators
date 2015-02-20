/*
  Copyright 2014, 2015 Dr. Markus Blatt - HPC-Simulation-Software & Services.
  Copyright 2014 Statoil AS
  Copyright 2015 NTNU

  This file is part of the Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "config.h"

#include <opm/autodiff/GridHelpers.hpp>
namespace Opm
{
namespace AutoDiffGrid
{

// Interface functions using Unstructured grid
/*    
int numCells(const UnstructuredGrid& grid)
{
    return grid.number_of_cells;
}

int numFaces(const UnstructuredGrid& grid)
{
    return grid.number_of_faces;
}
int dimensions(const UnstructuredGrid& grid)
{
    return grid.dimensions;
}
*/
Eigen::Array<int, Eigen::Dynamic, 2, Eigen::RowMajor>
faceCellsToEigen(const UnstructuredGrid& grid)
{
    typedef Eigen::Array<int, Eigen::Dynamic, 2, Eigen::RowMajor> TwoColInt;
    return Eigen::Map<TwoColInt>(grid.face_cells, grid.number_of_faces, 2);
}

Eigen::Array<double, Eigen::Dynamic, 1>
cellCentroidsZToEigen(const UnstructuredGrid& grid)
{
    return Eigen::Map<Eigen::Array<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> >
        (grid.cell_centroids, grid.number_of_cells, grid.dimensions).rightCols<1>();
}

const double*
cellCentroid(const UnstructuredGrid& grid, int cell_index)
{
    return grid.cell_centroids+(cell_index*grid.dimensions);
}

const double* faceCentroid(const UnstructuredGrid& grid, int face_index)
{
    return grid.face_centroids+(face_index*grid.dimensions);
}
/*
SparseTableView cell2Faces(const UnstructuredGrid& grid)
{
    return SparseTableView(grid.cell_faces, grid.cell_facepos, numCells(grid));
}
*/
double cellVolume(const UnstructuredGrid& grid, int cell_index)
{
    return grid.cell_volumes[cell_index];
}

const double* beginCellVolumes(const UnstructuredGrid& grid)
{
    return grid.cell_volumes;
}
const double* endCellVolumes(const UnstructuredGrid& grid)
{
    return grid.cell_volumes+numCells(grid);
}

void extractInternalFaces(const UnstructuredGrid& grid,
                          Eigen::Array<int, Eigen::Dynamic, 1>& internal_faces,
                          Eigen::Array<int, Eigen::Dynamic, 2, Eigen::RowMajor>& nbi)
{
    typedef Eigen::Array<bool, Eigen::Dynamic, 1> OneColBool;
    typedef Eigen::Array<int, Eigen::Dynamic, 2, Eigen::RowMajor> TwoColInt;
    typedef Eigen::Array<bool, Eigen::Dynamic, 2, Eigen::RowMajor> TwoColBool;
    TwoColInt nb = faceCellsToEigen(grid);
    // std::cout << "nb = \n" << nb << std::endl;
    // Extracts the internal faces of the grid.
    // These are stored in internal_faces.
    TwoColBool nbib = nb >= 0;
    OneColBool ifaces = nbib.rowwise().all();
    const int num_internal = ifaces.cast<int>().sum();
    // std::cout << num_internal << " internal faces." << std::endl;
    nbi.resize(num_internal, 2);
    internal_faces.resize(num_internal);
    int fi = 0;
    int nf = numFaces(grid);

    for (int f = 0; f < nf; ++f) {
        if (ifaces[f]) {
            internal_faces[fi] = f;
            nbi.row(fi) = nb.row(f);
            ++fi;
        }
    }
}
} // end namespace AutoDiffGrid

#ifdef HAVE_DUNE_CORNERPOINT
// Interface functions using CpGrid

namespace UgGridHelpers
{

int numCells(const Dune::CpGrid& grid)
{
    return grid.numCells();
}

int numFaces(const  Dune::CpGrid& grid)
{
    return grid.numFaces();
}

int dimensions(const Dune::CpGrid&)
{
    return Dune::CpGrid::dimension;
}

int numCellFaces(const Dune::CpGrid& grid)
{
    return grid.numCellFaces();    
}

const int* cartDims(const Dune::CpGrid& grid)
{
    return &(grid.logicalCartesianSize()[0]);
}

const int*  globalCell(const Dune::CpGrid& grid)
{
    return &(grid.globalCell()[0]);
}

CellCentroidTraits<Dune::CpGrid>::IteratorType
beginCellCentroids(const Dune::CpGrid& grid)
{
    return CellCentroidTraits<Dune::CpGrid>::IteratorType(grid, 0);
}

double cellCentroidCoordinate(const Dune::CpGrid& grid, int cell_index,
                              int coordinate)
{
    return grid.cellCentroid(cell_index)[coordinate];
}

FaceCentroidTraits<Dune::CpGrid>::IteratorType
beginFaceCentroids(const Dune::CpGrid& grid)
{
    return FaceCentroidTraits<Dune::CpGrid>::IteratorType(grid, 0);
}

FaceCentroidTraits<Dune::CpGrid>::ValueType
faceCentroid(const Dune::CpGrid& grid, int face_index)
{
    return grid.faceCentroid(face_index);
}

Opm::AutoDiffGrid::Cell2FacesContainer cell2Faces(const Dune::CpGrid& grid)
{
    return Opm::AutoDiffGrid::Cell2FacesContainer(&grid);
}

FaceCellTraits<Dune::CpGrid>::Type
faceCells(const Dune::CpGrid& grid)
{
    return Opm::AutoDiffGrid::FaceCellsContainerProxy(&grid);
}

Face2VerticesTraits<Dune::CpGrid>::Type
face2Vertices(const Dune::CpGrid& grid)
{
    return Opm::AutoDiffGrid::FaceVerticesContainerProxy(&grid);
}

const double* vertexCoordinates(const Dune::CpGrid& grid, int index)
{
    return &(grid.vertexPosition(index)[0]);
}

const double* faceNormal(const Dune::CpGrid& grid, int face_index)
{
    return &(grid.faceNormal(face_index)[0]);
}

double faceArea(const Dune::CpGrid& grid, int face_index)
{
    return grid.faceArea(face_index);
}

int faceTag(const Dune::CpGrid& grid,
            const Opm::AutoDiffGrid::Cell2FacesRow::iterator& cell_face)
{
    return grid.faceTag(cell_face);
}
} // end namespace UgGridHelpers

namespace AutoDiffGrid
{

ADFaceCellTraits<Dune::CpGrid>::Type
faceCellsToEigen(const Dune::CpGrid& grid)
{
    return Opm::AutoDiffGrid::FaceCellsContainerProxy(&grid);
}

Eigen::Array<double, Eigen::Dynamic, 1>
cellCentroidsZToEigen(const Dune::CpGrid& grid)
{
    // Create an Eigen array of appropriate size
    int rows=numCells(grid);
    Eigen::Array<double, Eigen::Dynamic, 1> array(rows);
    // Fill it with the z coordinate of the cell centroids.
    for (int i=0; i<rows; ++i)
        array[i]=cellCentroid(grid, i)[2];
    return array;
}

const double* cellCentroid(const Dune::CpGrid& grid, int cell_index)
{
    return &(grid.cellCentroid(cell_index)[0]);
}

const double* faceCentroid(const Dune::CpGrid& grid, int face_index)
{
    return &(grid.faceCentroid(face_index)[0]);
}

double cellVolume(const  Dune::CpGrid& grid, int cell_index)
{
    return grid.cellVolume(cell_index);
}

void extractInternalFaces(const Dune::CpGrid& grid,
                          Eigen::Array<int, Eigen::Dynamic, 1>& internal_faces,
                          Eigen::Array<int, Eigen::Dynamic, 2, Eigen::RowMajor>& nbi)
{
    // Extracts the internal faces of the grid.
    // These are stored in internal_faces.
    int nf=numFaces(grid);
    int num_internal=0;
    for(int f=0; f<nf; ++f)
    {
        if(grid.faceCell(f, 0)<0 || grid.faceCell(f, 1)<0)
            continue;
        ++num_internal;
    }
    // std::cout << num_internal << " internal faces." << std::endl;
    nbi.resize(num_internal, 2);
    internal_faces.resize(num_internal);
    int fi = 0;

    for (int f = 0; f < nf; ++f) {
        if(grid.faceCell(f, 0)>=0 && grid.faceCell(f, 1)>=0) {
            internal_faces[fi] = f;
            nbi(fi,0) = grid.faceCell(f, 0);
            nbi(fi,1) = grid.faceCell(f, 1);
            ++fi;
        }
    }
}


CellVolumeIterator beginCellVolumes(const Dune::CpGrid& grid)
{
    return CellVolumeIterator(grid, 0);
}

CellVolumeIterator endCellVolumes(const Dune::CpGrid& grid)
{
    return CellVolumeIterator(grid, numCells(grid));
}
}       // end namespace AutoDiffGrid
#endif  // HAVE_DUNE_CORNERPOINT
}       // end namespace Opm
