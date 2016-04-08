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

/*
SparseTableView cell2Faces(const UnstructuredGrid& grid)
{
    return SparseTableView(grid.cell_faces, grid.cell_facepos, numCells(grid));
}
*/

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

#ifdef HAVE_OPM_GRID

namespace AutoDiffGrid
{

ADFaceCellTraits<Dune::CpGrid>::Type
faceCellsToEigen(const Dune::CpGrid& grid)
{
    return Dune::cpgrid::FaceCellsContainerProxy(&grid);
}

Eigen::Array<double, Eigen::Dynamic, 1>
cellCentroidsZToEigen(const Dune::CpGrid& grid)
{
    // Create an Eigen array of appropriate size
    int rows=numCells(grid);
    Eigen::Array<double, Eigen::Dynamic, 1> array(rows);
    // Fill it with the z coordinate of the cell centroids.
    for (int i=0; i<rows; ++i)
        array[i]=Opm::UgGridHelpers::cellCentroid(grid, i)[2];
    return array;
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

}       // end namespace AutoDiffGrid
#endif  // HAVE_OPM_GRID
}       // end namespace Opm
