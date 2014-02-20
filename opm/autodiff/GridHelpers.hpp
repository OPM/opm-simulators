/*
  Copyright 2014 Dr. Markus Blatt - HPC-Simulation-Software & Services.
  Copyright 2014 Statoil AS

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
#ifndef OPM_GRIDHELPERS_HEADER_INCLUDED
#define OPM_GRIDHELPERS_HEADER_INCLUDED

#include <boost/range/iterator_range.hpp>
#include <opm/core/grid.h>
#include <opm/core/grid/GridHelpers.hpp>
#include <Eigen/Eigen>
#include <Eigen/Sparse>

#ifdef HAVE_DUNE_CORNERPOINT
#include <dune/grid/CpGrid.hpp>
#endif

namespace Opm
{

namespace AutoDiffGrid
{

using Opm::UgGridHelpers::SparseTableView;
using Opm::UgGridHelpers::numCells;
using Opm::UgGridHelpers::numFaces;
using Opm::UgGridHelpers::dimensions;
using Opm::UgGridHelpers::cartDims;
using Opm::UgGridHelpers::globalCell;
using Opm::UgGridHelpers::cell2Faces;
using Opm::UgGridHelpers::increment;
using Opm::UgGridHelpers::getCoordinate;

/// \brief Mapps a grid type to the corresponding face to cell mapping.
///
/// The value of the mapping is provided by the type Type.
template<class T>
struct ADFaceCellTraits
{
};

template<>
struct ADFaceCellTraits<UnstructuredGrid>
{
    typedef Eigen::Array<int, Eigen::Dynamic, 2, Eigen::RowMajor> Type;
};

/// \brief Get the face to cell mapping of a grid.
ADFaceCellTraits<UnstructuredGrid>::Type
faceCells(const UnstructuredGrid& grid);

/// \brief Get the z coordinates of the cell centroids of a grid.
Eigen::Array<double, Eigen::Dynamic, 1>
cellCentroidsZ(const UnstructuredGrid& grid);

/// \brief Get the centroid of a cell.
/// \param grid The grid whose cell centroid we query.
/// \param cell_index The index of the corresponding cell.
const double* cellCentroid(const UnstructuredGrid& grid, int cell_index);

/// \brief Get the cell centroid of a face.
/// \param grid The grid whose cell centroid we query.
/// \param face_index The index of the corresponding face.
const double* faceCentroid(const UnstructuredGrid& grid, int face_index);

/// \brief Mapping of the grid type to the type of the cell to faces mapping.
template<class T>
struct ADCell2FacesTraits
    : public Opm::UgGridHelpers::Cell2FacesTraits<T>
{
};

/// \brief Get the volume of a cell.
/// \param grid The grid the cell belongs to.
/// \param cell_index The index of the cell.
double cellVolume(const UnstructuredGrid& grid, int cell_index);

/// \brief The mapping of the grid type to type of the iterator over
/// the cell volumes.
///
/// The value of the mapping is stored in nested type IteratorType
/// \tparam T The type of the grid.
template<class T>
struct ADCellVolumesTraits
{
};

template<>
struct ADCellVolumesTraits<UnstructuredGrid>
{
    typedef const double* IteratorType;
};

/// \brief Get an iterator over the cell volumes of a grid positioned at the first cell.
const double* beginCellVolumes(const UnstructuredGrid& grid);

/// \brief Get an iterator over the cell volumes of a grid positioned after the last cell.
const double* endCellVolumes(const UnstructuredGrid& grid);

/// \brief extracts the internal faces of a grid.
/// \param[in] The grid whose internal faces we query.
/// \param[out] internal_faces The internal faces.
/// \param[out] nbi 
void extractInternalFaces(const UnstructuredGrid& grid,
                          Eigen::Array<int, Eigen::Dynamic, 1>& internal_faces,
                          Eigen::Array<int, Eigen::Dynamic, 2, Eigen::RowMajor>& nbi);

using Opm::UgGridHelpers::beginFaceCentroids;
using Opm::UgGridHelpers::beginCellCentroids;
}
}

#ifdef HAVE_DUNE_CORNERPOINT

#include <dune/common/iteratorfacades.hh>
namespace Opm
{

namespace AutoDiffGrid
{
/// \brief Get the number of cells of a grid.
int numCells(const Dune::CpGrid& grid);

/// \brief Get the number of faces of a grid.
int numFaces(const  Dune::CpGrid& grid);

/// \brief Get the dimensions of a grid
int dimensions(const Dune::CpGrid& grid);

/// \brief Get the cartesion dimension of the underlying structured grid.
const int* cartDims(const Dune::CpGrid& grid);

/// \brief Get the local to global index mapping.
///
/// The global index is the index of the active cell
/// in the underlying structured grid.
const int*  globalCell(const Dune::CpGrid&);

/// \brief A proxy class representing a row of FaceCellsContainer.
class FaceCellsProxy
{    
public:
    /// \brief Constructor.
    /// \param grid The grid whose face to cell mapping we represent.
    /// \param cell_index The index of the cell we repesent.
    FaceCellsProxy(const Dune::CpGrid* grid, int cell_index)
        : grid_(grid), cell_index_(cell_index)
    {}
    /// \brief Get the index of the cell associated with a local_index.
    int operator[](int local_index)
    {
        return grid_->faceCell(cell_index_, local_index);
    }
private:
    const Dune::CpGrid* grid_;
    int cell_index_;
};

/// \brief A class representing the face to cells mapping similar to the
/// way done in UnstructuredGrid.
class FaceCellsContainerProxy
{
public:
    /// \brief Constructor.
    /// \param grid The grid whose information we represent.
    FaceCellsContainerProxy(const Dune::CpGrid* grid)
        : grid_(grid)
    {}
    /// \brief Get the mapping for a cell.
    /// \param cell_index The index of the cell.
    FaceCellsProxy operator[](int cell_index) const
    {
        return FaceCellsProxy(grid_, cell_index);
    }
    /// \brief Get a face associated with a cell.
    /// \param cell_index The index of the cell.
    /// \param local_index The local index of the cell, either 0 or 1.
    /// \param The index of the face or -1 if it is not present because of
    /// a boundary.
    int operator()(int cell_index, int local_index) const
    {
        return grid_->faceCell(cell_index, local_index);
    }
private:
    const Dune::CpGrid* grid_;
};

template<>
struct ADFaceCellTraits<Dune::CpGrid>
{
    typedef FaceCellsContainerProxy Type;
};

/// \brief Get the face to cell mapping of a grid.
ADFaceCellTraits<Dune::CpGrid>::Type
faceCells(const Dune::CpGrid& grid);

/// \brief Get the z coordinates of the cell centroids of a grid.
Eigen::Array<double, Eigen::Dynamic, 1>
cellCentroidsZ(const  Dune::CpGrid& grid);

/// \brief Get the centroid of a cell.
/// \param grid The grid whose cell centroid we query.
/// \param cell_index The index of the corresponding cell.
const double* cellCentroid(const Dune::CpGrid& grid, int cell_index);

/// \brief Get the cell centroid of a face.
/// \param grid The grid whose cell centroid we query.
/// \param face_index The index of the corresponding face.
const double* faceCentroid(const Dune::CpGrid& grid, int face_index);

class Cell2FacesRow
{
public:
    class iterator
    {
    public:
        iterator(const Dune::cpgrid::OrientedEntityTable<0,1>::row_type* row,
                 int index)
            : row_(row), index_(index)
        {}

        iterator operator++()
        {
            ++index_;
            return *this;
        }
        iterator operator++(int)
        {
            iterator ret=*this;
            ++index_;
            return ret;
        }
        int operator*()
        {
            return row_->operator[](index_).index();
        }
    private:
        const Dune::cpgrid::OrientedEntityTable<0,1>::row_type* row_;
        int index_;
    };
    
    typedef iterator const_iterator;

    Cell2FacesRow(const Dune::cpgrid::OrientedEntityTable<0,1>::row_type row)
        : row_(row)
    {}

    const_iterator begin() const
    {
        return const_iterator(&row_, 0);
    }

    const_iterator end() const
    {
        return const_iterator(&row_, row_.size());
    }
    
private:
    const Dune::cpgrid::OrientedEntityTable<0,1>::row_type row_;
};

class Cell2FacesContainer
{
public:
    Cell2FacesContainer(const Dune::CpGrid* grid)
        : grid_(grid)
    {};
    
    Cell2FacesRow operator[](int cell_index)
    {
        return Cell2FacesRow(grid_->cellFaceRow(cell_index));
    }
    
private:
    const Dune::CpGrid* grid_;
};

template<>
struct ADCell2FacesTraits<Dune::CpGrid>
{
    typedef Cell2FacesContainer Type;
};

/// \brief Get the cell to faces mapping of a grid.
Cell2FacesContainer cell2Faces(const Dune::CpGrid& grid);

/// \brief Get the volume of a cell.
/// \param grid The grid the cell belongs to.
/// \param cell_index The index of the cell.
double cellVolume(const  Dune::CpGrid& grid, int cell_index);

/// \brief An iterator over the cell volumes.
class CellVolumeIterator
    : public Dune::RandomAccessIteratorFacade<CellVolumeIterator, double, double, int>
{
public:
    /// \brief Creates an iterator.
    /// \param grid The grid the iterator belongs to.
    /// \param cell_index The position of the iterator.
    CellVolumeIterator(const  Dune::CpGrid& grid, int cell_index)
        : grid_(&grid), cell_index_(cell_index)
    {}

    double dereference()
    {
        return grid_->cellVolume(cell_index_);
    }
    void increment()
    {
        ++cell_index_;
    }
    double elementAt(int n)
    {
        return grid_->cellVolume(n);
    }
    void adavance(int n)
    {
        cell_index_+=n;
    }
    void decrement()
    {
        --cell_index_;
    }
    int distanceTo(const CellVolumeIterator& o)
    {
        return o.cell_index_-cell_index_;
    }
    bool equals(const CellVolumeIterator& o)
    {
        return o.grid_==grid_ && o.cell_index_==cell_index_;
    }
    
private:
    const Dune::CpGrid* grid_;
    int cell_index_;
};

template<>
struct ADCellVolumesTraits<Dune::CpGrid>
{
    typedef CellVolumeIterator IteratorType;
};

/// \brief Get an iterator over the cell volumes of a grid positioned at the first cell.
CellVolumeIterator beginCellVolumes(const Dune::CpGrid& grid);

/// \brief Get an iterator over the cell volumes of a grid positioned one after the last cell.
CellVolumeIterator endCellVolumes(const Dune::CpGrid& grid);
} // end namespace AutoDiffGrid
} //end namespace OPM

#endif
#endif
