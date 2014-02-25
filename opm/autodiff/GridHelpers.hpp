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

#include <functional>

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

namespace UgGridHelpers
{

} //end namespace UgGridHelpers
namespace AutoDiffGrid
{

/// \brief Mapps a grid type to the corresponding face to cell mapping.
///
/// The value of the mapping is provided by the type Type.
template<class T>
struct ADFaceCellTraits
{
};

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

} // end namespace AutoDiffGrid
} // end namespace Opm

#ifdef HAVE_DUNE_CORNERPOINT

#include <dune/common/iteratorfacades.hh>
namespace Opm
{

namespace AutoDiffGrid
{
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
    typedef FaceCellsProxy row_type;
    
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

class Cell2FacesRow
{
public:
    class iterator
        : public Dune::RandomAccessIteratorFacade<iterator,int, int, int>
    {
    public:
        iterator(const Dune::cpgrid::OrientedEntityTable<0,1>::row_type* row,
                 int index)
            : row_(row), index_(index)
        {}

        void increment()
        {
            ++index_;
        }
        void decrement()
        {
            --index_;
        }
        int dereference() const
        {
            return row_->operator[](index_).index();
        }
        int elementAt(int n) const
        {
            return row_->operator[](n).index();
        }
        void advance(int n)
        {
            index_+=n;
        }
        int distanceTo(const iterator& o)const
        {
            return o.index_-index_;
        }
        bool equals(const iterator& o) const
        {
            return index_==o.index_;
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
    typedef  Cell2FacesRow row_type;
    
    Cell2FacesContainer(const Dune::CpGrid* grid)
        : grid_(grid)
    {};
    
    Cell2FacesRow operator[](int cell_index) const
    {
        return Cell2FacesRow(grid_->cellFaceRow(cell_index));
    }
    
        /// \brief Get the number of non-zero entries.
    std::size_t noEntries() const
    {
        return grid_->numCellFaces();
    }
private:
    const Dune::CpGrid* grid_;
};
}

namespace UgGridHelpers
{
template<>
struct Cell2FacesTraits<Dune::CpGrid>
{
    typedef Opm::AutoDiffGrid::Cell2FacesContainer Type;
};
/// \brief An iterator over the cell volumes.
template<const Dune::FieldVector<double, 3>& (Dune::CpGrid::*Method)(int)const>
class CpGridCentroidIterator
    : public Dune::RandomAccessIteratorFacade<CpGridCentroidIterator<Method>, Dune::FieldVector<double, 3>,
                                              const Dune::FieldVector<double, 3>&, int>
{
public:
    /// \brief Creates an iterator.
    /// \param grid The grid the iterator belongs to.
    /// \param cell_index The position of the iterator.
    CpGridCentroidIterator(const  Dune::CpGrid& grid, int cell_index)
        : grid_(&grid), cell_index_(cell_index)
    {}

    const Dune::FieldVector<double, 3>& dereference() const
    {
        return std::mem_fn(Method)(*grid_, cell_index_);
    }
    void increment()
    {
        ++cell_index_;
    }
    const Dune::FieldVector<double, 3>& elementAt(int n) const
    {
        return  std::mem_fn(Method)(*grid_, cell_index_);
    }
    void advance(int n)
    {
        cell_index_+=n;
    }
    void decrement()
    {
        --cell_index_;
    }
    int distanceTo(const CpGridCentroidIterator& o) const
    {
        return o.cell_index_-cell_index_;
    }
    bool equals(const CpGridCentroidIterator& o) const
    {
        return o.grid_==grid_ && o.cell_index_==cell_index_;
    }
    
private:
    const Dune::CpGrid* grid_;
    int cell_index_;
};

template<>
struct CellCentroidTraits<Dune::CpGrid>
{
    typedef CpGridCentroidIterator<&Dune::CpGrid::cellCentroid> IteratorType;
    typedef const double* ValueType;
};

/// \brief Get the number of cells of a grid.
int numCells(const Dune::CpGrid& grid);

/// \brief Get the number of faces of a grid.
int numFaces(const  Dune::CpGrid& grid);

/// \brief Get the dimensions of a grid
int dimensions(const Dune::CpGrid& grid);

/// \brief Get the number of faces, where each face counts as many times as there are adjacent faces
int numCellFaces(const Dune::CpGrid& grid);

/// \brief Get the cartesion dimension of the underlying structured grid.
const int* cartDims(const Dune::CpGrid& grid);

/// \brief Get the local to global index mapping.
///
/// The global index is the index of the active cell
/// in the underlying structured grid.
const int*  globalCell(const Dune::CpGrid&);

CellCentroidTraits<Dune::CpGrid>::IteratorType
beginCellCentroids(const Dune::CpGrid& grid);

/// \brief Get a coordinate of a specific cell centroid.
/// \brief grid The grid.
/// \brief cell_index The index of the specific cell.
/// \breif coordinate The coordinate index.
double cellCentroidCoordinate(const UnstructuredGrid& grid, int cell_index,
                                 int coordinate);

template<>
struct FaceCentroidTraits<Dune::CpGrid>
{
    typedef CpGridCentroidIterator<&Dune::CpGrid::faceCentroid> IteratorType;
    typedef const Dune::CpGrid::Vector ValueType;
};

/// \brief Get an iterator over the face centroids positioned at the first cell.
FaceCentroidTraits<Dune::CpGrid>::IteratorType
beginFaceCentroids(const Dune::CpGrid& grid);

/// \brief Get a coordinate of a specific face centroid.
/// \param grid The grid.
/// \param face_index The index of the specific face.
/// \param coordinate The coordinate index.
FaceCentroidTraits<Dune::CpGrid>::ValueType
faceCentroid(const Dune::CpGrid& grid, int face_index);

template<>
struct FaceCellTraits<Dune::CpGrid>
{
    typedef Opm::AutoDiffGrid::FaceCellsContainerProxy Type;
};
/// \brief Get the cell to faces mapping of a grid.
Opm::AutoDiffGrid::Cell2FacesContainer cell2Faces(const Dune::CpGrid& grid);

/// \brief Get the face to cell mapping of a grid.
FaceCellTraits<Dune::CpGrid>::Type
faceCells(const Dune::CpGrid& grid);

const double* faceNormal(const Dune::CpGrid& grid, int face_index);
} // end namespace UgGridHelperHelpers

namespace AutoDiffGrid
{

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

template<>
struct ADCell2FacesTraits<Dune::CpGrid>
{
    typedef Cell2FacesContainer Type;
};

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

    double dereference() const
    {
        return grid_->cellVolume(cell_index_);
    }
    void increment()
    {
        ++cell_index_;
    }
    double elementAt(int n) const
    {
        return grid_->cellVolume(n);
    }
    void advance(int n)
    {
        cell_index_+=n;
    }
    void decrement()
    {
        --cell_index_;
    }
    int distanceTo(const CellVolumeIterator& o) const
    {
        return o.cell_index_-cell_index_;
    }
    bool equals(const CellVolumeIterator& o) const
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

/// \brief extracts the internal faces of a grid.
/// \param[in] The grid whose internal faces we query.
/// \param[out] internal_faces The internal faces.
/// \param[out] nbi 
void extractInternalFaces(const Dune::CpGrid& grid,
                          Eigen::Array<int, Eigen::Dynamic, 1>& internal_faces,
                          Eigen::Array<int, Eigen::Dynamic, 2, Eigen::RowMajor>& nbi);

template<>
struct ADFaceCellTraits<Dune::CpGrid>
    : public Opm::UgGridHelpers::FaceCellTraits<Dune::CpGrid>
{};
/// \brief Get the face to cell mapping of a grid.
inline ADFaceCellTraits<Dune::CpGrid>::Type
faceCells(const Dune::CpGrid& grid)
{
    return Opm::UgGridHelpers::faceCells(grid);
}
} // end namespace AutoDiffGrid
} //end namespace OPM

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
using Opm::UgGridHelpers::numCellFaces;
using Opm::UgGridHelpers::beginFaceCentroids;
using Opm::UgGridHelpers::beginCellCentroids;

template<>
struct ADFaceCellTraits<UnstructuredGrid>
{
    typedef Eigen::Array<int, Eigen::Dynamic, 2, Eigen::RowMajor> Type;
};

/// \brief Get the face to cell mapping of a grid.
ADFaceCellTraits<UnstructuredGrid>::Type
faceCells(const UnstructuredGrid& grid);

}
}

#endif
