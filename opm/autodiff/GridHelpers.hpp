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
#ifndef OPM_AUTODIFF_GRIDHELPERS_HEADER_INCLUDED
#define OPM_AUTODIFF_GRIDHELPERS_HEADER_INCLUDED

#include <functional>

#include <boost/range/iterator_range.hpp>
#include <opm/core/grid.h>
#include <opm/core/grid/GridHelpers.hpp>

#include <opm/core/utility/platform_dependent/disable_warnings.h>

#include <Eigen/Eigen>
#include <Eigen/Sparse>

#ifdef HAVE_DUNE_CORNERPOINT
#include <dune/grid/CpGrid.hpp>
#include <dune/grid/cpgrid/GridHelpers.hpp>
#endif

#include <opm/core/utility/platform_dependent/reenable_warnings.h>


namespace Opm
{

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
cellCentroidsZToEigen(const UnstructuredGrid& grid);

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
}

namespace UgGridHelpers
{
} // end namespace UgGridHelperHelpers

namespace AutoDiffGrid
{

/// \brief Get the z coordinates of the cell centroids of a grid.
/// \return The z coordinates of the cell centroids in an Eigen array
Eigen::Array<double, Eigen::Dynamic, 1>
cellCentroidsZToEigen(const  Dune::CpGrid& grid);

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
    typedef Dune::cpgrid::Cell2FacesContainer Type;
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
ADFaceCellTraits<Dune::CpGrid>::Type
faceCellsToEigen(const Dune::CpGrid& grid);
} // end namespace AutoDiffGrid
} //end namespace OPM

#endif
namespace Opm
{
namespace AutoDiffGrid
{

using Opm::UgGridHelpers::SparseTableView;
using Opm::UgGridHelpers::numCells;
using Opm::UgGridHelpers::faceCells;
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
faceCellsToEigen(const UnstructuredGrid& grid);

} // end namespace AutoDiffGrid
} //end namespace OPM

#endif
