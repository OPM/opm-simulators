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
#include <dune/grid/polyhedralgrid.hh>
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

/// \brief Mapping of the grid type to the type of the cell to faces mapping.
template<class T>
struct ADCell2FacesTraits
    : public Opm::UgGridHelpers::Cell2FacesTraits<T>
{
};

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
namespace Opm
{

namespace AutoDiffGrid
{

/// \brief Get the z coordinates of the cell centroids of a grid.
/// \return The z coordinates of the cell centroids in an Eigen array
Eigen::Array<double, Eigen::Dynamic, 1>
cellCentroidsZToEigen(const  Dune::CpGrid& grid);

template<>
struct ADCell2FacesTraits<Dune::CpGrid>
{
    typedef Dune::cpgrid::Cell2FacesContainer Type;
};

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
using Opm::UgGridHelpers::cellCentroid;
using Opm::UgGridHelpers::faceCentroid;
using Opm::UgGridHelpers::beginCellVolumes;
using Opm::UgGridHelpers::cellVolume;

template<>
struct ADFaceCellTraits<UnstructuredGrid>
{
    typedef Eigen::Array<int, Eigen::Dynamic, 2, Eigen::RowMajor> Type;
};

#ifdef HAVE_DUNE_CORNERPOINT
// specialization for PolyhedralGrid as a fallback to UnstructuredGrid
template< int dim, int dimworld >
struct ADFaceCellTraits< Dune::PolyhedralGrid< dim, dimworld > >
 : public ADFaceCellTraits<UnstructuredGrid>
{
};
#endif


/// \brief Get the face to cell mapping of a grid.
ADFaceCellTraits<UnstructuredGrid>::Type
faceCellsToEigen(const UnstructuredGrid& grid);

} // end namespace AutoDiffGrid
} //end namespace OPM

#endif
