/*
  Copyright 2021 Total SE

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

#ifndef OPM_ASPINPARTITION_HEADER_INCLUDED
#define OPM_ASPINPARTITION_HEADER_INCLUDED

#include <functional>
#include <string>
#include <utility>
#include <vector>

namespace Opm {

class Well;

} // namespace Opm

namespace Opm {

/// Control parameters for on-rank subdomain partitioning using Zoltan library.
///
/// \tparam Element Grid view entity type
template <typename Element>
struct ZoltanPartitioningControl
{
    /// Partition imbalance, percentage units.  Assumed to be the run's
    /// value of ZoltanImbalanceTol.
    double domain_imbalance{1.03};

    /// Compute a locally unique, for this MPI process, ID of a local
    /// cell/element/entity.
    std::function<int(const Element&)> index;

    /// Compute a globally unique, across all MPI processes, ID for a local
    /// cell/element/entity.  Might for instance return the cell's
    /// linearised index
    ///
    ///    i + nx*(j + ny*k)
    ///
    /// of the Cartesian cell (i,j,k).
    std::function<int(int)> local_to_global;
};

/// Partition rank's interior cells.
///
/// \param[in] method Partitioning method.  Supported values are \c
///    "zoltan", \c "simple", or a filename with the extension \c
///    ".partition".  The \c "zoltan" method invokes the Zoltan graph
///    partitioning package and requires both MPI and an active Zoltan
///    library.  The \c "simple" method uses a one-dimensional load-balanced
///    approach, and the filename method will read a precomputed partition
///    vector from the named file.
///
/// \param[in] num_local_domains Number of subdomains.  Not used when
///    explicit partitioning is input from a file.
///
/// \param[in] comm MPI Communication object for exchanging globally unique
///    cell IDs and for communication within the Zoltan library.  Not used
///    unless \code method == "zoltan" \endcode.
///
/// \param[in] grid_view View of rank's cells, both interior and overlap
///    cells.  Not used unless \code method == "zoltan" \endcode.
///
/// \param[in] wells Collection of simulation model wells.  Not used unless
///    \code method == "zoltan" \endcode.
///
/// \param[in] zoltan_ctrl Control parameters for local Zoltan-based
///    partitioning.  Not used unless \code method == "zoltan" \endcode.
///
/// \return Partition vector--subdomain ID for each cell in \p grid_view
///    traversal order for its interior cells--and the number of subdomains
///    on current rank.
template <class GridView, class Element>
std::pair<std::vector<int>, int>
partitionCells(const std::string&                        method,
               const int                                 num_local_domains,
               const GridView&                           grid_view,
               const std::vector<Well>&                  wells,
               const ZoltanPartitioningControl<Element>& zoltan_ctrl);

/// Read a partitioning from file, assumed to contain one number per cell, its partition number.
/// \return pair containing a partition vector (partition number for each cell), and the number of partitions.
std::pair<std::vector<int>, int> partitionCellsFromFile(const std::string& filename, const int num_cells);

/// Trivially simple partitioner assigning partitions en bloc, consecutively by cell index.
/// \return pair containing a partition vector (partition number for each cell), and the number of partitions.
std::pair<std::vector<int>, int> partitionCellsSimple(const int num_cells, const int num_domains);

} // namespace Opm

#endif // OPM_ASPINPARTITION_HEADER_INCLUDED
