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

#include <string>
#include <utility>
#include <vector>

namespace Opm
{

    /// Read a partitioning from file, assumed to contain one number per cell, its partition number.
    /// \return pair containing a partition vector (partition number for each cell), and the number of partitions.
    std::pair<std::vector<int>, int> partitionCellsFromFile(const std::string& filename, const int num_cells);

    /// Trivially simple partitioner assigning partitions en bloc, consecutively by cell index.
    /// \return pair containing a partition vector (partition number for each cell), and the number of partitions.
    std::pair<std::vector<int>, int> partitionCellsSimple(const int num_cells, const int num_domains);

} // namespace Opm


#endif // OPM_ASPINPARTITION_HEADER_INCLUDED
