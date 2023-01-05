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

#include <opm/simulators/flow/aspinPartition.hpp>

#include <fmt/format.h>

#include <algorithm>
#include <cassert>
#include <fstream>
#include <iterator>
#include <numeric>
#include <string>

namespace Opm
{


namespace
{

    // Read from file, containing one number per cell, from [0, ... , num_domains - 1].
    std::pair<std::vector<int>, int> partitionCellsFromFile([[maybe_unused]] const int num_cells)
    {
        // TODO: refactor to make more flexible.
        // Read file into single vector.
        const std::string filename = "partition.txt";
        std::ifstream is(filename);
        const std::vector<int> cellpart{std::istream_iterator<int>(is), std::istream_iterator<int>()};
        if (cellpart.size() != size_t(num_cells)) {
            auto msg = fmt::format("Partition file contains {} entries, but there are {} cells.",
                                   cellpart.size(), num_cells);
            throw std::runtime_error(msg);
        }

        // Create and return the output domain vector.
        const int num_domains = (*std::max_element(cellpart.begin(), cellpart.end())) + 1;
        return { cellpart, num_domains };
    }


    // Trivially simple partitioner
    std::pair<std::vector<int>, int> partitionCellsSimple(const int num_cells, const int num_domains)
    {
        // Build the partitions.
        std::vector<int> bounds(num_domains + 1);
        bounds[0] = 0;
        const int dom_sz = num_cells / num_domains;
        for (int i = 1; i < num_domains; ++i) {
            bounds[i] = dom_sz * i;
        }
        bounds[num_domains] = num_cells;
        std::vector<int> part(num_cells);
        for (int i = 0; i < num_domains; ++i) {
            std::fill(part.begin() + bounds[i], part.begin() + bounds[i + 1], i);
        }
        return { part, num_domains };
    }

} // anonymous namespace


std::pair<std::vector<int>, int>
partitionCells(const int num_cells)
{
    // const std::string method = "simple";
    const std::string method = "file";
    if (method == "simple") {
        return partitionCellsSimple(num_cells, 1);
    } else if (method == "file") {
        return partitionCellsFromFile(num_cells);
    } else {
        return {};
    }
}

} // namespace Opm
