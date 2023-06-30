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

#include <config.h>
#include <opm/simulators/flow/partitionCells.hpp>

#include <opm/grid/CpGrid.hpp>
#include <opm/grid/polyhedralgrid.hh>

#include <opm/input/eclipse/Schedule/Well/Well.hpp>

#include <opm/simulators/flow/countGlobalCells.hpp>

#if HAVE_DUNE_ALUGRID
#include <dune/alugrid/grid.hh>
#endif // HAVE_DUNE_ALUGRID

#include <fmt/format.h>

#include <algorithm>
#include <cassert>
#include <fstream>
#include <iterator>
#include <numeric>
#include <stdexcept>
#include <string>
#include <type_traits>

namespace {

std::pair<std::vector<int>, int>
countDomains(std::vector<int> partition_vector)
{
    auto maxPos = std::max_element(partition_vector.begin(),
                                   partition_vector.end());

    const auto num_domains = (maxPos == partition_vector.end())
        ? 0 : *maxPos + 1;

    return { std::move(partition_vector), num_domains };
}

template <typename, class = void>
struct HasZoltanPartitioning : public std::false_type {};

template <typename GridType>
struct HasZoltanPartitioning<
    GridType,
    std::void_t<decltype(std::declval<const GridType&>().zoltanPartitionWithoutScatter
                         (std::declval<const std::vector<Opm::Well>*>(),
                          std::declval<const double*>(),
                          std::declval<const int>(),
                          std::declval<const double>()))>
    > : public std::true_type {};

} // anonymous namespace

namespace Opm {

template<class Grid>
std::pair<std::vector<int>, int> partitionCells(const Grid& grid,
                                                const std::vector<Well>& wells,
                                                const std::string& method,
                                                const int num_local_domains,
                                                const double partition_imbalance)
{
    if (method == "zoltan") {
        if constexpr (HasZoltanPartitioning<Grid>::value) {
            return partitionCellsZoltan(grid, wells, num_local_domains, partition_imbalance);
        } else {
            OPM_THROW(std::runtime_error, "Zoltan requested for local domain partitioning, "
                      "but is not available for the current grid type.");
        }
    } else if (method == "simple") {
        const int num_cells = detail::countLocalInteriorCells(grid);
        return partitionCellsSimple(num_cells, num_local_domains);
    } else if (method.size() > 10 && method.substr(method.size() - 10, 10) == ".partition") {
        // Method string ends with ".partition", interpret as filename for partitioning.
        const int num_cells = detail::countLocalInteriorCells(grid);
        return partitionCellsFromFile(method, num_cells);
    } else {
        OPM_THROW(std::runtime_error, "Unknown local domain partitioning method requested: " + method);
    }
}


// Read from file, containing one number per cell.
std::pair<std::vector<int>, int> partitionCellsFromFile(const std::string& filename, const int num_cells)
{
    // Read file into single vector.
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
    const int dom_sz = num_cells / num_domains;
    std::vector<int> bounds(num_domains + 1, dom_sz);
    bounds[0] = 0;
    for (int i = 0; i < num_cells % num_domains; ++i) {
        ++bounds[i + 1];
    }
    std::partial_sum(bounds.begin(), bounds.end(), bounds.begin());
    assert(bounds[num_domains] == num_cells);
    std::vector<int> part(num_cells);
    for (int i = 0; i < num_domains; ++i) {
        std::fill(part.begin() + bounds[i], part.begin() + bounds[i + 1], i);
    }
    return { part, num_domains };
}


template<class Grid>
std::pair<std::vector<int>, int> partitionCellsZoltan(const Grid& grid,
                                                      const std::vector<Well>& wells,
                                                      const int num_domains,
                                                      const double domain_imbalance)
{
    auto partition_vector = grid.zoltanPartitionWithoutScatter
        (&wells, nullptr, num_domains, domain_imbalance);

    return countDomains(std::move(partition_vector));
}

template std::pair<std::vector<int>,int>
partitionCells<Dune::CpGrid>(const Dune::CpGrid&,
                             const std::vector<Well>&,
                             const std::string&,
                             const int,
                             const double);

template std::pair<std::vector<int>,int>
partitionCells<Dune::PolyhedralGrid<3,3,double>>(const Dune::PolyhedralGrid<3,3,double>&,
                                                 const std::vector<Well>&,
                                                 const std::string&,
                                                 const int,
                                                 const double);

#if HAVE_DUNE_ALUGRID
#if HAVE_MPI
using ALUGrid3CN = Dune::ALUGrid<3, 3, Dune::cube, Dune::nonconforming, Dune::ALUGridMPIComm>;
#else
using ALUGrid3CN = Dune::ALUGrid<3, 3, Dune::cube, Dune::nonconforming, Dune::ALUGridNoComm>;
#endif //HAVE_MPI

template std::pair<std::vector<int>,int>
partitionCells<ALUGrid3CN>(const ALUGrid3CN&,
                           const std::vector<Well>&,
                           const std::string&,
                           const int,
                           const double);
#endif

} // namespace Opm
