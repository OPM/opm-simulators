/*
  Copyright 2023 Equinor ASA

  This file is part of the Open Porous Media Project (OPM).

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

#include <opm/simulators/utils/compressPartition.hpp>

#include <algorithm>
#include <tuple>
#include <utility>
#include <vector>

namespace {

    template <typename T>
    std::pair<T, T> valueRange(const std::vector<T>& x)
    {
        auto mmPos = std::minmax_element(x.begin(), x.end());

        return { *mmPos.first, *mmPos.second };
    }

    void compressAndCountPartitionIDs(std::vector<int>& partition,
                                      int&              num_domains)
    {
        const auto& [low, high] = valueRange(partition);

        auto seen = std::vector<bool>(high - low + 1, false);
        for (const auto& domain : partition) {
            seen[domain - low] = domain >= 0;
        }

        auto compressed = std::vector<int>(seen.size(), -1);
        for (auto i = 0*compressed.size(); i < compressed.size(); ++i) {
            if (seen[i]) {
                compressed[i] = num_domains++;
            }
        }

        for (auto& domain : partition) {
            if (domain >= 0) {
                domain = compressed[domain - low];
            }
        }
    }
} // Anonymous namespace

std::pair<std::vector<int>, int>
Opm::util::compressAndCountPartitionIDs(std::vector<int>&& parts0)
{
    auto parts = std::pair<std::vector<int>, int> { std::move(parts0), 0 };

    if (! parts.first.empty()) {
        ::compressAndCountPartitionIDs(parts.first, parts.second);
    }

    return parts;
}

std::vector<int> Opm::util::compressPartitionIDs(std::vector<int>&& parts0)
{
    return compressAndCountPartitionIDs(std::move(parts0)).first;
}

void Opm::util::compressPartitionIDs(std::vector<int>& parts0)
{
    [[maybe_unused]] auto num_domains = 0;
    ::compressAndCountPartitionIDs(parts0, num_domains);
}
