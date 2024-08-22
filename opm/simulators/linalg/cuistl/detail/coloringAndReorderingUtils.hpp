/*
  Copyright 2024 SINTEF AS

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
#ifndef OPM_COLORING_AND_REORDERING_UTILS_HPP
#define OPM_COLORING_AND_REORDERING_UTILS_HPP

#include <fmt/core.h>
#include <memory>
#include <opm/common/ErrorMacros.hpp>
#include <opm/grid/utility/SparseTable.hpp>
#include <opm/simulators/linalg/cuistl/detail/safe_conversion.hpp>
#include <tuple>
#include <vector>
/*
This file contains a collection of utility functions used in the GPU implementation of ILU and DILU
The functions deal with creating the mappings between reordered and natural indices, as well as
extracting sparsity structures from dune matrices and creating cusparsematrix indices
*/
namespace Opm::gpuistl::detail
{
inline std::vector<int>
createReorderedToNatural(const Opm::SparseTable<size_t>& levelSets)
{
    auto res = std::vector<int>(Opm::gpuistl::detail::to_size_t(levelSets.dataSize()));
    int globCnt = 0;
    for (auto row : levelSets) {
        for (auto col : row) {
            OPM_ERROR_IF(Opm::gpuistl::detail::to_size_t(globCnt) >= res.size(),
                         fmt::format("Internal error. globCnt = {}, res.size() = {}", globCnt, res.size()));
            res[globCnt++] = static_cast<int>(col);
        }
    }
    return res;
}

inline std::vector<int>
createNaturalToReordered(const Opm::SparseTable<size_t>& levelSets)
{
    auto res = std::vector<int>(Opm::gpuistl::detail::to_size_t(levelSets.dataSize()));
    int globCnt = 0;
    for (auto row : levelSets) {
        for (auto col : row) {
            OPM_ERROR_IF(Opm::gpuistl::detail::to_size_t(globCnt) >= res.size(),
                         fmt::format("Internal error. globCnt = {}, res.size() = {}", globCnt, res.size()));
            res[col] = globCnt++;
        }
    }
    return res;
}

template <class M, class field_type, class GPUM>
inline std::unique_ptr<GPUM>
createReorderedMatrix(const M& naturalMatrix, const std::vector<int>& reorderedToNatural)
{
    M reorderedMatrix(naturalMatrix.N(), naturalMatrix.N(), naturalMatrix.nonzeroes(), M::row_wise);
    for (auto dstRowIt = reorderedMatrix.createbegin(); dstRowIt != reorderedMatrix.createend(); ++dstRowIt) {
        auto srcRow = naturalMatrix.begin() + reorderedToNatural[dstRowIt.index()];
        // For elements in A
        for (auto elem = srcRow->begin(); elem != srcRow->end(); elem++) {
            dstRowIt.insert(elem.index());
        }
    }
    return std::unique_ptr<GPUM>(new auto(GPUM::fromMatrix(reorderedMatrix, true)));
}

template <class M, class field_type, class GPUM>
inline std::tuple<std::unique_ptr<GPUM>, std::unique_ptr<GPUM>>
extractLowerAndUpperMatrices(const M& naturalMatrix, const std::vector<int>& reorderedToNatural)
{
    const size_t new_nnz = (naturalMatrix.nonzeroes() - naturalMatrix.N()) / 2;

    M reorderedLower(naturalMatrix.N(), naturalMatrix.N(), new_nnz, M::row_wise);
    M reorderedUpper(naturalMatrix.N(), naturalMatrix.N(), new_nnz, M::row_wise);

    for (auto lowerIt = reorderedLower.createbegin(), upperIt = reorderedUpper.createbegin();
         lowerIt != reorderedLower.createend();
         ++lowerIt, ++upperIt) {

        auto srcRow = naturalMatrix.begin() + reorderedToNatural[lowerIt.index()];

        for (auto elem = srcRow->begin(); elem != srcRow->end(); ++elem) {
            if (elem.index() < srcRow.index()) { // add index to lower matrix if under the diagonal
                lowerIt.insert(elem.index());
            } else if (elem.index() > srcRow.index()) { // add element to upper matrix if above the diagonal
                upperIt.insert(elem.index());
            }
        }
    }

    return {std::unique_ptr<GPUM>(new auto(GPUM::fromMatrix(reorderedLower, true))),
            std::unique_ptr<GPUM>(new auto(GPUM::fromMatrix(reorderedUpper, true)))};
}
} // namespace Opm::gpuistl::detail

#endif
