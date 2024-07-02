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
#ifndef OPM_COLRING_AND_REORDERING_UTILS_HPP
#define OPM_COLRING_AND_REORDERING_UTILS_HPP

#include <dune/common/fmatrix.hh>
#include <dune/istl/bcrsmatrix.hh>
#include <fmt/core.h>
#include <opm/common/ErrorMacros.hpp>
#include <opm/simulators/linalg/cuistl/CuSparseMatrix.hpp>
#include <opm/simulators/linalg/cuistl/CuVector.hpp>
#include <opm/simulators/linalg/cuistl/detail/cuda_safe_call.hpp>
#include <opm/simulators/linalg/cuistl/detail/cusparse_matrix_operations.hpp>
#include <opm/simulators/linalg/GraphColoring.hpp>
#include <opm/simulators/linalg/matrixblock.hh>
#include <vector>
#include <config.h>
#include <chrono>
#include <tuple>
/*
This file contains a collection of utility functions used in the GPU implementation of ILU and DILU
The functions deal with creating the mappings between reordered and natural indices, as well as
extracting sparsity structures from dune matrices and creating cusparsematrix indices
*/
namespace Opm::cuistl::detail
{
    inline std::vector<int>
    createReorderedToNatural(Opm::SparseTable<size_t>& levelSets)
    {
        auto res = std::vector<int>(Opm::cuistl::detail::to_size_t(levelSets.dataSize()));
        int globCnt = 0;
        for (auto row : levelSets) {
            for (auto col : row) {
                OPM_ERROR_IF(Opm::cuistl::detail::to_size_t(globCnt) >= res.size(),
                            fmt::format("Internal error. globCnt = {}, res.size() = {}", globCnt, res.size()));
                res[globCnt++] = static_cast<int>(col);
            }
        }
        return res;
    }

    inline std::vector<int>
    createNaturalToReordered(Opm::SparseTable<size_t>& levelSets)
    {
        auto res = std::vector<int>(Opm::cuistl::detail::to_size_t(levelSets.dataSize()));
        int globCnt = 0;
        for (auto row : levelSets) {
            for (auto col : row) {
                OPM_ERROR_IF(Opm::cuistl::detail::to_size_t(globCnt) >= res.size(),
                            fmt::format("Internal error. globCnt = {}, res.size() = {}", globCnt, res.size()));
                res[col] = globCnt++;
            }
        }
        return res;
    }

    template <class M, class field_type, class GPUM>
    inline void
    createReorderedMatrix(const M& naturalMatrix,
                        std::vector<int>& reorderedToNatural,
                        std::unique_ptr<GPUM>& reorderedGpuMat)
    {
        M reorderedMatrix(naturalMatrix.N(), naturalMatrix.N(), naturalMatrix.nonzeroes(), M::row_wise);
        for (auto dstRowIt = reorderedMatrix.createbegin(); dstRowIt != reorderedMatrix.createend(); ++dstRowIt) {
            auto srcRow = naturalMatrix.begin() + reorderedToNatural[dstRowIt.index()];
            // For elements in A
            for (auto elem = srcRow->begin(); elem != srcRow->end(); elem++) {
                dstRowIt.insert(elem.index());
            }
        }

        reorderedGpuMat.reset(new auto (GPUM::fromMatrix(reorderedMatrix, true)));
    }

    template <class M, class field_type, class GPUM>
    inline void
    extractLowerAndUpperMatrices(const M& naturalMatrix,
                                std::vector<int>& reorderedToNatural,
                                std::unique_ptr<GPUM>& lower,
                                std::unique_ptr<GPUM>& upper)
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

        lower.reset(new auto (GPUM::fromMatrix(reorderedLower, true)));
        upper.reset(new auto (GPUM::fromMatrix(reorderedUpper, true)));
        return;
    }
} // END OF NAMESPACE

#endif
