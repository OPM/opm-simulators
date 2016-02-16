/*
  Copyright 2015 SINTEF ICT, Applied Mathematics.

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

#ifndef OPM_FASTSPARSESUM_HEADER_INCLUDED
#define OPM_FASTSPARSESUM_HEADER_INCLUDED

#include <opm/common/utility/platform_dependent/disable_warnings.h>

#include <Eigen/Eigen>
#include <Eigen/Sparse>

#include <opm/common/utility/platform_dependent/reenable_warnings.h>


namespace Opm
{

    /// Add Eigen sparse matrices.
    /// This function is faster than Eigen's own operator+ for these matrices,
    /// but no attempt has been made to make it generic (in terms of matrix
    /// ordering etc). Is is intended for use by the AutoDiffMatrix class.
    inline void fastSparseSum(const Eigen::SparseMatrix<double>& lhs,
                              const Eigen::SparseMatrix<double>& rhs,
                              Eigen::SparseMatrix<double>& sum)
    {
        typedef Eigen::SparseMatrix<double> Mat;
        typedef Mat::Scalar Scalar;
        typedef Mat::Index Index;

        sum = Mat(lhs.rows(), lhs.cols());
        if (lhs.nonZeros() == 0) {
            sum = rhs;
            return;
        }
        if (rhs.nonZeros() == 0) {
            sum = lhs;
            return;
        }
        sum.setZero();
        sum.reserve(lhs.nonZeros() + rhs.nonZeros());

        // we compute each column of the result, one after the other
        const Index cols = lhs.cols();
        using ColEntry = std::pair<int, double>; // Row index and value.
        std::vector<ColEntry> col_entries;
        for (Index col = 0; col < cols; ++col) {
            col_entries.clear();
            Mat::InnerIterator lhs_it(lhs, col);
            Mat::InnerIterator rhs_it(rhs, col);
            while (true) {
                if (!lhs_it) {
                    for (; rhs_it; ++rhs_it) {
                        col_entries.emplace_back(rhs_it.index(), rhs_it.value());
                    }
                    break;
                } else if (!rhs_it) {
                    for (; lhs_it; ++lhs_it) {
                        col_entries.emplace_back(lhs_it.index(), lhs_it.value());
                    }
                    break;
                }
                // Assumes both iterators are valid at this point.
                if (lhs_it.index() == rhs_it.index()) {
                    col_entries.emplace_back(lhs_it.index(), lhs_it.value() + rhs_it.value());
                    ++lhs_it;
                    ++rhs_it;
                } else if (lhs_it.index() < rhs_it.index()) {
                    col_entries.emplace_back(lhs_it.index(), lhs_it.value());
                    ++lhs_it;
                } else {
                    assert(lhs_it.index() > rhs_it.index());
                    col_entries.emplace_back(rhs_it.index(), rhs_it.value());
                    ++rhs_it;
                }
            }
            sum.startVec(col);
            for (const ColEntry& ce : col_entries) {
                sum.insertBackByOuterInnerUnordered(col, ce.first) = ce.second;
            }
        }
        sum.finalize();
        sum.makeCompressed();
        // Mat sumx = lhs + rhs;
        // assert(sum == sumx);
    }

} // namespace Opm



#endif // OPM_FASTSPARSESUM_HEADER_INCLUDED
