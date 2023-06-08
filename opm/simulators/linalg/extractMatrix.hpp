/*
  Copyright 2021 SINTEF Digital, Mathematics and Cybernetics.

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

#ifndef OPM_EXTRACT_MATRIX_HEADER_INCLUDED
#define OPM_EXTRACT_MATRIX_HEADER_INCLUDED

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <vector>

namespace Opm
{

namespace Details
{

    template <class Matrix>
    void copySubMatrix(const Matrix& A, const std::vector<int>& indices, Matrix& B)
    {
        // Copy elements so that B_{i, j} = A_{indices[i], indices[j]}.
        for (auto row = B.begin(); row != B.end(); ++row) {
            for (auto col = row->begin(); col != row->end(); ++col) {
                *col = A[indices[row.index()]][indices[col.index()]];
            }
        }
    }


    template <class Matrix>
    Matrix extractMatrix(const Matrix& m, const std::vector<int>& indices)
    {
        assert(std::is_sorted(indices.begin(), indices.end()));

        // Set up reverse index map.
        const std::size_t n = indices.size();
        std::size_t newrow = 0;
        enum { NotIncluded = -1 };
        std::vector<int> index_map(m.N(), NotIncluded);
        for (auto row = m.begin(); row != m.end(); ++row) {
            const int row_index = row.index();
            if (row_index == indices[newrow]) {
                index_map[row_index] = newrow;
                ++newrow;
                if (newrow == n) {
                    break;
                }
            }
        }
        assert(newrow == n);

        // Count nonzeroes.
        std::size_t nnz = 0;
        for (auto row = m.begin(); row != m.end(); ++row) {
            if (index_map[row.index()] != NotIncluded) {
                for (auto col = row->begin(); col != row->end(); ++col) {
                    if (index_map[col.index()] != NotIncluded) {
                        ++nnz;
                    }
                }
            }
        }

        // Create the matrix structure.
        Matrix res(n, n, nnz, Matrix::row_wise);
        auto from_row = m.begin();
        for (auto row = res.createbegin(); row != res.createend(); ++row) {
            // Move from_row to point to the row to be extracted.
            while (static_cast<int>(from_row.index()) < indices[row.index()]) {
                ++from_row;
            }
            assert(static_cast<int>(from_row.index()) == indices[row.index()]);
            // Insert nonzeros for row.
            for (auto from_col = from_row->begin(); from_col != from_row->end(); ++from_col) {
                const int new_col = index_map[from_col.index()];
                if (new_col != NotIncluded) {
                    row.insert(new_col);
                }
            }
        }

        copySubMatrix(m, indices, res);
        return res;
    }


    template <class Vector>
    Vector extractVector(const Vector& x, const std::vector<int>& indices)
    {
        Vector res(indices.size());
        for (std::size_t ii = 0; ii < indices.size(); ++ii) {
            res[ii] = x[indices[ii]];
        }
        return res;
    }


    template <class Vector>
    void setGlobal(const Vector& x, const std::vector<int>& indices, Vector& global_x)
    {
        for (std::size_t ii = 0; ii < indices.size(); ++ii) {
            global_x[indices[ii]] = x[ii];
        }
    }



    template <class Matrix>
    bool matrixEqual(const Matrix& m1, const Matrix& m2)
    {
        // Compare size and nonzeroes.
        if (m1.N() != m2.N()) return false;
        if (m1.M() != m2.M()) return false;
        if (m1.nonzeroes() != m2.nonzeroes()) return false;

        auto row1 = m1.begin();
        auto row2 = m2.begin();
        for (; row1 != m1.end(); ++row1, ++row2) {
            if (row2 == m2.end()) return false;
            if (row1.index() != row2.index()) return false;
            auto col1 = row1->begin();
            auto col2 = row2->begin();
            for (; col1 != row1->end(); ++col1, ++col2) {
                if (col2 == row2->end()) return false;
                if (col1.index() != col2.index()) return false;
                if (*col1 != *col2) return false;
            }
        }
        return true;
    }


} // namespace Details


} // namespace Opm

#endif // OPM_EXTRACT_MATRIX_HEADER_INCLUDED
