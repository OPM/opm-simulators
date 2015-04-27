// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2010-2013 by Andreas Lauser                               *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTBILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the            *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
/*!
 * \file
 * \copydetails Opm::TridiagonalMatrix
 */
#ifndef OPM_TRIDIAGONAL_MATRIX_HH
#define OPM_TRIDIAGONAL_MATRIX_HH

#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>

#include <assert.h>

namespace Opm {

/*!
 * \brief Provides a tridiagonal matrix that also supports non-zero
 *        entries in the upper right and lower left
 *
 * The entries in the lower left and upper right are supported to make
 * implementing periodic systems easy.
 *
 * The API of this class is designed to be close to the one used by
 * the DUNE matrix classes.
 */
template <class Scalar>
class TridiagonalMatrix
{
    struct TridiagRow_ {
        TridiagRow_(TridiagonalMatrix &m, size_t rowIdx)
            : matrix_(m)
            , rowIdx_(rowIdx)
        {};

        Scalar &operator[](size_t colIdx)
        { return matrix_.at(rowIdx_, colIdx); }

        Scalar operator[](size_t colIdx) const
        { return matrix_.at(rowIdx_, colIdx); }

        /*!
         * \brief Prefix increment operator
         */
        TridiagRow_ &operator++()
        { ++ rowIdx_; return *this; }

        /*!
         * \brief Prefix decrement operator
         */
        TridiagRow_ &operator--()
        { -- rowIdx_; return *this; }

        /*!
         * \brief Comparision operator
         */
        bool operator==(const TridiagRow_ &other) const
        { return other.rowIdx_ == rowIdx_ && &other.matrix_ == &matrix_; }

        /*!
         * \brief Comparision operator
         */
        bool operator!=(const TridiagRow_ &other) const
        { return !operator==(other); }

        /*!
         * \brief Dereference operator
         */
        TridiagRow_ &operator*()
        { return *this; }

        /*!
         * \brief Return the row index of the this row.
         *
         * 0 is the first row.
         */
        size_t index() const
        { return rowIdx_; }

    private:
        TridiagonalMatrix &matrix_;
        mutable size_t rowIdx_;
    };

public:
    typedef Scalar FieldType;
    typedef TridiagRow_ RowType;
    typedef size_t SizeType;
    typedef TridiagRow_ iterator;
    typedef TridiagRow_ const_iterator;

    explicit TridiagonalMatrix(int numRows = 0)
    {
        resize(numRows);
    };

    TridiagonalMatrix(int numRows, Scalar value)
    {
        resize(numRows);
        this->operator=(value);
    };

    /*!
     * \brief Copy constructor.
     */
    TridiagonalMatrix(const TridiagonalMatrix &source)
    { *this = source; };

    /*!
     * \brief Return the number of rows/columns of the matrix.
     */
    size_t size() const
    { return diag_[0].size(); }

    /*!
     * \brief Return the number of rows of the matrix.
     */
    size_t rows() const
    { return size(); }

    /*!
     * \brief Return the number of columns of the matrix.
     */
    size_t cols() const
    { return size(); }

    /*!
     * \brief Change the number of rows of the matrix.
     */
    void resize(size_t n)
    {
        if (n == size())
            return;

        for (int diagIdx = 0; diagIdx < 3; ++ diagIdx)
            diag_[diagIdx].resize(n);
    }

    /*!
     * \brief Access an entry.
     */
    Scalar &at(size_t rowIdx, size_t colIdx)
    {
        size_t n = size();

        // special cases
        if (n > 2) {
            if (rowIdx == 0 && colIdx == n - 1)
                return diag_[2][0];
            if (rowIdx == n - 1 && colIdx == 0)
                return diag_[0][n - 1];
        }

        int diagIdx = 1 + colIdx - rowIdx;
        // make sure that the requested column is in range
        assert(0 <= diagIdx && diagIdx < 3);
        return diag_[diagIdx][colIdx];
    }

    /*!
     * \brief Access an entry.
     */
    Scalar at(size_t rowIdx, size_t colIdx) const
    {
        int n = size();

        // special cases
        if (rowIdx == 0 && colIdx == n - 1)
            return diag_[2][0];
        if (rowIdx == n - 1 && colIdx == 0)
            return diag_[0][n - 1];

        int diagIdx = 1 + colIdx - rowIdx;
        // make sure that the requested column is in range
        assert(0 <= diagIdx && diagIdx < 3);
        return diag_[diagIdx][colIdx];
    }

    /*!
     * \brief Assignment operator from another tridiagonal matrix.
     */
    TridiagonalMatrix &operator=(const TridiagonalMatrix &source)
    {
        for (int diagIdx = 0; diagIdx < 3; ++ diagIdx)
            diag_[diagIdx] = source.diag_[diagIdx];

        return *this;
    }

    /*!
     * \brief Assignment operator from a Scalar.
     */
    TridiagonalMatrix &operator=(Scalar value)
    {
        for (int diagIdx = 0; diagIdx < 3; ++ diagIdx)
            diag_[diagIdx].assign(size(), value);

        return *this;
    }

    /*!
     * \begin Iterator for the first row
     */
    iterator begin()
    { return TridiagRow_(*this, 0); }

    /*!
     * \begin Const iterator for the first row
     */
    const_iterator begin() const
    { return TridiagRow_(const_cast<TridiagonalMatrix&>(*this), 0); }

    /*!
     * \begin Const iterator for the next-to-last row
     */
    const_iterator end() const
    { return TridiagRow_(const_cast<TridiagonalMatrix&>(*this), size()); }

    /*!
     * \brief Row access operator.
     */
    TridiagRow_ operator[](size_t rowIdx)
    { return TridiagRow_(*this, rowIdx); }

    /*!
     * \brief Row access operator.
     */
    const TridiagRow_ operator[](size_t rowIdx) const
    { return TridiagRow_(*this, rowIdx); }

    /*!
     * \brief Multiplication with a Scalar
     */
    TridiagonalMatrix &operator*=(Scalar alpha)
    {
        int n = size();
        for (int diagIdx = 0; diagIdx < 3; ++ diagIdx) {
            for (int i = 0; i < n; ++i) {
                diag_[diagIdx][i] *= alpha;
            }
        }

        return *this;
    }

    /*!
     * \brief Division by a Scalar
     */
    TridiagonalMatrix &operator/=(Scalar alpha)
    {
        alpha = 1.0/alpha;
        int n = size();
        for (int diagIdx = 0; diagIdx < 3; ++ diagIdx) {
            for (int i = 0; i < n; ++i) {
                diag_[diagIdx][i] *= alpha;
            }
        }

        return *this;
    }

    /*!
     * \brief Subtraction operator
     */
    TridiagonalMatrix &operator-=(const TridiagonalMatrix &other)
    { return axpy(-1.0, other); }

    /*!
     * \brief Addition operator
     */
    TridiagonalMatrix &operator+=(const TridiagonalMatrix &other)
    { return axpy(1.0, other); }


    /*!
     * \brief Multiply and add the matrix entries of another
     *        tridiagonal matrix.
     *
     * This means that
     * \code
     * A.axpy(alpha, B)
     * \endcode
     * is equivalent to
     * \code
     * A += alpha*C
     * \endcode
     */
    TridiagonalMatrix &axpy(Scalar alpha, const TridiagonalMatrix &other)
    {
        assert(size() == other.size());

        int n = size();
        for (int diagIdx = 0; diagIdx < 3; ++ diagIdx)
            for (int i = 0; i < n; ++ i)
                diag_[diagIdx][i] += alpha * other[diagIdx][i];

        return *this;
    }

    /*!
     * \brief Matrix-vector product
     *
     * This means that
     * \code
     * A.mv(x, y)
     * \endcode
     * is equivalent to
     * \code
     * y = A*x
     * \endcode
     */
    template<class Vector>
    void mv(const Vector &source, Vector &dest) const
    {
        assert(source.size() == size());
        assert(dest.size() == size());
        assert(size() > 1);

        // deal with rows 1 .. n-2
        int n = size();
        for (int i = 1; i < n - 1; ++ i) {
            dest[i] =
                diag_[0][i - 1]*source[i-1] +
                diag_[1][i]*source[i] +
                diag_[2][i + 1]*source[i + 1];
        }

        // rows 0 and n-1
        dest[0] =
            diag_[1][0]*source[0] +
            diag_[2][1]*source[1] +
            diag_[2][0]*source[n - 1];

        dest[n - 1] =
            diag_[0][n-1]*source[0] +
            diag_[0][n-2]*source[n-2] +
            diag_[1][n-1]*source[n-1];
    }

    /*!
     * \brief Additive matrix-vector product
     *
     * This means that
     * \code
     * A.umv(x, y)
     * \endcode
     * is equivalent to
     * \code
     * y += A*x
     * \endcode
     */
    template<class Vector>
    void umv(const Vector &source, Vector &dest) const
    {
        assert(source.size() == size());
        assert(dest.size() == size());
        assert(size() > 1);

        // deal with rows 1 .. n-2
        int n = size();
        for (int i = 1; i < n - 1; ++ i) {
            dest[i] +=
                diag_[0][i - 1]*source[i-1] +
                diag_[1][i]*source[i] +
                diag_[2][i + 1]*source[i + 1];
        }

        // rows 0 and n-1
        dest[0] +=
            diag_[1][0]*source[0] +
            diag_[2][1]*source[1] +
            diag_[2][0]*source[n - 1];

        dest[n - 1] +=
            diag_[0][n-1]*source[0] +
            diag_[0][n-2]*source[n-2] +
            diag_[1][n-1]*source[n-1];
    }

    /*!
     * \brief Subtractive matrix-vector product
     *
     * This means that
     * \code
     * A.mmv(x, y)
     * \endcode
     * is equivalent to
     * \code
     * y -= A*x
     * \endcode
     */
    template<class Vector>
    void mmv(const Vector &source, Vector &dest) const
    {
        assert(source.size() == size());
        assert(dest.size() == size());
        assert(size() > 1);

        // deal with rows 1 .. n-2
        int n = size();
        for (int i = 1; i < n - 1; ++ i) {
            dest[i] -=
                diag_[0][i - 1]*source[i-1] +
                diag_[1][i]*source[i] +
                diag_[2][i + 1]*source[i + 1];
        }

        // rows 0 and n-1
        dest[0] -=
            diag_[1][0]*source[0] +
            diag_[2][1]*source[1] +
            diag_[2][0]*source[n - 1];

        dest[n - 1] -=
            diag_[0][n-1]*source[0] +
            diag_[0][n-2]*source[n-2] +
            diag_[1][n-1]*source[n-1];
    }

    /*!
     * \brief Scaled additive matrix-vector product
     *
     * This means that
     * \code
     * A.usmv(x, y)
     * \endcode
     * is equivalent to
     * \code
     * y += alpha*(A*x)
     * \endcode
     */
    template<class Vector>
    void usmv(Scalar alpha, const Vector &source, Vector &dest) const
    {
        assert(source.size() == size());
        assert(dest.size() == size());
        assert(size() > 1);

        // deal with rows 1 .. n-2
        int n = size();
        for (int i = 1; i < n - 1; ++ i) {
            dest[i] +=
                alpha*(
                    diag_[0][i - 1]*source[i-1] +
                    diag_[1][i]*source[i] +
                    diag_[2][i + 1]*source[i + 1]);
        }

        // rows 0 and n-1
        dest[0] +=
            alpha*(
                diag_[1][0]*source[0] +
                diag_[2][1]*source[1] +
                diag_[2][0]*source[n - 1]);

        dest[n - 1] +=
            alpha*(
                diag_[0][n-1]*source[0] +
                diag_[0][n-2]*source[n-2] +
                diag_[1][n-1]*source[n-1]);
    }

    /*!
     * \brief Transposed matrix-vector product
     *
     * This means that
     * \code
     * A.mtv(x, y)
     * \endcode
     * is equivalent to
     * \code
     * y = A^T*x
     * \endcode
     */
    template<class Vector>
    void mtv(const Vector &source, Vector &dest) const
    {
        assert(source.size() == size());
        assert(dest.size() == size());
        assert(size() > 1);

        // deal with rows 1 .. n-2
        int n = size();
        for (int i = 1; i < n - 1; ++ i) {
            dest[i] =
                diag_[2][i + 1]*source[i-1] +
                diag_[1][i]*source[i] +
                diag_[0][i - 1]*source[i + 1];
        }

        // rows 0 and n-1
        dest[0] =
            diag_[1][0]*source[0] +
            diag_[0][1]*source[1] +
            diag_[0][n-1]*source[n - 1];

        dest[n - 1] =
            diag_[2][0]*source[0] +
            diag_[2][n-1]*source[n-2] +
            diag_[1][n-1]*source[n-1];
    }

    /*!
     * \brief Transposed additive matrix-vector product
     *
     * This means that
     * \code
     * A.umtv(x, y)
     * \endcode
     * is equivalent to
     * \code
     * y += A^T*x
     * \endcode
     */
    template<class Vector>
    void umtv(const Vector &source, Vector &dest) const
    {
        assert(source.size() == size());
        assert(dest.size() == size());
        assert(size() > 1);

        // deal with rows 1 .. n-2
        int n = size();
        for (int i = 1; i < n - 1; ++ i) {
            dest[i] +=
                diag_[2][i + 1]*source[i-1] +
                diag_[1][i]*source[i] +
                diag_[0][i - 1]*source[i + 1];
        }

        // rows 0 and n-1
        dest[0] +=
            diag_[1][0]*source[0] +
            diag_[0][1]*source[1] +
            diag_[0][n-1]*source[n - 1];

        dest[n - 1] +=
            diag_[2][0]*source[0] +
            diag_[2][n-1]*source[n-2] +
            diag_[1][n-1]*source[n-1];
    }

    /*!
     * \brief Transposed subtractive matrix-vector product
     *
     * This means that
     * \code
     * A.mmtv(x, y)
     * \endcode
     * is equivalent to
     * \code
     * y -= A^T*x
     * \endcode
     */
    template<class Vector>
    void mmtv (const Vector &source, Vector &dest) const
    {
        assert(source.size() == size());
        assert(dest.size() == size());
        assert(size() > 1);

        // deal with rows 1 .. n-2
        int n = size();
        for (int i = 1; i < n - 1; ++ i) {
            dest[i] -=
                diag_[2][i + 1]*source[i-1] +
                diag_[1][i]*source[i] +
                diag_[0][i - 1]*source[i + 1];
        }

        // rows 0 and n-1
        dest[0] -=
            diag_[1][0]*source[0] +
            diag_[0][1]*source[1] +
            diag_[0][n-1]*source[n - 1];

        dest[n - 1] -=
            diag_[2][0]*source[0] +
            diag_[2][n-1]*source[n-2] +
            diag_[1][n-1]*source[n-1];
    }

    /*!
     * \brief Transposed scaled additive matrix-vector product
     *
     * This means that
     * \code
     * A.umtv(alpha, x, y)
     * \endcode
     * is equivalent to
     * \code
     * y += alpha*A^T*x
     * \endcode
     */
    template<class Vector>
    void usmtv(Scalar alpha, const Vector &source, Vector &dest) const
    {
        assert(source.size() == size());
        assert(dest.size() == size());
        assert(size() > 1);

        // deal with rows 1 .. n-2
        int n = size();
        for (int i = 1; i < n - 1; ++ i) {
            dest[i] +=
                alpha*(
                    diag_[2][i + 1]*source[i-1] +
                    diag_[1][i]*source[i] +
                    diag_[0][i - 1]*source[i + 1]);
        }

        // rows 0 and n-1
        dest[0] +=
            alpha*(
                diag_[1][0]*source[0] +
                diag_[0][1]*source[1] +
                diag_[0][n-1]*source[n - 1]);

        dest[n - 1] +=
            alpha*(
                diag_[2][0]*source[0] +
                diag_[2][n-1]*source[n-2] +
                diag_[1][n-1]*source[n-1]);
    }

    /*!
     * \brief Calculate the frobenius norm
     *
     * i.e., the square root of the sum of all squared entries. This
     * corresponds to the euclidean norm for vectors.
     */
    Scalar frobeniusNorm() const
    { return std::sqrt(frobeniusNormSquared()); }

    /*!
     * \brief Calculate the squared frobenius norm
     *
     * i.e., the sum of all squared entries.
     */
    Scalar frobeniusNormSquared() const
    {
        Scalar result = 0;

        int n = size();
        for (int i = 0; i < n; ++ i)
            for (int diagIdx = 0; diagIdx < 3; ++ diagIdx)
                result += diag_[diagIdx][i];

        return result;
    }

    /*!
     * \brief Calculate the infinity norm
     *
     * i.e., the maximum of the sum of the absolute values of all rows.
     */
    Scalar infinityNorm() const
    {
        Scalar result=0;

        // deal with rows 1 .. n-2
        int n = size();
        for (int i = 1; i < n - 1; ++ i) {
            result = std::max(result,
                              std::abs(diag_[0][i - 1]) +
                              std::abs(diag_[1][i]) +
                              std::abs(diag_[2][i + 1]));
        }

        // rows 0 and n-1
        result = std::max(result,
                          std::abs(diag_[1][0]) +
                          std::abs(diag_[2][1]) +
                          std::abs(diag_[2][0]));


        result = std::max(result,
                          std::abs(diag_[0][n-1]) +
                          std::abs(diag_[0][n-2]) +
                          std::abs(diag_[1][n-2]));

        return result;
    }

    /*!
     * \brief Calculate the solution for a linear system of equations
     *
     * i.e., calculate x, so that it solves Ax = b, where A is a
     * tridiagonal matrix.
     */
    template <class XVector, class BVector>
    void solve(XVector &x, const BVector &b) const
    {
        if (size() > 2 && diag_[2][0] != 0)
            solveWithUpperRight_(x, b);
        else
            solveWithoutUpperRight_(x, b);
    }

    /*!
     * \brief Print the matrix to a given output stream.
     */
    void print(std::ostream &os = std::cout) const
    {
        int n = size();

        // row 0
        os << at(0, 0) << "\t"
           << at(0, 1) << "\t";

        if (n > 3)
            os << "\t";
        if (n > 2)
            os << at(0, n-1);
        os << "\n";

        // row 1 .. n - 2
        for (int rowIdx = 1; rowIdx < n-1; ++rowIdx) {
            if (rowIdx > 1)
                os << "\t";
            if (rowIdx == n - 2)
                os << "\t";

            os << at(rowIdx, rowIdx - 1) << "\t"
               << at(rowIdx, rowIdx) << "\t"
               << at(rowIdx, rowIdx + 1) << "\n";
        }

        // row n - 1
        if (n > 2)
            os << at(n-1, 0) << "\t";
        if (n > 3)
            os << "\t";
        if (n > 4)
            os << "\t";
        os << at(n-1, n-2) << "\t"
           << at(n-1, n-1) << "\n";
    }

private:
    template <class XVector, class BVector>
    void solveWithUpperRight_(XVector &x, const BVector &b) const
    {
        size_t n = size();

        std::vector<Scalar> lowerDiag(diag_[0]), mainDiag(diag_[1]), upperDiag(diag_[2]), lastColumn(n);
        std::vector<Scalar> bStar(n);
        std::copy(b.begin(), b.end(), bStar.begin());

        lastColumn[0] = upperDiag[0];

        // forward elimination
        for (size_t i = 1; i < n; ++i) {
            Scalar alpha = lowerDiag[i - 1]/mainDiag[i - 1];

            lowerDiag[i - 1] -= alpha * mainDiag[i - 1];
            mainDiag[i] -= alpha * upperDiag[i];

            bStar[i] -= alpha * bStar[i - 1];
        };

        // deal with the last row if the entry on the lower left is not zero
        if (lowerDiag[n - 1] != 0.0 && size() > 2) {
            Scalar lastRow = lowerDiag[n - 1];
            for (size_t i = 0; i < n - 1; ++i) {
                Scalar alpha = lastRow/mainDiag[i];
                lastRow = - alpha*upperDiag[i + 1];
                bStar[n - 1] -= alpha * bStar[i];
            }

            mainDiag[n-1] += lastRow;
        }

        // backward elimination
        x[n - 1] = bStar[n - 1]/mainDiag[n-1];
        for (int i = n - 2; i >= 0; --i) {
            x[i] = (bStar[i] - x[i + 1]*upperDiag[i+1] - x[n-1]*lastColumn[i])/mainDiag[i];
        }
    }

    template <class XVector, class BVector>
    void solveWithoutUpperRight_(XVector &x, const BVector &b) const
    {
        size_t n = size();

        std::vector<Scalar> lowerDiag(diag_[0]), mainDiag(diag_[1]), upperDiag(diag_[2]);
        std::vector<Scalar> bStar(n);
        std::copy(b.begin(), b.end(), bStar.begin());

        // forward elimination
        for (size_t i = 1; i < n; ++i) {
            Scalar alpha = lowerDiag[i - 1]/mainDiag[i - 1];

            lowerDiag[i - 1] -= alpha * mainDiag[i - 1];
            mainDiag[i] -= alpha * upperDiag[i];

            bStar[i] -= alpha * bStar[i - 1];
        };

        // deal with the last row if the entry on the lower left is not zero
        if (lowerDiag[n - 1] != 0.0 && size() > 2) {
            Scalar lastRow = lowerDiag[n - 1];
            for (size_t i = 0; i < n - 1; ++i) {
                Scalar alpha = lastRow/mainDiag[i];
                lastRow = - alpha*upperDiag[i + 1];
                bStar[n - 1] -= alpha * bStar[i];
            }

            mainDiag[n-1] += lastRow;
        }

        // backward elimination
        x[n - 1] = bStar[n - 1]/mainDiag[n-1];
        for (int i = n - 2; i >= 0; --i) {
            x[i] = (bStar[i] - x[i + 1]*upperDiag[i+1])/mainDiag[i];
        }
    }

    mutable std::vector<Scalar> diag_[3];
};

} // namespace Opm

template <class Scalar>
std::ostream &operator<<(std::ostream &os, const Opm::TridiagonalMatrix<Scalar> &mat)
{
    mat.print(os);
    return os;
}

#endif
