/*
  Copyright 2025 Equinor ASA

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

// Minimal fixed-size matrix class for CUDA kernels
#ifndef OPM_MINIMATRIX_HPP
#define OPM_MINIMATRIX_HPP
#include <array>
#include <cstddef>
#include <initializer_list>

#include <opm/common/utility/gpuDecorators.hpp>
#include <opm/simulators/linalg/gpuistl/MiniVector.hpp>
#include <opm/simulators/linalg/matrixblock.hh>

namespace Opm::gpuistl
{

/**
 * @brief A small fixed-size square matrix class for use in CUDA kernels.
 * @tparam T Element type.
 * @tparam dimension Number of rows and columns (matrix is dimension x dimension).
 */
template<class T, int dimension>
class MiniMatrix {
public:
	using value_type = T;
	using size_type = std::size_t;
	using array_type = std::array<T, dimension * dimension>;

	/**
	 * @brief Default constructor. Elements are default-initialized.
	 */
	OPM_HOST_DEVICE MiniMatrix() = default;

	/**
	 * @brief Constructor from initializer list.
	 * @param init Initializer list of elements (row-major order).
	 */
	OPM_HOST_DEVICE MiniMatrix(std::initializer_list<T> init) {
		size_type i = 0;
		for (auto it = init.begin(); it != init.end() && i < data_.size(); ++it, ++i)
        {
			data_[i] = *it;
        }
	}

    // We need a conversion from CPU based matrix block to these GPU-friendly MiniMatrices
    // Some copy-to-gpu routines will convert types using this ctor
    template <class Type, int n, int m>
    MiniMatrix(MatrixBlock<Type, n, m> mb) {
        static_assert(n == m, "MiniMatrix can only be constructed from square MatrixBlock");
        static_assert(n == dimension, "MiniMatrix dimension must match MatrixBlock dimension");
        for (size_type i = 0; i < dimension; ++i)
        {
            for (size_type j = 0; j < dimension; ++j)
            {
                (*this)(i, j) = mb(i, j);
            }
        }
    }

    OPM_HOST_DEVICE MiniMatrix(const T& value) {
        fill(value);
    }

	/**
	 * @brief Access element at (row, col).
	 * @param row Row index.
	 * @param col Column index.
	 * @return Reference to element.
	 */
	OPM_HOST_DEVICE T& operator()(size_type row, size_type col) {
		return data_[row * dimension + col];
	}

	/**
	 * @brief Access element at (row, col) (const version).
	 * @param row Row index.
	 * @param col Column index.
	 * @return Const reference to element.
	 */
	OPM_HOST_DEVICE const T& operator()(size_type row, size_type col) const {
		return data_[row * dimension + col];
	}

	/**
	 * @brief Get iterator to beginning of data.
	 */
	OPM_HOST_DEVICE auto begin() { return data_.begin(); }

	/**
	 * @brief Get iterator to end of data.
	 */
	OPM_HOST_DEVICE auto end() { return data_.end(); }

	/**
	 * @brief Get const iterator to beginning of data.
	 */
	OPM_HOST_DEVICE auto begin() const { return data_.begin(); }

	/**
	 * @brief Get const iterator to end of data.
	 */
	OPM_HOST_DEVICE auto end() const { return data_.end(); }

	/**
	 * @brief Get pointer to raw data.
	 */
	OPM_HOST_DEVICE T* data() { return data_.data(); }

	/**
	 * @brief Get const pointer to raw data.
	 */
	OPM_HOST_DEVICE const T* data() const { return data_.data(); }

	/**
	 * @brief Fill all elements with a value.
	 * @param value Value to fill.
	 */
	OPM_HOST_DEVICE void fill(const T& value) {
		for (auto& x : data_) x = value;
	}

	/**
	 * @brief Get matrix dimension.
	 */
	OPM_HOST_DEVICE static constexpr size_type size() { return dimension; }

    /**
     * @brief Add another matrix to this one (element-wise).
     * @param other Matrix to add.
     * @return Reference to this matrix.
     */
    OPM_HOST_DEVICE MiniMatrix& operator+=(const MiniMatrix& other) {
        for (size_type i = 0; i < data_.size(); ++i)
            data_[i] += other.data_[i];
        return *this;
    }

	/**
	 * @brief Add two matrices (element-wise).
	 * @param other Matrix to add.
	 * @return Resulting matrix.
	 */
    OPM_HOST_DEVICE MiniMatrix operator+(const MiniMatrix& other) const {
        MiniMatrix result = *this;
        result += other;
        return result;
    }

    /**
     * @brief Subtract another matrix from this one (element-wise).
     * @param other Matrix to subtract.
     * @return Reference to this matrix.
     */
    OPM_HOST_DEVICE MiniMatrix& operator-=(const MiniMatrix& other) {
        for (size_type i = 0; i < data_.size(); ++i)
            data_[i] -= other.data_[i];
        return *this;
    }

	/**
	 * @brief Subtract two matrices (element-wise).
	 * @param other Matrix to subtract.
	 * @return Resulting matrix.
	 */
	OPM_HOST_DEVICE MiniMatrix operator-(const MiniMatrix& other) const {
		MiniMatrix result = *this;
		result -= other;
		return result;
	}

    /**
     * @brief Matrix-vector multiplication: y = A * x
     * @param x Input vector (length = dimension)
     * @return Resulting vector (length = dimension)
     */
    OPM_HOST_DEVICE std::array<T, dimension> operator*(const std::array<T, dimension>& x) const {
        std::array<T, dimension> result{};
        for (size_type row = 0; row < dimension; ++row) {
            T sum = T{};
            for (size_type col = 0; col < dimension; ++col)
                sum += (*this)(row, col) * x[col];
            result[row] = sum;
        }
        return result;
    }

    /**
     * @brief Matrix-matrix multiplication: C = A * B
     * @param B Right-hand side matrix
     * @return Resulting matrix
     */
    OPM_HOST_DEVICE MiniMatrix operator*(const MiniMatrix& B) const {
        MiniMatrix C;
        for (size_type row = 0; row < dimension; ++row) {
            for (size_type col = 0; col < dimension; ++col) {
                T sum = T{};
                for (size_type k = 0; k < dimension; ++k)
                    sum += (*this)(row, k) * B(k, col);
                C(row, col) = sum;
            }
        }
        return C;
    }

	/**
	 * @brief Matrix-vector multiplication: y = A * x, with MiniVector.
	 * @param x Input MiniVector (length = dimension)
	 * @return Resulting MiniVector (length = dimension)
	 */
	OPM_HOST_DEVICE Opm::gpuistl::MiniVector<T, dimension> operator*(const Opm::gpuistl::MiniVector<T, dimension>& x) const {
		Opm::gpuistl::MiniVector<T, dimension> result{};
		for (size_type row = 0; row < dimension; ++row) {
			T sum = T{};
			for (size_type col = 0; col < dimension; ++col)
				sum += (*this)(row, col) * x[col];
			result[row] = sum;
		}
		return result;
	}

private:
	array_type data_;
};

} // namespace Opm::gpuistl

#endif // OPM_MINIMATRIX_HPP
