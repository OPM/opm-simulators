/*
  Copyright 2022-2023 SINTEF AS

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
#ifndef OPM_CUVECTOR_HEADER_HPP
#define OPM_CUVECTOR_HEADER_HPP
#include <dune/common/fvector.hh>
#include <dune/istl/bvector.hh>
#include <exception>
#include <fmt/core.h>
#include <opm/common/ErrorMacros.hpp>
#include <opm/simulators/linalg/cuistl/detail/CuBlasHandle.hpp>
#include <opm/simulators/linalg/cuistl/detail/safe_conversion.hpp>
#include <vector>
#include <string>


namespace Opm::cuistl
{

/**
 * @brief The CuVector class is a simple (arithmetic) vector class for the GPU.
 *
 * @note we currently only support simple raw primitives for T (double, float and int)
 *
 * @note We currently only support arithmetic operations on double and float.
 *
 * @note this vector has no notion of block size. The user is responsible for allocating
 *       the correct number of primitives (double or floats)
 *
 * Example usage:
 *
 * @code{.cpp}
 * #include <opm/simulators/linalg/cuistl/CuVector.hpp>
 *
 * void someFunction() {
 *     auto someDataOnCPU = std::vector<double>({1.0, 2.0, 42.0, 59.9451743, 10.7132692});
 *
 *     auto dataOnGPU = CuVector<double>(someDataOnCPU);
 *
 *     // Multiply by 4.0:
 *     dataOnGPU *= 4.0;
 *
 *     // Get data back on CPU in another vector:
 *     auto stdVectorOnCPU = dataOnGPU.asStdVector();
 * }
 *
 * @tparam T the type to store. Can be either float, double or int.
 */
template <typename T>
class CuVector
{
public:
    using field_type = T;
    using size_type = size_t;


    /**
     * @brief CuVector allocates new GPU memory of the same size as other and copies the content of the other vector to
     * this newly allocated memory.
     *
     * @note This does synchronous transfer.
     *
     * @param other the vector to copy from
     */
    CuVector(const CuVector<T>& other);

    /**
     * @brief CuVector allocates new GPU memory of the same size as data and copies the content of the data vector to
     * this newly allocated memory.
     *
     * @note This does CPU to GPU transfer.
     * @note This does synchronous transfer.
     *
     * @note For now data.size() needs to be within the limits of int due to restrctions of CuBlas.
     *
     * @param data the vector to copy from
     */
    explicit CuVector(const std::vector<T>& data);

    /**
     * @brief operator= copies the content of the data vector to the memory of this vector.
     *
     * @note This requires the two vectors to be of the same size.
     * @note This does synchronous transfer.
     *
     * @param other the vector to copy from
     */
    CuVector& operator=(const CuVector<T>& other);

    /**
     * @brief operator= sets the whole vector equal to the scalar value.
     *
     * @note This does asynchronous operations
     *
     * @param scalar the value all elements will be set to.
     */
    CuVector& operator=(T scalar);

    /**
     * @brief CuVector allocates new GPU memory of size numberOfElements * sizeof(T)
     *
     * @note For now numberOfElements needs to be within the limits of int due to restrictions in cublas
     *
     * @param numberOfElements number of T elements to allocate
     */
    explicit CuVector(const size_t numberOfElements);


    /**
     * @brief CuVector allocates new GPU memory of size numberOfElements * sizeof(T) and copies numberOfElements from
     * data
     *
     * @note This assumes the data is on the CPU.
     *
     * @param numberOfElements number of T elements to allocate
     * @param dataOnHost data on host/CPU
     *
     * @note For now numberOfElements needs to be within the limits of int due to restrictions in cublas
     */
    CuVector(const T* dataOnHost, const size_t numberOfElements);

    /**
     * @brief ~CuVector calls cudaFree
     */
    virtual ~CuVector();

    /**
     * @return the raw pointer to the GPU data
     */
    T* data();

    /**
     * @return the raw pointer to the GPU data
     */
    const T* data() const;

    /**
     * @brief copyFromHost copies data from a Dune::BlockVector
     * @param bvector the vector to copy from
     *
     * @note This does synchronous transfer.
     * @note This assumes that the size of this vector is equal to the dim of the input vector.
     */
    template <int BlockDimension>
    void copyFromHost(const Dune::BlockVector<Dune::FieldVector<T, BlockDimension>>& bvector)
    {
        // TODO: [perf] vector.dim() can be replaced by bvector.N() * BlockDimension
        if (detail::to_size_t(m_numberOfElements) != bvector.dim()) {
            OPM_THROW(std::runtime_error,
                      fmt::format("Given incompatible vector size. CuVector has size {}, \n"
                                  "however, BlockVector has N() = {}, and dim = {}.",
                                  m_numberOfElements,
                                  bvector.N(),
                                  bvector.dim()));
        }
        const auto dataPointer = static_cast<const T*>(&(bvector[0][0]));
        copyFromHost(dataPointer, m_numberOfElements);
    }

    /**
     * @brief copyToHost copies data to a Dune::BlockVector
     * @param bvector the vector to copy to
     *
     * @note This does synchronous transfer.
     * @note This assumes that the size of this vector is equal to the dim of the input vector.
     */
    template <int BlockDimension>
    void copyToHost(Dune::BlockVector<Dune::FieldVector<T, BlockDimension>>& bvector) const
    {
        // TODO: [perf] vector.dim() can be replaced by bvector.N() * BlockDimension
        if (detail::to_size_t(m_numberOfElements) != bvector.dim()) {
            OPM_THROW(std::runtime_error,
                      fmt::format("Given incompatible vector size. CuVector has size {},\n however, the BlockVector "
                                  "has has N() = {}, and dim() = {}.",
                                  m_numberOfElements,
                                  bvector.N(),
                                  bvector.dim()));
        }
        const auto dataPointer = static_cast<T*>(&(bvector[0][0]));
        copyToHost(dataPointer, m_numberOfElements);
    }

    /**
     * @brief copyFromHost copies numberOfElements from the CPU memory dataPointer
     * @param dataPointer raw pointer to CPU memory
     * @param numberOfElements number of elements to copy
     * @note This does synchronous transfer.
     * @note assumes that this vector has numberOfElements elements
     */
    void copyFromHost(const T* dataPointer, size_t numberOfElements);

    /**
     * @brief copyFromHost copies numberOfElements to the CPU memory dataPointer
     * @param dataPointer raw pointer to CPU memory
     * @param numberOfElements number of elements to copy
     * @note This does synchronous transfer.
     * @note assumes that this vector has numberOfElements elements
     */
    void copyToHost(T* dataPointer, size_t numberOfElements) const;

    /**
     * @brief copyToHost copies data from an std::vector
     * @param data the vector to copy from
     *
     * @note This does synchronous transfer.
     * @note This assumes that the size of this vector is equal to the size of the input vector.
     */
    void copyFromHost(const std::vector<T>& data);

    /**
     * @brief copyToHost copies data to an std::vector
     * @param data the vector to copy to
     *
     * @note This does synchronous transfer.
     * @note This assumes that the size of this vector is equal to the size of the input vector.
     */
    void copyToHost(std::vector<T>& data) const;

    /**
     * @brief operator *= multiplies every element by scalar
     * @param scalar the scalar to with which to multiply every element
     *
     * @note This operation is asynchronous.
     *
     * @note int is not supported
     */
    CuVector<T>& operator*=(const T& scalar);

    /**
     * @brief axpy sets this vector equal to this + alha * y
     * @param alpha the scalar with which to multiply y
     * @param y input vector of same size as this
     *
     * @note this will call CuBlas in the background
     * @note int is not supported
     */
    CuVector<T>& axpy(T alpha, const CuVector<T>& y);

    /**
     * @brief operator+= adds the other vector to this vector
     *
     * @note this will call CuBlas in the background
     * @note int is not supported
     */
    CuVector<T>& operator+=(const CuVector<T>& other);

    /**
     * @brief operator-= subtracts the other vector from this vector
     *
     * @note this will call CuBlas in the background
     * @note int is not supported
     */
    CuVector<T>& operator-=(const CuVector<T>& other);

    /**
     * @brief dot computes the dot product (standard inner product) against the other vector
     * @param other vector of same size as this
     * @note this will call CuBlas in the background
     * @note int is not supported
     *
     * @return the result on the inner product
     */
    T dot(const CuVector<T>& other) const;

    /**
     * @brief returns the l2 norm of the vector
     * @note this will call CuBlas in the background
     * @note int is not supported
     *
     * @return the l2 norm
     */
    T two_norm() const;

    /**
     * Computes the dot product sum_i this[indexSet[i]] * other[indexSet[i]]
     *
     * @note int is not supported
     */
    T dot(const CuVector<T>& other, const CuVector<int>& indexSet, CuVector<T>& buffer) const;

    /**
     * Computes the norm sqrt(sum_i this[indexSet[i]] * this[indexSet[i]])
     *
     * @note int is not supported
     */
    T two_norm(const CuVector<int>& indexSet, CuVector<T>& buffer) const;


    /**
     * Computes the dot product sum_i this[indexSet[i]] * other[indexSet[i]]
     *
     * @note int is not supported
     */
    T dot(const CuVector<T>& other, const CuVector<int>& indexSet) const;

    /**
     * Computes the norm sqrt(sum_i this[indexSet[i]] * this[indexSet[i]])
     *
     * @note int is not supported
     */
    T two_norm(const CuVector<int>& indexSet) const;


    /**
     * @brief dim returns the dimension (number of T elements) in the vector
     * @return number of elements
     */
    size_type dim() const;


    /**
     * @brief creates an std::vector of the same size and copies the GPU data to this std::vector
     * @return an std::vector containing the elements copied from the GPU.
     */
    std::vector<T> asStdVector() const;

    /**
     * @brief creates an std::vector of the same size and copies the GPU data to this std::vector
     * @return an std::vector containing the elements copied from the GPU.
     */
    template <int blockSize>
    Dune::BlockVector<Dune::FieldVector<T, blockSize>> asDuneBlockVector() const
    {
        OPM_ERROR_IF(dim() % blockSize != 0,
                     fmt::format("blockSize is not a multiple of dim(). Given blockSize = {}, and dim() = {}",
                                 blockSize,
                                 dim()));

        Dune::BlockVector<Dune::FieldVector<T, blockSize>> returnValue(dim() / blockSize);
        copyToHost(returnValue);
        return returnValue;
    }


    /**
     * @brief setZeroAtIndexSet for each element in indexSet, sets the index of this vector to be zero
     * @param indexSet the set of indices to set to zero
     *
     * @note Assumes all indices are within range
     *
     * This is supposed to do the same as the following code on the CPU:
     * @code{.cpp}
     * for (int index : indexSet) {
     *     this->data[index] = T(0.0);
     * }
     * @endcode
     */
    void setZeroAtIndexSet(const CuVector<int>& indexSet);

    // Slow method that creates a string representation of a CuVector for debug purposes
    std::string toDebugString()
    {
        std::vector<T> v = asStdVector();
        std::string res = "";
        for (T element : v){
            res += std::to_string(element) + " ";
        }
        res += std::to_string(v[v.size()-1]);
        return res;
    }

private:
    T* m_dataOnDevice = nullptr;

    // Note that we store this as int to make sure we are always cublas compatible.
    // This gives the added benefit that a size_t to int conversion error occurs during construction.
    const int m_numberOfElements;
    detail::CuBlasHandle& m_cuBlasHandle;

    void assertSameSize(const CuVector<T>& other) const;
    void assertSameSize(int size) const;

    void assertHasElements() const;
};

} // namespace Opm::cuistl
#endif
