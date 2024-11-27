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
#ifndef OPM_GPUBUFFER_HEADER_HPP
#define OPM_GPUBUFFER_HEADER_HPP
#include <dune/common/fvector.hh>
#include <dune/istl/bvector.hh>
#include <exception>
#include <fmt/core.h>
#include <opm/common/ErrorMacros.hpp>
#include <opm/simulators/linalg/gpuistl/detail/safe_conversion.hpp>
#include <opm/simulators/linalg/gpuistl/GpuView.hpp>
#include <vector>
#include <string>


namespace Opm::gpuistl
{

/**
 * @brief The GpuBuffer class is a simple container class for the GPU.
 *
 *
 * Example usage:
 *
 * @code{.cpp}
 * #include <opm/simulators/linalg/gpuistl/GpuBuffer.hpp>
 *
 * void someFunction() {
 *     auto someDataOnCPU = std::vector<double>({1.0, 2.0, 42.0, 59.9451743, 10.7132692});
 *
 *     auto dataOnGPU = GpuBuffer<double>(someDataOnCPU);
 *
 *     auto stdVectorOnCPU = dataOnGPU.asStdVector();
 * }
 *
 * @tparam T the type to store. Can be either float, double or int.
 */
template <typename T>
class GpuBuffer
{
public:
    using field_type = T;
    using size_type = size_t;
    using value_type = T;

    /**
     * @brief GpuBuffer allocates new GPU memory of the same size as other and copies the content of the other buffer to
     * this newly allocated memory.
     *
     * @note This does synchronous transfer.
     *
     * @param other the buffer to copy from
     */
    GpuBuffer(const GpuBuffer<T>& other);

    /**
     * @brief GpuBuffer allocates new GPU memory of the same size as data and copies the content of the data vector to
     * this newly allocated memory.
     *
     * @note This does CPU to GPU transfer.
     * @note This does synchronous transfer.
     *
     * @param data the vector to copy from
     */
    explicit GpuBuffer(const std::vector<T>& data);

    /**
     * @brief Default constructor not allocating any memory
     */
    GpuBuffer() = default;

    /**
     * @brief GpuBuffer allocates new GPU memory of size numberOfElements * sizeof(T)
     *
     * @param numberOfElements number of T elements to allocate
     */
    explicit GpuBuffer(const size_t numberOfElements);


    /**
     * @brief GpuBuffer allocates new GPU memory of size numberOfElements * sizeof(T) and copies numberOfElements from
     * data
     *
     * @note This assumes the data is on the CPU.
     *
     * @param numberOfElements number of T elements to allocate
     * @param dataOnHost data on host/CPU
     */
    GpuBuffer(const T* dataOnHost, const size_t numberOfElements);

    /**
     * @brief ~GpuBuffer calls cudaFree
     */
    virtual ~GpuBuffer();

    /**
     * @return the raw pointer to the GPU data
     */
    T* data();

    /**
     * @return the raw pointer to the GPU data
     */
    const T* data() const;

    /**
     * @return fetch the first element in a GpuBuffer
     */
    __host__ __device__ T& front()
    {
#ifndef NDEBUG
        assertHasElements();
#endif
    return m_dataOnDevice[0];
    }

    /**
     * @return fetch the last element in a GpuBuffer
     */
    __host__ __device__ T& back()
    {
#ifndef NDEBUG
        assertHasElements();
#endif
    return m_dataOnDevice[m_numberOfElements-1];
    }

    /**
     * @return fetch the first element in a GpuBuffer
     */
    __host__ __device__ T front() const
    {
#ifndef NDEBUG
        assertHasElements();
#endif
    return m_dataOnDevice[0];
    }

    /**
     * @return fetch the last element in a GpuBuffer
     */
    __host__ __device__ T back() const
    {
#ifndef NDEBUG
        assertHasElements();
#endif
    return m_dataOnDevice[m_numberOfElements-1];
    }

    /**
     * @brief copyFromHost copies data from a Dune::BlockVector
     * @param bvector the vector to copy from
     *
     * @note This does synchronous transfer.
     * @note This assumes that the size of this vector is equal to the size of the input vector.
     */
    template <int BlockDimension>
    void copyFromHost(const Dune::BlockVector<Dune::FieldVector<T, BlockDimension>>& bvector)
    {
        // TODO: [perf] vector.size() can be replaced by bvector.N() * BlockDimension
        if (m_numberOfElements != bvector.size()) {
            OPM_THROW(std::runtime_error,
                      fmt::format("Given incompatible vector size. GpuBuffer has size {}, \n"
                                  "however, BlockVector has N() = {}, and size = {}.",
                                  m_numberOfElements,
                                  bvector.N(),
                                  bvector.size()));
        }
        const auto dataPointer = static_cast<const T*>(&(bvector[0][0]));
        copyFromHost(dataPointer, m_numberOfElements);
    }

    /**
     * @brief copyToHost copies data to a Dune::BlockVector
     * @param bvector the vector to copy to
     *
     * @note This does synchronous transfer.
     * @note This assumes that the size of this vector is equal to the size of the input vector.
     */
    template <int BlockDimension>
    void copyToHost(Dune::BlockVector<Dune::FieldVector<T, BlockDimension>>& bvector) const
    {
        // TODO: [perf] vector.size() can be replaced by bvector.N() * BlockDimension
        if (m_numberOfElements != bvector.size()) {
            OPM_THROW(std::runtime_error,
                      fmt::format("Given incompatible vector size. GpuBuffer has size {},\n however, the BlockVector "
                                  "has has N() = {}, and size() = {}.",
                                  m_numberOfElements,
                                  bvector.N(),
                                  bvector.size()));
        }
        const auto dataPointer = static_cast<T*>(&(bvector[0][0]));
        copyToHost(dataPointer, m_numberOfElements);
    }

    /**
     * @brief copyFromHost copies numberOfElements from the CPU memory dataPointer
     * @param dataPointer raw pointer to CPU memory
     * @param numberOfElements number of elements to copy
     * @note This does synchronous transfer.
     * @note assumes that this buffer has numberOfElements elements
     */
    void copyFromHost(const T* dataPointer, size_t numberOfElements);

    /**
     * @brief copyFromHost copies numberOfElements to the CPU memory dataPointer
     * @param dataPointer raw pointer to CPU memory
     * @param numberOfElements number of elements to copy
     * @note This does synchronous transfer.
     * @note assumes that this buffer has numberOfElements elements
     */
    void copyToHost(T* dataPointer, size_t numberOfElements) const;

    /**
     * @brief copyToHost copies data from an std::vector
     * @param data the vector to copy from
     *
     * @note This does synchronous transfer.
     * @note This assumes that the size of this buffer is equal to the size of the input vector.
     */
    void copyFromHost(const std::vector<T>& data);

    /**
     * @brief copyToHost copies data to an std::vector
     * @param data the vector to copy to
     *
     * @note This does synchronous transfer.
     * @note This assumes that the size of this buffer is equal to the size of the input vector.
     */
    void copyToHost(std::vector<T>& data) const;

    /**
     * @brief size returns the size (number of T elements) in the buffer
     * @return number of elements
     */
    size_type size() const;

    /**
     * @brief resize the number of elements that fit in the buffer, shrinking it causes truncation
     * @param number of elements in the new buffer
     */
    void resize(size_t);

    /**
     * @brief creates an std::vector of the same size and copies the GPU data to this std::vector
     * @return an std::vector containing the elements copied from the GPU.
     */
    std::vector<T> asStdVector() const;

private:
    T* m_dataOnDevice = nullptr;
    size_t m_numberOfElements = 0;

    void assertSameSize(const GpuBuffer<T>& other) const;
    void assertSameSize(size_t size) const;

    void assertHasElements() const;
};

template <class T>
GpuView<const T> make_view(const GpuBuffer<T>&);

} // namespace Opm::gpuistl
#endif
