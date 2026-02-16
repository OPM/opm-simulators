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
#include <opm/simulators/linalg/gpuistl/detail/gpu_safe_call.hpp>
#include <opm/simulators/linalg/gpuistl/detail/safe_conversion.hpp>
#include <opm/simulators/linalg/gpuistl/GpuView.hpp>
#include <vector>
#include <string>
#include <cuda_runtime.h>


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
     * @brief Default constructor not allocating any memory
     */
    GpuBuffer() = default;

    /**
     * @brief GpuBuffer allocates new GPU memory of the same size as other and copies the content of the other buffer to
     * this newly allocated memory.
     *
     * @note This does synchronous transfer.
     *
     * @param other the buffer to copy from
     */
    GpuBuffer(const GpuBuffer<T>& other)
        : GpuBuffer(other.m_numberOfElements)
    {
        assertSameSize(other);
        if (m_numberOfElements == 0) {
            return;
        }
        OPM_GPU_SAFE_CALL(cudaMemcpy(m_dataOnDevice,
                                    other.m_dataOnDevice,
                                    m_numberOfElements * sizeof(T),
                                    cudaMemcpyDeviceToDevice));
    }

    /**
     * @brief GpuBuffer allocates new GPU memory of the same size as data and copies the content of the data vector to
     * this newly allocated memory.
     *
     * @note This does CPU to GPU transfer.
     * @note This does synchronous transfer.
     *
     * @param data the vector to copy from
     */
    explicit GpuBuffer(const std::vector<T>& data)
        : GpuBuffer(data.size())
    {
        copyFromHost(data);
    }

    /**
     * @brief GpuBuffer allocates new GPU memory of size numberOfElements * sizeof(T)
     *
     * @param numberOfElements number of T elements to allocate
     */
    explicit GpuBuffer(const size_t numberOfElements)
        : m_numberOfElements(numberOfElements)
    {
        OPM_GPU_SAFE_CALL(cudaMalloc(&m_dataOnDevice, sizeof(T) * m_numberOfElements));
    }


    /**
     * @brief GpuBuffer allocates new GPU memory of size numberOfElements * sizeof(T) and copies numberOfElements from
     * data
     *
     * @note This assumes the data is on the CPU.
     *
     * @param numberOfElements number of T elements to allocate
     * @param dataOnHost data on host/CPU
     */
    GpuBuffer(const T* dataOnHost, const size_t numberOfElements)
        : GpuBuffer(numberOfElements)
    {
        OPM_GPU_SAFE_CALL(cudaMemcpy(
            m_dataOnDevice, dataOnHost, m_numberOfElements * sizeof(T), cudaMemcpyHostToDevice));
    }


    /**
     * @brief ~GpuBuffer calls cudaFree
     */
    virtual ~GpuBuffer()
    {
        OPM_GPU_WARN_IF_ERROR(cudaFree(m_dataOnDevice));
    }

    /**
     * @return the raw pointer to the GPU data
     */
    T* data()
    {
        return m_dataOnDevice;
    }

    /**
     * @return the raw pointer to the GPU data
     */
    const T* data() const
    {
        return m_dataOnDevice;
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
    void copyFromHost(const T* dataPointer, size_t numberOfElements)
    {
        if (numberOfElements > size()) {
            OPM_THROW(std::runtime_error,
                    fmt::format(fmt::runtime("Requesting to copy too many elements. "
                                             "buffer has {} elements, while {} was requested."),
                                size(),
                                numberOfElements));
        }
        OPM_GPU_SAFE_CALL(cudaMemcpy(data(), dataPointer, numberOfElements * sizeof(T), cudaMemcpyHostToDevice));
    }

    /**
     * @brief copyFromHost copies numberOfElements to the CPU memory dataPointer
     * @param dataPointer raw pointer to CPU memory
     * @param numberOfElements number of elements to copy
     * @note This does synchronous transfer.
     * @note assumes that this buffer has numberOfElements elements
     */
    void copyToHost(T* dataPointer, size_t numberOfElements) const
    {
        assertSameSize(numberOfElements);
        OPM_GPU_SAFE_CALL(cudaMemcpy(dataPointer, data(), numberOfElements * sizeof(T), cudaMemcpyDeviceToHost));
    }

    /**
     * @brief copyToHost copies data from an std::vector
     * @param data the vector to copy from
     *
     * @note This does synchronous transfer.
     * @note This assumes that the size of this buffer is equal to the size of the input vector.
     */
    void copyFromHost(const std::vector<T>& data)
    {
        assertSameSize(data.size());

        if (data.empty()) {
            return;
        }

        if constexpr (std::is_same_v<T, bool>)
        {
            auto tmp = std::make_unique<bool[]>(data.size());
            for (size_t i = 0; i < data.size(); ++i) {
                tmp[i] = static_cast<bool>(data[i]);
            }
            copyFromHost(tmp.get(), data.size());
        }
        else {
            copyFromHost(data.data(), data.size());
        }
    }

    /**
     * @brief copyToHost copies data to an std::vector
     * @param data the vector to copy to
     *
     * @note This does synchronous transfer.
     * @note This assumes that the size of this buffer is equal to the size of the input vector.
     */
    void copyToHost(std::vector<T>& data) const
    {
        assertSameSize(data.size());

        if (data.empty()) {
            return;
        }

        if constexpr (std::is_same_v<T, bool>)
        {
            auto tmp = std::make_unique<bool[]>(data.size());
            copyToHost(tmp.get(), data.size());
            for (size_t i = 0; i < data.size(); ++i) {
                data[i] = static_cast<bool>(tmp[i]);
            }
            return;
        }
        else {
            copyToHost(data.data(), data.size());
        }
    }

    /**
     * @brief size returns the size (number of T elements) in the buffer
     * @return number of elements
     */
    size_type size() const
    {
        return m_numberOfElements;
    }

    /**
     * @brief resize the number of elements that fit in the buffer, shrinking it causes truncation
     * @param number of elements in the new buffer
     */
    void resize(size_t newSize)
    {
        if (newSize < 1) {
            OPM_THROW(std::invalid_argument, "Setting a GpuBuffer size to a non-positive number is not allowed");
        }

        if (m_numberOfElements == 0) {
            // We have no data, so we can just allocate new memory
            OPM_GPU_SAFE_CALL(cudaMalloc(&m_dataOnDevice, sizeof(T) * newSize));
        }
        else {
            // Allocate memory for temporary buffer
            T* tmpBuffer = nullptr;
            OPM_GPU_SAFE_CALL(cudaMalloc(&tmpBuffer, sizeof(T) * m_numberOfElements));

            // Move the data from the old to the new buffer with truncation
            size_t sizeOfMove = std::min({m_numberOfElements, newSize});
            OPM_GPU_SAFE_CALL(cudaMemcpy(tmpBuffer,
                                        m_dataOnDevice,
                                        sizeOfMove * sizeof(T),
                                        cudaMemcpyDeviceToDevice));

            // free the old buffer
            OPM_GPU_SAFE_CALL(cudaFree(m_dataOnDevice));

            // swap the buffers
            m_dataOnDevice = tmpBuffer;
        }

        // update size
        m_numberOfElements = newSize;
    }

    /**
     * @brief creates an std::vector of the same size and copies the GPU data to this std::vector
     * @return an std::vector containing the elements copied from the GPU.
     */
    std::vector<T> asStdVector() const
    {
        std::vector<T> temporary(m_numberOfElements);
        copyToHost(temporary);
        return temporary;
    }

private:
    T* m_dataOnDevice = nullptr;
    size_t m_numberOfElements = 0;

    void assertSameSize(const GpuBuffer<T>& other) const
    {
        assertSameSize(other.m_numberOfElements);
    }

    void assertSameSize(size_t size) const
    {
        if (size != m_numberOfElements) {
            OPM_THROW(std::invalid_argument,
                    fmt::format(fmt::runtime("Given buffer has {}, while we have {}."),
                                size, m_numberOfElements));
        }
    }

    void assertHasElements() const
    {
        if (m_numberOfElements <= 0) {
            OPM_THROW(std::invalid_argument, "We have 0 elements");
        }
    }
};

template <class T>
GpuView<T> make_view(GpuBuffer<T>& buf) {
    return GpuView<T>(buf.data(), buf.size());
}

template <class T>
GpuView<const T> make_view(const GpuBuffer<T>& buf) {
    return GpuView<const T>(buf.data(), buf.size());
}

} // namespace Opm::gpuistl
#endif
