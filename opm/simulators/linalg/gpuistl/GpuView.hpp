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
#ifndef OPM_GPUVIEW_HEADER_HPP
#define OPM_GPUVIEW_HEADER_HPP

#include <dune/common/fvector.hh>

#include <opm/common/ErrorMacros.hpp>

#include <opm/common/utility/gpuDecorators.hpp>
#include <opm/simulators/linalg/gpuistl/detail/safe_conversion.hpp>

#include <stdexcept>
#include <vector>

#include <fmt/core.h>

namespace Opm::gpuistl
{

/**
 * @brief The GpuView class is provides a view of some data allocated on the GPU
 * Essenstially is only stores a pointer and a size.
 *
 *  This class supports being used from inside a CUDA/HIP Kernel.
 *  Implementations are placed in this headerfile for functions that may be called
 *  inside a kernel to avoid expensive RDC (relocatable device code)
 *
 * The view will typically provide a view into a GpuBuffer and be able to
 * manipulate the data within it
 *
 * @param T Type of the data we store, typically int/float/double w/o const specifier
 *
 **/
template <typename T>
class GpuView
{
public:
    using value_type = T;
    /**
     * @brief Default constructor that will initialize cublas and allocate 0 bytes of memory
     */
    explicit GpuView() = default;

    //TODO: we probably dont need anything like this or is it useful to have views also be able to handle things on CPU?
    /// @brief constructor based on std::vectors, this will make a view on the CPU
    /// @param data std vector to pr
    GpuView(std::vector<T>& data);

    /**
     * @brief operator[] to retrieve a reference to an item in the buffer
     *
     * @param idx The index of the element
     */
    __host__ __device__ T& operator[](size_t idx){
#ifndef NDEBUG
        assertInRange(idx);
#endif
    return m_dataPtr[idx];
    }

    /**
     * @brief operator[] to retrieve a copy of an item in the buffer
     *
     * @param idx The index of the element
     */
    __host__ __device__ T operator[](size_t idx) const {
#ifndef NDEBUG
        assertInRange(idx);
#endif
    return m_dataPtr[idx];
    }


    /**
     * @brief GpuView allocates new GPU memory of size numberOfElements * sizeof(T) and copies numberOfElements from
     * data
     *
     * @note This assumes the data is on the CPU.
     *
     * @param numberOfElements number of T elements to allocate
     * @param dataOnHost data on host/CPU
     */
    GpuView(T* dataOnHost, size_t numberOfElements)
        : m_dataPtr(dataOnHost), m_numberOfElements(numberOfElements)
    {
    }

    /**
     * @brief ~GpuView calls cudaFree
     */
    ~GpuView() = default;

    /**
     * @return the raw pointer to the GPU data
     */
    __host__ __device__ T* data(){
        return m_dataPtr;
    }

    /**
     * @return the raw pointer to the GPU data
     */
    __host__ __device__ const T* data() const{
        return m_dataPtr;
    }

    /**
     * @return true if the view has no elements, false otherwise
     */
    __host__ __device__ bool empty() const {
        return m_numberOfElements == 0;
    }

    /**
     * @return fetch the first element in a GpuView
     */
    __host__ __device__ T& front()
    {
#ifndef NDEBUG
        assertHasElements();
#endif
        return m_dataPtr[0];
    }

    /**
     * @return fetch the last element in a GpuView
     */
    __host__ __device__ T& back()
    {
#ifndef NDEBUG
        assertHasElements();
#endif
        return m_dataPtr[m_numberOfElements-1];
    }

    /**
     * @return fetch the first element in a GpuView
     */
    __host__ __device__ T front() const
    {
#ifndef NDEBUG
        assertHasElements();
#endif
        return m_dataPtr[0];
    }

    /**
     * @return fetch the last element in a GpuView
     */
    __host__ __device__ T back() const
    {
#ifndef NDEBUG
        assertHasElements();
#endif
        return m_dataPtr[m_numberOfElements-1];
    }

    /**
     * @brief copyFromHost copies numberOfElements from the CPU memory dataPointer
     * @param dataPointer raw pointer to CPU memory
     * @param numberOfElements number of elements to copy
     * @note This does synchronous transfer.
     * @note assumes that this view has numberOfElements elements
     */
    void copyFromHost(const T* dataPointer, size_t numberOfElements);

    /**
     * @brief copyFromHost copies numberOfElements to the CPU memory dataPointer
     * @param dataPointer raw pointer to CPU memory
     * @param numberOfElements number of elements to copy
     * @note This does synchronous transfer.
     * @note assumes that this view has numberOfElements elements
     */
    void copyToHost(T* dataPointer, size_t numberOfElements) const;

    /**
     * @brief copyToHost copies data from an std::vector
     * @param data the vector to copy from
     *
     * @note This does synchronous transfer.
     * @note This assumes that the size of this view is equal to the size of the input vector.
     */
    void copyFromHost(const std::vector<T>& data);

    /**
     * @brief copyToHost copies data to an std::vector
     * @param data the vector to copy to
     *
     * @note This does synchronous transfer.
     * @note This assumes that the size of this view is equal to the size of the input vector.
     */
    void copyToHost(std::vector<T>& data) const;

    /**
     * @brief size returns the size (number of T elements) in the vector
     * @return number of elements
     */
    __host__ __device__ size_t size() const{
        return m_numberOfElements;
    }

    /**
     * @brief creates an std::vector of the same size and copies the GPU data to this std::vector
     * @return an std::vector containing the elements copied from the GPU.
     */
    std::vector<T> asStdVector() const;
    /// @brief Iterator class to make GpuViews more similar to std containers
    class iterator {
    public:
        // Iterator typedefs
        using iterator_category = std::forward_iterator_tag;
        using difference_type = std::ptrdiff_t;
        using value_type = T;
        using pointer = T*;
        using reference = T&;

        /// @brief Create iterator from a pointer
        /// @param ptr provided pointer that will become an iterator
        /// @return // the created iterator object
        __host__ __device__ iterator(T* ptr) : m_ptr(ptr) {}

        /// @brief Dereference operator
        /// @return retrieve what the iterator points at
        __host__ __device__ reference operator*() const {
            return *m_ptr;
        }

        /// @brief Pre-increment operator
        /// @return return the pointer after it is incremented
        __host__ __device__ iterator& operator++() {
            ++m_ptr;
            return *this;
        }

        /// @brief Post-increment operator
        /// @param no parameter, int is placeholder for c++ implementation to differentiate from pre-increment
        /// @return Iterator before it is incremented
        __host__ __device__ iterator operator++(int) {
            iterator tmp = *this;
            ++m_ptr;
            return tmp;
        }

        /// @brief Pre-decrement operator
        /// @return return the pointer after it is decremented
        __host__ __device__ iterator& operator--() {
            --m_ptr;
            return *this;
        }

        /// @brief Post-decrement operator
        /// @param no parameter, int is placeholder for c++ implementation to differentiate from pre-decrement
        /// @return Iterator before it is decremented
        __host__ __device__ iterator operator--(int) {
            iterator tmp = *this;
            --m_ptr;
            return tmp;
        }

        /// @brief Inequality comparison operator
        /// @return boolean value that is true if the pointers contains different addresses
        __host__ __device__ bool operator!=(const iterator& other) const {
            return !(m_ptr == other.m_ptr);
        }

        /// @brief Inequality comparison operator
        /// @return boolean value that is true if the pointers contains the same address
        __host__ __device__ bool operator==(const iterator& other) const {
            return m_ptr == other.m_ptr;
        }

        /// @brief subtraction operator
        /// @param other iterator to subtract
        /// @return diffptr that represents difference between the iterators
        __host__ __device__ difference_type operator-(const iterator& other) const {
            return std::distance(other.m_ptr, m_ptr);
        }

        /// @brief Subtraction of given number of elements from iterator
        /// @param n the number of elements to step backwards
        /// @return An iterator pointing to a location n steps behind
        __host__ __device__ iterator operator-(difference_type n) const {
            return iterator(m_ptr-n);
        }

        /// @brief Addition operator with diffptr
        /// @param n diffptr to add
        /// @return new iterator with diffptr added
        __host__ __device__ iterator operator+(difference_type n) const {
            return iterator(m_ptr + n);
        }

        /// @brief Less than comparison
        /// @param other iterator
        /// @return true if this objects iterator is less than the other iterator
        __host__ __device__ bool operator<(const iterator& other) const {
            return m_ptr < other.m_ptr;
        }

        /// @brief Greater than comparison
        /// @param other iterator
        /// @return true if this objects iterator is greater than than the other iterator
        __host__ __device__ bool operator>(const iterator& other) const {
            return m_ptr > other.m_ptr;
        }

    private:
        // Pointer to the current element
        T* m_ptr;
    };

    /**
     * @brief Get an iterator pointing to the first element of the buffer
     * @param iterator to traverse the buffer
     */
    __host__ __device__ iterator begin(){
        return iterator(m_dataPtr);
    }

    /**
     * @brief Get a const iterator pointing to the first element of the buffer
     * @param iterator to traverse the buffer
     */
    __host__ __device__ iterator begin() const {
        return iterator(m_dataPtr);
    }

    /**
     * @brief Get an iterator pointing to the address after the last element of the buffer
     * @param iterator pointing to the first value after the end of the buffer
     */
    __host__ __device__ iterator end(){
        return iterator(m_dataPtr + m_numberOfElements);
    }

    /**
     * @brief Get a const iterator pointing to the address after the last element of the buffer
     * @param iterator pointing to the first value after the end of the buffer
     */
    __host__ __device__ iterator end() const {
        return iterator(m_dataPtr + m_numberOfElements);
    }

private:
    T* m_dataPtr;
    size_t m_numberOfElements;

    /// @brief Helper function to assert if another view has the same size
    /// @param other view
    __host__ __device__ void assertSameSize(const GpuView<T>& other) const
    {
        assertSameSize(other.m_numberOfElements);
    }
    /// @brief Helper function to assert if the size of this view equal to a given value
    /// @param size The value to compare with the size of this view
    __host__ __device__ void assertSameSize(size_t size) const
    {
#if OPM_IS_INSIDE_DEVICE_FUNCTION
        // TODO: find a better way to handle exceptions in kernels, this will possibly be printed many times
        assert(size == m_numberOfElements && "Views did not have the same size");
#else
        if (size != m_numberOfElements) {
            OPM_THROW(std::invalid_argument,
                    fmt::format(fmt::runtime("Given view has {}, while this View has {}."), size, m_numberOfElements));
        }
#endif
    }

    /// @brief Helper function to assert that the view has at least one element
    __host__ __device__ void assertHasElements() const
    {
#if OPM_IS_INSIDE_DEVICE_FUNCTION
        // TODO: find a better way to handle exceptions in kernels, this will possibly be printed many times
        assert(m_numberOfElements > 0 && "View has 0 elements");
#else
        if (m_numberOfElements <= 0) {
            OPM_THROW(std::invalid_argument, "View has 0 elements");
        }
#endif
    }

    /// @brief Helper function to determine if an index is within the range of valid indexes in the view
    __host__ __device__ void assertInRange(size_t idx) const
    {
#if OPM_IS_INSIDE_DEVICE_FUNCTION
        // TODO: find a better way to handle exceptions in kernels, this will possibly be printed many times
        assert(idx < m_numberOfElements && "The index provided was not in the range [0, buffersize-1]");
#else
        if (idx >= m_numberOfElements) {
            OPM_THROW(std::invalid_argument,
                    fmt::format(fmt::runtime("The index provided was not in the range [0, buffersize-1]")));
        }
#endif
    }
};

} // namespace Opm::gpuistl

#endif
