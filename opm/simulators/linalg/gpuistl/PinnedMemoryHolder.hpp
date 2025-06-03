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

#ifndef OPM_SIMULATORS_LINALG_GPUISTL_PINNEDMEMORYHOLDER_HPP
#define OPM_SIMULATORS_LINALG_GPUISTL_PINNEDMEMORYHOLDER_HPP


#include "opm/simulators/linalg/gpuistl/detail/gpu_safe_call.hpp"

#include <cuda_runtime.h>

#include <fmt/format.h>

namespace Opm
{
namespace gpuistl
{

    /**
     * @brief RAII class for pinning host memory using cudaHostRegister.
     *
     * This class registers a given host memory region in its constructor
     * and unregisters it in its destructor. This is useful for speeding up
     * asynchronous memory transfers between host and GPU.
     *
     * @tparam T The type of data pointed to by the host pointer.
     */
    template <class T>
    class PinnedMemoryHolder
    {
    public:
        /**
         * @brief Constructs a PinnedMemoryHolder and registers the host memory.
         * @param ptr Pointer to the host memory to be pinned.
         * @param numberOfElements The number of elements of type T in the memory region.
         * @param flags Flags for cudaHostRegister. Defaults to cudaHostRegisterDefault.
         *
         * @throws std::runtime_error if cudaHostRegister fails.
         */
        PinnedMemoryHolder(T* ptr, std::size_t numberOfElements, unsigned int flags = cudaHostRegisterDefault)
            : m_ptr(ptr)
            , m_numberOfElements(numberOfElements)
        {
            if (m_ptr && m_numberOfElements > 0) {
                std::size_t sizeInBytes = numberOfElements * sizeof(T);
                OPM_GPU_SAFE_CALL(cudaHostRegister(static_cast<void*>(m_ptr), sizeInBytes, flags));
            }
        }

        /**
         * @brief Destructor. Unregisters the host memory.
         */
        ~PinnedMemoryHolder()
        {
            if (m_ptr && m_numberOfElements > 0)
                OPM_GPU_SAFE_CALL(cudaHostUnregister(static_cast<void*>(m_ptr)));
        }

        // Disable copy constructor and copy assignment operator
        PinnedMemoryHolder(const PinnedMemoryHolder&) = delete;
        PinnedMemoryHolder& operator=(const PinnedMemoryHolder&) = delete;

        // Enable move constructor and move assignment operator
        PinnedMemoryHolder(PinnedMemoryHolder&& other) noexcept
            : m_ptr(other.m_ptr)
            , m_numberOfElements(other.m_numberOfElements)
        {
            other.m_ptr = nullptr;
            other.m_numberOfElements = 0;
        }

        PinnedMemoryHolder& operator=(PinnedMemoryHolder&& other) noexcept
        {
            if (this != &other) {
                // Unregister existing memory, if any
                if (m_ptr && m_numberOfElements > 0) {
                    OPM_GPU_SAFE_CALL(cudaHostUnregister(static_cast<void*>(m_ptr)));
                }

                // Transfer ownership
                m_ptr = other.m_ptr;
                m_numberOfElements = other.m_numberOfElements;

                // Leave other in a valid but empty state
                other.m_ptr = nullptr;
                other.m_numberOfElements = 0;
            }
            return *this;
        }

        /**
         * @brief Gets the pointer to the pinned memory.
         * @return Pointer to the pinned memory, or nullptr if not valid.
         */
        T* get() const
        {
            return m_ptr;
        }

        /**
         * @brief Gets the number of elements in the pinned memory region.
         * @return Number of elements.
         */
        std::size_t numberOfElements() const
        {
            return m_numberOfElements;
        }

    private:
        T* m_ptr;
        std::size_t m_numberOfElements;
    };

} // namespace gpuistl
} // namespace Opm

#endif // OPM_SIMULATORS_LINALG_GPUISTL_PINNEDMEMORYHOLDER_HPP
