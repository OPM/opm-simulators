/*
  Copyright 2019 Equinor AS.

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
#ifndef MPI_PACKER_HPP
#define MPI_PACKER_HPP

#include <opm/common/utility/TimeService.hpp>
#include <opm/simulators/utils/ParallelCommunication.hpp>

#include <dune/common/parallel/mpitraits.hh>

#include <bitset>
#include <cstddef>
#include <string>

namespace Opm {
namespace Mpi {
namespace detail {

//! \brief Abstract struct for packing which is (partially) specialized for specific types.
template <bool pod, class T>
struct Packing
{
    static std::size_t packSize(const T&, Parallel::MPIComm);
    static void pack(const T&, std::vector<char>&, int&, Parallel::MPIComm);
    static void unpack(T&, std::vector<char>&, int&, Parallel::MPIComm);
};

//! \brief Packaging for pod data.
template<class T>
struct Packing<true,T>
{
    //! \brief Calculates the pack size for a POD.
    //! \param data The data to pack
    //! \param comm The communicator to use
    static std::size_t packSize(const T& data, Parallel::MPIComm comm)
    {
        return packSize(&data, 1, comm);
    }

    //! \brief Calculates the pack size for an array of POD.
    //! \param data The array to pack
    //! \param n Length of array
    //! \param comm The communicator to use
    static std::size_t packSize(const T*, std::size_t n, Parallel::MPIComm comm)
    {
        int size = 0;
        MPI_Pack_size(n, Dune::MPITraits<T>::getType(), comm, &size);
        return size;
    }

    //! \brief Pack a POD.
    //! \param data The variable to pack
    //! \param buffer Buffer to pack into
    //! \param position Position in buffer to use
    //! \param comm The communicator to use
    static void pack(const T& data,
                     std::vector<char>& buffer,
                     int& position,
                     Parallel::MPIComm comm)
    {
        pack(&data, 1, buffer, position, comm);
    }

    //! \brief Pack an array of POD.
    //! \param data The array to pack
    //! \param n Length of array
    //! \param buffer Buffer to pack into
    //! \param position Position in buffer to use
    //! \param comm The communicator to use
    static void pack(const T* data,
                     std::size_t n,
                     std::vector<char>& buffer,
                     int& position,
                     Parallel::MPIComm comm)
    {
        MPI_Pack(data, n, Dune::MPITraits<T>::getType(), buffer.data(),
                 buffer.size(), &position, comm);
    }

    //! \brief Unpack a POD.
    //! \param data The variable to unpack
    //! \param buffer Buffer to unpack from
    //! \param position Position in buffer to use
    //! \param comm The communicator to use
    static void unpack(T& data,
                       std::vector<char>& buffer,
                       int& position,
                       Parallel::MPIComm comm)
    {
        unpack(&data, 1, buffer, position, comm);
    }

    //! \brief Unpack an array of POD.
    //! \param data The array to unpack
    //! \param n Length of array
    //! \param buffer Buffer to unpack from
    //! \param position Position in buffer to use
    //! \param comm The communicator to use
    static void unpack(T* data,
                       std::size_t n,
                       std::vector<char>& buffer,
                       int& position,
                       Parallel::MPIComm comm)
    {
        MPI_Unpack(buffer.data(), buffer.size(), &position, data, n,
                   Dune::MPITraits<T>::getType(), comm);
    }
};

//! \brief Default handling for unsupported types.
template<class T>
struct Packing<false,T>
{
    static std::size_t packSize(const T&, Parallel::MPIComm)
    {
        static_assert(!std::is_same_v<T,T>, "Packing not supported for type");
        return 0;
    }

    static void pack(const T&, std::vector<char>&, int&,
                     Parallel::MPIComm)
    {
      static_assert(!std::is_same_v<T,T>, "Packing not supported for type");
    }

    static void unpack(T&, std::vector<char>&, int&,
                       Parallel::MPIComm)
    {
        static_assert(!std::is_same_v<T,T>, "Packing not supported for type");
    }
};

//! \brief Specialization for std::bitset
template <std::size_t Size>
struct Packing<false,std::bitset<Size>>
{
    static std::size_t packSize(const std::bitset<Size>&, Opm::Parallel::MPIComm);
    static void pack(const std::bitset<Size>&, std::vector<char>&, int&, Opm::Parallel::MPIComm);
    static void unpack(std::bitset<Size>&, std::vector<char>&, int&, Opm::Parallel::MPIComm);
};

#define ADD_PACK_SPECIALIZATION(T) \
    template<> \
    struct Packing<false,T> \
    { \
        static std::size_t packSize(const T&, Parallel::MPIComm); \
        static void pack(const T&, std::vector<char>&, int&, Parallel::MPIComm); \
        static void unpack(T&, std::vector<char>&, int&, Parallel::MPIComm); \
    };

ADD_PACK_SPECIALIZATION(std::string)
ADD_PACK_SPECIALIZATION(time_point)

#undef ADD_PACK_SPECIALIZATION

}

//! \brief Struct handling packing of serialization for MPI communication.
struct Packer {
    //! \brief Constructor.
    //! \param comm The communicator to use
    Packer(Parallel::Communication comm)
        : m_comm(comm)
    {}

    //! \brief Calculates the pack size for a variable.
    //! \tparam T The type of the data to be packed
    //! \param data The data to pack
    template<class T>
    std::size_t packSize(const T& data) const
    {
        return detail::Packing<std::is_pod_v<T>,T>::packSize(data, m_comm);
    }

    //! \brief Calculates the pack size for an array.
    //! \tparam T The type of the data to be packed
    //! \param data The array to pack
    //! \param n Length of array
    template<class T>
    std::size_t packSize(const T* data, std::size_t n) const
    {
        static_assert(std::is_pod_v<T>, "Array packing not supported for non-pod data");
        return detail::Packing<true,T>::packSize(data, n, m_comm);
    }

    //! \brief Pack a variable.
    //! \tparam T The type of the data to be packed
    //! \param data The variable to pack
    //! \param buffer Buffer to pack into
    //! \param position Position in buffer to use
    template<class T>
    void pack(const T& data,
              std::vector<char>& buffer,
              int& position) const
    {
        detail::Packing<std::is_pod_v<T>,T>::pack(data, buffer, position, m_comm);
    }

    //! \brief Pack an array.
    //! \tparam T The type of the data to be packed
    //! \param data The array to pack
    //! \param n Length of array
    //! \param buffer Buffer to pack into
    //! \param position Position in buffer to use
    template<class T>
    void pack(const T* data,
              std::size_t n,
              std::vector<char>& buffer,
              int& position) const
    {
        static_assert(std::is_pod_v<T>, "Array packing not supported for non-pod data");
        detail::Packing<true,T>::pack(data, n, buffer, position, m_comm);
    }

    //! \brief Unpack a variable.
    //! \tparam T The type of the data to be unpacked
    //! \param data The variable to unpack
    //! \param buffer Buffer to unpack from
    //! \param position Position in buffer to use
    template<class T>
    void unpack(T& data,
                std::vector<char>& buffer,
                int& position) const
    {
        detail::Packing<std::is_pod_v<T>,T>::unpack(data, buffer, position, m_comm);
    }

    //! \brief Unpack an array.
    //! \tparam T The type of the data to be unpacked
    //! \param data The array to unpack
    //! \param n Length of array
    //! \param buffer Buffer to unpack from
    //! \param position Position in buffer to use
    template<class T>
    void unpack(T* data,
                std::size_t n,
                std::vector<char>& buffer,
                int& position) const
    {
        static_assert(std::is_pod_v<T>, "Array packing not supported for non-pod data");
        detail::Packing<true,T>::unpack(data, n, buffer, position, m_comm);
    }

private:
    Parallel::Communication m_comm; //!< Communicator to use
};

} // end namespace Mpi
} // end namespace Opm

#endif // MPI_PACKER_HPP
