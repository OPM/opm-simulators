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
#ifndef MPI_SERIALIZER_HPP
#define MPI_SERIALIZER_HPP

#include <opm/common/ErrorMacros.hpp>
#include <opm/common/utility/TimeService.hpp>

#include <opm/simulators/utils/ParallelCommunication.hpp>

#include <bitset>
#include <cstddef>
#include <string>
#include <tuple>
#include <typeinfo>
#include <vector>

namespace Opm
{

namespace Mpi
{
template<class T>
std::size_t packSize(const T*, std::size_t, Opm::Parallel::MPIComm,
                     std::integral_constant<bool, false>);

template<class T>
std::size_t packSize(const T*, std::size_t l, Opm::Parallel::MPIComm comm,
                     std::integral_constant<bool, true>);

template<class T>
std::size_t packSize(const T* data, std::size_t l, Opm::Parallel::MPIComm comm);

template<class T>
std::size_t packSize(const T&, Opm::Parallel::MPIComm,
                     std::integral_constant<bool, false>)
{
    std::string msg = std::string{"Packing not (yet) supported for non-pod type: "} + typeid(T).name();
    OPM_THROW(std::logic_error, msg);
}

template<class T>
std::size_t packSize(const T&, Opm::Parallel::MPIComm comm,
                     std::integral_constant<bool, true>)
{
#if HAVE_MPI
    int size{};
    MPI_Pack_size(1, Dune::MPITraits<T>::getType(), comm, &size);
    return size;
#else
    (void) comm;
    return 0;
#endif
}

template<class T>
std::size_t packSize(const T& data, Opm::Parallel::MPIComm comm)
{
    return packSize(data, comm, typename std::is_pod<T>::type());
}

template<class T1, class T2>
std::size_t packSize(const std::pair<T1,T2>& data, Opm::Parallel::MPIComm comm);

template<class T, class A>
std::size_t packSize(const std::vector<T,A>& data, Opm::Parallel::MPIComm comm);

template<class A>
std::size_t packSize(const std::vector<bool,A>& data, Opm::Parallel::MPIComm comm);

template<class... Ts>
std::size_t packSize(const std::tuple<Ts...>& data, Opm::Parallel::MPIComm comm);

std::size_t packSize(const char* str, Opm::Parallel::MPIComm comm);

template<std::size_t Size>
std::size_t packSize(const std::bitset<Size>& data, Opm::Parallel::MPIComm comm);

////// pack routines

template<class T>
void pack(const T*, std::size_t, std::vector<char>&, int&,
          Opm::Parallel::MPIComm, std::integral_constant<bool, false>);

template<class T>
void pack(const T* data, std::size_t l, std::vector<char>& buffer, int& position,
          Opm::Parallel::MPIComm comm, std::integral_constant<bool, true>);

template<class T>
void pack(const T* data, std::size_t l, std::vector<char>& buffer, int& position,
          Opm::Parallel::MPIComm comm);

template<class T>
void pack(const T&, std::vector<char>&, int&,
          Opm::Parallel::MPIComm, std::integral_constant<bool, false>)
{
    OPM_THROW(std::logic_error, "Packing not (yet) supported for this non-pod type.");
}

template<class T>
void pack(const T& data, std::vector<char>& buffer, int& position,
          Opm::Parallel::MPIComm comm, std::integral_constant<bool, true>)
{
#if HAVE_MPI
    MPI_Pack(&data, 1, Dune::MPITraits<T>::getType(), buffer.data(),
             buffer.size(), &position, comm);
#else
    (void) data;
    (void) comm;
    (void) buffer;
    (void) position;
#endif
}

template<class T>
void pack(const T& data, std::vector<char>& buffer, int& position,
          Opm::Parallel::MPIComm comm)
{
    pack(data, buffer, position, comm, typename std::is_pod<T>::type());
}

template<class T1, class T2>
void pack(const std::pair<T1,T2>& data, std::vector<char>& buffer, int& position,
          Opm::Parallel::MPIComm comm);

template<class T, class A>
void pack(const std::vector<T,A>& data, std::vector<char>& buffer, int& position,
          Opm::Parallel::MPIComm comm);

template<class A>
void pack(const std::vector<bool,A>& data, std::vector<char>& buffer, int& position,
          Opm::Parallel::MPIComm comm);

template<class... Ts>
void pack(const std::tuple<Ts...>& data, std::vector<char>& buffer,
          int& position, Opm::Parallel::MPIComm comm);

void pack(const char* str, std::vector<char>& buffer, int& position,
          Opm::Parallel::MPIComm comm);

template<std::size_t Size>
void pack(const std::bitset<Size>& data, std::vector<char>& buffer, int& position,
          Opm::Parallel::MPIComm comm);

/// unpack routines

template<class T>
void unpack(T*, const std::size_t&, std::vector<char>&, int&,
            Opm::Parallel::MPIComm, std::integral_constant<bool, false>);

template<class T>
void unpack(T* data, const std::size_t& l, std::vector<char>& buffer, int& position,
            Opm::Parallel::MPIComm comm,
            std::integral_constant<bool, true>);

template<class T>
void unpack(T* data, const std::size_t& l, std::vector<char>& buffer, int& position,
            Opm::Parallel::MPIComm comm);

template<class T>
void unpack(T&, std::vector<char>&, int&,
            Opm::Parallel::MPIComm, std::integral_constant<bool, false>)
{
    OPM_THROW(std::logic_error, "Packing not (yet) supported for this non-pod type.");
}

template<class T>
void unpack(T& data, std::vector<char>& buffer, int& position,
            Opm::Parallel::MPIComm comm, std::integral_constant<bool, true>)
{
#if HAVE_MPI
    MPI_Unpack(buffer.data(), buffer.size(), &position, &data, 1,
               Dune::MPITraits<T>::getType(), comm);
#else
    (void) data;
    (void) comm;
    (void) buffer;
    (void) position;
#endif
}

template<class T>
void unpack(T& data, std::vector<char>& buffer, int& position,
            Opm::Parallel::MPIComm comm)
{
    unpack(data, buffer, position, comm, typename std::is_pod<T>::type());
}

template<class T1, class T2>
void unpack(std::pair<T1,T2>& data, std::vector<char>& buffer, int& position,
            Opm::Parallel::MPIComm comm);

template<class T, class A>
void unpack(std::vector<T,A>& data, std::vector<char>& buffer, int& position,
            Opm::Parallel::MPIComm comm);

template<class A>
void unpack(std::vector<bool,A>& data, std::vector<char>& buffer, int& position,
          Opm::Parallel::MPIComm comm);

template<class... Ts>
void unpack(std::tuple<Ts...>& data, std::vector<char>& buffer,
            int& position, Opm::Parallel::MPIComm comm);

void unpack(char* str, std::size_t length, std::vector<char>& buffer, int& position,
            Opm::Parallel::MPIComm comm);

template<std::size_t Size>
void unpack(std::bitset<Size>& data, std::vector<char>& buffer, int& position,
            Opm::Parallel::MPIComm comm);

/// prototypes for complex types

#define ADD_PACK_PROTOTYPES(T) \
  std::size_t packSize(const T& data, Opm::Parallel::MPIComm comm); \
  void pack(const T& data, std::vector<char>& buffer, int& position, \
          Opm::Parallel::MPIComm comm); \
  void unpack(T& data, std::vector<char>& buffer, int& position, \
              Opm::Parallel::MPIComm comm);

ADD_PACK_PROTOTYPES(std::string)
ADD_PACK_PROTOTYPES(time_point)

} // end namespace Mpi

} // end namespace Opm
#endif // MPI_SERIALIZER_HPP
