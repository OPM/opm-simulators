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
#include <config.h>
#include "MPIPacker.hpp"

#include <cstdint>
#include <cstring>
#include <ctime>

#include <dune/common/parallel/mpitraits.hh>
#if HAVE_MPI
#include <ebos/eclmpiserializer.hh>
#endif

namespace Opm
{
namespace Mpi
{
template<class T>
std::size_t packSize(const T*, std::size_t, Opm::Parallel::MPIComm,
                     std::integral_constant<bool, false>)
{
    OPM_THROW(std::logic_error, "Packing not (yet) supported for this non-pod type.");
}

template<class T>
std::size_t packSize(const T*, std::size_t l, Opm::Parallel::MPIComm comm,
                     std::integral_constant<bool, true>)
{
#if HAVE_MPI
    int size;
    MPI_Pack_size(1, Dune::MPITraits<std::size_t>::getType(), comm, &size);
    std::size_t totalSize = size;
    MPI_Pack_size(l, Dune::MPITraits<T>::getType(), comm, &size);
    return totalSize + size;
#else
    (void) comm;
    return l-l;
#endif
}

template<class T>
std::size_t packSize(const T* data, std::size_t l, Opm::Parallel::MPIComm comm)
{
    return packSize(data, l, comm, typename std::is_pod<T>::type());
}

std::size_t packSize(const char* str, Opm::Parallel::MPIComm comm)
{
#if HAVE_MPI
    int size;
    MPI_Pack_size(1, Dune::MPITraits<std::size_t>::getType(), comm, &size);
    int totalSize = size;
    MPI_Pack_size(strlen(str)+1, MPI_CHAR, comm, &size);
    return totalSize + size;
#else
    (void) str;
    (void) comm;
    return 0;
#endif
}

std::size_t packSize(const std::string& str, Opm::Parallel::MPIComm comm)
{
    return packSize(str.c_str(), comm);
}

template <class T>
struct Packing
{
};

template <std::size_t Size>
struct Packing<std::bitset<Size>>
{
    static std::size_t packSize(const std::bitset<Size>& data, Opm::Parallel::MPIComm comm)
    {
        return Mpi::packSize(data.to_ullong(), comm);
    }

    static void pack(const std::bitset<Size>& data, std::vector<char>& buffer, int& position, Opm::Parallel::MPIComm comm)
    {
        Mpi::pack(data.to_ullong(), buffer, position, comm);
    }

    static void unpack(std::bitset<Size>& data, std::vector<char>& buffer, int& position, Opm::Parallel::MPIComm comm)
    {
        unsigned long long d;
        Mpi::unpack(d, buffer, position, comm);
        data = std::bitset<Size>(d);
    }
};

template<std::size_t Size>
std::size_t packSize(const std::bitset<Size>& data, Opm::Parallel::MPIComm comm)
{
    return Packing<std::bitset<Size>>::packSize(data, comm);
}

std::size_t packSize(const Opm::time_point&, Opm::Parallel::MPIComm comm)
{
    std::time_t tp = 0;
    return packSize(tp, comm);
}


////// pack routines

template<class T>
void pack(const T*, std::size_t, std::vector<char>&, int&,
          Opm::Parallel::MPIComm, std::integral_constant<bool, false>)
{
    OPM_THROW(std::logic_error, "Packing not (yet) supported for this non-pod type.");
}

template<class T>
void pack(const T* data, std::size_t l, std::vector<char>& buffer, int& position,
          Opm::Parallel::MPIComm comm,
          std::integral_constant<bool, true>)
{
#if HAVE_MPI
    MPI_Pack(&l, 1, Dune::MPITraits<std::size_t>::getType(), buffer.data(),
             buffer.size(), &position, comm);
    MPI_Pack(data, l, Dune::MPITraits<T>::getType(), buffer.data(),
             buffer.size(), &position, comm);
#else
    (void) data;
    (void) comm;
    (void) l;
    (void) buffer;
    (void) position;
#endif
}

template<class T>
void pack(const T* data, std::size_t l, std::vector<char>& buffer, int& position,
          Opm::Parallel::MPIComm comm)
{
    pack(data, l, buffer, position, comm, typename std::is_pod<T>::type());
}

void pack(const char* str, std::vector<char>& buffer, int& position,
          Opm::Parallel::MPIComm comm)
{
#if HAVE_MPI
    std::size_t length = strlen(str)+1;
    MPI_Pack(&length, 1, Dune::MPITraits<std::size_t>::getType(), buffer.data(),
        buffer.size(), &position, comm);
    MPI_Pack(str, strlen(str)+1, MPI_CHAR, buffer.data(), buffer.size(),
         &position, comm);
#else
    (void) str;
    (void) comm;
    (void) buffer;
    (void) position;
#endif
}

void pack(const std::string& str, std::vector<char>& buffer, int& position,
          Opm::Parallel::MPIComm comm)
{
    pack(str.c_str(), buffer, position, comm);
}

template<std::size_t Size>
void pack(const std::bitset<Size>& data, std::vector<char>& buffer,
          int& position, Opm::Parallel::MPIComm comm)
{
    Packing<std::bitset<Size>>::pack(data, buffer, position, comm);
}

void pack(const Opm::time_point& data, std::vector<char>& buffer, int& position,
          Opm::Parallel::MPIComm comm)
{
    pack(Opm::TimeService::to_time_t(data), buffer, position, comm);
}


/// Mpi::unpack routines

template<class T>
void unpack(T*, const std::size_t&, std::vector<char>&, int&,
            Opm::Parallel::MPIComm, std::integral_constant<bool, false>)
{
    OPM_THROW(std::logic_error, "Packing not (yet) supported for this non-pod type.");
}

template<class T>
void unpack(T* data, const std::size_t& l, std::vector<char>& buffer, int& position,
            Opm::Parallel::MPIComm comm,
            std::integral_constant<bool, true>)
{
#if HAVE_MPI
    MPI_Unpack(buffer.data(), buffer.size(), &position, data, l,
               Dune::MPITraits<T>::getType(), comm);
#else
    (void) data;
    (void) comm;
    (void) l;
    (void) buffer;
    (void) position;
#endif
}

template<class T>
void unpack(T* data, const std::size_t& l, std::vector<char>& buffer, int& position,
            Opm::Parallel::MPIComm comm)
{
    unpack(data, l, buffer, position, comm, typename std::is_pod<T>::type());
}

void unpack(char* str, std::size_t length, std::vector<char>& buffer, int& position,
            Opm::Parallel::MPIComm comm)
{
#if HAVE_MPI
    MPI_Unpack(buffer.data(), buffer.size(), &position, const_cast<char*>(str), length, MPI_CHAR, comm);
#else
    (void) str;
    (void) comm;
    (void) length;
    (void) buffer;
    (void) position;
#endif
}

void unpack(std::string& str, std::vector<char>& buffer, int& position,
            Opm::Parallel::MPIComm comm)
{
    std::size_t length=0;
    unpack(length, buffer, position, comm);
    std::vector<char> cStr(length, '\0');
    unpack(cStr.data(), length, buffer, position, comm);
    str.clear();
    str.append(cStr.data());
}

template<std::size_t Size>
void unpack(std::bitset<Size>& data, std::vector<char>& buffer, int& position,
            Opm::Parallel::MPIComm comm)
{
    Packing<std::bitset<Size>>::unpack(data, buffer, position, comm);
}

void unpack([[maybe_unused]] Opm::time_point& data, std::vector<char>& buffer, int& position,
            Opm::Parallel::MPIComm comm)
{
    std::time_t tp;
    unpack(tp, buffer, position, comm);
#if HAVE_MPI
    data = Opm::TimeService::from_time_t(tp);
#endif
}

#define INSTANTIATE_PACK(...) \
template std::size_t packSize(const __VA_ARGS__& data, \
                              Opm::Parallel::MPIComm comm); \
template void pack(const __VA_ARGS__& data, \
                   std::vector<char>& buffer, int& position, \
                   Opm::Parallel::MPIComm comm); \
template void unpack(__VA_ARGS__& data, \
                     std::vector<char>& buffer, int& position, \
                     Opm::Parallel::MPIComm comm);

INSTANTIATE_PACK(float)
INSTANTIATE_PACK(double)
INSTANTIATE_PACK(bool)
INSTANTIATE_PACK(int)
INSTANTIATE_PACK(unsigned char)
INSTANTIATE_PACK(unsigned int)
INSTANTIATE_PACK(unsigned long int)
INSTANTIATE_PACK(unsigned long long int)
INSTANTIATE_PACK(std::bitset<4>)


#undef INSTANTIATE_PACK

} // end namespace Mpi

} // end namespace Opm
