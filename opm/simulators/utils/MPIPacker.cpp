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
#include <opm/simulators/utils/MPIPacker.hpp>

#include <opm/common/utility/TimeService.hpp>
#include <opm/input/eclipse/EclipseState/IOConfig/FIPConfig.hpp>

#include <bitset>
#include <cstdint>
#include <ctime>
#include <string>
#include <type_traits>


namespace Opm {
namespace Mpi {
namespace detail {

template<std::size_t Size>
std::size_t Packing<false,std::bitset<Size>>::
packSize(const std::bitset<Size>& data,
         Parallel::MPIComm comm)
{
    return Packing<true,unsigned long long>::packSize(data.to_ullong(), comm);
}

template<std::size_t Size>
void Packing<false,std::bitset<Size>>::
pack(const std::bitset<Size>& data,
     std::vector<char>& buffer,
     int& position,
     Parallel::MPIComm comm)
{
    Packing<true,unsigned long long>::pack(data.to_ullong(), buffer, position, comm);
}

template<std::size_t Size>
void Packing<false,std::bitset<Size>>::
unpack(std::bitset<Size>& data,
       std::vector<char>& buffer,
       int& position,
       Parallel::MPIComm comm)
{
    unsigned long long d;
    Packing<true,unsigned long long>::unpack(d, buffer, position, comm);
    data = std::bitset<Size>(d);
}

std::size_t Packing<false,std::string>::
packSize(const std::string& data, Parallel::MPIComm comm)
{
    int size;
    MPI_Pack_size(1, Dune::MPITraits<std::size_t>::getType(), comm, &size);
    int totalSize = size;
    MPI_Pack_size(data.size(), MPI_CHAR, comm, &size);
    return totalSize + size;
}

void Packing<false,std::string>::
pack(const std::string& data,
     std::vector<char>& buffer,
     int& position,
     Parallel::MPIComm comm)
{
    std::size_t length = data.size();
    MPI_Pack(&length, 1, Dune::MPITraits<std::size_t>::getType(), buffer.data(),
             buffer.size(), &position, comm);
    MPI_Pack(data.data(), length, MPI_CHAR, buffer.data(), buffer.size(),
             &position, comm);
}

void Packing<false,std::string>::
unpack(std::string& data,
       std::vector<char>& buffer,
       int& position,
       Opm::Parallel::MPIComm comm)
{
    std::size_t length = 0;
    MPI_Unpack(buffer.data(), buffer.size(), &position, &length, 1,
               Dune::MPITraits<std::size_t>::getType(), comm);
    std::vector<char> cStr(length+1, '\0');
    MPI_Unpack(buffer.data(), buffer.size(), &position, cStr.data(), length,
               MPI_CHAR, comm);
    data.clear();
    data.append(cStr.data(), length);
}

std::size_t Packing<false,time_point>::
packSize(const time_point&, Opm::Parallel::MPIComm comm)
{
    return Packing<true,std::time_t>::packSize(std::time_t(), comm);
}

void Packing<false,time_point>::
pack(const time_point& data,
     std::vector<char>& buffer,
     int& position,
     Parallel::MPIComm comm)
{
    Packing<true,std::time_t>::pack(TimeService::to_time_t(data),
                                    buffer, position, comm);
}

void Packing<false,time_point>::
unpack(time_point& data,
       std::vector<char>& buffer,
       int& position,
       Parallel::MPIComm comm)
{
    std::time_t res;
    Packing<true,std::time_t>::unpack(res, buffer, position, comm);
    data = TimeService::from_time_t(res);
}

template struct Packing<false,std::bitset<3>>;
template struct Packing<false,std::bitset<4>>;
template struct Packing<false,std::bitset<10>>;
constexpr int NumFip = static_cast<int>(FIPConfig::OutputField::NUM_FIP_REPORT);
template struct Packing<false,std::bitset<NumFip>>;

} // end namespace detail
} // end namespace Mpi
} // end namespace Opm
