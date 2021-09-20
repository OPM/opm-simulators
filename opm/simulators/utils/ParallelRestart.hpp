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
#ifndef PARALLEL_RESTART_HPP
#define PARALLEL_RESTART_HPP

#if HAVE_MPI
#include <mpi.h>
#endif

#include <opm/common/ErrorMacros.hpp>
#include <opm/common/utility/TimeService.hpp>

#include <dune/common/version.hh>
#include <dune/common/parallel/mpihelper.hh>

#include <chrono>
#include <optional>
#include <map>
#include <set>
#include <string>
#include <tuple>
#include <typeinfo>
#include <unordered_map>
#include <unordered_set>
#include <vector>

using MPIComm = typename Dune::MPIHelper::MPICommunicator;
#if DUNE_VERSION_NEWER(DUNE_COMMON, 2, 7)
    using Communication = Dune::Communication<MPIComm>; 
#else
    using Communication = Dune::CollectiveCommunication<MPIComm>;
#endif

namespace Opm
{

class EclipseIO;
class SummaryState;
class RestartKey;
class RestartValue;

namespace data
{
struct AquiferData;
struct CarterTracyData;
struct CellData;
struct Connection;
struct CurrentControl;
struct FetkovichData;
class GroupAndNetworkValues;
struct GroupConstraints;
struct GroupData;
struct GroupGuideRates;
class GuideRateValue;
struct NodeData;
struct NumericAquiferData;
class Rates;
struct Segment;
class Solution;
struct Well;
class Wells;
}

namespace Action
{
class State;
}

namespace Mpi
{
template<class T>
std::size_t packSize(const T*, std::size_t, MPIComm,
                     std::integral_constant<bool, false>);

template<class T>
std::size_t packSize(const T*, std::size_t l, MPIComm comm,
                     std::integral_constant<bool, true>);

template<class T>
std::size_t packSize(const T* data, std::size_t l, MPIComm comm);

template<class T>
std::size_t packSize(const T&, MPIComm,
                     std::integral_constant<bool, false>)
{
    std::string msg = std::string{"Packing not (yet) supported for non-pod type: "} + typeid(T).name();
    OPM_THROW(std::logic_error, msg);
}

template<class T>
std::size_t packSize(const T&, MPIComm comm,
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
std::size_t packSize(const T& data, MPIComm comm)
{
    return packSize(data, comm, typename std::is_pod<T>::type());
}

template<class T1, class T2>
std::size_t packSize(const std::pair<T1,T2>& data, MPIComm comm);

template<class T>
std::size_t packSize(const std::optional<T>& data, MPIComm comm);

template<class T, class A>
std::size_t packSize(const std::vector<T,A>& data, MPIComm comm);

template<class K, class C, class A>
std::size_t packSize(const std::set<K,C,A>& data,
                     MPIComm comm);

template<class T, class H, class KE, class A>
std::size_t packSize(const std::unordered_set<T,H,KE,A>& data,
                     MPIComm comm);

template<class A>
std::size_t packSize(const std::vector<bool,A>& data, MPIComm comm);

template<class... Ts>
std::size_t packSize(const std::tuple<Ts...>& data, MPIComm comm);

template<class T, std::size_t N>
std::size_t packSize(const std::array<T,N>& data, MPIComm comm);

std::size_t packSize(const char* str, MPIComm comm);

template<class T1, class T2, class C, class A>
std::size_t packSize(const std::map<T1,T2,C,A>& data, MPIComm comm);

template<class T1, class T2, class H, class P, class A>
std::size_t packSize(const std::unordered_map<T1,T2,H,P,A>& data, MPIComm comm);

////// pack routines

template<class T>
void pack(const T*, std::size_t, std::vector<char>&, int&,
          MPIComm, std::integral_constant<bool, false>);

template<class T>
void pack(const T* data, std::size_t l, std::vector<char>& buffer, int& position,
          MPIComm comm, std::integral_constant<bool, true>);

template<class T>
void pack(const T* data, std::size_t l, std::vector<char>& buffer, int& position,
          MPIComm comm);

template<class T>
void pack(const T&, std::vector<char>&, int&,
          MPIComm, std::integral_constant<bool, false>)
{
    OPM_THROW(std::logic_error, "Packing not (yet) supported for this non-pod type.");
}

template<class T>
void pack(const T& data, std::vector<char>& buffer, int& position,
          MPIComm comm, std::integral_constant<bool, true>)
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
          MPIComm comm)
{
    pack(data, buffer, position, comm, typename std::is_pod<T>::type());
}

template<class T1, class T2>
void pack(const std::pair<T1,T2>& data, std::vector<char>& buffer, int& position,
          MPIComm comm);

template<class T>
void pack(const std::optional<T>& data, std::vector<char>& buffer, int& position,
          MPIComm comm);

template<class T, class A>
void pack(const std::vector<T,A>& data, std::vector<char>& buffer, int& position,
          MPIComm comm);

template<class A>
void pack(const std::vector<bool,A>& data, std::vector<char>& buffer, int& position,
          MPIComm comm);

template<class... Ts>
void pack(const std::tuple<Ts...>& data, std::vector<char>& buffer,
          int& position, MPIComm comm);

template<class K, class C, class A>
void pack(const std::set<K,C,A>& data,
          std::vector<char>& buffer, int& position,
          MPIComm comm);

template<class T, class H, class KE, class A>
void pack(const std::unordered_set<T,H,KE,A>& data,
          std::vector<char>& buffer, int& position,
          MPIComm comm);

template<class T, size_t N>
void pack(const std::array<T,N>& data, std::vector<char>& buffer, int& position,
          MPIComm comm);

template<class T1, class T2, class C, class A>
void pack(const std::map<T1,T2,C,A>& data, std::vector<char>& buffer, int& position,
          MPIComm comm);

template<class T1, class T2, class H, class P, class A>
void pack(const std::unordered_map<T1,T2,H,P,A>& data, std::vector<char>& buffer, int& position,
          MPIComm comm);

void pack(const char* str, std::vector<char>& buffer, int& position,
          MPIComm comm);

/// unpack routines

template<class T>
void unpack(T*, const std::size_t&, std::vector<char>&, int&,
            MPIComm, std::integral_constant<bool, false>);

template<class T>
void unpack(T* data, const std::size_t& l, std::vector<char>& buffer, int& position,
            MPIComm comm,
            std::integral_constant<bool, true>);

template<class T>
void unpack(T* data, const std::size_t& l, std::vector<char>& buffer, int& position,
            MPIComm comm);

template<class T>
void unpack(T&, std::vector<char>&, int&,
            MPIComm, std::integral_constant<bool, false>)
{
    OPM_THROW(std::logic_error, "Packing not (yet) supported for this non-pod type.");
}

template<class T>
void unpack(T& data, std::vector<char>& buffer, int& position,
            MPIComm comm, std::integral_constant<bool, true>)
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
            MPIComm comm)
{
    unpack(data, buffer, position, comm, typename std::is_pod<T>::type());
}

template<class T1, class T2>
void unpack(std::pair<T1,T2>& data, std::vector<char>& buffer, int& position,
            MPIComm comm);

template<class T>
void unpack(std::optional<T>& data, std::vector<char>& buffer, int& position,
            MPIComm comm);

template<class T, class A>
void unpack(std::vector<T,A>& data, std::vector<char>& buffer, int& position,
            MPIComm comm);

template<class A>
void unpack(std::vector<bool,A>& data, std::vector<char>& buffer, int& position,
          MPIComm comm);

template<class... Ts>
void unpack(std::tuple<Ts...>& data, std::vector<char>& buffer,
            int& position, MPIComm comm);

template<class K, class C, class A>
void unpack(std::set<K,C,A>& data,
            std::vector<char>& buffer, int& position,
            MPIComm comm);

template<class T, class H, class KE, class A>
void unpack(std::unordered_set<T,H,KE,A>& data,
            std::vector<char>& buffer, int& position,
            MPIComm comm);

template<class T, size_t N>
void unpack(std::array<T,N>& data, std::vector<char>& buffer, int& position,
          MPIComm comm);

template<class T1, class T2, class C, class A>
void unpack(std::map<T1,T2,C,A>& data, std::vector<char>& buffer, int& position,
            MPIComm comm);

template<class T1, class T2, class H, class P, class A>
void unpack(std::unordered_map<T1,T2,H,P,A>& data, std::vector<char>& buffer, int& position,
            MPIComm comm);

void unpack(char* str, std::size_t length, std::vector<char>& buffer, int& position,
            MPIComm comm);

/// prototypes for complex types

#define ADD_PACK_PROTOTYPES(T) \
  std::size_t packSize(const T& data, MPIComm comm); \
  void pack(const T& data, std::vector<char>& buffer, int& position, \
          MPIComm comm); \
  void unpack(T& data, std::vector<char>& buffer, int& position, \
              MPIComm comm);

ADD_PACK_PROTOTYPES(data::AquiferData)
ADD_PACK_PROTOTYPES(data::CarterTracyData)
ADD_PACK_PROTOTYPES(data::CellData)
ADD_PACK_PROTOTYPES(data::Connection)
ADD_PACK_PROTOTYPES(data::CurrentControl)
ADD_PACK_PROTOTYPES(data::FetkovichData)
ADD_PACK_PROTOTYPES(data::Rates)
ADD_PACK_PROTOTYPES(data::Segment)
ADD_PACK_PROTOTYPES(data::Solution)
ADD_PACK_PROTOTYPES(data::GuideRateValue)
ADD_PACK_PROTOTYPES(data::GroupConstraints)
ADD_PACK_PROTOTYPES(data::GroupGuideRates)
ADD_PACK_PROTOTYPES(data::GroupData)
ADD_PACK_PROTOTYPES(data::NodeData)
ADD_PACK_PROTOTYPES(data::GroupAndNetworkValues)
ADD_PACK_PROTOTYPES(data::NumericAquiferData)
ADD_PACK_PROTOTYPES(data::Well)
ADD_PACK_PROTOTYPES(data::Wells)
ADD_PACK_PROTOTYPES(RestartKey)
ADD_PACK_PROTOTYPES(RestartValue)
ADD_PACK_PROTOTYPES(std::string)
ADD_PACK_PROTOTYPES(time_point)

} // end namespace Mpi

RestartValue loadParallelRestart(const EclipseIO* eclIO, Action::State& actionState, SummaryState& summaryState,
                                 const std::vector<RestartKey>& solutionKeys,
                                 const std::vector<RestartKey>& extraKeys,
                                 Communication comm);

} // end namespace Opm
#endif // PARALLEL_RESTART_HPP
