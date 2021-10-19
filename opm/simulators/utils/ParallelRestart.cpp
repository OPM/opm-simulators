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
#if HAVE_MPI
#include <mpi.h>
#endif

#include "ParallelRestart.hpp"
#include <cassert>
#include <cstdint>
#include <cstring>
#include <ctime>
#include <memory>
#include <dune/common/parallel/mpitraits.hh>
#include <opm/output/data/Aquifer.hpp>
#include <opm/output/data/Cells.hpp>
#include <opm/output/data/Groups.hpp>
#include <opm/output/data/GuideRateValue.hpp>
#include <opm/output/data/Solution.hpp>
#include <opm/output/data/Wells.hpp>
#include <opm/output/eclipse/EclipseIO.hpp>
#include <opm/output/eclipse/RestartValue.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/SummaryState.hpp>

#define HANDLE_AS_POD(T) \
  std::size_t packSize(const T& data, Opm::Parallel::MPIComm comm) \
  { \
      return packSize(data, comm, std::integral_constant<bool,true>()); \
  } \
  void pack(const T& data, std::vector<char>& buffer, int& position, \
            Opm::Parallel::MPIComm comm) \
  { \
      pack(data, buffer, position, comm, std::integral_constant<bool,true>()); \
  } \
  void unpack(T& data, std::vector<char>& buffer, int& position, \
              Opm::Parallel::MPIComm comm) \
  { \
      unpack(data, buffer, position, comm, std::integral_constant<bool,true>()); \
  }

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

template<class T1, class T2>
std::size_t packSize(const std::pair<T1,T2>& data, Opm::Parallel::MPIComm comm)
{
    return packSize(data.first, comm) + packSize(data.second, comm);
}

template<class T>
std::size_t packSize(const std::optional<T>& data, Opm::Parallel::MPIComm comm)
{
    bool has_value = data.has_value();
    std::size_t pack_size = packSize(has_value, comm);
    if (has_value)
        pack_size += packSize(*data, comm);
    return pack_size;
}


template<class T, class A>
std::size_t packSize(const std::vector<T,A>& data, Opm::Parallel::MPIComm comm)
{
    if (std::is_pod<T>::value)
        // size written automatically
        return packSize(data.data(), data.size(), comm);

    std::size_t size = packSize(data.size(), comm);

    for (const auto& entry: data)
        size += packSize(entry, comm);

    return size;
}

template<class A>
std::size_t packSize(const std::vector<bool,A>& data, Opm::Parallel::MPIComm comm)
{
    bool entry;
    return packSize(data.size(), comm) + data.size()*packSize(entry,comm);
}

template<std::size_t I = 0, typename Tuple>
typename std::enable_if<I == std::tuple_size<Tuple>::value, std::size_t>::type
pack_size_tuple_entry(const Tuple&, Opm::Parallel::MPIComm)
{
    return 0;
}

template<std::size_t I = 0, typename Tuple>
typename std::enable_if<I != std::tuple_size<Tuple>::value, std::size_t>::type
pack_size_tuple_entry(const Tuple& tuple, Opm::Parallel::MPIComm comm)
{
    return packSize(std::get<I>(tuple), comm) + pack_size_tuple_entry<I+1>(tuple, comm);
}

template<class... Ts>
std::size_t packSize(const std::tuple<Ts...>& data, Opm::Parallel::MPIComm comm)
{
    return pack_size_tuple_entry(data, comm);
}

template<class T, class H, class KE, class A>
std::size_t packSize(const std::unordered_set<T,H,KE,A>& data,
                     Opm::Parallel::MPIComm comm)
{
    std::size_t totalSize = packSize(data.size(), comm);
    for (const auto& entry : data)
    {
        totalSize += packSize(entry, comm);
    }
    return totalSize;
}

template<class K, class C, class A>
std::size_t packSize(const std::set<K,C,A>& data,
                     Opm::Parallel::MPIComm comm)
{
    std::size_t totalSize = packSize(data.size(), comm);
    for (const auto& entry : data)
    {
        totalSize += packSize(entry, comm);
    }
    return totalSize;
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

template<class T1, class T2, class C, class A>
std::size_t packSize(const std::map<T1,T2,C,A>& data, Opm::Parallel::MPIComm comm)
{
    std::size_t totalSize = packSize(data.size(), comm);
    for (const auto& entry: data)
    {
        totalSize += packSize(entry, comm);
    }
    return totalSize;
}

template<class T1, class T2, class H, class P, class A>
std::size_t packSize(const std::unordered_map<T1,T2,H,P,A>& data, Opm::Parallel::MPIComm comm)
{
    std::size_t totalSize = packSize(data.size(), comm);
    for (const auto& entry: data)
    {
        totalSize += packSize(entry, comm);
    }
    return totalSize;
}

template<class T, std::size_t N>
std::size_t packSize(const std::array<T,N>& data, Opm::Parallel::MPIComm comm)
{
    return N*packSize(data[0], comm);
}

HANDLE_AS_POD(data::CarterTracyData)
HANDLE_AS_POD(data::Connection)
HANDLE_AS_POD(data::CurrentControl)
HANDLE_AS_POD(data::FetkovichData)
HANDLE_AS_POD(data::GroupConstraints)
HANDLE_AS_POD(data::NodeData)
HANDLE_AS_POD(data::Rates)
HANDLE_AS_POD(data::Segment)

std::size_t packSize(const data::NumericAquiferData& data, Opm::Parallel::MPIComm comm)
{
    return packSize(data.initPressure, comm);
}

std::size_t packSize(const data::AquiferData& data, Opm::Parallel::MPIComm comm)
{
    const auto type = 0ull;

    const auto base = packSize(data.aquiferID, comm)
        + packSize(data.pressure, comm)
        + packSize(data.fluxRate, comm)
        + packSize(data.volume, comm)
        + packSize(data.initPressure, comm)
        + packSize(data.datumDepth, comm)
        + packSize(type, comm);

    if (auto const* aquFet = data.typeData.get<data::AquiferType::Fetkovich>();
        aquFet != nullptr)
    {
        return base + packSize(*aquFet, comm);
    }
    else if (auto const* aquCT = data.typeData.get<data::AquiferType::CarterTracy>();
             aquCT != nullptr)
    {
        return base + packSize(*aquCT, comm);
    }
    else if (auto const* aquNum = data.typeData.get<data::AquiferType::Numerical>();
             aquNum != nullptr)
    {
        return base + packSize(*aquNum, comm);
    }

    return base;
}

std::size_t packSize(const data::GuideRateValue&, Opm::Parallel::MPIComm comm)
{
    const auto nItem = static_cast<std::size_t>(data::GuideRateValue::Item::NumItems);

    return packSize(std::array<int   , nItem>{}, comm)
        +  packSize(std::array<double, nItem>{}, comm);
}

std::size_t packSize(const data::GroupGuideRates& data, Opm::Parallel::MPIComm comm)
{
    return packSize(data.production, comm)
        +  packSize(data.injection, comm);
}

std::size_t packSize(const data::GroupData& data, Opm::Parallel::MPIComm comm)
{
    return packSize(data.currentControl, comm)
        +  packSize(data.guideRates, comm);
}

std::size_t packSize(const data::Well& data, Opm::Parallel::MPIComm comm)
{
    std::size_t size = packSize(data.rates, comm);
    size += packSize(data.bhp, comm) + packSize(data.thp, comm);
    size += packSize(data.temperature, comm);
    size += packSize(data.control, comm);
    size += packSize(data.connections, comm);
    size += packSize(data.segments, comm);
    size += packSize(data.current_control, comm);
    size += packSize(data.guide_rates, comm);
    return size;
}

std::size_t packSize(const data::CellData& data, Opm::Parallel::MPIComm comm)
{
    return packSize(data.dim, comm) + packSize(data.data, comm) + packSize(data.target, comm);
}

std::size_t packSize(const RestartKey& data, Opm::Parallel::MPIComm comm)
{
    return packSize(data.key, comm) + packSize(data.dim, comm) + packSize(data.required, comm);
}

std::size_t packSize(const data::Solution& data, Opm::Parallel::MPIComm comm)
{
    // Needs explicit conversion to a supported base type holding the data
    // to prevent throwing.
    return packSize(static_cast<const std::map< std::string, data::CellData>&>(data), comm);
}

std::size_t packSize(const data::GroupAndNetworkValues& data, Opm::Parallel::MPIComm comm)
{
    return packSize(data.groupData, comm)
        +  packSize(data.nodeData, comm);
}

std::size_t packSize(const data::Wells& data, Opm::Parallel::MPIComm comm)
{
    // Needs explicit conversion to a supported base type holding the data
    // to prevent throwing.
    return packSize(static_cast<const std::map< std::string, data::Well>&>(data), comm);
}

std::size_t packSize(const RestartValue& data, Opm::Parallel::MPIComm comm)
{
    return packSize(data.solution, comm)
        +  packSize(data.wells, comm)
        +  packSize(data.grp_nwrk, comm)
        +  packSize(data.aquifer, comm)
        +  packSize(data.extra, comm);
}

std::size_t packSize(const Opm::time_point&, Opm::Parallel::MPIComm comm)
{
    std::time_t tp;
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

template<class T1, class T2>
void pack(const std::pair<T1,T2>& data, std::vector<char>& buffer, int& position,
          Opm::Parallel::MPIComm comm)
{
    pack(data.first, buffer, position, comm);
    pack(data.second, buffer, position, comm);
}

template<class T>
void pack(const std::optional<T>& data, std::vector<char>& buffer, int& position,
    Opm::Parallel::MPIComm comm)
{
    bool has_value = data.has_value();
    pack(has_value, buffer, position, comm);
    if (has_value)
        pack(*data, buffer, position, comm);
}


template<class T, class A>
void pack(const std::vector<T, A>& data, std::vector<char>& buffer, int& position,
          Opm::Parallel::MPIComm comm)
{
    if (std::is_pod<T>::value)
    {
        // size written automatically
        pack(data.data(), data.size(), buffer, position, comm);
        return;
    }

    pack(data.size(), buffer, position, comm);

    for (const auto& entry: data)
        pack(entry, buffer, position, comm);
}

template<class K, class C, class A>
void pack(const std::set<K,C,A>& data,
          std::vector<char>& buffer, int& position,
          Opm::Parallel::MPIComm comm)
{
    pack(data.size(), buffer, position, comm);

    for (const auto& entry : data)
    {
        pack(entry, buffer, position, comm);
    }
}

template<class T, class H, class KE, class A>
void pack(const std::unordered_set<T,H,KE,A>& data,
          std::vector<char>& buffer, int& position,
          Opm::Parallel::MPIComm comm)
{
    pack(data.size(), buffer, position, comm);

    for (const auto& entry : data)
    {
        pack(entry, buffer, position, comm);
    }
}

template<class T, size_t N>
void pack(const std::array<T,N>& data, std::vector<char>& buffer, int& position,
          Opm::Parallel::MPIComm comm)
{
    for (const T& entry : data)
        pack(entry, buffer, position, comm);
}

template<class A>
void pack(const std::vector<bool,A>& data, std::vector<char>& buffer, int& position,
          Opm::Parallel::MPIComm comm)
{
    pack(data.size(), buffer, position, comm);
    for (const auto entry : data) { // Not a reference: vector<bool> range
        bool b = entry;
        pack(b, buffer, position, comm);
    }
}

template<std::size_t I = 0, typename Tuple>
typename std::enable_if<I == std::tuple_size<Tuple>::value, void>::type
pack_tuple_entry(const Tuple&, std::vector<char>&, int&,
                      Opm::Parallel::MPIComm)
{
}

template<std::size_t I = 0, typename Tuple>
typename std::enable_if<I != std::tuple_size<Tuple>::value, void>::type
pack_tuple_entry(const Tuple& tuple, std::vector<char>& buffer,
                 int& position, Opm::Parallel::MPIComm comm)
{
    pack(std::get<I>(tuple), buffer, position, comm);
    pack_tuple_entry<I+1>(tuple, buffer, position, comm);
}

template<class... Ts>
void pack(const std::tuple<Ts...>& data, std::vector<char>& buffer,
          int& position, Opm::Parallel::MPIComm comm)
{
    pack_tuple_entry(data, buffer, position, comm);
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

template<class T1, class T2, class C, class A>
void pack(const std::map<T1,T2,C,A>& data, std::vector<char>& buffer, int& position,
          Opm::Parallel::MPIComm comm)
{
    pack(data.size(), buffer, position, comm);

    for (const auto& entry: data)
    {
        pack(entry, buffer, position, comm);
    }
}

template<class T1, class T2, class H, class P, class A>
void pack(const std::unordered_map<T1,T2,H,P,A>& data, std::vector<char>& buffer, int& position,
          Opm::Parallel::MPIComm comm)
{
    pack(data.size(), buffer, position, comm);

    for (const auto& entry: data)
    {
        pack(entry, buffer, position, comm);
    }
}

void pack(const data::NumericAquiferData& data, std::vector<char>& buffer, int& position,
          Opm::Parallel::MPIComm comm)
{
    pack(data.initPressure, buffer, position, comm);
}

void pack(const data::AquiferData& data, std::vector<char>& buffer, int& position,
          Opm::Parallel::MPIComm comm)
{
    const auto type =
          (data.typeData.is<data::AquiferType::Fetkovich>()   * (1ull << 0))
        + (data.typeData.is<data::AquiferType::CarterTracy>() * (1ull << 1))
        + (data.typeData.is<data::AquiferType::Numerical>()   * (1ull << 2));

    pack(data.aquiferID, buffer, position, comm);
    pack(data.pressure, buffer, position, comm);
    pack(data.fluxRate, buffer, position, comm);
    pack(data.volume, buffer, position, comm);
    pack(data.initPressure, buffer, position, comm);
    pack(data.datumDepth, buffer, position, comm);
    pack(type, buffer, position, comm);

    if (auto const* aquFet = data.typeData.get<data::AquiferType::Fetkovich>();
        aquFet != nullptr)
    {
        pack(*aquFet, buffer, position, comm);
    }
    else if (auto const* aquCT = data.typeData.get<data::AquiferType::CarterTracy>();
             aquCT != nullptr)
    {
        pack(*aquCT, buffer, position, comm);
    }
    else if (auto const* aquNum = data.typeData.get<data::AquiferType::Numerical>();
             aquNum != nullptr)
    {
        pack(*aquNum, buffer, position, comm);
    }
}

void pack(const data::GuideRateValue& data, std::vector<char>& buffer, int& position,
          Opm::Parallel::MPIComm comm)
{
    using Item = data::GuideRateValue::Item;
    const auto nItem = static_cast<std::size_t>(Item::NumItems);

    auto has = std::array<int   , nItem>{};  has.fill(0);
    auto val = std::array<double, nItem>{};  val.fill(0.0);

    for (auto itemID = 0*nItem; itemID < nItem; ++itemID) {
        const auto item = static_cast<Item>(itemID);

        if (data.has(item)) {
            has[itemID] = 1;
            val[itemID] = data.get(item);
        }
    }

    pack(has, buffer, position, comm);
    pack(val, buffer, position, comm);
}

void pack(const data::GroupGuideRates& data, std::vector<char>& buffer, int& position,
          Opm::Parallel::MPIComm comm)
{
    pack(data.production, buffer, position, comm);
    pack(data.injection, buffer, position, comm);
}

void pack(const data::GroupData& data, std::vector<char>& buffer, int& position,
          Opm::Parallel::MPIComm comm)
{
    pack(data.currentControl, buffer, position, comm);
    pack(data.guideRates, buffer, position, comm);
}

void pack(const data::Well& data, std::vector<char>& buffer, int& position,
          Opm::Parallel::MPIComm comm)
{
    pack(data.rates, buffer, position, comm);
    pack(data.bhp, buffer, position, comm);
    pack(data.thp, buffer, position, comm);
    pack(data.temperature, buffer, position, comm);
    pack(data.control, buffer, position, comm);
    pack(data.connections, buffer, position, comm);
    pack(data.segments, buffer, position, comm);
    pack(data.current_control, buffer, position, comm);
    pack(data.guide_rates, buffer, position, comm);
}

void pack(const RestartKey& data, std::vector<char>& buffer, int& position,
          Opm::Parallel::MPIComm comm)
{
    pack(data.key, buffer, position, comm);
    pack(data.dim, buffer, position, comm);
    pack(data.required, buffer, position, comm);
}

void pack(const data::CellData& data, std::vector<char>& buffer, int& position,
          Opm::Parallel::MPIComm comm)
{
    pack(data.dim, buffer, position, comm);
    pack(data.data, buffer, position, comm);
    pack(data.target, buffer, position, comm);
}

void pack(const data::Solution& data, std::vector<char>& buffer, int& position,
          Opm::Parallel::MPIComm comm)
{
    // Needs explicit conversion to a supported base type holding the data
    // to prevent throwing.
    pack(static_cast<const std::map< std::string, data::CellData>&>(data),
         buffer, position, comm);
}

void pack(const data::Wells& data, std::vector<char>& buffer, int& position,
          Opm::Parallel::MPIComm comm)
{
    // Needs explicit conversion to a supported base type holding the data
    // to prevent throwing.
    pack(static_cast<const std::map< std::string, data::Well>&>(data),
         buffer, position, comm);
}

void pack(const data::GroupAndNetworkValues& data, std::vector<char>& buffer, int& position,
          Opm::Parallel::MPIComm comm)
{
    pack(data.groupData, buffer, position, comm);
    pack(data.nodeData, buffer, position, comm);
}

void pack(const RestartValue& data, std::vector<char>& buffer, int& position,
          Opm::Parallel::MPIComm comm)
{
    pack(data.solution, buffer, position, comm);
    pack(data.wells, buffer, position, comm);
    pack(data.grp_nwrk, buffer, position, comm);
    pack(data.aquifer, buffer, position, comm);
    pack(data.extra, buffer, position, comm);
}

void pack(const Opm::time_point& data, std::vector<char>& buffer, int& position,
          Opm::Parallel::MPIComm comm)
{
    pack(Opm::TimeService::to_time_t(data), buffer, position, comm);
}


/// unpack routines

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

template<class T1, class T2>
void unpack(std::pair<T1,T2>& data, std::vector<char>& buffer, int& position,
            Opm::Parallel::MPIComm comm)
{
    unpack(data.first, buffer, position, comm);
    unpack(data.second, buffer, position, comm);
}

template<class T>
void unpack(std::optional<T>&data, std::vector<char>& buffer, int& position,
            Opm::Parallel::MPIComm comm)
{
    bool has_value;
    unpack(has_value, buffer, position, comm);
    if (has_value) {
        T val;
        unpack(val, buffer, position, comm);
        data = std::optional<T>(val);
    } else
        data.reset();
}


template<class T, class A>
void unpack(std::vector<T,A>& data, std::vector<char>& buffer, int& position,
            Opm::Parallel::MPIComm comm)
{
    std::size_t length = 0;
    unpack(length, buffer, position, comm);
    data.resize(length);

    if (std::is_pod<T>::value)
    {
        unpack(data.data(), data.size(), buffer, position, comm);
        return;
    }

    for (auto& entry: data)
        unpack(entry, buffer, position, comm);
}

template<class A>
void unpack(std::vector<bool,A>& data, std::vector<char>& buffer, int& position,
            Opm::Parallel::MPIComm comm)
{
    size_t size;
    unpack(size, buffer, position, comm);
    data.clear();
    data.reserve(size);
    for (size_t i = 0; i < size; ++i) {
        bool entry;
        unpack(entry, buffer, position, comm);
        data.push_back(entry);
    }
}

template<std::size_t I = 0, typename Tuple>
typename std::enable_if<I == std::tuple_size<Tuple>::value, void>::type
unpack_tuple_entry(Tuple&, std::vector<char>&, int&,
                   Opm::Parallel::MPIComm)
{
}

template<std::size_t I = 0, typename Tuple>
typename std::enable_if<I != std::tuple_size<Tuple>::value, void>::type
unpack_tuple_entry(Tuple& tuple, std::vector<char>& buffer,
                   int& position, Opm::Parallel::MPIComm comm)
{
    unpack(std::get<I>(tuple), buffer, position, comm);
    unpack_tuple_entry<I+1>(tuple, buffer, position, comm);
}

template<class... Ts>
void unpack(std::tuple<Ts...>& data, std::vector<char>& buffer,
            int& position, Opm::Parallel::MPIComm comm)
{
    unpack_tuple_entry(data, buffer, position, comm);
}

template<class K, class C, class A>
void unpack(std::set<K,C,A>& data,
            std::vector<char>& buffer, int& position,
            Opm::Parallel::MPIComm comm)
{
    std::size_t size = 0;
    unpack(size, buffer, position, comm);

    for (;size>0; size--)
    {
        K entry;
        unpack(entry, buffer, position, comm);
        data.insert(entry);
    }
}

template<class T, class H, class KE, class A>
void unpack(std::unordered_set<T,H,KE,A>& data,
            std::vector<char>& buffer, int& position,
            Opm::Parallel::MPIComm comm)
{
    std::size_t size=0;
    unpack(size, buffer, position, comm);

    for (;size>0; size--)
    {
        T entry;
        unpack(entry, buffer, position, comm);
        data.insert(entry);
    }
}

template<class T, size_t N>
void unpack(std::array<T,N>& data, std::vector<char>& buffer, int& position,
          Opm::Parallel::MPIComm comm)
{
    for (T& entry : data)
        unpack(entry, buffer, position, comm);
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

template<class T1, class T2, class C, class A>
void unpack(std::map<T1,T2,C,A>& data, std::vector<char>& buffer, int& position,
            Opm::Parallel::MPIComm comm)
{
    std::size_t size=0;
    unpack(size, buffer, position, comm);

    for (;size>0; size--)
    {
        std::pair<T1,T2> entry;
        unpack(entry, buffer, position, comm);
        data.insert(entry);
    }
}

template<class T1, class T2, class H, class P, class A>
void unpack(std::unordered_map<T1,T2,H,P,A>& data, std::vector<char>& buffer, int& position,
            Opm::Parallel::MPIComm comm)
{
    std::size_t size=0;
    unpack(size, buffer, position, comm);

    for (;size>0; size--)
    {
        std::pair<T1,T2> entry;
        unpack(entry, buffer, position, comm);
        data.insert(entry);
    }
}

void unpack(data::Well& data, std::vector<char>& buffer, int& position,
            Opm::Parallel::MPIComm comm)
{
    unpack(data.rates, buffer, position, comm);
    unpack(data.bhp, buffer, position, comm);
    unpack(data.thp, buffer, position, comm);
    unpack(data.temperature, buffer, position, comm);
    unpack(data.control, buffer, position, comm);
    unpack(data.connections, buffer, position, comm);
    unpack(data.segments, buffer, position, comm);
    unpack(data.current_control, buffer, position, comm);
    unpack(data.guide_rates, buffer, position, comm);
}

void unpack(data::NumericAquiferData& data, std::vector<char>& buffer, int& position,
            Opm::Parallel::MPIComm comm)
{
    unpack(data.initPressure, buffer, position, comm);
}

void unpack(data::AquiferData& data, std::vector<char>& buffer, int& position,
            Opm::Parallel::MPIComm comm)
{
    auto type = 0ull;

    unpack(data.aquiferID, buffer, position, comm);
    unpack(data.pressure, buffer, position, comm);
    unpack(data.fluxRate, buffer, position, comm);
    unpack(data.volume, buffer, position, comm);
    unpack(data.initPressure, buffer, position, comm);
    unpack(data.datumDepth, buffer, position, comm);
    unpack(type, buffer, position, comm);

    if (type == (1ull << 0)) {
        auto* aquFet = data.typeData.create<data::AquiferType::Fetkovich>();
        unpack(*aquFet, buffer, position, comm);
    }
    else if (type == (1ull << 1)) {
        auto* aquCT = data.typeData.create<data::AquiferType::CarterTracy>();
        unpack(*aquCT, buffer, position, comm);
    }
    else if (type == (1ull << 2)) {
        auto* aquNum = data.typeData.create<data::AquiferType::Numerical>();
        unpack(*aquNum, buffer, position, comm);
    }
}

void unpack(data::GuideRateValue& data, std::vector<char>& buffer, int& position,
            Opm::Parallel::MPIComm comm)
{
    using Item = data::GuideRateValue::Item;
    const auto nItem = static_cast<std::size_t>(Item::NumItems);

    auto has = std::array<int   , nItem>{};
    auto val = std::array<double, nItem>{};

    unpack(has, buffer, position, comm);
    unpack(val, buffer, position, comm);

    for (auto itemID = 0*nItem; itemID < nItem; ++itemID) {
        if (has[itemID] != 0) {
            data.set(static_cast<Item>(itemID), val[itemID]);
        }
    }
}

void unpack(data::GroupGuideRates& data, std::vector<char>& buffer, int& position,
            Opm::Parallel::MPIComm comm)
{
    unpack(data.production, buffer, position, comm);
    unpack(data.injection, buffer, position, comm);
}

void unpack(data::GroupData& data, std::vector<char>& buffer, int& position,
            Opm::Parallel::MPIComm comm)
{
    unpack(data.currentControl, buffer, position, comm);
    unpack(data.guideRates, buffer, position, comm);
}

void unpack(RestartKey& data, std::vector<char>& buffer, int& position,
            Opm::Parallel::MPIComm comm)
{
    unpack(data.key, buffer, position, comm);
    unpack(data.dim, buffer, position, comm);
    unpack(data.required, buffer, position, comm);
}

void unpack(data::CellData& data, std::vector<char>& buffer, int& position,
            Opm::Parallel::MPIComm comm)
{
    unpack(data.dim, buffer, position, comm);
    unpack(data.data, buffer, position, comm);
    unpack(data.target, buffer, position, comm);
}

void unpack(data::Solution& data, std::vector<char>& buffer, int& position,
            Opm::Parallel::MPIComm comm)
{
    // Needs explicit conversion to a supported base type holding the data
    // to prevent throwing.
    unpack(static_cast<std::map< std::string, data::CellData>&>(data),
           buffer, position, comm);
}

void unpack(data::Wells& data, std::vector<char>& buffer, int& position,
            Opm::Parallel::MPIComm comm)
{
    // Needs explicit conversion to a supported base type holding the data
    // to prevent throwing.
    unpack(static_cast<std::map< std::string, data::Well>&>(data),
           buffer, position, comm);
}

void unpack(data::GroupAndNetworkValues& data, std::vector<char>& buffer, int& position,
            Opm::Parallel::MPIComm comm)
{
    unpack(data.groupData, buffer, position, comm);
    unpack(data.nodeData, buffer, position, comm);
}

void unpack(RestartValue& data, std::vector<char>& buffer, int& position,
            Opm::Parallel::MPIComm comm)
{
    unpack(data.solution, buffer, position, comm);
    unpack(data.wells, buffer, position, comm);
    unpack(data.grp_nwrk, buffer, position, comm);
    unpack(data.aquifer, buffer, position, comm);
    unpack(data.extra, buffer, position, comm);
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


#define INSTANTIATE_PACK_VECTOR(...) \
template std::size_t packSize(const std::vector<__VA_ARGS__>& data, \
                              Opm::Parallel::MPIComm comm); \
template void pack(const std::vector<__VA_ARGS__>& data, \
                   std::vector<char>& buffer, int& position, \
                   Opm::Parallel::MPIComm comm); \
template void unpack(std::vector<__VA_ARGS__>& data, \
                     std::vector<char>& buffer, int& position, \
                     Opm::Parallel::MPIComm comm);

INSTANTIATE_PACK_VECTOR(float)
INSTANTIATE_PACK_VECTOR(double)
INSTANTIATE_PACK_VECTOR(std::vector<double>)
INSTANTIATE_PACK_VECTOR(bool)
INSTANTIATE_PACK_VECTOR(char)
INSTANTIATE_PACK_VECTOR(int)
INSTANTIATE_PACK_VECTOR(unsigned char)
INSTANTIATE_PACK_VECTOR(unsigned int)
INSTANTIATE_PACK_VECTOR(unsigned long int)
INSTANTIATE_PACK_VECTOR(unsigned long long int)
INSTANTIATE_PACK_VECTOR(std::time_t)
INSTANTIATE_PACK_VECTOR(std::array<double, 3>)
INSTANTIATE_PACK_VECTOR(std::pair<bool,double>)
INSTANTIATE_PACK_VECTOR(std::map<std::string,int>)
INSTANTIATE_PACK_VECTOR(std::pair<std::string,std::vector<size_t>>)
INSTANTIATE_PACK_VECTOR(std::pair<int,std::vector<int>>)
INSTANTIATE_PACK_VECTOR(std::pair<int,std::vector<size_t>>)
INSTANTIATE_PACK_VECTOR(std::string)

#undef INSTANTIATE_PACK_VECTOR

#undef INSTANTIATE_PACK_SET

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
INSTANTIATE_PACK(std::array<short,3>)
INSTANTIATE_PACK(std::array<bool,3>)
INSTANTIATE_PACK(std::array<int,3>)
INSTANTIATE_PACK(std::map<std::pair<int,int>,std::pair<bool,double>>)
INSTANTIATE_PACK(std::optional<double>)
INSTANTIATE_PACK(std::optional<std::string>)
INSTANTIATE_PACK(std::pair<double, double>)
INSTANTIATE_PACK(std::optional<std::pair<double,double>>)
INSTANTIATE_PACK(std::map<std::string,std::vector<int>>)
INSTANTIATE_PACK(std::map<std::string,std::map<std::pair<int,int>,int>>)
INSTANTIATE_PACK(std::map<std::string,int>)
INSTANTIATE_PACK(std::map<std::string,double>)
INSTANTIATE_PACK(std::map<int,int>)
INSTANTIATE_PACK(std::map<int,data::AquiferData>)
INSTANTIATE_PACK(std::unordered_map<std::string,size_t>)
INSTANTIATE_PACK(std::unordered_map<std::string,std::string>)
INSTANTIATE_PACK(std::unordered_set<std::string>)
INSTANTIATE_PACK(std::set<std::string>)

#undef INSTANTIATE_PACK

} // end namespace Mpi

RestartValue loadParallelRestart(const EclipseIO* eclIO, Action::State& actionState, SummaryState& summaryState,
                                 const std::vector<Opm::RestartKey>& solutionKeys,
                                 const std::vector<Opm::RestartKey>& extraKeys,
                                 Parallel::Communication comm)
{
#if HAVE_MPI
    RestartValue restartValues{};

    if (eclIO)
    {
        assert(comm.rank() == 0);
        restartValues = eclIO->loadRestart(actionState, summaryState, solutionKeys, extraKeys);
        int packedSize = Mpi::packSize(restartValues, comm);
        std::vector<char> buffer(packedSize);
        int position=0;
        Mpi::pack(restartValues, buffer, position, comm);
        comm.broadcast(&position, 1, 0);
        comm.broadcast(buffer.data(), position, 0);
        std::vector<char> buf2 = summaryState.serialize();
        int size = buf2.size();
        comm.broadcast(&size, 1, 0);
        comm.broadcast(buf2.data(), size, 0);
    }
    else
    {
        int bufferSize{};
        comm.broadcast(&bufferSize, 1, 0);
        std::vector<char> buffer(bufferSize);
        comm.broadcast(buffer.data(), bufferSize, 0);
        int position{};
        Mpi::unpack(restartValues, buffer, position, comm);
        comm.broadcast(&bufferSize, 1, 0);
        buffer.resize(bufferSize);
        comm.broadcast(buffer.data(), bufferSize, 0);
        summaryState.deserialize(buffer);
    }
    return restartValues;
#else
    (void) comm;
    return eclIO->loadRestart(actionState, summaryState, solutionKeys, extraKeys);
#endif
}

} // end namespace Opm
