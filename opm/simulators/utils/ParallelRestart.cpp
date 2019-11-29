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
#include <opm/parser/eclipse/EclipseState/Grid/NNC.hpp>
#include <opm/parser/eclipse/EclipseState/Edit/EDITNNC.hpp>
#include <opm/parser/eclipse/EclipseState/SimulationConfig/ThresholdPressure.hpp>
#include <opm/parser/eclipse/EclipseState/Tables/ColumnSchema.hpp>
#include <opm/parser/eclipse/EclipseState/Tables/Rock2dTable.hpp>
#include <opm/parser/eclipse/EclipseState/Tables/Rock2dtrTable.hpp>
#include <opm/parser/eclipse/EclipseState/Tables/TableColumn.hpp>
#include <opm/parser/eclipse/EclipseState/Tables/TableSchema.hpp>
#include <dune/common/parallel/mpitraits.hh>

#define HANDLE_AS_POD(T) \
  std::size_t packSize(const T& data, Dune::MPIHelper::MPICommunicator comm) \
  { \
      return packSize(data, comm, std::integral_constant<bool,true>()); \
  } \
  void pack(const T& data, std::vector<char>& buffer, int& position, \
            Dune::MPIHelper::MPICommunicator comm) \
  { \
      pack(data, buffer, position, comm, std::integral_constant<bool,true>()); \
  } \
  void unpack(T& data, std::vector<char>& buffer, int& position, \
              Dune::MPIHelper::MPICommunicator comm) \
  { \
      unpack(data, buffer, position, comm, std::integral_constant<bool,true>()); \
  }

namespace Opm
{
namespace Mpi
{
template<class T>
std::size_t packSize(const T*, std::size_t, Dune::MPIHelper::MPICommunicator,
                     std::integral_constant<bool, false>)
{
    OPM_THROW(std::logic_error, "Packing not (yet) supported for this non-pod type.");
}

template<class T>
std::size_t packSize(const T*, std::size_t l, Dune::MPIHelper::MPICommunicator comm,
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
std::size_t packSize(const T* data, std::size_t l, Dune::MPIHelper::MPICommunicator comm)
{
    return packSize(data, l, comm, typename std::is_pod<T>::type());
}

template<class T>
std::size_t packSize(const T&, Dune::MPIHelper::MPICommunicator,
                     std::integral_constant<bool, false>)
{
    OPM_THROW(std::logic_error, "Packing not (yet) supported for this non-pod type.");
}

template<class T>
std::size_t packSize(const T&, Dune::MPIHelper::MPICommunicator comm,
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
std::size_t packSize(const T& data, Dune::MPIHelper::MPICommunicator comm)
{
    return packSize(data, comm, typename std::is_pod<T>::type());
}

template<class T1, class T2>
std::size_t packSize(const std::pair<T1,T2>& data, Dune::MPIHelper::MPICommunicator comm)
{
    return packSize(data.first, comm) + packSize(data.second, comm);
}

template<class T, class A>
std::size_t packSize(const std::vector<T,A>& data, Dune::MPIHelper::MPICommunicator comm)
{
    if (std::is_pod<T>::value)
        // size written automatically
        return packSize(data.data(), data.size(), comm);

    std::size_t size = packSize(data.size(), comm);

    for (const auto& entry: data)
        size+=packSize(entry, comm);

    return size;
}

template<class A>
std::size_t packSize(const std::vector<bool,A>& data, Dune::MPIHelper::MPICommunicator comm)
{
    bool entry;
    return packSize(data.size(), comm) + data.size()*packSize(entry,comm);
}

template<class Key, class Value>
std::size_t packSize(const OrderedMap<Key,Value>& data, Dune::MPIHelper::MPICommunicator comm)
{
  return packSize(data.getIndex(), comm) + packSize(data.getStorage(), comm);
}


std::size_t packSize(const char* str, Dune::MPIHelper::MPICommunicator comm)
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

std::size_t packSize(const std::string& str, Dune::MPIHelper::MPICommunicator comm)
{
    return packSize(str.c_str(), comm);
}

template<class T1, class T2, class C, class A>
std::size_t packSize(const std::map<T1,T2,C,A>& data, Dune::MPIHelper::MPICommunicator comm)
{
    std::size_t totalSize = packSize(data.size(), comm);
    for (const auto& entry: data)
    {
        totalSize += packSize(entry, comm);
    }
    return totalSize;
}

template<class T1, class T2, class H, class P, class A>
std::size_t packSize(const std::unordered_map<T1,T2,H,P,A>& data, Dune::MPIHelper::MPICommunicator comm)
{
    std::size_t totalSize = packSize(data.size(), comm);
    for (const auto& entry: data)
    {
        totalSize += packSize(entry, comm);
    }
    return totalSize;
}

HANDLE_AS_POD(data::Connection)
HANDLE_AS_POD(data::Rates)
HANDLE_AS_POD(data::Segment)

std::size_t packSize(const data::Well& data, Dune::MPIHelper::MPICommunicator comm)
{
    std::size_t size = packSize(data.rates, comm);
    size += packSize(data.bhp, comm) + packSize(data.thp, comm);
    size += packSize(data.temperature, comm);
    size += packSize(data.control, comm);
    size += packSize(data.connections, comm);
    size += packSize(data.segments, comm);
    return size;
}

std::size_t packSize(const data::CellData& data, Dune::MPIHelper::MPICommunicator comm)
{
    return packSize(data.dim, comm) + packSize(data.data, comm) + packSize(data.target, comm);
}

std::size_t packSize(const RestartKey& data, Dune::MPIHelper::MPICommunicator comm)
{
    return packSize(data.key, comm) + packSize(data.dim, comm) + packSize(data.required, comm);
}

std::size_t packSize(const data::Solution& data, Dune::MPIHelper::MPICommunicator comm)
{
    // Needs explicit conversion to a supported base type holding the data
    // to prevent throwing.
    return packSize(static_cast<const std::map< std::string, data::CellData>&>(data), comm);
}

std::size_t packSize(const data::WellRates& data, Dune::MPIHelper::MPICommunicator comm)
{
    // Needs explicit conversion to a supported base type holding the data
    // to prevent throwing.
    return packSize(static_cast<const std::map< std::string, data::Well>&>(data), comm);
}

std::size_t packSize(const RestartValue& data, Dune::MPIHelper::MPICommunicator comm)
{
    return packSize(data.solution, comm) + packSize(data.wells, comm) + packSize(data.extra, comm);
}

std::size_t packSize(const ThresholdPressure& data, Dune::MPIHelper::MPICommunicator comm)
{
   return packSize(data.active(), comm) +
          packSize(data.restart(), comm) +
          packSize(data.thresholdPressureTable(), comm) +
          packSize(data.pressureTable(), comm);
}

std::size_t packSize(const NNC& data, Dune::MPIHelper::MPICommunicator comm)
{
   return packSize(data.data(), comm);
}

std::size_t packSize(const EDITNNC& data, Dune::MPIHelper::MPICommunicator comm)
{
   return packSize(data.data(), comm);
}

std::size_t packSize(const Rock2dTable& data, Dune::MPIHelper::MPICommunicator comm)
{
   return packSize(data.pvmultValues(), comm) +
          packSize(data.pressureValues(), comm);
}

std::size_t packSize(const Rock2dtrTable& data, Dune::MPIHelper::MPICommunicator comm)
{
   return packSize(data.transMultValues(), comm) +
          packSize(data.pressureValues(), comm);
}

std::size_t packSize(const ColumnSchema& data, Dune::MPIHelper::MPICommunicator comm)
{
    std::size_t res = packSize(data.name(), comm) +
                      packSize(data.order(), comm) +
                      packSize(data.getDefaultMode(), comm);
    if (data.getDefaultMode() == Table::DEFAULT_CONST) {
        res += packSize(data.getDefaultValue(), comm);
    }

    return res;
}

std::size_t packSize(const TableSchema& data, Dune::MPIHelper::MPICommunicator comm)
{
   return packSize(data.getColumns(), comm);
}

std::size_t packSize(const TableColumn& data, Dune::MPIHelper::MPICommunicator comm)
{
   return packSize(data.schema(), comm) +
          packSize(data.name(), comm) +
          packSize(data.values(), comm) +
          packSize(data.defaults(), comm) +
          packSize(data.defaultCount(), comm);
}

////// pack routines

template<class T>
void pack(const T*, std::size_t, std::vector<char>&, int&,
          Dune::MPIHelper::MPICommunicator, std::integral_constant<bool, false>)
{
    OPM_THROW(std::logic_error, "Packing not (yet) supported for this non-pod type.");
}

template<class T>
void pack(const T* data, std::size_t l, std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm,
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
          Dune::MPIHelper::MPICommunicator comm)
{
    pack(data, l, buffer, position, comm, typename std::is_pod<T>::type());
}

template<class T>
void pack(const T&, std::vector<char>&, int&,
          Dune::MPIHelper::MPICommunicator, std::integral_constant<bool, false>)
{
    OPM_THROW(std::logic_error, "Packing not (yet) supported for this non-pod type.");
}

template<class T>
void pack(const T& data, std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm, std::integral_constant<bool, true>)
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
          Dune::MPIHelper::MPICommunicator comm)
{
    pack(data, buffer, position, comm, typename std::is_pod<T>::type());
}

template<class T1, class T2>
void pack(const std::pair<T1,T2>& data, std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    pack(data.first, buffer, position, comm);
    pack(data.second, buffer, position, comm);
}

template<class T, class A>
void pack(const std::vector<T, A>& data, std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
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

template<class A>
void pack(const std::vector<bool,A>& data, std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    pack(data.size(), buffer, position, comm);
    for (const auto& entry : data) {
        bool b = entry;
        pack(b, buffer, position, comm);
    }
}

template<class Key, class Value>
void pack(const OrderedMap<Key, Value>& data, std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    pack(data.getIndex(), buffer, position, comm);
    pack(data.getStorage(), buffer, position, comm);
}

void pack(const char* str, std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
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
          Dune::MPIHelper::MPICommunicator comm)
{
    pack(str.c_str(), buffer, position, comm);
}

template<class T1, class T2, class C, class A>
void pack(const std::map<T1,T2,C,A>& data, std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    pack(data.size(), buffer, position, comm);

    for (const auto& entry: data)
    {
        pack(entry, buffer, position, comm);
    }
}

template<class T1, class T2, class H, class P, class A>
void pack(const std::unordered_map<T1,T2,H,P,A>& data, std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    pack(data.size(), buffer, position, comm);

    for (const auto& entry: data)
    {
        pack(entry, buffer, position, comm);
    }
}

void pack(const data::Well& data, std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    pack(data.rates, buffer, position, comm);
    pack(data.bhp, buffer, position, comm);
    pack(data.thp, buffer, position, comm);
    pack(data.temperature, buffer, position, comm);
    pack(data.control, buffer, position, comm);
    pack(data.connections, buffer, position, comm);
    pack(data.segments, buffer, position, comm);
}

void pack(const RestartKey& data, std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    pack(data.key, buffer, position, comm);
    pack(data.dim, buffer, position, comm);
    pack(data.required, buffer, position, comm);
}

void pack(const data::CellData& data, std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    pack(data.dim, buffer, position, comm);
    pack(data.data, buffer, position, comm);
    pack(data.target, buffer, position, comm);
}

void pack(const data::Solution& data, std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    // Needs explicit conversion to a supported base type holding the data
    // to prevent throwing.
    pack(static_cast<const std::map< std::string, data::CellData>&>(data),
         buffer, position, comm);
}

void pack(const data::WellRates& data, std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    // Needs explicit conversion to a supported base type holding the data
    // to prevent throwing.
    pack(static_cast<const std::map< std::string, data::Well>&>(data),
         buffer, position, comm);
}

void pack(const RestartValue& data, std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    pack(data.solution, buffer, position, comm);
    pack(data.wells, buffer, position, comm);
    pack(data.extra, buffer, position, comm);
}

void pack(const ThresholdPressure& data, std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    pack(data.active(), buffer, position, comm);
    pack(data.restart(), buffer, position, comm);
    pack(data.thresholdPressureTable(), buffer, position, comm);
    pack(data.pressureTable(), buffer, position, comm);
}

void pack(const NNC& data, std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    pack(data.data(), buffer, position, comm);
}

void pack(const EDITNNC& data, std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    pack(data.data(), buffer, position, comm);
}

void pack(const Rock2dTable& data, std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    pack(data.pvmultValues(), buffer, position, comm);
    pack(data.pressureValues(), buffer, position, comm);
}

void pack(const Rock2dtrTable& data, std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    pack(data.transMultValues(), buffer, position, comm);
    pack(data.pressureValues(), buffer, position, comm);
}

void pack(const ColumnSchema& data, std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    pack(data.name(), buffer, position, comm);
    pack(data.order(), buffer, position, comm);
    pack(data.getDefaultMode(), buffer, position, comm);
    if (data.getDefaultMode() == Table::DEFAULT_CONST)
        pack(data.getDefaultValue(), buffer, position, comm);
}

void pack(const TableSchema& data, std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    pack(data.getColumns(), buffer, position, comm);
}

void pack(const TableColumn& data, std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    pack(data.schema(), buffer, position, comm);
    pack(data.name(), buffer, position, comm);
    pack(data.values(), buffer, position, comm);
    pack(data.defaults(), buffer, position, comm);
    pack(data.defaultCount(), buffer, position, comm);
}

/// unpack routines

template<class T>
void unpack(T*, const std::size_t&, std::vector<char>&, int&,
            Dune::MPIHelper::MPICommunicator, std::integral_constant<bool, false>)
{
    OPM_THROW(std::logic_error, "Packing not (yet) supported for this non-pod type.");
}

template<class T>
void unpack(T* data, const std::size_t& l, std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm,
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
            Dune::MPIHelper::MPICommunicator comm)
{
    unpack(data, l, buffer, position, comm, typename std::is_pod<T>::type());
}

template<class T>
void unpack(T&, std::vector<char>&, int&,
            Dune::MPIHelper::MPICommunicator, std::integral_constant<bool, false>)
{
    OPM_THROW(std::logic_error, "Packing not (yet) supported for this non-pod type.");
}

template<class T>
void unpack(T& data, std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm, std::integral_constant<bool, true>)
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
            Dune::MPIHelper::MPICommunicator comm)
{
    unpack(data, buffer, position, comm, typename std::is_pod<T>::type());
}

template<class T1, class T2>
void unpack(std::pair<T1,T2>& data, std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm)
{
    unpack(data.first, buffer, position, comm);
    unpack(data.second, buffer, position, comm);
}

template<class T, class A>
void unpack(std::vector<T,A>& data, std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm)
{
    std::size_t length=0;
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
            Dune::MPIHelper::MPICommunicator comm)
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

template<class Key, class Value>
void unpack(OrderedMap<Key,Value>& data, std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm)
{
  typename OrderedMap<Key,Value>::index_type index;
  typename OrderedMap<Key,Value>::storage_type storage;
  unpack(index, buffer, position, comm);
  unpack(storage, buffer, position, comm);
  data = OrderedMap<Key,Value>(index, storage);
}

void unpack(char* str, std::size_t length, std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm)
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
            Dune::MPIHelper::MPICommunicator comm)
{
    std::size_t length=0;
    unpack(length, buffer, position, comm);
    std::vector<char> cStr(length, '\0');
    unpack(cStr.data(), length, buffer, position, comm);
    assert(str.empty());
    str.append(cStr.data());
}

template<class T1, class T2, class C, class A>
void unpack(std::map<T1,T2,C,A>& data, std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm)
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
            Dune::MPIHelper::MPICommunicator comm)
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
            Dune::MPIHelper::MPICommunicator comm)
{
    unpack(data.rates, buffer, position, comm);
    unpack(data.bhp, buffer, position, comm);
    unpack(data.thp, buffer, position, comm);
    unpack(data.temperature, buffer, position, comm);
    unpack(data.control, buffer, position, comm);
    unpack(data.connections, buffer, position, comm);
    unpack(data.segments, buffer, position, comm);
}

void unpack(RestartKey& data, std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm)
{
    unpack(data.key, buffer, position, comm);
    unpack(data.dim, buffer, position, comm);
    unpack(data.required, buffer, position, comm);
}

void unpack(data::CellData& data, std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm)
{
    unpack(data.dim, buffer, position, comm);
    unpack(data.data, buffer, position, comm);
    unpack(data.target, buffer, position, comm);
}

void unpack(data::Solution& data, std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm)
{
    // Needs explicit conversion to a supported base type holding the data
    // to prevent throwing.
    unpack(static_cast<std::map< std::string, data::CellData>&>(data),
           buffer, position, comm);
}

void unpack(data::WellRates& data, std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm)
{
    // Needs explicit conversion to a supported base type holding the data
    // to prevent throwing.
    unpack(static_cast<std::map< std::string, data::Well>&>(data),
           buffer, position, comm);
}

void unpack(RestartValue& data, std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm)
{
    unpack(data.solution, buffer, position, comm);
    unpack(data.wells, buffer, position, comm);
    unpack(data.extra, buffer, position, comm);
}

void unpack(ThresholdPressure& data, std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm)
{
    ThresholdPressure::ThresholdPressureTable thpTable;
    ThresholdPressure::PressureTable pTable;
    bool active, restart;
    unpack(active, buffer, position, comm);
    unpack(restart, buffer, position, comm);
    unpack(thpTable, buffer, position, comm);
    unpack(pTable, buffer, position, comm);

    data = ThresholdPressure(active, restart, thpTable, pTable);
}

void unpack(NNC& data, std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm)
{
    std::vector<NNCdata> res;
    unpack(res, buffer, position, comm);
    data = NNC(res);
}

void unpack(EDITNNC& data, std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm)
{
    std::vector<NNCdata> res;
    unpack(res, buffer, position, comm);
    data = EDITNNC(res);
}

void unpack(Rock2dTable& data, std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm)
{
    std::vector<std::vector<double>> pvmultValues;
    std::vector<double> pressureValues;
    unpack(pvmultValues, buffer, position, comm);
    unpack(pressureValues, buffer, position, comm);
    data = Rock2dTable(pvmultValues, pressureValues);
}

void unpack(Rock2dtrTable& data, std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm)
{
    std::vector<std::vector<double>> transMultValues;
    std::vector<double> pressureValues;
    unpack(transMultValues, buffer, position, comm);
    unpack(pressureValues, buffer, position, comm);
    data = Rock2dtrTable(transMultValues, pressureValues);
}

void unpack(ColumnSchema& data, std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm)
{
    std::string name;
    Table::ColumnOrderEnum order;
    Table::DefaultAction action;
    unpack(name, buffer, position, comm);
    unpack(order, buffer, position, comm);
    unpack(action, buffer, position, comm);
    if (action == Table::DEFAULT_CONST) {
        double value;
        unpack(value, buffer, position, comm);
        data = ColumnSchema(name, order, value);
    } else
        data = ColumnSchema(name, order, action);
}

void unpack(TableSchema& data, std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm)
{
    OrderedMap<std::string, ColumnSchema> columns;
    unpack(columns, buffer, position, comm);
    data = TableSchema(columns);
}

void unpack(TableColumn& data, std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm)
{
    ColumnSchema schema;
    std::string name;
    std::vector<double> values;
    std::vector<bool> defaults;
    size_t defaultCount;
    unpack(schema, buffer, position, comm);
    unpack(name, buffer, position, comm);
    unpack(values, buffer, position, comm);
    unpack(defaults, buffer, position, comm);
    unpack(defaultCount, buffer, position, comm);
    data = TableColumn(schema, name, values, defaults, defaultCount);
}

} // end namespace Mpi
RestartValue loadParallelRestart(const EclipseIO* eclIO, SummaryState& summaryState,
                                 const std::vector<Opm::RestartKey>& solutionKeys,
                                 const std::vector<Opm::RestartKey>& extraKeys,
                                 Dune::CollectiveCommunication<Dune::MPIHelper::MPICommunicator> comm)
{
#if HAVE_MPI
    data::Solution sol;
    data::Wells wells;
    RestartValue restartValues(sol, wells);

    if (eclIO)
    {
        assert(comm.rank() == 0);
        restartValues = eclIO->loadRestart(summaryState, solutionKeys, extraKeys);
        int packedSize = Mpi::packSize(restartValues, comm);
        std::vector<char> buffer(packedSize);
        int position=0;
        Mpi::pack(restartValues, buffer, position, comm);
        comm.broadcast(&position, 1, 0);
        comm.broadcast(buffer.data(), position, 0);
    }
    else
    {
        int bufferSize{};
        comm.broadcast(&bufferSize, 1, 0);
        std::vector<char> buffer(bufferSize);
        comm.broadcast(buffer.data(), bufferSize, 0);
        int position{};
        Mpi::unpack(restartValues, buffer, position, comm);
    }
    return restartValues;
#else
    (void) comm;
    return eclIO->loadRestart(summaryState, solutionKeys, extraKeys);
#endif
}
} // end namespace Opm
