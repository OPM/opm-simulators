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
#include <ebos/eclmpiserializer.hh>
#include <opm/output/data/Aquifer.hpp>
#include <opm/output/data/Cells.hpp>
#include <opm/output/data/Groups.hpp>
#include <opm/output/data/GuideRateValue.hpp>
#include <opm/output/data/Solution.hpp>
#include <opm/output/data/Wells.hpp>
#include <opm/input/eclipse/EclipseState/Util/OrderedMap.hpp>
#include <opm/output/eclipse/EclipseIO.hpp>
#include <opm/input/eclipse/EclipseState/Util/OrderedMap.hpp>
#include <opm/output/eclipse/RestartValue.hpp>
#include <opm/input/eclipse/Schedule/SummaryState.hpp>

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
    bool entry = false;
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
HANDLE_AS_POD(data::CurrentControl)
HANDLE_AS_POD(data::FetkovichData)
HANDLE_AS_POD(data::GroupConstraints)
HANDLE_AS_POD(data::NodeData)

template <class T>
struct Packing
{
};

template <>
struct Packing<data::Rates>
{
    static std::size_t packSize(const data::Rates& data, Opm::Parallel::MPIComm comm)
    {
        std::size_t totalSize = 0;
        totalSize += Mpi::packSize(data.mask, comm);
        totalSize += Mpi::packSize(data.wat, comm);
        totalSize += Mpi::packSize(data.oil, comm);
        totalSize += Mpi::packSize(data.gas, comm);
        totalSize += Mpi::packSize(data.polymer, comm);
        totalSize += Mpi::packSize(data.solvent, comm);
        totalSize += Mpi::packSize(data.energy, comm);
        totalSize += Mpi::packSize(data.dissolved_gas, comm);
        totalSize += Mpi::packSize(data.vaporized_oil, comm);
        totalSize += Mpi::packSize(data.reservoir_water, comm);
        totalSize += Mpi::packSize(data.reservoir_oil, comm);
        totalSize += Mpi::packSize(data.reservoir_gas, comm);
        totalSize += Mpi::packSize(data.productivity_index_water, comm);
        totalSize += Mpi::packSize(data.productivity_index_oil, comm);
        totalSize += Mpi::packSize(data.productivity_index_gas, comm);
        totalSize += Mpi::packSize(data.well_potential_water, comm);
        totalSize += Mpi::packSize(data.well_potential_oil, comm);
        totalSize += Mpi::packSize(data.well_potential_gas, comm);
        totalSize += Mpi::packSize(data.brine, comm);
        totalSize += Mpi::packSize(data.alq, comm);
        totalSize += Mpi::packSize(data.tracer, comm);
        totalSize += Mpi::packSize(data.micp, comm);
        totalSize += Mpi::packSize(data.vaporized_water, comm);
        return totalSize;
    }

    static void pack(const data::Rates& data, std::vector<char>& buffer, int& position, Opm::Parallel::MPIComm comm)
    {
        Mpi::pack(data.mask, buffer, position, comm);
        Mpi::pack(data.wat, buffer, position, comm);
        Mpi::pack(data.oil, buffer, position, comm);
        Mpi::pack(data.gas, buffer, position, comm);
        Mpi::pack(data.polymer, buffer, position, comm);
        Mpi::pack(data.solvent, buffer, position, comm);
        Mpi::pack(data.energy, buffer, position, comm);
        Mpi::pack(data.dissolved_gas, buffer, position, comm);
        Mpi::pack(data.vaporized_oil, buffer, position, comm);
        Mpi::pack(data.reservoir_water, buffer, position, comm);
        Mpi::pack(data.reservoir_oil, buffer, position, comm);
        Mpi::pack(data.reservoir_gas, buffer, position, comm);
        Mpi::pack(data.productivity_index_water, buffer, position, comm);
        Mpi::pack(data.productivity_index_oil, buffer, position, comm);
        Mpi::pack(data.productivity_index_gas, buffer, position, comm);
        Mpi::pack(data.well_potential_water, buffer, position, comm);
        Mpi::pack(data.well_potential_oil, buffer, position, comm);
        Mpi::pack(data.well_potential_gas, buffer, position, comm);
        Mpi::pack(data.brine, buffer, position, comm);
        Mpi::pack(data.alq, buffer, position, comm);
        Mpi::pack(data.tracer, buffer, position, comm);
        Mpi::pack(data.micp, buffer, position, comm);
        Mpi::pack(data.vaporized_water, buffer, position, comm);
    }

    static void unpack(data::Rates& data, std::vector<char>& buffer, int& position, Opm::Parallel::MPIComm comm)
    {
        Mpi::unpack(data.mask, buffer, position, comm);
        Mpi::unpack(data.wat, buffer, position, comm);
        Mpi::unpack(data.oil, buffer, position, comm);
        Mpi::unpack(data.gas, buffer, position, comm);
        Mpi::unpack(data.polymer, buffer, position, comm);
        Mpi::unpack(data.solvent, buffer, position, comm);
        Mpi::unpack(data.energy, buffer, position, comm);
        Mpi::unpack(data.dissolved_gas, buffer, position, comm);
        Mpi::unpack(data.vaporized_oil, buffer, position, comm);
        Mpi::unpack(data.reservoir_water, buffer, position, comm);
        Mpi::unpack(data.reservoir_oil, buffer, position, comm);
        Mpi::unpack(data.reservoir_gas, buffer, position, comm);
        Mpi::unpack(data.productivity_index_water, buffer, position, comm);
        Mpi::unpack(data.productivity_index_oil, buffer, position, comm);
        Mpi::unpack(data.productivity_index_gas, buffer, position, comm);
        Mpi::unpack(data.well_potential_water, buffer, position, comm);
        Mpi::unpack(data.well_potential_oil, buffer, position, comm);
        Mpi::unpack(data.well_potential_gas, buffer, position, comm);
        Mpi::unpack(data.brine, buffer, position, comm);
        Mpi::unpack(data.alq, buffer, position, comm);
        Mpi::unpack(data.tracer, buffer, position, comm);
        Mpi::unpack(data.micp, buffer, position, comm);
        Mpi::unpack(data.vaporized_water, buffer, position, comm);
    }
};

template <>
struct Packing<data::Connection>
{
    static std::size_t packSize(const data::Connection& data, Opm::Parallel::MPIComm comm)
    {
        std::size_t totalSize = 0;
        totalSize += Mpi::packSize(data.index, comm);
        totalSize += Mpi::packSize(data.rates, comm);
        totalSize += Mpi::packSize(data.pressure, comm);
        totalSize += Mpi::packSize(data.reservoir_rate, comm);
        totalSize += Mpi::packSize(data.cell_pressure, comm);
        totalSize += Mpi::packSize(data.cell_saturation_water, comm);
        totalSize += Mpi::packSize(data.cell_saturation_gas, comm);
        totalSize += Mpi::packSize(data.effective_Kh, comm);
        totalSize += Mpi::packSize(data.trans_factor, comm);
        return totalSize;
    }

    static void pack(const data::Connection& data, std::vector<char>& buffer, int& position, Opm::Parallel::MPIComm comm)
    {
        Mpi::pack(data.index, buffer, position, comm);
        Mpi::pack(data.rates, buffer, position, comm);
        Mpi::pack(data.pressure, buffer, position, comm);
        Mpi::pack(data.reservoir_rate, buffer, position, comm);
        Mpi::pack(data.cell_pressure, buffer, position, comm);
        Mpi::pack(data.cell_saturation_water, buffer, position, comm);
        Mpi::pack(data.cell_saturation_gas, buffer, position, comm);
        Mpi::pack(data.effective_Kh, buffer, position, comm);
        Mpi::pack(data.trans_factor, buffer, position, comm);
    }

    static void unpack(data::Connection& data, std::vector<char>& buffer, int& position, Opm::Parallel::MPIComm comm)
    {
        Mpi::unpack(data.index, buffer, position, comm);
        Mpi::unpack(data.rates, buffer, position, comm);
        Mpi::unpack(data.pressure, buffer, position, comm);
        Mpi::unpack(data.reservoir_rate, buffer, position, comm);
        Mpi::unpack(data.cell_pressure, buffer, position, comm);
        Mpi::unpack(data.cell_saturation_water, buffer, position, comm);
        Mpi::unpack(data.cell_saturation_gas, buffer, position, comm);
        Mpi::unpack(data.effective_Kh, buffer, position, comm);
        Mpi::unpack(data.trans_factor, buffer, position, comm);
    }
};

template <>
struct Packing<data::SegmentPressures>
{
    static std::size_t packSize(const data::SegmentPressures& data, Opm::Parallel::MPIComm comm)
    {
        return Mpi::packSize(data.values_, comm);
    }

    static void pack(const data::SegmentPressures& data, std::vector<char>& buffer, int& position, Opm::Parallel::MPIComm comm)
    {
        Mpi::pack(data.values_, buffer, position, comm);
    }

    static void unpack(data::SegmentPressures& data, std::vector<char>& buffer, int& position, Opm::Parallel::MPIComm comm)
    {
        Mpi::unpack(data.values_, buffer, position, comm);
    }
};

template <>
struct Packing<data::Segment>
{
    static std::size_t packSize(const data::Segment& data, Opm::Parallel::MPIComm comm)
    {
        std::size_t totalSize = 0;
        totalSize += Mpi::packSize(data.rates, comm);
        totalSize += Mpi::packSize(data.pressures, comm);
        totalSize += Mpi::packSize(data.segNumber, comm);
        return totalSize;
    }

    static void pack(const data::Segment& data, std::vector<char>& buffer, int& position, Opm::Parallel::MPIComm comm)
    {
        Mpi::pack(data.rates, buffer, position, comm);
        Mpi::pack(data.pressures, buffer, position, comm);
        Mpi::pack(data.segNumber, buffer, position, comm);
    }

    static void unpack(data::Segment& data, std::vector<char>& buffer, int& position, Opm::Parallel::MPIComm comm)
    {
        Mpi::unpack(data.rates, buffer, position, comm);
        Mpi::unpack(data.pressures, buffer, position, comm);
        Mpi::unpack(data.segNumber, buffer, position, comm);
    }
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

std::size_t packSize(const data::Rates& data, Opm::Parallel::MPIComm comm)
{
    return Packing<data::Rates>::packSize(data, comm);
}

std::size_t packSize(const data::Connection& data, Opm::Parallel::MPIComm comm)
{
    return Packing<data::Connection>::packSize(data, comm);
}

std::size_t packSize(const data::SegmentPressures& data, Opm::Parallel::MPIComm comm)
{
    return Packing<data::SegmentPressures>::packSize(data, comm);
}

std::size_t packSize(const data::Segment& data, Opm::Parallel::MPIComm comm)
{
    return Packing<data::Segment>::packSize(data, comm);
}

template<std::size_t Size>
std::size_t packSize(const std::bitset<Size>& data, Opm::Parallel::MPIComm comm)
{
    return Packing<std::bitset<Size>>::packSize(data, comm);
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

void pack(const data::Rates& data, std::vector<char>& buffer, int& position,
          Opm::Parallel::MPIComm comm)
{
    Packing<data::Rates>::pack(data, buffer, position, comm);
}

void pack(const data::Connection& data, std::vector<char>& buffer, int& position,
          Opm::Parallel::MPIComm comm)
{
    Packing<data::Connection>::pack(data, buffer, position, comm);
}

void pack(const data::SegmentPressures& data, std::vector<char>& buffer, int& position,
          Opm::Parallel::MPIComm comm)
{
    Packing<data::SegmentPressures>::pack(data, buffer, position, comm);
}

void pack(const data::Segment& data, std::vector<char>& buffer, int& position,
          Opm::Parallel::MPIComm comm)
{
    Packing<data::Segment>::pack(data, buffer, position, comm);
}

template<std::size_t Size>
void pack(const std::bitset<Size>& data, std::vector<char>& buffer,
          int& position, Opm::Parallel::MPIComm comm)
{
    Packing<std::bitset<Size>>::pack(data, buffer, position, comm);
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

void unpack(data::Rates& data, std::vector<char>& buffer, int& position,
            Opm::Parallel::MPIComm comm)
{
    Packing<data::Rates>::unpack(data, buffer, position, comm);
}

void unpack(data::Connection& data, std::vector<char>& buffer, int& position,
            Opm::Parallel::MPIComm comm)
{
    Packing<data::Connection>::unpack(data, buffer, position, comm);
}

void unpack(data::SegmentPressures& data, std::vector<char>& buffer, int& position,
            Opm::Parallel::MPIComm comm)
{
    Packing<data::SegmentPressures>::unpack(data, buffer, position, comm);
}

void unpack(data::Segment& data, std::vector<char>& buffer, int& position,
            Opm::Parallel::MPIComm comm)
{
    Packing<data::Segment>::unpack(data, buffer, position, comm);
}

template<std::size_t Size>
void unpack(std::bitset<Size>& data, std::vector<char>& buffer, int& position,
            Opm::Parallel::MPIComm comm)
{
    Packing<std::bitset<Size>>::unpack(data, buffer, position, comm);
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
INSTANTIATE_PACK(std::array<double,4>)
INSTANTIATE_PACK(std::array<double,5>)
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
INSTANTIATE_PACK(std::unordered_map<std::string,size_t,Opm::OrderedMapDetail::TruncatedStringHash<std::string::npos>, Opm::OrderedMapDetail::TruncatedStringEquals<std::string::npos>>)
INSTANTIATE_PACK(std::unordered_map<std::string,size_t,Opm::OrderedMapDetail::TruncatedStringHash<8>, Opm::OrderedMapDetail::TruncatedStringEquals<8>>)
INSTANTIATE_PACK(std::unordered_map<std::string,std::string>)
INSTANTIATE_PACK(std::unordered_set<std::string>)
INSTANTIATE_PACK(std::set<std::string>)
INSTANTIATE_PACK(std::bitset<4>)


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
    }

    EclMpiSerializer ser(comm);
    ser.broadcast(restartValues);
    ser.broadcast(summaryState);
    return restartValues;
#else
    (void) comm;
    return eclIO->loadRestart(actionState, summaryState, solutionKeys, extraKeys);
#endif
}

} // end namespace Opm
