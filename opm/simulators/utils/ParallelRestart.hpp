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

#include <opm/output/eclipse/RestartValue.hpp>
#include <opm/output/eclipse/EclipseIO.hpp>
#include <opm/output/eclipse/Summary.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/DynamicState.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/DynamicVector.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/Group/GConSale.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/Group/GConSump.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/Group/Group.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/Group/GuideRateConfig.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/TimeMap.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/Well/Well.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/Well/WellTestConfig.hpp>
#include <opm/parser/eclipse/EclipseState/Util/OrderedMap.hpp>

#include <dune/common/parallel/mpihelper.hh>

#include <set>
#include <tuple>
#include <vector>
#include <map>
#include <unordered_map>

namespace Opm
{

namespace Action {
    class Actions;
    class ActionX;
    class AST;
    class ASTNode;
    class Condition;
    class Quantity;
}

class Aqudims;
class Connection;
class DeckItem;
class DeckRecord;
class Dimension;
class Events;
template<class T> class IOrderSet;
class Location;
class MessageLimits;
class MLimits;
class OilVaporizationProperties;
class PlyvmhRecord;
class PlyvmhTable;
class Regdims;
class RestartConfig;
class RestartSchedule;
class RFTConfig;
class RocktabTable;
class Rock2dTable;
class Rock2dtrTable;
class Segment;
class SimpleTable;
class SkprpolyTable;
class SolventDensityTable;
class SpiralICD;
class StandardCond;
class Tabdims;
class TableColumn;
class TableContainer;
class TableManager;
class TableSchema;
class TimeStampUTC;
class Tuning;
class UDAValue;
class UnitSystem;
class Valve;
class VFPInjTable;
class VFPProdTable;
class WellConnections;
class WellEconProductionLimits;
class WellFoamProperties;
class WellPolymerProperties;
class WellSegments;
class WellTracerProperties;
class WList;
class WListManager;

namespace Mpi
{
template<class T>
std::size_t packSize(const T*, std::size_t, Dune::MPIHelper::MPICommunicator,
                     std::integral_constant<bool, false>);

template<class T>
std::size_t packSize(const T*, std::size_t l, Dune::MPIHelper::MPICommunicator comm,
                     std::integral_constant<bool, true>);

template<class T>
std::size_t packSize(const T* data, std::size_t l, Dune::MPIHelper::MPICommunicator comm);

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
std::size_t packSize(const std::pair<T1,T2>& data, Dune::MPIHelper::MPICommunicator comm);

template<class T, class A>
std::size_t packSize(const std::vector<T,A>& data, Dune::MPIHelper::MPICommunicator comm);

template<class K, class C, class A>
std::size_t packSize(const std::set<K,C,A>& data,
                     Dune::MPIHelper::MPICommunicator comm);

template<class T, class H, class KE, class A>
std::size_t packSize(const std::unordered_set<T,H,KE,A>& data,
                     Dune::MPIHelper::MPICommunicator comm);

template<class A>
std::size_t packSize(const std::vector<bool,A>& data, Dune::MPIHelper::MPICommunicator comm);

template<class... Ts>
std::size_t packSize(const std::tuple<Ts...>& data, Dune::MPIHelper::MPICommunicator comm);

template<class T>
std::size_t packSize(const std::shared_ptr<T>& data,
                     Dune::MPIHelper::MPICommunicator comm);

template<class T, std::size_t N>
std::size_t packSize(const std::array<T,N>& data, Dune::MPIHelper::MPICommunicator comm);

template<class T>
std::size_t packSize(const std::unique_ptr<T>& data,
                     Dune::MPIHelper::MPICommunicator comm);

std::size_t packSize(const char* str, Dune::MPIHelper::MPICommunicator comm);

std::size_t packSize(const std::string& str, Dune::MPIHelper::MPICommunicator comm);

template<class T1, class T2, class C, class A>
std::size_t packSize(const std::map<T1,T2,C,A>& data, Dune::MPIHelper::MPICommunicator comm);

template<class T1, class T2, class H, class P, class A>
std::size_t packSize(const std::unordered_map<T1,T2,H,P,A>& data, Dune::MPIHelper::MPICommunicator comm);

template<class Key, class Value>
std::size_t packSize(const OrderedMap<Key,Value>& data, Dune::MPIHelper::MPICommunicator comm);

template<class T>
std::size_t packSize(const DynamicVector<T>& data, Dune::MPIHelper::MPICommunicator comm);

template<class T>
std::size_t packSize(const DynamicState<T>& data, Dune::MPIHelper::MPICommunicator comm);

template<class T>
std::size_t packSize(const IOrderSet<T>& data, Dune::MPIHelper::MPICommunicator comm);

////// pack routines

template<class T>
void pack(const T*, std::size_t, std::vector<char>&, int&,
          Dune::MPIHelper::MPICommunicator, std::integral_constant<bool, false>);

template<class T>
void pack(const T* data, std::size_t l, std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm, std::integral_constant<bool, true>);

template<class T>
void pack(const T* data, std::size_t l, std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm);

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
          Dune::MPIHelper::MPICommunicator comm);

template<class T, class A>
void pack(const std::vector<T,A>& data, std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm);

template<class A>
void pack(const std::vector<bool,A>& data, std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm);

template<class... Ts>
void pack(const std::tuple<Ts...>& data, std::vector<char>& buffer,
          int& position, Dune::MPIHelper::MPICommunicator comm);

template<class K, class C, class A>
void pack(const std::set<K,C,A>& data,
          std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm);

template<class T, class H, class KE, class A>
void pack(const std::unordered_set<T,H,KE,A>& data,
          std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm);

template<class T>
void pack(const std::shared_ptr<T>& data, std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm);

template<class T, size_t N>
void pack(const std::array<T,N>& data, std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm);

template<class T>
void pack(const std::unique_ptr<T>& data, std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm);

template<class T1, class T2, class C, class A>
void pack(const std::map<T1,T2,C,A>& data, std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm);

template<class T1, class T2, class H, class P, class A>
void pack(const std::unordered_map<T1,T2,H,P,A>& data, std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm);

template<class Key, class Value>
void pack(const OrderedMap<Key,Value>& data, std::vector<char>& buffer,
          int& position, Dune::MPIHelper::MPICommunicator comm);

template<class T>
void pack(const DynamicState<T>& data, std::vector<char>& buffer,
          int& position, Dune::MPIHelper::MPICommunicator comm);

template<class T>
void pack(const DynamicVector<T>& data, std::vector<char>& buffer,
          int& position, Dune::MPIHelper::MPICommunicator comm);

template<class T>
void pack(const IOrderSet<T>& data, std::vector<char>& buffer,
          int& position, Dune::MPIHelper::MPICommunicator comm);

void pack(const char* str, std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm);

/// unpack routines

template<class T>
void unpack(T*, const std::size_t&, std::vector<char>&, int&,
            Dune::MPIHelper::MPICommunicator, std::integral_constant<bool, false>);

template<class T>
void unpack(T* data, const std::size_t& l, std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm,
            std::integral_constant<bool, true>);

template<class T>
void unpack(T* data, const std::size_t& l, std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm);

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
            Dune::MPIHelper::MPICommunicator comm);

template<class T, class A>
void unpack(std::vector<T,A>& data, std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm);

template<class A>
void unpack(std::vector<bool,A>& data, std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm);

template<class... Ts>
void unpack(std::tuple<Ts...>& data, std::vector<char>& buffer,
            int& position, Dune::MPIHelper::MPICommunicator comm);

template<class K, class C, class A>
void unpack(std::set<K,C,A>& data,
            std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm);

template<class T, class H, class KE, class A>
void unpack(std::unordered_set<T,H,KE,A>& data,
            std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm);

template<class T>
void unpack(std::shared_ptr<T>& data, std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm);

template<class T, size_t N>
void unpack(std::array<T,N>& data, std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm);

template<class T>
void unpack(std::unique_ptr<T>& data, std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm);

template<class T1, class T2, class C, class A>
void unpack(std::map<T1,T2,C,A>& data, std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm);

template<class T1, class T2, class H, class P, class A>
void unpack(std::unordered_map<T1,T2,H,P,A>& data, std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm);

template<class Key, class Value>
void unpack(OrderedMap<Key,Value>& data, std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm);

template<class T>
void unpack(DynamicState<T>& data, std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm);

template<class T>
void unpack(DynamicVector<T>& data, std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm);

template<class T>
void unpack(IOrderSet<T>& data, std::vector<char>& buffer,
            int& position, Dune::MPIHelper::MPICommunicator comm);

void unpack(char* str, std::size_t length, std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm);

/// prototypes for complex types

#define ADD_PACK_PROTOTYPES(T) \
  std::size_t packSize(const T& data, Dune::MPIHelper::MPICommunicator comm); \
  void pack(const T& data, std::vector<char>& buffer, int& position, \
          Dune::MPIHelper::MPICommunicator comm); \
  void unpack(T& data, std::vector<char>& buffer, int& position, \
              Dune::MPIHelper::MPICommunicator comm);

ADD_PACK_PROTOTYPES(Action::Actions)
ADD_PACK_PROTOTYPES(Action::ActionX)
ADD_PACK_PROTOTYPES(Action::AST)
ADD_PACK_PROTOTYPES(Action::ASTNode)
ADD_PACK_PROTOTYPES(Action::Condition)
ADD_PACK_PROTOTYPES(Action::Quantity)
ADD_PACK_PROTOTYPES(Aqudims)
ADD_PACK_PROTOTYPES(Connection)
ADD_PACK_PROTOTYPES(data::CellData)
ADD_PACK_PROTOTYPES(data::Connection)
ADD_PACK_PROTOTYPES(data::CurrentControl)
ADD_PACK_PROTOTYPES(data::Rates)
ADD_PACK_PROTOTYPES(data::Segment)
ADD_PACK_PROTOTYPES(data::Solution)
ADD_PACK_PROTOTYPES(data::Well)
ADD_PACK_PROTOTYPES(data::WellRates)
ADD_PACK_PROTOTYPES(Deck)
ADD_PACK_PROTOTYPES(DeckItem)
ADD_PACK_PROTOTYPES(DeckKeyword)
ADD_PACK_PROTOTYPES(DeckRecord)
ADD_PACK_PROTOTYPES(Dimension)
ADD_PACK_PROTOTYPES(Events)
ADD_PACK_PROTOTYPES(GConSale)
ADD_PACK_PROTOTYPES(GConSale::GCONSALEGroup)
ADD_PACK_PROTOTYPES(GConSump)
ADD_PACK_PROTOTYPES(GConSump::GCONSUMPGroup)
ADD_PACK_PROTOTYPES(GuideRateConfig)
ADD_PACK_PROTOTYPES(GuideRateConfig::GroupTarget)
ADD_PACK_PROTOTYPES(GuideRateConfig::WellTarget)
ADD_PACK_PROTOTYPES(GuideRateModel)
ADD_PACK_PROTOTYPES(Group)
ADD_PACK_PROTOTYPES(Group::GroupInjectionProperties)
ADD_PACK_PROTOTYPES(Group::GroupProductionProperties)
ADD_PACK_PROTOTYPES(Location)
ADD_PACK_PROTOTYPES(MessageLimits)
ADD_PACK_PROTOTYPES(MLimits)
ADD_PACK_PROTOTYPES(OilVaporizationProperties)
ADD_PACK_PROTOTYPES(PlyvmhRecord)
ADD_PACK_PROTOTYPES(PlyvmhTable)
ADD_PACK_PROTOTYPES(Regdims)
ADD_PACK_PROTOTYPES(RestartConfig)
ADD_PACK_PROTOTYPES(RestartKey)
ADD_PACK_PROTOTYPES(RestartSchedule)
ADD_PACK_PROTOTYPES(RestartValue)
ADD_PACK_PROTOTYPES(RFTConfig)
ADD_PACK_PROTOTYPES(Rock2dTable)
ADD_PACK_PROTOTYPES(Rock2dtrTable)
ADD_PACK_PROTOTYPES(RocktabTable)
ADD_PACK_PROTOTYPES(Segment)
ADD_PACK_PROTOTYPES(SimpleTable)
ADD_PACK_PROTOTYPES(SkprpolyTable)
ADD_PACK_PROTOTYPES(SolventDensityTable)
ADD_PACK_PROTOTYPES(SpiralICD)
ADD_PACK_PROTOTYPES(std::string)
ADD_PACK_PROTOTYPES(Tabdims)
ADD_PACK_PROTOTYPES(TableColumn)
ADD_PACK_PROTOTYPES(TableContainer)
ADD_PACK_PROTOTYPES(TableManager)
ADD_PACK_PROTOTYPES(TableSchema)
ADD_PACK_PROTOTYPES(TimeMap)
ADD_PACK_PROTOTYPES(TimeStampUTC)
ADD_PACK_PROTOTYPES(Tuning)
ADD_PACK_PROTOTYPES(UDAValue)
ADD_PACK_PROTOTYPES(UnitSystem)
ADD_PACK_PROTOTYPES(Valve)
ADD_PACK_PROTOTYPES(VFPInjTable)
ADD_PACK_PROTOTYPES(VFPProdTable)
ADD_PACK_PROTOTYPES(Well)
ADD_PACK_PROTOTYPES(WellType)
ADD_PACK_PROTOTYPES(Well::WellGuideRate)
ADD_PACK_PROTOTYPES(Well::WellInjectionProperties)
ADD_PACK_PROTOTYPES(Well::WellProductionProperties)
ADD_PACK_PROTOTYPES(WellBrineProperties)
ADD_PACK_PROTOTYPES(WellConnections)
ADD_PACK_PROTOTYPES(WellEconProductionLimits)
ADD_PACK_PROTOTYPES(WellFoamProperties)
ADD_PACK_PROTOTYPES(WellPolymerProperties)
ADD_PACK_PROTOTYPES(WellSegments)
ADD_PACK_PROTOTYPES(WellTestConfig)
ADD_PACK_PROTOTYPES(WellTestConfig::WTESTWell)
ADD_PACK_PROTOTYPES(WellTracerProperties)
ADD_PACK_PROTOTYPES(WList)
ADD_PACK_PROTOTYPES(WListManager)

} // end namespace Mpi

RestartValue loadParallelRestart(const EclipseIO* eclIO, SummaryState& summaryState,
                                 const std::vector<Opm::RestartKey>& solutionKeys,
                                 const std::vector<Opm::RestartKey>& extraKeys,
                                 Dune::CollectiveCommunication<Dune::MPIHelper::MPICommunicator> comm);

} // end namespace Opm
#endif // PARALLEL_RESTART_HPP
