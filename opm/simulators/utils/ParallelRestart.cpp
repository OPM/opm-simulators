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
#include <opm/parser/eclipse/EclipseState/Runspec.hpp>
#include <opm/parser/eclipse/EclipseState/Grid/NNC.hpp>
#include <opm/parser/eclipse/EclipseState/InitConfig/Equil.hpp>
#include <opm/parser/eclipse/EclipseState/InitConfig/FoamConfig.hpp>
#include <opm/parser/eclipse/EclipseState/InitConfig/InitConfig.hpp>
#include <opm/parser/eclipse/EclipseState/IOConfig/IOConfig.hpp>
#include <opm/parser/eclipse/EclipseState/IOConfig/RestartConfig.hpp>
#include <opm/parser/eclipse/EclipseState/Edit/EDITNNC.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/TimeMap.hpp>
#include <opm/parser/eclipse/EclipseState/SimulationConfig/SimulationConfig.hpp>
#include <opm/parser/eclipse/EclipseState/SimulationConfig/ThresholdPressure.hpp>
#include <opm/parser/eclipse/EclipseState/Tables/Aqudims.hpp>
#include <opm/parser/eclipse/EclipseState/Tables/ColumnSchema.hpp>
#include <opm/parser/eclipse/EclipseState/Tables/Eqldims.hpp>
#include <opm/parser/eclipse/EclipseState/Tables/FlatTable.hpp>
#include <opm/parser/eclipse/EclipseState/Tables/JFunc.hpp>
#include <opm/parser/eclipse/EclipseState/Tables/PlymwinjTable.hpp>
#include <opm/parser/eclipse/EclipseState/Tables/PolyInjTable.hpp>
#include <opm/parser/eclipse/EclipseState/Tables/PvtgTable.hpp>
#include <opm/parser/eclipse/EclipseState/Tables/PvtoTable.hpp>
#include <opm/parser/eclipse/EclipseState/Tables/Regdims.hpp>
#include <opm/parser/eclipse/EclipseState/Tables/Rock2dTable.hpp>
#include <opm/parser/eclipse/EclipseState/Tables/Rock2dtrTable.hpp>
#include <opm/parser/eclipse/EclipseState/Tables/SimpleTable.hpp>
#include <opm/parser/eclipse/EclipseState/Tables/SkprpolyTable.hpp>
#include <opm/parser/eclipse/EclipseState/Tables/SkprwatTable.hpp>
#include <opm/parser/eclipse/EclipseState/Tables/TableColumn.hpp>
#include <opm/parser/eclipse/EclipseState/Tables/TableContainer.hpp>
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

template<class T>
std::size_t packSize(const DynamicState<T>& data, Dune::MPIHelper::MPICommunicator comm)
{
    return packSize(data.data(), comm) + packSize(data.initialRange(), comm);
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

HANDLE_AS_POD(Actdims)
HANDLE_AS_POD(Aqudims)
HANDLE_AS_POD(data::Connection)
HANDLE_AS_POD(data::Rates)
HANDLE_AS_POD(data::Segment)
HANDLE_AS_POD(DENSITYRecord)
HANDLE_AS_POD(EclHysterConfig)
HANDLE_AS_POD(Eqldims)
HANDLE_AS_POD(EquilRecord)
HANDLE_AS_POD(FoamData)
HANDLE_AS_POD(JFunc)
HANDLE_AS_POD(PVTWRecord)
HANDLE_AS_POD(PVCDORecord)
HANDLE_AS_POD(Regdims)
HANDLE_AS_POD(RestartSchedule)
HANDLE_AS_POD(ROCKRecord)
HANDLE_AS_POD(Tabdims)
HANDLE_AS_POD(TimeMap::StepData)
HANDLE_AS_POD(VISCREFRecord)
HANDLE_AS_POD(WATDENTRecord)
HANDLE_AS_POD(Welldims)
HANDLE_AS_POD(WellSegmentDims)

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

std::size_t packSize(const SimpleTable& data, Dune::MPIHelper::MPICommunicator comm)
{
   return packSize(data.schema(), comm) +
          packSize(data.columns(), comm) +
          packSize(data.jfunc(), comm);
}

std::size_t packSize(const TableContainer& data, Dune::MPIHelper::MPICommunicator comm)
{
    size_t res = 2*packSize(data.max(), comm);
    for (const auto& it : data.tables()) {
        if (it.second) {
            res += packSize(it.first, comm) + packSize(*it.second, comm);
        }
    }

    return res;
}

std::size_t packSize(const Equil& data, Dune::MPIHelper::MPICommunicator comm)
{
    return packSize(data.records(), comm);
}

std::size_t packSize(const FoamConfig& data, Dune::MPIHelper::MPICommunicator comm)
{
    return packSize(data.records(), comm);
}

std::size_t packSize(const InitConfig& data, Dune::MPIHelper::MPICommunicator comm)
{
    return packSize(data.getEquil(), comm) +
           packSize(data.getFoamConfig(), comm) +
           packSize(data.filleps(), comm) +
           packSize(data.restartRequested(), comm) +
           packSize(data.getRestartStep(), comm) +
           packSize(data.getRestartRootName(), comm);
}

std::size_t packSize(const SimulationConfig& data, Dune::MPIHelper::MPICommunicator comm)
{
    return packSize(data.getThresholdPressure(), comm) +
           packSize(data.useCPR(), comm) +
           packSize(data.hasDISGAS(), comm) +
           packSize(data.hasVAPOIL(), comm) +
           packSize(data.isThermal(), comm);
}

std::size_t packSize(const TimeMap& data, Dune::MPIHelper::MPICommunicator comm)
{
    return packSize(data.timeList(), comm) +
           packSize(data.firstTimeStepMonths(), comm) +
           packSize(data.firstTimeStepYears(), comm);
}

std::size_t packSize(const RestartConfig& data, Dune::MPIHelper::MPICommunicator comm)
{
    return packSize(data.timeMap(), comm) +
           packSize(data.getFirstRestartStep(), comm) +
           packSize(data.writeInitialRst(), comm) +
           packSize(data.restartSchedule(), comm) +
           packSize(data.restartKeywords(), comm) +
           packSize(data.saveKeywords(), comm);
}

std::size_t packSize(const IOConfig& data, Dune::MPIHelper::MPICommunicator comm)
{
    return packSize(data.getWriteINITFile(), comm) +
           packSize(data.getWriteEGRIDFile(), comm) +
           packSize(data.getUNIFIN(), comm) +
           packSize(data.getUNIFOUT(), comm) +
           packSize(data.getFMTIN(), comm) +
           packSize(data.getFMTOUT(), comm) +
           packSize(data.getFirstRestartStep(), comm) +
           packSize(data.getDeckFileName(), comm) +
           packSize(data.getOutputEnabled(), comm) +
           packSize(data.getOutputDir(), comm) +
           packSize(data.getNoSim(), comm) +
           packSize(data.getBaseName(), comm) +
           packSize(data.getEclCompatibleRST(), comm);
}

std::size_t packSize(const Phases& data, Dune::MPIHelper::MPICommunicator comm)
{
    return packSize(data.getBits(), comm);
}

std::size_t packSize(const EndpointScaling& data, Dune::MPIHelper::MPICommunicator comm)
{
    return packSize(data.getBits(), comm);
}

std::size_t packSize(const UDQParams& data, Dune::MPIHelper::MPICommunicator comm)
{
    return packSize(data.reseed(), comm) +
           packSize(data.rand_seed(), comm) +
           packSize(data.range(), comm) +
           packSize(data.undefinedValue(), comm) +
           packSize(data.cmpEpsilon(), comm);
}

std::size_t packSize(const Runspec& data, Dune::MPIHelper::MPICommunicator comm)
{
    return packSize(data.phases(), comm) +
           packSize(data.tabdims(), comm) +
           packSize(data.endpointScaling(), comm) +
           packSize(data.wellDimensions(), comm) +
           packSize(data.wellSegmentDimensions(), comm) +
           packSize(data.udqParams(), comm) +
           packSize(data.hysterPar(), comm) +
           packSize(data.actdims(), comm);
}

std::size_t packSize(const PvtxTable& data, Dune::MPIHelper::MPICommunicator comm)
{
    return packSize(data.getOuterColumnSchema(), comm) +
           packSize(data.getOuterColumn(), comm) +
           packSize(data.getUnderSaturatedSchema(), comm) +
           packSize(data.getSaturatedSchema(), comm) +
           packSize(data.getUnderSaturatedTables(), comm) +
           packSize(data.getSaturatedTable(), comm);
}

std::size_t packSize(const PvtgTable& data, Dune::MPIHelper::MPICommunicator comm)
{
    return packSize(static_cast<const PvtxTable&>(data), comm);
}

std::size_t packSize(const PvtoTable& data, Dune::MPIHelper::MPICommunicator comm)
{
    return packSize(static_cast<const PvtxTable&>(data), comm);
}

std::size_t packSize(const PvtwTable& data, Dune::MPIHelper::MPICommunicator comm)
{
    return packSize(static_cast<const std::vector<PVTWRecord>&>(data), comm);
}

std::size_t packSize(const PvcdoTable& data, Dune::MPIHelper::MPICommunicator comm)
{
    return packSize(static_cast<const std::vector<PVCDORecord>&>(data), comm);
}

std::size_t packSize(const DensityTable& data, Dune::MPIHelper::MPICommunicator comm)
{
    return packSize(static_cast<const std::vector<DENSITYRecord>&>(data), comm);
}

std::size_t packSize(const ViscrefTable& data, Dune::MPIHelper::MPICommunicator comm)
{
    return packSize(static_cast<const std::vector<VISCREFRecord>&>(data), comm);
}

std::size_t packSize(const WatdentTable& data, Dune::MPIHelper::MPICommunicator comm)
{
    return packSize(static_cast<const std::vector<WATDENTRecord>&>(data), comm);
}

std::size_t packSize(const PolyInjTable& data, Dune::MPIHelper::MPICommunicator comm)
{
    return packSize(data.getThroughputs(), comm) +
           packSize(data.getVelocities(), comm) +
           packSize(data.getTableNumber(), comm) +
           packSize(data.getTableData(), comm);
}

std::size_t packSize(const PlymwinjTable& data, Dune::MPIHelper::MPICommunicator comm)
{
    return packSize(static_cast<const PolyInjTable&>(data), comm);
}

std::size_t packSize(const SkprpolyTable& data, Dune::MPIHelper::MPICommunicator comm)
{
    return packSize(static_cast<const PolyInjTable&>(data), comm) +
           packSize(data.referenceConcentration(), comm);
}

std::size_t packSize(const SkprwatTable& data, Dune::MPIHelper::MPICommunicator comm)
{
    return packSize(static_cast<const PolyInjTable&>(data), comm);
}

std::size_t packSize(const RockTable& data, Dune::MPIHelper::MPICommunicator comm)
{
    return packSize(static_cast<const std::vector<ROCKRecord>&>(data), comm);
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

template<class T>
void pack(const DynamicState<T>& data, std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    pack(data.data(), buffer, position, comm);
    pack(data.initialRange(), buffer, position, comm);
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

void pack(const SimpleTable& data, std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    pack(data.schema(), buffer, position, comm);
    pack(data.columns(), buffer, position, comm);
    pack(data.jfunc(), buffer, position, comm);
}

void pack(const TableContainer& data, std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    pack(data.max(), buffer, position, comm);
    size_t entries = 0;
    for (const auto& it : data.tables()) {
        if (it.second) {
          ++entries;
        }
    }
    pack(entries, buffer, position, comm);
    for (const auto& it : data.tables()) {
        if (it.second) {
          pack(it.first, buffer, position, comm);
          pack(*it.second, buffer, position, comm);
        }
    }
}

void pack(const Equil& data, std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    pack(data.records(), buffer, position, comm);
}

void pack(const FoamConfig& data, std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    pack(data.records(), buffer, position, comm);
}

void pack(const InitConfig& data, std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    pack(data.getEquil(), buffer, position, comm);
    pack(data.getFoamConfig(), buffer, position, comm);
    pack(data.filleps(), buffer, position, comm);
    pack(data.restartRequested(), buffer, position, comm);
    pack(data.getRestartStep(), buffer, position, comm);
    pack(data.getRestartRootName(), buffer, position, comm);
}

void pack(const SimulationConfig& data, std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    pack(data.getThresholdPressure(), buffer, position, comm);
    pack(data.useCPR(), buffer, position, comm);
    pack(data.hasDISGAS(), buffer, position, comm);
    pack(data.hasVAPOIL(), buffer, position, comm);
    pack(data.isThermal(), buffer, position, comm);
}

void pack(const TimeMap& data, std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    pack(data.timeList(), buffer, position, comm);
    pack(data.firstTimeStepMonths(), buffer, position, comm);
    pack(data.firstTimeStepYears(), buffer, position, comm);
}

void pack(const RestartConfig& data, std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    pack(data.timeMap(), buffer, position, comm);
    pack(data.getFirstRestartStep(), buffer, position, comm);
    pack(data.writeInitialRst(), buffer, position, comm);
    pack(data.restartSchedule(), buffer, position, comm);
    pack(data.restartKeywords(), buffer, position, comm);
    pack(data.saveKeywords(), buffer, position, comm);
}

void pack(const IOConfig& data, std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    pack(data.getWriteINITFile(), buffer, position, comm);
    pack(data.getWriteEGRIDFile(), buffer, position, comm);
    pack(data.getUNIFIN(), buffer, position, comm);
    pack(data.getUNIFOUT(), buffer, position, comm);
    pack(data.getFMTIN(), buffer, position, comm);
    pack(data.getFMTOUT(), buffer, position, comm);
    pack(data.getFirstRestartStep(), buffer, position, comm);
    pack(data.getDeckFileName(), buffer, position, comm);
    pack(data.getOutputEnabled(), buffer, position, comm);
    pack(data.getOutputDir(), buffer, position, comm);
    pack(data.getNoSim(), buffer, position, comm);
    pack(data.getBaseName(), buffer, position, comm);
    pack(data.getEclCompatibleRST(), buffer, position, comm);
}

void pack(const Phases& data, std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    pack(data.getBits(), buffer, position, comm);
}

void pack(const EndpointScaling& data, std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    pack(data.getBits(), buffer, position, comm);
}

void pack(const UDQParams& data, std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    pack(data.reseed(), buffer, position, comm);
    pack(data.rand_seed(), buffer, position, comm);
    pack(data.range(), buffer, position, comm);
    pack(data.undefinedValue(), buffer, position, comm);
    pack(data.cmpEpsilon(), buffer, position, comm);
}

void pack(const Runspec& data, std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    pack(data.phases(), buffer, position, comm);
    pack(data.tabdims(), buffer, position, comm);
    pack(data.endpointScaling(), buffer, position, comm);
    pack(data.wellDimensions(), buffer, position, comm);
    pack(data.wellSegmentDimensions(), buffer, position, comm);
    pack(data.udqParams(), buffer, position, comm);
    pack(data.hysterPar(), buffer, position, comm);
    pack(data.actdims(), buffer, position, comm);
}

void pack(const PvtxTable& data, std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    pack(data.getOuterColumnSchema(), buffer, position, comm);
    pack(data.getOuterColumn(), buffer, position, comm);
    pack(data.getUnderSaturatedSchema(), buffer, position, comm);
    pack(data.getSaturatedSchema(), buffer, position, comm);
    pack(data.getUnderSaturatedTables(), buffer, position, comm);
    pack(data.getSaturatedTable(), buffer, position, comm);
}

void pack(const PvtgTable& data, std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    pack(static_cast<const PvtxTable&>(data), buffer, position, comm);
}

void pack(const PvtoTable& data, std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    pack(static_cast<const PvtxTable&>(data), buffer, position, comm);
}

void pack(const PvtwTable& data, std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    pack(static_cast<const std::vector<PVTWRecord>&>(data), buffer, position, comm);
}

void pack(const PvcdoTable& data, std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    pack(static_cast<const std::vector<PVCDORecord>&>(data), buffer, position, comm);
}

void pack(const DensityTable& data, std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    pack(static_cast<const std::vector<DENSITYRecord>&>(data), buffer, position, comm);
}

void pack(const ViscrefTable& data, std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    pack(static_cast<const std::vector<VISCREFRecord>&>(data), buffer, position, comm);
}

void pack(const WatdentTable& data, std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    pack(static_cast<const std::vector<WATDENTRecord>&>(data), buffer, position, comm);
}

void pack(const PolyInjTable& data, std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    pack(data.getThroughputs(), buffer, position, comm);
    pack(data.getVelocities(), buffer, position, comm);
    pack(data.getTableNumber(), buffer, position, comm);
    pack(data.getTableData(), buffer, position, comm);
}

void pack(const PlymwinjTable& data, std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    pack(static_cast<const PolyInjTable&>(data), buffer, position, comm);
}

void pack(const SkprpolyTable& data, std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    pack(static_cast<const PolyInjTable&>(data), buffer, position, comm);
    pack(data.referenceConcentration(), buffer, position, comm);
}

void pack(const SkprwatTable& data, std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    pack(static_cast<const PolyInjTable&>(data), buffer, position, comm);
}

void pack(const RockTable& data, std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    pack(static_cast<const std::vector<ROCKRecord>&>(data), buffer, position, comm);
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

template<class T>
void unpack(DynamicState<T>& data, std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm)
{
    std::vector<T> ddata;
    size_t initial_range;
    unpack(ddata, buffer, position, comm);
    unpack(initial_range, buffer, position, comm);
    data = DynamicState<T>(ddata, initial_range);
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

void unpack(SimpleTable& data, std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm)
{
    TableSchema schema;
    OrderedMap<std::string, TableColumn> columns;
    bool jf;
    unpack(schema, buffer, position, comm);
    unpack(columns, buffer, position, comm);
    unpack(jf, buffer, position, comm);
    data = SimpleTable(schema, columns, jf);
}

void unpack(TableContainer& data, std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm)
{
    size_t max;
    unpack(max, buffer, position, comm);
    data = TableContainer(max);
    size_t entries;
    unpack(entries, buffer, position, comm);
    for (size_t i = 0; i < entries; ++i) {
        size_t id;
        unpack(id, buffer, position, comm);
        SimpleTable table;
        unpack(table, buffer, position, comm);
        data.addTable(id, std::make_shared<const SimpleTable>(table));
    }
}

void unpack(Equil& data, std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm)
{
    std::vector<EquilRecord> records;
    unpack(records, buffer, position, comm);
    data = Equil(records);
}

void unpack(FoamConfig& data, std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm)
{
    std::vector<FoamData> records;
    unpack(records, buffer, position, comm);
    data = FoamConfig(records);
}

void unpack(InitConfig& data, std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm)
{
    Equil equil;
    FoamConfig foam;
    bool filleps, restartRequested;
    int restartStep;
    std::string restartRootName;
    unpack(equil, buffer, position, comm);
    unpack(foam, buffer, position, comm);
    unpack(filleps, buffer, position, comm);
    unpack(restartRequested, buffer, position, comm);
    unpack(restartStep, buffer, position, comm);
    unpack(restartRootName, buffer, position, comm);
    data = InitConfig(equil, foam, filleps, restartRequested,
                      restartStep, restartRootName);
}

void unpack(SimulationConfig& data, std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm)
{
    ThresholdPressure thresholdPressure;
    bool useCPR, DISGAS, VAPOIL, isThermal;
    unpack(thresholdPressure, buffer, position, comm);
    unpack(useCPR, buffer, position, comm);
    unpack(DISGAS, buffer, position, comm);
    unpack(VAPOIL, buffer, position, comm);
    unpack(isThermal, buffer, position, comm);
    data = SimulationConfig(thresholdPressure, useCPR, DISGAS, VAPOIL, isThermal);
}

void unpack(TimeMap& data, std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm)
{
    std::vector<std::time_t> timeList;
    std::vector<TimeMap::StepData> firstStepMonths;
    std::vector<TimeMap::StepData> firstStepYears;
    unpack(timeList, buffer, position, comm);
    unpack(firstStepMonths, buffer, position, comm);
    unpack(firstStepYears, buffer, position, comm);

    data = TimeMap(timeList, firstStepMonths, firstStepYears);
}

void unpack(RestartConfig& data, std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm)
{
    TimeMap timemap;
    int firstRstStep;
    bool writeInitialRst;
    DynamicState<RestartSchedule> restart_sched;
    DynamicState<std::map<std::string,int>> restart_keyw;
    std::vector<bool> save_keyw;
    unpack(timemap, buffer, position, comm);
    unpack(firstRstStep, buffer, position, comm);
    unpack(writeInitialRst, buffer, position, comm);
    unpack(restart_sched, buffer, position, comm);
    unpack(restart_keyw, buffer, position, comm);
    unpack(save_keyw, buffer, position, comm);
    data = RestartConfig(timemap, firstRstStep, writeInitialRst, restart_sched,
                         restart_keyw, save_keyw);
}

void unpack(IOConfig& data, std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm)
{
    bool write_init, write_egrid, unifin, unifout, fmtin, fmtout;
    int firstRestartStep;
    std::string deck_name, output_dir, base_name;
    bool output_enabled, no_sim, ecl_compatible_rst;

    unpack(write_init, buffer, position, comm);
    unpack(write_egrid, buffer, position, comm);
    unpack(unifin, buffer, position, comm);
    unpack(unifout, buffer, position, comm);
    unpack(fmtin, buffer, position, comm);
    unpack(fmtout, buffer, position, comm);
    unpack(firstRestartStep, buffer, position, comm);
    unpack(deck_name, buffer, position, comm);
    unpack(output_enabled, buffer, position, comm);
    unpack(output_dir, buffer, position, comm);
    unpack(no_sim, buffer, position, comm);
    unpack(base_name, buffer, position, comm);
    unpack(ecl_compatible_rst, buffer, position, comm);
    data = IOConfig(write_init, write_egrid, unifin, unifout, fmtin, fmtout,
                    firstRestartStep, deck_name, output_enabled, output_dir,
                    no_sim, base_name, ecl_compatible_rst);
}

void unpack(Phases& data, std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm)
{
    unsigned long bits;
    unpack(bits, buffer, position, comm);
    data = Phases(std::bitset<NUM_PHASES_IN_ENUM>(bits));
}

void unpack(EndpointScaling& data, std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm)
{
    unsigned long bits;
    unpack(bits, buffer, position, comm);
    data = EndpointScaling(std::bitset<4>(bits));
}

void unpack(UDQParams& data, std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm)
{
    bool reseed;
    int rand_seed;
    double range, undefVal, cmp_eps;

    unpack(reseed, buffer, position, comm);
    unpack(rand_seed, buffer, position, comm);
    unpack(range, buffer, position, comm);
    unpack(undefVal, buffer, position, comm);
    unpack(cmp_eps, buffer, position, comm);
    data = UDQParams(reseed, rand_seed, range, undefVal, cmp_eps);
}

void unpack(Runspec& data, std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm)
{
    Phases phases;
    Tabdims tabdims;
    EndpointScaling endScale;
    Welldims wellDims;
    WellSegmentDims wsegDims;
    UDQParams udqparams;
    EclHysterConfig hystPar;
    Actdims actdims;
    unpack(phases, buffer, position, comm);
    unpack(tabdims, buffer, position, comm);
    unpack(endScale, buffer, position, comm);
    unpack(wellDims, buffer, position, comm);
    unpack(wsegDims, buffer, position, comm);
    unpack(udqparams, buffer, position, comm);
    unpack(hystPar, buffer, position, comm);
    unpack(actdims, buffer, position, comm);
    data = Runspec(phases, tabdims, endScale, wellDims, wsegDims,
                   udqparams, hystPar, actdims);
}

template<class PVTType>
void unpack_pvt(PVTType& data, std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm)
{
    ColumnSchema outer_schema;
    TableColumn outer_column;
    TableSchema undersat_schema, sat_schema;
    std::vector<SimpleTable> undersat_tables;
    SimpleTable sat_table;
    unpack(outer_schema, buffer, position, comm);
    unpack(outer_column, buffer, position, comm);
    unpack(undersat_schema, buffer, position, comm);
    unpack(sat_schema, buffer, position, comm);
    unpack(undersat_tables, buffer, position, comm);
    unpack(sat_table, buffer, position, comm);
    data = PVTType(outer_schema, outer_column, undersat_schema, sat_schema,
                   undersat_tables, sat_table);
}

void unpack(PvtgTable& data, std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm)
{
    unpack_pvt(data, buffer, position, comm);
}

void unpack(PvtoTable& data, std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm)
{
    unpack_pvt(data, buffer, position, comm);
}

void unpack(PvtwTable& data, std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm)
{
    std::vector<PVTWRecord> pdata;
    unpack(pdata, buffer, position, comm);
    data = PvtwTable(pdata);
}

void unpack(PvcdoTable& data, std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm)
{
    std::vector<PVCDORecord> pdata;
    unpack(pdata, buffer, position, comm);
    data = PvcdoTable(pdata);
}

void unpack(DensityTable& data, std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm)
{
    std::vector<DENSITYRecord> pdata;
    unpack(pdata, buffer, position, comm);
    data = DensityTable(pdata);
}

void unpack(ViscrefTable& data, std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm)
{
    std::vector<VISCREFRecord> pdata;
    unpack(pdata, buffer, position, comm);
    data = ViscrefTable(pdata);
}

void unpack(WatdentTable& data, std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm)
{
    std::vector<WATDENTRecord> pdata;
    unpack(pdata, buffer, position, comm);
    data = WatdentTable(pdata);
}

void unpack(PolyInjTable& data, std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm)
{
    std::vector<double> throughputs, velocities;
    int tableNumber;
    std::vector<std::vector<double>> tableData;
    unpack(throughputs, buffer, position, comm);
    unpack(velocities, buffer, position, comm);
    unpack(tableNumber, buffer, position, comm);
    unpack(tableData, buffer, position, comm);
    data = PolyInjTable(throughputs, velocities, tableNumber, tableData);
}

void unpack(PlymwinjTable& data, std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm)
{
    unpack(static_cast<PolyInjTable&>(data), buffer, position, comm);
}

void unpack(SkprpolyTable& data, std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm)
{
    unpack(static_cast<PolyInjTable&>(data), buffer, position, comm);
    double refConcentration;
    unpack(refConcentration, buffer, position, comm);
    data.setReferenceConcentration(refConcentration);
}

void unpack(SkprwatTable& data, std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm)
{
    unpack(static_cast<PolyInjTable&>(data), buffer, position, comm);
}

void unpack(RockTable& data, std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm)
{
    std::vector<ROCKRecord> pdata;
    unpack(pdata, buffer, position, comm);
    data = RockTable(pdata);
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
