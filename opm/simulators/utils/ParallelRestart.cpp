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
#include <opm/parser/eclipse/EclipseState/Schedule/OilVaporizationProperties.hpp>
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
#include <opm/parser/eclipse/EclipseState/Tables/TableManager.hpp>
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

template<std::size_t I = 0, typename Tuple>
typename std::enable_if<I == std::tuple_size<Tuple>::value, std::size_t>::type
pack_size_tuple_entry(const Tuple&, Dune::MPIHelper::MPICommunicator)
{
    return 0;
}

template<std::size_t I = 0, typename Tuple>
typename std::enable_if<I != std::tuple_size<Tuple>::value, std::size_t>::type
pack_size_tuple_entry(const Tuple& tuple, Dune::MPIHelper::MPICommunicator comm)
{
    return packSize(std::get<I>(tuple), comm) + pack_size_tuple_entry<I+1>(tuple, comm);
}

template<class... Ts>
std::size_t packSize(const std::tuple<Ts...>& data, Dune::MPIHelper::MPICommunicator comm)
{
    return pack_size_tuple_entry(data, comm);
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

std::size_t packSize(const TableManager& data, Dune::MPIHelper::MPICommunicator comm)
{
    return packSize(data.getSimpleTables(), comm) +
           packSize(data.getPvtgTables(), comm) +
           packSize(data.getPvtoTables(), comm) +
           packSize(data.getRock2dTables(), comm) +
           packSize(data.getRock2dtrTables(), comm) +
           packSize(data.getPvtwTable(), comm) +
           packSize(data.getPvcdoTable(), comm) +
           packSize(data.getDensityTable(), comm) +
           packSize(data.getRockTable(), comm) +
           packSize(data.getViscrefTable(), comm) +
           packSize(data.getWatdentTable(), comm) +
           packSize(data.getPlymwinjTables(), comm) +
           packSize(data.getSkprwatTables(), comm) +
           packSize(data.getSkprpolyTables(), comm) +
           packSize(data.getTabdims(), comm) +
           packSize(data.getRegdims(), comm) +
           packSize(data.getEqldims(), comm) +
           packSize(data.getAqudims(), comm) +
           packSize(data.useImptvd(), comm) +
           packSize(data.useEnptvd(), comm) +
           packSize(data.useEqlnum(), comm) +
           packSize(data.useJFunc(), comm) +
           (data.useJFunc() ? packSize(data.getJFunc(), comm) : 0) +
           packSize(data.rtemp(), comm);
}

template<class Scalar>
std::size_t packSize(const Tabulated1DFunction<Scalar>& data,
                     Dune::MPIHelper::MPICommunicator comm)
{
    return packSize(data.xValues(), comm) +
           packSize(data.yValues(), comm);
}

template std::size_t packSize(const Tabulated1DFunction<double>& data,
                              Dune::MPIHelper::MPICommunicator comm);

template<class Scalar>
std::size_t packSize(const SolventPvt<Scalar>& data,
                     Dune::MPIHelper::MPICommunicator comm)
{
    return packSize(data.solventReferenceDensity(), comm) +
           packSize(data.inverseSolventB(), comm) +
           packSize(data.solventMu(), comm) +
           packSize(data.inverseSolventBMu(), comm);
}

template std::size_t packSize(const SolventPvt<double>& data,
                              Dune::MPIHelper::MPICommunicator comm);

template<class Scalar>
std::size_t packSize(const IntervalTabulated2DFunction<Scalar>& data,
                     Dune::MPIHelper::MPICommunicator comm)
{
    return packSize(data.xPos(), comm) +
           packSize(data.yPos(), comm) +
           packSize(data.samples(), comm) +
           packSize(data.xExtrapolate(), comm) +
           packSize(data.yExtrapolate(), comm);
}

template std::size_t packSize(const IntervalTabulated2DFunction<double>& data,
                              Dune::MPIHelper::MPICommunicator comm);

template<class Scalar>
std::size_t packSize(const UniformXTabulated2DFunction<Scalar>& data,
                     Dune::MPIHelper::MPICommunicator comm)
{
    return packSize(data.xPos(), comm) +
           packSize(data.yPos(), comm) +
           packSize(data.samples(), comm) +
           packSize(data.interpolationGuide(), comm);
}

template std::size_t packSize(const UniformXTabulated2DFunction<double>& data,
                              Dune::MPIHelper::MPICommunicator comm);

template<class Scalar>
std::size_t packSize(const DryGasPvt<Scalar>& data,
                     Dune::MPIHelper::MPICommunicator comm)
{
    return packSize(data.gasReferenceDensity(), comm) +
           packSize(data.inverseGasB(), comm) +
           packSize(data.gasMu(), comm) +
           packSize(data.inverseGasBMu(), comm);
}

template std::size_t packSize(const DryGasPvt<double>& data,
                              Dune::MPIHelper::MPICommunicator comm);

template<class Scalar>
std::size_t packSize(const GasPvtThermal<Scalar>& data,
                     Dune::MPIHelper::MPICommunicator comm)
{
    std::size_t size = packSize(data.gasvisctCurves(), comm) +
                       packSize(data.gasdentRefTemp(), comm) +
                       packSize(data.gasdentCT1(), comm) +
                       packSize(data.gasdentCT2(), comm) +
                       packSize(data.internalEnergyCurves(), comm) +
                       packSize(data.enableThermalDensity(), comm) +
                       packSize(data.enableThermalViscosity(), comm) +
                       packSize(data.enableInternalEnergy(), comm);
    size += packSize(bool(), comm);
    if (data.isoThermalPvt())
        size += packSize(*data.isoThermalPvt(), comm);

    return size;
}

template std::size_t packSize(const GasPvtThermal<double>& data,
                              Dune::MPIHelper::MPICommunicator comm);

template<class Scalar, bool enableThermal>
std::size_t packSize(const GasPvtMultiplexer<Scalar,enableThermal>& data,
                     Dune::MPIHelper::MPICommunicator comm)
{
    std::size_t size = packSize(data.gasPvtApproach(), comm);
    const void* realGasPvt = data.realGasPvt();
    using PvtApproach = GasPvtMultiplexer<Scalar,enableThermal>;
    if (data.gasPvtApproach() == PvtApproach::DryGasPvt) {
        const auto& pvt = *static_cast<const DryGasPvt<Scalar>*>(realGasPvt);
        size += packSize(pvt, comm);
    } else if (data.gasPvtApproach() == PvtApproach::WetGasPvt) {
        const auto& pvt = *static_cast<const WetGasPvt<Scalar>*>(realGasPvt);
        size += packSize(pvt, comm);
    } else if (data.gasPvtApproach() == PvtApproach::ThermalGasPvt) {
        const auto& pvt = *static_cast<const GasPvtThermal<Scalar>*>(realGasPvt);
        size += packSize(pvt, comm);
    }

    return size;
}

template std::size_t packSize(const GasPvtMultiplexer<double,true>& data,
                              Dune::MPIHelper::MPICommunicator comm);
template std::size_t packSize(const GasPvtMultiplexer<double,false>& data,
                              Dune::MPIHelper::MPICommunicator comm);

template<class Scalar>
std::size_t packSize(const WetGasPvt<Scalar>& data,
                     Dune::MPIHelper::MPICommunicator comm)
{
    return packSize(data.gasReferenceDensity(), comm) +
           packSize(data.oilReferenceDensity(), comm) +
           packSize(data.inverseGasB(), comm) +
           packSize(data.inverseSaturatedGasB(), comm) +
           packSize(data.gasMu(), comm) +
           packSize(data.inverseGasBMu(), comm) +
           packSize(data.inverseSaturatedGasBMu(), comm) +
           packSize(data.saturatedOilVaporizationFactorTable(), comm) +
           packSize(data.saturationPressure(), comm) +
           packSize(data.vapPar1(), comm);
}

template std::size_t packSize(const WetGasPvt<double>& data,
                              Dune::MPIHelper::MPICommunicator comm);

template<class Scalar, bool enableThermal>
std::size_t packSize(const OilPvtMultiplexer<Scalar,enableThermal>& data,
                     Dune::MPIHelper::MPICommunicator comm)
{
    std::size_t size = packSize(data.approach(), comm);
    const void* realOilPvt = data.realOilPvt();
    using PvtApproach = OilPvtMultiplexer<Scalar,enableThermal>;
    if (data.approach() == PvtApproach::ConstantCompressibilityOilPvt) {
        const auto& pvt = *static_cast<const ConstantCompressibilityOilPvt<Scalar>*>(realOilPvt);
        size += packSize(pvt, comm);
    } else if (data.approach() == PvtApproach::DeadOilPvt) {
        const auto& pvt = *static_cast<const DeadOilPvt<Scalar>*>(realOilPvt);
        size += packSize(pvt, comm);
    } else if (data.approach() == PvtApproach::LiveOilPvt) {
        const auto& pvt = *static_cast<const LiveOilPvt<Scalar>*>(realOilPvt);
        size += packSize(pvt, comm);
    } else if (data.approach() == PvtApproach::ThermalOilPvt) {
        const auto& pvt = *static_cast<const OilPvtThermal<Scalar>*>(realOilPvt);
        size += packSize(pvt, comm);
    }

    return size;
}

template std::size_t packSize(const OilPvtMultiplexer<double,true>& data,
                              Dune::MPIHelper::MPICommunicator comm);
template std::size_t packSize(const OilPvtMultiplexer<double,false>& data,
                              Dune::MPIHelper::MPICommunicator comm);

template<class Scalar>
std::size_t packSize(const ConstantCompressibilityOilPvt<Scalar>& data,
                     Dune::MPIHelper::MPICommunicator comm)
{
    return packSize(data.oilReferenceDensity(), comm) +
           packSize(data.oilReferencePressure(), comm) +
           packSize(data.oilReferenceFormationVolumeFactor(), comm) +
           packSize(data.oilCompressibility(), comm) +
           packSize(data.oilViscosity(), comm) +
           packSize(data.oilViscosibility(), comm);
}

template std::size_t packSize(const ConstantCompressibilityOilPvt<double>& data,
                              Dune::MPIHelper::MPICommunicator comm);

template<class Scalar>
std::size_t packSize(const DeadOilPvt<Scalar>& data,
                     Dune::MPIHelper::MPICommunicator comm)
{
    return packSize(data.oilReferenceDensity(), comm) +
           packSize(data.inverseOilB(), comm) +
           packSize(data.oilMu(), comm) +
           packSize(data.inverseOilBMu(), comm);
}

template std::size_t packSize(const DeadOilPvt<double>& data,
                              Dune::MPIHelper::MPICommunicator comm);

template<class Scalar>
std::size_t packSize(const LiveOilPvt<Scalar>& data,
                     Dune::MPIHelper::MPICommunicator comm)
{
    return packSize(data.gasReferenceDensity(), comm) +
           packSize(data.oilReferenceDensity(), comm) +
           packSize(data.inverseOilBTable(), comm) +
           packSize(data.oilMuTable(), comm) +
           packSize(data.inverseOilBMuTable(), comm) +
           packSize(data.saturatedOilMuTable(), comm) +
           packSize(data.inverseSaturatedOilBTable(), comm) +
           packSize(data.inverseSaturatedOilBMuTable(), comm) +
           packSize(data.saturatedGasDissolutionFactorTable(), comm) +
           packSize(data.saturationPressure(), comm) +
           packSize(data.vapPar2(), comm);
}

template std::size_t packSize(const LiveOilPvt<double>& data,
                              Dune::MPIHelper::MPICommunicator comm);

template<class Scalar>
std::size_t packSize(const OilPvtThermal<Scalar>& data,
                     Dune::MPIHelper::MPICommunicator comm)
{
    std::size_t size = packSize(data.oilvisctCurves(), comm) +
                       packSize(data.viscrefPress(), comm) +
                       packSize(data.viscrefRs(), comm) +
                       packSize(data.viscRef(), comm) +
                       packSize(data.oildentRefTemp(), comm) +
                       packSize(data.oildentCT1(), comm) +
                       packSize(data.oildentCT2(), comm) +
                       packSize(data.internalEnergyCurves(), comm) +
                       packSize(data.enableThermalDensity(), comm) +
                       packSize(data.enableThermalViscosity(), comm) +
                       packSize(data.enableInternalEnergy(), comm);
    size += packSize(bool(), comm);
    if (data.isoThermalPvt())
        size += packSize(*data.isoThermalPvt(), comm);

    return size;
}

template std::size_t packSize(const OilPvtThermal<double>& data,
                              Dune::MPIHelper::MPICommunicator comm);

template<class Scalar, bool enableThermal>
std::size_t packSize(const WaterPvtMultiplexer<Scalar,enableThermal>& data,
                     Dune::MPIHelper::MPICommunicator comm)
{
    std::size_t size = packSize(data.approach(), comm);
    const void* realWaterPvt = data.realWaterPvt();
    using PvtApproach = WaterPvtMultiplexer<Scalar,enableThermal>;
    if (data.approach() == PvtApproach::ConstantCompressibilityWaterPvt) {
        const auto& pvt = *static_cast<const ConstantCompressibilityWaterPvt<Scalar>*>(realWaterPvt);
        size += packSize(pvt, comm);
    } else if (data.approach() == PvtApproach::ThermalWaterPvt) {
        const auto& pvt = *static_cast<const WaterPvtThermal<Scalar>*>(realWaterPvt);
        size += packSize(pvt, comm);
    }

    return size;
}

template std::size_t packSize(const WaterPvtMultiplexer<double,true>& data,
                              Dune::MPIHelper::MPICommunicator comm);
template std::size_t packSize(const WaterPvtMultiplexer<double,false>& data,
                              Dune::MPIHelper::MPICommunicator comm);

template<class Scalar>
std::size_t packSize(const ConstantCompressibilityWaterPvt<Scalar>& data,
                     Dune::MPIHelper::MPICommunicator comm)
{
    return packSize(data.waterReferenceDensity(), comm) +
           packSize(data.waterReferencePressure(), comm) +
           packSize(data.waterReferenceFormationVolumeFactor(), comm) +
           packSize(data.waterCompressibility(), comm) +
           packSize(data.waterViscosity(), comm) +
           packSize(data.waterViscosibility(), comm);
}

template std::size_t packSize(const ConstantCompressibilityWaterPvt<double>& data,
                              Dune::MPIHelper::MPICommunicator comm);

template<class Scalar>
std::size_t packSize(const WaterPvtThermal<Scalar>& data,
                     Dune::MPIHelper::MPICommunicator comm)
{
    std::size_t size = packSize(data.viscrefPress(), comm) +
                       packSize(data.watdentRefTemp(), comm) +
                       packSize(data.watdentCT1(), comm) +
                       packSize(data.watdentCT2(), comm) +
                       packSize(data.pvtwRefPress(), comm) +
                       packSize(data.pvtwRefB(), comm) +
                       packSize(data.pvtwCompressibility(), comm) +
                       packSize(data.pvtwViscosity(), comm) +
                       packSize(data.pvtwViscosibility(), comm) +
                       packSize(data.watvisctCurves(), comm) +
                       packSize(data.internalEnergyCurves(), comm) +
                       packSize(data.enableThermalDensity(), comm) +
                       packSize(data.enableThermalViscosity(), comm) +
                       packSize(data.enableInternalEnergy(), comm);
    size += packSize(bool(), comm);
    if (data.isoThermalPvt())
        size += packSize(*data.isoThermalPvt(), comm);

    return size;
}

template std::size_t packSize(const WaterPvtThermal<double>& data,
                              Dune::MPIHelper::MPICommunicator comm);

std::size_t packSize(const OilVaporizationProperties& data,
                     Dune::MPIHelper::MPICommunicator comm)
{
    return packSize(data.getType(), comm) +
           packSize(data.vap1(), comm) +
           packSize(data.vap2(), comm) +
           packSize(data.maxDRSDT(), comm) +
           packSize(data.maxDRSDT_allCells(), comm) +
           packSize(data.maxDRVDT(), comm);
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

template<std::size_t I = 0, typename Tuple>
typename std::enable_if<I == std::tuple_size<Tuple>::value, void>::type
pack_tuple_entry(const Tuple&, std::vector<char>&, int&,
                      Dune::MPIHelper::MPICommunicator)
{
}

template<std::size_t I = 0, typename Tuple>
typename std::enable_if<I != std::tuple_size<Tuple>::value, void>::type
pack_tuple_entry(const Tuple& tuple, std::vector<char>& buffer,
                 int& position, Dune::MPIHelper::MPICommunicator comm)
{
    pack(std::get<I>(tuple), buffer, position, comm);
    pack_tuple_entry<I+1>(tuple, buffer, position, comm);
}

template<class... Ts>
void pack(const std::tuple<Ts...>& data, std::vector<char>& buffer,
          int& position, Dune::MPIHelper::MPICommunicator comm)
{
    pack_tuple_entry(data, buffer, position, comm);
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

void pack(const TableManager& data, std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    pack(data.getSimpleTables(), buffer, position, comm);
    pack(data.getPvtgTables(), buffer, position, comm);
    pack(data.getPvtoTables(), buffer, position, comm);
    pack(data.getRock2dTables(), buffer, position, comm);
    pack(data.getRock2dtrTables(), buffer, position, comm);
    pack(data.getPvtwTable(), buffer, position, comm);
    pack(data.getPvcdoTable(), buffer, position, comm);
    pack(data.getDensityTable(), buffer, position, comm);
    pack(data.getRockTable(), buffer, position, comm);
    pack(data.getViscrefTable(), buffer, position, comm);
    pack(data.getWatdentTable(), buffer, position, comm);
    pack(data.getPlymwinjTables(), buffer, position, comm);
    pack(data.getSkprwatTables(), buffer, position, comm);
    pack(data.getSkprpolyTables(), buffer, position, comm);
    pack(data.getTabdims(), buffer, position, comm);
    pack(data.getRegdims(), buffer, position, comm);
    pack(data.getEqldims(), buffer, position, comm);
    pack(data.getAqudims(), buffer, position, comm);
    pack(data.useImptvd(), buffer, position, comm);
    pack(data.useEnptvd(), buffer, position, comm);
    pack(data.useEqlnum(), buffer, position, comm);
    pack(data.useJFunc(), buffer, position, comm);
    if (data.useJFunc())
        pack(data.getJFunc(), buffer, position, comm);
    pack(data.rtemp(), buffer, position, comm);
}

template<class Scalar>
void pack(const Tabulated1DFunction<Scalar>& data, std::vector<char>& buffer,
          int& position, Dune::MPIHelper::MPICommunicator comm)
{
    pack(data.xValues(), buffer, position, comm);
    pack(data.yValues(), buffer, position, comm);
}

template void pack(const Tabulated1DFunction<double>& data, std::vector<char>& buffer,
                   int& position, Dune::MPIHelper::MPICommunicator comm);

template<class Scalar>
void pack(const IntervalTabulated2DFunction<Scalar>& data, std::vector<char>& buffer,
          int& position, Dune::MPIHelper::MPICommunicator comm)
{
    pack(data.xPos(), buffer, position, comm);
    pack(data.yPos(), buffer, position, comm);
    pack(data.samples(), buffer, position, comm);
    pack(data.xExtrapolate(), buffer, position, comm);
    pack(data.yExtrapolate(), buffer, position, comm);
}

template void pack(const IntervalTabulated2DFunction<double>& data,
                   std::vector<char>& buffer,
                   int& position, Dune::MPIHelper::MPICommunicator comm);

template
void pack(const std::vector<IntervalTabulated2DFunction<double>>& data,
          std::vector<char>& buffer,
          int& position, Dune::MPIHelper::MPICommunicator comm);

template
void pack(const std::map<int,IntervalTabulated2DFunction<double>>& data,
          std::vector<char>& buffer,
          int& position, Dune::MPIHelper::MPICommunicator comm);

template<class Scalar>
void pack(const UniformXTabulated2DFunction<Scalar>& data, std::vector<char>& buffer,
          int& position, Dune::MPIHelper::MPICommunicator comm)
{
    pack(data.xPos(), buffer, position, comm);
    pack(data.yPos(), buffer, position, comm);
    pack(data.samples(), buffer, position, comm);
    pack(data.interpolationGuide(), buffer, position, comm);
}

template void pack(const UniformXTabulated2DFunction<double>& data,
                   std::vector<char>& buffer,
                   int& position, Dune::MPIHelper::MPICommunicator comm);

template<class Scalar>
void pack(const SolventPvt<Scalar>& data, std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    pack(data.solventReferenceDensity(), buffer, position, comm);
    pack(data.inverseSolventB(), buffer, position, comm);
    pack(data.solventMu(), buffer, position, comm);
    pack(data.inverseSolventBMu(), buffer, position, comm);
}

template void pack(const SolventPvt<double>& data,
                   std::vector<char>& buffer, int& position,
                   Dune::MPIHelper::MPICommunicator comm);

template<class Scalar, bool enableThermal>
void pack(const GasPvtMultiplexer<Scalar,enableThermal>& data,
          std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    pack(data.gasPvtApproach(), buffer, position, comm);
    const void* realGasPvt = data.realGasPvt();
    using PvtApproach = GasPvtMultiplexer<Scalar,enableThermal>;
    if (data.gasPvtApproach() == PvtApproach::DryGasPvt) {
        const auto& pvt = *static_cast<const DryGasPvt<Scalar>*>(realGasPvt);
        pack(pvt, buffer, position, comm);
    } else if (data.gasPvtApproach() == PvtApproach::WetGasPvt) {
        const auto& pvt = *static_cast<const WetGasPvt<Scalar>*>(realGasPvt);
        pack(pvt, buffer, position, comm);
    } else if (data.gasPvtApproach() == PvtApproach::ThermalGasPvt) {
        const auto& pvt = *static_cast<const GasPvtThermal<Scalar>*>(realGasPvt);
        pack(pvt, buffer, position, comm);
    }
}

template void pack(const GasPvtMultiplexer<double,true>& data,
                   std::vector<char>& buffer, int& position,
                   Dune::MPIHelper::MPICommunicator comm);
template void pack(const GasPvtMultiplexer<double,false>& data,
                   std::vector<char>& buffer, int& position,
                   Dune::MPIHelper::MPICommunicator comm);

template<class Scalar>
void pack(const DryGasPvt<Scalar>& data, std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    pack(data.gasReferenceDensity(), buffer, position, comm);
    pack(data.inverseGasB(), buffer, position, comm);
    pack(data.gasMu(), buffer, position, comm);
    pack(data.inverseGasBMu(), buffer, position, comm);
}

template void pack(const DryGasPvt<double>& data,
                   std::vector<char>& buffer, int& position,
                   Dune::MPIHelper::MPICommunicator comm);

template<class Scalar>
void pack(const GasPvtThermal<Scalar>& data,
          std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    pack(data.gasvisctCurves(), buffer, position, comm);
    pack(data.gasdentRefTemp(), buffer, position, comm);
    pack(data.gasdentCT1(), buffer, position, comm);
    pack(data.gasdentCT2(), buffer, position, comm);
    pack(data.internalEnergyCurves(), buffer, position, comm);
    pack(data.enableThermalDensity(), buffer, position, comm);
    pack(data.enableThermalViscosity(), buffer, position, comm);
    pack(data.enableInternalEnergy(), buffer, position, comm);
    pack(data.isoThermalPvt() != nullptr, buffer, position, comm);
    if (data.isoThermalPvt())
        pack(*data.isoThermalPvt(), buffer, position, comm);
}

template void pack(const GasPvtThermal<double>& data,
                   std::vector<char>& buffer, int& position,
                   Dune::MPIHelper::MPICommunicator comm);

template<class Scalar>
void pack(const WetGasPvt<Scalar>& data, std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    pack(data.gasReferenceDensity(), buffer, position, comm);
    pack(data.oilReferenceDensity(), buffer, position, comm);
    pack(data.inverseGasB(), buffer, position, comm);
    pack(data.inverseSaturatedGasB(), buffer, position, comm);
    pack(data.gasMu(), buffer, position, comm);
    pack(data.inverseGasBMu(), buffer, position, comm);
    pack(data.inverseSaturatedGasBMu(), buffer, position, comm);
    pack(data.saturatedOilVaporizationFactorTable(), buffer, position, comm);
    pack(data.saturationPressure(), buffer, position, comm);
    pack(data.vapPar1(), buffer, position, comm);
}

template void pack(const WetGasPvt<double>& data,
                   std::vector<char>& buffer, int& position,
                   Dune::MPIHelper::MPICommunicator comm);

template<class Scalar, bool enableThermal>
void pack(const OilPvtMultiplexer<Scalar,enableThermal>& data,
          std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    pack(data.approach(), buffer, position, comm);
    const void* realOilPvt = data.realOilPvt();
    using PvtApproach = OilPvtMultiplexer<Scalar,enableThermal>;
    if (data.approach() == PvtApproach::ConstantCompressibilityOilPvt) {
        const auto& pvt = *static_cast<const ConstantCompressibilityOilPvt<Scalar>*>(realOilPvt);
        pack(pvt, buffer, position, comm);
    } else if (data.approach() == PvtApproach::DeadOilPvt) {
        const auto& pvt = *static_cast<const DeadOilPvt<Scalar>*>(realOilPvt);
        pack(pvt, buffer, position, comm);
    } else if (data.approach() == PvtApproach::LiveOilPvt) {
        const auto& pvt = *static_cast<const LiveOilPvt<Scalar>*>(realOilPvt);
        pack(pvt, buffer, position, comm);
    } else if (data.approach() == PvtApproach::ThermalOilPvt) {
        const auto& pvt = *static_cast<const OilPvtThermal<Scalar>*>(realOilPvt);
        pack(pvt, buffer, position, comm);
    }
}

template void pack(const OilPvtMultiplexer<double,true>& data,
                   std::vector<char>& buffer, int& position,
                   Dune::MPIHelper::MPICommunicator comm);
template void pack(const OilPvtMultiplexer<double,false>& data,
                   std::vector<char>& buffer, int& position,
                   Dune::MPIHelper::MPICommunicator comm);

template<class Scalar>
void pack(const ConstantCompressibilityOilPvt<Scalar>& data,
          std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    pack(data.oilReferenceDensity(), buffer, position, comm);
    pack(data.oilReferencePressure(), buffer, position, comm);
    pack(data.oilReferenceFormationVolumeFactor(), buffer, position, comm);
    pack(data.oilCompressibility(), buffer, position, comm);
    pack(data.oilViscosity(), buffer, position, comm);
    pack(data.oilViscosibility(), buffer, position, comm);
}

template void pack(const ConstantCompressibilityOilPvt<double>& data,
                   std::vector<char>& buffer, int& position,
                   Dune::MPIHelper::MPICommunicator comm);

template<class Scalar>
void pack(const DeadOilPvt<Scalar>& data,
          std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    pack(data.oilReferenceDensity(), buffer, position, comm);
    pack(data.inverseOilB(), buffer, position, comm);
    pack(data.oilMu(), buffer, position, comm);
    pack(data.inverseOilBMu(), buffer, position, comm);
}

template void pack(const DeadOilPvt<double>& data,
                   std::vector<char>& buffer, int& position,
                   Dune::MPIHelper::MPICommunicator comm);

template<class Scalar>
void pack(const LiveOilPvt<Scalar>& data,
          std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    pack(data.gasReferenceDensity(), buffer, position, comm);
    pack(data.oilReferenceDensity(), buffer, position, comm);
    pack(data.inverseOilBTable(), buffer, position, comm);
    pack(data.oilMuTable(), buffer, position, comm);
    pack(data.inverseOilBMuTable(), buffer, position, comm);
    pack(data.saturatedOilMuTable(), buffer, position, comm);
    pack(data.inverseSaturatedOilBTable(), buffer, position, comm);
    pack(data.inverseSaturatedOilBMuTable(), buffer, position, comm);
    pack(data.saturatedGasDissolutionFactorTable(), buffer, position, comm);
    pack(data.saturationPressure(), buffer, position, comm);
    pack(data.vapPar2(), buffer, position, comm);
}

template void pack(const LiveOilPvt<double>& data,
                   std::vector<char>& buffer, int& position,
                   Dune::MPIHelper::MPICommunicator comm);

template<class Scalar>
void pack(const OilPvtThermal<Scalar>& data,
          std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    pack(data.oilvisctCurves(), buffer, position, comm);
    pack(data.viscrefPress(), buffer, position, comm);
    pack(data.viscrefRs(), buffer, position, comm);
    pack(data.viscRef(), buffer, position, comm);
    pack(data.oildentRefTemp(), buffer, position, comm);
    pack(data.oildentCT1(), buffer, position, comm);
    pack(data.oildentCT2(), buffer, position, comm);
    pack(data.internalEnergyCurves(), buffer, position, comm);
    pack(data.enableThermalDensity(), buffer, position, comm);
    pack(data.enableThermalViscosity(), buffer, position, comm);
    pack(data.enableInternalEnergy(), buffer, position, comm);
    pack(data.isoThermalPvt() != nullptr, buffer, position, comm);
    if (data.isoThermalPvt())
        pack(*data.isoThermalPvt(), buffer, position, comm);
}

template void pack(const OilPvtThermal<double>& data,
                   std::vector<char>& buffer, int& position,
                   Dune::MPIHelper::MPICommunicator comm);

template<class Scalar, bool enableThermal>
void pack(const WaterPvtMultiplexer<Scalar,enableThermal>& data,
          std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    pack(data.approach(), buffer, position, comm);
    const void* realWaterPvt = data.realWaterPvt();
    using PvtApproach = WaterPvtMultiplexer<Scalar,enableThermal>;
    if (data.approach() == PvtApproach::ConstantCompressibilityWaterPvt) {
        const auto& pvt = *static_cast<const ConstantCompressibilityWaterPvt<Scalar>*>(realWaterPvt);
        pack(pvt, buffer, position, comm);
    } else if (data.approach() == PvtApproach::ThermalWaterPvt) {
        const auto& pvt = *static_cast<const WaterPvtThermal<Scalar>*>(realWaterPvt);
        pack(pvt, buffer, position, comm);
    }
}

template void pack(const WaterPvtMultiplexer<double,true>& data,
                   std::vector<char>& buffer, int& position,
                   Dune::MPIHelper::MPICommunicator comm);
template void pack(const WaterPvtMultiplexer<double,false>& data,
                   std::vector<char>& buffer, int& position,
                   Dune::MPIHelper::MPICommunicator comm);

template<class Scalar>
void pack(const ConstantCompressibilityWaterPvt<Scalar>& data,
          std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    pack(data.waterReferenceDensity(), buffer, position, comm);
    pack(data.waterReferencePressure(), buffer, position, comm);
    pack(data.waterReferenceFormationVolumeFactor(), buffer, position, comm);
    pack(data.waterCompressibility(), buffer, position, comm);
    pack(data.waterViscosity(), buffer, position, comm);
    pack(data.waterViscosibility(), buffer, position, comm);
}

template void pack(const ConstantCompressibilityWaterPvt<double>& data,
                   std::vector<char>& buffer, int& position,
                   Dune::MPIHelper::MPICommunicator comm);

template<class Scalar>
void pack(const WaterPvtThermal<Scalar>& data,
          std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    pack(data.viscrefPress(), buffer, position, comm);
    pack(data.watdentRefTemp(), buffer, position, comm);
    pack(data.watdentCT1(), buffer, position, comm);
    pack(data.watdentCT2(), buffer, position, comm);
    pack(data.pvtwRefPress(), buffer, position, comm);
    pack(data.pvtwRefB(), buffer, position, comm);
    pack(data.pvtwCompressibility(), buffer, position, comm);
    pack(data.pvtwViscosity(), buffer, position, comm);
    pack(data.pvtwViscosibility(), buffer, position, comm);
    pack(data.watvisctCurves(), buffer, position, comm);
    pack(data.internalEnergyCurves(), buffer, position, comm);
    pack(data.enableThermalDensity(), buffer, position, comm);
    pack(data.enableThermalViscosity(), buffer, position, comm);
    pack(data.enableInternalEnergy(), buffer, position, comm);
    pack(data.isoThermalPvt() != nullptr, buffer, position, comm);
    if (data.isoThermalPvt())
        pack(*data.isoThermalPvt(), buffer, position, comm);
}

template void pack(const WaterPvtThermal<double>& data,
                   std::vector<char>& buffer, int& position,
                   Dune::MPIHelper::MPICommunicator comm);

void pack(const OilVaporizationProperties& data,
          std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    pack(data.getType(), buffer, position, comm);
    pack(data.vap1(), buffer, position, comm);
    pack(data.vap2(), buffer, position, comm);
    pack(data.maxDRSDT(), buffer, position, comm);
    pack(data.maxDRSDT_allCells(), buffer, position, comm);
    pack(data.maxDRVDT(), buffer, position, comm);
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

template<std::size_t I = 0, typename Tuple>
typename std::enable_if<I == std::tuple_size<Tuple>::value, void>::type
unpack_tuple_entry(Tuple&, std::vector<char>&, int&,
                   Dune::MPIHelper::MPICommunicator)
{
}

template<std::size_t I = 0, typename Tuple>
typename std::enable_if<I != std::tuple_size<Tuple>::value, void>::type
unpack_tuple_entry(Tuple& tuple, std::vector<char>& buffer,
                   int& position, Dune::MPIHelper::MPICommunicator comm)
{
    unpack(std::get<I>(tuple), buffer, position, comm);
    unpack_tuple_entry<I+1>(tuple, buffer, position, comm);
}

template<class... Ts>
void unpack(std::tuple<Ts...>& data, std::vector<char>& buffer,
            int& position, Dune::MPIHelper::MPICommunicator comm)
{
    unpack_tuple_entry(data, buffer, position, comm);
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

void unpack(TableManager& data, std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    std::map<std::string, TableContainer> simpleTables;
    std::vector<PvtgTable> pvtgTables;
    std::vector<PvtoTable> pvtoTables;
    std::vector<Rock2dTable> rock2dTables;
    std::vector<Rock2dtrTable> rock2dtrTables;
    PvtwTable pvtwTable;
    PvcdoTable pvcdoTable;
    DensityTable densityTable;
    RockTable rockTable;
    ViscrefTable viscrefTable;
    WatdentTable watdentTable;
    std::map<int, PlymwinjTable> plymwinjTables;
    std::map<int, SkprwatTable> skprwatTables;
    std::map<int, SkprpolyTable> skprpolyTables;
    Tabdims tabdims;
    Regdims regdims;
    Eqldims eqldims;
    Aqudims aqudims;
    bool hasImptvd;
    bool hasEntpvd;
    bool hasEqlnum;
    std::shared_ptr<JFunc> jfunc;
    double rtemp;
    unpack(simpleTables, buffer, position, comm);
    unpack(pvtgTables, buffer, position, comm);
    unpack(pvtoTables, buffer, position, comm);
    unpack(rock2dTables, buffer, position, comm);
    unpack(rock2dtrTables, buffer, position, comm);
    unpack(pvtwTable, buffer, position, comm);
    unpack(pvcdoTable, buffer, position, comm);
    unpack(densityTable, buffer, position, comm);
    unpack(rockTable, buffer, position, comm);
    unpack(viscrefTable, buffer, position, comm);
    unpack(watdentTable, buffer, position, comm);
    unpack(plymwinjTables, buffer, position, comm);
    unpack(skprwatTables, buffer, position, comm);
    unpack(skprpolyTables, buffer, position, comm);
    unpack(tabdims, buffer, position, comm);
    unpack(regdims, buffer, position, comm);
    unpack(eqldims, buffer, position, comm);
    unpack(aqudims, buffer, position, comm);
    unpack(hasImptvd, buffer, position, comm);
    unpack(hasEntpvd, buffer, position, comm);
    unpack(hasEqlnum, buffer, position, comm);
    bool hasJf;
    unpack(hasJf, buffer, position, comm);
    if (hasJf) {
        jfunc = std::make_shared<JFunc>();
        unpack(*jfunc, buffer, position, comm);
    }
    unpack(rtemp, buffer, position, comm);
    data = TableManager(simpleTables, pvtgTables, pvtoTables, rock2dTables,
                        rock2dtrTables, pvtwTable, pvcdoTable, densityTable,
                        rockTable, viscrefTable, watdentTable, plymwinjTables,
                        skprwatTables, skprpolyTables, tabdims, regdims, eqldims,
                        aqudims, hasImptvd, hasEntpvd, hasEqlnum, jfunc, rtemp);
}

template<class Scalar>
void unpack(Tabulated1DFunction<Scalar>& data, std::vector<char>& buffer,
            int& position, Dune::MPIHelper::MPICommunicator comm)
{
    std::vector<Scalar> xValues, yValues;
    unpack(xValues, buffer, position, comm);
    unpack(yValues, buffer, position, comm);
    data = Tabulated1DFunction<Scalar>(xValues, yValues, false);
}

template void unpack(Tabulated1DFunction<double>& data, std::vector<char>& buffer,
                     int& position, Dune::MPIHelper::MPICommunicator comm);

template<class Scalar>
void unpack(IntervalTabulated2DFunction<Scalar>& data, std::vector<char>& buffer,
            int& position, Dune::MPIHelper::MPICommunicator comm)
{
    std::vector<Scalar> xPos, yPos;
    std::vector<std::vector<Scalar>> samples;
    bool xExtrapolate, yExtrapolate;
    unpack(xPos, buffer, position, comm);
    unpack(yPos, buffer, position, comm);
    unpack(samples, buffer, position, comm);
    unpack(xExtrapolate, buffer, position, comm);
    unpack(yExtrapolate, buffer, position, comm);
    data = IntervalTabulated2DFunction<Scalar>(xPos, yPos, samples,
                                               xExtrapolate, yExtrapolate);
}

template void unpack(IntervalTabulated2DFunction<double>& data,
                     std::vector<char>& buffer,
                     int& position, Dune::MPIHelper::MPICommunicator comm);

template<class Scalar>
void unpack(UniformXTabulated2DFunction<Scalar>& data, std::vector<char>& buffer,
            int& position, Dune::MPIHelper::MPICommunicator comm)
{
    std::vector<Scalar> xPos, yPos;
    std::vector<std::vector<typename UniformXTabulated2DFunction<Scalar>::SamplePoint>> samples;
    typename UniformXTabulated2DFunction<Scalar>::InterpolationPolicy interpolationGuide;
    unpack(xPos, buffer, position, comm);
    unpack(yPos, buffer, position, comm);
    unpack(samples, buffer, position, comm);
    unpack(interpolationGuide, buffer, position, comm);
    data = UniformXTabulated2DFunction<Scalar>(xPos, yPos, samples,
                                               interpolationGuide);
}

template void unpack(UniformXTabulated2DFunction<double>& data,
                     std::vector<char>& buffer,
                     int& position, Dune::MPIHelper::MPICommunicator comm);

template<class Scalar>
void unpack(SolventPvt<Scalar>& data, std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm)
{
    std::vector<Scalar> solventReferenceDensity;
    std::vector<typename SolventPvt<Scalar>::TabulatedOneDFunction> inverseSolventB;
    std::vector<typename SolventPvt<Scalar>::TabulatedOneDFunction> solventMu;
    std::vector<typename SolventPvt<Scalar>::TabulatedOneDFunction> inverseSolventBMu;
    unpack(solventReferenceDensity, buffer, position, comm);
    unpack(inverseSolventB, buffer, position, comm);
    unpack(solventMu, buffer, position, comm);
    unpack(inverseSolventBMu, buffer, position, comm);
    data = SolventPvt<Scalar>(solventReferenceDensity, inverseSolventB,
                              solventMu, inverseSolventBMu);
}

template void unpack(SolventPvt<double>& data,
                     std::vector<char>& buffer, int& position,
                     Dune::MPIHelper::MPICommunicator comm);

template<class Scalar, bool enableThermal>
void unpack(GasPvtMultiplexer<Scalar,enableThermal>& data,
            std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm)
{
    typename GasPvtMultiplexer<Scalar,enableThermal>::GasPvtApproach approach;
    unpack(approach, buffer, position, comm);
    using PvtApproach = GasPvtMultiplexer<Scalar,enableThermal>;
    void* pvt = nullptr;
    if (approach == PvtApproach::DryGasPvt) {
        DryGasPvt<Scalar>* realPvt = new DryGasPvt<Scalar>;
        unpack(*realPvt, buffer, position, comm);
        pvt = realPvt;
    } else if (data.gasPvtApproach() == PvtApproach::WetGasPvt) {
        WetGasPvt<Scalar>* realPvt = new WetGasPvt<Scalar>;
        unpack(*realPvt, buffer, position, comm);
        pvt = realPvt;
    } else if (data.gasPvtApproach() == PvtApproach::ThermalGasPvt) {
        GasPvtThermal<Scalar>* realPvt = new GasPvtThermal<Scalar>;
        unpack(*realPvt, buffer, position, comm);
        pvt = realPvt;
    }
    data = GasPvtMultiplexer<Scalar,enableThermal>(approach, pvt);
}

template void unpack(GasPvtMultiplexer<double,true>& data,
                     std::vector<char>& buffer, int& position,
                     Dune::MPIHelper::MPICommunicator comm);
template void unpack(GasPvtMultiplexer<double,false>& data,
                     std::vector<char>& buffer, int& position,
                     Dune::MPIHelper::MPICommunicator comm);

template<class Scalar>
void unpack(DryGasPvt<Scalar>& data, std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm)
{
    std::vector<Scalar> gasReferenceDensity;
    std::vector<typename DryGasPvt<Scalar>::TabulatedOneDFunction> inverseGasB;
    std::vector<typename DryGasPvt<Scalar>::TabulatedOneDFunction> gasMu;
    std::vector<typename DryGasPvt<Scalar>::TabulatedOneDFunction> inverseGasBMu;
    unpack(gasReferenceDensity, buffer, position, comm);
    unpack(inverseGasB, buffer, position, comm);
    unpack(gasMu, buffer, position, comm);
    unpack(inverseGasBMu, buffer, position, comm);
    data = DryGasPvt<Scalar>(gasReferenceDensity, inverseGasB,
                                gasMu, inverseGasBMu);
}

template void unpack(DryGasPvt<double>& data,
                     std::vector<char>& buffer, int& position,
                     Dune::MPIHelper::MPICommunicator comm);

template<class Scalar>
void unpack(GasPvtThermal<Scalar>& data,
            std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm)
{
    std::vector<typename GasPvtThermal<Scalar>::TabulatedOneDFunction> gasvisctCurves;
    std::vector<Scalar> gasdentRefTemp, gasdentCT1, gasdentCT2;
    std::vector<typename GasPvtThermal<Scalar>::TabulatedOneDFunction> internalEnergyCurves;
    bool enableThermalDensity, enableThermalViscosity, enableInternalEnergy;
    unpack(gasvisctCurves, buffer, position, comm);
    unpack(gasdentRefTemp, buffer, position, comm);
    unpack(gasdentCT1, buffer, position, comm);
    unpack(gasdentCT2, buffer, position, comm);
    unpack(internalEnergyCurves, buffer, position, comm);
    unpack(enableThermalDensity, buffer, position, comm);
    unpack(enableThermalViscosity, buffer, position, comm);
    unpack(enableInternalEnergy, buffer, position, comm);
    bool isothermal;
    unpack(isothermal, buffer, position, comm);
    typename GasPvtThermal<Scalar>::IsothermalPvt* pvt = nullptr;
    if (isothermal) {
        pvt = new typename GasPvtThermal<Scalar>::IsothermalPvt;
        unpack(*pvt, buffer, position, comm);
    }
    data = GasPvtThermal<Scalar>(pvt, gasvisctCurves, gasdentRefTemp,
                                 gasdentCT1, gasdentCT2,
                                 internalEnergyCurves,
                                 enableThermalDensity,
                                 enableThermalViscosity,
                                 enableInternalEnergy);
}

template void unpack(GasPvtThermal<double>& data,
                     std::vector<char>& buffer, int& position,
                     Dune::MPIHelper::MPICommunicator comm);

template<class Scalar>
void unpack(WetGasPvt<Scalar>& data, std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm)
{
    std::vector<Scalar> gasReferenceDensity, oilReferenceDensity;
    std::vector<typename WetGasPvt<Scalar>::TabulatedTwoDFunction> inverseGasB;
    std::vector<typename WetGasPvt<Scalar>::TabulatedOneDFunction> inverseSaturatedGasB;
    std::vector<typename WetGasPvt<Scalar>::TabulatedTwoDFunction> gasMu;
    std::vector<typename WetGasPvt<Scalar>::TabulatedTwoDFunction> inverseGasBMu;
    std::vector<typename WetGasPvt<Scalar>::TabulatedOneDFunction> inverseSaturatedGasBMu;
    std::vector<typename WetGasPvt<Scalar>::TabulatedOneDFunction> satOilVapFacTable;
    std::vector<typename WetGasPvt<Scalar>::TabulatedOneDFunction> saturationPressure;
    Scalar vapPar1;
    unpack(gasReferenceDensity, buffer, position, comm);
    unpack(oilReferenceDensity, buffer, position, comm);
    unpack(inverseGasB, buffer, position, comm);
    unpack(inverseSaturatedGasB, buffer, position, comm);
    unpack(gasMu, buffer, position, comm);
    unpack(inverseGasBMu, buffer, position, comm);
    unpack(inverseSaturatedGasBMu, buffer, position, comm);
    unpack(satOilVapFacTable, buffer, position, comm);
    unpack(saturationPressure, buffer, position, comm);
    unpack(vapPar1, buffer, position, comm);
    data = WetGasPvt<Scalar>(gasReferenceDensity, oilReferenceDensity, inverseGasB,
                             inverseSaturatedGasB, gasMu, inverseGasBMu,
                             inverseSaturatedGasBMu, satOilVapFacTable,
                             saturationPressure, vapPar1);
}

template void unpack(WetGasPvt<double>& data,
                     std::vector<char>& buffer, int& position,
                     Dune::MPIHelper::MPICommunicator comm);

template<class Scalar, bool enableThermal>
void unpack(OilPvtMultiplexer<Scalar,enableThermal>& data,
            std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm)
{
    typename OilPvtMultiplexer<Scalar,enableThermal>::OilPvtApproach approach;
    unpack(approach, buffer, position, comm);
    using PvtApproach = OilPvtMultiplexer<Scalar,enableThermal>;
    void* pvt = nullptr;
    if (approach == PvtApproach::ConstantCompressibilityOilPvt) {
        auto* realPvt = new ConstantCompressibilityOilPvt<Scalar>;
        unpack(*realPvt, buffer, position, comm);
        pvt = realPvt;
    } else if (approach == PvtApproach::DeadOilPvt) {
        auto* realPvt = new DeadOilPvt<Scalar>;
        unpack(*realPvt, buffer, position, comm);
        pvt = realPvt;
    } else if (approach == PvtApproach::LiveOilPvt) {
        auto* realPvt = new LiveOilPvt<Scalar>;
        unpack(*realPvt, buffer, position, comm);
        pvt = realPvt;
    } else if (approach == PvtApproach::ThermalOilPvt) {
        auto* realPvt = new OilPvtThermal<Scalar>;
        unpack(*realPvt, buffer, position, comm);
        pvt = realPvt;
    }
    data = OilPvtMultiplexer<Scalar,enableThermal>(approach, pvt);
}

template void unpack(OilPvtMultiplexer<double,true>& data,
                     std::vector<char>& buffer, int& position,
                     Dune::MPIHelper::MPICommunicator comm);
template void unpack(OilPvtMultiplexer<double,false>& data,
                     std::vector<char>& buffer, int& position,
                     Dune::MPIHelper::MPICommunicator comm);

template<class Scalar>
void unpack(ConstantCompressibilityOilPvt<Scalar>& data,
            std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm)
{
    std::vector<Scalar> oilReferenceDensity, oilReferencePressure,
                        oilReferenceFormationVolumeFactor,
                        oilCompressibility, oilViscosity, oilViscosibility;

    unpack(oilReferenceDensity, buffer, position, comm);
    unpack(oilReferencePressure, buffer, position, comm);
    unpack(oilReferenceFormationVolumeFactor, buffer, position, comm);
    unpack(oilCompressibility, buffer, position, comm);
    unpack(oilViscosity, buffer, position, comm);
    unpack(oilViscosibility, buffer, position, comm);
    data = ConstantCompressibilityOilPvt<Scalar>(oilReferenceDensity,
                                                 oilReferencePressure,
                                                 oilReferenceFormationVolumeFactor,
                                                 oilCompressibility,
                                                 oilViscosity,
                                                 oilViscosibility);
}

template void unpack(ConstantCompressibilityOilPvt<double>& data,
                     std::vector<char>& buffer, int& position,
                     Dune::MPIHelper::MPICommunicator comm);

template<class Scalar>
void unpack(DeadOilPvt<Scalar>& data,
            std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm)
{
    std::vector<Scalar> oilReferenceDensity;
    using FuncVec = std::vector<typename DeadOilPvt<Scalar>::TabulatedOneDFunction>;
    FuncVec inverseOilB, oilMu, inverseOilBMu;

    unpack(oilReferenceDensity, buffer, position, comm);
    unpack(inverseOilB, buffer, position, comm);
    unpack(oilMu, buffer, position, comm);
    unpack(inverseOilBMu, buffer, position, comm);
    data = DeadOilPvt<Scalar>(oilReferenceDensity, inverseOilB,
                              oilMu, inverseOilBMu);
}

template void unpack(DeadOilPvt<double>& data,
                     std::vector<char>& buffer, int& position,
                     Dune::MPIHelper::MPICommunicator comm);

template<class Scalar>
void unpack(LiveOilPvt<Scalar>& data,
          std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    std::vector<Scalar> gasReferenceDensity, oilReferenceDensity;
    using FuncVec = std::vector<typename LiveOilPvt<Scalar>::TabulatedOneDFunction>;
    using FuncVec2 = std::vector<typename LiveOilPvt<Scalar>::TabulatedTwoDFunction>;
    FuncVec2 inverseOilBTable, oilMuTable, inverseOilBMuTable;
    FuncVec saturatedOilMuTable, inverseSaturatedOilBTable,
            inverseSaturatedOilBMuTable,
            saturatedGasDissolutionFactorTable, saturationPressure;
    Scalar vapPar2;

    unpack(gasReferenceDensity, buffer, position, comm);
    unpack(oilReferenceDensity, buffer, position, comm);
    unpack(inverseOilBTable, buffer, position, comm);
    unpack(oilMuTable, buffer, position, comm);
    unpack(inverseOilBMuTable, buffer, position, comm);
    unpack(saturatedOilMuTable, buffer, position, comm);
    unpack(inverseSaturatedOilBTable, buffer, position, comm);
    unpack(inverseSaturatedOilBMuTable, buffer, position, comm);
    unpack(saturatedGasDissolutionFactorTable, buffer, position, comm);
    unpack(saturationPressure, buffer, position, comm);
    unpack(vapPar2, buffer, position, comm);
    data = LiveOilPvt<Scalar>(gasReferenceDensity,
                              oilReferenceDensity,
                              inverseOilBTable,
                              oilMuTable,
                              inverseOilBMuTable,
                              saturatedOilMuTable,
                              inverseSaturatedOilBTable,
                              inverseSaturatedOilBMuTable,
                              saturatedGasDissolutionFactorTable,
                              saturationPressure,
                              vapPar2);
}

template void unpack(LiveOilPvt<double>& data,
                     std::vector<char>& buffer, int& position,
                     Dune::MPIHelper::MPICommunicator comm);

template<class Scalar>
void unpack(OilPvtThermal<Scalar>& data,
            std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm)
{
    std::vector<typename OilPvtThermal<Scalar>::TabulatedOneDFunction> oilvisctCurves;
    std::vector<Scalar> oildentRefTemp, oildentCT1, oildentCT2, viscrefPress, viscrefRs, viscRef;
    std::vector<typename OilPvtThermal<Scalar>::TabulatedOneDFunction> internalEnergyCurves;
    bool enableThermalDensity, enableThermalViscosity, enableInternalEnergy;
    unpack(oilvisctCurves, buffer, position, comm);
    unpack(viscrefPress, buffer, position, comm);
    unpack(viscrefRs, buffer, position, comm);
    unpack(viscRef, buffer, position, comm);
    unpack(oildentRefTemp, buffer, position, comm);
    unpack(oildentCT1, buffer, position, comm);
    unpack(oildentCT2, buffer, position, comm);
    unpack(internalEnergyCurves, buffer, position, comm);
    unpack(enableThermalDensity, buffer, position, comm);
    unpack(enableThermalViscosity, buffer, position, comm);
    unpack(enableInternalEnergy, buffer, position, comm);
    bool isothermal;
    unpack(isothermal, buffer, position, comm);
    typename OilPvtThermal<Scalar>::IsothermalPvt* pvt = nullptr;
    if (isothermal) {
        pvt = new typename OilPvtThermal<Scalar>::IsothermalPvt;
        unpack(*pvt, buffer, position, comm);
    }
    data = OilPvtThermal<Scalar>(pvt, oilvisctCurves,
                                 viscrefPress, viscrefRs, viscRef,
                                 oildentRefTemp,
                                 oildentCT1, oildentCT2,
                                 internalEnergyCurves,
                                 enableThermalDensity,
                                 enableThermalViscosity,
                                 enableInternalEnergy);
}

template void unpack(OilPvtThermal<double>& data,
                     std::vector<char>& buffer, int& position,
                     Dune::MPIHelper::MPICommunicator comm);

template<class Scalar, bool enableThermal>
void unpack(WaterPvtMultiplexer<Scalar,enableThermal>& data,
            std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm)
{
    typename WaterPvtMultiplexer<Scalar,enableThermal>::WaterPvtApproach approach;
    unpack(approach, buffer, position, comm);
    using PvtApproach = WaterPvtMultiplexer<Scalar,enableThermal>;
    void* pvt = nullptr;
    if (approach == PvtApproach::ConstantCompressibilityWaterPvt) {
        auto* realPvt = new ConstantCompressibilityWaterPvt<Scalar>;
        unpack(*realPvt, buffer, position, comm);
        pvt = realPvt;
    } else if (data.approach() == PvtApproach::ThermalWaterPvt) {
        auto* realPvt = new WaterPvtThermal<Scalar>;
        unpack(*realPvt, buffer, position, comm);
        pvt = realPvt;
    }
    data = WaterPvtMultiplexer<Scalar,enableThermal>(approach, pvt);
}

template void unpack(WaterPvtMultiplexer<double,true>& data,
                     std::vector<char>& buffer, int& position,
                     Dune::MPIHelper::MPICommunicator comm);
template void unpack(WaterPvtMultiplexer<double,false>& data,
                     std::vector<char>& buffer, int& position,
                     Dune::MPIHelper::MPICommunicator comm);

template<class Scalar>
void unpack(ConstantCompressibilityWaterPvt<Scalar>& data,
            std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm)
{
    std::vector<Scalar> waterReferenceDensity, waterReferencePressure,
                        waterReferenceFormationVolumeFactor,
                        waterCompressibility, waterViscosity, waterViscosibility;

    unpack(waterReferenceDensity, buffer, position, comm);
    unpack(waterReferencePressure, buffer, position, comm);
    unpack(waterReferenceFormationVolumeFactor, buffer, position, comm);
    unpack(waterCompressibility, buffer, position, comm);
    unpack(waterViscosity, buffer, position, comm);
    unpack(waterViscosibility, buffer, position, comm);
    data = ConstantCompressibilityWaterPvt<Scalar>(waterReferenceDensity,
                                                   waterReferencePressure,
                                                   waterReferenceFormationVolumeFactor,
                                                   waterCompressibility,
                                                   waterViscosity,
                                                   waterViscosibility);
}

template void unpack(ConstantCompressibilityWaterPvt<double>& data,
                     std::vector<char>& buffer, int& position,
                     Dune::MPIHelper::MPICommunicator comm);

template<class Scalar>
void unpack(WaterPvtThermal<Scalar>& data,
            std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm)
{
    std::vector<Scalar> viscrefPress, watdentRefTemp, watdentCT1, watdentCT2,
                        pvtwRefPress, pvtwRefB, pvtwCompressibility,
                        pvtwViscosity, pvtwViscosibility;
    std::vector<typename WaterPvtThermal<Scalar>::TabulatedOneDFunction> watvisctCurves;
    std::vector<typename WaterPvtThermal<Scalar>::TabulatedOneDFunction> internalEnergyCurves;
    bool enableThermalDensity, enableThermalViscosity, enableInternalEnergy;
    unpack(viscrefPress, buffer, position, comm);
    unpack(watdentRefTemp, buffer, position, comm);
    unpack(watdentCT1, buffer, position, comm);
    unpack(watdentCT2, buffer, position, comm);
    unpack(pvtwRefPress, buffer, position, comm);
    unpack(pvtwRefB, buffer, position, comm);
    unpack(pvtwCompressibility, buffer, position, comm);
    unpack(pvtwViscosity, buffer, position, comm);
    unpack(pvtwViscosibility, buffer, position, comm);
    unpack(watvisctCurves, buffer, position, comm);
    unpack(internalEnergyCurves, buffer, position, comm);
    unpack(enableThermalDensity, buffer, position, comm);
    unpack(enableThermalViscosity, buffer, position, comm);
    unpack(enableInternalEnergy, buffer, position, comm);
    bool isothermal;
    unpack(isothermal, buffer, position, comm);
    typename WaterPvtThermal<Scalar>::IsothermalPvt* pvt = nullptr;
    if (isothermal) {
        pvt = new typename WaterPvtThermal<Scalar>::IsothermalPvt;
        unpack(*pvt, buffer, position, comm);
    }
    data = WaterPvtThermal<Scalar>(pvt, viscrefPress, watdentRefTemp,
                                   watdentCT1, watdentCT2,
                                   pvtwRefPress, pvtwRefB,
                                   pvtwCompressibility,
                                   pvtwViscosity,
                                   pvtwViscosibility,
                                   watvisctCurves,
                                   internalEnergyCurves,
                                   enableThermalDensity,
                                   enableThermalViscosity,
                                   enableInternalEnergy);
}

template void unpack(WaterPvtThermal<double>& data,
                     std::vector<char>& buffer, int& position,
                     Dune::MPIHelper::MPICommunicator comm);

void unpack(OilVaporizationProperties& data,
          std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    OilVaporizationProperties::OilVaporization type;
    std::vector<double> vap1, vap2, maxDRSDT, maxDRVDT;
    std::vector<bool> maxDRSDT_allCells;
    unpack(type, buffer, position, comm);
    unpack(vap1, buffer, position, comm);
    unpack(vap2, buffer, position, comm);
    unpack(maxDRSDT, buffer, position, comm);
    unpack(maxDRSDT_allCells, buffer, position, comm);
    unpack(maxDRVDT, buffer, position, comm);
    data = OilVaporizationProperties(type, vap1, vap2, maxDRSDT,
                                     maxDRSDT_allCells, maxDRVDT);
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
