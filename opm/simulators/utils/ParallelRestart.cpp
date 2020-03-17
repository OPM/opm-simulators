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
#include <opm/common/OpmLog/Location.hpp>
#include <opm/parser/eclipse/EclipseState/Grid/FaceDir.hpp>
#include <opm/parser/eclipse/EclipseState/IOConfig/RestartConfig.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/Action/ActionAST.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/Action/Actions.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/Action/ActionX.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/Action/ASTNode.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/Action/Condition.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/Events.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/Group/GuideRateConfig.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/MessageLimits.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/MSW/icd.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/MSW/SpiralICD.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/MSW/Valve.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/OilVaporizationProperties.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/RFTConfig.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/ScheduleTypes.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/TimeMap.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/Tuning.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/UDQ/UDQFunction.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/UDQ/UDQFunctionTable.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/VFPInjTable.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/VFPProdTable.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/Well/Connection.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/Well/WellConnections.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/Well/WellFoamProperties.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/Well/WellPolymerProperties.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/Well/WellTracerProperties.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/Well/WList.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/Well/WListManager.hpp>
#include <opm/parser/eclipse/EclipseState/Tables/Aqudims.hpp>
#include <opm/parser/eclipse/EclipseState/Tables/FlatTable.hpp>
#include <opm/parser/eclipse/EclipseState/Tables/TableManager.hpp>
#include <opm/parser/eclipse/EclipseState/Util/IOrderSet.hpp>
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

namespace
{
template<class Type>
std::pair<std::vector<Type>, std::vector<int>>
splitDynState(const Opm::DynamicState<Type>& state)
{
    std::vector<Type> unique;
    for (const auto& w : state.data()) {
        if (std::find(unique.begin(), unique.end(), w) == unique.end())
            unique.push_back(w);
    }
    std::vector<int> idxVec;
    idxVec.reserve(state.data().size()+1);
    for (const auto& w : state.data()) {
        auto uIt = std::find(unique.begin(), unique.end(), w);
        idxVec.push_back(uIt-unique.begin());
    }
    idxVec.push_back(state.initialRange());

    return std::make_pair(unique, idxVec);
}

template<class Type>
void reconstructDynState(const std::vector<Type>& unique,
                         const std::vector<int>& idxVec,
                         Opm::DynamicState<Type>& result)
{
    std::vector<Type> ptrData;
    for (size_t i = 0; i < idxVec.size()-1; ++i) {
        ptrData.push_back(unique[idxVec[i]]);
    }
    result = Opm::DynamicState<Type>(ptrData, idxVec.back());
}

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
        size += packSize(entry, comm);

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

template<class T, class H, class KE, class A>
std::size_t packSize(const std::unordered_set<T,H,KE,A>& data,
                     Dune::MPIHelper::MPICommunicator comm)
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
                     Dune::MPIHelper::MPICommunicator comm)
{
    std::size_t totalSize = packSize(data.size(), comm);
    for (const auto& entry : data)
    {
        totalSize += packSize(entry, comm);
    }
    return totalSize;
}

template<class Key, class Value>
std::size_t packSize(const OrderedMap<Key,Value>& data, Dune::MPIHelper::MPICommunicator comm)
{
  return packSize(data.getIndex(), comm) + packSize(data.getStorage(), comm);
}

template<class T>
std::size_t packSize(const DynamicState<T>& data, Dune::MPIHelper::MPICommunicator comm)
{

    auto split = splitDynState(data);
    return packSize(split.first, comm) + packSize(split.second, comm);
}

template<class T>
std::size_t packSize(const DynamicVector<T>& data, Dune::MPIHelper::MPICommunicator comm)
{
    return packSize(data.data(), comm);
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

template<class T, std::size_t N>
std::size_t packSize(const std::array<T,N>& data, Dune::MPIHelper::MPICommunicator comm)
{
    return N*packSize(data[0], comm);
}

HANDLE_AS_POD(Aqudims)
HANDLE_AS_POD(data::Connection)
HANDLE_AS_POD(data::CurrentControl)
HANDLE_AS_POD(data::Rates)
HANDLE_AS_POD(data::Segment)
HANDLE_AS_POD(MLimits)
HANDLE_AS_POD(PlyvmhRecord)
HANDLE_AS_POD(Tabdims)
HANDLE_AS_POD(TimeStampUTC::YMD)
HANDLE_AS_POD(Tuning)
HANDLE_AS_POD(WellBrineProperties)
HANDLE_AS_POD(WellFoamProperties)

std::size_t packSize(const data::Well& data, Dune::MPIHelper::MPICommunicator comm)
{
    std::size_t size = packSize(data.rates, comm);
    size += packSize(data.bhp, comm) + packSize(data.thp, comm);
    size += packSize(data.temperature, comm);
    size += packSize(data.control, comm);
    size += packSize(data.connections, comm);
    size += packSize(data.segments, comm);
    size += packSize(data.current_control, comm);
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

std::size_t packSize(const WellType& data, Dune::MPIHelper::MPICommunicator comm)
{
    return packSize(data.producer(), comm) +
           packSize(data.preferred_phase(), comm);
}

std::size_t packSize(const TimeMap& data, Dune::MPIHelper::MPICommunicator comm)
{
    return packSize(data.timeList(), comm);
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

namespace {

struct SplitSimpleTables {
    size_t plyshMax = 0;
    size_t rockMax = 0;
    std::map<size_t, std::shared_ptr<PlyshlogTable>> plyshMap;
    std::map<size_t, std::shared_ptr<RocktabTable>> rockMap;
};

SplitSimpleTables
splitSimpleTable(std::map<std::string, TableContainer>& simpleTables)
{
    SplitSimpleTables result;

    // PlyshlogTable need special treatment
    auto it = simpleTables.find("PLYSHLOG");
    if (it != simpleTables.end()) {
        result.plyshMax = it->second.max();
        for (const auto& mapIt : it->second.tables()) {
            auto ptr = std::static_pointer_cast<PlyshlogTable>(mapIt.second);
            result.plyshMap.insert(std::make_pair(mapIt.first, ptr));
        }
        simpleTables.erase(it);
    }

    // RocktabTable need special treatment
    it = simpleTables.find("ROCKMAP");
    if (it != simpleTables.end()) {
        result.rockMax = it->second.max();
        for (const auto& mapIt : it->second.tables()) {
            auto ptr = std::static_pointer_cast<RocktabTable>(mapIt.second);
            result.rockMap.insert(std::make_pair(mapIt.first,  ptr));
        }
        simpleTables.erase(it);
    }

    return result;
}

}

std::size_t packSize(const TableManager& data, Dune::MPIHelper::MPICommunicator comm)
{
    auto simpleTables = data.getSimpleTables();
    auto splitTab = splitSimpleTable(simpleTables);

    return packSize(simpleTables, comm) +
           packSize(splitTab.plyshMax, comm) +
           packSize(splitTab.plyshMap, comm) +
           packSize(splitTab.rockMax, comm) +
           packSize(splitTab.rockMap, comm) +
           packSize(data.getPvtgTables(), comm) +
           packSize(data.getPvtoTables(), comm) +
           packSize(data.getRock2dTables(), comm) +
           packSize(data.getRock2dtrTables(), comm) +
           packSize(data.getPvtwTable(), comm) +
           packSize(data.getPvcdoTable(), comm) +
           packSize(data.getDensityTable(), comm) +
           packSize(data.getPlyvmhTable(), comm) +
           packSize(data.getRockTable(), comm) +
           packSize(data.getPlmixparTable(), comm) +
           packSize(data.getShrateTable(), comm) +
           packSize(data.getStone1exTable(), comm) +
           packSize(data.getTlmixparTable(), comm) +
           packSize(data.getViscrefTable(), comm) +
           packSize(data.getWatdentTable(), comm) +
           packSize(data.getPvtwSaltTables(), comm) +
           packSize(data.getBrineDensityTables(), comm) +
           packSize(data.getSolventDensityTables(), comm) +
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
           packSize(data.useShrate(), comm) +
           packSize(data.useJFunc(), comm) +
          (data.useJFunc() ? packSize(data.getJFunc(), comm) : 0) +
           packSize(data.OilDenT(), comm) +
           packSize(data.GasDenT(), comm) +
           packSize(data.WatDenT(), comm) +
           packSize(data.stCond(), comm) +
           packSize(data.gas_comp_index(), comm) +
           packSize(data.rtemp(), comm);
}

template
std::size_t packSize(const std::map<Phase,Group::GroupInjectionProperties>& data,
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

std::size_t packSize(const Events& data,
                     Dune::MPIHelper::MPICommunicator comm)
{
    return packSize(data.events(), comm);
}

std::size_t packSize(const MessageLimits& data,
                     Dune::MPIHelper::MPICommunicator comm)
{
    return packSize(data.getLimits(), comm);
}

std::size_t packSize(const VFPInjTable& data,
                     Dune::MPIHelper::MPICommunicator comm)
{
    return packSize(data.getTableNum(), comm) +
           packSize(data.getDatumDepth(), comm) +
           packSize(data.getFloType(), comm) +
           packSize(data.getFloAxis(), comm) +
           packSize(data.getTHPAxis(), comm) +
           packSize(data.getTable(), comm);
}

std::size_t packSize(const VFPProdTable& data,
                     Dune::MPIHelper::MPICommunicator comm)
{
    return packSize(data.getTableNum(), comm) +
           packSize(data.getDatumDepth(), comm) +
           packSize(data.getFloType(), comm) +
           packSize(data.getWFRType(), comm) +
           packSize(data.getGFRType(), comm) +
           packSize(data.getALQType(), comm) +
           packSize(data.getFloAxis(), comm) +
           packSize(data.getTHPAxis(), comm) +
           packSize(data.getWFRAxis(), comm) +
           packSize(data.getGFRAxis(), comm) +
           packSize(data.getALQAxis(), comm) +
           packSize(data.getTable(), comm);
}

std::size_t packSize(const WellTestConfig::WTESTWell& data,
                     Dune::MPIHelper::MPICommunicator comm)
{
    return packSize(data.name, comm) +
           packSize(data.shut_reason, comm) +
           packSize(data.test_interval, comm) +
           packSize(data.num_test, comm) +
           packSize(data.startup_time, comm) +
           packSize(data.begin_report_step, comm);
}

std::size_t packSize(const WellTestConfig& data,
                     Dune::MPIHelper::MPICommunicator comm)
{
    return packSize(data.getWells(), comm);
}

std::size_t packSize(const WellTracerProperties& data,
                     Dune::MPIHelper::MPICommunicator comm)
{
    return packSize(data.getConcentrations(), comm);
}

std::size_t packSize(const UDAValue& data,
                     Dune::MPIHelper::MPICommunicator comm)
{
    return packSize(data.get_dim(), comm) +
           packSize(data.is<double>(), comm) +
           (data.is<double>() ? packSize(data.get<double>(), comm) :
                                packSize(data.get<std::string>(), comm));
}

std::size_t packSize(const Connection& data,
                     Dune::MPIHelper::MPICommunicator comm)
{
    return packSize(data.dir(), comm) +
           packSize(data.depth(), comm) +
           packSize(data.state(), comm) +
           packSize(data.satTableId(), comm) +
           packSize(data.complnum(), comm) +
           packSize(data.CF(), comm) +
           packSize(data.Kh(), comm) +
           packSize(data.rw(), comm) +
           packSize(data.r0(), comm) +
           packSize(data.skinFactor(), comm) +
           packSize(data.getI(), comm) +
           packSize(data.getJ(), comm) +
           packSize(data.getK(), comm) +
           packSize(data.kind(), comm) +
           packSize(data.getSeqIndex(), comm) +
           packSize(data.getSegDistStart(), comm) +
           packSize(data.getSegDistEnd(), comm) +
           packSize(data.getDefaultSatTabId(), comm) +
           packSize(data.getCompSegSeqIndex(), comm) +
           packSize(data.segment(), comm) +
           packSize(data.wellPi(), comm);
}

std::size_t packSize(const Well::WellInjectionProperties& data,
                     Dune::MPIHelper::MPICommunicator comm)
{
    return packSize(data.name, comm) +
           packSize(data.surfaceInjectionRate, comm) +
           packSize(data.reservoirInjectionRate, comm) +
           packSize(data.BHPTarget, comm) +
           packSize(data.THPTarget, comm) +
           packSize(data.bhp_hist_limit, comm) +
           packSize(data.thp_hist_limit, comm) +
           packSize(data.temperature, comm) +
           packSize(data.BHPH, comm) +
           packSize(data.THPH, comm) +
           packSize(data.VFPTableNumber, comm) +
           packSize(data.predictionMode, comm) +
           packSize(data.injectionControls, comm) +
           packSize(data.injectorType, comm) +
           packSize(data.controlMode, comm);
}

std::size_t packSize(const WellEconProductionLimits& data,
                     Dune::MPIHelper::MPICommunicator comm)
{
    return packSize(data.minOilRate(), comm) +
           packSize(data.minGasRate(), comm) +
           packSize(data.maxWaterCut(), comm) +
           packSize(data.maxGasOilRatio(), comm) +
           packSize(data.maxWaterGasRatio(), comm) +
           packSize(data.workover(), comm) +
           packSize(data.endRun(), comm) +
           packSize(data.followonWell(), comm) +
           packSize(data.quantityLimit(), comm) +
           packSize(data.maxSecondaryMaxWaterCut(), comm) +
           packSize(data.workoverSecondary(), comm) +
           packSize(data.maxGasLiquidRatio(), comm) +
           packSize(data.minLiquidRate(), comm) +
           packSize(data.maxTemperature(), comm) +
           packSize(data.minReservoirFluidRate(), comm);
}

std::size_t packSize(const WellConnections& data,
                     Dune::MPIHelper::MPICommunicator comm)
{
    return packSize(data.getHeadI(), comm) +
           packSize(data.getHeadJ(), comm) +
           packSize(data.getNumRemoved(), comm) +
           packSize(data.getConnections(), comm);
}

std::size_t packSize(const Well::WellProductionProperties& data,
                     Dune::MPIHelper::MPICommunicator comm)
{
    return packSize(data.name, comm) +
           packSize(data.OilRate, comm) +
           packSize(data.WaterRate, comm) +
           packSize(data.GasRate, comm) +
           packSize(data.LiquidRate, comm) +
           packSize(data.ResVRate, comm) +
           packSize(data.BHPTarget, comm) +
           packSize(data.THPTarget, comm) +
           packSize(data.bhp_hist_limit, comm) +
           packSize(data.thp_hist_limit, comm) +
           packSize(data.BHPH, comm) +
           packSize(data.THPH, comm) +
           packSize(data.VFPTableNumber, comm) +
           packSize(data.ALQValue, comm) +
           packSize(data.predictionMode, comm) +
           packSize(data.controlMode, comm) +
           packSize(data.whistctl_cmode, comm) +
           packSize(data.getNumProductionControls(), comm);
}

std::size_t packSize(const SpiralICD& data,
                     Dune::MPIHelper::MPICommunicator comm)
{
    return packSize(data.strength(), comm) +
           packSize(data.length(), comm) +
           packSize(data.densityCalibration(), comm) +
           packSize(data.viscosityCalibration(), comm) +
           packSize(data.criticalValue(), comm) +
           packSize(data.widthTransitionRegion(), comm) +
           packSize(data.maxViscosityRatio(), comm) +
           packSize(data.methodFlowScaling(), comm) +
           packSize(data.maxAbsoluteRate(), comm) +
           packSize(data.status(), comm) +
           packSize(data.scalingFactor(), comm);
}

std::size_t packSize(const Valve& data,
                     Dune::MPIHelper::MPICommunicator comm)
{
    return packSize(data.conFlowCoefficient(), comm) +
           packSize(data.conCrossArea(), comm) +
           packSize(data.conMaxCrossArea(), comm) +
           packSize(data.pipeAdditionalLength(), comm) +
           packSize(data.pipeDiameter(), comm) +
           packSize(data.pipeRoughness(), comm) +
           packSize(data.pipeCrossArea(), comm) +
           packSize(data.status(), comm);
}

std::size_t packSize(const Segment& data,
                     Dune::MPIHelper::MPICommunicator comm)
{
    return packSize(data.segmentNumber(), comm) +
           packSize(data.branchNumber(), comm) +
           packSize(data.outletSegment(), comm) +
           packSize(data.inletSegments(), comm) +
           packSize(data.totalLength(), comm) +
           packSize(data.depth(), comm) +
           packSize(data.internalDiameter(), comm) +
           packSize(data.roughness(), comm) +
           packSize(data.crossArea(), comm) +
           packSize(data.volume(), comm) +
           packSize(data.dataReady(), comm) +
           packSize(data.segmentType(), comm) +
           packSize(data.spiralICD(), comm) +
           packSize(data.getValve(), comm);
}

template<class T>
std::size_t packSize(const std::shared_ptr<T>& data,
                     Dune::MPIHelper::MPICommunicator comm)
{
    std::size_t size = packSize(bool(), comm);
    if (data)
         size += packSize(*data, comm);

    return size;
}

template<class T>
std::size_t packSize(const std::unique_ptr<T>& data,
                     Dune::MPIHelper::MPICommunicator comm)
{
    std::size_t size = packSize(bool(), comm);
    if (data)
         size += packSize(*data, comm);

    return size;
}

std::size_t packSize(const Dimension& data,
                     Dune::MPIHelper::MPICommunicator comm)
{
    return packSize(data.getSIScalingRaw(), comm) +
           packSize(data.getSIOffset(), comm);
}

std::size_t packSize(const UnitSystem& data,
                     Dune::MPIHelper::MPICommunicator comm)
{
    return packSize(data.getName(), comm) +
           packSize(data.getType(), comm) +
           packSize(data.getDimensions(), comm) +
           packSize(data.use_count(), comm);
}

std::size_t packSize(const WellSegments& data,
                     Dune::MPIHelper::MPICommunicator comm)
{
    return packSize(data.compPressureDrop(), comm) +
           packSize(data.segments(), comm);
}

std::size_t packSize(const Well& data,
                     Dune::MPIHelper::MPICommunicator comm)
{
    std::size_t size = packSize(data.name(), comm) +
                       packSize(data.groupName(), comm) +
                       packSize(data.firstTimeStep(), comm) +
                       packSize(data.seqIndex(), comm) +
                       packSize(data.getHeadI(), comm) +
                       packSize(data.getHeadJ(), comm) +
                       packSize(data.getRefDepth(), comm) +
                       packSize(data.wellType(), comm) +
                       packSize(data.getWellConnectionOrdering(), comm) +
                       packSize(data.units(), comm) +
                       packSize(data.udqUndefined(), comm) +
                       packSize(data.getStatus(), comm) +
                       packSize(data.getDrainageRadius(), comm) +
                       packSize(data.getAllowCrossFlow(), comm) +
                       packSize(data.getAutomaticShutIn(), comm) +
                       packSize(data.wellGuideRate(), comm) +
                       packSize(data.getEfficiencyFactor(), comm) +
                       packSize(data.getSolventFraction(), comm) +
                       packSize(data.predictionMode(), comm) +
                       packSize(data.getEconLimits(), comm) +
                       packSize(data.getFoamProperties(), comm) +
                       packSize(data.getPolymerProperties(), comm) +
                       packSize(data.getBrineProperties(), comm) +
                       packSize(data.getTracerProperties(), comm) +
                       packSize(data.getConnections(), comm) +
                       packSize(data.getProductionProperties(), comm) +
                       packSize(data.getInjectionProperties(), comm) +
                       packSize(data.hasSegments(), comm);
    if (data.hasSegments())
        size += packSize(data.getSegments(), comm);

    return size;
}

template<class T>
std::size_t packSize(const IOrderSet<T>& data,
                     Dune::MPIHelper::MPICommunicator comm)
{
    return packSize(data.index(), comm) +
           packSize(data.data(), comm);
}

std::size_t packSize(const Group::GroupInjectionProperties& data,
                     Dune::MPIHelper::MPICommunicator comm)
{
    return packSize(data.phase, comm) +
           packSize(data.cmode, comm) +
           packSize(data.surface_max_rate, comm) +
           packSize(data.resv_max_rate, comm) +
           packSize(data.target_reinj_fraction, comm) +
           packSize(data.target_void_fraction, comm) +
           packSize(data.reinj_group, comm) +
           packSize(data.voidage_group, comm) +
           packSize(data.injection_controls, comm);
}

std::size_t packSize(const Group::GroupProductionProperties& data,
                     Dune::MPIHelper::MPICommunicator comm)
{
    return packSize(data.cmode, comm) +
           packSize(data.exceed_action, comm) +
           packSize(data.oil_target, comm) +
           packSize(data.water_target, comm) +
           packSize(data.gas_target, comm) +
           packSize(data.liquid_target, comm) +
           packSize(data.guide_rate, comm) +
           packSize(data.guide_rate_def, comm) +
           packSize(data.resv_target, comm) +
           packSize(data.production_controls, comm);
}

std::size_t packSize(const Group& data,
                     Dune::MPIHelper::MPICommunicator comm)
{
    return packSize(data.name(), comm) +
           packSize(data.insert_index(), comm) +
           packSize(data.initStep(), comm) +
           packSize(data.udqUndefined(), comm) +
           packSize(data.units(), comm) +
           packSize(data.type(), comm) +
           packSize(data.getGroupEfficiencyFactor(), comm) +
           packSize(data.getTransferGroupEfficiencyFactor(), comm) +
           packSize(data.isAvailableForGroupControl(), comm) +
           packSize(data.getGroupNetVFPTable(), comm) +
           packSize(data.parent(), comm) +
           packSize(data.iwells(), comm) +
           packSize(data.igroups(), comm) +
           packSize(data.injectionProperties(), comm) +
           packSize(data.productionProperties(), comm);
}

std::size_t packSize(const WList& data,
                     Dune::MPIHelper::MPICommunicator comm)
{
    return packSize(data.wellList(), comm);
}

std::size_t packSize(const WListManager& data,
                     Dune::MPIHelper::MPICommunicator comm)
{
    return packSize(data.lists(), comm);
}

std::size_t packSize(const GuideRateModel& data,
                     Dune::MPIHelper::MPICommunicator comm)
{
    return packSize(data.timeInterval(), comm) +
           packSize(data.target(), comm) +
           packSize(data.coefs(), comm) +
           packSize(data.allow_increase(), comm) +
           packSize(data.damping_factor(), comm) +
           packSize(data.free_gas(), comm) +
           packSize(data.defaultModel(), comm) +
           packSize(data.udaCoefs(), comm);
}

std::size_t packSize(const GuideRateConfig& data,
                     Dune::MPIHelper::MPICommunicator comm)
{
    return packSize(data.getModel(), comm) +
           packSize(data.getWells(), comm) +
           packSize(data.getGroups(), comm);
}

std::size_t packSize(const GConSale::GCONSALEGroup& data,
                     Dune::MPIHelper::MPICommunicator comm)
{
    return packSize(data.sales_target, comm) +
           packSize(data.max_sales_rate, comm) +
           packSize(data.min_sales_rate, comm) +
           packSize(data.max_proc, comm) +
           packSize(data.udq_undefined, comm) +
           packSize(data.unit_system, comm);
}

std::size_t packSize(const GConSale& data,
                     Dune::MPIHelper::MPICommunicator comm)
{
    return packSize(data.getGroups(), comm);
}

std::size_t packSize(const GConSump::GCONSUMPGroup& data,
                     Dune::MPIHelper::MPICommunicator comm)
{
    return packSize(data.consumption_rate, comm) +
           packSize(data.import_rate, comm) +
           packSize(data.network_node, comm) +
           packSize(data.udq_undefined, comm) +
           packSize(data.unit_system, comm);
}

std::size_t packSize(const GConSump& data,
                     Dune::MPIHelper::MPICommunicator comm)
{
    return packSize(data.getGroups(), comm);
}

std::size_t packSize(const RFTConfig& data,
                     Dune::MPIHelper::MPICommunicator comm)
{
    return packSize(data.timeMap(), comm) +
           packSize(data.firstRFTOutput(), comm) +
           packSize(data.wellOpenRftTime(), comm) +
           packSize(data.wellOpenRftName(), comm) +
           packSize(data.wellOpen(), comm) +
           packSize(data.rftConfig(), comm) +
           packSize(data.pltConfig(), comm);
}

std::size_t packSize(const DeckItem& data,
                     Dune::MPIHelper::MPICommunicator comm)
{
    return packSize(data.dVal(), comm) +
           packSize(data.iVal(), comm) +
           packSize(data.sVal(), comm) +
           packSize(data.uVal(), comm) +
           packSize(data.getType(), comm) +
           packSize(data.name(), comm) +
           packSize(data.valueStatus(), comm) +
           packSize(data.rawData(), comm) +
           packSize(data.activeDimensions(), comm) +
           packSize(data.defaultDimensions(), comm);
}

std::size_t packSize(const DeckRecord& data,
                     Dune::MPIHelper::MPICommunicator comm)
{
    return packSize(data.getItems(), comm);
}

std::size_t packSize(const Location& data,
                     Dune::MPIHelper::MPICommunicator comm)
{
    return packSize(data.filename, comm) +
           packSize(data.lineno, comm);
}

std::size_t packSize(const DeckKeyword& data,
                     Dune::MPIHelper::MPICommunicator comm)
{
    return packSize(data.name(), comm) +
           packSize(data.location(), comm) +
           packSize(data.records(), comm) +
           packSize(data.isDataKeyword(), comm) +
           packSize(data.isSlashTerminated(), comm);
}


std::size_t packSize(const Deck& data,
                     Dune::MPIHelper::MPICommunicator comm)
{
    return packSize(data.keywords(), comm) +
           packSize(data.getDefaultUnitSystem(), comm) +
           packSize(data.activeUnitSystem(), comm) +
           packSize(data.getDataFile(), comm) +
           packSize(data.getInputPath(), comm) +
           packSize(data.unitSystemAccessCount(), comm);
}

std::size_t packSize(const Action::ASTNode& data,
                     Dune::MPIHelper::MPICommunicator comm)
{
    return packSize(data.type, comm) +
           packSize(data.func_type, comm) +
           packSize(data.func, comm) +
           packSize(data.argList(), comm) +
           packSize(data.getNumber(), comm) +
           packSize(data.childrens(), comm);
}

std::size_t packSize(const Action::AST& data,
                     Dune::MPIHelper::MPICommunicator comm)
{
    return packSize(data.getCondition(), comm);
}

std::size_t packSize(const Action::Quantity& data,
                     Dune::MPIHelper::MPICommunicator comm)
{
    return packSize(data.quantity, comm) +
           packSize(data.args, comm);
}

std::size_t packSize(const Action::Condition& data,
                     Dune::MPIHelper::MPICommunicator comm)
{
    return packSize(data.lhs, comm) +
           packSize(data.rhs, comm) +
           packSize(data.logic, comm) +
           packSize(data.cmp, comm) +
           packSize(data.cmp_string, comm);
}

std::size_t packSize(const Action::ActionX& data,
                     Dune::MPIHelper::MPICommunicator comm)
{
    return packSize(data.name(), comm) +
           packSize(data.max_run(), comm) +
           packSize(data.min_wait(), comm) +
           packSize(data.start_time(), comm) +
           packSize(data.getKeywords(), comm) +
           packSize(data.getCondition(), comm) +
           packSize(data.conditions(), comm) +
           packSize(data.getRunCount(), comm) +
           packSize(data.getLastRun(), comm);
}

std::size_t packSize(const Action::Actions& data,
                     Dune::MPIHelper::MPICommunicator comm)
{
    return packSize(data.getActions(), comm);
}

std::size_t packSize(const RestartSchedule& data,
                     Dune::MPIHelper::MPICommunicator comm)
{
    return packSize(data.timestep, comm) +
           packSize(data.basic, comm) +
           packSize(data.frequency, comm) +
           packSize(data.rptsched_restart_set, comm) +
           packSize(data.rptsched_restart, comm);
}

std::size_t packSize(const TimeStampUTC& data,
                     Dune::MPIHelper::MPICommunicator comm)
{
    return packSize(data.ymd(), comm) +
           packSize(data.hour(), comm) +
           packSize(data.minutes(), comm) +
           packSize(data.seconds(), comm) +
           packSize(data.microseconds(), comm);
}

std::size_t packSize(const WellPolymerProperties& data,
                     Dune::MPIHelper::MPICommunicator comm)
{
    return packSize(data.m_polymerConcentration, comm) +
           packSize(data.m_saltConcentration, comm) +
           packSize(data.m_plymwinjtable, comm) +
           packSize(data.m_skprwattable, comm) +
           packSize(data.m_skprpolytable, comm);
}

std::size_t packSize(const Well::WellGuideRate& data,
                     Dune::MPIHelper::MPICommunicator comm)
{
    return packSize(data.available, comm) +
           packSize(data.guide_rate, comm) +
           packSize(data.guide_phase, comm) +
           packSize(data.scale_factor, comm);
}

std::size_t packSize(const GuideRateConfig::WellTarget& data,
                     Dune::MPIHelper::MPICommunicator comm)
{
    return packSize(data.guide_rate, comm) +
           packSize(data.target, comm) +
           packSize(data.scaling_factor, comm);
}

std::size_t packSize(const GuideRateConfig::GroupTarget& data,
                     Dune::MPIHelper::MPICommunicator comm)
{
    return packSize(data.guide_rate, comm) +
           packSize(data.target, comm);
}

std::size_t packSize(const PlyvmhTable& data, Dune::MPIHelper::MPICommunicator comm)
{
    return packSize(static_cast<const std::vector<PlyvmhRecord>&>(data), comm);
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

template<class K, class C, class A>
void pack(const std::set<K,C,A>& data,
          std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
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
          Dune::MPIHelper::MPICommunicator comm)
{
    pack(data.size(), buffer, position, comm);

    for (const auto& entry : data)
    {
        pack(entry, buffer, position, comm);
    }
}

template<class T, size_t N>
void pack(const std::array<T,N>& data, std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    for (const T& entry : data)
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
    auto split = splitDynState(data);
    pack(split.first, buffer, position, comm);
    pack(split.second, buffer, position, comm);
}

template<class T>
void pack(const DynamicVector<T>& data, std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    pack(data.data(), buffer, position, comm);
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


template void pack(const std::map<Phase, Group::GroupInjectionProperties>& data,
          std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm);



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
    pack(data.current_control, buffer, position, comm);
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

void pack(const WellType& data, std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm) {
    pack(data.producer(), buffer, position, comm);
    pack(data.preferred_phase(), buffer, position, comm);
}

void pack(const TimeMap& data, std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    pack(data.timeList(), buffer, position, comm);
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

void pack(const TableManager& data, std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    auto simpleTables = data.getSimpleTables();
    auto splitTab = splitSimpleTable(simpleTables);

    pack(simpleTables, buffer, position, comm);
    pack(splitTab.plyshMax, buffer, position, comm);
    pack(splitTab.plyshMap, buffer, position, comm);
    pack(splitTab.rockMax, buffer, position, comm);
    pack(splitTab.rockMap, buffer, position, comm);
    pack(data.getPvtgTables(), buffer, position, comm);
    pack(data.getPvtoTables(), buffer, position, comm);
    pack(data.getRock2dTables(), buffer, position, comm);
    pack(data.getRock2dtrTables(), buffer, position, comm);
    pack(data.getPvtwTable(), buffer, position, comm);
    pack(data.getPvcdoTable(), buffer, position, comm);
    pack(data.getDensityTable(), buffer, position, comm);
    pack(data.getPlyvmhTable(), buffer, position, comm);
    pack(data.getRockTable(), buffer, position, comm);
    pack(data.getPlmixparTable(), buffer, position, comm);
    pack(data.getShrateTable(), buffer, position, comm);
    pack(data.getStone1exTable(), buffer, position, comm);
    pack(data.getTlmixparTable(), buffer, position, comm);
    pack(data.getViscrefTable(), buffer, position, comm);
    pack(data.getWatdentTable(), buffer, position, comm);
    pack(data.getPvtwSaltTables(), buffer, position, comm);
    pack(data.getBrineDensityTables(), buffer, position, comm);
    pack(data.getSolventDensityTables(), buffer, position, comm);
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
    pack(data.useShrate(), buffer, position, comm);
    pack(data.useJFunc(), buffer, position, comm);
    if (data.useJFunc())
        pack(data.getJFunc(), buffer, position, comm);
    pack(data.OilDenT(), buffer, position, comm);
    pack(data.GasDenT(), buffer, position, comm);
    pack(data.WatDenT(), buffer, position, comm);
    pack(data.stCond(), buffer, position, comm);
    pack(data.gas_comp_index(), buffer, position, comm);
    pack(data.rtemp(), buffer, position, comm);
}

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

void pack(const Events& data,
          std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    pack(data.events(), buffer, position, comm);
}

void pack(const MessageLimits& data,
          std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    pack(data.getLimits(), buffer, position, comm);
}
void pack(const VFPInjTable& data,
          std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    pack(data.getTableNum(), buffer, position, comm);
    pack(data.getDatumDepth(), buffer, position, comm);
    pack(data.getFloType(), buffer, position, comm);
    pack(data.getFloAxis(), buffer, position, comm);
    pack(data.getTHPAxis(), buffer, position, comm);
    pack(data.getTable(), buffer, position, comm);
}

void pack(const VFPProdTable& data,
          std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    pack(data.getTableNum(), buffer, position, comm);
    pack(data.getDatumDepth(), buffer, position, comm);
    pack(data.getFloType(), buffer, position, comm);
    pack(data.getWFRType(), buffer, position, comm);
    pack(data.getGFRType(), buffer, position, comm);
    pack(data.getALQType(), buffer, position, comm);
    pack(data.getFloAxis(), buffer, position, comm);
    pack(data.getTHPAxis(), buffer, position, comm);
    pack(data.getWFRAxis(), buffer, position, comm);
    pack(data.getGFRAxis(), buffer, position, comm);
    pack(data.getALQAxis(), buffer, position, comm);
    pack(data.getTable(), buffer, position, comm);
}

void pack(const WellTestConfig::WTESTWell& data,
          std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    pack(data.name, buffer, position, comm);
    pack(data.shut_reason, buffer, position, comm);
    pack(data.test_interval, buffer, position, comm);
    pack(data.num_test, buffer, position, comm);
    pack(data.startup_time, buffer, position, comm);
    pack(data.begin_report_step, buffer, position, comm);
}

void pack(const WellTestConfig& data,
          std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    pack(data.getWells(), buffer, position, comm);
}

void pack(const WellTracerProperties& data,
          std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    pack(data.getConcentrations(), buffer, position, comm);
}

void pack(const UDAValue& data,
          std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    pack(data.get_dim(), buffer, position, comm);
    pack(data.is<double>(), buffer, position, comm);
    if (data.is<double>())
        pack(data.get<double>(), buffer, position, comm);
    else
        pack(data.get<std::string>(), buffer, position, comm);
}

void pack(const Connection& data,
          std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    pack(data.dir(), buffer, position, comm);
    pack(data.depth(), buffer, position, comm);
    pack(data.state(), buffer, position, comm);
    pack(data.satTableId(), buffer, position, comm);
    pack(data.complnum(), buffer, position, comm);
    pack(data.CF(), buffer, position, comm);
    pack(data.Kh(), buffer, position, comm);
    pack(data.rw(), buffer, position, comm);
    pack(data.r0(), buffer, position, comm);
    pack(data.skinFactor(), buffer, position, comm);
    pack(data.getI(), buffer, position, comm);
    pack(data.getJ(), buffer, position, comm);
    pack(data.getK(), buffer, position, comm);
    pack(data.kind(), buffer, position, comm);
    pack(data.getSeqIndex(), buffer, position, comm);
    pack(data.getSegDistStart(), buffer, position, comm);
    pack(data.getSegDistEnd(), buffer, position, comm);
    pack(data.getDefaultSatTabId(), buffer, position, comm);
    pack(data.getCompSegSeqIndex(), buffer, position, comm);
    pack(data.segment(), buffer, position, comm);
    pack(data.wellPi(), buffer, position, comm);
}

void pack(const Well::WellInjectionProperties& data,
          std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    pack(data.name, buffer, position, comm);
    pack(data.surfaceInjectionRate, buffer, position, comm);
    pack(data.reservoirInjectionRate, buffer, position, comm);
    pack(data.BHPTarget, buffer, position, comm);
    pack(data.THPTarget, buffer, position, comm);
    pack(data.bhp_hist_limit, buffer, position, comm);
    pack(data.thp_hist_limit, buffer, position, comm);
    pack(data.temperature, buffer, position, comm);
    pack(data.BHPH, buffer, position, comm);
    pack(data.THPH, buffer, position, comm);
    pack(data.VFPTableNumber, buffer, position, comm);
    pack(data.predictionMode, buffer, position, comm);
    pack(data.injectionControls, buffer, position, comm);
    pack(data.injectorType, buffer, position, comm);
    pack(data.controlMode, buffer, position, comm);
}

void pack(const WellEconProductionLimits& data,
          std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    pack(data.minOilRate(), buffer, position, comm);
    pack(data.minGasRate(), buffer, position, comm);
    pack(data.maxWaterCut(), buffer, position, comm);
    pack(data.maxGasOilRatio(), buffer, position, comm);
    pack(data.maxWaterGasRatio(), buffer, position, comm);
    pack(data.workover(), buffer, position, comm);
    pack(data.endRun(), buffer, position, comm);
    pack(data.followonWell(), buffer, position, comm);
    pack(data.quantityLimit(), buffer, position, comm);
    pack(data.maxSecondaryMaxWaterCut(), buffer, position, comm);
    pack(data.workoverSecondary(), buffer, position, comm);
    pack(data.maxGasLiquidRatio(), buffer, position, comm);
    pack(data.minLiquidRate(), buffer, position, comm);
    pack(data.maxTemperature(), buffer, position, comm);
    pack(data.minReservoirFluidRate(), buffer, position, comm);
}

void pack(const WellConnections& data,
          std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    pack(data.getHeadI(), buffer, position, comm);
    pack(data.getHeadJ(), buffer, position, comm);
    pack(data.getNumRemoved(), buffer, position, comm);
    pack(data.getConnections(), buffer, position, comm);
}

void pack(const Well::WellProductionProperties& data,
          std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    pack(data.name, buffer, position, comm);
    pack(data.OilRate, buffer, position, comm);
    pack(data.WaterRate, buffer, position, comm);
    pack(data.GasRate, buffer, position, comm);
    pack(data.LiquidRate, buffer, position, comm);
    pack(data.ResVRate, buffer, position, comm);
    pack(data.BHPTarget, buffer, position, comm);
    pack(data.THPTarget, buffer, position, comm);
    pack(data.bhp_hist_limit, buffer, position, comm);
    pack(data.thp_hist_limit, buffer, position, comm);
    pack(data.BHPH, buffer, position, comm);
    pack(data.THPH, buffer, position, comm);
    pack(data.VFPTableNumber, buffer, position, comm);
    pack(data.ALQValue, buffer, position, comm);
    pack(data.predictionMode, buffer, position, comm);
    pack(data.controlMode, buffer, position, comm);
    pack(data.whistctl_cmode, buffer, position, comm);
    pack(data.getNumProductionControls(), buffer, position, comm);
}

void pack(const SpiralICD& data,
          std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    pack(data.strength(), buffer, position, comm);
    pack(data.length(), buffer, position, comm);
    pack(data.densityCalibration(), buffer, position, comm);
    pack(data.viscosityCalibration(), buffer, position, comm);
    pack(data.criticalValue(), buffer, position, comm);
    pack(data.widthTransitionRegion(), buffer, position, comm);
    pack(data.maxViscosityRatio(), buffer, position, comm);
    pack(data.methodFlowScaling(), buffer, position, comm);
    pack(data.maxAbsoluteRate(), buffer, position, comm);
    pack(data.status(), buffer, position, comm);
    pack(data.scalingFactor(), buffer, position, comm);
}

void pack(const Valve& data,
          std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    pack(data.conFlowCoefficient(), buffer, position, comm);
    pack(data.conCrossArea(), buffer, position, comm);
    pack(data.conMaxCrossArea(), buffer, position, comm);
    pack(data.pipeAdditionalLength(), buffer, position, comm);
    pack(data.pipeDiameter(), buffer, position, comm);
    pack(data.pipeRoughness(), buffer, position, comm);
    pack(data.pipeCrossArea(), buffer, position, comm);
    pack(data.status(), buffer, position, comm);
}

void pack(const Segment& data,
          std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    pack(data.segmentNumber(), buffer, position, comm);
    pack(data.branchNumber(), buffer, position, comm);
    pack(data.outletSegment(), buffer, position, comm);
    pack(data.inletSegments(), buffer, position, comm);
    pack(data.totalLength(), buffer, position, comm);
    pack(data.depth(), buffer, position, comm);
    pack(data.internalDiameter(), buffer, position, comm);
    pack(data.roughness(), buffer, position, comm);
    pack(data.crossArea(), buffer, position, comm);
    pack(data.volume(), buffer, position, comm);
    pack(data.dataReady(), buffer, position, comm);
    pack(data.segmentType(), buffer, position, comm);
    pack(data.spiralICD(), buffer, position, comm);
    pack(data.getValve(), buffer, position, comm);
}

template<class T>
void pack(const std::shared_ptr<T>& data, std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    pack(data != nullptr, buffer, position, comm);
    if (data)
        pack(*data, buffer, position, comm);
}

template<class T>
void pack(const std::unique_ptr<T>& data, std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    pack(data != nullptr, buffer, position, comm);
    if (data)
        pack(*data, buffer, position, comm);
}

void pack(const Dimension& data,
          std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    pack(data.getSIScalingRaw(), buffer, position, comm);
    pack(data.getSIOffset(), buffer, position, comm);
}

void pack(const UnitSystem& data,
          std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    pack(data.getName(), buffer, position, comm);
    pack(data.getType(), buffer, position, comm);
    pack(data.getDimensions(), buffer, position, comm);
    pack(data.use_count(), buffer, position, comm);
}

void pack(const WellSegments& data,
          std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    pack(data.compPressureDrop(), buffer, position, comm);
    pack(data.segments(), buffer, position, comm);
}

void pack(const Well& data,
          std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    pack(data.name(), buffer, position, comm);
    pack(data.groupName(), buffer, position, comm);
    pack(data.firstTimeStep(), buffer, position, comm);
    pack(data.seqIndex(), buffer, position, comm);
    pack(data.getHeadI(), buffer, position, comm);
    pack(data.getHeadJ(), buffer, position, comm);
    pack(data.getRefDepth(), buffer, position, comm);
    pack(data.wellType(), buffer, position, comm);
    pack(data.getWellConnectionOrdering(), buffer, position, comm);
    pack(data.units(), buffer, position, comm);
    pack(data.udqUndefined(), buffer, position, comm);
    pack(data.getStatus(), buffer, position, comm);
    pack(data.getDrainageRadius(), buffer, position, comm);
    pack(data.getAllowCrossFlow(), buffer, position, comm);
    pack(data.getAutomaticShutIn(), buffer, position, comm);
    pack(data.wellGuideRate(), buffer, position, comm);
    pack(data.getEfficiencyFactor(), buffer, position, comm);
    pack(data.getSolventFraction(), buffer, position, comm);
    pack(data.predictionMode(), buffer, position, comm);
    pack(data.getEconLimits(), buffer, position, comm);
    pack(data.getFoamProperties(), buffer, position, comm);
    pack(data.getPolymerProperties(), buffer, position, comm);
    pack(data.getBrineProperties(), buffer, position, comm);
    pack(data.getTracerProperties(), buffer, position, comm);
    pack(data.getConnections(), buffer, position, comm);
    pack(data.getProductionProperties(), buffer, position, comm);
    pack(data.getInjectionProperties(), buffer, position, comm);
    pack(data.hasSegments(), buffer, position, comm);
    if (data.hasSegments())
        pack(data.getSegments(), buffer, position, comm);
}

template<class T>
void pack(const IOrderSet<T>& data, std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    pack(data.index(), buffer, position, comm);
    pack(data.data(), buffer, position, comm);
}

void pack(const Group::GroupInjectionProperties& data,
          std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    pack(data.phase, buffer, position, comm);
    pack(data.cmode, buffer, position, comm);
    pack(data.surface_max_rate, buffer, position, comm);
    pack(data.resv_max_rate, buffer, position, comm);
    pack(data.target_reinj_fraction, buffer, position, comm);
    pack(data.target_void_fraction, buffer, position, comm);
    pack(data.reinj_group, buffer, position, comm);
    pack(data.voidage_group, buffer, position, comm);
    pack(data.injection_controls, buffer, position, comm);
}

void pack(const Group::GroupProductionProperties& data,
          std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    pack(data.cmode, buffer, position, comm);
    pack(data.exceed_action, buffer, position, comm);
    pack(data.oil_target, buffer, position, comm);
    pack(data.water_target, buffer, position, comm);
    pack(data.gas_target, buffer, position, comm);
    pack(data.liquid_target, buffer, position, comm);
    pack(data.guide_rate, buffer, position, comm);
    pack(data.guide_rate_def, buffer, position, comm);
    pack(data.resv_target, buffer, position, comm);
    pack(data.production_controls, buffer, position, comm);
}

void pack(const Group& data,
          std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    pack(data.name(), buffer, position, comm);
    pack(data.insert_index(), buffer, position, comm);
    pack(data.initStep(), buffer, position, comm);
    pack(data.udqUndefined(), buffer, position, comm);
    pack(data.units(), buffer, position, comm);
    pack(data.type(), buffer, position, comm);
    pack(data.getGroupEfficiencyFactor(), buffer, position, comm);
    pack(data.getTransferGroupEfficiencyFactor(), buffer, position, comm);
    pack(data.isAvailableForGroupControl(), buffer, position, comm);
    pack(data.getGroupNetVFPTable(), buffer, position, comm);
    pack(data.parent(), buffer, position, comm);
    pack(data.iwells(), buffer, position, comm);
    pack(data.igroups(), buffer, position, comm);
    pack(data.injectionProperties(), buffer, position, comm);
    pack(data.productionProperties(), buffer, position, comm);
}

void pack(const WList& data,
          std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    pack(data.wellList(), buffer, position, comm);
}

void pack(const WListManager& data,
          std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    pack(data.lists(), buffer, position, comm);
}

void pack(const GuideRateModel& data,
          std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    pack(data.timeInterval(), buffer, position, comm);
    pack(data.target(), buffer, position, comm);
    pack(data.coefs(), buffer, position, comm);
    pack(data.allow_increase(), buffer, position, comm);
    pack(data.damping_factor(), buffer, position, comm);
    pack(data.free_gas(), buffer, position, comm);
    pack(data.defaultModel(), buffer, position, comm);
    pack(data.udaCoefs(), buffer, position, comm);
}

void pack(const GuideRateConfig& data,
          std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    pack(data.getModel(), buffer, position, comm);
    pack(data.getWells(), buffer, position, comm);
    pack(data.getGroups(), buffer, position, comm);
}

void pack(const GConSale::GCONSALEGroup& data,
          std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    pack(data.sales_target, buffer, position, comm);
    pack(data.max_sales_rate, buffer, position, comm);
    pack(data.min_sales_rate, buffer, position, comm);
    pack(data.max_proc, buffer, position, comm);
    pack(data.udq_undefined, buffer, position, comm);
    pack(data.unit_system, buffer, position, comm);
}

void pack(const GConSale& data,
          std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    pack(data.getGroups(), buffer, position, comm);
}

void pack(const GConSump::GCONSUMPGroup& data,
          std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    pack(data.consumption_rate, buffer, position, comm);
    pack(data.import_rate, buffer, position, comm);
    pack(data.network_node, buffer, position, comm);
    pack(data.udq_undefined, buffer, position, comm);
    pack(data.unit_system, buffer, position, comm);
}

void pack(const GConSump& data,
          std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    pack(data.getGroups(), buffer, position, comm);
}

void pack(const RFTConfig& data,
          std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    pack(data.timeMap(), buffer, position, comm);
    pack(data.firstRFTOutput(), buffer, position, comm);
    pack(data.wellOpenRftTime(), buffer, position, comm);
    pack(data.wellOpenRftName(), buffer, position, comm);
    pack(data.wellOpen(), buffer, position, comm);
    pack(data.rftConfig(), buffer, position, comm);
    pack(data.pltConfig(), buffer, position, comm);
}

void pack(const DeckItem& data,
          std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    pack(data.dVal(), buffer, position, comm);
    pack(data.iVal(), buffer, position, comm);
    pack(data.sVal(), buffer, position, comm);
    pack(data.uVal(), buffer, position, comm);
    pack(data.getType(), buffer, position, comm);
    pack(data.name(), buffer, position, comm);
    pack(data.valueStatus(), buffer, position, comm);
    pack(data.rawData(), buffer, position, comm);
    pack(data.activeDimensions(), buffer, position, comm);
    pack(data.defaultDimensions(), buffer, position, comm);
}

void pack(const DeckRecord& data,
          std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    pack(data.getItems(), buffer, position, comm);
}

void pack(const Location& data,
          std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    pack(data.filename, buffer, position, comm);
    pack(data.lineno, buffer, position, comm);
}

void pack(const DeckKeyword& data,
          std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    pack(data.name(), buffer, position, comm);
    pack(data.location(), buffer, position, comm);
    pack(data.records(), buffer, position, comm);
    pack(data.isDataKeyword(), buffer, position, comm);
    pack(data.isSlashTerminated(), buffer, position, comm);
}

void pack(const Deck& data,
          std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    pack(data.keywords(), buffer, position, comm);
    pack(data.getDefaultUnitSystem(), buffer, position, comm);
    pack(data.activeUnitSystem(), buffer, position, comm);
    pack(data.getDataFile(), buffer, position, comm);
    pack(data.getInputPath(), buffer, position, comm);
    pack(data.unitSystemAccessCount(), buffer, position, comm);
}

void pack(const Action::ASTNode& data,
          std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    pack(data.type, buffer, position, comm);
    pack(data.func_type, buffer, position, comm);
    pack(data.func, buffer, position, comm);
    pack(data.argList(), buffer, position, comm);
    pack(data.getNumber(), buffer, position, comm);
    pack(data.childrens(), buffer, position, comm);
}

void pack(const Action::AST& data,
          std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    pack(data.getCondition(), buffer, position, comm);
}

void pack(const Action::Quantity& data,
          std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    pack(data.quantity, buffer, position, comm);
    pack(data.args, buffer, position, comm);
}

void pack(const Action::Condition& data,
          std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    pack(data.lhs, buffer, position, comm);
    pack(data.rhs, buffer, position, comm);
    pack(data.logic, buffer, position, comm);
    pack(data.cmp, buffer, position, comm);
    pack(data.cmp_string, buffer, position, comm);
}

void pack(const Action::ActionX& data,
          std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    pack(data.name(), buffer, position, comm);
    pack(data.max_run(), buffer, position, comm);
    pack(data.min_wait(), buffer, position, comm);
    pack(data.start_time(), buffer, position, comm);
    pack(data.getKeywords(), buffer, position, comm);
    pack(data.getCondition(), buffer, position, comm);
    pack(data.conditions(), buffer, position, comm);
    pack(data.getRunCount(), buffer, position, comm);
    pack(data.getLastRun(), buffer, position, comm);
}

void pack(const Action::Actions& data,
          std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    pack(data.getActions(), buffer, position, comm);
}

void pack(const RestartSchedule& data,
          std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    pack(data.timestep, buffer, position, comm);
    pack(data.basic, buffer, position, comm);
    pack(data.frequency, buffer, position, comm);
    pack(data.rptsched_restart_set, buffer, position, comm);
    pack(data.rptsched_restart, buffer, position, comm);
}

void pack(const TimeStampUTC& data,
          std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    pack(data.ymd(), buffer, position, comm);
    pack(data.hour(), buffer, position, comm);
    pack(data.minutes(), buffer, position, comm);
    pack(data.seconds(), buffer, position, comm);
    pack(data.microseconds(), buffer, position, comm);
}

void pack(const WellPolymerProperties& data,
          std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    pack(data.m_polymerConcentration, buffer, position, comm);
    pack(data.m_saltConcentration, buffer, position, comm);
    pack(data.m_plymwinjtable, buffer, position, comm);
    pack(data.m_skprwattable, buffer, position, comm);
    pack(data.m_skprpolytable, buffer, position, comm);
}

void pack(const Well::WellGuideRate& data,
          std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    pack(data.available, buffer, position, comm);
    pack(data.guide_rate, buffer, position, comm);
    pack(data.guide_phase, buffer, position, comm);
    pack(data.scale_factor, buffer, position, comm);
}

void pack(const GuideRateConfig::WellTarget& data,
          std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    pack(data.guide_rate, buffer, position, comm);
    pack(data.target, buffer, position, comm);
    pack(data.scaling_factor, buffer, position, comm);
}

void pack(const GuideRateConfig::GroupTarget& data,
          std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    pack(data.guide_rate, buffer, position, comm);
    pack(data.target, buffer, position, comm);
}

void pack(const PlyvmhTable& data, std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    pack(static_cast<const std::vector<PlyvmhRecord>&>(data), buffer, position, comm);
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

template<class K, class C, class A>
void unpack(std::set<K,C,A>& data,
            std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm)
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
            Dune::MPIHelper::MPICommunicator comm)
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
          Dune::MPIHelper::MPICommunicator comm)
{
    for (T& entry : data)
        unpack(entry, buffer, position, comm);
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
    std::vector<T> unique;
    std::vector<int> indices;
    Opm::Mpi::unpack(unique, buffer, position, comm);
    Opm::Mpi::unpack(indices, buffer, position, comm);
    reconstructDynState(unique, indices, data);
}

template<class T>
void unpack(DynamicVector<T>& data, std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm)
{
    std::vector<T> ddata;
    unpack(ddata, buffer, position, comm);
    data = DynamicVector<T>(ddata);
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
    str.clear();
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
    unpack(data.current_control, buffer, position, comm);
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

void unpack(WellType& data, std::vector<char>& buffer, int& position, Dune::MPIHelper::MPICommunicator comm)
{
    Phase preferred_phase;
    bool producer;
    unpack(producer, buffer, position, comm);
    unpack(preferred_phase, buffer, position, comm);
    data = WellType( producer, preferred_phase );
}

void unpack(TimeMap& data, std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm)
{
    std::vector<std::time_t> timeList;
    unpack(timeList, buffer, position, comm);

    data = TimeMap(timeList);
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

void unpack(TableManager& data, std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    std::map<std::string, TableContainer> simpleTables;
    SplitSimpleTables split;
    std::vector<PvtgTable> pvtgTables;
    std::vector<PvtoTable> pvtoTables;
    std::vector<Rock2dTable> rock2dTables;
    std::vector<Rock2dtrTable> rock2dtrTables;
    PvtwTable pvtwTable;
    PvcdoTable pvcdoTable;
    DensityTable densityTable;
    PlyvmhTable plyvmhTable;
    RockTable rockTable;
    ViscrefTable viscrefTable;
    PlmixparTable plmixparTable;
    ShrateTable shrateTable;
    Stone1exTable stone1exTable;
    TlmixparTable tlmixparTable;
    WatdentTable watdentTable;
    std::vector<PvtwsaltTable> pvtwsaltTables;
    std::vector<BrineDensityTable> bdensityTables;
    std::vector<SolventDensityTable> sdensityTables;
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
    bool hasShrate;
    DenT oilDenT, gasDenT, watDenT;
    StandardCond stcond;
    std::size_t gas_comp_index;
    std::shared_ptr<JFunc> jfunc;
    double rtemp;

    unpack(simpleTables, buffer, position, comm);
    unpack(split.plyshMax, buffer, position, comm);
    unpack(split.plyshMap, buffer, position, comm);
    unpack(split.rockMax, buffer, position, comm);
    unpack(split.rockMap, buffer, position, comm);
    unpack(pvtgTables, buffer, position, comm);
    unpack(pvtoTables, buffer, position, comm);
    unpack(rock2dTables, buffer, position, comm);
    unpack(rock2dtrTables, buffer, position, comm);
    unpack(pvtwTable, buffer, position, comm);
    unpack(pvcdoTable, buffer, position, comm);
    unpack(densityTable, buffer, position, comm);
    unpack(plyvmhTable, buffer, position, comm);
    unpack(rockTable, buffer, position, comm);
    unpack(plmixparTable, buffer, position, comm);
    unpack(shrateTable, buffer, position, comm);
    unpack(stone1exTable, buffer, position, comm);
    unpack(tlmixparTable, buffer, position, comm);
    unpack(viscrefTable, buffer, position, comm);
    unpack(watdentTable, buffer, position, comm);
    unpack(pvtwsaltTables, buffer, position, comm);
    unpack(bdensityTables, buffer, position, comm);
    unpack(sdensityTables, buffer, position, comm);
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
    unpack(hasShrate, buffer, position, comm);
    bool hasJf;
    unpack(hasJf, buffer, position, comm);
    if (hasJf) {
        jfunc = std::make_shared<JFunc>();
        unpack(*jfunc, buffer, position, comm);
    }
    unpack(oilDenT, buffer, position, comm);
    unpack(gasDenT, buffer, position, comm);
    unpack(watDenT, buffer, position, comm);
    unpack(stcond, buffer, position, comm);
    unpack(gas_comp_index, buffer, position, comm);
    unpack(rtemp, buffer, position, comm);

    if (split.plyshMax > 0) {
        TableContainer container(split.plyshMax);
        for (const auto& it : split.plyshMap) {
            container.addTable(it.first, it.second);
        }
        simpleTables.insert(std::make_pair("PLYSHLOG", container));
    }
    if (split.rockMax > 0) {
        TableContainer container(split.rockMax);
        for (const auto& it : split.rockMap) {
            container.addTable(it.first, it.second);
        }
        simpleTables.insert(std::make_pair("ROCKTAB", container));
    }

    data = TableManager(simpleTables, pvtgTables, pvtoTables, rock2dTables,
                        rock2dtrTables, pvtwTable, pvcdoTable, densityTable,
                        plyvmhTable, rockTable, plmixparTable, shrateTable, stone1exTable,
                        tlmixparTable, viscrefTable, watdentTable, pvtwsaltTables,
                        bdensityTables, sdensityTables, plymwinjTables, skprwatTables,
                        skprpolyTables, tabdims, regdims, eqldims, aqudims, hasImptvd,
                        hasEntpvd, hasEqlnum, hasShrate, jfunc, oilDenT, gasDenT,
                        watDenT, stcond, gas_comp_index, rtemp);
}

void unpack(OilVaporizationProperties& data,
          std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    OilVaporizationProperties::OilVaporization type;
    double vap1, vap2;
    std::vector<double> maxDRSDT, maxDRVDT;
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

void unpack(Events& data,
          std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    DynamicVector<uint64_t> events;
    unpack(events, buffer, position, comm);
    data = Events(events);
}

void unpack(MessageLimits& data,
          std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    DynamicState<MLimits> limits;
    unpack(limits, buffer, position, comm);
    data = MessageLimits(limits);
}

void unpack(VFPInjTable& data,
            std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm)
{
    int tableNum;
    double datumDepth;
    VFPInjTable::FLO_TYPE floType;
    std::vector<double> floAxis, thpAxis;
    VFPInjTable::array_type table;

    unpack(tableNum, buffer, position, comm);
    unpack(datumDepth, buffer, position, comm);
    unpack(floType, buffer, position, comm);
    unpack(floAxis, buffer, position, comm);
    unpack(thpAxis, buffer, position, comm);
    unpack(table, buffer, position, comm);

    data = VFPInjTable(tableNum, datumDepth, floType,
                       floAxis, thpAxis, table);
}

void unpack(VFPProdTable& data,
            std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm)
{
    int tableNum;
    double datumDepth;
    VFPProdTable::FLO_TYPE floType;
    VFPProdTable::WFR_TYPE wfrType;
    VFPProdTable::GFR_TYPE gfrType;
    VFPProdTable::ALQ_TYPE alqType;
    std::vector<double> floAxis, thpAxis, wfrAxis, gfrAxis, alqAxis;
    VFPProdTable::array_type table;

    unpack(tableNum, buffer, position, comm);
    unpack(datumDepth, buffer, position, comm);
    unpack(floType, buffer, position, comm);
    unpack(wfrType, buffer, position, comm);
    unpack(gfrType, buffer, position, comm);
    unpack(alqType, buffer, position, comm);
    unpack(floAxis, buffer, position, comm);
    unpack(thpAxis, buffer, position, comm);
    unpack(wfrAxis, buffer, position, comm);
    unpack(gfrAxis, buffer, position, comm);
    unpack(alqAxis, buffer, position, comm);
    unpack(table, buffer, position, comm);

    data = VFPProdTable(tableNum, datumDepth, floType, wfrType,
                        gfrType, alqType, floAxis, thpAxis,
                        wfrAxis, gfrAxis, alqAxis, table);
}

void unpack(WellTestConfig::WTESTWell& data,
            std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm)
{
    unpack(data.name, buffer, position, comm);
    unpack(data.shut_reason, buffer, position, comm);
    unpack(data.test_interval, buffer, position, comm);
    unpack(data.num_test, buffer, position, comm);
    unpack(data.startup_time, buffer, position, comm);
    unpack(data.begin_report_step, buffer, position, comm);
}

void unpack(WellTestConfig& data,
            std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm)
{
    std::vector<WellTestConfig::WTESTWell> ddata;
    unpack(ddata, buffer, position, comm);
    data = WellTestConfig(ddata);
}

void unpack(WellTracerProperties& data,
            std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm)
{
    WellTracerProperties::ConcentrationMap ddata;
    unpack(ddata, buffer, position, comm);
    data = WellTracerProperties(ddata);
}

void unpack(UDAValue& data,
            std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm)
{
    bool isDouble;
    Dimension dim;
    unpack(dim, buffer, position, comm);
    unpack(isDouble, buffer, position, comm);
    if (isDouble) {
        double val;
        unpack(val, buffer, position, comm);
        data = UDAValue(val, dim);
    } else {
        std::string val;
        unpack(val, buffer, position, comm);
        data = UDAValue(val, dim);
    }
}

void unpack(Connection& data,
            std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm)
{
    Connection::Direction dir;
    double depth;
    Connection::State state;
    int satTableId, complnum;
    double CF, Kh, rw, r0, skinFactor;
    int I, J, K;
    size_t seqIndex;
    double segDistStart, segDistEnd;
    bool defaultSatTabId;
    size_t compSegSeqIndex;
    int segment;
    double wellPi;
    Connection::CTFKind kind;

    unpack(dir, buffer, position, comm);
    unpack(depth, buffer, position, comm);
    unpack(state, buffer, position, comm);
    unpack(satTableId, buffer, position, comm);
    unpack(complnum, buffer, position, comm);
    unpack(CF, buffer, position, comm);
    unpack(Kh, buffer, position, comm);
    unpack(rw, buffer, position, comm);
    unpack(r0, buffer, position, comm);
    unpack(skinFactor, buffer, position, comm);
    unpack(I, buffer, position, comm);
    unpack(J, buffer, position, comm);
    unpack(K, buffer, position, comm);
    unpack(kind, buffer, position, comm);
    unpack(seqIndex, buffer, position, comm);
    unpack(segDistStart, buffer, position, comm);
    unpack(segDistEnd, buffer, position, comm);
    unpack(defaultSatTabId, buffer, position, comm);
    unpack(compSegSeqIndex, buffer, position, comm);
    unpack(segment, buffer, position, comm);
    unpack(wellPi, buffer, position, comm);

    data = Connection(dir, depth, state, satTableId,
                      complnum, CF, Kh, rw, r0,
                      skinFactor, {I,J,K}, kind, seqIndex,
                      segDistStart, segDistEnd,
                      defaultSatTabId, compSegSeqIndex,
                      segment, wellPi);
}

void unpack(Well::WellInjectionProperties& data,
            std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm)
{
    unpack(data.name, buffer, position, comm);
    unpack(data.surfaceInjectionRate, buffer, position, comm);
    unpack(data.reservoirInjectionRate, buffer, position, comm);
    unpack(data.BHPTarget, buffer, position, comm);
    unpack(data.THPTarget, buffer, position, comm);
    unpack(data.bhp_hist_limit, buffer, position, comm);
    unpack(data.thp_hist_limit, buffer, position, comm);
    unpack(data.temperature, buffer, position, comm);
    unpack(data.BHPH, buffer, position, comm);
    unpack(data.THPH, buffer, position, comm);
    unpack(data.VFPTableNumber, buffer, position, comm);
    unpack(data.predictionMode, buffer, position, comm);
    unpack(data.injectionControls, buffer, position, comm);
    unpack(data.injectorType, buffer, position, comm);
    unpack(data.controlMode, buffer, position, comm);
}

void unpack(WellEconProductionLimits& data,
            std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm)
{
    double minOilRate, minGasRate, maxWaterCut, maxGasOilRatio, maxWaterGasRatio;
    WellEconProductionLimits::EconWorkover workover, workoverSecondary;
    bool endRun;
    std::string followonWell;
    WellEconProductionLimits::QuantityLimit quantityLimit;
    double secondaryMaxWaterCut, maxGasLiquidRatio, minLiquidRate,
           maxTemperature, minReservoirFluidRate;
    unpack(minOilRate, buffer, position, comm);
    unpack(minGasRate, buffer, position, comm);
    unpack(maxWaterCut, buffer, position, comm);
    unpack(maxGasOilRatio, buffer, position, comm);
    unpack(maxWaterGasRatio, buffer, position, comm);
    unpack(workover, buffer, position, comm);
    unpack(endRun, buffer, position, comm);
    unpack(followonWell, buffer, position, comm);
    unpack(quantityLimit, buffer, position, comm);
    unpack(secondaryMaxWaterCut, buffer, position, comm);
    unpack(workoverSecondary, buffer, position, comm);
    unpack(maxGasLiquidRatio, buffer, position, comm);
    unpack(minLiquidRate, buffer, position, comm);
    unpack(maxTemperature, buffer, position, comm);
    unpack(minReservoirFluidRate, buffer, position, comm);
    data = WellEconProductionLimits(minOilRate, minGasRate, maxWaterCut,
                                    maxGasOilRatio, maxWaterGasRatio,
                                    workover, endRun, followonWell,
                                    quantityLimit, secondaryMaxWaterCut,
                                    workoverSecondary, maxGasLiquidRatio,
                                    minLiquidRate, maxTemperature,
                                    minReservoirFluidRate);
}

void unpack(WellConnections& data,
            std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm)
{
    int headI, headJ;
    size_t numRemoved;
    std::vector<Connection> connections;

    unpack(headI, buffer, position, comm),
    unpack(headJ, buffer, position, comm),
    unpack(numRemoved, buffer, position, comm),
    unpack(connections, buffer, position, comm),

    data = WellConnections(headI, headJ, numRemoved, connections);
}

void unpack(Well::WellProductionProperties& data,
            std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm)
{
    std::string name;
    UDAValue OilRate, WaterRate, GasRate, LiquidRate, ResVRate;
    UDAValue BHPTarget, THPTarget;
    double bhp_hist_limit, thp_hist_limit;
    double BHPH, THPH;
    int VFPTableNumber;
    double ALQValue;
    bool predictionMode;
    Well::ProducerCMode controlMode, whistctl_cmode;
    int prodCtrls;

    unpack(name, buffer, position, comm);
    unpack(OilRate, buffer, position, comm);
    unpack(WaterRate, buffer, position, comm);
    unpack(GasRate, buffer, position, comm);
    unpack(LiquidRate, buffer, position, comm);
    unpack(ResVRate, buffer, position, comm);
    unpack(BHPTarget, buffer, position, comm);
    unpack(THPTarget, buffer, position, comm);
    unpack(bhp_hist_limit, buffer, position, comm);
    unpack(thp_hist_limit, buffer, position, comm);
    unpack(BHPH, buffer, position, comm);
    unpack(THPH, buffer, position, comm);
    unpack(VFPTableNumber, buffer, position, comm);
    unpack(ALQValue, buffer, position, comm);
    unpack(predictionMode, buffer, position, comm);
    unpack(controlMode, buffer, position, comm);
    unpack(whistctl_cmode, buffer, position, comm);
    unpack(prodCtrls, buffer, position, comm);
    data = Well::WellProductionProperties(name, OilRate, WaterRate, GasRate,
                                          LiquidRate, ResVRate, BHPTarget,
                                          THPTarget, bhp_hist_limit, thp_hist_limit,
                                          BHPH, THPH, VFPTableNumber,
                                          ALQValue, predictionMode, controlMode,
                                          whistctl_cmode, prodCtrls);
}

void unpack(SpiralICD& data,
            std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm)
{
    double strength, length, densityCalibration,
           viscosityCalibration, criticalValue,
           widthTransitionRegion, maxViscosityRatio;
    int methodFlowScaling;
    double maxAbsoluteRate;
    ICDStatus status;
    double scalingFactor;

    unpack(strength, buffer, position, comm);
    unpack(length, buffer, position, comm);
    unpack(densityCalibration, buffer, position, comm);
    unpack(viscosityCalibration, buffer, position, comm);
    unpack(criticalValue, buffer, position, comm);
    unpack(widthTransitionRegion, buffer, position, comm);
    unpack(maxViscosityRatio, buffer, position, comm);
    unpack(methodFlowScaling, buffer, position, comm);
    unpack(maxAbsoluteRate, buffer, position, comm);
    unpack(status, buffer, position, comm);
    unpack(scalingFactor, buffer, position, comm);

    data = SpiralICD(strength, length, densityCalibration,
                     viscosityCalibration, criticalValue,
                     widthTransitionRegion, maxViscosityRatio,
                     methodFlowScaling, maxAbsoluteRate,
                     status, scalingFactor);
}

void unpack(Valve& data,
            std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm)
{
    double conFlowCoefficient;
    double conCrossArea;
    double conMaxCrossArea;
    double pipeAdditionalLength;
    double pipeDiameter;
    double pipeRoughness;
    double pipeCrossArea;
    ICDStatus status;

    unpack(conFlowCoefficient, buffer, position, comm);
    unpack(conCrossArea, buffer, position, comm);
    unpack(conMaxCrossArea, buffer, position, comm);
    unpack(pipeAdditionalLength, buffer, position, comm);
    unpack(pipeDiameter, buffer, position, comm);
    unpack(pipeRoughness, buffer, position, comm);
    unpack(pipeCrossArea, buffer, position, comm);
    unpack(status, buffer, position, comm);
    data = Valve(conFlowCoefficient, conCrossArea, conMaxCrossArea,
                 pipeAdditionalLength, pipeDiameter, pipeRoughness,
                 pipeCrossArea, status);
}

void unpack(Segment& data,
            std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm)
{
    int segmentNumber, branchNumber, outletSegment;
    std::vector<int> inletSegments;
    double totalLength, depth, internalDiameter, roughness, crossArea, volume;
    bool dataReady;
    Segment::SegmentType segmentType;
    std::shared_ptr<SpiralICD> spiralICD;
    std::shared_ptr<Valve> valve;

    unpack(segmentNumber, buffer, position, comm);
    unpack(branchNumber, buffer, position, comm);
    unpack(outletSegment, buffer, position, comm);
    unpack(inletSegments, buffer, position, comm);
    unpack(totalLength, buffer, position, comm);
    unpack(depth, buffer, position, comm);
    unpack(internalDiameter, buffer, position, comm);
    unpack(roughness, buffer, position, comm);
    unpack(crossArea, buffer, position, comm);
    unpack(volume, buffer, position, comm);
    unpack(dataReady, buffer, position, comm);
    unpack(segmentType, buffer, position, comm);
    unpack(spiralICD, buffer, position, comm);
    unpack(valve, buffer, position, comm);
    data = Segment(segmentNumber, branchNumber, outletSegment,
                   inletSegments, totalLength, depth,
                   internalDiameter, roughness, crossArea,
                   volume, dataReady, segmentType, spiralICD, valve);
}

template<class T>
void unpack(std::shared_ptr<T>& data, std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    bool hasVal;
    unpack(hasVal, buffer, position, comm);
    if (hasVal) {
        data = std::make_shared<T>();
        unpack(*data, buffer, position, comm);
    }
}

template<class T>
void unpack(std::unique_ptr<T>& data, std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    bool hasVal;
    unpack(hasVal, buffer, position, comm);
    if (hasVal) {
        data.reset(new T);
        unpack(*data, buffer, position, comm);
    }
}

void unpack(Dimension& data,
            std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm)
{
    double siScaling, siOffset;

    unpack(siScaling, buffer, position, comm);
    unpack(siOffset, buffer, position, comm);
    data = Dimension(siScaling, siOffset);
}

void unpack(UnitSystem& data,
            std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm)
{
    std::string name;
    UnitSystem::UnitType type;
    std::map<std::string, Dimension> dimensions;
    size_t use_count;
    unpack(name, buffer, position, comm);
    unpack(type, buffer, position, comm);
    unpack(dimensions, buffer, position, comm);
    unpack(use_count, buffer, position, comm);

    data = UnitSystem(name, type, dimensions, use_count);
}

void unpack(WellSegments& data,
            std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm)
{
    WellSegments::CompPressureDrop compPressureDrop;
    std::vector<Segment> segments;

    unpack(compPressureDrop, buffer, position, comm);
    unpack(segments, buffer, position, comm);

    data = WellSegments(compPressureDrop, segments);
}

void unpack(Well& data,
            std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm)
{
    std::string name, groupName;
    std::size_t firstTimeStep, seqIndex;
    int headI, headJ;
    double ref_depth;
    WellType wtype;
    Connection::Order ordering;
    UnitSystem units;
    double udq_undefined;
    Well::Status status;
    double drainageRadius;
    bool allowCrossFlow, automaticShutIn;
    Well::WellGuideRate guideRate;
    double efficiencyFactor;
    double solventFraction;
    bool prediction_mode;
    auto econLimits = std::make_shared<WellEconProductionLimits>();
    auto foamProperties = std::make_shared<WellFoamProperties>();
    auto polymerProperties = std::make_shared<WellPolymerProperties>();
    auto brineProperties = std::make_shared<WellBrineProperties>();
    auto tracerProperties = std::make_shared<WellTracerProperties>();
    auto connection = std::make_shared<WellConnections>();
    auto production = std::make_shared<Well::WellProductionProperties>();
    auto injection = std::make_shared<Well::WellInjectionProperties>();
    std::shared_ptr<WellSegments> segments;

    unpack(name, buffer, position, comm);
    unpack(groupName, buffer, position, comm);
    unpack(firstTimeStep, buffer, position, comm);
    unpack(seqIndex, buffer, position, comm);
    unpack(headI, buffer, position, comm);
    unpack(headJ, buffer, position, comm);
    unpack(ref_depth, buffer, position, comm);
    unpack(wtype, buffer, position, comm);
    unpack(ordering, buffer, position, comm);
    unpack(units, buffer, position, comm);
    unpack(udq_undefined, buffer, position, comm);
    unpack(status, buffer, position, comm);
    unpack(drainageRadius, buffer, position, comm);
    unpack(allowCrossFlow, buffer, position, comm);
    unpack(automaticShutIn, buffer, position, comm);
    unpack(guideRate, buffer, position, comm);
    unpack(efficiencyFactor, buffer, position, comm);
    unpack(solventFraction, buffer, position, comm);
    unpack(prediction_mode, buffer, position, comm);
    unpack(*econLimits, buffer, position, comm);
    unpack(*foamProperties, buffer, position, comm);
    unpack(*polymerProperties, buffer, position, comm);
    unpack(*brineProperties, buffer, position, comm);
    unpack(*tracerProperties, buffer, position, comm);
    unpack(*connection, buffer, position, comm);
    unpack(*production, buffer, position, comm);
    unpack(*injection, buffer, position, comm);
    bool hasSegments;
    unpack(hasSegments, buffer, position, comm);
    if (hasSegments) {
        segments = std::make_shared<WellSegments>();
        unpack(*segments, buffer, position, comm);
    }
    data = Well(name, groupName, firstTimeStep, seqIndex, headI, headJ,
                ref_depth, wtype, ordering, units, udq_undefined, status,
                drainageRadius, allowCrossFlow, automaticShutIn,
                guideRate, efficiencyFactor, solventFraction, prediction_mode,
                econLimits, foamProperties, polymerProperties, brineProperties,
                tracerProperties, connection, production, injection, segments);
}

template<class T>
void unpack(IOrderSet<T>& data, std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm)
{
    typename IOrderSet<T>::index_type index;
    typename IOrderSet<T>::storage_type storage;
    unpack(index, buffer, position, comm);
    unpack(storage, buffer, position, comm);
    data = IOrderSet<T>(index, storage);
}

template void unpack(std::map<Phase,Group::GroupInjectionProperties>& data,
            std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm);

void unpack(Group::GroupInjectionProperties& data,
            std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm)
{
    unpack(data.phase, buffer, position, comm);
    unpack(data.cmode, buffer, position, comm);
    unpack(data.surface_max_rate, buffer, position, comm);
    unpack(data.resv_max_rate, buffer, position, comm);
    unpack(data.target_reinj_fraction, buffer, position, comm);
    unpack(data.target_void_fraction, buffer, position, comm);
    unpack(data.reinj_group, buffer, position, comm);
    unpack(data.voidage_group, buffer, position, comm);
    unpack(data.injection_controls, buffer, position, comm);
}


void unpack(Group::GroupProductionProperties& data,
            std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm)
{
    unpack(data.cmode, buffer, position, comm);
    unpack(data.exceed_action, buffer, position, comm);
    unpack(data.oil_target, buffer, position, comm);
    unpack(data.water_target, buffer, position, comm);
    unpack(data.gas_target, buffer, position, comm);
    unpack(data.liquid_target, buffer, position, comm);
    unpack(data.guide_rate, buffer, position, comm);
    unpack(data.guide_rate_def, buffer, position, comm);
    unpack(data.resv_target, buffer, position, comm);
    unpack(data.production_controls, buffer, position, comm);
}

void unpack(Group& data,
            std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm)
{
    std::string name;
    std::size_t insert_index, initStep;
    double udqUndefined;
    UnitSystem units;
    Group::GroupType type;
    double groupEfficiencyFactor;
    bool transferGroupEfficiencyFactor;
    bool availableForGroupControl;
    int groupNetVFPTable;
    std::string parent;
    IOrderSet<std::string> wells, groups;
    std::map<Phase, Group::GroupInjectionProperties> injection;
    Group::GroupProductionProperties production;

    unpack(name, buffer, position, comm);
    unpack(insert_index, buffer, position, comm);
    unpack(initStep, buffer, position, comm);
    unpack(udqUndefined, buffer, position, comm);
    unpack(units, buffer, position, comm);
    unpack(type, buffer, position, comm);
    unpack(groupEfficiencyFactor, buffer, position, comm);
    unpack(transferGroupEfficiencyFactor, buffer, position, comm);
    unpack(availableForGroupControl, buffer, position, comm);
    unpack(groupNetVFPTable, buffer, position, comm);
    unpack(parent, buffer, position, comm);
    unpack(wells, buffer, position, comm);
    unpack(groups, buffer, position, comm);
    unpack(injection, buffer, position, comm);
    unpack(production, buffer, position, comm);
    data = Group(name, insert_index, initStep, udqUndefined,
                 units, type, groupEfficiencyFactor,
                 transferGroupEfficiencyFactor,
                 availableForGroupControl,
                 groupNetVFPTable, parent, wells, groups,
                 injection, production);
}

void unpack(WList& data,
            std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm)
{
    WList::storage ddata;
    unpack(ddata, buffer, position, comm);
    data = WList(ddata);
}

void unpack(WListManager& data,
            std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm)
{
    std::map<std::string,WList> lists;
    unpack(lists, buffer, position, comm);
    data = WListManager(lists);
}

void unpack(GuideRateModel& data,
            std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm)
{
    double timeInterval;
    GuideRateModel::Target target;
    std::array<double,6> coefs;
    bool allow_increase, free_gas, defaultModel;
    double damping_factor;
    std::array<UDAValue,3> udaCoefs;

    unpack(timeInterval, buffer, position, comm);
    unpack(target, buffer, position, comm);
    unpack(coefs, buffer, position, comm);
    unpack(allow_increase, buffer, position, comm);
    unpack(damping_factor, buffer, position, comm);
    unpack(free_gas, buffer, position, comm);
    unpack(defaultModel, buffer, position, comm);
    unpack(udaCoefs, buffer, position, comm);
    data = GuideRateModel(timeInterval, target, coefs, allow_increase,
                          damping_factor, free_gas, defaultModel, udaCoefs);
}

void unpack(GuideRateConfig& data,
            std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm)
{
    std::shared_ptr<GuideRateModel> model;
    std::unordered_map<std::string, GuideRateConfig::WellTarget> wells;
    std::unordered_map<std::string, GuideRateConfig::GroupTarget> groups;

    unpack(model, buffer, position, comm);
    unpack(wells, buffer, position, comm);
    unpack(groups, buffer, position, comm);
    data = GuideRateConfig(model, wells, groups);
}

void unpack(GConSale::GCONSALEGroup& data,
            std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm)
{
    unpack(data.sales_target, buffer, position, comm);
    unpack(data.max_sales_rate, buffer, position, comm);
    unpack(data.min_sales_rate, buffer, position, comm);
    unpack(data.max_proc, buffer, position, comm);
    unpack(data.udq_undefined, buffer, position, comm);
    unpack(data.unit_system, buffer, position, comm);
}

void unpack(GConSale& data,
            std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm)
{
    std::map<std::string,GConSale::GCONSALEGroup> groups;
    unpack(groups, buffer, position, comm);
    data = GConSale(groups);
}

void unpack(GConSump::GCONSUMPGroup& data,
            std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm)
{
    unpack(data.consumption_rate, buffer, position, comm);
    unpack(data.import_rate, buffer, position, comm);
    unpack(data.network_node, buffer, position, comm);
    unpack(data.udq_undefined, buffer, position, comm);
    unpack(data.unit_system, buffer, position, comm);
}

void unpack(GConSump& data,
            std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm)
{
    std::map<std::string,GConSump::GCONSUMPGroup> groups;
    unpack(groups, buffer, position, comm);
    data = GConSump(groups);
}

void unpack(RFTConfig& data,
            std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm)
{
    TimeMap timeMap;
    std::size_t first_rft;
    std::pair<bool, std::size_t> wellOpenRftTime;
    RFTConfig::WellOpenTimeMap wellOpenRftName;
    RFTConfig::WellOpenTimeMap wellOpen;
    RFTConfig::RFTMap rftConfig;
    RFTConfig::PLTMap pltConfig;

    unpack(timeMap, buffer, position, comm);
    unpack(first_rft, buffer, position, comm);
    unpack(wellOpenRftTime, buffer, position, comm);
    unpack(wellOpenRftName, buffer, position, comm);
    unpack(wellOpen, buffer, position, comm);
    unpack(rftConfig, buffer, position, comm);
    unpack(pltConfig, buffer, position, comm);
    data = RFTConfig(timeMap, first_rft, wellOpenRftTime, wellOpenRftName,
                     wellOpen, rftConfig, pltConfig);
}


void unpack(DeckItem& data,
            std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm)
{
    std::vector<double> dVal;
    std::vector<int> iVal;
    std::vector<std::string> sVal;
    std::vector<UDAValue> uVal;
    type_tag type;
    std::string name;
    std::vector<value::status> valueStatus;
    bool rawData;
    std::vector<Dimension> activeDimensions, defaultDimensions;

    unpack(dVal, buffer, position, comm);
    unpack(iVal, buffer, position, comm);
    unpack(sVal, buffer, position, comm);
    unpack(uVal, buffer, position, comm);
    unpack(type, buffer, position, comm);
    unpack(name, buffer, position, comm);
    unpack(valueStatus, buffer, position, comm);
    unpack(rawData, buffer, position, comm);
    unpack(activeDimensions, buffer, position, comm);
    unpack(defaultDimensions, buffer, position, comm);
    data = DeckItem(dVal, iVal, sVal, uVal, type, name,
                    valueStatus, rawData, activeDimensions, defaultDimensions);
}

void unpack(DeckRecord& data,
            std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm)
{
    std::vector<DeckItem> items;
    unpack(items, buffer, position, comm);
    data = DeckRecord(std::move(items));
}

void unpack(Location& data,
            std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm)
{
    data.filename.clear();
    unpack(data.filename, buffer, position, comm);
    unpack(data.lineno, buffer, position, comm);
}

void unpack(DeckKeyword& data,
            std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm)
{
    std::string name;
    Location location;
    std::vector<DeckRecord> records;
    bool isDataKeyword, isSlashTerminated;

    unpack(name, buffer, position, comm);
    unpack(location, buffer, position, comm);
    unpack(records, buffer, position, comm);
    unpack(isDataKeyword, buffer, position, comm);
    unpack(isSlashTerminated, buffer, position, comm);
    data = DeckKeyword(name, location, records,
                       isDataKeyword, isSlashTerminated);
}

void unpack(Deck& data, std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm)
{
    std::vector<DeckKeyword> keywords;
    UnitSystem defaultUnitSystem;
    std::unique_ptr<UnitSystem> activeUnitSystem;
    std::string dataFile, inputPath;
    size_t accessCount;

    unpack(keywords, buffer, position, comm);
    unpack(defaultUnitSystem, buffer, position, comm);
    unpack(activeUnitSystem, buffer, position, comm);
    unpack(dataFile, buffer, position, comm);
    unpack(inputPath, buffer, position, comm);
    unpack(accessCount, buffer, position, comm);
    data = Deck(keywords, defaultUnitSystem,
                activeUnitSystem.get(), dataFile, inputPath, accessCount);
}

void unpack(Action::ASTNode& data, std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm)
{
    TokenType token;
    FuncType func_type;
    std::string func;
    std::vector<std::string> argList;
    double number;
    std::vector<Action::ASTNode> children;

    unpack(token, buffer, position, comm);
    unpack(func_type, buffer, position, comm);
    unpack(func, buffer, position, comm);
    unpack(argList, buffer, position, comm);
    unpack(number, buffer, position, comm);
    unpack(children, buffer, position, comm);
    data = Action::ASTNode(token, func_type, func, argList, number, children);
}

void unpack(Action::AST& data, std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm)
{
    std::shared_ptr<Action::ASTNode> condition;
    unpack(condition, buffer, position, comm);
    data = Action::AST(condition);
}

void unpack(Action::Quantity& data, std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm)
{
    unpack(data.quantity, buffer, position, comm);
    unpack(data.args, buffer, position, comm);
}

void unpack(Action::Condition& data, std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm)
{
    unpack(data.lhs, buffer, position, comm);
    unpack(data.rhs, buffer, position, comm);
    unpack(data.logic, buffer, position, comm);
    unpack(data.cmp, buffer, position, comm);
    unpack(data.cmp_string, buffer, position, comm);
}

void unpack(Action::ActionX& data, std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm)
{
    std::string name;
    size_t max_run;
    double min_wait;
    std::time_t start_time;
    std::vector<DeckKeyword> keywords;
    Action::AST condition;
    std::vector<Action::Condition> conditions;
    size_t run_count;
    std::time_t last_run;

    unpack(name, buffer, position, comm);
    unpack(max_run, buffer, position, comm);
    unpack(min_wait, buffer, position, comm);
    unpack(start_time, buffer, position, comm);
    unpack(keywords, buffer, position, comm);
    unpack(condition, buffer, position, comm);
    unpack(conditions, buffer, position, comm);
    unpack(run_count, buffer, position, comm);
    unpack(last_run, buffer, position, comm);
    data = Action::ActionX(name, max_run, min_wait, start_time, keywords,
                           condition, conditions, run_count, last_run);
}

void unpack(Action::Actions& data, std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm)
{
    std::vector<Action::ActionX> actions;
    unpack(actions, buffer, position, comm);
    data = Action::Actions(actions);
}

void unpack(RestartSchedule& data,
            std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm)
{
    unpack(data.timestep, buffer, position, comm);
    unpack(data.basic, buffer, position, comm);
    unpack(data.frequency, buffer, position, comm);
    unpack(data.rptsched_restart_set, buffer, position, comm);
    unpack(data.rptsched_restart, buffer, position, comm);
}

void unpack(TimeStampUTC& data,
            std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm)
{
    TimeStampUTC::YMD ymd;
    int hour, minutes, seconds, usec;

    unpack(ymd, buffer, position, comm);
    unpack(hour, buffer, position, comm);
    unpack(minutes, buffer, position, comm);
    unpack(seconds, buffer, position, comm);
    unpack(usec, buffer, position, comm);
    data = TimeStampUTC(ymd, hour, minutes, seconds, usec);
}

void unpack(WellPolymerProperties& data,
            std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm)
{
    unpack(data.m_polymerConcentration, buffer, position, comm);
    unpack(data.m_saltConcentration, buffer, position, comm);
    unpack(data.m_plymwinjtable, buffer, position, comm);
    unpack(data.m_skprwattable, buffer, position, comm);
    unpack(data.m_skprpolytable, buffer, position, comm);
}

void unpack(Well::WellGuideRate& data,
            std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm)
{
    unpack(data.available, buffer, position, comm);
    unpack(data.guide_rate, buffer, position, comm);
    unpack(data.guide_phase, buffer, position, comm);
    unpack(data.scale_factor, buffer, position, comm);
}

void unpack(GuideRateConfig::WellTarget& data,
            std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm)
{
    unpack(data.guide_rate, buffer, position, comm);
    unpack(data.target, buffer, position, comm);
    unpack(data.scaling_factor, buffer, position, comm);
}

void unpack(GuideRateConfig::GroupTarget& data,
            std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm)
{
    unpack(data.guide_rate, buffer, position, comm);
    unpack(data.target, buffer, position, comm);
}

void unpack(PlyvmhTable& data, std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm)
{
    std::vector<PlyvmhRecord> pdata;
    unpack(pdata, buffer, position, comm);
    data = PlyvmhTable(pdata);
}

#define INSTANTIATE_PACK_VECTOR(...) \
template std::size_t packSize(const std::vector<__VA_ARGS__>& data, \
                              Dune::MPIHelper::MPICommunicator comm); \
template void pack(const std::vector<__VA_ARGS__>& data, \
                   std::vector<char>& buffer, int& position, \
                   Dune::MPIHelper::MPICommunicator comm); \
template void unpack(std::vector<__VA_ARGS__>& data, \
                     std::vector<char>& buffer, int& position, \
                     Dune::MPIHelper::MPICommunicator comm);

INSTANTIATE_PACK_VECTOR(double)
INSTANTIATE_PACK_VECTOR(std::vector<double>)
INSTANTIATE_PACK_VECTOR(bool)
INSTANTIATE_PACK_VECTOR(char)
INSTANTIATE_PACK_VECTOR(int)
INSTANTIATE_PACK_VECTOR(size_t)
INSTANTIATE_PACK_VECTOR(std::array<double, 3>)
INSTANTIATE_PACK_VECTOR(std::pair<bool,double>)
INSTANTIATE_PACK_VECTOR(std::shared_ptr<Group>)
INSTANTIATE_PACK_VECTOR(std::shared_ptr<VFPInjTable>)
INSTANTIATE_PACK_VECTOR(std::shared_ptr<VFPProdTable>)
INSTANTIATE_PACK_VECTOR(std::shared_ptr<Well>)
INSTANTIATE_PACK_VECTOR(std::pair<std::string,std::vector<int>>)
INSTANTIATE_PACK_VECTOR(std::pair<int,std::vector<int>>)

#undef INSTANTIATE_PACK_VECTOR

#define INSTANTIATE_PACK_SET(...) \
template std::size_t packSize(const std::set<__VA_ARGS__>& data, \
                              Dune::MPIHelper::MPICommunicator comm); \
template void pack(const std::set<__VA_ARGS__>& data, \
                   std::vector<char>& buffer, int& position, \
                   Dune::MPIHelper::MPICommunicator comm); \
template void unpack(std::set<__VA_ARGS__>& data, \
                     std::vector<char>& buffer, int& position, \
                     Dune::MPIHelper::MPICommunicator comm);

INSTANTIATE_PACK_SET(std::string)

#undef INSTANTIATE_PACK_SET

#define INSTANTIATE_PACK_SHARED_PTR(...) \
template std::size_t packSize(const std::shared_ptr<__VA_ARGS__>& data, \
                              Dune::MPIHelper::MPICommunicator comm); \
template void pack(const std::shared_ptr<__VA_ARGS__>& data, \
                   std::vector<char>& buffer, int& position, \
                   Dune::MPIHelper::MPICommunicator comm); \
template void unpack(std::shared_ptr<__VA_ARGS__>& data, \
                     std::vector<char>& buffer, int& position, \
                     Dune::MPIHelper::MPICommunicator comm);

INSTANTIATE_PACK_SHARED_PTR(SpiralICD)
INSTANTIATE_PACK_SHARED_PTR(VFPInjTable)
INSTANTIATE_PACK_SHARED_PTR(Well)
INSTANTIATE_PACK_SHARED_PTR(WellTestConfig)
#undef INSTANTIATE_PACK_SHARED_PTR

#define INSTANTIATE_PACK(...) \
template std::size_t packSize(const __VA_ARGS__& data, \
                              Dune::MPIHelper::MPICommunicator comm); \
template void pack(const __VA_ARGS__& data, \
                   std::vector<char>& buffer, int& position, \
                   Dune::MPIHelper::MPICommunicator comm); \
template void unpack(__VA_ARGS__& data, \
                     std::vector<char>& buffer, int& position, \
                     Dune::MPIHelper::MPICommunicator comm);

INSTANTIATE_PACK(double)
INSTANTIATE_PACK(std::size_t)
INSTANTIATE_PACK(bool)
INSTANTIATE_PACK(int)
INSTANTIATE_PACK(std::array<short,3>)
INSTANTIATE_PACK(std::array<bool,3>)
INSTANTIATE_PACK(unsigned char)
INSTANTIATE_PACK(std::map<std::pair<int,int>,std::pair<bool,double>>)
INSTANTIATE_PACK(std::map<FaceDir::DirEnum,std::string>)
INSTANTIATE_PACK(std::map<FaceDir::DirEnum,std::vector<double>>)
INSTANTIATE_PACK(std::map<std::string,Events>)
INSTANTIATE_PACK(std::map<std::string,std::vector<int>>)
INSTANTIATE_PACK(std::map<std::string,std::map<std::pair<int,int>,int>>)
INSTANTIATE_PACK(std::map<UDQVarType,std::size_t>)
INSTANTIATE_PACK(std::unordered_map<std::string,size_t>)
INSTANTIATE_PACK(std::unordered_map<std::string,std::string>)
INSTANTIATE_PACK(std::pair<bool,double>)
INSTANTIATE_PACK(DynamicState<int>)
INSTANTIATE_PACK(DynamicState<OilVaporizationProperties>)
INSTANTIATE_PACK(DynamicState<std::shared_ptr<Action::Actions>>)
INSTANTIATE_PACK(DynamicState<std::shared_ptr<GConSale>>)
INSTANTIATE_PACK(DynamicState<std::shared_ptr<GConSump>>)
INSTANTIATE_PACK(DynamicState<std::shared_ptr<GuideRateConfig>>)
INSTANTIATE_PACK(DynamicState<Tuning>)
INSTANTIATE_PACK(DynamicState<Well::ProducerCMode>)
INSTANTIATE_PACK(DynamicState<std::shared_ptr<WellTestConfig>>)
INSTANTIATE_PACK(DynamicState<std::shared_ptr<WListManager>>)
INSTANTIATE_PACK(DynamicVector<Deck>)

#undef INSTANTIATE_PACK

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
