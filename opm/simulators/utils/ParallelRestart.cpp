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
#include <opm/parser/eclipse/EclipseState/Grid/FaceDir.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/MSW/icd.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/MSW/SpiralICD.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/MSW/Valve.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/RFTConfig.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/ScheduleTypes.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/UDQ/UDQFunction.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/UDQ/UDQFunctionTable.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/VFPInjTable.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/VFPProdTable.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/Well/Connection.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/Well/WellConnections.hpp>
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

HANDLE_AS_POD(data::Connection)
HANDLE_AS_POD(data::CurrentControl)
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

template
std::size_t packSize(const std::map<Phase,Group::GroupInjectionProperties>& data,
                     Dune::MPIHelper::MPICommunicator comm);

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

std::size_t packSize(const WellConnections& data,
                     Dune::MPIHelper::MPICommunicator comm)
{
    return packSize(data.ordering(), comm) +
           packSize(data.getHeadI(), comm) +
           packSize(data.getHeadJ(), comm) +
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

void pack(const WellConnections& data,
          std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm)
{
    pack(data.ordering(), buffer, position, comm);
    pack(data.getHeadI(), buffer, position, comm);
    pack(data.getHeadJ(), buffer, position, comm);
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

void unpack(WellConnections& data,
            std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm)
{
    int headI, headJ;
    Connection::Order ordering;
    std::vector<Connection> connections;

    unpack(ordering, buffer, position, comm),
    unpack(headI, buffer, position, comm),
    unpack(headJ, buffer, position, comm),
    unpack(connections, buffer, position, comm),

    data = WellConnections(ordering, headI, headJ, connections);
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
                ref_depth, wtype, units, udq_undefined, status,
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
INSTANTIATE_PACK_VECTOR(std::time_t)
INSTANTIATE_PACK_VECTOR(std::array<double, 3>)
INSTANTIATE_PACK_VECTOR(std::pair<bool,double>)
INSTANTIATE_PACK_VECTOR(std::shared_ptr<Group>)
INSTANTIATE_PACK_VECTOR(std::shared_ptr<VFPInjTable>)
INSTANTIATE_PACK_VECTOR(std::shared_ptr<VFPProdTable>)
INSTANTIATE_PACK_VECTOR(std::shared_ptr<Well>)
INSTANTIATE_PACK_VECTOR(std::map<std::string,int>)
INSTANTIATE_PACK_VECTOR(std::pair<std::string,std::vector<size_t>>)
INSTANTIATE_PACK_VECTOR(std::pair<int,std::vector<int>>)
INSTANTIATE_PACK_VECTOR(std::pair<int,std::vector<size_t>>)
INSTANTIATE_PACK_VECTOR(std::pair<RFTConfig::RFT,std::size_t>)
INSTANTIATE_PACK_VECTOR(std::pair<RFTConfig::PLT,std::size_t>)
INSTANTIATE_PACK_VECTOR(std::string)

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
INSTANTIATE_PACK(std::map<std::string,std::vector<int>>)
INSTANTIATE_PACK(std::map<std::string,std::map<std::pair<int,int>,int>>)
INSTANTIATE_PACK(std::map<std::string,int>)
INSTANTIATE_PACK(std::map<UDQVarType,std::size_t>)
INSTANTIATE_PACK(std::unordered_map<std::string,size_t>)
INSTANTIATE_PACK(std::unordered_map<std::string,std::string>)
INSTANTIATE_PACK(std::unordered_set<std::string>)
INSTANTIATE_PACK(std::pair<bool,double>)
INSTANTIATE_PACK(std::pair<bool,std::size_t>)

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
