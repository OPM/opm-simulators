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

#include <opm/material/common/Tabulated1DFunction.hpp>
#include <opm/material/common/IntervalTabulated2DFunction.hpp>
#include <opm/material/common/UniformXTabulated2DFunction.hpp>
#include <opm/material/fluidsystems/blackoilpvt/ConstantCompressibilityOilPvt.hpp>
#include <opm/material/fluidsystems/blackoilpvt/ConstantCompressibilityWaterPvt.hpp>
#include <opm/material/fluidsystems/blackoilpvt/DeadOilPvt.hpp>
#include <opm/material/fluidsystems/blackoilpvt/DryGasPvt.hpp>
#include <opm/material/fluidsystems/blackoilpvt/GasPvtMultiplexer.hpp>
#include <opm/material/fluidsystems/blackoilpvt/LiveOilPvt.hpp>
#include <opm/material/fluidsystems/blackoilpvt/OilPvtMultiplexer.hpp>
#include <opm/material/fluidsystems/blackoilpvt/SolventPvt.hpp>
#include <opm/material/fluidsystems/blackoilpvt/WaterPvtMultiplexer.hpp>
#include <opm/material/fluidsystems/blackoilpvt/WetGasPvt.hpp>
#include <opm/output/eclipse/RestartValue.hpp>
#include <opm/output/eclipse/EclipseIO.hpp>
#include <opm/output/eclipse/Summary.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/DynamicState.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/DynamicVector.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/TimeMap.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/Well/WellTestConfig.hpp>
#include <opm/parser/eclipse/EclipseState/Util/OrderedMap.hpp>

#include <dune/common/parallel/mpihelper.hh>

#include <tuple>
#include <vector>
#include <map>
#include <unordered_map>

namespace Opm
{

class Actdims;
class Aqudims;
class ColumnSchema;
class DENSITYRecord;
class DensityTable;
class EclHysterConfig;
class Eqldims;
class EDITNNC;
class EndpointScaling;
class Equil;
class EquilRecord;
class Events;
class FoamConfig;
class FoamData;
class InitConfig;
class IOConfig;
class JFunc;
class MessageLimits;
class MLimits;
class NNC;
struct NNCdata;
class OilVaporizationProperties;
class Phases;
class PlymwinjTable;
class PolyInjTable;
class PVCDORecord;
class PvcdoTable;
class PvtgTable;
class PvtoTable;
class PVTWRecord;
class PvtwTable;
class Regdims;
class RestartConfig;
class RestartSchedule;
class ROCKRecord;
class RockTable;
class Rock2dTable;
class Rock2dtrTable;
class Runspec;
class SimulationConfig;
class SimpleTable;
class SkprpolyTable;
class SkprwatTable;
class Tabdims;
class TableColumn;
class TableContainer;
class TableManager;
class TableSchema;
class ThresholdPressure;
class UDQParams;
class VFPInjTable;
class VFPProdTable;
class VISCREFRecord;
class ViscrefTable;
class WATDENTRecord;
class WatdentTable;
class Welldims;
class WellSegmentDims;

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
                     std::integral_constant<bool, false>);

template<class T>
std::size_t packSize(const T&, Dune::MPIHelper::MPICommunicator comm,
                     std::integral_constant<bool, true>);

template<class T>
std::size_t packSize(const T& data, Dune::MPIHelper::MPICommunicator comm);

template<class T1, class T2>
std::size_t packSize(const std::pair<T1,T2>& data, Dune::MPIHelper::MPICommunicator comm);

template<class T, class A>
std::size_t packSize(const std::vector<T,A>& data, Dune::MPIHelper::MPICommunicator comm);

template<class A>
std::size_t packSize(const std::vector<bool,A>& data, Dune::MPIHelper::MPICommunicator comm);

template<class... Ts>
std::size_t packSize(const std::tuple<Ts...>& data, Dune::MPIHelper::MPICommunicator comm);

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

template<class Scalar>
std::size_t packSize(const Tabulated1DFunction<Scalar>& data, Dune::MPIHelper::MPICommunicator comm);

template<class Scalar>
std::size_t packSize(const IntervalTabulated2DFunction<Scalar>& data, Dune::MPIHelper::MPICommunicator comm);

template<class Scalar>
std::size_t packSize(const UniformXTabulated2DFunction<Scalar>& data, Dune::MPIHelper::MPICommunicator comm);

template<class Scalar>
std::size_t packSize(const SolventPvt<Scalar>& data, Dune::MPIHelper::MPICommunicator comm);

template<class Scalar, bool enableThermal>
std::size_t packSize(const GasPvtMultiplexer<Scalar,enableThermal>& data,
                     Dune::MPIHelper::MPICommunicator comm);

template<class Scalar>
std::size_t packSize(const DryGasPvt<Scalar>& data, Dune::MPIHelper::MPICommunicator comm);

template<class Scalar>
std::size_t packSize(const GasPvtThermal<Scalar>& data, Dune::MPIHelper::MPICommunicator comm);

template<class Scalar>
std::size_t packSize(const WetGasPvt<Scalar>& data, Dune::MPIHelper::MPICommunicator comm);

template<class Scalar, bool enableThermal>
std::size_t packSize(const OilPvtMultiplexer<Scalar,enableThermal>& data,
                     Dune::MPIHelper::MPICommunicator comm);

template<class Scalar>
std::size_t packSize(const ConstantCompressibilityOilPvt<Scalar>& data,
                     Dune::MPIHelper::MPICommunicator comm);

template<class Scalar>
std::size_t packSize(const DeadOilPvt<Scalar>& data,
                     Dune::MPIHelper::MPICommunicator comm);

template<class Scalar>
std::size_t packSize(const LiveOilPvt<Scalar>& data,
                     Dune::MPIHelper::MPICommunicator comm);

template<class Scalar>
std::size_t packSize(const OilPvtThermal<Scalar>& data, Dune::MPIHelper::MPICommunicator comm);

template<class Scalar, bool enableThermal>
std::size_t packSize(const WaterPvtMultiplexer<Scalar,enableThermal>& data,
                     Dune::MPIHelper::MPICommunicator comm);

template<class Scalar>
std::size_t packSize(const ConstantCompressibilityWaterPvt<Scalar>& data,
                     Dune::MPIHelper::MPICommunicator comm);

template<class Scalar>
std::size_t packSize(const WaterPvtThermal<Scalar>& data, Dune::MPIHelper::MPICommunicator comm);

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
          Dune::MPIHelper::MPICommunicator, std::integral_constant<bool, false>);

template<class T>
void pack(const T& data, std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm, std::integral_constant<bool, true>);


template<class T>
void pack(const T& data, std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm);

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

template<class Scalar>
void pack(const Tabulated1DFunction<Scalar>& data, std::vector<char>& buffer,
          int& position, Dune::MPIHelper::MPICommunicator comm);

template<class Scalar>
void pack(const IntervalTabulated2DFunction<Scalar>& data, std::vector<char>& buffer,
          int& position, Dune::MPIHelper::MPICommunicator comm);

template<class Scalar>
void pack(const UniformXTabulated2DFunction<Scalar>& data, std::vector<char>& buffer,
          int& position, Dune::MPIHelper::MPICommunicator comm);

template<class Scalar>
void pack(const SolventPvt<Scalar>& data, std::vector<char>& buffer,
          int& position, Dune::MPIHelper::MPICommunicator comm);

template<class Scalar, bool enableThermal>
void pack(const GasPvtMultiplexer<Scalar,enableThermal>& data,
          const std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm);

template<class Scalar>
void pack(const DryGasPvt<Scalar>& data, std::vector<char>& buffer,
          int& position, Dune::MPIHelper::MPICommunicator comm);

template<class Scalar>
void pack(const GasPvtThermal<Scalar>& data, std::vector<char>& buffer,
          int& position, Dune::MPIHelper::MPICommunicator comm);

template<class Scalar>
void pack(const WetGasPvt<Scalar>& data, std::vector<char>& buffer,
          int& position, Dune::MPIHelper::MPICommunicator comm);

template<class Scalar, bool enableThermal>
void pack(const OilPvtMultiplexer<Scalar,enableThermal>& data,
          const std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm);

template<class Scalar>
void pack(const ConstantCompressibilityOilPvt<Scalar>& data,
          std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm);

template<class Scalar>
void pack(const DeadOilPvt<Scalar>& data,
          std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm);

template<class Scalar>
void pack(const LiveOilPvt<Scalar>& data,
          std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm);

template<class Scalar>
void pack(const OilPvtThermal<Scalar>& data, std::vector<char>& buffer,
          int& position, Dune::MPIHelper::MPICommunicator comm);

template<class Scalar, bool enableThermal>
void pack(const WaterPvtMultiplexer<Scalar,enableThermal>& data,
          const std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm);

template<class Scalar>
void pack(const ConstantCompressibilityWaterPvt<Scalar>& data,
          std::vector<char>& buffer, int& position,
          Dune::MPIHelper::MPICommunicator comm);

template<class Scalar>
void pack(const WaterPvtThermal<Scalar>& data, std::vector<char>& buffer,
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
            Dune::MPIHelper::MPICommunicator, std::integral_constant<bool, false>);

template<class T>
void unpack(T& data, std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm, std::integral_constant<bool, true>);

template<class T>
void unpack(T& data, std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm);

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

template<class Scalar>
void unpack(Tabulated1DFunction<Scalar>& data, std::vector<char>& buffer,
            int& position, Dune::MPIHelper::MPICommunicator comm);

template<class Scalar>
void unpack(IntervalTabulated2DFunction<Scalar>& data, std::vector<char>& buffer,
            int& position, Dune::MPIHelper::MPICommunicator comm);

template<class Scalar>
void unpack(UniformXTabulated2DFunction<Scalar>& data, std::vector<char>& buffer,
            int& position, Dune::MPIHelper::MPICommunicator comm);

template<class Scalar>
void unpack(SolventPvt<Scalar>& data, std::vector<char>& buffer,
            int& position, Dune::MPIHelper::MPICommunicator comm);

template<class Scalar, bool enableThermal>
void unpack(GasPvtMultiplexer<Scalar,enableThermal>& data,
            const std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm);

template<class Scalar>
void unpack(DryGasPvt<Scalar>& data, std::vector<char>& buffer,
            int& position, Dune::MPIHelper::MPICommunicator comm);

template<class Scalar>
void unpack(GasPvtThermal<Scalar>& data, std::vector<char>& buffer,
            int& position, Dune::MPIHelper::MPICommunicator comm);

template<class Scalar>
void unpack(WetGasPvt<Scalar>& data, std::vector<char>& buffer,
            int& position, Dune::MPIHelper::MPICommunicator comm);

template<class Scalar, bool enableThermal>
void unpack(OilPvtMultiplexer<Scalar,enableThermal>& data,
            const std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm);

template<class Scalar>
void unpack(ConstantCompressibilityOilPvt<Scalar>& data, std::vector<char>& buffer,
            int& position, Dune::MPIHelper::MPICommunicator comm);

template<class Scalar>
void unpack(DeadOilPvt<Scalar>& data, std::vector<char>& buffer,
            int& position, Dune::MPIHelper::MPICommunicator comm);

template<class Scalar>
void unpack(LiveOilPvt<Scalar>& data, std::vector<char>& buffer,
            int& position, Dune::MPIHelper::MPICommunicator comm);

template<class Scalar>
void unpack(OilPvtThermal<Scalar>& data, std::vector<char>& buffer,
            int& position, Dune::MPIHelper::MPICommunicator comm);

template<class Scalar, bool enableThermal>
void unpack(WaterPvtMultiplexer<Scalar,enableThermal>& data,
            const std::vector<char>& buffer, int& position,
            Dune::MPIHelper::MPICommunicator comm);

template<class Scalar>
void unpack(WaterPvtThermal<Scalar>& data, std::vector<char>& buffer,
            int& position, Dune::MPIHelper::MPICommunicator comm);

template<class Scalar>
void unpack(ConstantCompressibilityWaterPvt<Scalar>& data, std::vector<char>& buffer,
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

ADD_PACK_PROTOTYPES(Actdims)
ADD_PACK_PROTOTYPES(Aqudims)
ADD_PACK_PROTOTYPES(ColumnSchema)
ADD_PACK_PROTOTYPES(data::CellData)
ADD_PACK_PROTOTYPES(data::Connection)
ADD_PACK_PROTOTYPES(data::Rates)
ADD_PACK_PROTOTYPES(data::Segment)
ADD_PACK_PROTOTYPES(data::Solution)
ADD_PACK_PROTOTYPES(data::Well)
ADD_PACK_PROTOTYPES(data::WellRates)
ADD_PACK_PROTOTYPES(DENSITYRecord)
ADD_PACK_PROTOTYPES(DensityTable)
ADD_PACK_PROTOTYPES(EDITNNC)
ADD_PACK_PROTOTYPES(EndpointScaling)
ADD_PACK_PROTOTYPES(Equil)
ADD_PACK_PROTOTYPES(EquilRecord)
ADD_PACK_PROTOTYPES(Events)
ADD_PACK_PROTOTYPES(FoamConfig)
ADD_PACK_PROTOTYPES(FoamData)
ADD_PACK_PROTOTYPES(EclHysterConfig)
ADD_PACK_PROTOTYPES(Eqldims)
ADD_PACK_PROTOTYPES(InitConfig)
ADD_PACK_PROTOTYPES(IOConfig)
ADD_PACK_PROTOTYPES(JFunc)
ADD_PACK_PROTOTYPES(MessageLimits)
ADD_PACK_PROTOTYPES(MLimits)
ADD_PACK_PROTOTYPES(NNC)
ADD_PACK_PROTOTYPES(NNCdata)
ADD_PACK_PROTOTYPES(OilVaporizationProperties)
ADD_PACK_PROTOTYPES(Phases)
ADD_PACK_PROTOTYPES(PlymwinjTable)
ADD_PACK_PROTOTYPES(PolyInjTable)
ADD_PACK_PROTOTYPES(PVCDORecord)
ADD_PACK_PROTOTYPES(PvcdoTable)
ADD_PACK_PROTOTYPES(PvtgTable)
ADD_PACK_PROTOTYPES(PvtoTable)
ADD_PACK_PROTOTYPES(PVTWRecord)
ADD_PACK_PROTOTYPES(PvtwTable)
ADD_PACK_PROTOTYPES(Regdims)
ADD_PACK_PROTOTYPES(RestartConfig)
ADD_PACK_PROTOTYPES(RestartKey)
ADD_PACK_PROTOTYPES(RestartSchedule)
ADD_PACK_PROTOTYPES(RestartValue)
ADD_PACK_PROTOTYPES(ROCKRecord)
ADD_PACK_PROTOTYPES(RockTable)
ADD_PACK_PROTOTYPES(Rock2dTable)
ADD_PACK_PROTOTYPES(Rock2dtrTable)
ADD_PACK_PROTOTYPES(Runspec)
ADD_PACK_PROTOTYPES(std::string)
ADD_PACK_PROTOTYPES(SimulationConfig)
ADD_PACK_PROTOTYPES(SimpleTable)
ADD_PACK_PROTOTYPES(SkprpolyTable)
ADD_PACK_PROTOTYPES(SkprwatTable)
ADD_PACK_PROTOTYPES(Tabdims)
ADD_PACK_PROTOTYPES(TableColumn)
ADD_PACK_PROTOTYPES(TableContainer)
ADD_PACK_PROTOTYPES(TableManager)
ADD_PACK_PROTOTYPES(TableSchema)
ADD_PACK_PROTOTYPES(ThresholdPressure)
ADD_PACK_PROTOTYPES(TimeMap)
ADD_PACK_PROTOTYPES(TimeMap::StepData)
ADD_PACK_PROTOTYPES(UDQParams)
ADD_PACK_PROTOTYPES(VFPInjTable)
ADD_PACK_PROTOTYPES(VFPProdTable)
ADD_PACK_PROTOTYPES(VISCREFRecord)
ADD_PACK_PROTOTYPES(ViscrefTable)
ADD_PACK_PROTOTYPES(WATDENTRecord)
ADD_PACK_PROTOTYPES(WatdentTable)
ADD_PACK_PROTOTYPES(Welldims)
ADD_PACK_PROTOTYPES(WellSegmentDims)
ADD_PACK_PROTOTYPES(WellTestConfig::WTESTWell)

} // end namespace Mpi
RestartValue loadParallelRestart(const EclipseIO* eclIO, SummaryState& summaryState,
                                 const std::vector<Opm::RestartKey>& solutionKeys,
                                 const std::vector<Opm::RestartKey>& extraKeys,
                                 Dune::CollectiveCommunication<Dune::MPIHelper::MPICommunicator> comm);

} // end namespace Opm
#endif // PARALLEL_RESTART_HPP
