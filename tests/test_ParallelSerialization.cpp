/*
  Copyright 2019 SINTEF Digital, Mathematics and Cybernetics.

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

#define BOOST_TEST_MODULE TestParallelSerialization
#define BOOST_TEST_NO_MAIN

#include <boost/test/unit_test.hpp>

#include <memory>
#include <tuple>
#include <utility>

#include <opm/common/OpmLog/KeywordLocation.hpp>
#include <opm/input/eclipse/Deck/Deck.hpp>
#include <opm/input/eclipse/Deck/DeckItem.hpp>
#include <opm/input/eclipse/EclipseState/Aquifer/Aquancon.hpp>
#include <opm/input/eclipse/EclipseState/Aquifer/AquiferCT.hpp>
#include <opm/input/eclipse/EclipseState/Aquifer/AquiferConfig.hpp>
#include <opm/input/eclipse/EclipseState/Aquifer/Aquifetp.hpp>
#include <opm/input/eclipse/EclipseState/EclipseConfig.hpp>
#include <opm/input/eclipse/EclipseState/Runspec.hpp>
#include <opm/input/eclipse/EclipseState/TracerConfig.hpp>
#include <opm/input/eclipse/EclipseState/Grid/FaceDir.hpp>
#include <opm/input/eclipse/EclipseState/Grid/Fault.hpp>
#include <opm/input/eclipse/EclipseState/Grid/FaultCollection.hpp>
#include <opm/input/eclipse/EclipseState/Grid/FaultFace.hpp>
#include <opm/input/eclipse/EclipseState/Grid/MULTREGTScanner.hpp>
#include <opm/input/eclipse/EclipseState/Grid/NNC.hpp>
#include <opm/input/eclipse/EclipseState/Grid/TranCalculator.hpp>
#include <opm/input/eclipse/EclipseState/Grid/TransMult.hpp>
#include <opm/input/eclipse/EclipseState/InitConfig/Equil.hpp>
#include <opm/input/eclipse/EclipseState/InitConfig/FoamConfig.hpp>
#include <opm/input/eclipse/EclipseState/InitConfig/InitConfig.hpp>
#include <opm/input/eclipse/EclipseState/IOConfig/IOConfig.hpp>
#include <opm/input/eclipse/Schedule/RPTConfig.hpp>
#include <opm/input/eclipse/Schedule/RSTConfig.hpp>
#include <opm/input/eclipse/Schedule/Source.hpp>
#include <opm/input/eclipse/Schedule/SummaryState.hpp>
#include <opm/input/eclipse/Schedule/Action/ActionResult.hpp>
#include <opm/input/eclipse/Schedule/Action/State.hpp>
#include <opm/input/eclipse/Schedule/Action/ActionAST.hpp>
#include <opm/input/eclipse/Schedule/Action/PyAction.hpp>
#include <opm/input/eclipse/Schedule/Action/Actions.hpp>
#include <opm/input/eclipse/Schedule/Action/ActionX.hpp>
#include <opm/input/eclipse/Schedule/Action/ASTNode.hpp>
#include <opm/input/eclipse/Schedule/Action/Condition.hpp>
#include <opm/input/eclipse/Schedule/Events.hpp>
#include <opm/input/eclipse/Schedule/GasLiftOpt.hpp>
#include <opm/input/eclipse/Schedule/Group/GroupEconProductionLimits.hpp>
#include <opm/input/eclipse/Schedule/Group/GConSale.hpp>
#include <opm/input/eclipse/Schedule/Group/GConSump.hpp>
#include <opm/input/eclipse/Schedule/Group/Group.hpp>
#include <opm/input/eclipse/Schedule/Group/GuideRateConfig.hpp>
#include <opm/input/eclipse/Schedule/Group/GuideRateModel.hpp>
#include <opm/input/eclipse/Schedule/MessageLimits.hpp>
#include <opm/input/eclipse/Schedule/MSW/icd.hpp>
#include <opm/input/eclipse/Schedule/MSW/AICD.hpp>
#include <opm/input/eclipse/Schedule/MSW/SICD.hpp>
#include <opm/input/eclipse/Schedule/MSW/Valve.hpp>
#include <opm/input/eclipse/Schedule/MSW/WellSegments.hpp>
#include <opm/input/eclipse/Schedule/Network/Balance.hpp>
#include <opm/input/eclipse/Schedule/Network/ExtNetwork.hpp>
#include <opm/input/eclipse/Schedule/Network/Node.hpp>
#include <opm/input/eclipse/Schedule/OilVaporizationProperties.hpp>
#include <opm/input/eclipse/Schedule/RFTConfig.hpp>
#include <opm/input/eclipse/Schedule/Schedule.hpp>
#include <opm/input/eclipse/Schedule/ScheduleTypes.hpp>
#include <opm/input/eclipse/Schedule/Tuning.hpp>
#include <opm/input/eclipse/Schedule/UDQ/UDQActive.hpp>
#include <opm/input/eclipse/Schedule/UDQ/UDQAssign.hpp>
#include <opm/input/eclipse/Schedule/UDQ/UDQASTNode.hpp>
#include <opm/input/eclipse/Schedule/UDQ/UDQConfig.hpp>
#include <opm/input/eclipse/Schedule/UDQ/UDQDefine.hpp>
#include <opm/input/eclipse/Schedule/UDQ/UDQFunction.hpp>
#include <opm/input/eclipse/Schedule/UDQ/UDQFunctionTable.hpp>
#include <opm/input/eclipse/Schedule/UDQ/UDQInput.hpp>
#include <opm/input/eclipse/Schedule/UDQ/UDQState.hpp>
#include <opm/input/eclipse/Schedule/VFPInjTable.hpp>
#include <opm/input/eclipse/Schedule/VFPProdTable.hpp>
#include <opm/input/eclipse/Schedule/Well/Connection.hpp>
#include <opm/input/eclipse/Schedule/Well/NameOrder.hpp>
#include <opm/input/eclipse/Schedule/Well/PAvg.hpp>
#include <opm/input/eclipse/Schedule/Well/Well.hpp>
#include <opm/input/eclipse/Schedule/Well/WellBrineProperties.hpp>
#include <opm/input/eclipse/Schedule/Well/WellConnections.hpp>
#include <opm/input/eclipse/Schedule/Well/WellEconProductionLimits.hpp>
#include <opm/input/eclipse/Schedule/Well/WellFoamProperties.hpp>
#include <opm/input/eclipse/Schedule/Well/WellMICPProperties.hpp>
#include <opm/input/eclipse/Schedule/Well/WellPolymerProperties.hpp>
#include <opm/input/eclipse/Schedule/Well/WellTestConfig.hpp>
#include <opm/input/eclipse/Schedule/Well/WellTestState.hpp>
#include <opm/input/eclipse/Schedule/Well/WellTracerProperties.hpp>
#include <opm/input/eclipse/Schedule/Well/WList.hpp>
#include <opm/input/eclipse/Schedule/Well/WListManager.hpp>
#include <opm/input/eclipse/Schedule/Well/WVFPDP.hpp>
#include <opm/input/eclipse/Schedule/Well/WVFPEXP.hpp>
#include <opm/input/eclipse/Schedule/Well/WDFAC.hpp>
#include <opm/input/eclipse/Schedule/WriteRestartFileEvents.hpp>
#include <opm/input/eclipse/EclipseState/SimulationConfig/BCConfig.hpp>
#include <opm/input/eclipse/EclipseState/SimulationConfig/RockConfig.hpp>
#include <opm/input/eclipse/EclipseState/SimulationConfig/SimulationConfig.hpp>
#include <opm/input/eclipse/EclipseState/SimulationConfig/ThresholdPressure.hpp>
#include <opm/input/eclipse/EclipseState/SummaryConfig/SummaryConfig.hpp>
#include <opm/input/eclipse/EclipseState/Tables/Aqudims.hpp>
#include <opm/input/eclipse/EclipseState/Tables/ColumnSchema.hpp>
#include <opm/input/eclipse/EclipseState/Tables/Eqldims.hpp>
#include <opm/input/eclipse/EclipseState/Tables/FlatTable.hpp>
#include <opm/input/eclipse/EclipseState/Tables/DenT.hpp>
#include <opm/input/eclipse/EclipseState/Tables/JFunc.hpp>
#include <opm/input/eclipse/EclipseState/Tables/PlymwinjTable.hpp>
#include <opm/input/eclipse/EclipseState/Tables/PlyshlogTable.hpp>
#include <opm/input/eclipse/EclipseState/Tables/PvtgTable.hpp>
#include <opm/input/eclipse/EclipseState/Tables/PvtoTable.hpp>
#include <opm/input/eclipse/EclipseState/Tables/Regdims.hpp>
#include <opm/input/eclipse/EclipseState/Tables/Rock2dTable.hpp>
#include <opm/input/eclipse/EclipseState/Tables/Rock2dtrTable.hpp>
#include <opm/input/eclipse/EclipseState/Tables/RocktabTable.hpp>
#include <opm/input/eclipse/EclipseState/Tables/SimpleTable.hpp>
#include <opm/input/eclipse/EclipseState/Tables/SkprpolyTable.hpp>
#include <opm/input/eclipse/EclipseState/Tables/SkprwatTable.hpp>
#include <opm/input/eclipse/EclipseState/Tables/Tabdims.hpp>
#include <opm/input/eclipse/EclipseState/Tables/TableColumn.hpp>
#include <opm/input/eclipse/EclipseState/Tables/TableContainer.hpp>
#include <opm/input/eclipse/EclipseState/Tables/TableManager.hpp>
#include <opm/input/eclipse/EclipseState/Tables/TableSchema.hpp>
#include <opm/output/data/Aquifer.hpp>
#include <opm/output/eclipse/RestartValue.hpp>
#include <ebos/eclmpiserializer.hh>

template<class T>
std::tuple<T,int,int> PackUnpack(T& in)
{
    const auto& comm = Dune::MPIHelper::getCommunication();

    Opm::EclMpiSerializer ser(comm);
    ser.pack(in);
    const size_t pos1 = ser.position();
    T out;
    ser.unpack(out);
    const size_t pos2 = ser.position();

    return std::make_tuple(out, pos1, pos2);
}


#define TEST_FOR_TYPE_NAMED_OBJ(TYPE, NAME, OBJ) \
BOOST_AUTO_TEST_CASE(NAME) \
{ \
    auto val1 = Opm::TYPE::OBJ(); \
    auto val2 = PackUnpack(val1); \
    BOOST_CHECK_MESSAGE(std::get<1>(val2) == std::get<2>(val2), "Packed size differ from unpack size for " #TYPE); \
    BOOST_CHECK_MESSAGE(val1 == std::get<0>(val2), "Deserialized " #TYPE " differ"); \
}

#define TEST_FOR_TYPE_NAMED(TYPE, NAME) \
    TEST_FOR_TYPE_NAMED_OBJ(TYPE, NAME, serializationTestObject)

#define TEST_FOR_TYPE(TYPE) \
    TEST_FOR_TYPE_NAMED(TYPE, TYPE)


TEST_FOR_TYPE(Actdims)
TEST_FOR_TYPE(Aqudims)
TEST_FOR_TYPE(Aquancon)
TEST_FOR_TYPE(AquiferConfig)
TEST_FOR_TYPE(AquiferCT)
TEST_FOR_TYPE(Aquifetp)
TEST_FOR_TYPE(AutoICD)
TEST_FOR_TYPE_NAMED(Action::Actions, Actions)
TEST_FOR_TYPE_NAMED(Action::ActionX, ActionX)
TEST_FOR_TYPE_NAMED(Action::AST, ActionAST)
TEST_FOR_TYPE_NAMED(Action::ASTNode, ActionASTNode)
TEST_FOR_TYPE_NAMED(Action::State, ActionState)
TEST_FOR_TYPE(BCConfig)
TEST_FOR_TYPE(BrineDensityTable)
TEST_FOR_TYPE(ColumnSchema)
TEST_FOR_TYPE(Connection)
TEST_FOR_TYPE_NAMED_OBJ(data::AquiferData, AquiferData_CarterTracy, serializationTestObjectC)
TEST_FOR_TYPE_NAMED_OBJ(data::AquiferData, AquiferData_Fetkovich, serializationTestObjectF)
TEST_FOR_TYPE_NAMED_OBJ(data::AquiferData, AquiferData_Numeric, serializationTestObjectN)
TEST_FOR_TYPE_NAMED(data::CarterTracyData, CarterTracyData)
TEST_FOR_TYPE_NAMED(data::CellData, CellData)
TEST_FOR_TYPE_NAMED(data::Connection, dataConnection)
TEST_FOR_TYPE_NAMED(data::CurrentControl, CurrentControl)
TEST_FOR_TYPE_NAMED(data::FetkovichData, FetkovichData)
TEST_FOR_TYPE_NAMED(data::GroupAndNetworkValues, GroupAndNetworkValues)
TEST_FOR_TYPE_NAMED(data::GroupConstraints, GroupConstraints)
TEST_FOR_TYPE_NAMED(data::GroupData, GroupData)
TEST_FOR_TYPE_NAMED(data::GroupGuideRates, GroupGuideRates)
TEST_FOR_TYPE_NAMED(data::GuideRateValue, GuideRateValue)
TEST_FOR_TYPE_NAMED(data::NodeData, NodeData)
TEST_FOR_TYPE_NAMED(data::NumericAquiferData, NumericAquiferData)
TEST_FOR_TYPE_NAMED(data::Rates, Rates)
TEST_FOR_TYPE_NAMED(data::Segment, dataSegment)
TEST_FOR_TYPE_NAMED(data::SegmentPressures, SegmentPressures)
TEST_FOR_TYPE_NAMED(data::Solution, Solution)
TEST_FOR_TYPE_NAMED(data::Well, dataWell)
TEST_FOR_TYPE_NAMED(data::Wells, Wells)
TEST_FOR_TYPE_NAMED(data::WellBlockAvgPress, dataWBPObject)
TEST_FOR_TYPE_NAMED(data::WellBlockAveragePressures, dataWBPCollection)
TEST_FOR_TYPE(Deck)
TEST_FOR_TYPE(DeckItem)
TEST_FOR_TYPE(DeckKeyword)
TEST_FOR_TYPE(DeckRecord)
TEST_FOR_TYPE(DensityTable)
TEST_FOR_TYPE(DenT)
TEST_FOR_TYPE(Dimension)
TEST_FOR_TYPE(EclHysterConfig)
TEST_FOR_TYPE(EclipseConfig)
TEST_FOR_TYPE(EndpointScaling)
TEST_FOR_TYPE(Eqldims)
TEST_FOR_TYPE(Equil)
TEST_FOR_TYPE(TLMixpar)
TEST_FOR_TYPE(Events)
TEST_FOR_TYPE(Fault)
TEST_FOR_TYPE(FaultCollection)
TEST_FOR_TYPE(FaultFace)
TEST_FOR_TYPE_NAMED(Fieldprops::TranCalculator, TranCalculator)
TEST_FOR_TYPE(FoamConfig)
TEST_FOR_TYPE(FoamData)
TEST_FOR_TYPE(GConSale)
TEST_FOR_TYPE(GConSump)
TEST_FOR_TYPE(GridDims)
TEST_FOR_TYPE(Group)
TEST_FOR_TYPE_NAMED(Group::GroupInjectionProperties, GroupInjectionProperties)
TEST_FOR_TYPE_NAMED(Group::GroupProductionProperties, GroupProductionProperties)
TEST_FOR_TYPE(GuideRateConfig)
TEST_FOR_TYPE(GuideRateModel)
TEST_FOR_TYPE(InitConfig)
TEST_FOR_TYPE(IOConfig)
TEST_FOR_TYPE(JFunc)
TEST_FOR_TYPE(KeywordLocation)
TEST_FOR_TYPE(MessageLimits)
TEST_FOR_TYPE(MULTREGTScanner)
TEST_FOR_TYPE(NNC)
TEST_FOR_TYPE_NAMED(Network::Node, NetworkNode)
TEST_FOR_TYPE(OilVaporizationProperties)
TEST_FOR_TYPE(PAvg)
TEST_FOR_TYPE(Phases)
TEST_FOR_TYPE(PlymwinjTable)
TEST_FOR_TYPE(PlyshlogTable)
TEST_FOR_TYPE(PvcdoTable)
TEST_FOR_TYPE(PvtgTable)
TEST_FOR_TYPE(PvtoTable)
TEST_FOR_TYPE(PvtwsaltTable)
TEST_FOR_TYPE(PvtwTable)
TEST_FOR_TYPE(Regdims)
TEST_FOR_TYPE(RestartKey)
TEST_FOR_TYPE(RestartValue)
TEST_FOR_TYPE(RSTConfig)
TEST_FOR_TYPE(RFTConfig)
TEST_FOR_TYPE(RockConfig)
TEST_FOR_TYPE(RockTable)
TEST_FOR_TYPE(RocktabTable)
TEST_FOR_TYPE(Rock2dtrTable)
TEST_FOR_TYPE(Rock2dTable)
TEST_FOR_TYPE(Runspec)
TEST_FOR_TYPE(Schedule)
TEST_FOR_TYPE(ScheduleDeck)
TEST_FOR_TYPE(Segment)
TEST_FOR_TYPE(SimpleTable)
TEST_FOR_TYPE(SimulationConfig)
TEST_FOR_TYPE(SkprpolyTable)
TEST_FOR_TYPE(SkprwatTable)
TEST_FOR_TYPE(SICD)
TEST_FOR_TYPE(SolventDensityTable)
TEST_FOR_TYPE(Source)
TEST_FOR_TYPE(SummaryConfig)
TEST_FOR_TYPE(SummaryConfigNode)
TEST_FOR_TYPE(SummaryState)
TEST_FOR_TYPE(Tabdims)
TEST_FOR_TYPE(TableColumn)
TEST_FOR_TYPE(TableContainer)
TEST_FOR_TYPE(TableManager)
TEST_FOR_TYPE(TableSchema)
TEST_FOR_TYPE(ThresholdPressure)
TEST_FOR_TYPE(TracerConfig)
TEST_FOR_TYPE(TransMult)
TEST_FOR_TYPE(Tuning)
TEST_FOR_TYPE(UDAValue)
TEST_FOR_TYPE(UDQAssign)
TEST_FOR_TYPE(UDQActive)
TEST_FOR_TYPE(UDQASTNode)
TEST_FOR_TYPE(UDQConfig)
TEST_FOR_TYPE(UDQDefine)
TEST_FOR_TYPE(UDQIndex)
TEST_FOR_TYPE(UDQParams)
TEST_FOR_TYPE(UDQState)
TEST_FOR_TYPE(UnitSystem)
TEST_FOR_TYPE(Valve)
TEST_FOR_TYPE(VFPInjTable)
TEST_FOR_TYPE(VFPProdTable)
TEST_FOR_TYPE(ViscrefTable)
TEST_FOR_TYPE(WatdentTable)
TEST_FOR_TYPE(Well)
TEST_FOR_TYPE(Welldims)
TEST_FOR_TYPE(WellBrineProperties)
TEST_FOR_TYPE(WellConnections)
TEST_FOR_TYPE(WellEconProductionLimits)
TEST_FOR_TYPE(WellFoamProperties)
TEST_FOR_TYPE_NAMED(Well::WellGuideRate, WellGuideRate)
TEST_FOR_TYPE_NAMED(Well::WellInjectionProperties, WellInjectionProperties)
TEST_FOR_TYPE(WellPolymerProperties)
TEST_FOR_TYPE_NAMED(Well::WellProductionProperties, WellProductionProperties)
TEST_FOR_TYPE(WellTracerProperties)
TEST_FOR_TYPE(WellSegmentDims)
TEST_FOR_TYPE(WellSegments)
TEST_FOR_TYPE(WellTestConfig)
TEST_FOR_TYPE(WellTestState)
TEST_FOR_TYPE(WellType)
TEST_FOR_TYPE(WListManager)
TEST_FOR_TYPE(WriteRestartFileEvents)


bool init_unit_test_func()
{
    return true;
}


int main(int argc, char** argv)
{
    Dune::MPIHelper::instance(argc, argv);
    return boost::unit_test::unit_test_main(&init_unit_test_func, argc, argv);
}
