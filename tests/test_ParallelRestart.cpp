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

#define BOOST_TEST_MODULE TestParallelRestart
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
#include <opm/input/eclipse/EclipseState/Grid/TransMult.hpp>
#include <opm/input/eclipse/EclipseState/InitConfig/Equil.hpp>
#include <opm/input/eclipse/EclipseState/InitConfig/FoamConfig.hpp>
#include <opm/input/eclipse/EclipseState/InitConfig/InitConfig.hpp>
#include <opm/input/eclipse/EclipseState/IOConfig/IOConfig.hpp>
#include <opm/input/eclipse/Schedule/RSTConfig.hpp>
#include <opm/input/eclipse/Schedule/Action/ActionResult.hpp>
#include <opm/input/eclipse/Schedule/Action/State.hpp>
#include <opm/input/eclipse/Schedule/Action/ActionAST.hpp>
#include <opm/input/eclipse/Schedule/Action/PyAction.hpp>
#include <opm/input/eclipse/Schedule/Action/Actions.hpp>
#include <opm/input/eclipse/Schedule/Action/ActionX.hpp>
#include <opm/input/eclipse/Schedule/Action/ASTNode.hpp>
#include <opm/input/eclipse/Schedule/Action/Condition.hpp>
#include <opm/input/eclipse/Schedule/Events.hpp>
#include <opm/input/eclipse/Schedule/Group/GConSale.hpp>
#include <opm/input/eclipse/Schedule/Group/Group.hpp>
#include <opm/input/eclipse/Schedule/Group/GuideRateModel.hpp>
#include <opm/input/eclipse/Schedule/MessageLimits.hpp>
#include <opm/input/eclipse/Schedule/MSW/icd.hpp>
#include <opm/input/eclipse/Schedule/MSW/AICD.hpp>
#include <opm/input/eclipse/Schedule/MSW/SICD.hpp>
#include <opm/input/eclipse/Schedule/MSW/Valve.hpp>
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
#include <opm/input/eclipse/Schedule/Well/PAvg.hpp>
#include <opm/input/eclipse/Schedule/Well/Connection.hpp>
#include <opm/input/eclipse/Schedule/Well/WellFoamProperties.hpp>
#include <opm/input/eclipse/Schedule/Well/WellPolymerProperties.hpp>
#include <opm/input/eclipse/Schedule/Well/WellTracerProperties.hpp>
#include <opm/input/eclipse/Schedule/Well/Well.hpp>
#include <opm/input/eclipse/Schedule/Well/WList.hpp>
#include <opm/input/eclipse/Schedule/Well/WListManager.hpp>
#include <opm/input/eclipse/Schedule/WriteRestartFileEvents.hpp>
#include <opm/input/eclipse/Schedule/Well/WellTestConfig.hpp>
#include <opm/input/eclipse/Schedule/Well/WellTestState.hpp>
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
#include <opm/simulators/utils/ParallelRestart.hpp>
#include <ebos/eclmpiserializer.hh>


namespace {


Opm::data::Solution getSolution()
{
    Opm::data::Solution sol1;
    sol1.insert("testdata", Opm::UnitSystem::measure::length,
                {1.0, 2.0, 3.0}, Opm::data::TargetType::RESTART_SOLUTION);
    return sol1;
}


Opm::data::Rates getRates()
{
    Opm::data::Rates rat1;
    rat1.set(Opm::data::Rates::opt::wat, 1.0);
    rat1.set(Opm::data::Rates::opt::oil, 2.0);
    rat1.set(Opm::data::Rates::opt::gas, 3.0);
    rat1.set(Opm::data::Rates::opt::polymer, 4.0);
    rat1.set(Opm::data::Rates::opt::solvent, 5.0);
    rat1.set(Opm::data::Rates::opt::energy, 6.0);
    rat1.set(Opm::data::Rates::opt::dissolved_gas, 7.0);
    rat1.set(Opm::data::Rates::opt::vaporized_oil, 8.0);
    rat1.set(Opm::data::Rates::opt::reservoir_water, 9.0);
    rat1.set(Opm::data::Rates::opt::reservoir_oil, 10.0);
    rat1.set(Opm::data::Rates::opt::reservoir_gas, 11.0);
    rat1.set(Opm::data::Rates::opt::productivity_index_water, 12.0);
    rat1.set(Opm::data::Rates::opt::productivity_index_oil, 13.0);
    rat1.set(Opm::data::Rates::opt::productivity_index_gas, 14.0);
    rat1.set(Opm::data::Rates::opt::well_potential_water, 15.0);
    rat1.set(Opm::data::Rates::opt::well_potential_oil, 16.0);
    rat1.set(Opm::data::Rates::opt::well_potential_gas, 17.0);

    return rat1;
}


Opm::data::Connection getConnection()
{
    Opm::data::Connection con1;
    con1.rates = getRates();
    con1.index = 1;
    con1.pressure = 2.0;
    con1.reservoir_rate = 3.0;
    con1.cell_pressure = 4.0;
    con1.cell_saturation_water = 5.0;
    con1.cell_saturation_gas = 6.0;
    con1.effective_Kh = 7.0;
    con1.trans_factor = 8.0;

    return con1;
}


Opm::data::Segment getSegment()
{
    Opm::data::Segment seg1;
    seg1.rates = getRates();
    seg1.segNumber = 1;
    const auto pres_idx = Opm::data::SegmentPressures::Value::Pressure;
    seg1.pressures[pres_idx] = 2.0;
    return seg1;
}


Opm::data::CurrentControl getCurrentControl()
{
    Opm::data::CurrentControl curr;
    curr.isProducer = true;
    curr.prod = ::Opm::Well::ProducerCMode::CRAT;
    return curr;
}

Opm::data::GuideRateValue getWellGuideRate()
{
    using Item = Opm::data::GuideRateValue::Item;

    return Opm::data::GuideRateValue{}.set(Item::Oil  , 1.23)
                                      .set(Item::Gas  , 2.34)
                                      .set(Item::Water, 3.45)
                                      .set(Item::ResV , 4.56);
}

Opm::data::Well getWell()
{
    Opm::data::Well well1;
    well1.rates = getRates();
    well1.bhp = 1.0;
    well1.thp = 2.0;
    well1.temperature = 3.0;
    well1.control = 4;
    well1.connections.push_back(getConnection());
    well1.segments.insert({0, getSegment()});
    well1.current_control = getCurrentControl();
    well1.guide_rates = getWellGuideRate();
    return well1;
}

Opm::data::GroupGuideRates getGroupGuideRates()
{
    using Item = Opm::data::GuideRateValue::Item;

    auto gr = Opm::data::GroupGuideRates{};

    gr.production.set(Item::Oil ,   999.888)
                 .set(Item::Gas ,  8888.777)
                 .set(Item::ResV, 12345.678);

    gr.injection.set(Item::Gas  , 9876.543)
                .set(Item::Water, 2345.987);

    return gr;
}

Opm::data::GroupConstraints getGroupConstraints()
{
    using PMode = ::Opm::Group::ProductionCMode;
    using IMode = ::Opm::Group::InjectionCMode;

    return Opm::data::GroupConstraints{}
    .set(PMode::ORAT,           // Production
         IMode::VREP,           // Gas Injection
         IMode::NONE);          // Water Injection
}

Opm::data::GroupData getGroupData()
{
    return Opm::data::GroupData {
        getGroupConstraints(),
        getGroupGuideRates()
    };
}

Opm::data::NodeData getNodeData()
{
    return Opm::data::NodeData {
        123.457
    };
}

Opm::data::AquiferData getFetkovichAquifer(const int aquiferID = 1)
{
    auto aquifer = Opm::data::AquiferData {
        aquiferID, 123.456, 56.78, 9.0e10, 290.0, 2515.5
    };

    auto* aquFet = aquifer.typeData.create<Opm::data::AquiferType::Fetkovich>();

    aquFet->initVolume = 1.23;
    aquFet->prodIndex = 45.67;
    aquFet->timeConstant = 890.123;

    return aquifer;
}

Opm::data::AquiferData getCarterTracyAquifer(const int aquiferID = 5)
{
    auto aquifer = Opm::data::AquiferData {
        aquiferID, 123.456, 56.78, 9.0e10, 290.0, 2515.5
    };

    auto* aquCT = aquifer.typeData.create<Opm::data::AquiferType::CarterTracy>();

    aquCT->timeConstant = 987.65;
    aquCT->influxConstant = 43.21;
    aquCT->waterDensity = 1014.5;
    aquCT->waterViscosity = 0.00318;
    aquCT->dimensionless_time = 42.0;
    aquCT->dimensionless_pressure = 2.34;

    return aquifer;
}

Opm::data::AquiferData getNumericalAquifer(const int aquiferID = 2)
{
    auto aquifer = Opm::data::AquiferData {
        aquiferID, 123.456, 56.78, 9.0e10, 290.0, 2515.5
    };

    auto* aquNum = aquifer.typeData.create<Opm::data::AquiferType::Numerical>();

    aquNum->initPressure.push_back(1.234);
    aquNum->initPressure.push_back(2.345);
    aquNum->initPressure.push_back(3.4);
    aquNum->initPressure.push_back(9.876);

    return aquifer;
}
}


template<class T>
std::tuple<T,int,int> PackUnpack(const T& in)
{
    auto comm = Dune::MPIHelper::getCollectiveCommunication();
    std::size_t packSize = Opm::Mpi::packSize(in, comm);
    std::vector<char> buffer(packSize);
    int pos1 = 0;
    Opm::Mpi::pack(in, buffer, pos1, comm);
    int pos2 = 0;
    T out;
    Opm::Mpi::unpack(out, buffer, pos2, comm);

    return std::make_tuple(out, pos1, pos2);
}

template<class T>
std::tuple<T,int,int> PackUnpack2(T& in)
{
    auto comm = Dune::MPIHelper::getCollectiveCommunication();
    Opm::EclMpiSerializer ser(comm);
    ser.pack(in);
    size_t pos1 = ser.position();
    T out;
    ser.unpack(out);
    size_t pos2 = ser.position();

    return std::make_tuple(out, pos1, pos2);
}


#define DO_CHECKS(TYPE_NAME) \
    BOOST_CHECK_MESSAGE(std::get<1>(val2) == std::get<2>(val2), "Packed size differ from unpack size for " #TYPE_NAME);  \
    BOOST_CHECK_MESSAGE(val1 == std::get<0>(val2), "Deserialized " #TYPE_NAME " differ");


BOOST_AUTO_TEST_CASE(Solution_)
{
    Opm::data::Solution val1 = getSolution();
    auto val2 = PackUnpack(val1);
    DO_CHECKS(data::Solution)
}


BOOST_AUTO_TEST_CASE(Rates_)
{
    Opm::data::Rates val1 = getRates();
    auto val2 = PackUnpack(val1);
    DO_CHECKS(data::Rates)
}

BOOST_AUTO_TEST_CASE(dataFetkovichData)
{
    const auto val1 = getFetkovichAquifer();
    const auto val2 = PackUnpack(val1);

    DO_CHECKS(data::FetkovichData)
}

BOOST_AUTO_TEST_CASE(dataCarterTracyData)
{
    const auto val1 = getCarterTracyAquifer();
    const auto val2 = PackUnpack(val1);

    DO_CHECKS(data::CarterTracyData)
}

BOOST_AUTO_TEST_CASE(dataNumericAquiferData)
{
    const auto val1 = getNumericalAquifer();
    const auto val2 = PackUnpack(val1);

    DO_CHECKS(data::NumericAquiferData)
}

BOOST_AUTO_TEST_CASE(dataAquifers)
{
    const auto val1 = Opm::data::Aquifers {
        { 1, getFetkovichAquifer(1) },
        { 4, getFetkovichAquifer(4) },
        { 5, getCarterTracyAquifer(5) },
        { 9, getNumericalAquifer(9) },
    };

    const auto val2 = PackUnpack(val1);

    DO_CHECKS(data::Aquifers)
}

BOOST_AUTO_TEST_CASE(dataGuideRateValue)
{
    using Item = Opm::data::GuideRateValue::Item;

    const auto val1 = Opm::data::GuideRateValue{}
    .set(Item::Oil ,   999.888)
    .set(Item::Gas ,  8888.777)
    .set(Item::ResV, 12345.678);

    const auto val2 = PackUnpack(val1);

    BOOST_CHECK_MESSAGE(! std::get<0>(val2).has(Item::Water),
                        "Water Must Not Appear From "
                        "Serializing GuideRateValues");

    DO_CHECKS(data::GuideRateValue)
}

BOOST_AUTO_TEST_CASE(dataConnection_)
{
    Opm::data::Connection val1 = getConnection();
    auto val2 = PackUnpack(val1);
    DO_CHECKS(data::Connection)
}


BOOST_AUTO_TEST_CASE(dataCurrentControl)
{
    Opm::data::CurrentControl val1 = getCurrentControl();
    auto val2 = PackUnpack(val1);
    DO_CHECKS(data::CurrentControl)
}


BOOST_AUTO_TEST_CASE(dataSegment_)
{
    Opm::data::Segment val1 = getSegment();
    auto val2 = PackUnpack(val1);
    DO_CHECKS(data::Segment)
}


BOOST_AUTO_TEST_CASE(dataWell_)
{
    Opm::data::Well val1 = getWell();
    auto val2 = PackUnpack(val1);
    DO_CHECKS(data::Well)
}


BOOST_AUTO_TEST_CASE(Wells_)
{
    Opm::data::Wells val1;
    val1.insert({"test_well", getWell()});
    auto val2 = PackUnpack(val1);
    DO_CHECKS(data::Wells)
}

BOOST_AUTO_TEST_CASE(dataGroupConstraints)
{
    const auto val1 = getGroupConstraints();
    const auto val2 = PackUnpack(val1);

    DO_CHECKS(data::GroupConstraints)
}

BOOST_AUTO_TEST_CASE(dataGroupGuideRates)
{
    const auto val1 = getGroupData().guideRates;
    const auto val2 = PackUnpack(val1);

    DO_CHECKS(data::GroupGuideRates)
}

BOOST_AUTO_TEST_CASE(dataGroupData)
{
    const auto val1 = getGroupData();
    const auto val2 = PackUnpack(val1);

    DO_CHECKS(data::GroupData)
}

BOOST_AUTO_TEST_CASE(dataNodeData)
{
    const auto val1 = getNodeData();
    const auto val2 = PackUnpack(val1);

    DO_CHECKS(data::NodeData)
}

BOOST_AUTO_TEST_CASE(CellData_)
{
    Opm::data::CellData val1;
    val1.dim = Opm::UnitSystem::measure::length;
    val1.data = {1.0, 2.0, 3.0};
    val1.target = Opm::data::TargetType::RESTART_SOLUTION;
    auto val2 = PackUnpack(val1);
    DO_CHECKS(data::cellData)
}


BOOST_AUTO_TEST_CASE(RestartKey_)
{
    Opm::RestartKey val1("key", Opm::UnitSystem::measure::length, true);
    auto val2 = PackUnpack(val1);
    DO_CHECKS(RestartKey)
}


BOOST_AUTO_TEST_CASE(RestartValue)
{
    auto wells1 = Opm::data::Wells {{
        { "test_well", getWell() },
    }};
    auto grp_nwrk_1 = Opm::data::GroupAndNetworkValues {
        {                       // .groupData
            { "test_group1", getGroupData() },
        },
        {                       // .nodeData
            { "test_node1", getNodeData() },
        }
    };

    auto aquifers = Opm::data::Aquifers {
        { 2, getCarterTracyAquifer(2) },
        {11, getFetkovichAquifer(11) },
    };

    const auto val1 = Opm::RestartValue {
        getSolution(), std::move(wells1), std::move(grp_nwrk_1), std::move(aquifers)
    };
    const auto val2 = PackUnpack(val1);

    DO_CHECKS(RestartValue)
}

#define TEST_FOR_TYPE_NAMED(TYPE, NAME) \
BOOST_AUTO_TEST_CASE(NAME) \
{ \
    auto val1 = Opm::TYPE::serializeObject(); \
    auto val2 = PackUnpack2(val1); \
    DO_CHECKS(TYPE) \
}

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
TEST_FOR_TYPE(SummaryConfig)
TEST_FOR_TYPE(SummaryConfigNode)
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
