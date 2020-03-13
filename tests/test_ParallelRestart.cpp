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

#include <opm/common/OpmLog/Location.hpp>
#include <opm/parser/eclipse/Deck/Deck.hpp>
#include <opm/parser/eclipse/Deck/DeckItem.hpp>
#include <opm/parser/eclipse/EclipseState/Aquancon.hpp>
#include <opm/parser/eclipse/EclipseState/AquiferCT.hpp>
#include <opm/parser/eclipse/EclipseState/Aquifetp.hpp>
#include <opm/parser/eclipse/EclipseState/EclipseConfig.hpp>
#include <opm/parser/eclipse/EclipseState/Runspec.hpp>
#include <opm/parser/eclipse/EclipseState/Edit/EDITNNC.hpp>
#include <opm/parser/eclipse/EclipseState/Grid/FaceDir.hpp>
#include <opm/parser/eclipse/EclipseState/Grid/Fault.hpp>
#include <opm/parser/eclipse/EclipseState/Grid/FaultCollection.hpp>
#include <opm/parser/eclipse/EclipseState/Grid/FaultFace.hpp>
#include <opm/parser/eclipse/EclipseState/Grid/MULTREGTScanner.hpp>
#include <opm/parser/eclipse/EclipseState/Grid/NNC.hpp>
#include <opm/parser/eclipse/EclipseState/Grid/TransMult.hpp>
#include <opm/parser/eclipse/EclipseState/InitConfig/Equil.hpp>
#include <opm/parser/eclipse/EclipseState/InitConfig/FoamConfig.hpp>
#include <opm/parser/eclipse/EclipseState/InitConfig/InitConfig.hpp>
#include <opm/parser/eclipse/EclipseState/IOConfig/IOConfig.hpp>
#include <opm/parser/eclipse/EclipseState/IOConfig/RestartConfig.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/Action/ActionAST.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/Action/Actions.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/Action/ActionX.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/Action/ASTNode.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/Action/Condition.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/Events.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/Group/GConSale.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/Group/Group.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/Group/GuideRateModel.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/MessageLimits.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/MSW/icd.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/MSW/SpiralICD.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/MSW/Valve.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/OilVaporizationProperties.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/RFTConfig.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/Schedule.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/ScheduleTypes.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/TimeMap.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/Tuning.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/UDQ/UDQActive.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/UDQ/UDQAssign.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/UDQ/UDQASTNode.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/UDQ/UDQConfig.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/UDQ/UDQDefine.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/UDQ/UDQFunction.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/UDQ/UDQFunctionTable.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/UDQ/UDQInput.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/VFPInjTable.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/VFPProdTable.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/Well/Connection.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/Well/WellFoamProperties.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/Well/WellPolymerProperties.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/Well/WellTracerProperties.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/Well/Well.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/Well/WList.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/Well/WListManager.hpp>
#include <opm/parser/eclipse/EclipseState/SimulationConfig/BCConfig.hpp>
#include <opm/parser/eclipse/EclipseState/SimulationConfig/RockConfig.hpp>
#include <opm/parser/eclipse/EclipseState/SimulationConfig/SimulationConfig.hpp>
#include <opm/parser/eclipse/EclipseState/SimulationConfig/ThresholdPressure.hpp>
#include <opm/parser/eclipse/EclipseState/SummaryConfig/SummaryConfig.hpp>
#include <opm/parser/eclipse/EclipseState/Tables/Aqudims.hpp>
#include <opm/parser/eclipse/EclipseState/Tables/ColumnSchema.hpp>
#include <opm/parser/eclipse/EclipseState/Tables/Eqldims.hpp>
#include <opm/parser/eclipse/EclipseState/Tables/FlatTable.hpp>
#include <opm/parser/eclipse/EclipseState/Tables/DenT.hpp>
#include <opm/parser/eclipse/EclipseState/Tables/JFunc.hpp>
#include <opm/parser/eclipse/EclipseState/Tables/PlymwinjTable.hpp>
#include <opm/parser/eclipse/EclipseState/Tables/PlyshlogTable.hpp>
#include <opm/parser/eclipse/EclipseState/Tables/PvtgTable.hpp>
#include <opm/parser/eclipse/EclipseState/Tables/PvtoTable.hpp>
#include <opm/parser/eclipse/EclipseState/Tables/Regdims.hpp>
#include <opm/parser/eclipse/EclipseState/Tables/Rock2dTable.hpp>
#include <opm/parser/eclipse/EclipseState/Tables/Rock2dtrTable.hpp>
#include <opm/parser/eclipse/EclipseState/Tables/RocktabTable.hpp>
#include <opm/parser/eclipse/EclipseState/Tables/SimpleTable.hpp>
#include <opm/parser/eclipse/EclipseState/Tables/SkprpolyTable.hpp>
#include <opm/parser/eclipse/EclipseState/Tables/SkprwatTable.hpp>
#include <opm/parser/eclipse/EclipseState/Tables/Tabdims.hpp>
#include <opm/parser/eclipse/EclipseState/Tables/TableColumn.hpp>
#include <opm/parser/eclipse/EclipseState/Tables/TableContainer.hpp>
#include <opm/parser/eclipse/EclipseState/Tables/TableManager.hpp>
#include <opm/parser/eclipse/EclipseState/Tables/TableSchema.hpp>
#include <opm/output/eclipse/RestartValue.hpp>
#include <opm/simulators/utils/ParallelRestart.hpp>
#include <ebos/eclmpiserializer.hh>


namespace {


#if HAVE_MPI
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
    return con1;
}


Opm::data::Segment getSegment()
{
    Opm::data::Segment seg1;
    seg1.rates = getRates();
    seg1.segNumber = 1;
    seg1.pressure = 2.0;
    return seg1;
}


Opm::data::CurrentControl getCurrentControl()
{
    Opm::data::CurrentControl curr;
    curr.isProducer = true;
    curr.prod = ::Opm::Well::ProducerCMode::CRAT;
    return curr;
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
    return well1;
}


Opm::ThresholdPressure getThresholdPressure()
{
    return Opm::ThresholdPressure(false, true, {{true, 1.0}, {false, 2.0}},
                                  {{{1,2},{false,3.0}},{{2,3},{true,4.0}}});
}


Opm::RockConfig getRockConfig()
{
    return Opm::RockConfig(true, {{100, 0.25}, {200, 0.30}}, "ROCKNUM", 10, false, Opm::RockConfig::Hysteresis::HYSTER);
}


Opm::BCConfig getBCConfig()
{
    return Opm::BCConfig({{10,11,12,13,14,15,Opm::BCType::RATE,Opm::FaceDir::XPlus, Opm::BCComponent::GAS, 100.0}});
}


Opm::TableSchema getTableSchema()
{
    Opm::OrderedMap<std::string, Opm::ColumnSchema> data;
    data.insert({"test1", Opm::ColumnSchema("test1", Opm::Table::INCREASING,
                                                     Opm::Table::DEFAULT_LINEAR)});
    data.insert({"test2", Opm::ColumnSchema("test2", Opm::Table::INCREASING, 1.0)});
    return Opm::TableSchema(data);
}


Opm::TableColumn getTableColumn()
{
    return Opm::TableColumn(Opm::ColumnSchema("test1", Opm::Table::INCREASING,
                                              Opm::Table::DEFAULT_LINEAR),
                            "test2", {1.0, 2.0}, {false, true}, 2);
}


Opm::DenT getDenT()
{
    std::vector<Opm::DenT::entry> records;
    records.emplace_back(1,2,3);
    records.emplace_back(4,5,6);
    return Opm::DenT(records);
}

Opm::SimpleTable getSimpleTable()
{
    Opm::OrderedMap<std::string, Opm::TableColumn> data;
    data.insert({"test3", getTableColumn()});
    return Opm::SimpleTable(getTableSchema(), data, true);
}


Opm::EquilRecord getEquilRecord()
{
    return Opm::EquilRecord(1.0, 2.0, 3.0, 4.0, 5.0, 6.0, true, false, 1);
}


Opm::FoamData getFoamData()
{
    return Opm::FoamData(1.0, 2.0, 3.0, true, 4.0);
}


Opm::FoamConfig getFoamConfig()
{
    return Opm::FoamConfig({getFoamData(), getFoamData()},
                           Opm::Phase::GAS,
                           Opm::FoamConfig::MobilityModel::TAB);
}


Opm::TimeMap getTimeMap()
{
    return Opm::TimeMap({123});
}


Opm::RestartConfig getRestartConfig()
{
    Opm::DynamicState<Opm::RestartSchedule> rsched({Opm::RestartSchedule(1, 2, 3)}, 2);
    Opm::DynamicState<std::map<std::string,int>> rkw({{{"test",3}}}, 3);
    Opm::IOConfig io(true, false, true, false, false, true, "test1", true,
                     "test2", true, "test3", false);
    return Opm::RestartConfig(getTimeMap(), 1, true, rsched, rkw, {false, true});
}


Opm::PvtgTable getPvtgTable()
{
    return Opm::PvtgTable(Opm::ColumnSchema("test1", Opm::Table::INCREASING,
                                            Opm::Table::DEFAULT_LINEAR),
                          getTableColumn(),
                          getTableSchema(),
                          getTableSchema(),
                          {getSimpleTable()},
                          getSimpleTable());
}


Opm::PvtoTable getPvtoTable()
{
    return Opm::PvtoTable(Opm::ColumnSchema("test1", Opm::Table::INCREASING,
                                            Opm::Table::DEFAULT_LINEAR),
                          getTableColumn(),
                          getTableSchema(),
                          getTableSchema(),
                          {getSimpleTable()},
                          getSimpleTable());
}


Opm::TableContainer getTableContainer()
{
    Opm::OrderedMap<std::string, Opm::TableColumn> data;
    data.insert({"test3", getTableColumn()});
    Opm::SimpleTable tab1(getTableSchema(), data, true);
    Opm::TableContainer result(2);
    result.addTable(0, std::make_shared<Opm::SimpleTable>(tab1));
    result.addTable(1, std::make_shared<Opm::SimpleTable>(tab1));
    return result;
}


Opm::Well getFullWell()
{
    Opm::UnitSystem unitSystem;
    return Opm::Well("test1", "test2", 1, 2, 3, 4, 5.0,
                     Opm::WellType(Opm::Phase::WATER), Opm::Connection::Order::DEPTH,
                     unitSystem, 6.0, Opm::Well::Status::SHUT,
                     7.0, true, false,
                     Opm::Well::WellGuideRate{true, 1.0, Opm::Well::GuideRateTarget::COMB, 2.0},
                     8.0, 9.0, false,
                     std::make_shared<Opm::WellEconProductionLimits>(),
                     std::make_shared<Opm::WellFoamProperties>(),
                     std::make_shared<Opm::WellPolymerProperties>(),
                     std::make_shared<Opm::WellBrineProperties>(),
                     std::make_shared<Opm::WellTracerProperties>(),
                     std::make_shared<Opm::WellConnections>(),
                     std::make_shared<Opm::Well::WellProductionProperties>(),
                     std::make_shared<Opm::Well::WellInjectionProperties>(),
                     std::make_shared<Opm::WellSegments>());
}


Opm::VFPInjTable getVFPInjTable()
{
    Opm::VFPInjTable::array_type table;
    table.resize(3*2);
    std::iota(table.begin(), table.end(), 1.0);

    return Opm::VFPInjTable(1, 2.0, Opm::VFPInjTable::FLO_WAT, {1.0, 2.0},
                            {3.0, 4.0, 5.0}, table);
}


Opm::VFPProdTable getVFPProdTable()
{
    Opm::VFPProdTable::array_type table;
    table.resize(1*2*3*4*5);
    std::iota(table.begin(), table.end(), 1.0);

    return Opm::VFPProdTable(1, 2.0, Opm::VFPProdTable::FLO_OIL,
                             Opm::VFPProdTable::WFR_WOR,
                             Opm::VFPProdTable::GFR_GLR,
                             Opm::VFPProdTable::ALQ_TGLR,
                             {1.0, 2.0, 3.0, 4.0, 5.0},
                             {1.0},
                             {1.0, 2.0},
                             {1.0, 2.0, 3.0},
                             {1.0, 2.0, 3.0, 4.0}, table);
}


Opm::UDQConfig getUDQConfig()
{
    Opm::UDQParams params(true, 1, 2.0, 3.0, 4.0);
    std::shared_ptr<Opm::UDQASTNode> n0;
    Opm::UDQASTNode n1(Opm::UDQVarType::NONE,
                       Opm::UDQTokenType::error,
                       "test", 1.0, {"test1", "test2"}, n0, n0);
    Opm::UDQDefine def("test", std::make_shared<Opm::UDQASTNode>(n1),
                       Opm::UDQVarType::NONE, "test2");
    Opm::UDQAssign ass("test", Opm::UDQVarType::NONE,
                       {Opm::UDQAssign::AssignRecord{{"test1"}, 1.0},
                        Opm::UDQAssign::AssignRecord{{"test2"}, 2.0}});
    Opm::OrderedMap<std::string, Opm::UDQIndex> omap;
    omap.insert({"test8", Opm::UDQIndex(1, 2, Opm::UDQAction::ASSIGN,
                                        Opm::UDQVarType::WELL_VAR)});
    omap.insert({"test9", Opm::UDQIndex(3, 4, Opm::UDQAction::ASSIGN,
                                        Opm::UDQVarType::WELL_VAR)});
    return Opm::UDQConfig(params,
                          Opm::UDQFunctionTable(params),
                          {{"test1", def}, {"test2", def}},
                          {{"test3", ass}, {"test4", ass}},
                          {{"test5", "test6"}, {"test7", "test8"}},
                          omap,
                          {{Opm::UDQVarType::SCALAR, 5}, {Opm::UDQVarType::WELL_VAR, 6}});
}


Opm::GuideRateModel getGuideRateModel()
{
    return Opm::GuideRateModel(1.0, Opm::GuideRateModel::Target::WAT,
                               {2.0, 3.0, 4.0, 5.0, 6.0, 7.0},
                               true, 8.0, false, false,
                               {Opm::UDAValue(9.0),
                               Opm::UDAValue(10.0),
                               Opm::UDAValue(11.0)});
}


Opm::GuideRateConfig::GroupTarget getGuideRateConfigGroup()
{
    return Opm::GuideRateConfig::GroupTarget{1.0, Opm::Group::GuideRateTarget::COMB};
}


Opm::GuideRateConfig::WellTarget getGuideRateConfigWell()
{
    return Opm::GuideRateConfig::WellTarget{1.0, Opm::Well::GuideRateTarget::COMB, 2.0};
}


Opm::DeckRecord getDeckRecord()
{
    Opm::DeckItem item1({1.0}, {2}, {"test3"}, {Opm::UDAValue(4)},
                       Opm::type_tag::string, "test5",
                       {Opm::value::status::deck_value},
                       true,
                       {Opm::Dimension(7.0, 8.0)},
                       {Opm::Dimension(10.0, 11.0)});

    Opm::DeckItem item2({1.0}, {2}, {"test3"}, {Opm::UDAValue(4)},
                       Opm::type_tag::string, "test6",
                       {Opm::value::status::deck_value},
                       true,
                       {Opm::Dimension(7.0, 8.0)},
                       {Opm::Dimension(10.0, 11.0)});

    return Opm::DeckRecord({item1, item2});
}


Opm::Tuning getTuning()
{
    return Opm::Tuning{1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, true,
                       11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0, 19.0,
                       20.0, 21.0, false, 22.0, 3, 4, 5, 6, 7, 8, 9, 23.0, 24.0,
                       25.0, 26.0, true};
}


Opm::Action::Condition getCondition()
{
    Opm::Action::Quantity q;
    q.quantity = "test1";
    q.args = {"test2", "test3"};
    Opm::Action::Condition val1;
    val1.lhs = val1.rhs = q;
    val1.logic = Opm::Action::Condition::Logical::OR;
    val1.cmp = Opm::Action::Condition::Comparator::LESS;
    val1.cmp_string = "test";
    return val1;
}


Opm::Action::ActionX getActionX()
{
    std::shared_ptr<Opm::Action::ASTNode> node;
    node.reset(new Opm::Action::ASTNode(number, FuncType::field,
                                        "test1", {"test2"}, 1.0, {}));
    Opm::Action::AST ast(node);
    return Opm::Action::ActionX("test", 1, 2.0, 3,
                                {Opm::DeckKeyword("test", {"test",1},
                                                  {getDeckRecord(), getDeckRecord()},
                                                  true, false)},
                                ast, {getCondition()}, 4, 5);
}


Opm::AquiferCT getAquiferCT() {
    Opm::AquiferCT::AQUCT_data data;
    data.aquiferID = 1;
    data.inftableID = 2;
    data.pvttableID = 3;
    data.phi_aq = 100;
    data.d0 = 1;
    data.C_t = 10;
    data.r_o = 1.5;
    data.k_a = 100;
    data.c1 = 0.78;
    data.h = 1;
    data.c2 = 45;
    data.p0 = std::make_pair(true, 98);
    data.td = {1,2,3};
    data.pi = {4,5,6};
    data.cell_id = {0,10,100};

    return Opm::AquiferCT( { data } );
}

Opm::Aquifetp getAquifetp() {
    Opm::Aquifetp::AQUFETP_data data;

    data.aquiferID = 1;
    data.pvttableID = 3;
    data.C_t = 10;
    data.p0 = std::make_pair(true, 98);
    data.V0 = 0;
    data.d0 = 0;

    return Opm::Aquifetp( { data } );
}



Opm::Aquancon getAquancon() {
    Opm::Aquancon::AquancCell cell(1, 100, std::make_pair(false, 0), 100, Opm::FaceDir::XPlus);
    return Opm::Aquancon( std::unordered_map<int, std::vector<Opm::Aquancon::AquancCell>>{{1, {cell}}});
}

#endif


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


BOOST_AUTO_TEST_CASE(Solution)
{
#if HAVE_MPI
    Opm::data::Solution val1 = getSolution();
    auto val2 = PackUnpack(val1);
    DO_CHECKS(data::Solution)
#endif
}


BOOST_AUTO_TEST_CASE(Rates)
{
#if HAVE_MPI
    Opm::data::Rates val1 = getRates();
    auto val2 = PackUnpack(val1);
    DO_CHECKS(data::Rates)
#endif
}


BOOST_AUTO_TEST_CASE(dataConnection)
{
#if HAVE_MPI
    Opm::data::Connection val1 = getConnection();
    auto val2 = PackUnpack(val1);
    DO_CHECKS(data::Connection)
#endif
}


BOOST_AUTO_TEST_CASE(dataCurrentControl)
{
#if HAVE_MPI
    Opm::data::CurrentControl cur1 = getCurrentControl();
    auto cur2 = PackUnpack(cur1);
    BOOST_CHECK(std::get<1>(cur2) == std::get<2>(cur2));
    BOOST_CHECK(cur1 == std::get<0>(cur2));
#endif
}


BOOST_AUTO_TEST_CASE(dataSegment)
{
#if HAVE_MPI
    Opm::data::Segment val1 = getSegment();
    auto val2 = PackUnpack(val1);
    DO_CHECKS(data::Segment)
#endif
}


BOOST_AUTO_TEST_CASE(dataWell)
{
#if HAVE_MPI
    Opm::data::Well val1 = getWell();
    auto val2 = PackUnpack(val1);
    DO_CHECKS(data::Well)
#endif
}


BOOST_AUTO_TEST_CASE(WellRates)
{
#if HAVE_MPI
    Opm::data::WellRates val1;
    val1.insert({"test_well", getWell()});
    auto val2 = PackUnpack(val1);
    DO_CHECKS(data::WellRates)
#endif
}


BOOST_AUTO_TEST_CASE(CellData)
{
#if HAVE_MPI
    Opm::data::CellData val1;
    val1.dim = Opm::UnitSystem::measure::length;
    val1.data = {1.0, 2.0, 3.0};
    val1.target = Opm::data::TargetType::RESTART_SOLUTION;
    auto val2 = PackUnpack(val1);
    DO_CHECKS(data::cellData)
#endif
}


BOOST_AUTO_TEST_CASE(RestartKey)
{
#if HAVE_MPI
    Opm::RestartKey val1("key", Opm::UnitSystem::measure::length, true);
    auto val2 = PackUnpack(val1);
    DO_CHECKS(RestartKey)
#endif
}


BOOST_AUTO_TEST_CASE(RestartValue)
{
#if HAVE_MPI
    Opm::data::WellRates wells1;
    wells1.insert({"test_well", getWell()});
    Opm::RestartValue val1(getSolution(), wells1);
    auto val2 = PackUnpack(val1);
    DO_CHECKS(RestartValue)
#endif
}


BOOST_AUTO_TEST_CASE(ThresholdPressure)
{
#if HAVE_MPI
    Opm::ThresholdPressure val1 = getThresholdPressure();
    auto val2 = PackUnpack2(val1);
    DO_CHECKS(ThresholdPressure)
#endif
}

BOOST_AUTO_TEST_CASE(RockConfig)
{
#if HAVE_MPI
    Opm::RockConfig val1 = getRockConfig();
    auto val2 = PackUnpack2(val1);
    DO_CHECKS(RockConfig)
#endif
}


BOOST_AUTO_TEST_CASE(EDITNNC)
{
#if HAVE_MPI
    Opm::EDITNNC val1({{1,2,1.0},{2,3,2.0}});
    auto val2 = PackUnpack2(val1);
    DO_CHECKS(EDITNNC)
#endif
}


BOOST_AUTO_TEST_CASE(NNC)
{
#if HAVE_MPI
    Opm::NNC val1({{1,2,1.0},{2,3,2.0}});
    auto val2 = PackUnpack2(val1);
    DO_CHECKS(NNC)
#endif
}


BOOST_AUTO_TEST_CASE(Rock2dTable)
{
#if HAVE_MPI
    Opm::Rock2dTable val1({{1.0,2.0},{3.0,4.0}}, {1.0, 2.0, 3.0});
    auto val2 = PackUnpack(val1);
    DO_CHECKS(Rock2dTable)
#endif
}


BOOST_AUTO_TEST_CASE(Rock2dtrTable)
{
#if HAVE_MPI
    Opm::Rock2dtrTable val1({{1.0,2.0},{3.0,4.0}}, {1.0, 2.0, 3.0});
    auto val2 = PackUnpack(val1);
    DO_CHECKS(Rock2dtrTable)
#endif
}


BOOST_AUTO_TEST_CASE(ColumnSchema)
{
#if HAVE_MPI
    Opm::ColumnSchema val1("test1", Opm::Table::INCREASING,
                           Opm::Table::DEFAULT_LINEAR);
    auto val2 = PackUnpack(val1);
    DO_CHECKS(ColumnSchema)
    val1 = Opm::ColumnSchema("test2", Opm::Table::DECREASING, 1.0);
    val2 = PackUnpack(val1);
    DO_CHECKS(ColumnSchema)
#endif
}


BOOST_AUTO_TEST_CASE(TableSchema)
{
#if HAVE_MPI
    Opm::TableSchema val1 = getTableSchema();
    auto val2 = PackUnpack(val1);
    DO_CHECKS(TableSchema)
#endif
}


BOOST_AUTO_TEST_CASE(TableColumn)
{
#if HAVE_MPI
    Opm::TableColumn val1 = getTableColumn();
    auto val2 = PackUnpack(val1);
    DO_CHECKS(TableColumn)
#endif
}


BOOST_AUTO_TEST_CASE(SimpleTable)
{
#if HAVE_MPI
    Opm::SimpleTable val1 = getSimpleTable();
    auto val2 = PackUnpack(val1);
    DO_CHECKS(SimpleTable)
#endif
}


BOOST_AUTO_TEST_CASE(TableContainer)
{
#if HAVE_MPI
    Opm::OrderedMap<std::string, Opm::TableColumn> data;
    data.insert({"test3", getTableColumn()});
    Opm::SimpleTable tab1(getTableSchema(), data, true);
    Opm::TableContainer val1(2);
    val1.addTable(0, std::make_shared<Opm::SimpleTable>(tab1));
    val1.addTable(1, std::make_shared<Opm::SimpleTable>(tab1));
    auto val2 = PackUnpack(val1);
    DO_CHECKS(TableContainer)
#endif
}


BOOST_AUTO_TEST_CASE(EquilRecord)
{
#if HAVE_MPI
    Opm::EquilRecord val1 = getEquilRecord();
    auto val2 = PackUnpack2(val1);
    DO_CHECKS(EquilRecord)
#endif
}


BOOST_AUTO_TEST_CASE(Equil)
{
#if HAVE_MPI
    Opm::Equil val1({getEquilRecord(), getEquilRecord()});
    auto val2 = PackUnpack2(val1);
    DO_CHECKS(Equil)
#endif
}


BOOST_AUTO_TEST_CASE(FoamData)
{
#if HAVE_MPI
    Opm::FoamData val1 = getFoamData();
    auto val2 = PackUnpack2(val1);
    DO_CHECKS(FoamData)
#endif
}


BOOST_AUTO_TEST_CASE(FoamConfig)
{
#if HAVE_MPI
    Opm::FoamConfig val1 = getFoamConfig();
    auto val2 = PackUnpack2(val1);
    DO_CHECKS(FoamConfig)
#endif
}


BOOST_AUTO_TEST_CASE(InitConfig)
{
#if HAVE_MPI
    Opm::InitConfig val1(Opm::Equil({getEquilRecord(), getEquilRecord()}),
                         getFoamConfig(),
                         true, true, true, 20, "test1");
    auto val2 = PackUnpack2(val1);
    DO_CHECKS(InitConfig)
#endif
}


BOOST_AUTO_TEST_CASE(SimulationConfig)
{
#if HAVE_MPI
    Opm::SimulationConfig val1(getThresholdPressure(), getBCConfig(), getRockConfig(), false, true, false, true);
    auto val2 = PackUnpack2(val1);
    DO_CHECKS(SimulationConfig)
#endif
}


BOOST_AUTO_TEST_CASE(BCConfig)
{
#if HAVE_MPI
    Opm::BCConfig val1({{10,11,12,13,14,15,Opm::BCType::RATE, Opm::FaceDir::XPlus, Opm::BCComponent::GAS, 100}});
    auto val2 = PackUnpack2(val1);
    DO_CHECKS(BCConfig)
#endif
}


BOOST_AUTO_TEST_CASE(RestartSchedule)
{
#if HAVE_MPI
    Opm::RestartSchedule val1(1, 2, 3);
    auto val2 = PackUnpack(val1);
    DO_CHECKS(RestartSchedule)
#endif
}



BOOST_AUTO_TEST_CASE(TimeMap)
{
#if HAVE_MPI
    Opm::TimeMap val1 = getTimeMap();
    auto val2 = PackUnpack(val1);
    DO_CHECKS(TimeMap)
#endif
}


BOOST_AUTO_TEST_CASE(RestartConfig)
{
#if HAVE_MPI
    Opm::DynamicState<Opm::RestartSchedule> rsched({Opm::RestartSchedule(1, 2, 3)}, 2);
    Opm::DynamicState<std::map<std::string,int>> rkw({{{"test",3}}}, 3);
    Opm::IOConfig io(true, false, true, false, false, true, "test1", true,
                     "test2", true, "test3", false);
    Opm::RestartConfig val1(getTimeMap(), 1, true, rsched, rkw, {false, true});
    auto val2 = PackUnpack(val1);
    DO_CHECKS(RestartConfig)
#endif
}


BOOST_AUTO_TEST_CASE(IOConfig)
{
#if HAVE_MPI
    Opm::IOConfig val1(true, false, true, false, false, true, "test1", true,
                       "test2", true, "test3", false);
    auto val2 = PackUnpack2(val1);
    DO_CHECKS(IOConfig)
#endif
}


BOOST_AUTO_TEST_CASE(Phases)
{
#if HAVE_MPI
    Opm::Phases val1(true, true, true, false, true, false, true, false);
    auto val2 = PackUnpack(val1);
    DO_CHECKS(Phases)
#endif
}


BOOST_AUTO_TEST_CASE(Tabdims)
{
#if HAVE_MPI
    Opm::Tabdims val1(1,2,3,4,5,6);
    auto val2 = PackUnpack(val1);
    DO_CHECKS(Tabdims)
#endif
}


BOOST_AUTO_TEST_CASE(EndpointScaling)
{
#if HAVE_MPI
    Opm::EndpointScaling val1(std::bitset<4>(13));
    auto val2 = PackUnpack(val1);
    DO_CHECKS(EndpointScaling)
#endif
}


BOOST_AUTO_TEST_CASE(Welldims)
{
#if HAVE_MPI
    Opm::Welldims val1(1,2,3,4);
    auto val2 = PackUnpack(val1);
    DO_CHECKS(Welldims)
#endif
}


BOOST_AUTO_TEST_CASE(WellSegmentDims)
{
#if HAVE_MPI
    Opm::WellSegmentDims val1(1,2,3);
    auto val2 = PackUnpack(val1);
    DO_CHECKS(WellSegmentDims)
#endif
}


BOOST_AUTO_TEST_CASE(UDQParams)
{
#if HAVE_MPI
    Opm::UDQParams val1(true, 1, 2.0, 3.0, 4.0);
    auto val2 = PackUnpack(val1);
    DO_CHECKS(UDQParams)
#endif
}


BOOST_AUTO_TEST_CASE(EclHysterConfig)
{
#if HAVE_MPI
    Opm::EclHysterConfig val1(true, 1, 2);
    auto val2 = PackUnpack(val1);
    DO_CHECKS(EclHysterConfig)
#endif
}


BOOST_AUTO_TEST_CASE(Actdims)
{
#if HAVE_MPI
    Opm::Actdims val1(1,2,3,4);
    auto val2 = PackUnpack(val1);
    DO_CHECKS(Actdims)
#endif
}


BOOST_AUTO_TEST_CASE(Runspec)
{
#if HAVE_MPI
    Opm::Runspec val1(Opm::Phases(true, true, true, false, true, false, true, false),
                      Opm::Tabdims(1,2,3,4,5,6),
                      Opm::EndpointScaling(std::bitset<4>(13)),
                      Opm::Welldims(1,2,3,4),
                      Opm::WellSegmentDims(1,2,3),
                      Opm::UDQParams(true, 1, 2.0, 3.0, 4.0),
                      Opm::EclHysterConfig(true, 1, 2),
                      Opm::Actdims(1,2,3,4),
                      Opm::SatFuncControls(5.0e-7, Opm::SatFuncControls::ThreePhaseOilKrModel::Stone2));

    auto val2 = PackUnpack(val1);
    DO_CHECKS(Runspec)
#endif
}


BOOST_AUTO_TEST_CASE(PvtgTable)
{
#if HAVE_MPI
    Opm::PvtgTable val1 = getPvtgTable();
    auto val2 = PackUnpack(val1);
    DO_CHECKS(PvtgTable)
#endif
}


BOOST_AUTO_TEST_CASE(PvtoTable)
{
#if HAVE_MPI
    Opm::PvtoTable val1 = getPvtoTable();
    auto val2 = PackUnpack(val1);
    DO_CHECKS(PvtoTable)
#endif
}


BOOST_AUTO_TEST_CASE(JFunc)
{
#if HAVE_MPI
    Opm::JFunc val1(Opm::JFunc::Flag::BOTH, 1.0, 2.0,
                    3.0, 4.0, Opm::JFunc::Direction::XY);
    auto val2 = PackUnpack(val1);
    DO_CHECKS(JFunc)
#endif
}


BOOST_AUTO_TEST_CASE(PVTWRecord)
{
#if HAVE_MPI
    Opm::PVTWRecord val1{1.0, 2.0, 3.0, 4.0, 5.0};
    auto val2 = PackUnpack(val1);
    DO_CHECKS(PVTWRecord)
#endif
}


BOOST_AUTO_TEST_CASE(PvtwTable)
{
#if HAVE_MPI
    Opm::PvtwTable val1({Opm::PVTWRecord{1.0, 2.0, 3.0, 4.0, 5.0}});
    auto val2 = PackUnpack(val1);
    DO_CHECKS(PvtwTable)
#endif
}


BOOST_AUTO_TEST_CASE(PVCDORecord)
{
#if HAVE_MPI
    Opm::PVTWRecord val1{1.0, 2.0, 3.0, 4.0, 5.0};
    auto val2 = PackUnpack(val1);
    DO_CHECKS(PVTWRecord)
#endif
}


BOOST_AUTO_TEST_CASE(PvcdoTable)
{
#if HAVE_MPI
    Opm::PvcdoTable val1({Opm::PVCDORecord{1.0, 2.0, 3.0, 4.0, 5.0}});
    auto val2 = PackUnpack(val1);
    DO_CHECKS(PvcdoTable)
#endif
}


BOOST_AUTO_TEST_CASE(DENSITYRecord)
{
#if HAVE_MPI
    Opm::DENSITYRecord val1{1.0, 2.0, 3.0};
    auto val2 = PackUnpack(val1);
    DO_CHECKS(DENSITYRecord)
#endif
}


BOOST_AUTO_TEST_CASE(DensityTable)
{
#if HAVE_MPI
    Opm::DensityTable val1({Opm::DENSITYRecord{1.0, 2.0, 3.0}});
    auto val2 = PackUnpack(val1);
    DO_CHECKS(DensityTable)
#endif
}


BOOST_AUTO_TEST_CASE(VISCREFRecord)
{
#if HAVE_MPI
    Opm::VISCREFRecord val1{1.0, 2.0};
    auto val2 = PackUnpack(val1);
    DO_CHECKS(VISCREFRecord)
#endif
}


BOOST_AUTO_TEST_CASE(ViscrefTable)
{
#if HAVE_MPI
    Opm::ViscrefTable val1({Opm::VISCREFRecord{1.0, 2.0}});
    auto val2 = PackUnpack(val1);
    DO_CHECKS(ViscrefTable)
#endif
}


BOOST_AUTO_TEST_CASE(WATDENTRecord)
{
#if HAVE_MPI
    Opm::WATDENTRecord val1{1.0, 2.0, 3.0};
    auto val2 = PackUnpack(val1);
    DO_CHECKS(WATDENTRecord)
#endif
}


BOOST_AUTO_TEST_CASE(WatdentTable)
{
#if HAVE_MPI
    Opm::WatdentTable val1({Opm::WATDENTRecord{1.0, 2.0, 3.0}});
    auto val2 = PackUnpack(val1);
    DO_CHECKS(WatdentTable)
#endif
}


BOOST_AUTO_TEST_CASE(PlymwinjTable)
{
#if HAVE_MPI
    Opm::PlymwinjTable val1({1.0}, {2.0}, 1, {{1.0}, {2.0}});
    auto val2 = PackUnpack(val1);
    DO_CHECKS(PlymwinjTable)
#endif
}


BOOST_AUTO_TEST_CASE(SkprpolyTable)
{
#if HAVE_MPI
    Opm::SkprpolyTable val1({1.0}, {2.0}, 1, {{1.0}, {2.0}}, 3.0);
    auto val2 = PackUnpack(val1);
    DO_CHECKS(SkprpolyTable)
#endif
}


BOOST_AUTO_TEST_CASE(SkprwatTable)
{
#if HAVE_MPI
    Opm::SkprwatTable val1({1.0}, {2.0}, 1, {{1.0}, {2.0}});
    auto val2 = PackUnpack(val1);
    DO_CHECKS(SkprwatTable)
#endif
}


BOOST_AUTO_TEST_CASE(Regdims)
{
#if HAVE_MPI
    Opm::Regdims val1(1,2,3,4,5);
    auto val2 = PackUnpack(val1);
    DO_CHECKS(Regdims)
#endif
}


BOOST_AUTO_TEST_CASE(Eqldims)
{
#if HAVE_MPI
    Opm::Eqldims val1(1,2,3,4,5);
    auto val2 = PackUnpack(val1);
    DO_CHECKS(Eqldims)
#endif
}


BOOST_AUTO_TEST_CASE(Aqudims)
{
#if HAVE_MPI
    Opm::Aqudims val1(1,2,3,4,5,6,7,8);
    auto val2 = PackUnpack(val1);
    DO_CHECKS(Aqudims)
#endif
}


BOOST_AUTO_TEST_CASE(ROCKRecord)
{
#if HAVE_MPI
    Opm::ROCKRecord val1{1.0,2.0};
    auto val2 = PackUnpack(val1);
    DO_CHECKS(ROCKRecord)
#endif
}


BOOST_AUTO_TEST_CASE(RockTable)
{
#if HAVE_MPI
    Opm::RockTable val1({Opm::ROCKRecord{1.0,2.0}});
    auto val2 = PackUnpack(val1);
    DO_CHECKS(RockTable)
#endif
}


BOOST_AUTO_TEST_CASE(TableManager)
{
#if HAVE_MPI
    auto jfunc = std::make_shared<Opm::JFunc>(Opm::JFunc::Flag::BOTH,
                                              1.0, 2.0, 3.0, 4.0,
                                              Opm::JFunc::Direction::XY);
    Opm::TableManager val1({{"test", getTableContainer()}},
                           {getPvtgTable()},
                           {getPvtoTable()},
                           {Opm::Rock2dTable({{1.0,2.0},{3.0,4.0}}, {1.0, 2.0, 3.0})},
                           {Opm::Rock2dtrTable({{1.0,2.0},{3.0,4.0}}, {1.0, 2.0, 3.0})},
                           Opm::PvtwTable({Opm::PVTWRecord{1.0, 2.0, 3.0, 4.0, 5.0}}),
                           Opm::PvcdoTable({Opm::PVCDORecord{1.0, 2.0, 3.0, 4.0, 5.0}}),
                           Opm::DensityTable({Opm::DENSITYRecord{1.0, 2.0, 3.0}}),
                           Opm::PlyvmhTable({Opm::PlyvmhRecord{1.0, 2.0, 3.0, 4.0}}),
                           Opm::RockTable({Opm::ROCKRecord{1.0,2.0}}),
                           Opm::PlmixparTable({Opm::PlmixparRecord{1.0}}),
                           Opm::ShrateTable({Opm::ShrateRecord{1.0}}),
                           Opm::Stone1exTable({Opm::Stone1exRecord{1.0}}),
                           Opm::TlmixparTable({Opm::TlmixparRecord{1.0, 2.0}}),
                           Opm::ViscrefTable({Opm::VISCREFRecord{1.0, 2.0}}),
                           Opm::WatdentTable({Opm::WATDENTRecord{1.0, 2.0, 3.0}}),
                           {{1.0, 2.0, {1.0, 2.0, 3.0}}},
                           {{{1.0, 2.0, 3.0}}},
                           {{{4.0, 5.0, 6.0}}},
                           {{1, Opm::PlymwinjTable({1.0}, {2.0}, 1, {{1.0}, {2.0}})}},
                           {{2, Opm::SkprwatTable({1.0}, {2.0}, 1, {{1.0}, {2.0}})}},
                           {{3, Opm::SkprpolyTable({1.0}, {2.0}, 1, {{1.0}, {2.0}}, 3.0)}},
                           Opm::Tabdims(1,2,3,4,5,6),
                           Opm::Regdims(1,2,3,4,5),
                           Opm::Eqldims(1,2,3,4,5),
                           Opm::Aqudims(1,2,3,4,5,6,7,8),
                           true,
                           true,
                           true,
                           true,
                           jfunc,
                           getDenT(),
                           getDenT(),
                           getDenT(),
                           {7.0, 8.0},
                           77,
                           1.0);
    auto val2 = PackUnpack(val1);
    DO_CHECKS(TableManager)
#endif
}


BOOST_AUTO_TEST_CASE(OilVaporizationProperties)
{
#ifdef HAVE_MPI
    using VapType = Opm::OilVaporizationProperties::OilVaporization;
    Opm::OilVaporizationProperties val1(VapType::VAPPARS,
                                        1.0, 2.0, {5.0, 6.0},
                                        {false, true}, {7.0, 8.0});
    auto val2 = PackUnpack(val1);
    DO_CHECKS(OilVaporizationProperties)
    val1 = Opm::OilVaporizationProperties(VapType::DRDT,
                                          1.0, 2.0, {5.0, 6.0},
                                          {false, true}, {7.0, 8.0});
    val2 = PackUnpack(val1);
    DO_CHECKS(OilVaporizationProperties)
#endif
}


BOOST_AUTO_TEST_CASE(Events)
{
#ifdef HAVE_MPI
    Opm::Events val1(Opm::DynamicVector<uint64_t>({1,2,3,4,5}));
    auto val2 = PackUnpack(val1);
    DO_CHECKS(Events)
#endif
}


BOOST_AUTO_TEST_CASE(MLimits)
{
#ifdef HAVE_MPI
    Opm::MLimits val1{1,2,3,4,5,6,7,8,9,10,11,12};
    auto val2 = PackUnpack(val1);
    DO_CHECKS(MLimits)
#endif
}


BOOST_AUTO_TEST_CASE(MessageLimits)
{
#ifdef HAVE_MPI
    std::vector<Opm::MLimits> limits{Opm::MLimits{1,2,3,4,5,6,7,8,9,10,11,12}};
    Opm::MessageLimits val1(Opm::DynamicState<Opm::MLimits>(limits,2));
    auto val2 = PackUnpack(val1);
    DO_CHECKS(MessageLimits)
#endif
}


BOOST_AUTO_TEST_CASE(VFPInjTable)
{
#ifdef HAVE_MPI
    Opm::VFPInjTable val1 = getVFPInjTable();
    auto val2 = PackUnpack(val1);
    DO_CHECKS(VFPInjTable)
#endif
}


BOOST_AUTO_TEST_CASE(VFPProdTable)
{
#ifdef HAVE_MPI
    Opm::VFPProdTable val1 = getVFPProdTable();
    auto val2 = PackUnpack(val1);
    DO_CHECKS(VFPProdTable)
#endif
}


BOOST_AUTO_TEST_CASE(WTESTWell)
{
#ifdef HAVE_MPI
    Opm::WellTestConfig::WTESTWell val1{"test", Opm::WellTestConfig::ECONOMIC,
                                         1.0, 2, 3.0, 4};
    auto val2 = PackUnpack(val1);
    DO_CHECKS(WellTestConfig::WTESTWell)
#endif
}


BOOST_AUTO_TEST_CASE(WellTestConfig)
{
#ifdef HAVE_MPI
    Opm::WellTestConfig::WTESTWell tw{"test", Opm::WellTestConfig::ECONOMIC,
                                         1.0, 2, 3.0, 4};
    Opm::WellTestConfig val1({tw, tw, tw});
    auto val2 = PackUnpack(val1);
    DO_CHECKS(WellTestConfig)
#endif
}


BOOST_AUTO_TEST_CASE(WellPolymerProperties)
{
#ifdef HAVE_MPI
    Opm::WellPolymerProperties val1{1.0, 2.0, 3, 4, 5};
    auto val2 = PackUnpack(val1);
    DO_CHECKS(WellPolymerProperties)
#endif
}


BOOST_AUTO_TEST_CASE(WellFoamProperties)
{
#ifdef HAVE_MPI
    Opm::WellFoamProperties val1{1.0};
    auto val2 = PackUnpack(val1);
    DO_CHECKS(WellFoamProperties)
#endif
}


BOOST_AUTO_TEST_CASE(WellTracerProperties)
{
#ifdef HAVE_MPI
    Opm::WellTracerProperties val1({{"test", 1.0}, {"test2", 2.0}});
    auto val2 = PackUnpack(val1);
    DO_CHECKS(WellTracerProperties)
#endif
}


BOOST_AUTO_TEST_CASE(UDAValue)
{
#ifdef HAVE_MPI
    Opm::UDAValue val1("test");
    auto val2 = PackUnpack(val1);
    DO_CHECKS(UDAValue)
    val1 = Opm::UDAValue(1.0);
    val2 = PackUnpack(val1);
    DO_CHECKS(UDAValue)
#endif
}


BOOST_AUTO_TEST_CASE(Connection)
{
#ifdef HAVE_MPI
    Opm::Connection val1(Opm::Connection::Direction::Y,
                         1.0, Opm::Connection::State::SHUT,
                         2, 3, 4.0, 5.0, 6.0, 7.0, 8.0,
                         {9, 10, 11}, Opm::Connection::CTFKind::Defaulted,
                         12, 13.0, 14.0, true,
                         15, 16, 17.0);
    auto val2 = PackUnpack(val1);
    DO_CHECKS(Connection)
#endif
}


BOOST_AUTO_TEST_CASE(WellInjectionProperties)
{
#ifdef HAVE_MPI
    Opm::Well::WellInjectionProperties val1("test",
                                            Opm::UDAValue(1.0),
                                            Opm::UDAValue("test"),
                                            Opm::UDAValue(2.0),
                                            Opm::UDAValue(3.0),
                                            2.0, 3.0,
                                            4.0, 5.0, 6.0,
                                            7,
                                            true,
                                            8,
                                            Opm::InjectorType::OIL,
                                            Opm::Well::InjectorCMode::BHP);
    auto val2 = PackUnpack(val1);
    DO_CHECKS(Well::WellInjectionProperties)
#endif
}


BOOST_AUTO_TEST_CASE(WellEconProductionLimits)
{
#ifdef HAVE_MPI
    Opm::WellEconProductionLimits val1(1.0, 2.0, 3.0, 4.0, 5.0,
                                       Opm::WellEconProductionLimits::EconWorkover::CONP,
                                       true, "test",
                                       Opm::WellEconProductionLimits::QuantityLimit::POTN,
                                       6.0,
                                       Opm::WellEconProductionLimits::EconWorkover::WELL,
                                       7.0, 8.0, 9.0, 10.0);
    auto val2 = PackUnpack(val1);
    DO_CHECKS(WellEconProductionLimits)
#endif
}


BOOST_AUTO_TEST_CASE(WellGuideRate)
{
#ifdef HAVE_MPI
    Opm::Well::WellGuideRate val1{true, 1.0, Opm::Well::GuideRateTarget::COMB, 2.0};
    auto val2 = PackUnpack(val1);
    DO_CHECKS(Well::WellGuideRate)
#endif
}


BOOST_AUTO_TEST_CASE(WellConnections)
{
#ifdef HAVE_MPI
    Opm::Connection conn(Opm::Connection::Direction::Y,
                         1.0, Opm::Connection::State::SHUT,
                         2, 3, 4.0, 5.0, 6.0, 7.0, 8.0,
                         {9, 10, 11}, Opm::Connection::CTFKind::Defaulted,
                         12, 13.0, 14.0, true,
                         15, 16, 17.0);
    Opm::WellConnections val1(1, 2, 3, {conn, conn});
    auto val2 = PackUnpack(val1);
    DO_CHECKS(WellConnections)
#endif
}


BOOST_AUTO_TEST_CASE(WellProductionProperties)
{
#ifdef HAVE_MPI
    Opm::Well::WellProductionProperties val1("test",
                                             Opm::UDAValue(1.0),
                                             Opm::UDAValue("test"),
                                             Opm::UDAValue(2.0),
                                             Opm::UDAValue(3.0),
                                             Opm::UDAValue(4.0),
                                             Opm::UDAValue(5.0),
                                             Opm::UDAValue(6.0),
                                             5.0, 6.0,
                                             7.0, 8.0,
                                             9,
                                             10.0,
                                             true,
                                             Opm::Well::ProducerCMode::CRAT,
                                             Opm::Well::ProducerCMode::BHP, 11);
    auto val2 = PackUnpack(val1);
    DO_CHECKS(Well::WellProductionProperties)
#endif
}


BOOST_AUTO_TEST_CASE(SpiralICD)
{
#ifdef HAVE_MPI
    Opm::SpiralICD val1(1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8, 9.0,
                        Opm::ICDStatus::OPEN, 10.0);
    auto val2 = PackUnpack(val1);
    DO_CHECKS(SpiralICD)
#endif
}


BOOST_AUTO_TEST_CASE(Valve)
{
#ifdef HAVE_MPI
    Opm::Valve val1(1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, Opm::ICDStatus::OPEN);
    auto val2 = PackUnpack(val1);
    DO_CHECKS(Valve)
#endif
}


BOOST_AUTO_TEST_CASE(Segment)
{
#ifdef HAVE_MPI
    Opm::Segment val1(1, 2, 3, {1, 2}, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, false,
                      Opm::Segment::SegmentType::SICD,
                      std::make_shared<Opm::SpiralICD>(),
                      std::make_shared<Opm::Valve>());
    auto val2 = PackUnpack(val1);
    DO_CHECKS(Segment)
#endif
}


BOOST_AUTO_TEST_CASE(Dimension)
{
#ifdef HAVE_MPI
    Opm::Dimension val1(1.0, 2.0);
    auto val2 = PackUnpack(val1);
    DO_CHECKS(Dimension)
#endif
}


BOOST_AUTO_TEST_CASE(UnitSystem)
{
#ifdef HAVE_MPI
    Opm::UnitSystem val1(Opm::UnitSystem::UnitType::UNIT_TYPE_METRIC);
    auto val2 = PackUnpack(val1);
    DO_CHECKS(UnitSystem)
#endif
}


BOOST_AUTO_TEST_CASE(WellSegments)
{
#ifdef HAVE_MPI
    Opm::Segment seg(1, 2, 3, {1, 2}, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, false,
                     Opm::Segment::SegmentType::SICD,
                     std::make_shared<Opm::SpiralICD>(),
                     std::make_shared<Opm::Valve>());
    Opm::WellSegments val1(Opm::WellSegments::CompPressureDrop::HF_,
                           {seg, seg});

    auto val2 = PackUnpack(val1);
    DO_CHECKS(WellSegments)
#endif
}


BOOST_AUTO_TEST_CASE(Well)
{
#ifdef HAVE_MPI
    Opm::Well val1 = getFullWell();
    auto val2 = PackUnpack(val1);
    DO_CHECKS(Well)
#endif
}


BOOST_AUTO_TEST_CASE(GroupInjectionProperties)
{
#ifdef HAVE_MPI
    Opm::Group::GroupInjectionProperties val1{Opm::Phase::WATER,
                                              Opm::Group::InjectionCMode::REIN,
                                              Opm::UDAValue(1.0),
                                              Opm::UDAValue(2.0),
                                              Opm::UDAValue(3.0),
                                              Opm::UDAValue(4.0),
                                              "test1", "test2", 5};

    auto val2 = PackUnpack(val1);
    DO_CHECKS(Group::GroupInjectionProperties)
#endif
}


BOOST_AUTO_TEST_CASE(GroupProductionProperties)
{
#ifdef HAVE_MPI
    Opm::Group::GroupProductionProperties val1{Opm::Group::ProductionCMode::PRBL,
                                               Opm::Group::ExceedAction::WELL,
                                               Opm::UDAValue(1.0),
                                               Opm::UDAValue(2.0),
                                               Opm::UDAValue(3.0),
                                               Opm::UDAValue(4.0),
                                               5.0, Opm::Group::GuideRateTarget::COMB,
                                               6.0, 7};

    auto val2 = PackUnpack(val1);
    DO_CHECKS(Group::GroupProductionProperties)
#endif
}


BOOST_AUTO_TEST_CASE(Group)
{
#ifdef HAVE_MPI
    Opm::UnitSystem unitSystem;

    std::map<Opm::Phase, Opm::Group::GroupInjectionProperties> injection;
    Opm::Group val1("test1", 1, 2, 3.0, unitSystem,
                    Opm::Group::GroupType::PRODUCTION,
                    4.0, true, false, 5, "test2",
                    Opm::IOrderSet<std::string>({"test3", "test4"}, {"test3","test4"}),
                    Opm::IOrderSet<std::string>({"test5", "test6"}, {"test5","test6"}),
                    injection,
                    Opm::Group::GroupProductionProperties());

    auto val2 = PackUnpack(val1);
    DO_CHECKS(Group)
#endif
}


BOOST_AUTO_TEST_CASE(WList)
{
#ifdef HAVE_MPI
    Opm::WList val1({"test1", "test2", "test3"});
    auto val2 = PackUnpack(val1);
    DO_CHECKS(WList)
#endif
}


BOOST_AUTO_TEST_CASE(WListManager)
{
#ifdef HAVE_MPI
    Opm::WList wl({"test1", "test2", "test3"});
    std::map<std::string,Opm::WList> data{{"test", wl}, {"test2", wl}};
    Opm::WListManager val1(data);
    auto val2 = PackUnpack(val1);
    DO_CHECKS(WListManager)
#endif
}


BOOST_AUTO_TEST_CASE(UDQASTNode)
{
#ifdef HAVE_MPI
  std::shared_ptr<Opm::UDQASTNode> n0;
  std::shared_ptr<Opm::UDQASTNode> n1(new Opm::UDQASTNode(Opm::UDQVarType::NONE,
                                                          Opm::UDQTokenType::error,
                                                          "test1", 1.0, {"test2"},
                                                          n0, n0));
    Opm::UDQASTNode val1(Opm::UDQVarType::NONE,
                         Opm::UDQTokenType::error,
                         "test", 1.0, {"test3"}, n1, n1);
    auto val2 = PackUnpack(val1);
    DO_CHECKS(UDQASTNode)
#endif
}

BOOST_AUTO_TEST_CASE(UDQDefine)
{
#ifdef HAVE_MPI
    std::shared_ptr<Opm::UDQASTNode> n0;
    Opm::UDQASTNode n1(Opm::UDQVarType::NONE,
                       Opm::UDQTokenType::error,
                       "test", 1.0, {"test1", "test2"}, n0, n0);
    Opm::UDQDefine val1("test", std::make_shared<Opm::UDQASTNode>(n1),
                        Opm::UDQVarType::NONE, "test2");
    auto val2 = PackUnpack(val1);
    DO_CHECKS(UDQDefine)
#endif
}


BOOST_AUTO_TEST_CASE(UDQAssign)
{
#ifdef HAVE_MPI
    Opm::UDQAssign val1("test", Opm::UDQVarType::NONE,
                        {Opm::UDQAssign::AssignRecord{{"test1"}, 1.0},
                         Opm::UDQAssign::AssignRecord{{"test2"}, 2.0}});
    auto val2 = PackUnpack(val1);
    DO_CHECKS(UDQAssign)
#endif
}


BOOST_AUTO_TEST_CASE(UDQIndex)
{
#ifdef HAVE_MPI
    Opm::UDQIndex val1(1, 2, Opm::UDQAction::ASSIGN, Opm::UDQVarType::WELL_VAR);
    auto val2 = PackUnpack(val1);
    DO_CHECKS(UDQIndex)
#endif
}


BOOST_AUTO_TEST_CASE(UDQConfig)
{
#ifdef HAVE_MPI
    Opm::UDQConfig val1 = getUDQConfig();
    auto val2 = PackUnpack(val1);
    DO_CHECKS(UDQConfig)
#endif
}


BOOST_AUTO_TEST_CASE(UDQActiveInputRecord)
{
#ifdef HAVE_MPI
    Opm::UDQActive::InputRecord val1(1, "test1", "test2",
                                     Opm::UDAControl::WCONPROD_ORAT);
    auto val2 = PackUnpack(val1);
    DO_CHECKS(UDQActive::InputRecord)
#endif
}


BOOST_AUTO_TEST_CASE(UDQActiveRecord)
{
#ifdef HAVE_MPI
    Opm::UDQActive::Record val1("test1", 1, 2, "test2",
                                Opm::UDAControl::WCONPROD_ORAT);
    auto val2 = PackUnpack(val1);
    DO_CHECKS(UDQActive::Record)
#endif
}


BOOST_AUTO_TEST_CASE(UDQActive)
{
#ifdef HAVE_MPI
    Opm::UDQActive val1({Opm::UDQActive::InputRecord(1, "test1", "test2",
                                                     Opm::UDAControl::WCONPROD_ORAT)},
                        {Opm::UDQActive::Record("test1", 1, 2, "test2",
                                                  Opm::UDAControl::WCONPROD_ORAT)},
                        {{"test1", 1}}, {{"test2", 2}});
    auto val2 = PackUnpack(val1);
    DO_CHECKS(UDQActive)
#endif
}


BOOST_AUTO_TEST_CASE(AquiferCT)
{
#ifdef HAVE_MPI
    Opm::AquiferCT val1 = getAquiferCT();
    auto val2 = PackUnpack(val1);
    DO_CHECKS(AquiferCT);
#endif
}


BOOST_AUTO_TEST_CASE(Aquifetp)
{
#ifdef HAVE_MPI
    Opm::Aquifetp val1 = getAquifetp();
    auto val2 = PackUnpack(val1);
    DO_CHECKS(Aquifetp);
#endif
}


BOOST_AUTO_TEST_CASE(Aquancon)
{
#ifdef HAVE_MPI
    Opm::Aquancon val1 = getAquancon();
    auto val2 = PackUnpack(val1);
    DO_CHECKS(Aquancon);
#endif
}

BOOST_AUTO_TEST_CASE(AquferConfig)
{
#ifdef HAVE_MPI
    Opm::Aquifetp fetp = getAquifetp();
    Opm::AquiferCT ct = getAquiferCT();
    Opm::Aquancon conn = getAquancon();
    Opm::AquiferConfig val1(fetp, ct, conn);
    auto val2 = PackUnpack(val1);
    DO_CHECKS(AquiferConfig);
#endif
}




BOOST_AUTO_TEST_CASE(GuideRateModel)
{
#ifdef HAVE_MPI
    Opm::GuideRateModel val1 = getGuideRateModel();
    auto val2 = PackUnpack(val1);
    DO_CHECKS(GuideRateModel)
#endif
}


BOOST_AUTO_TEST_CASE(GuideRateConfigGroup)
{
#ifdef HAVE_MPI
    Opm::GuideRateConfig::GroupTarget val1 = getGuideRateConfigGroup();
    auto val2 = PackUnpack(val1);
    DO_CHECKS(GuideRateConfig::GroupTarget)
#endif
}


BOOST_AUTO_TEST_CASE(GuideRateConfigWell)
{
#ifdef HAVE_MPI
    Opm::GuideRateConfig::WellTarget val1 = getGuideRateConfigWell();
    auto val2 = PackUnpack(val1);
    DO_CHECKS(GuideRateConfig::WellTarget)
#endif
}


BOOST_AUTO_TEST_CASE(GuideRateConfig)
{
#ifdef HAVE_MPI
    auto model = std::make_shared<Opm::GuideRateModel>(getGuideRateModel());
    Opm::GuideRateConfig val1(model,
                              {{"test1", getGuideRateConfigWell()}},
                              {{"test2", getGuideRateConfigGroup()}});
    auto val2 = PackUnpack(val1);
    DO_CHECKS(GuideRateConfig)
#endif
}


BOOST_AUTO_TEST_CASE(GConSaleGroup)
{
#ifdef HAVE_MPI
    Opm::GConSale::GCONSALEGroup val1{Opm::UDAValue(1.0),
                                      Opm::UDAValue(2.0),
                                      Opm::UDAValue(3.0),
                                      Opm::GConSale::MaxProcedure::PLUG,
                                      4.0, Opm::UnitSystem()};
    auto val2 = PackUnpack(val1);
    DO_CHECKS(GConSale::GCONSALEGroup)
#endif
}


BOOST_AUTO_TEST_CASE(GConSale)
{
#ifdef HAVE_MPI
    Opm::GConSale::GCONSALEGroup group{Opm::UDAValue(1.0),
                                       Opm::UDAValue(2.0),
                                       Opm::UDAValue(3.0),
                                       Opm::GConSale::MaxProcedure::PLUG,
                                       4.0, Opm::UnitSystem()};
    Opm::GConSale val1({{"test1", group}, {"test2", group}});
    auto val2 = PackUnpack(val1);
    DO_CHECKS(GConSale)
#endif
}


BOOST_AUTO_TEST_CASE(GConSumpGroup)
{
#ifdef HAVE_MPI
    Opm::GConSump::GCONSUMPGroup val1{Opm::UDAValue(1.0),
                                      Opm::UDAValue(2.0),
                                      "test",
                                      3.0, Opm::UnitSystem()};
    auto val2 = PackUnpack(val1);
    DO_CHECKS(GConSump::GCONSUMPGroup)
#endif
}


BOOST_AUTO_TEST_CASE(GConSump)
{
#ifdef HAVE_MPI
    Opm::GConSump::GCONSUMPGroup group{Opm::UDAValue(1.0),
                                       Opm::UDAValue(2.0),
                                       "test",
                                       3.0, Opm::UnitSystem()};
    Opm::GConSump val1({{"test1", group}, {"test2", group}});
    auto val2 = PackUnpack(val1);
    DO_CHECKS(GConSump)
#endif
}


BOOST_AUTO_TEST_CASE(RFTConfig)
{
#ifdef HAVE_MPI
    Opm::RFTConfig val1(getTimeMap(),
                        std::size_t{1729},
                        {true, 1},
                        {{"test1", 2}, {"test2", 3}},
                        {{"test3", 2}},
                        {{"test1", {{{Opm::RFTConfig::RFT::TIMESTEP, 3}}, 4}}},
                        {{"test2", {{{Opm::RFTConfig::PLT::REPT, 5}}, 6}}});
    auto val2 = PackUnpack(val1);
    DO_CHECKS(RFTConfig)
#endif
}


BOOST_AUTO_TEST_CASE(DeckItem)
{
#ifdef HAVE_MPI
    Opm::DeckItem val1({1.0}, {2}, {"test3"}, {Opm::UDAValue(4)},
                       Opm::type_tag::string, "test5",
                       {Opm::value::status::deck_value},
                       true,
                       {Opm::Dimension(7.0, 8.0)},
                       {Opm::Dimension(10.0, 11.0)});

    auto val2 = PackUnpack(val1);
    DO_CHECKS(DeckItem)
#endif
}


BOOST_AUTO_TEST_CASE(DeckRecord)
{
#ifdef HAVE_MPI
    Opm::DeckRecord val1 = getDeckRecord();
    auto val2 = PackUnpack(val1);
    DO_CHECKS(DeckRecord)
#endif
}


BOOST_AUTO_TEST_CASE(Location)
{
#ifdef HAVE_MPI
    Opm::Location val1{"test", 1};
    auto val2 = PackUnpack(val1);
    DO_CHECKS(Location)
#endif
}


BOOST_AUTO_TEST_CASE(DeckKeyword)
{
#ifdef HAVE_MPI
    Opm::DeckKeyword val1("test", {"test",1},
                          {getDeckRecord(), getDeckRecord()}, true, false);
    auto val2 = PackUnpack(val1);
    DO_CHECKS(DeckKeyword)
#endif
}


BOOST_AUTO_TEST_CASE(Deck)
{
#ifdef HAVE_MPI
    std::unique_ptr<Opm::UnitSystem> unitSys(new Opm::UnitSystem);
    Opm::Deck val1({Opm::DeckKeyword("test", {"test",1},
                                     {getDeckRecord(), getDeckRecord()}, true, false)},
                   Opm::UnitSystem(), unitSys.get(),
                   "test2", "test3", 2);
    auto val2 = PackUnpack(val1);
    DO_CHECKS(Deck)
#endif
}


BOOST_AUTO_TEST_CASE(Tuning)
{
#ifdef HAVE_MPI
    Opm::Tuning val1 = getTuning();
    auto val2 = PackUnpack(val1);
    DO_CHECKS(Tuning)
#endif
}


BOOST_AUTO_TEST_CASE(ASTNode)
{
#ifdef HAVE_MPI
    Opm::Action::ASTNode child(number, FuncType::field, "test3", {"test2"}, 2.0, {});
    Opm::Action::ASTNode val1(number, FuncType::field, "test1", {"test2"}, 1.0, {child});
    auto val2 = PackUnpack(val1);
    DO_CHECKS(Action::ASTNode)
#endif
}


BOOST_AUTO_TEST_CASE(AST)
{
#ifdef HAVE_MPI
    std::shared_ptr<Opm::Action::ASTNode> node;
    node.reset(new Opm::Action::ASTNode(number, FuncType::field,
                                        "test1", {"test2"}, 1.0, {}));
    Opm::Action::AST val1(node);
    auto val2 = PackUnpack(val1);
    DO_CHECKS(Action::AST)
#endif
}


BOOST_AUTO_TEST_CASE(Quantity)
{
#ifdef HAVE_MPI
    Opm::Action::Quantity val1;
    val1.quantity = "test1";
    val1.args = {"test2", "test3"};
    auto val2 = PackUnpack(val1);
    DO_CHECKS(Action::Quantity)
#endif
}


BOOST_AUTO_TEST_CASE(Condition)
{
#ifdef HAVE_MPI
    Opm::Action::Condition val1 = getCondition();
    auto val2 = PackUnpack(val1);
    DO_CHECKS(Action::Condition)
#endif
}


BOOST_AUTO_TEST_CASE(ActionX)
{
#ifdef HAVE_MPI
    Opm::Action::ActionX val1 = getActionX();
    auto val2 = PackUnpack(val1);
    DO_CHECKS(Action::ActionX)
#endif
}


BOOST_AUTO_TEST_CASE(Actions)
{
#ifdef HAVE_MPI
    Opm::Action::Actions val1({getActionX()});
    auto val2 = PackUnpack(val1);
    DO_CHECKS(Action::Actions)
#endif
}


BOOST_AUTO_TEST_CASE(Schedule)
{
#ifdef HAVE_MPI
    Opm::UnitSystem unitSystem;
    Opm::Schedule::WellMap wells;
    wells.insert({"test", {{std::make_shared<Opm::Well>(getFullWell())},1}});
    Opm::Schedule::GroupMap groups;
    std::map<Opm::Phase, Opm::Group::GroupInjectionProperties> injection;
    groups.insert({"test", {{std::make_shared<Opm::Group>("test1", 1, 2, 3.0, unitSystem,
                                                          Opm::Group::GroupType::PRODUCTION,
                                                          4.0, true, false, 5, "test2",
                                                          Opm::IOrderSet<std::string>({"test3", "test4"}, {"test3","test4"}),
                                                          Opm::IOrderSet<std::string>({"test5", "test6"}, {"test5","test6"}),
                                                          injection,
                                                          Opm::Group::GroupProductionProperties())},1}});
    using VapType = Opm::OilVaporizationProperties::OilVaporization;
    Opm::DynamicState<Opm::OilVaporizationProperties> oilvap{{Opm::OilVaporizationProperties(VapType::VAPPARS,
                                                                                   1.0, 2.0, {5.0, 6.0},
                                                                                   {false, true}, {7.0, 8.0})},1};
    Opm::Events events(Opm::DynamicVector<uint64_t>({1,2,3,4,5}));
    std::unique_ptr<Opm::UnitSystem> unitSys(new Opm::UnitSystem);
    Opm::Deck sdeck({Opm::DeckKeyword("test", {"test",1},
                                     {getDeckRecord(), getDeckRecord()}, true, false)},
                   Opm::UnitSystem(), unitSys.get(),
                   "test2", "test3", 2);
    Opm::DynamicVector<Opm::Deck> modifierDeck({sdeck});
    std::vector<Opm::MLimits> limits{Opm::MLimits{1,2,3,4,5,6,7,8,9,10,11,12}};
    Opm::MessageLimits messageLimits(Opm::DynamicState<Opm::MLimits>(limits,2));
    Opm::Runspec runspec(Opm::Phases(true, true, true, false, true, false, true, false),
                         Opm::Tabdims(1,2,3,4,5,6),
                         Opm::EndpointScaling(std::bitset<4>(13)),
                         Opm::Welldims(1,2,3,4),
                         Opm::WellSegmentDims(1,2,3),
                         Opm::UDQParams(true, 1, 2.0, 3.0, 4.0),
                         Opm::EclHysterConfig(true, 1, 2),
                         Opm::Actdims(1,2,3,4),
                         Opm::SatFuncControls(5.6e-7, Opm::SatFuncControls::ThreePhaseOilKrModel::Stone1));
    Opm::Schedule::VFPProdMap vfpProd {{1, {{std::make_shared<Opm::VFPProdTable>(getVFPProdTable())},1}}};
    Opm::Schedule::VFPInjMap vfpIn{{1, {{std::make_shared<Opm::VFPInjTable>(getVFPInjTable())},1}}};
    Opm::WellTestConfig::WTESTWell tw{"test", Opm::WellTestConfig::ECONOMIC,
                                         1.0, 2, 3.0, 4};
    Opm::WellTestConfig wtc({tw, tw, tw});

    Opm::WList wl({"test1", "test2", "test3"});
    std::map<std::string,Opm::WList> data{{"test", wl}, {"test2", wl}};
    Opm::WListManager wlm(data);

    Opm::UDQActive udqa({Opm::UDQActive::InputRecord(1, "test1", "test2",
                                                     Opm::UDAControl::WCONPROD_ORAT)},
                        {Opm::UDQActive::Record("test1", 1, 2, "test2",
                                                  Opm::UDAControl::WCONPROD_ORAT)},
                        {{"test1", 1}}, {{"test2", 2}});

    auto model = std::make_shared<Opm::GuideRateModel>(getGuideRateModel());
    Opm::GuideRateConfig grc(model,
                             {{"test1", getGuideRateConfigWell()}},
                             {{"test2", getGuideRateConfigGroup()}});

    Opm::GConSale::GCONSALEGroup group{Opm::UDAValue(1.0),
                                       Opm::UDAValue(2.0),
                                       Opm::UDAValue(3.0),
                                       Opm::GConSale::MaxProcedure::PLUG,
                                       4.0, Opm::UnitSystem()};
    Opm::GConSale gcs({{"test1", group}, {"test2", group}});

    Opm::GConSump::GCONSUMPGroup grp{Opm::UDAValue(1.0),
                                     Opm::UDAValue(2.0),
                                     "test",
                                     3.0, Opm::UnitSystem()};
    Opm::GConSump gcm({{"test1", grp}, {"test2", grp}});

    Opm::Action::Actions acnts({getActionX()});

    Opm::RFTConfig rftc(getTimeMap(),
                        std::size_t{1729},
                        {true, 1},
                        {{"test1", 2}, {"test2", 3}},
                        {{"test3", 2}},
                        {{"test1", {{{Opm::RFTConfig::RFT::TIMESTEP, 3}}, 4}}},
                        {{"test2", {{{Opm::RFTConfig::PLT::REPT, 5}}, 6}}});

    Opm::Schedule val1(getTimeMap(),
                       wells,
                       groups,
                       oilvap,
                       events,
                       modifierDeck,
                       Opm::DynamicState<Opm::Tuning>({getTuning()}, 1),
                       messageLimits,
                       runspec,
                       vfpProd,
                       vfpIn,
                       {{std::make_shared<Opm::WellTestConfig>(wtc)}, 1},
                       {{std::make_shared<Opm::WListManager>(wlm)}, 1},
                       {{std::make_shared<Opm::UDQConfig>(getUDQConfig())}, 1},
                       {{std::make_shared<Opm::UDQActive>(udqa)}, 1},
                       {{std::make_shared<Opm::GuideRateConfig>(grc)}, 1},
                       {{std::make_shared<Opm::GConSale>(gcs)}, 1},
                       {{std::make_shared<Opm::GConSump>(gcm)}, 1},
                       {{Opm::Well::ProducerCMode::CRAT}, 1},
                       {{std::make_shared<Opm::Action::Actions>(acnts)}, 1},
                       rftc,
                       {std::vector<int>{1}, 1},
                       getRestartConfig(),
                       {{"test", events}});

    auto val2 = PackUnpack2(val1);
    DO_CHECKS(Schedule)
#endif
}

BOOST_AUTO_TEST_CASE(BrineDensityTable)
{
#ifdef HAVE_MPI
    Opm::BrineDensityTable val1({1.0, 2.0, 3.0});
    auto val2 = PackUnpack(val1);
    DO_CHECKS(BrineDensityTable)
#endif
}

BOOST_AUTO_TEST_CASE(SummaryConfigNode)
{
#ifdef HAVE_MPI
    auto val1 = Opm::SummaryConfigNode{"test1", Opm::SummaryConfigNode::Category::Region,
                                 Opm::Location{"test2", 1}}
                                 .parameterType(Opm::SummaryConfigNode::Type::Pressure)
                                 .namedEntity("test3")
                                 .number(2)
                                 .isUserDefined(true);

    auto val2 = PackUnpack2(val1);
    DO_CHECKS(SummaryConfigNode)
#endif
}


BOOST_AUTO_TEST_CASE(SummaryConfig)
{
#ifdef HAVE_MPI
    auto node = Opm::SummaryConfigNode{"test1", Opm::SummaryConfigNode::Category::Region,
                                 Opm::Location{"test2", 1}}
                                 .parameterType(Opm::SummaryConfigNode::Type::Pressure)
                                 .namedEntity("test3")
                                 .number(2)
                                 .isUserDefined(true);
    Opm::SummaryConfig val1({node}, {"test1", "test2"}, {"test3", "test4"});

    auto val2 = PackUnpack2(val1);
    DO_CHECKS(SummaryConfig)
#endif
}


BOOST_AUTO_TEST_CASE(PvtwsaltTable)
{
#ifdef HAVE_MPI
    Opm::PvtwsaltTable val1(1.0, 2.0, {3.0, 4.0, 5.0});
    auto val2 = PackUnpack(val1);
    DO_CHECKS(PvtwsaltTable)
#endif
}


BOOST_AUTO_TEST_CASE(WellBrineProperties)
{
#ifdef HAVE_MPI
    Opm::WellBrineProperties val1{1.0};
    auto val2 = PackUnpack(val1);
    DO_CHECKS(WellBrineProperties)
#endif
}


BOOST_AUTO_TEST_CASE(MULTREGTRecord)
{
#ifdef HAVE_MPI
    Opm::MULTREGTRecord val1{1, 2, 3.0, 4, Opm::MULTREGT::ALL, "test"};
    auto val2 = PackUnpack2(val1);
    DO_CHECKS(MULTREGTRecord)
#endif
}


BOOST_AUTO_TEST_CASE(MULTREGTScanner)
{
#ifdef HAVE_MPI
    std::vector<Opm::MULTREGTRecord> records{{1, 2, 3.0, 4, Opm::MULTREGT::ALL, "test1"}};
    std::map<std::pair<int, int>, int> searchRecord{{{5,6},0}};
    Opm::MULTREGTScanner::ExternalSearchMap searchMap;
    searchMap.insert({"test2", searchRecord});
    Opm::MULTREGTScanner val1({1, 2, 3},
                              records,
                              searchMap,
                              {{"test3", {7,8}}},
                              "test4");

    auto val2 = PackUnpack2(val1);
    DO_CHECKS(MULTREGTScanner)
#endif
}


BOOST_AUTO_TEST_CASE(EclipseConfig)
{
#ifdef HAVE_MPI
    Opm::IOConfig io(true, false, true, false, false, true, "test1", true,
                     "test2", true, "test3", false);
    Opm::InitConfig init(Opm::Equil({getEquilRecord(), getEquilRecord()}),
                         getFoamConfig(),
                         true, true, true, 20, "test1");
    Opm::EclipseConfig val1{init, io};

    auto val2 = PackUnpack(val1);
    DO_CHECKS(EclipseConfig)
#endif
}


BOOST_AUTO_TEST_CASE(TransMult)
{
#ifdef HAVE_MPI
    std::vector<Opm::MULTREGTRecord> records{{1, 2, 3.0, 4, Opm::MULTREGT::ALL, "test1"}};
    std::map<std::pair<int, int>, int> searchRecord{{{5,6},0}};
    Opm::MULTREGTScanner::ExternalSearchMap searchMap;
    searchMap.insert({"test2", searchRecord});
    Opm::MULTREGTScanner scanner({1, 2, 3},
                                 records,
                                 searchMap,
                                 {{"test3", {7,8}}},
                                 "test4");

    Opm::TransMult val1({1, 2, 3},
                        {{Opm::FaceDir::YPlus, {4.0, 5.0}}},
                        {{Opm::FaceDir::ZPlus, "test1"}},
                        scanner);
    auto val2 = PackUnpack2(val1);
    DO_CHECKS(TransMult)
#endif
}


BOOST_AUTO_TEST_CASE(FaultFace)
{
#ifdef HAVE_MPI
    Opm::FaultFace val1({1,2,3,4,5,6}, Opm::FaceDir::YPlus);
    auto val2 = PackUnpack2(val1);
    DO_CHECKS(FaultFace)
#endif
}


BOOST_AUTO_TEST_CASE(Fault)
{
#ifdef HAVE_MPI
    Opm::Fault val1("test", 1.0, {{{1,2,3,4,5,6}, Opm::FaceDir::YPlus}});
    auto val2 = PackUnpack2(val1);
    DO_CHECKS(Fault)
#endif
}

BOOST_AUTO_TEST_CASE(WellType)
{
#ifdef HAVE_MPI
    Opm::WellType val1(true, Opm::Phase::OIL);
    auto val2 = PackUnpack(val1);
    DO_CHECKS(WellType)
#endif
}

BOOST_AUTO_TEST_CASE(DenT)
{
#ifdef HAVE_MPI
    Opm::DenT val1 = getDenT();
    auto val2 = PackUnpack(val1);
    DO_CHECKS(DenT)
#endif
}

BOOST_AUTO_TEST_CASE(FaultCollection)
{
#ifdef HAVE_MPI
    Opm::Fault fault("test", 1.0, {{{1,2,3,4,5,6}, Opm::FaceDir::YPlus}});
    Opm::OrderedMap<std::string, Opm::Fault> faults;
    faults.insert({"test2", fault});
    Opm::FaultCollection val1(faults);
    auto val2 = PackUnpack2(val1);
    DO_CHECKS(FaultCollection)
#endif
}


BOOST_AUTO_TEST_CASE(SolventDensityTable)
{
#ifdef HAVE_MPI
    Opm::SolventDensityTable val1({1.0, 2.0, 3.0});
    auto val2 = PackUnpack(val1);
    DO_CHECKS(SolventDensityTable)
#endif
}


BOOST_AUTO_TEST_CASE(GridDims)
{
#ifdef HAVE_MPI
    Opm::GridDims val1{ 1,  2,  3};
    auto val2 = PackUnpack2(val1);
    DO_CHECKS(GridDims)
#endif
}


BOOST_AUTO_TEST_CASE(PlyshlogTable)
{
#ifdef HAVE_MPI
    Opm::OrderedMap<std::string, Opm::TableColumn> data;
    data.insert({"test3", getTableColumn()});
    Opm::PlyshlogTable val1(getTableSchema(), data, true, 1.0, 2.0, 3.0, true, true);
    auto val2 = PackUnpack(val1);
    DO_CHECKS(PlyshlogTable)
#endif
}


BOOST_AUTO_TEST_CASE(RocktabTable)
{
#ifdef HAVE_MPI
    Opm::OrderedMap<std::string, Opm::TableColumn> data;
    data.insert({"test3", getTableColumn()});
    Opm::RocktabTable val1(getTableSchema(), data, true, true);
    auto val2 = PackUnpack(val1);
    DO_CHECKS(RocktabTable)
#endif
}


bool init_unit_test_func()
{
    return true;
}


int main(int argc, char** argv)
{
    Dune::MPIHelper::instance(argc, argv);
    return boost::unit_test::unit_test_main(&init_unit_test_func, argc, argv);
}
