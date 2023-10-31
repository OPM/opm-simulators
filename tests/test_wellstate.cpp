/*
  Copyright 2018 Equinor ASA.

  This file is part of the Open Porous Media Project (OPM).

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

#include <algorithm>
#include <config.h>
#include <functional>
#include <vector>

#define BOOST_TEST_MODULE WellStateFIBOTest

#include "MpiFixture.hpp"
#include <opm/common/ErrorMacros.hpp>
#include <opm/simulators/wells/GlobalWellInfo.hpp>
#include <opm/simulators/wells/ParallelWellInfo.hpp>
#include <opm/simulators/wells/WellState.hpp>
#include <opm/simulators/wells/SingleWellState.hpp>
#include <opm/simulators/wells/SegmentState.hpp>
#include <opm/simulators/wells/WellContainer.hpp>
#include <opm/simulators/wells/PerfData.hpp>
#include <opm/input/eclipse/Python/Python.hpp>

#include <boost/test/unit_test.hpp>

#include <opm/input/eclipse/Deck/Deck.hpp>
#include <opm/input/eclipse/EclipseState/EclipseState.hpp>
#include <opm/input/eclipse/Parser/Parser.hpp>
#include <opm/input/eclipse/Parser/ParseContext.hpp>
#include <opm/input/eclipse/Schedule/MSW/WellSegments.hpp>
#include <opm/input/eclipse/Schedule/Schedule.hpp>
#include <opm/input/eclipse/Schedule/SummaryState.hpp>
#include <opm/input/eclipse/Schedule/Well/Well.hpp>
#include <opm/input/eclipse/Schedule/Well/WellConnections.hpp>
#include <opm/input/eclipse/Units/Units.hpp>
#include <opm/common/utility/TimeService.hpp>

#include <opm/grid/GridHelpers.hpp>

#include <opm/core/props/BlackoilPhases.hpp>
#include <opm/core/props/phaseUsageFromDeck.hpp>

#include <opm/grid/GridManager.hpp>

#include <chrono>
#include <cstddef>
#include <string>

BOOST_GLOBAL_FIXTURE(MPIFixture);

struct Setup
{
    Setup(const std::string& filename)
        : Setup(Opm::Parser{}.parseFile(filename))
    {}

    Setup(const Opm::Deck& deck)
        : es   (deck)
        , pu   (Opm::phaseUsageFromDeck(es))
        , grid (es.getInputGrid())
        , python( std::make_shared<Opm::Python>() )
        , sched(deck, es, python)
        , st(Opm::TimeService::from_time_t(sched.getStartTime()))
    {
        initWellPerfData();
    }

    void initWellPerfData()
    {
        const auto& wells = sched.getWells(0);
        const auto& cartDims = Opm::UgGridHelpers::cartDims(*grid.c_grid());
        const int* compressed_to_cartesian = Opm::UgGridHelpers::globalCell(*grid.c_grid());
        std::vector<int> cartesian_to_compressed(cartDims[0] * cartDims[1] * cartDims[2], -1);
        for (int ii = 0; ii < Opm::UgGridHelpers::numCells(*grid.c_grid()); ++ii) {
            cartesian_to_compressed[compressed_to_cartesian[ii]] = ii;
        }
        well_perf_data.resize(wells.size());
        int well_index = 0;
        for (const auto& well : wells) {
            well_perf_data[well_index].clear();
            well_perf_data[well_index].reserve(well.getConnections().size());
            for (const auto& completion : well.getConnections()) {
                if (completion.state() == Opm::Connection::State::OPEN) {
                    const int i = completion.getI();
                    const int j = completion.getJ();
                    const int k = completion.getK();
                    const int cart_grid_indx = i + cartDims[0] * (j + cartDims[1] * k);
                    const int active_index = cartesian_to_compressed[cart_grid_indx];
                    if (active_index < 0) {
                        const std::string msg
                            = ("Cell with i,j,k indices " + std::to_string(i) + " " + std::to_string(j) + " "
                               + std::to_string(k) + " not found in grid (well = " + well.name() + ").");
                        OPM_THROW(std::runtime_error, msg);
                    } else {
                        Opm::PerforationData pd;
                        pd.cell_index = active_index;
                        pd.connection_transmissibility_factor = completion.CF();
                        pd.connection_d_factor = completion.dFactor();
                        pd.satnum_id = completion.satTableId();
                        well_perf_data[well_index].push_back(pd);
                    }
                } else {
                    if (completion.state() != Opm::Connection::State::SHUT) {
                        OPM_THROW(std::runtime_error,
                                  "Completion state: " +
                                  Opm::Connection::State2String(completion.state()) +
                                  " not handled");
                    }
                }
            }
            ++well_index;
        }
    }

    Opm::EclipseState es;
    Opm::PhaseUsage   pu;
    Opm::GridManager  grid;
    std::shared_ptr<Opm::Python> python;
    Opm::Schedule     sched;
    Opm::SummaryState st;
    std::vector<std::vector<Opm::PerforationData>> well_perf_data;
};

namespace {
    Opm::WellState
    buildWellState(const Setup& setup, const std::size_t timeStep,
                   std::vector<Opm::ParallelWellInfo>& pinfos)
    {
        auto state  = Opm::WellState{setup.pu};

        const auto cpress =
            std::vector<double>(setup.grid.c_grid()->number_of_cells,
                                100.0*Opm::unit::barsa);

        auto wells = setup.sched.getWells(timeStep);
        pinfos.resize(wells.size());
        std::vector<std::reference_wrapper<Opm::ParallelWellInfo>> ppinfos;
        auto pw = pinfos.begin();

        for (const auto& well : wells)
        {
            *pw = {well.name()};
            ppinfos.push_back(std::ref(*pw));
            pw->communicateFirstPerforation(true);
            ++pw;
        }

        state.init(cpress, setup.sched,
                   wells, ppinfos,
                   timeStep, nullptr, setup.well_perf_data, setup.st);

        state.initWellStateMSWell(setup.sched.getWells(timeStep),
                                  nullptr);

        return state;
    }


    void setSegPress(const std::vector<Opm::Well>& wells,
                     Opm::WellState& wstate)
    {
        const auto nWell = wells.size();


        for (auto wellID = 0*nWell; wellID < nWell; ++wellID) {
            const auto& well     = wells[wellID];
            if (! well.isMultiSegment()) {
                continue;
            }

            const auto  pressTop = 100.0 * wellID;
            auto& ws = wstate.well(wellID);
            auto& segments = ws.segments;
            segments.pressure[0] = pressTop;

            const auto& segSet = well.getSegments();
            const auto  nSeg   = segSet.size();

            for (auto segID = 0*nSeg + 1; segID < nSeg; ++segID) {
                // One-based numbering scheme for segments.
                const auto segNo = segSet[segID].segmentNumber();
                segments.pressure[segNo - 1] = pressTop + 1.0*(segNo - 1);
            }
        }
    }


  void setSegRates(const std::vector<Opm::Well>& wells,
                     const Opm::PhaseUsage&               pu,
                     Opm::WellState& wstate)
    {
        const auto wat = pu.phase_used[Opm::BlackoilPhases::Aqua];
        const auto iw  = wat ? pu.phase_pos[Opm::BlackoilPhases::Aqua] : -1;

        const auto oil = pu.phase_used[Opm::BlackoilPhases::Liquid];
        const auto io  = oil ? pu.phase_pos[Opm::BlackoilPhases::Liquid] : -1;

        const auto gas = pu.phase_used[Opm::BlackoilPhases::Vapour];
        const auto ig  = gas ? pu.phase_pos[Opm::BlackoilPhases::Vapour] : -1;

        const auto np = wstate.numPhases();

        const auto nWell = wells.size();


        for (auto wellID = 0*nWell; wellID < nWell; ++wellID) {
            const auto& well     = wells[wellID];
            if (! well.isMultiSegment()) {
                continue;
            }

            const auto  rateTop  = 1000.0 * wellID;
            auto& ws = wstate.well(wellID);
            auto& segments = ws.segments;
            auto& segRates = segments.rates;

            if (wat) { segRates[iw] = rateTop; }
            if (oil) { segRates[io] = rateTop; }
            if (gas) { segRates[ig] = rateTop; }

            const auto& segSet = well.getSegments();
            const auto  nSeg   = segSet.size();

            for (auto segID = 0*nSeg + 1; segID < nSeg; ++segID) {
                // One-based numbering scheme for segments.
                const auto segNo = segSet[segID].segmentNumber();

                auto* rates = &segRates[(segNo - 1) * np];

                if (wat) { rates[iw] = rateTop + 100.0*(segNo - 1); }
                if (oil) { rates[io] = rateTop + 200.0*(segNo - 1); }
                if (gas) { rates[ig] = rateTop + 400.0*(segNo - 1); }
            }
        }
    }
} // Anonymous

BOOST_AUTO_TEST_SUITE(Segment)

// ---------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(Linearisation)
{
    const Setup setup{ "msw.data" };
    const auto tstep = std::size_t{0};

    std::vector<Opm::ParallelWellInfo> pinfos;
    const auto wstate = buildWellState(setup, tstep, pinfos);
    const auto& ws = wstate.well("PROD01");

    BOOST_CHECK_EQUAL(ws.segments.size(), 6);

    const auto& wells = setup.sched.getWellsatEnd();
    BOOST_CHECK_EQUAL(wells.size(), 2);
}

// ---------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(Pressure)
{
    const Setup setup{ "msw.data" };
    const auto tstep = std::size_t{0};

    std::vector<Opm::ParallelWellInfo> pinfos;
    auto wstate = buildWellState(setup, tstep, pinfos);

    const auto& wells = setup.sched.getWells(tstep);
    const auto prod01_first = wells[0].name() == "PROD01";

    setSegPress(wells, wstate);

    const auto rpt = wstate.report(setup.grid.c_grid()->global_cell, [](const int){return false;});

    {
        const auto expect_nSeg = 6;
        const auto& xw = rpt.at("PROD01");

        BOOST_CHECK_EQUAL(xw.segments.size(), expect_nSeg);

        const auto pressTop = prod01_first ? 0.0 : 100.0;

        for (auto segID = 0; segID < expect_nSeg; ++segID) {
            const auto& xseg = xw.segments.at(segID + 1);

            BOOST_CHECK_EQUAL(xseg.segNumber, segID + 1);
            const auto pres_idx = Opm::data::SegmentPressures::Value::Pressure;
            BOOST_CHECK_CLOSE(xseg.pressures[pres_idx], pressTop + 1.0*segID, 1.0e-10);
        }
    }

    const auto& ws = wstate.well("PROD01");
    const auto& perf_data = ws.perf_data;
    (void) perf_data;
}

// ---------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(Rates)
{
    const Setup setup{ "msw.data" };
    const auto tstep = std::size_t{0};

    std::vector<Opm::ParallelWellInfo> pinfos;
    auto wstate = buildWellState(setup, tstep, pinfos);

    const auto wells = setup.sched.getWells(tstep);
    const auto prod01_first = wells[0].name() == "PROD01";

    const auto& pu = setup.pu;

    setSegRates(wells, pu,  wstate);

    const auto rpt = wstate.report(setup.grid.c_grid()->global_cell, [](const int){return false;});

    const auto wat = pu.phase_used[Opm::BlackoilPhases::Aqua];
    const auto oil = pu.phase_used[Opm::BlackoilPhases::Liquid];
    const auto gas = pu.phase_used[Opm::BlackoilPhases::Vapour];

    BOOST_CHECK(wat && oil && gas);

    {
        const auto expect_nSeg = 6;
        const auto& xw = rpt.at("PROD01");

        BOOST_CHECK_EQUAL(xw.segments.size(), expect_nSeg);

        const auto rateTop = prod01_first ? 0.0 : 1000.0;

        for (auto segNum = 1; segNum <= expect_nSeg; ++segNum) {
            const auto& xseg = xw.segments.at(segNum);

            BOOST_CHECK_EQUAL(xseg.segNumber, segNum);

            BOOST_CHECK_CLOSE(xseg.rates.get(Opm::data::Rates::opt::wat),
                              rateTop + 100.0*(segNum - 1), 1.0e-10);

            BOOST_CHECK_CLOSE(xseg.rates.get(Opm::data::Rates::opt::oil),
                              rateTop + 200.0*(segNum - 1), 1.0e-10);

            BOOST_CHECK_CLOSE(xseg.rates.get(Opm::data::Rates::opt::gas),
                              rateTop + 400.0*(segNum - 1), 1.0e-10);
        }
    }
}

// ---------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(STOP_well)
{
    /*
      This test verifies that the perforation pressures is correctly initialized
      also for wells in the STOP state.
    */
    const Setup setup{ "wells_manager_data_wellSTOP.data" };

    std::vector<Opm::ParallelWellInfo> pinfos;
    auto wstate = buildWellState(setup, 0, pinfos);
    for (std::size_t well_index = 0; well_index < setup.sched.numWells(0); well_index++) {
        const auto& ws = wstate.well(well_index);
        const auto& perf_data = ws.perf_data;
        for (const auto& p : perf_data.pressure)
            BOOST_CHECK(p > 0);
    }
}


// ---------------------------------------------------------------------

//BOOST_AUTO_TEST_CASE(GlobalWellInfo_TEST) {
//    const Setup setup{ "msw.data" };
//    std::vector<Opm::Well> local_wells = { setup.sched.getWell("PROD01", 1) };
//    Opm::GlobalWellInfo gwi(setup.sched, 1, local_wells);
//    Opm::WellContainer<Opm::Well::Status> status({{"PROD01", Opm::Well::Status::OPEN}});
//
//    BOOST_CHECK(!gwi.in_injecting_group("INJE01"));
//    BOOST_CHECK(!gwi.in_injecting_group("PROD01"));
//    BOOST_CHECK(!gwi.in_producing_group("INJE01"));
//    BOOST_CHECK(!gwi.in_producing_group("PROD01"));
//
//    BOOST_CHECK_EQUAL( gwi.well_name(0), "INJE01");
//    BOOST_CHECK_EQUAL( gwi.well_name(1), "PROD01");
//    BOOST_CHECK_EQUAL( gwi.well_index("PROD01"), 1);
//
//    BOOST_CHECK_THROW( gwi.update_group( {}, {}, {} ), std::exception);
//
//
//    Opm::WellContainer<Opm::Well::InjectorCMode> inj_cmode({{"PROD01", Opm::Well::InjectorCMode::CMODE_UNDEFINED}});
//    {
//        Opm::WellContainer<Opm::Well::ProducerCMode> prod_cmode({{"PROD01", Opm::Well::ProducerCMode::GRUP}});
//        gwi.update_group(status, inj_cmode, prod_cmode);
//    }
//    BOOST_CHECK(!gwi.in_producing_group("INJE01"));
//    BOOST_CHECK(gwi.in_producing_group("PROD01"));
//
//    {
//        Opm::WellContainer<Opm::Well::ProducerCMode> prod_cmode(
//            {{"PROD01", Opm::Well::ProducerCMode::CMODE_UNDEFINED}});
//        gwi.update_group(status, inj_cmode, prod_cmode);
//    }
//
//    {
//        Opm::WellContainer<Opm::Well::ProducerCMode> prod_cmode({{"PROD01", Opm::Well::ProducerCMode::GRUP}});
//        gwi.update_group(status, inj_cmode, prod_cmode);
//    }
//    BOOST_CHECK(!gwi.in_producing_group("INJE01"));
//    BOOST_CHECK(gwi.in_producing_group("PROD01"));
//
//    {
//        Opm::WellContainer<Opm::Well::ProducerCMode> prod_cmode({{"PROD01", Opm::Well::ProducerCMode::NONE}});
//        gwi.update_group(status, inj_cmode, prod_cmode);
//    }
//    BOOST_CHECK(!gwi.in_producing_group("INJE01"));
//    BOOST_CHECK(!gwi.in_producing_group("PROD01"));
//}

BOOST_AUTO_TEST_CASE(TESTWellContainer) {
    Opm::WellContainer<int> wc;
    BOOST_CHECK_EQUAL(wc.size(), 0);

    wc.add("W1", 1);
    wc.add("W2", 2);

    BOOST_CHECK_EQUAL(wc.size(), 2);

    BOOST_CHECK_THROW(wc.add("W1", 1), std::exception);

    BOOST_CHECK_THROW(wc[10], std::exception);
    BOOST_CHECK_EQUAL(wc[0], 1);
    BOOST_CHECK_EQUAL(wc[1], 2);

    BOOST_CHECK_THROW(wc["INVALID_WELL"], std::exception);
    BOOST_CHECK_EQUAL(wc["W1"], 1);
    BOOST_CHECK_EQUAL(wc["W2"], 2);
    BOOST_CHECK_EQUAL(wc.well_name(0), "W1");
    BOOST_CHECK_EQUAL(wc.well_name(1), "W2");
    BOOST_CHECK_THROW(wc.well_name(10), std::exception);

    const auto& wells = wc.wells();
    std::vector<std::string> expected = {"W1", "W2"};
    BOOST_CHECK( std::is_permutation( wells.begin(), wells.end(), expected.begin(), expected.end()) );

    Opm::WellContainer<int> wc2;
    wc2.copy_welldata(wc);
    BOOST_CHECK_EQUAL(wc2.size() , 0);

    wc2.add("W1", 100);
    BOOST_CHECK_EQUAL(wc2["W1"], 100);
    wc2.copy_welldata(wc);
    BOOST_CHECK_EQUAL(wc2["W1"], 1);

    Opm::WellContainer<int> wc3;
    wc3.add("W2", 100);
    wc3.copy_welldata(wc);
    BOOST_CHECK_EQUAL(wc3["W2"], 2);
    BOOST_CHECK_EQUAL(wc3[0], 2);

    wc3["W2"] = 200;
    wc3.add("W3", 300);
    wc3.copy_welldata(wc);
    BOOST_CHECK_EQUAL(wc3["W2"], 2);
    BOOST_CHECK_EQUAL(wc3[0], 2);
    BOOST_CHECK_EQUAL(wc3["W3"], 300);
    BOOST_CHECK_EQUAL(wc3[1], 300);

    BOOST_CHECK_THROW(wc3.copy_welldata(wc, "W1"), std::exception);
    BOOST_CHECK_THROW(wc3.copy_welldata(wc, "W3"), std::exception);

    wc.clear();
    BOOST_CHECK_EQUAL(wc.size(), 0);


    BOOST_CHECK(wc3.has("W2"));
    BOOST_CHECK(!wc3.has("NO_SUCH_WELL"));


    std::vector<int> vec_copy(wc3.begin(), wc3.end());
    BOOST_CHECK_EQUAL(vec_copy.size(), wc3.size());
    for (std::size_t i = 0; i < wc3.size(); i++)
        BOOST_CHECK_EQUAL(vec_copy[i], wc3[i]);


    Opm::WellContainer<int> wci({{"W1", 1}, {"W2", 2}, {"W3", 3}});
    BOOST_CHECK_EQUAL(wci.size(), 3);
    BOOST_CHECK(wci.has("W1"));
    BOOST_CHECK_EQUAL(wci[1], 2);
    BOOST_CHECK_EQUAL(wci["W3"], 3);

    auto w3 = wci.well_index("W3");
    BOOST_CHECK(w3.has_value());
    BOOST_CHECK_EQUAL(w3.value(), 2);

    auto wx = wci.well_index("WX");
    BOOST_CHECK(!wx.has_value());
}

BOOST_AUTO_TEST_CASE(TESTSegmentState) {
    const Setup setup{ "msw.data" };
    const auto& well = setup.sched.getWell("PROD01", 0);
    const auto& segments = well.getSegments();
    Opm::SegmentState ss1(3, segments);
    Opm::SegmentState ss2;


    ss1.pressure_drop_hydrostatic[0] = 1;
    ss1.pressure_drop_friction[0] = 2;
    ss1.pressure_drop_accel[0] = 3;

    BOOST_CHECK_EQUAL(ss1.pressure_drop(0), 6);
}

BOOST_AUTO_TEST_CASE(TESTSegmentState2) {
    const Setup setup{ "msw.data" };
    std::vector<Opm::ParallelWellInfo> pinfo;
    const auto wstate = buildWellState(setup, 0, pinfo);
    const auto& well = setup.sched.getWell("PROD01", 0);
    const auto& ws = wstate.well("PROD01");

    auto segments = ws.segments;
    BOOST_CHECK_EQUAL(segments.pressure.size(), well.getSegments().size());

    segments.pressure_drop_friction[0] = 1;
    segments.pressure_drop_accel[0] = 2;
    segments.pressure_drop_hydrostatic[0] = 4;

    BOOST_CHECK_EQUAL(segments.pressure_drop(0), 7);


    for (std::size_t i=0; i < segments.pressure.size(); i++)
        segments.pressure[i] = (i + 1);

    const double bhp = 2.0;
    segments.scale_pressure(bhp);

    for (std::size_t i=0; i < segments.pressure.size(); i++)
        BOOST_CHECK_EQUAL(segments.pressure[i], 1.+(i+1));

    BOOST_CHECK_EQUAL( segments.size(), well.getSegments().size() );
}


BOOST_AUTO_TEST_CASE(TESTPerfData) {
    Opm::PerfData pd1(3, 100, true, 3);
    Opm::PerfData pd2(3, 100, true, 3);
    Opm::PerfData pd3(2, 100, true, 3);
    Opm::PerfData pd4(3, 100, false, 3);


    for (std::size_t i = 0; i < 3; i++) {
        pd1.pressure[i] = i+1;
        pd2.pressure[i] = 10*(i+1);
    }

    BOOST_CHECK(pd1.try_assign(pd2));
    for (std::size_t i = 0; i < 3; i++) {
        BOOST_CHECK(pd2.pressure[i] == 10*(i+1));
        BOOST_CHECK(pd1.pressure[i] == 10*(i+1));
    }

    BOOST_CHECK(!pd1.try_assign(pd3));
    for (std::size_t i = 0; i < 3; i++) {
        BOOST_CHECK(pd1.pressure[i] == 10*(i+1));
    }

    BOOST_CHECK(!pd1.try_assign(pd4));
}


BOOST_AUTO_TEST_CASE(TestSingleWellState) {
    Opm::ParallelWellInfo pinfo;
    std::vector<Opm::PerforationData> connections = {{0,1,1,0,0},{1,1,1,0,1},{2,1,1,0,2}};
    Opm::PhaseUsage pu;

    // This is totally bonkers, but the pu needs a complete deck to initialize properly
    pu.num_phases = 3;

    Opm::SingleWellState ws1("W1", pinfo, true,  100, connections, pu, 1);
    Opm::SingleWellState ws2("W2", pinfo, true,  100, connections, pu, 2);
    Opm::SingleWellState ws3("W3", pinfo, false, 100, connections, pu, 3);

    ws1.bhp = 100;
    ws1.thp = 200;

    ws2.init_timestep(ws1);
    BOOST_CHECK_EQUAL(ws2.bhp, ws1.bhp);
    BOOST_CHECK_EQUAL(ws2.thp, ws1.thp);

    ws3.bhp = ws1.bhp * 2;
    ws3.thp = ws1.thp * 2;
    ws3.init_timestep(ws1);
    BOOST_CHECK(ws3.bhp != ws1.bhp);
    BOOST_CHECK(ws3.thp != ws1.thp);
}



BOOST_AUTO_TEST_CASE(TestPU) {
    Opm::PhaseUsage pu({Opm::BlackoilPhases::Polymer, Opm::BlackoilPhases::Solvent, Opm::BlackoilPhases::Aqua, Opm::BlackoilPhases::ZFraction});

    BOOST_CHECK(pu.phase_used[Opm::BlackoilPhases::Aqua]);      BOOST_CHECK_EQUAL(pu.phase_pos[Opm::BlackoilPhases::Aqua], 0);
    BOOST_CHECK(pu.phase_used[Opm::BlackoilPhases::Solvent]);   BOOST_CHECK_EQUAL(pu.phase_pos[Opm::BlackoilPhases::Solvent], 1);
    BOOST_CHECK(pu.phase_used[Opm::BlackoilPhases::Polymer]);   BOOST_CHECK_EQUAL(pu.phase_pos[Opm::BlackoilPhases::Polymer], 2);
    BOOST_CHECK(pu.phase_used[Opm::BlackoilPhases::ZFraction]); BOOST_CHECK_EQUAL(pu.phase_pos[Opm::BlackoilPhases::ZFraction], 3);

    BOOST_CHECK(!pu.phase_used[Opm::BlackoilPhases::Liquid]);   BOOST_CHECK_EQUAL(pu.phase_pos[Opm::BlackoilPhases::Liquid], -1);
    BOOST_CHECK(!pu.phase_used[Opm::BlackoilPhases::Energy]);   BOOST_CHECK_EQUAL(pu.phase_pos[Opm::BlackoilPhases::Energy], -1);

    BOOST_CHECK_EQUAL(pu.num_phases, 1);  // Only Aqua counts as a phase.
    BOOST_CHECK(pu.has_polymer);
    BOOST_CHECK(pu.has_solvent);
    BOOST_CHECK(pu.has_zFraction);

    BOOST_CHECK(!pu.has_energy);
    BOOST_CHECK(!pu.has_polymermw);
    BOOST_CHECK(!pu.has_foam);
    BOOST_CHECK(!pu.has_brine);
}


BOOST_AUTO_TEST_SUITE_END()
