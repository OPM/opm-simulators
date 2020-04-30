/*
  Copyright 2018, 2020 Equinor ASA.

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

#include <config.h>

#define BOOST_TEST_MODULE WellAndGroupStatesTest

#include <boost/test/unit_test.hpp>

#include <opm/simulators/wells/WellAndGroupStates.hpp>

#include <opm/simulators/wells/PerforationData.hpp>

#include <opm/parser/eclipse/Parser/Parser.hpp>
#include <opm/parser/eclipse/Parser/ParseContext.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/Schedule.hpp>
#include <opm/parser/eclipse/Units/Units.hpp>

#include <opm/grid/GridHelpers.hpp>

#include <opm/core/props/BlackoilPhases.hpp>
#include <opm/core/props/phaseUsageFromDeck.hpp>

#include <opm/grid/GridManager.hpp>

#include <chrono>
#include <cstddef>
#include <string>

struct Setup
{
    Setup(const std::string& filename)
        : Setup(Opm::Parser{}.parseFile(filename))
    {}

    Setup(const Opm::Deck& deck)
        : es   (deck)
        , pu   (Opm::phaseUsageFromDeck(es))
        , grid (es.getInputGrid())
        , sched(deck, es)
        , st(std::chrono::system_clock::from_time_t(sched.getStartTime()))
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
                        pd.satnum_id = completion.satTableId();
                        well_perf_data[well_index].push_back(pd);
                    }
                } else {
                    if (completion.state() != Opm::Connection::State::SHUT) {
                        OPM_THROW(std::runtime_error,
                                  "Completion state: " << Opm::Connection::State2String(completion.state()) << " not handled");
                    }
                }
            }
            ++well_index;
        }
    }

    Opm::EclipseState es;
    Opm::PhaseUsage   pu;
    Opm::GridManager  grid;
    Opm::Schedule     sched;
    Opm::SummaryState st;
    std::vector<std::vector<Opm::PerforationData>> well_perf_data;
};

namespace {
    Opm::WellAndGroupStates<3>
    buildWellState(const Setup& setup, const std::size_t timeStep)
    {
        auto state  = Opm::WellAndGroupStates<3>{};

        const auto cpress =
            std::vector<double>(setup.grid.c_grid()->number_of_cells,
                                100.0*Opm::unit::barsa);

        state.init(cpress, setup.sched,
                   setup.sched.getWells(timeStep),
                   timeStep, nullptr, setup.pu, setup.well_perf_data, setup.st);

        // TODO: ensure this is obsolete
        // state.initWellStateMSWell(setup.sched.getWells(timeStep),
        //                           setup.pu, nullptr);

        return state;
    }


    void setSegPressWell(const Opm::Well& well,
                         const int wellID,
                         Opm::SingleWellState<3>& wstate)
    {
        if (! well.isMultiSegment()) {
            return;
        }

        const auto pressTop = 100.0 * wellID;
        wstate.segments[0].pressure = pressTop;

        const auto& segSet = well.getSegments();
        const auto  nSeg   = segSet.size();
        for (auto segID = 0*nSeg + 1; segID < nSeg; ++segID) {
            // One-based numbering scheme for segments.
            const auto segNo = segSet[segID].segmentNumber();
            wstate.segments[segNo - 1].pressure = pressTop + 1.0*(segNo - 1);
        }
    }


    void setSegPress(const std::vector<Opm::Well>& wells,
                     Opm::WellAndGroupStates<3>& wstate)
    {
        const auto nWell = wells.size();

        for (auto wellID = 0*nWell; wellID < nWell; ++wellID) {
            const auto& well = wells[wellID];
            auto& state = wstate.wellStates()[wellID];
            setSegPressWell(well, wellID, state);
        }
    }


    void setSegRates(const std::vector<Opm::Well>& wells,
                     const Opm::PhaseUsage& pu,
                     Opm::WellAndGroupStates<3>& wstate)
    {
        // This is not a very instructive example, as the well and group
        // state objects should no longer use the PhaseUsage-derived orderings.
        const auto wat = pu.phase_used[Opm::BlackoilPhases::Aqua];
        const auto iw  = wat ? pu.phase_pos[Opm::BlackoilPhases::Aqua] : -1;

        const auto oil = pu.phase_used[Opm::BlackoilPhases::Liquid];
        const auto io  = oil ? pu.phase_pos[Opm::BlackoilPhases::Liquid] : -1;

        const auto gas = pu.phase_used[Opm::BlackoilPhases::Vapour];
        const auto ig  = gas ? pu.phase_pos[Opm::BlackoilPhases::Vapour] : -1;

        const auto nWell = wells.size();
        for (auto wellID = 0*nWell; wellID < nWell; ++wellID) {
            const auto& well = wells[wellID];
            if (! well.isMultiSegment()) {
                continue;
            }

            auto& state = wstate.wellStates()[wellID];
            const auto rateTop = 1000.0 * wellID;
            {
                auto& rates = state.segments[0].surface_rates;
                if (wat) { rates[iw] = rateTop; }
                if (oil) { rates[io] = rateTop; }
                if (gas) { rates[ig] = rateTop; }
            }

            const auto& segSet = well.getSegments();
            const auto  nSeg   = segSet.size();
            for (auto segID = 0*nSeg + 1; segID < nSeg; ++segID) {
                // One-based numbering scheme for segments.
                const auto segNo = segSet[segID].segmentNumber();
                auto& rates = state.segments[segNo - 1].surface_rates;
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

    const auto wstate = buildWellState(setup, tstep);

    BOOST_CHECK_EQUAL(wstate.wellStates()[0].segments.size(), 0);
    BOOST_CHECK_EQUAL(wstate.wellStates()[1].segments.size(), 6);

    const auto& wells = setup.sched.getWellsatEnd();
    BOOST_CHECK_EQUAL(wells.size(), 2);
}

// ---------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(Pressure)
{
    const Setup setup{ "msw.data" };
    const auto tstep = std::size_t{0};

    auto wstate = buildWellState(setup, tstep);

    const auto& wells = setup.sched.getWells(tstep);
    const auto prod01_first = wells[0].name() == "PROD01";

    setSegPress(wells, wstate);

    const auto rpt = wstate.report(setup.pu, setup.grid.c_grid()->global_cell);

    {
        const auto& xw = rpt.at("INJE01");

        BOOST_CHECK(xw.segments.empty()); // Not a multisegment well.
    }

    {
        const auto expect_nSeg = 6;
        const auto& xw = rpt.at("PROD01");

        BOOST_CHECK_EQUAL(xw.segments.size(), expect_nSeg);

        const auto pressTop = prod01_first ? 0.0 : 100.0;

        for (auto segID = 0; segID < expect_nSeg; ++segID) {
            const auto& xseg = xw.segments.at(segID + 1);

            BOOST_CHECK_EQUAL(xseg.segNumber, segID + 1);
            BOOST_CHECK_CLOSE(xseg.pressures[Opm::data::SegmentPressures::Value::Pressure],
                              pressTop + 1.0*segID, 1.0e-10);
        }
    }
}

// ---------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(Rates)
{
    const Setup setup{ "msw.data" };
    const auto tstep = std::size_t{0};

    auto wstate = buildWellState(setup, tstep);

    const auto wells = setup.sched.getWells(tstep);
    const auto prod01_first = wells[0].name() == "PROD01";

    const auto& pu = setup.pu;

    setSegRates(wells, pu, wstate);

    const auto rpt = wstate.report(pu, setup.grid.c_grid()->global_cell);

    const auto wat = pu.phase_used[Opm::BlackoilPhases::Aqua];
    const auto oil = pu.phase_used[Opm::BlackoilPhases::Liquid];
    const auto gas = pu.phase_used[Opm::BlackoilPhases::Vapour];

    BOOST_CHECK(wat && oil && gas);

    {
        const auto rateTop = prod01_first ? 1000.0 : 0.0;

        const auto& xw = rpt.at("INJE01");

        BOOST_CHECK(xw.segments.empty()); // Not a multisegment well.
    }

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
    auto wstate = buildWellState(setup, 0);
    for (const auto& ws : wstate.wellStates()) {
        for (const auto& conn : ws.connections) {
            BOOST_CHECK(conn.pressure > 0.0);
        }
    }
}

BOOST_AUTO_TEST_SUITE_END()
