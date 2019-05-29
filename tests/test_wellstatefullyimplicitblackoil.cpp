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

#include <config.h>

#define BOOST_TEST_MODULE WellStateFIBOTest

#include <opm/simulators/wells/WellStateFullyImplicitBlackoil.hpp>

#include <boost/test/unit_test.hpp>

#include <opm/parser/eclipse/Parser/Parser.hpp>
#include <opm/parser/eclipse/Parser/ParseContext.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/Schedule.hpp>
#include <opm/parser/eclipse/Units/Units.hpp>

#include <opm/core/props/BlackoilPhases.hpp>
#include <opm/core/props/phaseUsageFromDeck.hpp>
#include <opm/core/wells/WellsManager.hpp>
#include <opm/core/wells.h>

#include <opm/grid/GridManager.hpp>

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
        , st()
    {}

    Opm::EclipseState es;
    Opm::PhaseUsage   pu;
    Opm::GridManager  grid;
    Opm::Schedule     sched;
    Opm::SummaryState st;
};

namespace {
    Opm::WellStateFullyImplicitBlackoil
    buildWellState(const Setup& setup, const std::size_t timeStep)
    {
        auto state  = Opm::WellStateFullyImplicitBlackoil{};

        const auto cpress =
            std::vector<double>(setup.grid.c_grid()->number_of_cells,
                                100.0*Opm::unit::barsa);

        const Opm::WellsManager wmgr{setup.es, setup.sched, setup.st, timeStep, *setup.grid.c_grid()};

        state.init(wmgr.c_wells(), cpress, setup.sched,
                   setup.sched.getWells2(timeStep),
                   timeStep, nullptr, setup.pu);

        state.initWellStateMSWell(wmgr.c_wells(),
                                  setup.sched.getWells2(timeStep),
                                  timeStep, setup.pu, nullptr);

        return state;
    }


    void setSegPress(const std::vector<Opm::Well2>& wells,
                     const std::size_t                    tstep,
                     Opm::WellStateFullyImplicitBlackoil& wstate)
    {
        const auto nWell = wells.size();

        auto& segPress = wstate.segPress();

        for (auto wellID = 0*nWell; wellID < nWell; ++wellID) {
            const auto& well     = wells[wellID];
            const auto  topSegIx = wstate.topSegmentIndex(wellID);
            const auto  pressTop = 100.0 * wellID;

            auto* press = &segPress[topSegIx];

            press[0] = pressTop;

            if (! well.isMultiSegment()) {
                continue;
            }

            const auto& segSet = well.getSegments();
            const auto  nSeg   = segSet.size();

            for (auto segID = 0*nSeg + 1; segID < nSeg; ++segID) {
                // One-based numbering scheme for segments.
                const auto segNo = segSet[segID].segmentNumber();
                press[segNo - 1] = pressTop + 1.0*(segNo - 1);
            }
        }
    }


  void setSegRates(const std::vector<Opm::Well2>& wells,
                     const std::size_t                    tstep,
                     const Opm::PhaseUsage&               pu,
                     Opm::WellStateFullyImplicitBlackoil& wstate)
    {
        const auto wat = pu.phase_used[Opm::BlackoilPhases::Aqua];
        const auto iw  = wat ? pu.phase_pos[Opm::BlackoilPhases::Aqua] : -1;

        const auto oil = pu.phase_used[Opm::BlackoilPhases::Liquid];
        const auto io  = oil ? pu.phase_pos[Opm::BlackoilPhases::Liquid] : -1;

        const auto gas = pu.phase_used[Opm::BlackoilPhases::Vapour];
        const auto ig  = gas ? pu.phase_pos[Opm::BlackoilPhases::Vapour] : -1;

        const auto np = wstate.numPhases();

        const auto nWell = wells.size();

        auto& segRates = wstate.segRates();

        for (auto wellID = 0*nWell; wellID < nWell; ++wellID) {
            const auto& well     = wells[wellID];
            const auto  topSegIx = wstate.topSegmentIndex(wellID);
            const auto  rateTop  = 1000.0 * wellID;

            if (wat) { segRates[np*topSegIx + iw] = rateTop; }
            if (oil) { segRates[np*topSegIx + io] = rateTop; }
            if (gas) { segRates[np*topSegIx + ig] = rateTop; }

            if (! well.isMultiSegment()) {
                continue;
            }

            const auto& segSet = well.getSegments();
            const auto  nSeg   = segSet.size();

            for (auto segID = 0*nSeg + 1; segID < nSeg; ++segID) {
                // One-based numbering scheme for segments.
                const auto segNo = segSet[segID].segmentNumber();

                auto* rates = &segRates[(topSegIx + segNo - 1) * np];

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

    BOOST_CHECK_EQUAL(wstate.numSegment(), 6 + 1);

    const auto& wells = setup.sched.getWells2atEnd();
    BOOST_CHECK_EQUAL(wells.size(), 2);

    const auto prod01_first = wells[0].name() == "PROD01";

    BOOST_CHECK_EQUAL(wstate.topSegmentIndex(0), 0);
    BOOST_CHECK_EQUAL(wstate.topSegmentIndex(1),
                      prod01_first ? 6 : 1);
}

// ---------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(Pressure)
{
    const Setup setup{ "msw.data" };
    const auto tstep = std::size_t{0};

    auto wstate = buildWellState(setup, tstep);

    const auto& wells = setup.sched.getWells2(tstep);
    const auto prod01_first = wells[0].name() == "PROD01";

    setSegPress(wells, tstep, wstate);

    const auto rpt = wstate.report(setup.pu, setup.grid.c_grid()->global_cell);

    {
        const auto& xw = rpt.at("INJE01");

        BOOST_CHECK_EQUAL(xw.segments.size(), 1); // Top Segment

        const auto& xseg = xw.segments.at(1);

        BOOST_CHECK_EQUAL(xseg.segNumber, 1);
        BOOST_CHECK_CLOSE(xseg.pressure, prod01_first ? 100.0 : 0.0, 1.0e-10);
    }

    {
        const auto expect_nSeg = 6;
        const auto& xw = rpt.at("PROD01");

        BOOST_CHECK_EQUAL(xw.segments.size(), expect_nSeg);

        const auto pressTop = prod01_first ? 0.0 : 100.0;

        for (auto segID = 0; segID < expect_nSeg; ++segID) {
            const auto& xseg = xw.segments.at(segID + 1);

            BOOST_CHECK_EQUAL(xseg.segNumber, segID + 1);
            BOOST_CHECK_CLOSE(xseg.pressure, pressTop + 1.0*segID, 1.0e-10);
        }
    }
}

// ---------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(Rates)
{
    const Setup setup{ "msw.data" };
    const auto tstep = std::size_t{0};

    auto wstate = buildWellState(setup, tstep);

    const auto wells = setup.sched.getWells2(tstep);
    const auto prod01_first = wells[0].name() == "PROD01";

    const auto& pu = setup.pu;

    setSegRates(wells, tstep, pu, wstate);

    const auto rpt = wstate.report(pu, setup.grid.c_grid()->global_cell);

    const auto wat = pu.phase_used[Opm::BlackoilPhases::Aqua];
    const auto oil = pu.phase_used[Opm::BlackoilPhases::Liquid];
    const auto gas = pu.phase_used[Opm::BlackoilPhases::Vapour];

    BOOST_CHECK(wat && oil && gas);

    {
        const auto rateTop = prod01_first ? 1000.0 : 0.0;

        const auto& xw = rpt.at("INJE01");

        BOOST_CHECK_EQUAL(xw.segments.size(), 1); // Top Segment

        const auto& xseg = xw.segments.at(1);

        BOOST_CHECK_EQUAL(xseg.segNumber, 1);
        BOOST_CHECK_CLOSE(xseg.rates.get(Opm::data::Rates::opt::wat),
                          rateTop, 1.0e-10);

        BOOST_CHECK_CLOSE(xseg.rates.get(Opm::data::Rates::opt::oil),
                          rateTop, 1.0e-10);

        BOOST_CHECK_CLOSE(xseg.rates.get(Opm::data::Rates::opt::gas),
                          rateTop, 1.0e-10);
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

BOOST_AUTO_TEST_SUITE_END()
