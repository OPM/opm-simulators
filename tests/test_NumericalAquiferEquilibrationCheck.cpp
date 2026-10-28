/*
  Copyright 2026 Equinor ASA.

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

#define BOOST_TEST_MODULE TestNumericalAquiferEquilibrationCheck

#define BOOST_TEST_NO_MAIN

#include <boost/test/unit_test.hpp>

#include <opm/simulators/flow/NumericalAquiferEquilibrationCheck.hpp>

#include <opm/simulators/utils/ParallelCommunication.hpp>

#include <opm/common/OpmLog/OpmLog.hpp>
#include <opm/common/OpmLog/StreamLog.hpp>

#include <opm/input/eclipse/Deck/Deck.hpp>
#include <opm/input/eclipse/EclipseState/Aquifer/AquiferConfig.hpp>
#include <opm/input/eclipse/EclipseState/Aquifer/NumericalAquifer/NumericalAquifers.hpp>
#include <opm/input/eclipse/EclipseState/EclipseState.hpp>
#include <opm/input/eclipse/Parser/Parser.hpp>

#include <dune/common/parallel/mpihelper.hh>

#include <cstddef>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

// Aquifer 1, in EQLNUM region 2, connects to cells of regions 1 (K=1-2)
// and 2 (K=4-6).  Row J=1 and layer K=3 are inactive.
Opm::EclipseState makeEclipseState()
{
    const auto deck = Opm::Parser{}.parseString(R"(
RUNSPEC
DIMENS
 4 3 6 /
OIL
WATER
METRIC
EQLDIMS
 2 /
GRID
DX
 72*100 /
DY
 72*100 /
DZ
 72*5 /
TOPS
 12*2000 /
PORO
 72*0.25 /
PERMX
 72*100 /
PERMY
 72*100 /
PERMZ
 72*10 /
ACTNUM
 4*0 8*1  4*0 8*1  12*0  4*0 8*1  4*0 8*1  4*0 8*1 /
AQUNUM
 1 1 1 5 1.0E6 1000 0.25 500 /
/
AQUCON
 1 1 4 2 2 1 6 'J-' /
/
REGIONS
EQLNUM
 24*1 48*2 /
)");

    auto eclState = Opm::EclipseState { deck };
    const auto& grid = eclState.getInputGrid();
    eclState.aquifer().mutableNumericalAquifers().postProcessConnections(grid, grid.getACTNUM());

    return eclState;
}

Opm::Parallel::Communication worldComm()
{
    return Opm::Parallel::Communication { Dune::MPIHelper::getCommunicator() };
}

// Captures messages logged through OpmLog while in scope.
class LogCapture
{
public:
    LogCapture()
    {
        Opm::OpmLog::addBackend(name_, std::make_shared<Opm::StreamLog>(stream_, Opm::Log::DefaultMessageTypes));
    }

    ~LogCapture()
    {
        Opm::OpmLog::removeBackend(name_);
    }

    LogCapture(const LogCapture&) = delete;
    LogCapture& operator=(const LogCapture&) = delete;

    std::string text() const
    {
        return stream_.str();
    }

private:
    static constexpr const char* name_ = "NumAquEquilCheck";
    std::ostringstream stream_;
};

// Mismatches reported by the last process only, i.e. not by the I/O
// process whenever there is more than one process.
std::vector<Opm::NumericalAquiferEquilibrationMismatch> lastRankMismatches()
{
    const auto comm = worldComm();
    if (comm.rank() != comm.size() - 1) {
        return {};
    }

    return {
        {1, 2, 1, 2.0e6},
        {1, 2, 1, -5.0e5},
    };
}

bool init_unit_test_func()
{
    return true;
}

} // Anonymous namespace

BOOST_AUTO_TEST_CASE(EmptyInput)
{
    const auto eclState = makeEclipseState();
    const LogCapture log;

    std::size_t count = 1;
    BOOST_CHECK_NO_THROW(count = Opm::checkNumericalAquiferEquilibration({}, eclState, worldComm(),
                                                                         /*allowMismatch=*/false));
    BOOST_CHECK_EQUAL(count, std::size_t{0});
    BOOST_CHECK(log.text().empty());
}

BOOST_AUTO_TEST_CASE(MismatchOnLastRankAllowed)
{
    const auto eclState = makeEclipseState();
    const auto comm = worldComm();
    const LogCapture log;

    std::size_t count = 0;
    BOOST_CHECK_NO_THROW(count = Opm::checkNumericalAquiferEquilibration(lastRankMismatches(), eclState,
                                                                         comm, /*allowMismatch=*/true));

    // Same global count on all processes.
    BOOST_CHECK_EQUAL(count, std::size_t{2});
    BOOST_CHECK_EQUAL(comm.min(count), comm.max(count));

    const auto text = log.text();
    if (comm.rank() == 0) {
        BOOST_CHECK(text.find("Aquifer 1 (cells in EQLNUM region 2): 2 connection(s) "
                              "to cells in EQLNUM region 1") != std::string::npos);
        BOOST_CHECK(text.find("up to 20 BARS") != std::string::npos);

        // Deck changes: new aquifer 2 for the region 1 connections.
        BOOST_CHECK(text.find(" 2 1 4 2 2 1 2 'J-' 1 0 'NO' 1 1 /") != std::string::npos);
        BOOST_CHECK(text.find(" 1 1 4 2 2 4 6 'J-' 1 0 'NO' 1 1 /") != std::string::npos);
        BOOST_CHECK(text.find("'EQLNUM' 1") != std::string::npos);

        BOOST_CHECK(text.find("--allow-inconsistent-aquifer-equilibration=true") != std::string::npos);
        BOOST_CHECK(text.find("This option is set, so the simulation continues.") != std::string::npos);
    }
    else {
        BOOST_CHECK(text.empty());
    }
}

BOOST_AUTO_TEST_CASE(MismatchOnLastRankStopsAllRanks)
{
    const auto eclState = makeEclipseState();
    const auto comm = worldComm();
    const LogCapture log;

    // Every process must throw, also those which found no mismatch.
    BOOST_CHECK_THROW(Opm::checkNumericalAquiferEquilibration(lastRankMismatches(), eclState,
                                                              comm, /*allowMismatch=*/false),
                      std::domain_error);

    if (comm.rank() == 0) {
        const auto text = log.text();
        BOOST_CHECK(text.find("Numerical aquifer(s) 1 connect cells") != std::string::npos);
        BOOST_CHECK(text.find("--allow-inconsistent-aquifer-equilibration=true") != std::string::npos);
        BOOST_CHECK(text.find("This option is set") == std::string::npos);
    }
}

int main(int argc, char** argv)
{
    Dune::MPIHelper::instance(argc, argv);

    return boost::unit_test::unit_test_main(&init_unit_test_func, argc, argv);
}
