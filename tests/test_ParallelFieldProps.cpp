/*
  Copyright 2026 SINTEF Digital

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

#define BOOST_TEST_MODULE ParallelFieldProps
#include <boost/test/unit_test.hpp>

#include "MpiFixture.hpp"

#include <opm/simulators/utils/ParallelCommunication.hpp>
#include <opm/simulators/utils/ParallelEclipseState.hpp>

#include <opm/input/eclipse/Deck/Deck.hpp>
#include <opm/input/eclipse/EclipseState/Grid/EclipseGrid.hpp>
#include <opm/input/eclipse/EclipseState/Grid/FieldData.hpp>
#include <opm/input/eclipse/EclipseState/Grid/FieldPropsManager.hpp>
#include <opm/input/eclipse/EclipseState/Runspec.hpp>
#include <opm/input/eclipse/EclipseState/Tables/TableManager.hpp>
#include <opm/input/eclipse/Parser/Parser.hpp>

#include <stdexcept>
#include <string>
#include <vector>

BOOST_GLOBAL_FIXTURE(MPIFixture);

namespace
{

// Use reordered local indices, including the last Cartesian cell. Cell 5
// appears on both ranks to model overlap.
struct LocalMapper {
    std::vector<int> cells;

    int compressedSize(int) const
    {
        return static_cast<int>(cells.size());
    }
    int cartesianIndex(int localIdx, int) const
    {
        return cells[localIdx];
    }
};

struct Setup {
    explicit Setup(const bool inactiveCells)
        : local(global, comm)
        , mapper {comm.rank() == 0 ? std::vector<int> {5, 0, 3} : std::vector<int> {2, 5}}
    {
        if (comm.rank() == 0) {
            const auto deck = Opm::Parser {}.parseString(R"(
RUNSPEC
DIMENS
6 1 1 /
COMPS
3 /
OIL
GAS
GRID
PORO
0.10 0.11 0.12 0.13 0.14 0.15 /
SOLUTION
ZMF
0.10 0.11 0.12 0.13 0.14 0.15
0.20 0.21 0.22 0.23 0.24 0.25
0.70 0.68 0.66 0.64 0.62 0.60 /
XMF
0.20 0.21 0.22 0.23 0.24 0.25
0.30 0.31 0.32 0.33 0.34 0.35
0.50 0.48 0.46 0.44 0.42 0.40 /
YMF
0.40 0.41 0.42 0.43 0.44 0.45
0.50 0.49 0.48 0.47 0.46 0.45
0.10 0.10 0.10 0.10 0.10 0.10 /
)");
            global
                = Opm::FieldPropsManager {deck, Opm::Phases {true, true, false}, grid, tables, 3};
            if (inactiveCells) {
                // Simulate cells removed by grid processing after parsing.
                grid.resetACTNUM({1, 0, 1, 1, 0, 1});
                global.reset_actnum(grid.getACTNUM());
            }
        }

        // Leave the local cache empty so each property's first get_double()
        // call fetches it from rank zero.
        local.resetCartesianMapper(&mapper);
    }

    Opm::Parallel::Communication comm;
    Opm::EclipseGrid grid {6, 1, 1};
    Opm::TableManager tables;
    Opm::FieldPropsManager global;
    Opm::ParallelFieldPropsManager local;
    LocalMapper mapper;
};

} // namespace

BOOST_AUTO_TEST_CASE(RejectUndistributedMultiValuedProperties)
{
    Setup setup(false);
    for (const std::string keyword : {"ZMF", "XMF", "YMF"}) {
        BOOST_CHECK_THROW(setup.local.get_global_double(keyword), std::runtime_error);
        BOOST_CHECK_THROW(setup.local.get_double(keyword), std::runtime_error);
        BOOST_CHECK(!setup.local.has_double(keyword));
    }
}

BOOST_AUTO_TEST_CASE(ScalarAndMissingProperties)
{
    Setup setup(true);
    const std::vector<double> expectedGlobal {0.10, 0.0, 0.12, 0.13, 0.0, 0.15};
    const auto global = setup.local.get_global_double("PORO");
    BOOST_CHECK_EQUAL_COLLECTIONS(
        global.begin(), global.end(), expectedGlobal.begin(), expectedGlobal.end());
    const auto& local = setup.local.get_double("PORO");
    const std::vector<double> expectedLocal = setup.comm.rank() == 0
        ? std::vector<double> {0.15, 0.10, 0.13}
        : std::vector<double> {0.12, 0.15};
    BOOST_CHECK_EQUAL_COLLECTIONS(
        local.begin(), local.end(), expectedLocal.begin(), expectedLocal.end());
    const auto& ntg = setup.local.get_double("NTG");
    const std::vector<double> expectedNtg(setup.mapper.cells.size(), 1.0);
    BOOST_CHECK_EQUAL_COLLECTIONS(ntg.begin(), ntg.end(), expectedNtg.begin(), expectedNtg.end());

    // Ranks without parsed field data must receive the same failure.
    BOOST_CHECK_THROW(setup.local.get_double("NO_SUCH_PROPERTY"), std::exception);
    BOOST_CHECK(!setup.local.has_double("NO_SUCH_PROPERTY"));
}
