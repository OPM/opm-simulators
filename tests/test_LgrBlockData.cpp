/*
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

#define BOOST_TEST_MODULE LgrBlockDataGather

#include <boost/test/unit_test.hpp>

#include <opm/simulators/flow/CollectDataOnIORank_impl.hpp>

#include <map>
#include <string>
#include <tuple>

namespace {
    // Key = (summary keyword, LGR level, level-local linearised Cartesian index),
    // matching opm-common's Summary::DynamicSimulatorState::LgrBlockValues.
    using Key = std::tuple<std::string, int, int>;
    using LgrBlockMap = std::map<Key, double>;
}

// Opm::PackUnPackLgrBlockData is the gather handle for the tuple-keyed LB* per-LGR-cell
// summary map. Its IO-rank constructor self-packs the local map into a message buffer and
// unpacks it back into the supplied global map -- which is exactly how the IO rank folds in
// each rank's contribution during collect(). Constructing it with isIORank == true therefore
// exercises the pack -> unpack round-trip, and constructing it twice into the same global map
// exercises the disjoint union of two ranks' contributions -- all without needing MPI.

BOOST_AUTO_TEST_CASE(RoundTripPreservesEntries)
{
    const LgrBlockMap local = {
        { Key{"LBPR",   1, 0}, 250.0 },
        { Key{"LBPR",   1, 3}, 251.5 },
        { Key{"LBOSAT", 2, 0}, 0.30  },
    };

    LgrBlockMap gathered;
    Opm::PackUnPackLgrBlockData{ local, gathered, /*isIORank=*/true };

    BOOST_CHECK(gathered == local);
}

BOOST_AUTO_TEST_CASE(DisjointMergeAccumulatesBothRanks)
{
    const LgrBlockMap rank0 = {
        { Key{"LBPR", 1, 0}, 250.0 },
        { Key{"LBPR", 1, 3}, 251.5 },
    };
    const LgrBlockMap rank1 = {
        { Key{"LBPR",   1, 7}, 252.0 },
        { Key{"LBOSAT", 2, 0}, 0.30  },
    };

    LgrBlockMap gathered;
    Opm::PackUnPackLgrBlockData{ rank0, gathered, /*isIORank=*/true };
    Opm::PackUnPackLgrBlockData{ rank1, gathered, /*isIORank=*/true };

    LgrBlockMap expected = rank0;
    expected.insert(rank1.begin(), rank1.end());

    BOOST_CHECK(gathered == expected);
    BOOST_CHECK_EQUAL(gathered.size(), 4u);
    BOOST_CHECK_EQUAL(gathered.at(Key{"LBPR",   1, 0}), 250.0);
    BOOST_CHECK_EQUAL(gathered.at(Key{"LBPR",   1, 7}), 252.0);
    BOOST_CHECK_EQUAL(gathered.at(Key{"LBOSAT", 2, 0}), 0.30);
}

BOOST_AUTO_TEST_CASE(EmptyMapGathersNothing)
{
    const LgrBlockMap empty;

    LgrBlockMap gathered;
    Opm::PackUnPackLgrBlockData{ empty, gathered, /*isIORank=*/true };

    BOOST_CHECK(gathered.empty());
}
