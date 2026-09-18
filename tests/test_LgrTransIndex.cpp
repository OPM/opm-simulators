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

#define BOOST_TEST_MODULE LgrTransIndex

#include <boost/test/unit_test.hpp>

#include <opm/simulators/flow/LgrOutputTransGather.hpp>

#include <array>
#include <utility>
#include <vector>

// Unit tests for Opm::LgrTransIndex<N> -- the flat open-addressing table that holds the
// gathered LGR output transmissibilities on the I/O rank (see ADR-0001). The parallel INIT
// regression exercises it end to end; these tests pin the container in isolation: empty
// state (the non-I/O-rank case), exact value round-trip, absent-key termination (no infinite
// probe), N=3 and N=4 key widths, negative/large key components, and sizes straddling the
// power-of-two capacity boundary.

using Opm::LgrTransIndex;

namespace {
    using Key3 = std::array<int, 3>;
    using Rec3 = std::pair<Key3, double>;
    using Key4 = std::array<int, 4>;
    using Rec4 = std::pair<Key4, double>;
}

BOOST_AUTO_TEST_CASE(DefaultConstructedIsEmpty)
{
    // The state on every non-I/O rank: no records, every lookup misses (and must not crash).
    const LgrTransIndex<3> idx;
    BOOST_CHECK(idx.find(Key3{0, 0, 0}) == nullptr);
    BOOST_CHECK(idx.find(Key3{1, 2, 3}) == nullptr);
}

BOOST_AUTO_TEST_CASE(EmptyRecords)
{
    const LgrTransIndex<3> idx{std::vector<Rec3>{}};
    BOOST_CHECK(idx.find(Key3{0, 0, 0}) == nullptr);
}

BOOST_AUTO_TEST_CASE(SingleRecord)
{
    const std::vector<Rec3> recs{ {Key3{1, 4, 5}, 42.5} };
    const LgrTransIndex<3> idx{recs};

    const double* v = idx.find(Key3{1, 4, 5});
    BOOST_REQUIRE(v != nullptr);
    BOOST_CHECK_EQUAL(*v, 42.5);

    BOOST_CHECK(idx.find(Key3{1, 4, 6}) == nullptr);  // same level, absent
    BOOST_CHECK(idx.find(Key3{0, 4, 5}) == nullptr);  // different level, absent
}

BOOST_AUTO_TEST_CASE(ManyRecordsFoundWithExactValues)
{
    // Real-shaped keys: +X/+Y/+Z neighbour connections of a small refined grid, several levels.
    std::vector<Rec3> recs;
    const int NX = 6, NY = 5, NZ = 4;
    double v = 1.0;
    for (int level = 1; level <= 3; ++level) {
        for (int iz = 0; iz < NZ; ++iz) {
            for (int iy = 0; iy < NY; ++iy) {
                for (int ix = 0; ix < NX; ++ix) {
                    const int c = ix + iy * NX + iz * NX * NY;
                    if (ix + 1 < NX) recs.push_back({Key3{level, c, c + 1},          v += 1.0});
                    if (iy + 1 < NY) recs.push_back({Key3{level, c, c + NX},         v += 1.0});
                    if (iz + 1 < NZ) recs.push_back({Key3{level, c, c + NX * NY},    v += 1.0});
                }
            }
        }
    }

    const LgrTransIndex<3> idx{recs};

    for (const auto& r : recs) {
        const double* got = idx.find(r.first);
        BOOST_REQUIRE_MESSAGE(got != nullptr, "inserted key not found");
        BOOST_CHECK_EQUAL(*got, r.second);
    }
    // A clearly-absent key must terminate (the table always keeps at least one empty slot).
    BOOST_CHECK(idx.find(Key3{9, 999999, 1000000}) == nullptr);
}

BOOST_AUTO_TEST_CASE(CrossLevelN4)
{
    std::vector<Rec4> recs;
    double v = 100.0;
    for (int i = 0; i < 500; ++i) {
        recs.push_back({Key4{0, i, 1, i + 3}, v += 0.25});
    }
    const LgrTransIndex<4> idx{recs};

    for (const auto& r : recs) {
        const double* got = idx.find(r.first);
        BOOST_REQUIRE(got != nullptr);
        BOOST_CHECK_EQUAL(*got, r.second);
    }
    BOOST_CHECK(idx.find(Key4{0, 12345, 1, 67890}) == nullptr);
}

BOOST_AUTO_TEST_CASE(NegativeAndLargeKeyComponents)
{
    // hash mixes the ints through unsigned; negative and large components must round-trip.
    const std::vector<Rec3> recs{
        {Key3{-1, -2, -3},                  1.0},
        {Key3{0, 2000000000, 2000000001},   2.0},
        {Key3{-5, 7, 7},                    3.0},  // equal trailing components
    };
    const LgrTransIndex<3> idx{recs};

    BOOST_REQUIRE(idx.find(Key3{-1, -2, -3}) != nullptr);
    BOOST_CHECK_EQUAL(*idx.find(Key3{-1, -2, -3}),                1.0);
    BOOST_CHECK_EQUAL(*idx.find(Key3{0, 2000000000, 2000000001}), 2.0);
    BOOST_CHECK_EQUAL(*idx.find(Key3{-5, 7, 7}),                  3.0);
    BOOST_CHECK(idx.find(Key3{-1, -2, -4}) == nullptr);
}

BOOST_AUTO_TEST_CASE(LoadFactorBoundaryAllFound)
{
    // Sizes straddling the power-of-two capacity boundaries (load approaching the 0.7 cap):
    // every key must still be found and every miss must terminate.
    for (const int n : {1, 2, 3, 7, 8, 9, 100, 127, 128, 129, 1000}) {
        std::vector<Rec3> recs;
        recs.reserve(n);
        for (int i = 0; i < n; ++i) {
            recs.push_back({Key3{i % 4, i, i + 1}, static_cast<double>(i)});
        }
        const LgrTransIndex<3> idx{recs};

        for (const auto& r : recs) {
            BOOST_REQUIRE(idx.find(r.first) != nullptr);
        }
        BOOST_CHECK(idx.find(Key3{0, -1, -1}) == nullptr);  // absent must not loop forever
    }
}
