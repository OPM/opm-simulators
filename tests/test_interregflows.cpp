/*
  Copyright 2022 Equinor.

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

#define BOOST_TEST_MODULE TestEclInterRegFlows

#ifndef HAVE_MPI
// Suppress GCC diagnostics of the form
//
//   warning: "HAVE_MPI" is not defined, evaluates to 0
//
// when compiling with "-Wundef".
#define HAVE_MPI 0
#endif

#include <boost/test/unit_test.hpp>

#include <opm/output/data/InterRegFlowMap.hpp>

#include <opm/grid/common/p2pcommunicator.hh>

#include <opm/simulators/flow/InterRegFlows.hpp>

#include <cmath>
#include <cstddef>
#include <functional>
#include <string>
#include <type_traits>
#include <vector>

// +-----------> I
// |
// |   +-----+-----+
// |   |     |     |
// |   |  0  |  1  |
// |   |     |     |
// |   +-----+-----+
// |   |     |     |
// |   |  2  |  3  |
// |   |     |     |
// |   +-----+-----+
// |
// v J

namespace {
    Opm::data::InterRegFlowMap::FlowRates conn_01()
    {
        using Component = Opm::data::InterRegFlowMap::Component;

        auto rate = Opm::data::InterRegFlowMap::FlowRates{};

        rate[Component::Oil] = 1.0;
        rate[Component::Gas] = 2.0;
        rate[Component::Water] = 3.0;
        rate[Component::Disgas] = 4.0;
        rate[Component::Vapoil] = 5.0;

        return rate;
    }

    Opm::data::InterRegFlowMap::FlowRates conn_02()
    {
        using Component = Opm::data::InterRegFlowMap::Component;

        auto rate = Opm::data::InterRegFlowMap::FlowRates{};

        rate[Component::Oil] = 0.1;
        rate[Component::Gas] = 0.2;
        rate[Component::Water] = 0.3;
        rate[Component::Disgas] = 0.4;
        rate[Component::Vapoil] = 0.5;

        return rate;
    }

    Opm::data::InterRegFlowMap::FlowRates conn_13()
    {
        using Component = Opm::data::InterRegFlowMap::Component;

        auto rate = Opm::data::InterRegFlowMap::FlowRates{};

        rate[Component::Oil] = -0.2;
        rate[Component::Gas] = -0.4;
        rate[Component::Water] = -0.6;
        rate[Component::Disgas] = -0.8;
        rate[Component::Vapoil] = -1.0;

        return rate;
    }

    Opm::data::InterRegFlowMap::FlowRates conn_23()
    {
        using Component = Opm::data::InterRegFlowMap::Component;

        auto rate = Opm::data::InterRegFlowMap::FlowRates{};

        rate[Component::Oil] = -0.12;
        rate[Component::Gas] = -0.24;
        rate[Component::Water] = -0.36;
        rate[Component::Disgas] = -0.48;
        rate[Component::Vapoil] = -1.25;

        return rate;
    }

    std::vector<int> all_same_region()
    {
        return {
            1, 1,
            1, 1,
        };
    }

    std::vector<int> left_right_split_region()
    {
        return {
            1, 2,
            1, 2,
        };
    }

    std::vector<int> checker_board_region()
    {
        return {
            1, 2,
            2, 1,
        };
    }

    std::vector<int> all_separate_region()
    {
        return {
            4, 3,
            2, 1,
        };
    }

    namespace TwoProc {
        namespace P1 {
            bool isInterior(const int cellID)
            {
                return (cellID % 2) == 1;
            }
        } // namespace P1

        namespace P2 {
            bool isInterior(const int cellID)
            {
                return (cellID % 2) == 0;
            }
        } // namespace P2
    } // namespace TwoProc

    namespace FourProc {
        namespace P1 {
            bool isInterior(const int cellID)
            {
                return cellID == 1;
            }
        } // namespace P1

        namespace P2 {
            bool isInterior(const int cellID)
            {
                return cellID == 2;
            }
        } // namespace P2

        namespace P3 {
            bool isInterior(const int cellID)
            {
                return cellID == 3;
            }
        } // namespace P3

        namespace P4 {
            bool isInterior(const int cellID)
            {
                return cellID == 0;
            }
        } // namespace P4
    } // namespace FourProc
}

// =====================================================================

BOOST_AUTO_TEST_SUITE(Single_FIP_Region)

BOOST_AUTO_TEST_CASE(AllInternal)
{
    auto flows = Opm::InterRegFlowMapSingleFIP{ all_same_region() };
    BOOST_CHECK_EQUAL(flows.getLocalMaxRegionID(), std::size_t{1});

    BOOST_CHECK_MESSAGE(flows.assignGlobalMaxRegionID(5),
                        "Assigning large global maximum "
                        "Region ID must succeed");

    flows.addConnection({ 0, 0, true }, { 1, 1, true }, conn_01());
    flows.addConnection({ 0, 0, true }, { 2, 2, true }, conn_02());
    flows.addConnection({ 1, 1, true }, { 3, 3, true }, conn_13());
    flows.addConnection({ 2, 2, true }, { 3, 3, true }, conn_23());

    flows.compress();

    const auto& iregFlows = flows.getInterRegFlows();
    BOOST_CHECK_EQUAL(iregFlows.numRegions(), 5);

    for (auto r1 = 0*iregFlows.numRegions(); r1 < iregFlows.numRegions(); ++r1) {
        for (auto r2 = r1 + 1; r2 < iregFlows.numRegions(); ++r2) {
            auto flow = iregFlows.getInterRegFlows(r1, r2);
            BOOST_CHECK_MESSAGE(! flow.has_value(),
                                "There must not be inter-regional flow "
                                "between regions " << r1 << " and " << r2);
        }
    }
}

BOOST_AUTO_TEST_CASE(LeftInternal_RightOther)
{
    auto flows = Opm::InterRegFlowMapSingleFIP{ all_same_region() };
    BOOST_CHECK_EQUAL(flows.getLocalMaxRegionID(), std::size_t{1});

    BOOST_CHECK_MESSAGE(flows.assignGlobalMaxRegionID(5),
                        "Assigning large global maximum "
                        "Region ID must succeed");

    flows.addConnection({ 0, 0, true } , { 1, 1, false }, conn_01());
    flows.addConnection({ 0, 0, true } , { 2, 2, true } , conn_02());
    flows.addConnection({ 1, 1, false }, { 3, 3, false }, conn_13());
    flows.addConnection({ 2, 2, true } , { 3, 3, false }, conn_23());

    flows.compress();

    const auto& iregFlows = flows.getInterRegFlows();
    BOOST_CHECK_EQUAL(iregFlows.numRegions(), 5);

    for (auto r1 = 0*iregFlows.numRegions(); r1 < iregFlows.numRegions(); ++r1) {
        for (auto r2 = r1 + 1; r2 < iregFlows.numRegions(); ++r2) {
            auto flow = iregFlows.getInterRegFlows(r1, r2);
            BOOST_CHECK_MESSAGE(! flow.has_value(),
                                "There must not be inter-regional flow "
                                "between regions " << r1 << " and " << r2);
        }
    }
}

BOOST_AUTO_TEST_CASE(TopInternal_BottomOther)
{
    auto flows = Opm::InterRegFlowMapSingleFIP{ all_same_region() };
    BOOST_CHECK_EQUAL(flows.getLocalMaxRegionID(), std::size_t{1});

    BOOST_CHECK_MESSAGE(flows.assignGlobalMaxRegionID(5),
                        "Assigning large global maximum "
                        "Region ID must succeed");

    flows.addConnection({ 0, 0, true } , { 1, 1, true } , conn_01());
    flows.addConnection({ 0, 0, true } , { 2, 2, false }, conn_02());
    flows.addConnection({ 1, 1, true } , { 3, 3, false }, conn_13());
    flows.addConnection({ 2, 2, false }, { 3, 3, false }, conn_23());

    flows.compress();

    const auto& iregFlows = flows.getInterRegFlows();
    BOOST_CHECK_EQUAL(iregFlows.numRegions(), 5);

    for (auto r1 = 0*iregFlows.numRegions(); r1 < iregFlows.numRegions(); ++r1) {
        for (auto r2 = r1 + 1; r2 < iregFlows.numRegions(); ++r2) {
            auto flow = iregFlows.getInterRegFlows(r1, r2);
            BOOST_CHECK_MESSAGE(! flow.has_value(),
                                "There must not be inter-regional flow "
                                "between regions " << r1 << " and " << r2);
        }
    }
}

BOOST_AUTO_TEST_SUITE_END() // Single_FIP_Region

// =====================================================================

BOOST_AUTO_TEST_SUITE(Left_Right_Split_Region)

BOOST_AUTO_TEST_CASE(AllInternal)
{
    auto flows = Opm::InterRegFlowMapSingleFIP{ left_right_split_region() };
    BOOST_CHECK_EQUAL(flows.getLocalMaxRegionID(), std::size_t{2});

    flows.addConnection({ 0, 0, true }, { 1, 1, true }, conn_01());
    flows.addConnection({ 0, 0, true }, { 2, 2, true }, conn_02());
    flows.addConnection({ 1, 1, true }, { 3, 3, true }, conn_13());
    flows.addConnection({ 2, 2, true }, { 3, 3, true }, conn_23());

    flows.compress();

    const auto& iregFlows = flows.getInterRegFlows();
    BOOST_CHECK_EQUAL(iregFlows.numRegions(), 2);

    {
        const auto q12 = iregFlows.getInterRegFlows(0, 1);
        BOOST_REQUIRE_MESSAGE(q12.has_value(),
                              "There must be inter-region flows "
                              "between regions 1 and 2");

        const auto& [rate, sign] = q12.value();

        BOOST_CHECK_EQUAL(sign, 1.0);

        using FlowView = std::remove_cv_t<std::remove_reference_t<
            decltype(rate)>>;

        using Component = FlowView::Component;
        using Direction = FlowView::Direction;

        BOOST_CHECK_CLOSE(rate.flow(Component::Oil), 0.88f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Gas), 1.76f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Water), 2.64f, 1.0e-5);
        BOOST_CHECK_CLOSE(rate.flow(Component::Disgas), 3.52f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Vapoil), 3.75f, 1.0e-6);

        BOOST_CHECK_CLOSE(rate.flow(Component::Oil, Direction::Positive), 1.0f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Gas, Direction::Positive), 2.0f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Water, Direction::Positive), 3.0f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Disgas, Direction::Positive), 4.0f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Vapoil, Direction::Positive), 5.0f, 1.0e-6);

        BOOST_CHECK_CLOSE(rate.flow(Component::Oil, Direction::Negative), -0.12f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Gas, Direction::Negative), -0.24f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Water, Direction::Negative), -0.36f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Disgas, Direction::Negative), -0.48f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Vapoil, Direction::Negative), -1.25f, 1.0e-6);
    }

    {
        const auto q21 = iregFlows.getInterRegFlows(1, 0);
        BOOST_REQUIRE_MESSAGE(q21.has_value(),
                              "There must be inter-region flows "
                              "between regions 2 and 1");

        const auto& [rate, sign] = q21.value();

        BOOST_CHECK_EQUAL(sign, -1.0);

        using FlowView = std::remove_cv_t<std::remove_reference_t<
            decltype(rate)>>;

        using Component = FlowView::Component;
        using Direction = FlowView::Direction;

        BOOST_CHECK_CLOSE(rate.flow(Component::Oil), 0.88f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Gas), 1.76f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Water), 2.64f, 1.0e-5);
        BOOST_CHECK_CLOSE(rate.flow(Component::Disgas), 3.52f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Vapoil), 3.75f, 1.0e-6);

        BOOST_CHECK_CLOSE(rate.flow(Component::Oil, Direction::Positive), 1.0f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Gas, Direction::Positive), 2.0f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Water, Direction::Positive), 3.0f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Disgas, Direction::Positive), 4.0f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Vapoil, Direction::Positive), 5.0f, 1.0e-6);

        BOOST_CHECK_CLOSE(rate.flow(Component::Oil, Direction::Negative), -0.12f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Gas, Direction::Negative), -0.24f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Water, Direction::Negative), -0.36f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Disgas, Direction::Negative), -0.48f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Vapoil, Direction::Negative), -1.25f, 1.0e-6);
    }
}

BOOST_AUTO_TEST_CASE(LeftInternal_RightOther)
{
    auto flows = Opm::InterRegFlowMapSingleFIP{ left_right_split_region() };
    BOOST_CHECK_EQUAL(flows.getLocalMaxRegionID(), std::size_t{2});

    flows.addConnection({ 0, 0, true } , { 1, 1, false }, conn_01());
    flows.addConnection({ 0, 0, true } , { 2, 2, true } , conn_02());
    flows.addConnection({ 1, 1, false }, { 3, 3, false }, conn_13());
    flows.addConnection({ 2, 2, true } , { 3, 3, false }, conn_23());

    flows.compress();

    const auto& iregFlows = flows.getInterRegFlows();
    BOOST_CHECK_EQUAL(iregFlows.numRegions(), 2);

    {
        const auto q12 = iregFlows.getInterRegFlows(0, 1);
        BOOST_REQUIRE_MESSAGE(q12.has_value(),
                              "There must be inter-region flows "
                              "between regions 1 and 2");

        const auto& [rate, sign] = q12.value();

        BOOST_CHECK_EQUAL(sign, 1.0);

        using FlowView = std::remove_cv_t<std::remove_reference_t<
            decltype(rate)>>;

        using Component = FlowView::Component;
        using Direction = FlowView::Direction;

        BOOST_CHECK_CLOSE(rate.flow(Component::Oil), 0.88f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Gas), 1.76f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Water), 2.64f, 1.0e-5);
        BOOST_CHECK_CLOSE(rate.flow(Component::Disgas), 3.52f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Vapoil), 3.75f, 1.0e-6);

        BOOST_CHECK_CLOSE(rate.flow(Component::Oil, Direction::Positive), 1.0f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Gas, Direction::Positive), 2.0f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Water, Direction::Positive), 3.0f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Disgas, Direction::Positive), 4.0f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Vapoil, Direction::Positive), 5.0f, 1.0e-6);

        BOOST_CHECK_CLOSE(rate.flow(Component::Oil, Direction::Negative), -0.12f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Gas, Direction::Negative), -0.24f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Water, Direction::Negative), -0.36f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Disgas, Direction::Negative), -0.48f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Vapoil, Direction::Negative), -1.25f, 1.0e-6);
    }

    {
        const auto q21 = iregFlows.getInterRegFlows(1, 0);
        BOOST_REQUIRE_MESSAGE(q21.has_value(),
                              "There must be inter-region flows "
                              "between regions 2 and 1");

        const auto& [rate, sign] = q21.value();

        BOOST_CHECK_EQUAL(sign, -1.0);

        using FlowView = std::remove_cv_t<std::remove_reference_t<
            decltype(rate)>>;

        using Component = FlowView::Component;
        using Direction = FlowView::Direction;

        BOOST_CHECK_CLOSE(rate.flow(Component::Oil), 0.88f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Gas), 1.76f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Water), 2.64f, 1.0e-5);
        BOOST_CHECK_CLOSE(rate.flow(Component::Disgas), 3.52f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Vapoil), 3.75f, 1.0e-6);

        BOOST_CHECK_CLOSE(rate.flow(Component::Oil, Direction::Positive), 1.0f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Gas, Direction::Positive), 2.0f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Water, Direction::Positive), 3.0f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Disgas, Direction::Positive), 4.0f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Vapoil, Direction::Positive), 5.0f, 1.0e-6);

        BOOST_CHECK_CLOSE(rate.flow(Component::Oil, Direction::Negative), -0.12f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Gas, Direction::Negative), -0.24f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Water, Direction::Negative), -0.36f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Disgas, Direction::Negative), -0.48f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Vapoil, Direction::Negative), -1.25f, 1.0e-6);
    }
}

BOOST_AUTO_TEST_CASE(LeftOther_RightInternal)
{
    auto flows = Opm::InterRegFlowMapSingleFIP{ left_right_split_region() };
    BOOST_CHECK_EQUAL(flows.getLocalMaxRegionID(), std::size_t{2});

    flows.addConnection({ 0, 0, false }, { 1, 1, true } , conn_01());
    flows.addConnection({ 0, 0, false }, { 2, 2, false }, conn_02());
    flows.addConnection({ 1, 1, true } , { 3, 3, true } , conn_13());
    flows.addConnection({ 2, 2, false }, { 3, 3, true } , conn_23());

    flows.compress();

    const auto& iregFlows = flows.getInterRegFlows();
    BOOST_CHECK_EQUAL(iregFlows.numRegions(), 2);

    {
        const auto q12 = iregFlows.getInterRegFlows(0, 1);
        BOOST_REQUIRE_MESSAGE(! q12.has_value(),
                              "There must NOT be inter-region "
                              "flows between regions 1 and 2");
    }

    {
        const auto q21 = iregFlows.getInterRegFlows(1, 0);
        BOOST_REQUIRE_MESSAGE(! q21.has_value(),
                              "There must NOT be inter-region "
                              "flows between regions 1 and 2");
    }
}

BOOST_AUTO_TEST_CASE(TopOther_BottomInternal)
{
    auto flows = Opm::InterRegFlowMapSingleFIP{ left_right_split_region() };
    BOOST_CHECK_EQUAL(flows.getLocalMaxRegionID(), std::size_t{2});

    flows.addConnection({ 0, 0, true } , { 1, 1, true } , conn_01());
    flows.addConnection({ 0, 0, true } , { 2, 2, false }, conn_02());
    flows.addConnection({ 1, 1, true } , { 3, 3, false }, conn_13());
    flows.addConnection({ 2, 2, false }, { 3, 3, false }, conn_23());

    flows.compress();

    const auto& iregFlows = flows.getInterRegFlows();
    BOOST_CHECK_EQUAL(iregFlows.numRegions(), 2);

    {
        const auto q12 = iregFlows.getInterRegFlows(0, 1);
        BOOST_REQUIRE_MESSAGE(q12.has_value(),
                              "There must be inter-region flows "
                              "between regions 1 and 2");

        const auto& [rate, sign] = q12.value();

        BOOST_CHECK_EQUAL(sign, 1.0);

        using FlowView = std::remove_cv_t<std::remove_reference_t<
            decltype(rate)>>;

        using Component = FlowView::Component;
        using Direction = FlowView::Direction;

        BOOST_CHECK_CLOSE(rate.flow(Component::Oil), 1.0f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Gas), 2.0f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Water), 3.0f, 1.0e-5);
        BOOST_CHECK_CLOSE(rate.flow(Component::Disgas), 4.0f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Vapoil), 5.0f, 1.0e-6);

        BOOST_CHECK_CLOSE(rate.flow(Component::Oil, Direction::Positive), 1.0f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Gas, Direction::Positive), 2.0f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Water, Direction::Positive), 3.0f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Disgas, Direction::Positive), 4.0f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Vapoil, Direction::Positive), 5.0f, 1.0e-6);

        BOOST_CHECK_CLOSE(rate.flow(Component::Oil, Direction::Negative), 0.0f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Gas, Direction::Negative), 0.0f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Water, Direction::Negative), 0.0f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Disgas, Direction::Negative), 0.0f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Vapoil, Direction::Negative), 0.0f, 1.0e-6);
    }
}

BOOST_AUTO_TEST_CASE(CheckerBoard_Internal)
{
    auto flows = Opm::InterRegFlowMapSingleFIP{ left_right_split_region() };
    BOOST_CHECK_EQUAL(flows.getLocalMaxRegionID(), std::size_t{2});

    flows.addConnection({ 0, 0, false }, { 1, 1, true } , conn_01());
    flows.addConnection({ 0, 0, false }, { 2, 2, true } , conn_02());
    flows.addConnection({ 1, 1, true } , { 3, 3, false }, conn_13());
    flows.addConnection({ 2, 2, true } , { 3, 3, false }, conn_23());

    flows.compress();

    const auto& iregFlows = flows.getInterRegFlows();
    BOOST_CHECK_EQUAL(iregFlows.numRegions(), 2);

    {
        const auto q12 = iregFlows.getInterRegFlows(0, 1);
        BOOST_REQUIRE_MESSAGE(q12.has_value(),
                              "There must be inter-region flows "
                              "between regions 1 and 2");

        const auto& [rate, sign] = q12.value();

        BOOST_CHECK_EQUAL(sign, 1.0);

        using FlowView = std::remove_cv_t<std::remove_reference_t<
            decltype(rate)>>;

        using Component = FlowView::Component;
        using Direction = FlowView::Direction;

        BOOST_CHECK_CLOSE(rate.flow(Component::Oil), -0.12f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Gas), -0.24f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Water), -0.36f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Disgas), -0.48f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Vapoil), -1.25f, 1.0e-6);

        BOOST_CHECK_CLOSE(rate.flow(Component::Oil, Direction::Positive), 0.0f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Gas, Direction::Positive), 0.0f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Water, Direction::Positive), 0.0f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Disgas, Direction::Positive), 0.0f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Vapoil, Direction::Positive), 0.0f, 1.0e-6);

        BOOST_CHECK_CLOSE(rate.flow(Component::Oil, Direction::Negative), -0.12f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Gas, Direction::Negative), -0.24f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Water, Direction::Negative), -0.36f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Disgas, Direction::Negative), -0.48f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Vapoil, Direction::Negative), -1.25f, 1.0e-6);
    }
}

BOOST_AUTO_TEST_SUITE_END() // Left_Right_Split_Region

// =====================================================================

BOOST_AUTO_TEST_SUITE(Checker_Board_Region)

BOOST_AUTO_TEST_CASE(AllInternal)
{
    auto flows = Opm::InterRegFlowMapSingleFIP{ checker_board_region() };
    BOOST_CHECK_EQUAL(flows.getLocalMaxRegionID(), std::size_t{2});

    flows.addConnection({ 0, 0, true }, { 1, 1, true }, conn_01());
    flows.addConnection({ 0, 0, true }, { 2, 2, true }, conn_02());
    flows.addConnection({ 1, 1, true }, { 3, 3, true }, conn_13());
    flows.addConnection({ 2, 2, true }, { 3, 3, true }, conn_23());

    flows.compress();

    const auto& iregFlows = flows.getInterRegFlows();
    BOOST_CHECK_EQUAL(iregFlows.numRegions(), 2);

    {
        const auto q12 = iregFlows.getInterRegFlows(0, 1);
        BOOST_REQUIRE_MESSAGE(q12.has_value(),
                              "There must be inter-region flows "
                              "between regions 1 and 2");

        const auto& [rate, sign] = q12.value();

        BOOST_CHECK_EQUAL(sign, 1.0);

        using FlowView = std::remove_cv_t<std::remove_reference_t<
            decltype(rate)>>;

        using Component = FlowView::Component;
        using Direction = FlowView::Direction;

        BOOST_CHECK_CLOSE(rate.flow(Component::Oil), 1.42f, 1.0e-5);
        BOOST_CHECK_CLOSE(rate.flow(Component::Gas), 2.84f, 1.0e-5);
        BOOST_CHECK_CLOSE(rate.flow(Component::Water), 4.26f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Disgas), 5.68f, 1.0e-5);
        BOOST_CHECK_CLOSE(rate.flow(Component::Vapoil), 7.75f, 1.0e-6);

        BOOST_CHECK_CLOSE(rate.flow(Component::Oil, Direction::Positive), 1.42f, 1.0e-5);
        BOOST_CHECK_CLOSE(rate.flow(Component::Gas, Direction::Positive), 2.84f, 1.0e-5);
        BOOST_CHECK_CLOSE(rate.flow(Component::Water, Direction::Positive), 4.26f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Disgas, Direction::Positive), 5.68f, 1.0e-5);
        BOOST_CHECK_CLOSE(rate.flow(Component::Vapoil, Direction::Positive), 7.75f, 1.0e-6);

        BOOST_CHECK_CLOSE(rate.flow(Component::Oil, Direction::Negative), 0.0f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Gas, Direction::Negative), 0.0f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Water, Direction::Negative), 0.0f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Disgas, Direction::Negative), 0.0f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Vapoil, Direction::Negative), 0.0f, 1.0e-6);
    }
}

BOOST_AUTO_TEST_CASE(LeftInternal_RightOther)
{
    auto flows = Opm::InterRegFlowMapSingleFIP{ checker_board_region() };
    BOOST_CHECK_EQUAL(flows.getLocalMaxRegionID(), std::size_t{2});

    flows.addConnection({ 0, 0, true } , { 1, 1, false }, conn_01());
    flows.addConnection({ 0, 0, true } , { 2, 2, true } , conn_02());
    flows.addConnection({ 1, 1, false }, { 3, 3, false }, conn_13());
    flows.addConnection({ 2, 2, true } , { 3, 3, false }, conn_23());

    flows.compress();

    const auto& iregFlows = flows.getInterRegFlows();
    BOOST_CHECK_EQUAL(iregFlows.numRegions(), 2);

    {
        const auto q12 = iregFlows.getInterRegFlows(0, 1);
        BOOST_REQUIRE_MESSAGE(q12.has_value(),
                              "There must be inter-region flows "
                              "between regions 1 and 2");

        const auto& [rate, sign] = q12.value();

        BOOST_CHECK_EQUAL(sign, 1.0);

        using FlowView = std::remove_cv_t<std::remove_reference_t<
            decltype(rate)>>;

        using Component = FlowView::Component;
        using Direction = FlowView::Direction;

        BOOST_CHECK_CLOSE(rate.flow(Component::Oil), 1.22f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Gas), 2.44f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Water), 3.66f, 1.0e-5);
        BOOST_CHECK_CLOSE(rate.flow(Component::Disgas), 4.88f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Vapoil), 6.75f, 1.0e-6);

        BOOST_CHECK_CLOSE(rate.flow(Component::Oil, Direction::Positive), 1.22f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Gas, Direction::Positive), 2.44f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Water, Direction::Positive), 3.66f, 1.0e-5);
        BOOST_CHECK_CLOSE(rate.flow(Component::Disgas, Direction::Positive), 4.88f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Vapoil, Direction::Positive), 6.75f, 1.0e-6);

        BOOST_CHECK_CLOSE(rate.flow(Component::Oil, Direction::Negative), 0.0f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Gas, Direction::Negative), 0.0f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Water, Direction::Negative), 0.0f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Disgas, Direction::Negative), 0.0f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Vapoil, Direction::Negative), 0.0f, 1.0e-6);
    }
}

BOOST_AUTO_TEST_CASE(LeftOther_RightInternal)
{
    auto flows = Opm::InterRegFlowMapSingleFIP{ checker_board_region() };
    BOOST_CHECK_EQUAL(flows.getLocalMaxRegionID(), std::size_t{2});

    flows.addConnection({ 0, 0, false }, { 1, 1, true } , conn_01());
    flows.addConnection({ 0, 0, false }, { 2, 2, false }, conn_02());
    flows.addConnection({ 1, 1, true } , { 3, 3, true } , conn_13());
    flows.addConnection({ 2, 2, false }, { 3, 3, true } , conn_23());

    flows.compress();

    const auto& iregFlows = flows.getInterRegFlows();
    BOOST_CHECK_EQUAL(iregFlows.numRegions(), 2);

    {
        const auto q12 = iregFlows.getInterRegFlows(0, 1);
        BOOST_REQUIRE_MESSAGE(q12.has_value(),
                              "There must be inter-region flows "
                              "between regions 1 and 2");

        const auto& [rate, sign] = q12.value();

        BOOST_CHECK_EQUAL(sign, 1.0);

        using FlowView = std::remove_cv_t<std::remove_reference_t<
            decltype(rate)>>;

        using Component = FlowView::Component;
        using Direction = FlowView::Direction;

        BOOST_CHECK_CLOSE(rate.flow(Component::Oil), 0.2f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Gas), 0.4f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Water), 0.6f, 7.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Disgas), 0.8f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Vapoil), 1.0f, 1.0e-6);

        BOOST_CHECK_CLOSE(rate.flow(Component::Oil, Direction::Positive), 0.2f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Gas, Direction::Positive), 0.4f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Water, Direction::Positive), 0.6f, 7.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Disgas, Direction::Positive), 0.8f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Vapoil, Direction::Positive), 1.0f, 1.0e-6);

        BOOST_CHECK_CLOSE(rate.flow(Component::Oil, Direction::Negative), 0.0f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Gas, Direction::Negative), 0.0f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Water, Direction::Negative), 0.0f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Disgas, Direction::Negative), 0.0f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Vapoil, Direction::Negative), 0.0f, 1.0e-6);
    }
}

BOOST_AUTO_TEST_CASE(TopOther_BottomInternal)
{
    auto flows = Opm::InterRegFlowMapSingleFIP{ checker_board_region() };
    BOOST_CHECK_EQUAL(flows.getLocalMaxRegionID(), std::size_t{2});

    flows.addConnection({ 0, 0, false }, { 1, 1, false }, conn_01());
    flows.addConnection({ 0, 0, false }, { 2, 2, true } , conn_02());
    flows.addConnection({ 1, 1, false }, { 3, 3, true } , conn_13());
    flows.addConnection({ 2, 2, true } , { 3, 3, true } , conn_23());

    flows.compress();

    const auto& iregFlows = flows.getInterRegFlows();
    BOOST_CHECK_EQUAL(iregFlows.numRegions(), 2);

    {
        const auto q12 = iregFlows.getInterRegFlows(0, 1);
        BOOST_REQUIRE_MESSAGE(q12.has_value(),
                              "There must be inter-region flows "
                              "between regions 1 and 2");

        const auto& [rate, sign] = q12.value();

        BOOST_CHECK_EQUAL(sign, 1.0);

        using FlowView = std::remove_cv_t<std::remove_reference_t<
            decltype(rate)>>;

        using Component = FlowView::Component;
        using Direction = FlowView::Direction;

        BOOST_CHECK_CLOSE(rate.flow(Component::Oil), 0.12f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Gas), 0.24f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Water), 0.36f, 7.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Disgas), 0.48f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Vapoil), 1.25f, 1.0e-6);

        BOOST_CHECK_CLOSE(rate.flow(Component::Oil, Direction::Positive), 0.12f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Gas, Direction::Positive), 0.24f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Water, Direction::Positive), 0.36f, 7.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Disgas, Direction::Positive), 0.48f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Vapoil, Direction::Positive), 1.25f, 1.0e-6);

        BOOST_CHECK_CLOSE(rate.flow(Component::Oil, Direction::Negative), 0.0f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Gas, Direction::Negative), 0.0f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Water, Direction::Negative), 0.0f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Disgas, Direction::Negative), 0.0f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Vapoil, Direction::Negative), 0.0f, 1.0e-6);
    }
}

BOOST_AUTO_TEST_CASE(CheckerBoard_Internal)
{
    auto flows = Opm::InterRegFlowMapSingleFIP{ checker_board_region() };
    BOOST_CHECK_EQUAL(flows.getLocalMaxRegionID(), std::size_t{2});

    flows.addConnection({ 0, 0, false }, { 1, 1, true } , conn_01());
    flows.addConnection({ 0, 0, false }, { 2, 2, true } , conn_02());
    flows.addConnection({ 1, 1, true } , { 3, 3, false }, conn_13());
    flows.addConnection({ 2, 2, true } , { 3, 3, false }, conn_23());

    flows.compress();

    const auto& iregFlows = flows.getInterRegFlows();
    BOOST_CHECK_EQUAL(iregFlows.numRegions(), 2);

    {
        const auto q12 = iregFlows.getInterRegFlows(0, 1);
        BOOST_REQUIRE_MESSAGE(q12.has_value(),
                              "There must be inter-region flows "
                              "between regions 1 and 2");

        const auto& [rate, sign] = q12.value();

        BOOST_CHECK_EQUAL(sign, 1.0);

        using FlowView = std::remove_cv_t<std::remove_reference_t<
            decltype(rate)>>;

        using Component = FlowView::Component;
        using Direction = FlowView::Direction;

        BOOST_CHECK_CLOSE(rate.flow(Component::Oil), 0.32f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Gas), 0.64f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Water), 0.96f, 7.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Disgas), 1.28f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Vapoil), 2.25f, 1.0e-6);

        BOOST_CHECK_CLOSE(rate.flow(Component::Oil, Direction::Positive), 0.32f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Gas, Direction::Positive), 0.64f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Water, Direction::Positive), 0.96f, 7.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Disgas, Direction::Positive), 1.28f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Vapoil, Direction::Positive), 2.25f, 1.0e-6);

        BOOST_CHECK_CLOSE(rate.flow(Component::Oil, Direction::Negative), 0.0f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Gas, Direction::Negative), 0.0f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Water, Direction::Negative), 0.0f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Disgas, Direction::Negative), 0.0f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Vapoil, Direction::Negative), 0.0f, 1.0e-6);
    }
}

BOOST_AUTO_TEST_CASE(Reverse_CheckerBoard_Internal)
{
    auto flows = Opm::InterRegFlowMapSingleFIP{ checker_board_region() };
    BOOST_CHECK_EQUAL(flows.getLocalMaxRegionID(), std::size_t{2});

    flows.addConnection({ 0, 0, true } , { 1, 1, false }, conn_01());
    flows.addConnection({ 0, 0, true } , { 2, 2, false }, conn_02());
    flows.addConnection({ 1, 1, false }, { 3, 3, true } , conn_13());
    flows.addConnection({ 2, 2, false }, { 3, 3, true } , conn_23());

    flows.compress();

    const auto& iregFlows = flows.getInterRegFlows();
    BOOST_CHECK_EQUAL(iregFlows.numRegions(), 2);

    {
        const auto q12 = iregFlows.getInterRegFlows(0, 1);
        BOOST_REQUIRE_MESSAGE(q12.has_value(),
                              "There must be inter-region flows "
                              "between regions 1 and 2");

        const auto& [rate, sign] = q12.value();

        BOOST_CHECK_EQUAL(sign, 1.0);

        using FlowView = std::remove_cv_t<std::remove_reference_t<
            decltype(rate)>>;

        using Component = FlowView::Component;
        using Direction = FlowView::Direction;

        BOOST_CHECK_CLOSE(rate.flow(Component::Oil), 1.1f, 1.0e-5);
        BOOST_CHECK_CLOSE(rate.flow(Component::Gas), 2.2f, 1.0e-5);
        BOOST_CHECK_CLOSE(rate.flow(Component::Water), 3.3f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Disgas), 4.4f, 1.0e-5);
        BOOST_CHECK_CLOSE(rate.flow(Component::Vapoil), 5.5f, 1.0e-6);

        BOOST_CHECK_CLOSE(rate.flow(Component::Oil, Direction::Positive), 1.1f, 1.0e-5);
        BOOST_CHECK_CLOSE(rate.flow(Component::Gas, Direction::Positive), 2.2f, 1.0e-5);
        BOOST_CHECK_CLOSE(rate.flow(Component::Water, Direction::Positive), 3.3f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Disgas, Direction::Positive), 4.4f, 1.0e-5);
        BOOST_CHECK_CLOSE(rate.flow(Component::Vapoil, Direction::Positive), 5.5f, 1.0e-6);

        BOOST_CHECK_CLOSE(rate.flow(Component::Oil, Direction::Negative), 0.0f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Gas, Direction::Negative), 0.0f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Water, Direction::Negative), 0.0f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Disgas, Direction::Negative), 0.0f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Vapoil, Direction::Negative), 0.0f, 1.0e-6);
    }
}

BOOST_AUTO_TEST_SUITE_END() // Checker_Board_Region

// =====================================================================

BOOST_AUTO_TEST_SUITE(All_Separate_Region)

BOOST_AUTO_TEST_CASE(AllInternal)
{
    auto flows = Opm::InterRegFlowMapSingleFIP{ all_separate_region() };
    BOOST_CHECK_EQUAL(flows.getLocalMaxRegionID(), std::size_t{4});

    BOOST_CHECK_MESSAGE(! flows.assignGlobalMaxRegionID(2),
                        "Assigning small global maximum "
                        "Region ID must NOT succeed");

    flows.addConnection({ 0, 0, true }, { 1, 1, true }, conn_01());
    flows.addConnection({ 0, 0, true }, { 2, 2, true }, conn_02());
    flows.addConnection({ 1, 1, true }, { 3, 3, true }, conn_13());
    flows.addConnection({ 2, 2, true }, { 3, 3, true }, conn_23());

    flows.compress();

    const auto& iregFlows = flows.getInterRegFlows();
    BOOST_CHECK_EQUAL(iregFlows.numRegions(), 4);

    {
        const auto q12 = iregFlows.getInterRegFlows(0, 1);
        BOOST_REQUIRE_MESSAGE(q12.has_value(),
                              "There must be inter-region flows "
                              "between regions 1 and 2");

        const auto& [rate, sign] = q12.value();

        BOOST_CHECK_EQUAL(sign, 1.0);

        using FlowView = std::remove_cv_t<std::remove_reference_t<
            decltype(rate)>>;

        using Component = FlowView::Component;
        using Direction = FlowView::Direction;

        BOOST_CHECK_CLOSE(rate.flow(Component::Oil), 0.12f, 1.0e-5);
        BOOST_CHECK_CLOSE(rate.flow(Component::Gas), 0.24f, 1.0e-5);
        BOOST_CHECK_CLOSE(rate.flow(Component::Water), 0.36f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Disgas), 0.48f, 1.0e-5);
        BOOST_CHECK_CLOSE(rate.flow(Component::Vapoil), 1.25f, 1.0e-6);

        BOOST_CHECK_CLOSE(rate.flow(Component::Oil, Direction::Positive), 0.12f, 1.0e-5);
        BOOST_CHECK_CLOSE(rate.flow(Component::Gas, Direction::Positive), 0.24f, 1.0e-5);
        BOOST_CHECK_CLOSE(rate.flow(Component::Water, Direction::Positive), 0.36f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Disgas, Direction::Positive), 0.48f, 1.0e-5);
        BOOST_CHECK_CLOSE(rate.flow(Component::Vapoil, Direction::Positive), 1.25f, 1.0e-6);

        BOOST_CHECK_CLOSE(rate.flow(Component::Oil, Direction::Negative), 0.0f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Gas, Direction::Negative), 0.0f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Water, Direction::Negative), 0.0f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Disgas, Direction::Negative), 0.0f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Vapoil, Direction::Negative), 0.0f, 1.0e-6);
    }

    {
        const auto q13 = iregFlows.getInterRegFlows(0, 2);
        BOOST_REQUIRE_MESSAGE(q13.has_value(),
                              "There must be inter-region flows "
                              "between regions 1 and 3");

        const auto& [rate, sign] = q13.value();

        BOOST_CHECK_EQUAL(sign, 1.0);

        using FlowView = std::remove_cv_t<std::remove_reference_t<
            decltype(rate)>>;

        using Component = FlowView::Component;
        using Direction = FlowView::Direction;

        BOOST_CHECK_CLOSE(rate.flow(Component::Oil), 0.2f, 1.0e-5);
        BOOST_CHECK_CLOSE(rate.flow(Component::Gas), 0.4f, 1.0e-5);
        BOOST_CHECK_CLOSE(rate.flow(Component::Water), 0.6f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Disgas), 0.8f, 1.0e-5);
        BOOST_CHECK_CLOSE(rate.flow(Component::Vapoil), 1.0f, 1.0e-6);

        BOOST_CHECK_CLOSE(rate.flow(Component::Oil, Direction::Positive), 0.2f, 1.0e-5);
        BOOST_CHECK_CLOSE(rate.flow(Component::Gas, Direction::Positive), 0.4f, 1.0e-5);
        BOOST_CHECK_CLOSE(rate.flow(Component::Water, Direction::Positive), 0.6f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Disgas, Direction::Positive), 0.8f, 1.0e-5);
        BOOST_CHECK_CLOSE(rate.flow(Component::Vapoil, Direction::Positive), 1.0f, 1.0e-6);

        BOOST_CHECK_CLOSE(rate.flow(Component::Oil, Direction::Negative), 0.0f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Gas, Direction::Negative), 0.0f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Water, Direction::Negative), 0.0f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Disgas, Direction::Negative), 0.0f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Vapoil, Direction::Negative), 0.0f, 1.0e-6);
    }

    {
        const auto q14 = iregFlows.getInterRegFlows(0, 3);
        BOOST_CHECK_MESSAGE(! q14.has_value(),
                            "There must NOT be inter-region "
                            "flows between regions 1 and 4");
    }

    {
        const auto q23 = iregFlows.getInterRegFlows(1, 2);
        BOOST_CHECK_MESSAGE(! q23.has_value(),
                            "There must NOT be inter-region "
                            "flows between regions 2 and 3");
    }

    {
        const auto q24 = iregFlows.getInterRegFlows(1, 3);
        BOOST_REQUIRE_MESSAGE(q24.has_value(),
                              "There must be inter-region flows "
                              "between regions 2 and 4");

        const auto& [rate, sign] = q24.value();

        BOOST_CHECK_EQUAL(sign, 1.0);

        using FlowView = std::remove_cv_t<std::remove_reference_t<
            decltype(rate)>>;

        using Component = FlowView::Component;
        using Direction = FlowView::Direction;

        BOOST_CHECK_CLOSE(rate.flow(Component::Oil), -0.1f, 1.0e-5);
        BOOST_CHECK_CLOSE(rate.flow(Component::Gas), -0.2f, 1.0e-5);
        BOOST_CHECK_CLOSE(rate.flow(Component::Water), -0.3f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Disgas), -0.4f, 1.0e-5);
        BOOST_CHECK_CLOSE(rate.flow(Component::Vapoil), -0.5f, 1.0e-6);

        BOOST_CHECK_CLOSE(rate.flow(Component::Oil, Direction::Positive), 0.0f, 1.0e-5);
        BOOST_CHECK_CLOSE(rate.flow(Component::Gas, Direction::Positive), 0.0f, 1.0e-5);
        BOOST_CHECK_CLOSE(rate.flow(Component::Water, Direction::Positive), 0.0f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Disgas, Direction::Positive), 0.0f, 1.0e-5);
        BOOST_CHECK_CLOSE(rate.flow(Component::Vapoil, Direction::Positive), 0.0f, 1.0e-6);

        BOOST_CHECK_CLOSE(rate.flow(Component::Oil, Direction::Negative), -0.1f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Gas, Direction::Negative), -0.2f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Water, Direction::Negative), -0.3f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Disgas, Direction::Negative), -0.4f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Vapoil, Direction::Negative), -0.5f, 1.0e-6);
    }
}

BOOST_AUTO_TEST_CASE(LeftInternal_RightOther)
{
    auto flows = Opm::InterRegFlowMapSingleFIP{ all_separate_region() };
    BOOST_CHECK_EQUAL(flows.getLocalMaxRegionID(), std::size_t{4});

    flows.addConnection({ 0, 0, true } , { 1, 1, false }, conn_01());
    flows.addConnection({ 0, 0, true } , { 2, 2, true } , conn_02());
    flows.addConnection({ 1, 1, false }, { 3, 3, false }, conn_13());
    flows.addConnection({ 2, 2, true } , { 3, 3, false }, conn_23());

    flows.compress();

    const auto& iregFlows = flows.getInterRegFlows();
    BOOST_CHECK_EQUAL(iregFlows.numRegions(), 4);

    {
        const auto q12 = iregFlows.getInterRegFlows(0, 1);
        BOOST_REQUIRE_MESSAGE(q12.has_value(),
                              "There must be inter-region flows "
                              "between regions 1 and 2");

        const auto& [rate, sign] = q12.value();

        BOOST_CHECK_EQUAL(sign, 1.0);

        using FlowView = std::remove_cv_t<std::remove_reference_t<
            decltype(rate)>>;

        using Component = FlowView::Component;
        using Direction = FlowView::Direction;

        BOOST_CHECK_CLOSE(rate.flow(Component::Oil), 0.12f, 1.0e-5);
        BOOST_CHECK_CLOSE(rate.flow(Component::Gas), 0.24f, 1.0e-5);
        BOOST_CHECK_CLOSE(rate.flow(Component::Water), 0.36f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Disgas), 0.48f, 1.0e-5);
        BOOST_CHECK_CLOSE(rate.flow(Component::Vapoil), 1.25f, 1.0e-6);

        BOOST_CHECK_CLOSE(rate.flow(Component::Oil, Direction::Positive), 0.12f, 1.0e-5);
        BOOST_CHECK_CLOSE(rate.flow(Component::Gas, Direction::Positive), 0.24f, 1.0e-5);
        BOOST_CHECK_CLOSE(rate.flow(Component::Water, Direction::Positive), 0.36f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Disgas, Direction::Positive), 0.48f, 1.0e-5);
        BOOST_CHECK_CLOSE(rate.flow(Component::Vapoil, Direction::Positive), 1.25f, 1.0e-6);

        BOOST_CHECK_CLOSE(rate.flow(Component::Oil, Direction::Negative), 0.0f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Gas, Direction::Negative), 0.0f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Water, Direction::Negative), 0.0f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Disgas, Direction::Negative), 0.0f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Vapoil, Direction::Negative), 0.0f, 1.0e-6);
    }

    {
        const auto q13 = iregFlows.getInterRegFlows(0, 2);
        BOOST_CHECK_MESSAGE(! q13.has_value(),
                            "There must NOT be inter-region "
                            "flows between regions 1 and 3");
    }

    {
        const auto q14 = iregFlows.getInterRegFlows(0, 3);
        BOOST_CHECK_MESSAGE(! q14.has_value(),
                            "There must NOT be inter-region "
                            "flows between regions 1 and 4");
    }

    {
        const auto q23 = iregFlows.getInterRegFlows(1, 2);
        BOOST_CHECK_MESSAGE(! q23.has_value(),
                            "There must NOT be inter-region "
                            "flows between regions 2 and 3");
    }

    {
        const auto q24 = iregFlows.getInterRegFlows(1, 3);
        BOOST_REQUIRE_MESSAGE(q24.has_value(),
                              "There must be inter-region flows "
                              "between regions 2 and 4");

        const auto& [rate, sign] = q24.value();

        BOOST_CHECK_EQUAL(sign, 1.0);

        using FlowView = std::remove_cv_t<std::remove_reference_t<
            decltype(rate)>>;

        using Component = FlowView::Component;
        using Direction = FlowView::Direction;

        BOOST_CHECK_CLOSE(rate.flow(Component::Oil), -0.1f, 1.0e-5);
        BOOST_CHECK_CLOSE(rate.flow(Component::Gas), -0.2f, 1.0e-5);
        BOOST_CHECK_CLOSE(rate.flow(Component::Water), -0.3f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Disgas), -0.4f, 1.0e-5);
        BOOST_CHECK_CLOSE(rate.flow(Component::Vapoil), -0.5f, 1.0e-6);

        BOOST_CHECK_CLOSE(rate.flow(Component::Oil, Direction::Positive), 0.0f, 1.0e-5);
        BOOST_CHECK_CLOSE(rate.flow(Component::Gas, Direction::Positive), 0.0f, 1.0e-5);
        BOOST_CHECK_CLOSE(rate.flow(Component::Water, Direction::Positive), 0.0f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Disgas, Direction::Positive), 0.0f, 1.0e-5);
        BOOST_CHECK_CLOSE(rate.flow(Component::Vapoil, Direction::Positive), 0.0f, 1.0e-6);

        BOOST_CHECK_CLOSE(rate.flow(Component::Oil, Direction::Negative), -0.1f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Gas, Direction::Negative), -0.2f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Water, Direction::Negative), -0.3f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Disgas, Direction::Negative), -0.4f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Vapoil, Direction::Negative), -0.5f, 1.0e-6);
    }

    {
        const auto q34 = iregFlows.getInterRegFlows(2, 3);
        BOOST_REQUIRE_MESSAGE(q34.has_value(),
                              "There must be inter-region flows "
                              "between regions 3 and 4");

        const auto& [rate, sign] = q34.value();

        BOOST_CHECK_EQUAL(sign, 1.0);

        using FlowView = std::remove_cv_t<std::remove_reference_t<
            decltype(rate)>>;

        using Component = FlowView::Component;
        using Direction = FlowView::Direction;

        BOOST_CHECK_CLOSE(rate.flow(Component::Oil), -1.0f, 1.0e-5);
        BOOST_CHECK_CLOSE(rate.flow(Component::Gas), -2.0f, 1.0e-5);
        BOOST_CHECK_CLOSE(rate.flow(Component::Water), -3.0f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Disgas), -4.0f, 1.0e-5);
        BOOST_CHECK_CLOSE(rate.flow(Component::Vapoil), -5.0f, 1.0e-6);

        BOOST_CHECK_CLOSE(rate.flow(Component::Oil, Direction::Positive), 0.0f, 1.0e-5);
        BOOST_CHECK_CLOSE(rate.flow(Component::Gas, Direction::Positive), 0.0f, 1.0e-5);
        BOOST_CHECK_CLOSE(rate.flow(Component::Water, Direction::Positive), 0.0f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Disgas, Direction::Positive), 0.0f, 1.0e-5);
        BOOST_CHECK_CLOSE(rate.flow(Component::Vapoil, Direction::Positive), 0.0f, 1.0e-6);

        BOOST_CHECK_CLOSE(rate.flow(Component::Oil, Direction::Negative), -1.0f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Gas, Direction::Negative), -2.0f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Water, Direction::Negative), -3.0f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Disgas, Direction::Negative), -4.0f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Vapoil, Direction::Negative), -5.0f, 1.0e-6);
    }
}

BOOST_AUTO_TEST_CASE(LeftOther_RightInternal)
{
    auto flows = Opm::InterRegFlowMapSingleFIP{ all_separate_region() };
    BOOST_CHECK_EQUAL(flows.getLocalMaxRegionID(), std::size_t{4});

    flows.addConnection({ 0, 0, false }, { 1, 1, true } , conn_01());
    flows.addConnection({ 0, 0, false }, { 2, 2, false }, conn_02());
    flows.addConnection({ 1, 1, true } , { 3, 3, true } , conn_13());
    flows.addConnection({ 2, 2, false }, { 3, 3, true } , conn_23());

    flows.compress();

    const auto& iregFlows = flows.getInterRegFlows();
    BOOST_CHECK_EQUAL(iregFlows.numRegions(), 4);

    {
        const auto q12 = iregFlows.getInterRegFlows(0, 1);
        BOOST_REQUIRE_MESSAGE(! q12.has_value(),
                              "There must NOT be inter-region flows "
                              "between regions 1 and 2");
    }

    {
        const auto q13 = iregFlows.getInterRegFlows(0, 2);
        BOOST_REQUIRE_MESSAGE(q13.has_value(),
                              "There must be inter-region flows "
                              "between regions 1 and 3");

        const auto& [rate, sign] = q13.value();

        BOOST_CHECK_EQUAL(sign, 1.0);

        using FlowView = std::remove_cv_t<std::remove_reference_t<
            decltype(rate)>>;

        using Component = FlowView::Component;
        using Direction = FlowView::Direction;

        BOOST_CHECK_CLOSE(rate.flow(Component::Oil), 0.2f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Gas), 0.4f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Water), 0.6f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Disgas), 0.8f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Vapoil), 1.0f, 1.0e-6);

        BOOST_CHECK_CLOSE(rate.flow(Component::Oil, Direction::Positive), 0.2f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Gas, Direction::Positive), 0.4f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Water, Direction::Positive), 0.6f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Disgas, Direction::Positive), 0.8f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Vapoil, Direction::Positive), 1.0f, 1.0e-6);

        BOOST_CHECK_CLOSE(rate.flow(Component::Oil, Direction::Negative), 0.0f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Gas, Direction::Negative), 0.0f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Water, Direction::Negative), 0.0f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Disgas, Direction::Negative), 0.0f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Vapoil, Direction::Negative), 0.0f, 1.0e-6);
    }

    {
        const auto q14 = iregFlows.getInterRegFlows(0, 3);
        BOOST_CHECK_MESSAGE(! q14.has_value(),
                            "There must NOT be inter-region "
                            "flows between regions 1 and 4");
    }

    {
        const auto q23 = iregFlows.getInterRegFlows(1, 2);
        BOOST_CHECK_MESSAGE(! q23.has_value(),
                            "There must NOT be inter-region "
                            "flows between regions 2 and 3");
    }

    {
        const auto q24 = iregFlows.getInterRegFlows(1, 3);
        BOOST_CHECK_MESSAGE(! q24.has_value(),
                            "There must NOT be inter-region "
                            "flows between regions 2 and 4");
    }

    {
        const auto q34 = iregFlows.getInterRegFlows(2, 3);
        BOOST_REQUIRE_MESSAGE(! q34.has_value(),
                              "There must NOT be inter-region flows "
                              "between regions 3 and 4");
    }
}

BOOST_AUTO_TEST_CASE(TopOther_BottomInternal)
{
    auto flows = Opm::InterRegFlowMapSingleFIP{ all_separate_region() };
    BOOST_CHECK_EQUAL(flows.getLocalMaxRegionID(), std::size_t{4});

    flows.addConnection({ 0, 0, false }, { 1, 1, false }, conn_01());
    flows.addConnection({ 0, 0, false }, { 2, 2, true } , conn_02());
    flows.addConnection({ 1, 1, false }, { 3, 3, true } , conn_13());
    flows.addConnection({ 2, 2, true } , { 3, 3, true } , conn_23());

    flows.compress();

    const auto& iregFlows = flows.getInterRegFlows();
    BOOST_CHECK_EQUAL(iregFlows.numRegions(), 4);

    {
        const auto q12 = iregFlows.getInterRegFlows(0, 1);
        BOOST_REQUIRE_MESSAGE(q12.has_value(),
                              "There must be inter-region flows "
                              "between regions 1 and 2");

        const auto& [rate, sign] = q12.value();

        BOOST_CHECK_EQUAL(sign, 1.0);

        using FlowView = std::remove_cv_t<std::remove_reference_t<
            decltype(rate)>>;

        using Component = FlowView::Component;
        using Direction = FlowView::Direction;

        BOOST_CHECK_CLOSE(rate.flow(Component::Oil), 0.12f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Gas), 0.24f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Water), 0.36f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Disgas), 0.48f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Vapoil), 1.25f, 1.0e-6);

        BOOST_CHECK_CLOSE(rate.flow(Component::Oil, Direction::Positive), 0.12f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Gas, Direction::Positive), 0.24f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Water, Direction::Positive), 0.36f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Disgas, Direction::Positive), 0.48f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Vapoil, Direction::Positive), 1.25f, 1.0e-6);

        BOOST_CHECK_CLOSE(rate.flow(Component::Oil, Direction::Negative), 0.0f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Gas, Direction::Negative), 0.0f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Water, Direction::Negative), 0.0f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Disgas, Direction::Negative), 0.0f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Vapoil, Direction::Negative), 0.0f, 1.0e-6);
    }

    {
        const auto q13 = iregFlows.getInterRegFlows(0, 2);
        BOOST_REQUIRE_MESSAGE(! q13.has_value(),
                              "There must NOT be inter-region flows "
                              "between regions 1 and 3");
    }

    {
        const auto q14 = iregFlows.getInterRegFlows(0, 3);
        BOOST_CHECK_MESSAGE(! q14.has_value(),
                            "There must NOT be inter-region "
                            "flows between regions 1 and 4");
    }

    {
        const auto q23 = iregFlows.getInterRegFlows(1, 2);
        BOOST_CHECK_MESSAGE(! q23.has_value(),
                            "There must NOT be inter-region "
                            "flows between regions 2 and 3");
    }

    {
        const auto q24 = iregFlows.getInterRegFlows(1, 3);
        BOOST_REQUIRE_MESSAGE(! q24.has_value(),
                              "There must NOT be inter-region flows "
                              "between regions 2 and 4");
    }

    {
        const auto q34 = iregFlows.getInterRegFlows(2, 3);
        BOOST_REQUIRE_MESSAGE(! q34.has_value(),
                              "There must NOT be inter-region "
                              "flows between regions 3 and 4");
    }
}

BOOST_AUTO_TEST_CASE(Checker_Board_Region)
{
    auto flows = Opm::InterRegFlowMapSingleFIP{ all_separate_region() };
    BOOST_CHECK_EQUAL(flows.getLocalMaxRegionID(), std::size_t{4});

    flows.addConnection({ 0, 0, false }, { 1, 1, true } , conn_01());
    flows.addConnection({ 0, 0, false }, { 2, 2, true } , conn_02());
    flows.addConnection({ 1, 1, true } , { 3, 3, false }, conn_13());
    flows.addConnection({ 2, 2, true } , { 3, 3, false }, conn_23());

    flows.compress();

    const auto& iregFlows = flows.getInterRegFlows();
    BOOST_CHECK_EQUAL(iregFlows.numRegions(), 4);

    {
        const auto q12 = iregFlows.getInterRegFlows(0, 1);
        BOOST_REQUIRE_MESSAGE(q12.has_value(),
                              "There must be inter-region flows "
                              "between regions 1 and 2");

        const auto& [rate, sign] = q12.value();

        BOOST_CHECK_EQUAL(sign, 1.0);

        using FlowView = std::remove_cv_t<std::remove_reference_t<
            decltype(rate)>>;

        using Component = FlowView::Component;
        using Direction = FlowView::Direction;

        BOOST_CHECK_CLOSE(rate.flow(Component::Oil), 0.12f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Gas), 0.24f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Water), 0.36f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Disgas), 0.48f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Vapoil), 1.25f, 1.0e-6);

        BOOST_CHECK_CLOSE(rate.flow(Component::Oil, Direction::Positive), 0.12f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Gas, Direction::Positive), 0.24f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Water, Direction::Positive), 0.36f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Disgas, Direction::Positive), 0.48f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Vapoil, Direction::Positive), 1.25f, 1.0e-6);

        BOOST_CHECK_CLOSE(rate.flow(Component::Oil, Direction::Negative), 0.0f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Gas, Direction::Negative), 0.0f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Water, Direction::Negative), 0.0f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Disgas, Direction::Negative), 0.0f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Vapoil, Direction::Negative), 0.0f, 1.0e-6);
    }

    {
        const auto q13 = iregFlows.getInterRegFlows(0, 2);
        BOOST_REQUIRE_MESSAGE(q13.has_value(),
                              "There must be inter-region flows "
                              "between regions 1 and 3");

        const auto& [rate, sign] = q13.value();

        BOOST_CHECK_EQUAL(sign, 1.0);

        using FlowView = std::remove_cv_t<std::remove_reference_t<
            decltype(rate)>>;

        using Component = FlowView::Component;
        using Direction = FlowView::Direction;

        BOOST_CHECK_CLOSE(rate.flow(Component::Oil), 0.2f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Gas), 0.4f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Water), 0.6f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Disgas), 0.8f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Vapoil), 1.0f, 1.0e-6);

        BOOST_CHECK_CLOSE(rate.flow(Component::Oil, Direction::Positive), 0.2f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Gas, Direction::Positive), 0.4f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Water, Direction::Positive), 0.6f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Disgas, Direction::Positive), 0.8f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Vapoil, Direction::Positive), 1.0f, 1.0e-6);

        BOOST_CHECK_CLOSE(rate.flow(Component::Oil, Direction::Negative), 0.0f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Gas, Direction::Negative), 0.0f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Water, Direction::Negative), 0.0f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Disgas, Direction::Negative), 0.0f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Vapoil, Direction::Negative), 0.0f, 1.0e-6);
    }

    {
        const auto q14 = iregFlows.getInterRegFlows(0, 3);
        BOOST_CHECK_MESSAGE(! q14.has_value(),
                            "There must NOT be inter-region "
                            "flows between regions 1 and 4");
    }

    {
        const auto q23 = iregFlows.getInterRegFlows(1, 2);
        BOOST_CHECK_MESSAGE(! q23.has_value(),
                            "There must NOT be inter-region "
                            "flows between regions 2 and 3");
    }

    {
        const auto q24 = iregFlows.getInterRegFlows(1, 3);
        BOOST_REQUIRE_MESSAGE(! q24.has_value(),
                              "There must NOT be inter-region "
                              "flows between regions 2 and 4");
    }

    {
        const auto q34 = iregFlows.getInterRegFlows(2, 3);
        BOOST_REQUIRE_MESSAGE(! q34.has_value(),
                              "There must NOT be inter-region "
                              "flows between regions 3 and 4");
    }
}

BOOST_AUTO_TEST_CASE(Reverse_CheckerBoard_Internal)
{
    auto flows = Opm::InterRegFlowMapSingleFIP{ all_separate_region() };
    BOOST_CHECK_EQUAL(flows.getLocalMaxRegionID(), std::size_t{4});

    flows.addConnection({ 0, 0, true } , { 1, 1, false }, conn_01());
    flows.addConnection({ 0, 0, true } , { 2, 2, false }, conn_02());
    flows.addConnection({ 1, 1, false }, { 3, 3, true } , conn_13());
    flows.addConnection({ 2, 2, false }, { 3, 3, true } , conn_23());

    flows.compress();

    const auto& iregFlows = flows.getInterRegFlows();
    BOOST_CHECK_EQUAL(iregFlows.numRegions(), 4);

    {
        const auto q12 = iregFlows.getInterRegFlows(0, 1);
        BOOST_REQUIRE_MESSAGE(! q12.has_value(),
                              "There must NOT be inter-region "
                              "flows between regions 1 and 2");
    }

    {
        const auto q13 = iregFlows.getInterRegFlows(0, 2);
        BOOST_REQUIRE_MESSAGE(! q13.has_value(),
                              "There must NOT be inter-region "
                              "flows between regions 1 and 3");
    }

    {
        const auto q14 = iregFlows.getInterRegFlows(0, 3);
        BOOST_CHECK_MESSAGE(! q14.has_value(),
                            "There must NOT be inter-region "
                            "flows between regions 1 and 4");
    }

    {
        const auto q23 = iregFlows.getInterRegFlows(1, 2);
        BOOST_CHECK_MESSAGE(! q23.has_value(),
                            "There must NOT be inter-region "
                            "flows between regions 2 and 3");
    }

    {
        const auto q24 = iregFlows.getInterRegFlows(1, 3);
        BOOST_CHECK_MESSAGE(q24.has_value(),
                            "There must be inter-region "
                            "flows between regions 2 and 4");

        const auto& [rate, sign] = q24.value();

        BOOST_CHECK_EQUAL(sign, 1.0);

        using FlowView = std::remove_cv_t<std::remove_reference_t<
            decltype(rate)>>;

        using Component = FlowView::Component;
        using Direction = FlowView::Direction;

        BOOST_CHECK_CLOSE(rate.flow(Component::Oil), -0.1f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Gas), -0.2f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Water), -0.3f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Disgas), -0.4f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Vapoil), -0.5f, 1.0e-6);

        BOOST_CHECK_CLOSE(rate.flow(Component::Oil, Direction::Positive), 0.0f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Gas, Direction::Positive), 0.0f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Water, Direction::Positive), 0.0f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Disgas, Direction::Positive), 0.0f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Vapoil, Direction::Positive), 0.0f, 1.0e-6);

        BOOST_CHECK_CLOSE(rate.flow(Component::Oil, Direction::Negative), -0.1f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Gas, Direction::Negative), -0.2f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Water, Direction::Negative), -0.3f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Disgas, Direction::Negative), -0.4f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Vapoil, Direction::Negative), -0.5f, 1.0e-6);
    }

    {
        const auto q34 = iregFlows.getInterRegFlows(2, 3);
        BOOST_REQUIRE_MESSAGE(q34.has_value(),
                              "There must be inter-region "
                              "flows between regions 3 and 4");

        const auto& [rate, sign] = q34.value();

        BOOST_CHECK_EQUAL(sign, 1.0);

        using FlowView = std::remove_cv_t<std::remove_reference_t<
            decltype(rate)>>;

        using Component = FlowView::Component;
        using Direction = FlowView::Direction;

        BOOST_CHECK_CLOSE(rate.flow(Component::Oil), -1.0f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Gas), -2.0f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Water), -3.0f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Disgas), -4.0f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Vapoil), -5.0f, 1.0e-6);

        BOOST_CHECK_CLOSE(rate.flow(Component::Oil, Direction::Positive), 0.0f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Gas, Direction::Positive), 0.0f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Water, Direction::Positive), 0.0f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Disgas, Direction::Positive), 0.0f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Vapoil, Direction::Positive), 0.0f, 1.0e-6);

        BOOST_CHECK_CLOSE(rate.flow(Component::Oil, Direction::Negative), -1.0f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Gas, Direction::Negative), -2.0f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Water, Direction::Negative), -3.0f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Disgas, Direction::Negative), -4.0f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Vapoil, Direction::Negative), -5.0f, 1.0e-6);
    }
}

BOOST_AUTO_TEST_SUITE_END() // All_Separate_Region

// =====================================================================

BOOST_AUTO_TEST_SUITE(MultiRank_ReadWriteMerge)

namespace {
    auto makeMessageBuffer()
    {
        return Dune::Point2PointCommunicator<Dune::SimpleMessageBuffer>
            ::MessageBufferType{};
    }

    template <typename Predicate>
    void addConnections(Predicate&&                    isInterior,
                        Opm::InterRegFlowMapSingleFIP& rank)
    {
        rank.addConnection({ 0, 0, isInterior(0) }, { 1, 1, isInterior(1) }, conn_01());
        rank.addConnection({ 0, 0, isInterior(0) }, { 2, 2, isInterior(2) }, conn_02());
        rank.addConnection({ 1, 1, isInterior(1) }, { 3, 3, isInterior(3) }, conn_13());
        rank.addConnection({ 2, 2, isInterior(2) }, { 3, 3, isInterior(3) }, conn_23());
    }
}

BOOST_AUTO_TEST_SUITE(Two_Processes)

BOOST_AUTO_TEST_CASE(Empty)
{
    using Map = Opm::InterRegFlowMapSingleFIP;

    auto rank = std::vector<Map>(2, Map { all_same_region() });

    BOOST_CHECK_EQUAL(rank[0].getLocalMaxRegionID(), std::size_t{1});
    BOOST_CHECK_EQUAL(rank[1].getLocalMaxRegionID(), std::size_t{1});

    rank[0].assignGlobalMaxRegionID(3);
    rank[1].assignGlobalMaxRegionID(5);

    for (auto& r : rank) {
        r.compress();
    }

    auto buffer = makeMessageBuffer();

    rank[1].write(buffer);
    rank[0].read(buffer);

    BOOST_CHECK_THROW(rank[0].addConnection({ 0, TwoProc::P1::isInterior(0) },
                                            { 1, TwoProc::P1::isInterior(1) }, conn_01()),
                      std::logic_error);

    rank[0].compress();

    const auto& iregFlows = rank[0].getInterRegFlows();
    BOOST_CHECK_EQUAL(iregFlows.numRegions(), 5);

    for (auto r1 = 0*iregFlows.numRegions(); r1 < iregFlows.numRegions(); ++r1) {
        for (auto r2 = r1 + 1; r2 < iregFlows.numRegions(); ++r2) {
            auto flow = iregFlows.getInterRegFlows(r1, r2);
            BOOST_CHECK_MESSAGE(! flow.has_value(),
                                "There must not be inter-regional flow "
                                "between regions " << r1 << " and " << r2);
        }
    }
}

BOOST_AUTO_TEST_CASE(Single_FIP_Region)
{
    using Map = Opm::InterRegFlowMapSingleFIP;

    auto rank = std::vector<Map>(2, Map { all_same_region() });

    BOOST_CHECK_EQUAL(rank[0].getLocalMaxRegionID(), std::size_t{1});
    BOOST_CHECK_EQUAL(rank[1].getLocalMaxRegionID(), std::size_t{1});

    addConnections(TwoProc::P1::isInterior, rank[0]);
    rank[0].assignGlobalMaxRegionID(3);

    addConnections(TwoProc::P2::isInterior, rank[1]);
    rank[1].assignGlobalMaxRegionID(5);

    for (auto& r : rank) {
        r.compress();
    }

    auto buffer = makeMessageBuffer();

    rank[1].write(buffer);
    rank[0].read(buffer);

    BOOST_CHECK_THROW(rank[0].addConnection({ 0, TwoProc::P1::isInterior(0) },
                                            { 1, TwoProc::P1::isInterior(1) }, conn_01()),
                      std::logic_error);

    rank[0].compress();

    const auto& iregFlows = rank[0].getInterRegFlows();
    BOOST_CHECK_EQUAL(iregFlows.numRegions(), 5);

    for (auto r1 = 0*iregFlows.numRegions(); r1 < iregFlows.numRegions(); ++r1) {
        for (auto r2 = r1 + 1; r2 < iregFlows.numRegions(); ++r2) {
            auto flow = iregFlows.getInterRegFlows(r1, r2);
            BOOST_CHECK_MESSAGE(! flow.has_value(),
                                "There must not be inter-regional flow "
                                "between regions " << r1 << " and " << r2);
        }
    }
}

BOOST_AUTO_TEST_CASE(Left_Right_Split_Region)
{
    using Map = Opm::InterRegFlowMapSingleFIP;

    auto rank = std::vector<Map>(2, Map { left_right_split_region() });

    BOOST_CHECK_EQUAL(rank[0].getLocalMaxRegionID(), std::size_t{2});
    BOOST_CHECK_EQUAL(rank[1].getLocalMaxRegionID(), std::size_t{2});

    addConnections(TwoProc::P1::isInterior, rank[0]);
    rank[0].assignGlobalMaxRegionID(3);

    addConnections(TwoProc::P2::isInterior, rank[1]);
    rank[1].assignGlobalMaxRegionID(5);

    for (auto& r : rank) {
        r.compress();
    }

    auto buffer = makeMessageBuffer();

    rank[1].write(buffer);
    rank[0].read(buffer);

    rank[0].compress();

    const auto& iregFlows = rank[0].getInterRegFlows();
    BOOST_CHECK_EQUAL(iregFlows.numRegions(), 5);
    {
        const auto q12 = iregFlows.getInterRegFlows(0, 1);
        BOOST_REQUIRE_MESSAGE(q12.has_value(),
                              "There must be inter-region flows "
                              "between regions 1 and 2");

        const auto& [rate, sign] = q12.value();

        BOOST_CHECK_EQUAL(sign, 1.0);

        using FlowView = std::remove_cv_t<std::remove_reference_t<
            decltype(rate)>>;

        using Component = FlowView::Component;
        using Direction = FlowView::Direction;

        BOOST_CHECK_CLOSE(rate.flow(Component::Oil), 0.88f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Gas), 1.76f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Water), 2.64f, 1.0e-5);
        BOOST_CHECK_CLOSE(rate.flow(Component::Disgas), 3.52f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Vapoil), 3.75f, 1.0e-6);

        BOOST_CHECK_CLOSE(rate.flow(Component::Oil, Direction::Positive), 1.0f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Gas, Direction::Positive), 2.0f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Water, Direction::Positive), 3.0f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Disgas, Direction::Positive), 4.0f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Vapoil, Direction::Positive), 5.0f, 1.0e-6);

        BOOST_CHECK_CLOSE(rate.flow(Component::Oil, Direction::Negative), -0.12f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Gas, Direction::Negative), -0.24f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Water, Direction::Negative), -0.36f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Disgas, Direction::Negative), -0.48f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Vapoil, Direction::Negative), -1.25f, 1.0e-6);
    }
}

BOOST_AUTO_TEST_CASE(Checker_Board_Region)
{
    using Map = Opm::InterRegFlowMapSingleFIP;

    auto rank = std::vector<Map>(2, Map { checker_board_region() });

    BOOST_CHECK_EQUAL(rank[0].getLocalMaxRegionID(), std::size_t{2});
    BOOST_CHECK_EQUAL(rank[1].getLocalMaxRegionID(), std::size_t{2});

    addConnections(TwoProc::P1::isInterior, rank[0]);
    rank[0].assignGlobalMaxRegionID(3);

    addConnections(TwoProc::P2::isInterior, rank[1]);
    rank[1].assignGlobalMaxRegionID(5);

    for (auto& r : rank) {
        r.compress();
    }

    auto buffer = makeMessageBuffer();

    rank[1].write(buffer);
    rank[0].read(buffer);

    rank[0].compress();

    const auto& iregFlows = rank[0].getInterRegFlows();
    BOOST_CHECK_EQUAL(iregFlows.numRegions(), 5);

    {
        const auto q12 = iregFlows.getInterRegFlows(0, 1);
        BOOST_REQUIRE_MESSAGE(q12.has_value(),
                              "There must be inter-region flows "
                              "between regions 1 and 2");

        const auto& [rate, sign] = q12.value();

        BOOST_CHECK_EQUAL(sign, 1.0);

        using FlowView = std::remove_cv_t<std::remove_reference_t<
            decltype(rate)>>;

        using Component = FlowView::Component;
        using Direction = FlowView::Direction;

        BOOST_CHECK_CLOSE(rate.flow(Component::Oil), 1.42f, 1.0e-5);
        BOOST_CHECK_CLOSE(rate.flow(Component::Gas), 2.84f, 1.0e-5);
        BOOST_CHECK_CLOSE(rate.flow(Component::Water), 4.26f, 5.0e-5);
        BOOST_CHECK_CLOSE(rate.flow(Component::Disgas), 5.68f, 1.0e-5);
        BOOST_CHECK_CLOSE(rate.flow(Component::Vapoil), 7.75f, 1.0e-6);

        BOOST_CHECK_CLOSE(rate.flow(Component::Oil, Direction::Positive), 1.42f, 1.0e-5);
        BOOST_CHECK_CLOSE(rate.flow(Component::Gas, Direction::Positive), 2.84f, 1.0e-5);
        BOOST_CHECK_CLOSE(rate.flow(Component::Water, Direction::Positive), 4.26f, 5.0e-5);
        BOOST_CHECK_CLOSE(rate.flow(Component::Disgas, Direction::Positive), 5.68f, 1.0e-5);
        BOOST_CHECK_CLOSE(rate.flow(Component::Vapoil, Direction::Positive), 7.75f, 1.0e-6);

        BOOST_CHECK_CLOSE(rate.flow(Component::Oil, Direction::Negative), 0.0f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Gas, Direction::Negative), 0.0f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Water, Direction::Negative), 0.0f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Disgas, Direction::Negative), 0.0f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Vapoil, Direction::Negative), 0.0f, 1.0e-6);
    }
}

BOOST_AUTO_TEST_CASE(All_Separate_Region)
{
    using Map = Opm::InterRegFlowMapSingleFIP;

    auto rank = std::vector<Map>(2, Map { all_separate_region() });

    BOOST_CHECK_EQUAL(rank[0].getLocalMaxRegionID(), std::size_t{4});
    BOOST_CHECK_EQUAL(rank[1].getLocalMaxRegionID(), std::size_t{4});

    addConnections(TwoProc::P1::isInterior, rank[0]);
    BOOST_CHECK_MESSAGE(! rank[0].assignGlobalMaxRegionID(3),
                        "Assigning small maximum global region ID must NOT succeed");

    addConnections(TwoProc::P2::isInterior, rank[1]);
    rank[1].assignGlobalMaxRegionID(5);

    for (auto& r : rank) {
        r.compress();
    }

    auto buffer = makeMessageBuffer();

    rank[1].write(buffer);
    rank[0].read(buffer);

    rank[0].compress();

    const auto& iregFlows = rank[0].getInterRegFlows();
    BOOST_CHECK_EQUAL(iregFlows.numRegions(), 5);

    {
        const auto q12 = iregFlows.getInterRegFlows(0, 1);
        BOOST_REQUIRE_MESSAGE(q12.has_value(),
                              "There must be inter-region flows "
                              "between regions 1 and 2");

        const auto& [rate, sign] = q12.value();

        BOOST_CHECK_EQUAL(sign, 1.0);

        using FlowView = std::remove_cv_t<std::remove_reference_t<
            decltype(rate)>>;

        using Component = FlowView::Component;
        using Direction = FlowView::Direction;

        BOOST_CHECK_CLOSE(rate.flow(Component::Oil), 0.12f, 1.0e-5);
        BOOST_CHECK_CLOSE(rate.flow(Component::Gas), 0.24f, 1.0e-5);
        BOOST_CHECK_CLOSE(rate.flow(Component::Water), 0.36f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Disgas), 0.48f, 1.0e-5);
        BOOST_CHECK_CLOSE(rate.flow(Component::Vapoil), 1.25f, 1.0e-6);

        BOOST_CHECK_CLOSE(rate.flow(Component::Oil, Direction::Positive), 0.12f, 1.0e-5);
        BOOST_CHECK_CLOSE(rate.flow(Component::Gas, Direction::Positive), 0.24f, 1.0e-5);
        BOOST_CHECK_CLOSE(rate.flow(Component::Water, Direction::Positive), 0.36f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Disgas, Direction::Positive), 0.48f, 1.0e-5);
        BOOST_CHECK_CLOSE(rate.flow(Component::Vapoil, Direction::Positive), 1.25f, 1.0e-6);

        BOOST_CHECK_CLOSE(rate.flow(Component::Oil, Direction::Negative), 0.0f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Gas, Direction::Negative), 0.0f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Water, Direction::Negative), 0.0f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Disgas, Direction::Negative), 0.0f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Vapoil, Direction::Negative), 0.0f, 1.0e-6);
    }

    {
        const auto q13 = iregFlows.getInterRegFlows(0, 2);
        BOOST_REQUIRE_MESSAGE(q13.has_value(),
                              "There must be inter-region flows "
                              "between regions 1 and 3");

        const auto& [rate, sign] = q13.value();

        BOOST_CHECK_EQUAL(sign, 1.0);

        using FlowView = std::remove_cv_t<std::remove_reference_t<
            decltype(rate)>>;

        using Component = FlowView::Component;
        using Direction = FlowView::Direction;

        BOOST_CHECK_CLOSE(rate.flow(Component::Oil), 0.2f, 1.0e-5);
        BOOST_CHECK_CLOSE(rate.flow(Component::Gas), 0.4f, 1.0e-5);
        BOOST_CHECK_CLOSE(rate.flow(Component::Water), 0.6f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Disgas), 0.8f, 1.0e-5);
        BOOST_CHECK_CLOSE(rate.flow(Component::Vapoil), 1.0f, 1.0e-6);

        BOOST_CHECK_CLOSE(rate.flow(Component::Oil, Direction::Positive), 0.2f, 1.0e-5);
        BOOST_CHECK_CLOSE(rate.flow(Component::Gas, Direction::Positive), 0.4f, 1.0e-5);
        BOOST_CHECK_CLOSE(rate.flow(Component::Water, Direction::Positive), 0.6f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Disgas, Direction::Positive), 0.8f, 1.0e-5);
        BOOST_CHECK_CLOSE(rate.flow(Component::Vapoil, Direction::Positive), 1.0f, 1.0e-6);

        BOOST_CHECK_CLOSE(rate.flow(Component::Oil, Direction::Negative), 0.0f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Gas, Direction::Negative), 0.0f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Water, Direction::Negative), 0.0f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Disgas, Direction::Negative), 0.0f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Vapoil, Direction::Negative), 0.0f, 1.0e-6);
    }

    {
        const auto q14 = iregFlows.getInterRegFlows(0, 3);
        BOOST_CHECK_MESSAGE(! q14.has_value(),
                            "There must NOT be inter-region "
                            "flows between regions 1 and 4");
    }

    {
        const auto q23 = iregFlows.getInterRegFlows(1, 2);
        BOOST_CHECK_MESSAGE(! q23.has_value(),
                            "There must NOT be inter-region "
                            "flows between regions 2 and 3");
    }

    {
        const auto q24 = iregFlows.getInterRegFlows(1, 3);
        BOOST_REQUIRE_MESSAGE(q24.has_value(),
                              "There must be inter-region flows "
                              "between regions 2 and 4");

        const auto& [rate, sign] = q24.value();

        BOOST_CHECK_EQUAL(sign, 1.0);

        using FlowView = std::remove_cv_t<std::remove_reference_t<
            decltype(rate)>>;

        using Component = FlowView::Component;
        using Direction = FlowView::Direction;

        BOOST_CHECK_CLOSE(rate.flow(Component::Oil), -0.1f, 1.0e-5);
        BOOST_CHECK_CLOSE(rate.flow(Component::Gas), -0.2f, 1.0e-5);
        BOOST_CHECK_CLOSE(rate.flow(Component::Water), -0.3f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Disgas), -0.4f, 1.0e-5);
        BOOST_CHECK_CLOSE(rate.flow(Component::Vapoil), -0.5f, 1.0e-6);

        BOOST_CHECK_CLOSE(rate.flow(Component::Oil, Direction::Positive), 0.0f, 1.0e-5);
        BOOST_CHECK_CLOSE(rate.flow(Component::Gas, Direction::Positive), 0.0f, 1.0e-5);
        BOOST_CHECK_CLOSE(rate.flow(Component::Water, Direction::Positive), 0.0f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Disgas, Direction::Positive), 0.0f, 1.0e-5);
        BOOST_CHECK_CLOSE(rate.flow(Component::Vapoil, Direction::Positive), 0.0f, 1.0e-6);

        BOOST_CHECK_CLOSE(rate.flow(Component::Oil, Direction::Negative), -0.1f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Gas, Direction::Negative), -0.2f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Water, Direction::Negative), -0.3f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Disgas, Direction::Negative), -0.4f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Vapoil, Direction::Negative), -0.5f, 1.0e-6);
    }
}

BOOST_AUTO_TEST_SUITE_END() // Two_Processes

// ---------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE(Four_Processes)

BOOST_AUTO_TEST_CASE(Single_FIP_Region)
{
    using Map = Opm::InterRegFlowMapSingleFIP;

    auto rank = std::vector<Map>(4, Map { all_same_region() });
    BOOST_CHECK_EQUAL(rank[0].getLocalMaxRegionID(), std::size_t{1});
    BOOST_CHECK_EQUAL(rank[1].getLocalMaxRegionID(), std::size_t{1});
    BOOST_CHECK_EQUAL(rank[2].getLocalMaxRegionID(), std::size_t{1});
    BOOST_CHECK_EQUAL(rank[3].getLocalMaxRegionID(), std::size_t{1});

    addConnections(FourProc::P1::isInterior, rank[0]);
    rank[0].assignGlobalMaxRegionID(3);

    addConnections(FourProc::P2::isInterior, rank[1]);
    rank[1].assignGlobalMaxRegionID(5);

    addConnections(FourProc::P3::isInterior, rank[2]);
    rank[2].assignGlobalMaxRegionID(2);

    addConnections(FourProc::P4::isInterior, rank[3]);

    for (auto& r : rank) {
        r.compress();
    }
    auto buffer = makeMessageBuffer();

    rank[1].write(buffer);
    rank[3].write(buffer);
    rank[2].write(buffer);
    rank[0].read(buffer);
    rank[0].read(buffer);
    rank[0].read(buffer);

    BOOST_CHECK_THROW(rank[0].addConnection({ 0, FourProc::P1::isInterior(0) },
                                            { 1, FourProc::P1::isInterior(1) }, conn_01()),
                      std::logic_error);

    rank[0].compress();

    const auto& iregFlows = rank[0].getInterRegFlows();
    BOOST_CHECK_EQUAL(iregFlows.numRegions(), 5);

    for (auto r1 = 0*iregFlows.numRegions(); r1 < iregFlows.numRegions(); ++r1) {
        for (auto r2 = r1 + 1; r2 < iregFlows.numRegions(); ++r2) {
            auto flow = iregFlows.getInterRegFlows(r1, r2);
            BOOST_CHECK_MESSAGE(! flow.has_value(),
                                "There must not be inter-regional flow "
                                "between regions " << r1 << " and " << r2);
        }
    }
}

BOOST_AUTO_TEST_CASE(Left_Right_Split_Region)
{
    using Map = Opm::InterRegFlowMapSingleFIP;

    auto rank = std::vector<Map>(4, Map { left_right_split_region() });
    BOOST_CHECK_EQUAL(rank[0].getLocalMaxRegionID(), std::size_t{2});
    BOOST_CHECK_EQUAL(rank[1].getLocalMaxRegionID(), std::size_t{2});
    BOOST_CHECK_EQUAL(rank[2].getLocalMaxRegionID(), std::size_t{2});
    BOOST_CHECK_EQUAL(rank[3].getLocalMaxRegionID(), std::size_t{2});

    addConnections(FourProc::P1::isInterior, rank[0]);
    rank[0].assignGlobalMaxRegionID(3);

    addConnections(FourProc::P2::isInterior, rank[1]);
    rank[1].assignGlobalMaxRegionID(5);

    addConnections(FourProc::P3::isInterior, rank[2]);
    rank[2].assignGlobalMaxRegionID(2);

    addConnections(FourProc::P4::isInterior, rank[3]);

    for (auto& r : rank) {
        r.compress();
    }

    auto buffer = makeMessageBuffer();

    rank[0].write(buffer);
    rank[3].write(buffer);
    rank[2].write(buffer);
    rank[1].read(buffer);
    rank[1].read(buffer);
    rank[1].read(buffer);

    rank[1].compress();

    const auto& iregFlows = rank[1].getInterRegFlows();
    BOOST_CHECK_EQUAL(iregFlows.numRegions(), 5);
    {
        const auto q12 = iregFlows.getInterRegFlows(0, 1);
        BOOST_REQUIRE_MESSAGE(q12.has_value(),
                              "There must be inter-region flows "
                              "between regions 1 and 2");

        const auto& [rate, sign] = q12.value();

        BOOST_CHECK_EQUAL(sign, 1.0);

        using FlowView = std::remove_cv_t<std::remove_reference_t<
            decltype(rate)>>;

        using Component = FlowView::Component;
        using Direction = FlowView::Direction;

        BOOST_CHECK_CLOSE(rate.flow(Component::Oil), 0.88f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Gas), 1.76f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Water), 2.64f, 1.0e-5);
        BOOST_CHECK_CLOSE(rate.flow(Component::Disgas), 3.52f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Vapoil), 3.75f, 1.0e-6);

        BOOST_CHECK_CLOSE(rate.flow(Component::Oil, Direction::Positive), 1.0f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Gas, Direction::Positive), 2.0f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Water, Direction::Positive), 3.0f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Disgas, Direction::Positive), 4.0f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Vapoil, Direction::Positive), 5.0f, 1.0e-6);

        BOOST_CHECK_CLOSE(rate.flow(Component::Oil, Direction::Negative), -0.12f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Gas, Direction::Negative), -0.24f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Water, Direction::Negative), -0.36f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Disgas, Direction::Negative), -0.48f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Vapoil, Direction::Negative), -1.25f, 1.0e-6);
    }
}

BOOST_AUTO_TEST_CASE(Checker_Board_Region)
{
    using Map = Opm::InterRegFlowMapSingleFIP;

    auto rank = std::vector<Map>(4, Map { checker_board_region() });
    BOOST_CHECK_EQUAL(rank[0].getLocalMaxRegionID(), std::size_t{2});
    BOOST_CHECK_EQUAL(rank[1].getLocalMaxRegionID(), std::size_t{2});
    BOOST_CHECK_EQUAL(rank[2].getLocalMaxRegionID(), std::size_t{2});
    BOOST_CHECK_EQUAL(rank[3].getLocalMaxRegionID(), std::size_t{2});

    addConnections(FourProc::P1::isInterior, rank[0]);
    rank[0].assignGlobalMaxRegionID(3);

    addConnections(FourProc::P2::isInterior, rank[1]);
    rank[1].assignGlobalMaxRegionID(5);

    addConnections(FourProc::P3::isInterior, rank[2]);
    rank[2].assignGlobalMaxRegionID(2);

    addConnections(FourProc::P4::isInterior, rank[3]);

    for (auto& r : rank) {
        r.compress();
    }

    auto buffer = makeMessageBuffer();

    rank[2].write(buffer);
    rank[0].write(buffer);
    rank[1].write(buffer);
    rank[3].read(buffer);
    rank[3].read(buffer);
    rank[3].read(buffer);

    rank[3].compress();

    const auto& iregFlows = rank[3].getInterRegFlows();
    BOOST_CHECK_EQUAL(iregFlows.numRegions(), 5);

    {
        const auto q12 = iregFlows.getInterRegFlows(0, 1);
        BOOST_REQUIRE_MESSAGE(q12.has_value(),
                              "There must be inter-region flows "
                              "between regions 1 and 2");

        const auto& [rate, sign] = q12.value();

        BOOST_CHECK_EQUAL(sign, 1.0);

        using FlowView = std::remove_cv_t<std::remove_reference_t<
            decltype(rate)>>;

        using Component = FlowView::Component;
        using Direction = FlowView::Direction;

        BOOST_CHECK_CLOSE(rate.flow(Component::Oil), 1.42f, 1.0e-5);
        BOOST_CHECK_CLOSE(rate.flow(Component::Gas), 2.84f, 1.0e-5);
        BOOST_CHECK_CLOSE(rate.flow(Component::Water), 4.26f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Disgas), 5.68f, 1.0e-5);
        BOOST_CHECK_CLOSE(rate.flow(Component::Vapoil), 7.75f, 1.0e-6);

        BOOST_CHECK_CLOSE(rate.flow(Component::Oil, Direction::Positive), 1.42f, 1.0e-5);
        BOOST_CHECK_CLOSE(rate.flow(Component::Gas, Direction::Positive), 2.84f, 1.0e-5);
        BOOST_CHECK_CLOSE(rate.flow(Component::Water, Direction::Positive), 4.26f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Disgas, Direction::Positive), 5.68f, 1.0e-5);
        BOOST_CHECK_CLOSE(rate.flow(Component::Vapoil, Direction::Positive), 7.75f, 1.0e-6);

        BOOST_CHECK_CLOSE(rate.flow(Component::Oil, Direction::Negative), 0.0f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Gas, Direction::Negative), 0.0f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Water, Direction::Negative), 0.0f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Disgas, Direction::Negative), 0.0f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Vapoil, Direction::Negative), 0.0f, 1.0e-6);
    }
}

BOOST_AUTO_TEST_CASE(All_Separate_Region)
{
    using Map = Opm::InterRegFlowMapSingleFIP;

    auto rank = std::vector<Map>(4, Map { all_separate_region() });
    BOOST_CHECK_EQUAL(rank[0].getLocalMaxRegionID(), std::size_t{4});
    BOOST_CHECK_EQUAL(rank[1].getLocalMaxRegionID(), std::size_t{4});
    BOOST_CHECK_EQUAL(rank[2].getLocalMaxRegionID(), std::size_t{4});
    BOOST_CHECK_EQUAL(rank[3].getLocalMaxRegionID(), std::size_t{4});

    addConnections(FourProc::P1::isInterior, rank[0]);
    BOOST_CHECK_MESSAGE(! rank[0].assignGlobalMaxRegionID(3),
                        "Assgning small maximum global region ID must NOT succeed");

    addConnections(FourProc::P2::isInterior, rank[1]);
    rank[1].assignGlobalMaxRegionID(5);

    addConnections(FourProc::P3::isInterior, rank[2]);
    BOOST_CHECK_MESSAGE(! rank[2].assignGlobalMaxRegionID(2),
                        "Assgning small maximum global region ID must NOT succeed");

    addConnections(FourProc::P4::isInterior, rank[3]);

    for (auto& r : rank) {
        r.compress();
    }

    auto buffer = makeMessageBuffer();

    rank[2].write(buffer);
    rank[3].write(buffer);
    rank[1].write(buffer);
    rank[0].read(buffer);
    rank[0].read(buffer);
    rank[0].read(buffer);

    rank[0].compress();

    const auto& iregFlows = rank[0].getInterRegFlows();
    BOOST_CHECK_EQUAL(iregFlows.numRegions(), 5);

    {
        const auto q12 = iregFlows.getInterRegFlows(0, 1);
        BOOST_REQUIRE_MESSAGE(q12.has_value(),
                              "There must be inter-region flows "
                              "between regions 1 and 2");

        const auto& [rate, sign] = q12.value();

        BOOST_CHECK_EQUAL(sign, 1.0);

        using FlowView = std::remove_cv_t<std::remove_reference_t<
            decltype(rate)>>;

        using Component = FlowView::Component;
        using Direction = FlowView::Direction;

        BOOST_CHECK_CLOSE(rate.flow(Component::Oil), 0.12f, 1.0e-5);
        BOOST_CHECK_CLOSE(rate.flow(Component::Gas), 0.24f, 1.0e-5);
        BOOST_CHECK_CLOSE(rate.flow(Component::Water), 0.36f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Disgas), 0.48f, 1.0e-5);
        BOOST_CHECK_CLOSE(rate.flow(Component::Vapoil), 1.25f, 1.0e-6);

        BOOST_CHECK_CLOSE(rate.flow(Component::Oil, Direction::Positive), 0.12f, 1.0e-5);
        BOOST_CHECK_CLOSE(rate.flow(Component::Gas, Direction::Positive), 0.24f, 1.0e-5);
        BOOST_CHECK_CLOSE(rate.flow(Component::Water, Direction::Positive), 0.36f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Disgas, Direction::Positive), 0.48f, 1.0e-5);
        BOOST_CHECK_CLOSE(rate.flow(Component::Vapoil, Direction::Positive), 1.25f, 1.0e-6);

        BOOST_CHECK_CLOSE(rate.flow(Component::Oil, Direction::Negative), 0.0f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Gas, Direction::Negative), 0.0f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Water, Direction::Negative), 0.0f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Disgas, Direction::Negative), 0.0f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Vapoil, Direction::Negative), 0.0f, 1.0e-6);
    }

    {
        const auto q13 = iregFlows.getInterRegFlows(0, 2);
        BOOST_REQUIRE_MESSAGE(q13.has_value(),
                              "There must be inter-region flows "
                              "between regions 1 and 3");

        const auto& [rate, sign] = q13.value();

        BOOST_CHECK_EQUAL(sign, 1.0);

        using FlowView = std::remove_cv_t<std::remove_reference_t<
            decltype(rate)>>;

        using Component = FlowView::Component;
        using Direction = FlowView::Direction;

        BOOST_CHECK_CLOSE(rate.flow(Component::Oil), 0.2f, 1.0e-5);
        BOOST_CHECK_CLOSE(rate.flow(Component::Gas), 0.4f, 1.0e-5);
        BOOST_CHECK_CLOSE(rate.flow(Component::Water), 0.6f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Disgas), 0.8f, 1.0e-5);
        BOOST_CHECK_CLOSE(rate.flow(Component::Vapoil), 1.0f, 1.0e-6);

        BOOST_CHECK_CLOSE(rate.flow(Component::Oil, Direction::Positive), 0.2f, 1.0e-5);
        BOOST_CHECK_CLOSE(rate.flow(Component::Gas, Direction::Positive), 0.4f, 1.0e-5);
        BOOST_CHECK_CLOSE(rate.flow(Component::Water, Direction::Positive), 0.6f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Disgas, Direction::Positive), 0.8f, 1.0e-5);
        BOOST_CHECK_CLOSE(rate.flow(Component::Vapoil, Direction::Positive), 1.0f, 1.0e-6);

        BOOST_CHECK_CLOSE(rate.flow(Component::Oil, Direction::Negative), 0.0f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Gas, Direction::Negative), 0.0f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Water, Direction::Negative), 0.0f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Disgas, Direction::Negative), 0.0f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Vapoil, Direction::Negative), 0.0f, 1.0e-6);
    }

    {
        const auto q14 = iregFlows.getInterRegFlows(0, 3);
        BOOST_CHECK_MESSAGE(! q14.has_value(),
                            "There must NOT be inter-region "
                            "flows between regions 1 and 4");
    }

    {
        const auto q23 = iregFlows.getInterRegFlows(1, 2);
        BOOST_CHECK_MESSAGE(! q23.has_value(),
                            "There must NOT be inter-region "
                            "flows between regions 2 and 3");
    }

    {
        const auto q24 = iregFlows.getInterRegFlows(1, 3);
        BOOST_REQUIRE_MESSAGE(q24.has_value(),
                              "There must be inter-region flows "
                              "between regions 2 and 4");

        const auto& [rate, sign] = q24.value();

        BOOST_CHECK_EQUAL(sign, 1.0);

        using FlowView = std::remove_cv_t<std::remove_reference_t<
            decltype(rate)>>;

        using Component = FlowView::Component;
        using Direction = FlowView::Direction;

        BOOST_CHECK_CLOSE(rate.flow(Component::Oil), -0.1f, 1.0e-5);
        BOOST_CHECK_CLOSE(rate.flow(Component::Gas), -0.2f, 1.0e-5);
        BOOST_CHECK_CLOSE(rate.flow(Component::Water), -0.3f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Disgas), -0.4f, 1.0e-5);
        BOOST_CHECK_CLOSE(rate.flow(Component::Vapoil), -0.5f, 1.0e-6);

        BOOST_CHECK_CLOSE(rate.flow(Component::Oil, Direction::Positive), 0.0f, 1.0e-5);
        BOOST_CHECK_CLOSE(rate.flow(Component::Gas, Direction::Positive), 0.0f, 1.0e-5);
        BOOST_CHECK_CLOSE(rate.flow(Component::Water, Direction::Positive), 0.0f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Disgas, Direction::Positive), 0.0f, 1.0e-5);
        BOOST_CHECK_CLOSE(rate.flow(Component::Vapoil, Direction::Positive), 0.0f, 1.0e-6);

        BOOST_CHECK_CLOSE(rate.flow(Component::Oil, Direction::Negative), -0.1f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Gas, Direction::Negative), -0.2f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Water, Direction::Negative), -0.3f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Disgas, Direction::Negative), -0.4f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Vapoil, Direction::Negative), -0.5f, 1.0e-6);
    }
}

BOOST_AUTO_TEST_SUITE_END() // Four_Processes

BOOST_AUTO_TEST_SUITE_END() // MultiRank_ReadWriteMerge

// =====================================================================

BOOST_AUTO_TEST_SUITE(MultiArray_Wrapper)

namespace {
    auto makeMessageBuffer()
    {
        return Dune::Point2PointCommunicator<Dune::SimpleMessageBuffer>
            ::MessageBufferType{};
    }

    template <typename Predicate>
    void addConnections(Predicate&&           isInterior,
                        Opm::InterRegFlowMap& rank)
    {
        rank.addConnection({ 0, 0, isInterior(0) }, { 1, 1, isInterior(1) }, conn_01());
        rank.addConnection({ 0, 0, isInterior(0) }, { 2, 2, isInterior(2) }, conn_02());
        rank.addConnection({ 1, 1, isInterior(1) }, { 3, 3, isInterior(3) }, conn_13());
        rank.addConnection({ 2, 2, isInterior(2) }, { 3, 3, isInterior(3) }, conn_23());
    }
}

BOOST_AUTO_TEST_CASE(Single_Process)
{
    const auto fipnum = all_same_region();
    const auto fipspl = left_right_split_region();
    const auto fipchk = checker_board_region();
    const auto fipsep = all_separate_region();

    auto flows = Opm::InterRegFlowMap {
        fipnum.size(),
        {
            { "FIPNUM", std::cref(fipnum) },
            { "FIPSPL", std::cref(fipspl) },
            { "FIPCHK", std::cref(fipchk) },
            { "FIPSEP", std::cref(fipsep) },
        }
    };

    {
        const auto names = flows.names();
        const auto expect = std::vector<std::string> {
            "FIPNUM", "FIPSPL", "FIPCHK", "FIPSEP",
        };

        BOOST_CHECK_MESSAGE(names == expect,
                            "Region array names don't match expected");
    }

    flows.addConnection({ 0, true }, { 1, true }, conn_01());
    flows.addConnection({ 0, true }, { 2, true }, conn_02());
    flows.addConnection({ 1, true }, { 3, true }, conn_13());
    flows.addConnection({ 2, true }, { 3, true }, conn_23());

    {
        const auto maxLocalRegID = flows.getLocalMaxRegionID();
        const auto expect = std::vector<std::size_t> {
            1, 2, 2, 4,
        };

        BOOST_CHECK_MESSAGE(maxLocalRegID == expect,
                            "Maximum local region IDs must match expected");

        BOOST_CHECK_MESSAGE(! flows.assignGlobalMaxRegionID({ 5, 1, 4, 5, }),
                            "Assigning small global maximum "
                            "region IDs must NOT succeed");

        BOOST_CHECK_MESSAGE(flows.assignGlobalMaxRegionID({ 5, 4, 4, 5, }),
                            "Assigning global maximum region IDs must succeed");
    }

    flows.compress();

    const auto iregFlows = flows.getInterRegFlows();
    BOOST_REQUIRE_EQUAL(iregFlows.size(), std::size_t{4});

    // FIPNUM
    {
        const auto& map = iregFlows[0];
        BOOST_CHECK_EQUAL(map.numRegions(), 5);

        for (auto r1 = 0*map.numRegions(); r1 < map.numRegions(); ++r1) {
            for (auto r2 = r1 + 1; r2 < map.numRegions(); ++r2) {
                auto flow = map.getInterRegFlows(r1, r2);
                BOOST_CHECK_MESSAGE(! flow.has_value(),
                                    "There must not be inter-regional flow "
                                    "between regions " << r1 << " and " << r2);
            }
        }
    }

    // FIPSPL
    {
        const auto& map = iregFlows[1];
        BOOST_CHECK_EQUAL(map.numRegions(), 4);

        {
            const auto q12 = map.getInterRegFlows(0, 1);
            BOOST_REQUIRE_MESSAGE(q12.has_value(),
                                  "There must be inter-region flows "
                                  "between regions 1 and 2");

            const auto& [rate, sign] = q12.value();

            BOOST_CHECK_EQUAL(sign, 1.0);

            using FlowView = std::remove_cv_t<std::remove_reference_t<
                decltype(rate)>>;

            using Component = FlowView::Component;
            using Direction = FlowView::Direction;

            BOOST_CHECK_CLOSE(rate.flow(Component::Oil), 0.88f, 1.0e-6);
            BOOST_CHECK_CLOSE(rate.flow(Component::Gas), 1.76f, 1.0e-6);
            BOOST_CHECK_CLOSE(rate.flow(Component::Water), 2.64f, 1.0e-5);
            BOOST_CHECK_CLOSE(rate.flow(Component::Disgas), 3.52f, 1.0e-6);
            BOOST_CHECK_CLOSE(rate.flow(Component::Vapoil), 3.75f, 1.0e-6);

            BOOST_CHECK_CLOSE(rate.flow(Component::Oil, Direction::Positive), 1.0f, 1.0e-6);
            BOOST_CHECK_CLOSE(rate.flow(Component::Gas, Direction::Positive), 2.0f, 1.0e-6);
            BOOST_CHECK_CLOSE(rate.flow(Component::Water, Direction::Positive), 3.0f, 1.0e-6);
            BOOST_CHECK_CLOSE(rate.flow(Component::Disgas, Direction::Positive), 4.0f, 1.0e-6);
            BOOST_CHECK_CLOSE(rate.flow(Component::Vapoil, Direction::Positive), 5.0f, 1.0e-6);

            BOOST_CHECK_CLOSE(rate.flow(Component::Oil, Direction::Negative), -0.12f, 1.0e-6);
            BOOST_CHECK_CLOSE(rate.flow(Component::Gas, Direction::Negative), -0.24f, 1.0e-6);
            BOOST_CHECK_CLOSE(rate.flow(Component::Water, Direction::Negative), -0.36f, 1.0e-6);
            BOOST_CHECK_CLOSE(rate.flow(Component::Disgas, Direction::Negative), -0.48f, 1.0e-6);
            BOOST_CHECK_CLOSE(rate.flow(Component::Vapoil, Direction::Negative), -1.25f, 1.0e-6);
        }

        {
            const auto q21 = map.getInterRegFlows(1, 0);
            BOOST_REQUIRE_MESSAGE(q21.has_value(),
                                  "There must be inter-region flows "
                                  "between regions 2 and 1");

            const auto& [rate, sign] = q21.value();

            BOOST_CHECK_EQUAL(sign, -1.0);

            using FlowView = std::remove_cv_t<std::remove_reference_t<
                decltype(rate)>>;

            using Component = FlowView::Component;
            using Direction = FlowView::Direction;

            BOOST_CHECK_CLOSE(rate.flow(Component::Oil), 0.88f, 1.0e-6);
            BOOST_CHECK_CLOSE(rate.flow(Component::Gas), 1.76f, 1.0e-6);
            BOOST_CHECK_CLOSE(rate.flow(Component::Water), 2.64f, 1.0e-5);
            BOOST_CHECK_CLOSE(rate.flow(Component::Disgas), 3.52f, 1.0e-6);
            BOOST_CHECK_CLOSE(rate.flow(Component::Vapoil), 3.75f, 1.0e-6);

            BOOST_CHECK_CLOSE(rate.flow(Component::Oil, Direction::Positive), 1.0f, 1.0e-6);
            BOOST_CHECK_CLOSE(rate.flow(Component::Gas, Direction::Positive), 2.0f, 1.0e-6);
            BOOST_CHECK_CLOSE(rate.flow(Component::Water, Direction::Positive), 3.0f, 1.0e-6);
            BOOST_CHECK_CLOSE(rate.flow(Component::Disgas, Direction::Positive), 4.0f, 1.0e-6);
            BOOST_CHECK_CLOSE(rate.flow(Component::Vapoil, Direction::Positive), 5.0f, 1.0e-6);

            BOOST_CHECK_CLOSE(rate.flow(Component::Oil, Direction::Negative), -0.12f, 1.0e-6);
            BOOST_CHECK_CLOSE(rate.flow(Component::Gas, Direction::Negative), -0.24f, 1.0e-6);
            BOOST_CHECK_CLOSE(rate.flow(Component::Water, Direction::Negative), -0.36f, 1.0e-6);
            BOOST_CHECK_CLOSE(rate.flow(Component::Disgas, Direction::Negative), -0.48f, 1.0e-6);
            BOOST_CHECK_CLOSE(rate.flow(Component::Vapoil, Direction::Negative), -1.25f, 1.0e-6);
        }
    }

    // FIPCHK
    {
        const auto& map = iregFlows[2];
        BOOST_CHECK_EQUAL(map.numRegions(), 4);

        const auto q12 = map.getInterRegFlows(0, 1);
        BOOST_REQUIRE_MESSAGE(q12.has_value(),
                              "There must be inter-region flows "
                              "between regions 1 and 2");

        const auto& [rate, sign] = q12.value();

        BOOST_CHECK_EQUAL(sign, 1.0);

        using FlowView = std::remove_cv_t<std::remove_reference_t<
            decltype(rate)>>;

        using Component = FlowView::Component;
        using Direction = FlowView::Direction;

        BOOST_CHECK_CLOSE(rate.flow(Component::Oil), 1.42f, 1.0e-5);
        BOOST_CHECK_CLOSE(rate.flow(Component::Gas), 2.84f, 1.0e-5);
        BOOST_CHECK_CLOSE(rate.flow(Component::Water), 4.26f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Disgas), 5.68f, 1.0e-5);
        BOOST_CHECK_CLOSE(rate.flow(Component::Vapoil), 7.75f, 1.0e-6);

        BOOST_CHECK_CLOSE(rate.flow(Component::Oil, Direction::Positive), 1.42f, 1.0e-5);
        BOOST_CHECK_CLOSE(rate.flow(Component::Gas, Direction::Positive), 2.84f, 1.0e-5);
        BOOST_CHECK_CLOSE(rate.flow(Component::Water, Direction::Positive), 4.26f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Disgas, Direction::Positive), 5.68f, 1.0e-5);
        BOOST_CHECK_CLOSE(rate.flow(Component::Vapoil, Direction::Positive), 7.75f, 1.0e-6);

        BOOST_CHECK_CLOSE(rate.flow(Component::Oil, Direction::Negative), 0.0f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Gas, Direction::Negative), 0.0f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Water, Direction::Negative), 0.0f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Disgas, Direction::Negative), 0.0f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Vapoil, Direction::Negative), 0.0f, 1.0e-6);
    }

    // FIPSEP
    {
        const auto& map = iregFlows[3];
        BOOST_CHECK_EQUAL(map.numRegions(), 5);

        {
            const auto q12 = map.getInterRegFlows(0, 1);
            BOOST_REQUIRE_MESSAGE(q12.has_value(),
                                  "There must be inter-region flows "
                                  "between regions 1 and 2");

            const auto& [rate, sign] = q12.value();

            BOOST_CHECK_EQUAL(sign, 1.0);

            using FlowView = std::remove_cv_t<std::remove_reference_t<
                decltype(rate)>>;

            using Component = FlowView::Component;
            using Direction = FlowView::Direction;

            BOOST_CHECK_CLOSE(rate.flow(Component::Oil), 0.12f, 1.0e-5);
            BOOST_CHECK_CLOSE(rate.flow(Component::Gas), 0.24f, 1.0e-5);
            BOOST_CHECK_CLOSE(rate.flow(Component::Water), 0.36f, 1.0e-6);
            BOOST_CHECK_CLOSE(rate.flow(Component::Disgas), 0.48f, 1.0e-5);
            BOOST_CHECK_CLOSE(rate.flow(Component::Vapoil), 1.25f, 1.0e-6);

            BOOST_CHECK_CLOSE(rate.flow(Component::Oil, Direction::Positive), 0.12f, 1.0e-5);
            BOOST_CHECK_CLOSE(rate.flow(Component::Gas, Direction::Positive), 0.24f, 1.0e-5);
            BOOST_CHECK_CLOSE(rate.flow(Component::Water, Direction::Positive), 0.36f, 1.0e-6);
            BOOST_CHECK_CLOSE(rate.flow(Component::Disgas, Direction::Positive), 0.48f, 1.0e-5);
            BOOST_CHECK_CLOSE(rate.flow(Component::Vapoil, Direction::Positive), 1.25f, 1.0e-6);

            BOOST_CHECK_CLOSE(rate.flow(Component::Oil, Direction::Negative), 0.0f, 1.0e-6);
            BOOST_CHECK_CLOSE(rate.flow(Component::Gas, Direction::Negative), 0.0f, 1.0e-6);
            BOOST_CHECK_CLOSE(rate.flow(Component::Water, Direction::Negative), 0.0f, 1.0e-6);
            BOOST_CHECK_CLOSE(rate.flow(Component::Disgas, Direction::Negative), 0.0f, 1.0e-6);
            BOOST_CHECK_CLOSE(rate.flow(Component::Vapoil, Direction::Negative), 0.0f, 1.0e-6);
        }

        {
            const auto q13 = map.getInterRegFlows(0, 2);
            BOOST_REQUIRE_MESSAGE(q13.has_value(),
                                  "There must be inter-region flows "
                                  "between regions 1 and 3");

            const auto& [rate, sign] = q13.value();

            BOOST_CHECK_EQUAL(sign, 1.0);

            using FlowView = std::remove_cv_t<std::remove_reference_t<
                decltype(rate)>>;

            using Component = FlowView::Component;
            using Direction = FlowView::Direction;

            BOOST_CHECK_CLOSE(rate.flow(Component::Oil), 0.2f, 1.0e-5);
            BOOST_CHECK_CLOSE(rate.flow(Component::Gas), 0.4f, 1.0e-5);
            BOOST_CHECK_CLOSE(rate.flow(Component::Water), 0.6f, 1.0e-6);
            BOOST_CHECK_CLOSE(rate.flow(Component::Disgas), 0.8f, 1.0e-5);
            BOOST_CHECK_CLOSE(rate.flow(Component::Vapoil), 1.0f, 1.0e-6);

            BOOST_CHECK_CLOSE(rate.flow(Component::Oil, Direction::Positive), 0.2f, 1.0e-5);
            BOOST_CHECK_CLOSE(rate.flow(Component::Gas, Direction::Positive), 0.4f, 1.0e-5);
            BOOST_CHECK_CLOSE(rate.flow(Component::Water, Direction::Positive), 0.6f, 1.0e-6);
            BOOST_CHECK_CLOSE(rate.flow(Component::Disgas, Direction::Positive), 0.8f, 1.0e-5);
            BOOST_CHECK_CLOSE(rate.flow(Component::Vapoil, Direction::Positive), 1.0f, 1.0e-6);

            BOOST_CHECK_CLOSE(rate.flow(Component::Oil, Direction::Negative), 0.0f, 1.0e-6);
            BOOST_CHECK_CLOSE(rate.flow(Component::Gas, Direction::Negative), 0.0f, 1.0e-6);
            BOOST_CHECK_CLOSE(rate.flow(Component::Water, Direction::Negative), 0.0f, 1.0e-6);
            BOOST_CHECK_CLOSE(rate.flow(Component::Disgas, Direction::Negative), 0.0f, 1.0e-6);
            BOOST_CHECK_CLOSE(rate.flow(Component::Vapoil, Direction::Negative), 0.0f, 1.0e-6);
        }

        {
            const auto q14 = map.getInterRegFlows(0, 3);
            BOOST_CHECK_MESSAGE(! q14.has_value(),
                                "There must NOT be inter-region "
                                "flows between regions 1 and 4");
        }

        {
            const auto q23 = map.getInterRegFlows(1, 2);
            BOOST_CHECK_MESSAGE(! q23.has_value(),
                                "There must NOT be inter-region "
                                "flows between regions 2 and 3");
        }

        {
            const auto q24 = map.getInterRegFlows(1, 3);
            BOOST_REQUIRE_MESSAGE(q24.has_value(),
                                  "There must be inter-region flows "
                                  "between regions 2 and 4");

            const auto& [rate, sign] = q24.value();

            BOOST_CHECK_EQUAL(sign, 1.0);

            using FlowView = std::remove_cv_t<std::remove_reference_t<
                decltype(rate)>>;

            using Component = FlowView::Component;
            using Direction = FlowView::Direction;

            BOOST_CHECK_CLOSE(rate.flow(Component::Oil), -0.1f, 1.0e-5);
            BOOST_CHECK_CLOSE(rate.flow(Component::Gas), -0.2f, 1.0e-5);
            BOOST_CHECK_CLOSE(rate.flow(Component::Water), -0.3f, 1.0e-6);
            BOOST_CHECK_CLOSE(rate.flow(Component::Disgas), -0.4f, 1.0e-5);
            BOOST_CHECK_CLOSE(rate.flow(Component::Vapoil), -0.5f, 1.0e-6);

            BOOST_CHECK_CLOSE(rate.flow(Component::Oil, Direction::Positive), 0.0f, 1.0e-5);
            BOOST_CHECK_CLOSE(rate.flow(Component::Gas, Direction::Positive), 0.0f, 1.0e-5);
            BOOST_CHECK_CLOSE(rate.flow(Component::Water, Direction::Positive), 0.0f, 1.0e-6);
            BOOST_CHECK_CLOSE(rate.flow(Component::Disgas, Direction::Positive), 0.0f, 1.0e-5);
            BOOST_CHECK_CLOSE(rate.flow(Component::Vapoil, Direction::Positive), 0.0f, 1.0e-6);

            BOOST_CHECK_CLOSE(rate.flow(Component::Oil, Direction::Negative), -0.1f, 1.0e-6);
            BOOST_CHECK_CLOSE(rate.flow(Component::Gas, Direction::Negative), -0.2f, 1.0e-6);
            BOOST_CHECK_CLOSE(rate.flow(Component::Water, Direction::Negative), -0.3f, 1.0e-6);
            BOOST_CHECK_CLOSE(rate.flow(Component::Disgas, Direction::Negative), -0.4f, 1.0e-6);
            BOOST_CHECK_CLOSE(rate.flow(Component::Vapoil, Direction::Negative), -0.5f, 1.0e-6);
        }
    }
}

BOOST_AUTO_TEST_CASE(Two_Processes)
{
    using Map = Opm::InterRegFlowMap;

    const auto fipnum = all_same_region();
    const auto fipspl = left_right_split_region();
    const auto fipchk = checker_board_region();
    const auto fipsep = all_separate_region();

    auto rank = std::vector<Map>(2, Map {
            fipnum.size(), {
                { "FIPNUM", std::cref(fipnum) },
                { "FIPSPL", std::cref(fipspl) },
                { "FIPCHK", std::cref(fipchk) },
                { "FIPSEP", std::cref(fipsep) },
            }});

    {
        const auto expect = std::vector<std::size_t>{ 1, 2, 2, 4 };

        BOOST_CHECK_MESSAGE(rank[0].getLocalMaxRegionID() == expect,
                            "Local maximum region IDs must match expected on rank 0");

        BOOST_CHECK_MESSAGE(rank[1].getLocalMaxRegionID() == expect,
                            "Local maximum region IDs must match expected on rank 1");
    }

    addConnections(TwoProc::P1::isInterior, rank[0]);
    rank[0].assignGlobalMaxRegionID({3, 2, 2, 4});

    addConnections(TwoProc::P2::isInterior, rank[1]);
    rank[1].assignGlobalMaxRegionID({5, 2, 2, 4});

    for (auto& r : rank) {
        r.compress();
    }

    auto buffer = makeMessageBuffer();

    rank[1].write(buffer);
    rank[0].read(buffer);

    BOOST_CHECK_MESSAGE(rank[0].readIsConsistent(),
                        "Consistent write() must yield consistent read()");

    BOOST_CHECK_THROW(rank[0].addConnection({ 0, TwoProc::P1::isInterior(0) },
                                            { 1, TwoProc::P1::isInterior(1) }, conn_01()),
                      std::logic_error);

    rank[0].compress();

    {
        const auto names = rank[0].names();
        const auto expect = std::vector<std::string> {
            "FIPNUM", "FIPSPL", "FIPCHK", "FIPSEP",
        };

        BOOST_CHECK_MESSAGE(names == expect,
                            "Region array names don't match expected");
    }

    const auto iregFlows = rank[0].getInterRegFlows();
    BOOST_REQUIRE_EQUAL(iregFlows.size(), std::size_t{4});

    // FIPNUM
    {
        const auto& map = iregFlows[0];
        BOOST_CHECK_EQUAL(map.numRegions(), 5);

        for (auto r1 = 0*map.numRegions(); r1 < map.numRegions(); ++r1) {
            for (auto r2 = r1 + 1; r2 < map.numRegions(); ++r2) {
                auto flow = map.getInterRegFlows(r1, r2);
                BOOST_CHECK_MESSAGE(! flow.has_value(),
                                    "There must not be inter-regional flow "
                                    "between regions " << r1 << " and " << r2);
            }
        }
    }

    // FIPSPL
    {
        const auto& map = iregFlows[1];
        BOOST_CHECK_EQUAL(map.numRegions(), 2);

        {
            const auto q12 = map.getInterRegFlows(0, 1);
            BOOST_REQUIRE_MESSAGE(q12.has_value(),
                                  "There must be inter-region flows "
                                  "between regions 1 and 2");

            const auto& [rate, sign] = q12.value();

            BOOST_CHECK_EQUAL(sign, 1.0);

            using FlowView = std::remove_cv_t<std::remove_reference_t<
                decltype(rate)>>;

            using Component = FlowView::Component;
            using Direction = FlowView::Direction;

            BOOST_CHECK_CLOSE(rate.flow(Component::Oil), 0.88f, 1.0e-6);
            BOOST_CHECK_CLOSE(rate.flow(Component::Gas), 1.76f, 1.0e-6);
            BOOST_CHECK_CLOSE(rate.flow(Component::Water), 2.64f, 1.0e-5);
            BOOST_CHECK_CLOSE(rate.flow(Component::Disgas), 3.52f, 1.0e-6);
            BOOST_CHECK_CLOSE(rate.flow(Component::Vapoil), 3.75f, 1.0e-6);

            BOOST_CHECK_CLOSE(rate.flow(Component::Oil, Direction::Positive), 1.0f, 1.0e-6);
            BOOST_CHECK_CLOSE(rate.flow(Component::Gas, Direction::Positive), 2.0f, 1.0e-6);
            BOOST_CHECK_CLOSE(rate.flow(Component::Water, Direction::Positive), 3.0f, 1.0e-6);
            BOOST_CHECK_CLOSE(rate.flow(Component::Disgas, Direction::Positive), 4.0f, 1.0e-6);
            BOOST_CHECK_CLOSE(rate.flow(Component::Vapoil, Direction::Positive), 5.0f, 1.0e-6);

            BOOST_CHECK_CLOSE(rate.flow(Component::Oil, Direction::Negative), -0.12f, 1.0e-6);
            BOOST_CHECK_CLOSE(rate.flow(Component::Gas, Direction::Negative), -0.24f, 1.0e-6);
            BOOST_CHECK_CLOSE(rate.flow(Component::Water, Direction::Negative), -0.36f, 1.0e-6);
            BOOST_CHECK_CLOSE(rate.flow(Component::Disgas, Direction::Negative), -0.48f, 1.0e-6);
            BOOST_CHECK_CLOSE(rate.flow(Component::Vapoil, Direction::Negative), -1.25f, 1.0e-6);
        }

        {
            const auto q21 = map.getInterRegFlows(1, 0);
            BOOST_REQUIRE_MESSAGE(q21.has_value(),
                                  "There must be inter-region flows "
                                  "between regions 2 and 1");

            const auto& [rate, sign] = q21.value();

            BOOST_CHECK_EQUAL(sign, -1.0);

            using FlowView = std::remove_cv_t<std::remove_reference_t<
                decltype(rate)>>;

            using Component = FlowView::Component;
            using Direction = FlowView::Direction;

            BOOST_CHECK_CLOSE(rate.flow(Component::Oil), 0.88f, 1.0e-6);
            BOOST_CHECK_CLOSE(rate.flow(Component::Gas), 1.76f, 1.0e-6);
            BOOST_CHECK_CLOSE(rate.flow(Component::Water), 2.64f, 1.0e-5);
            BOOST_CHECK_CLOSE(rate.flow(Component::Disgas), 3.52f, 1.0e-6);
            BOOST_CHECK_CLOSE(rate.flow(Component::Vapoil), 3.75f, 1.0e-6);

            BOOST_CHECK_CLOSE(rate.flow(Component::Oil, Direction::Positive), 1.0f, 1.0e-6);
            BOOST_CHECK_CLOSE(rate.flow(Component::Gas, Direction::Positive), 2.0f, 1.0e-6);
            BOOST_CHECK_CLOSE(rate.flow(Component::Water, Direction::Positive), 3.0f, 1.0e-6);
            BOOST_CHECK_CLOSE(rate.flow(Component::Disgas, Direction::Positive), 4.0f, 1.0e-6);
            BOOST_CHECK_CLOSE(rate.flow(Component::Vapoil, Direction::Positive), 5.0f, 1.0e-6);

            BOOST_CHECK_CLOSE(rate.flow(Component::Oil, Direction::Negative), -0.12f, 1.0e-6);
            BOOST_CHECK_CLOSE(rate.flow(Component::Gas, Direction::Negative), -0.24f, 1.0e-6);
            BOOST_CHECK_CLOSE(rate.flow(Component::Water, Direction::Negative), -0.36f, 1.0e-6);
            BOOST_CHECK_CLOSE(rate.flow(Component::Disgas, Direction::Negative), -0.48f, 1.0e-6);
            BOOST_CHECK_CLOSE(rate.flow(Component::Vapoil, Direction::Negative), -1.25f, 1.0e-6);
        }
    }

    // FIPCHK
    {
        const auto& map = iregFlows[2];
        BOOST_CHECK_EQUAL(map.numRegions(), 2);

        const auto q12 = map.getInterRegFlows(0, 1);
        BOOST_REQUIRE_MESSAGE(q12.has_value(),
                              "There must be inter-region flows "
                              "between regions 1 and 2");

        const auto& [rate, sign] = q12.value();

        BOOST_CHECK_EQUAL(sign, 1.0);

        using FlowView = std::remove_cv_t<std::remove_reference_t<
            decltype(rate)>>;

        using Component = FlowView::Component;
        using Direction = FlowView::Direction;

        BOOST_CHECK_CLOSE(rate.flow(Component::Oil), 1.42f, 1.0e-5);
        BOOST_CHECK_CLOSE(rate.flow(Component::Gas), 2.84f, 1.0e-5);
        BOOST_CHECK_CLOSE(rate.flow(Component::Water), 4.26f, 2.0e-5);
        BOOST_CHECK_CLOSE(rate.flow(Component::Disgas), 5.68f, 1.0e-5);
        BOOST_CHECK_CLOSE(rate.flow(Component::Vapoil), 7.75f, 1.0e-6);

        BOOST_CHECK_CLOSE(rate.flow(Component::Oil, Direction::Positive), 1.42f, 1.0e-5);
        BOOST_CHECK_CLOSE(rate.flow(Component::Gas, Direction::Positive), 2.84f, 1.0e-5);
        BOOST_CHECK_CLOSE(rate.flow(Component::Water, Direction::Positive), 4.26f, 2.0e-5);
        BOOST_CHECK_CLOSE(rate.flow(Component::Disgas, Direction::Positive), 5.68f, 1.0e-5);
        BOOST_CHECK_CLOSE(rate.flow(Component::Vapoil, Direction::Positive), 7.75f, 1.0e-6);

        BOOST_CHECK_CLOSE(rate.flow(Component::Oil, Direction::Negative), 0.0f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Gas, Direction::Negative), 0.0f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Water, Direction::Negative), 0.0f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Disgas, Direction::Negative), 0.0f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Vapoil, Direction::Negative), 0.0f, 1.0e-6);
    }

    // FIPSEP
    {
        const auto& map = iregFlows[3];
        BOOST_CHECK_EQUAL(map.numRegions(), 4);

        {
            const auto q12 = map.getInterRegFlows(0, 1);
            BOOST_REQUIRE_MESSAGE(q12.has_value(),
                                  "There must be inter-region flows "
                                  "between regions 1 and 2");

            const auto& [rate, sign] = q12.value();

            BOOST_CHECK_EQUAL(sign, 1.0);

            using FlowView = std::remove_cv_t<std::remove_reference_t<
                decltype(rate)>>;

            using Component = FlowView::Component;
            using Direction = FlowView::Direction;

            BOOST_CHECK_CLOSE(rate.flow(Component::Oil), 0.12f, 1.0e-5);
            BOOST_CHECK_CLOSE(rate.flow(Component::Gas), 0.24f, 1.0e-5);
            BOOST_CHECK_CLOSE(rate.flow(Component::Water), 0.36f, 1.0e-6);
            BOOST_CHECK_CLOSE(rate.flow(Component::Disgas), 0.48f, 1.0e-5);
            BOOST_CHECK_CLOSE(rate.flow(Component::Vapoil), 1.25f, 1.0e-6);

            BOOST_CHECK_CLOSE(rate.flow(Component::Oil, Direction::Positive), 0.12f, 1.0e-5);
            BOOST_CHECK_CLOSE(rate.flow(Component::Gas, Direction::Positive), 0.24f, 1.0e-5);
            BOOST_CHECK_CLOSE(rate.flow(Component::Water, Direction::Positive), 0.36f, 1.0e-6);
            BOOST_CHECK_CLOSE(rate.flow(Component::Disgas, Direction::Positive), 0.48f, 1.0e-5);
            BOOST_CHECK_CLOSE(rate.flow(Component::Vapoil, Direction::Positive), 1.25f, 1.0e-6);

            BOOST_CHECK_CLOSE(rate.flow(Component::Oil, Direction::Negative), 0.0f, 1.0e-6);
            BOOST_CHECK_CLOSE(rate.flow(Component::Gas, Direction::Negative), 0.0f, 1.0e-6);
            BOOST_CHECK_CLOSE(rate.flow(Component::Water, Direction::Negative), 0.0f, 1.0e-6);
            BOOST_CHECK_CLOSE(rate.flow(Component::Disgas, Direction::Negative), 0.0f, 1.0e-6);
            BOOST_CHECK_CLOSE(rate.flow(Component::Vapoil, Direction::Negative), 0.0f, 1.0e-6);
        }

        {
            const auto q13 = map.getInterRegFlows(0, 2);
            BOOST_REQUIRE_MESSAGE(q13.has_value(),
                                  "There must be inter-region flows "
                                  "between regions 1 and 3");

            const auto& [rate, sign] = q13.value();

            BOOST_CHECK_EQUAL(sign, 1.0);

            using FlowView = std::remove_cv_t<std::remove_reference_t<
                decltype(rate)>>;

            using Component = FlowView::Component;
            using Direction = FlowView::Direction;

            BOOST_CHECK_CLOSE(rate.flow(Component::Oil), 0.2f, 1.0e-5);
            BOOST_CHECK_CLOSE(rate.flow(Component::Gas), 0.4f, 1.0e-5);
            BOOST_CHECK_CLOSE(rate.flow(Component::Water), 0.6f, 1.0e-6);
            BOOST_CHECK_CLOSE(rate.flow(Component::Disgas), 0.8f, 1.0e-5);
            BOOST_CHECK_CLOSE(rate.flow(Component::Vapoil), 1.0f, 1.0e-6);

            BOOST_CHECK_CLOSE(rate.flow(Component::Oil, Direction::Positive), 0.2f, 1.0e-5);
            BOOST_CHECK_CLOSE(rate.flow(Component::Gas, Direction::Positive), 0.4f, 1.0e-5);
            BOOST_CHECK_CLOSE(rate.flow(Component::Water, Direction::Positive), 0.6f, 1.0e-6);
            BOOST_CHECK_CLOSE(rate.flow(Component::Disgas, Direction::Positive), 0.8f, 1.0e-5);
            BOOST_CHECK_CLOSE(rate.flow(Component::Vapoil, Direction::Positive), 1.0f, 1.0e-6);

            BOOST_CHECK_CLOSE(rate.flow(Component::Oil, Direction::Negative), 0.0f, 1.0e-6);
            BOOST_CHECK_CLOSE(rate.flow(Component::Gas, Direction::Negative), 0.0f, 1.0e-6);
            BOOST_CHECK_CLOSE(rate.flow(Component::Water, Direction::Negative), 0.0f, 1.0e-6);
            BOOST_CHECK_CLOSE(rate.flow(Component::Disgas, Direction::Negative), 0.0f, 1.0e-6);
            BOOST_CHECK_CLOSE(rate.flow(Component::Vapoil, Direction::Negative), 0.0f, 1.0e-6);
        }

        {
            const auto q14 = map.getInterRegFlows(0, 3);
            BOOST_CHECK_MESSAGE(! q14.has_value(),
                                "There must NOT be inter-region "
                                "flows between regions 1 and 4");
        }

        {
            const auto q23 = map.getInterRegFlows(1, 2);
            BOOST_CHECK_MESSAGE(! q23.has_value(),
                                "There must NOT be inter-region "
                                "flows between regions 2 and 3");
        }

        {
            const auto q24 = map.getInterRegFlows(1, 3);
            BOOST_REQUIRE_MESSAGE(q24.has_value(),
                                  "There must be inter-region flows "
                                  "between regions 2 and 4");

            const auto& [rate, sign] = q24.value();

            BOOST_CHECK_EQUAL(sign, 1.0);

            using FlowView = std::remove_cv_t<std::remove_reference_t<
                decltype(rate)>>;

            using Component = FlowView::Component;
            using Direction = FlowView::Direction;

            BOOST_CHECK_CLOSE(rate.flow(Component::Oil), -0.1f, 1.0e-5);
            BOOST_CHECK_CLOSE(rate.flow(Component::Gas), -0.2f, 1.0e-5);
            BOOST_CHECK_CLOSE(rate.flow(Component::Water), -0.3f, 1.0e-6);
            BOOST_CHECK_CLOSE(rate.flow(Component::Disgas), -0.4f, 1.0e-5);
            BOOST_CHECK_CLOSE(rate.flow(Component::Vapoil), -0.5f, 1.0e-6);

            BOOST_CHECK_CLOSE(rate.flow(Component::Oil, Direction::Positive), 0.0f, 1.0e-5);
            BOOST_CHECK_CLOSE(rate.flow(Component::Gas, Direction::Positive), 0.0f, 1.0e-5);
            BOOST_CHECK_CLOSE(rate.flow(Component::Water, Direction::Positive), 0.0f, 1.0e-6);
            BOOST_CHECK_CLOSE(rate.flow(Component::Disgas, Direction::Positive), 0.0f, 1.0e-5);
            BOOST_CHECK_CLOSE(rate.flow(Component::Vapoil, Direction::Positive), 0.0f, 1.0e-6);

            BOOST_CHECK_CLOSE(rate.flow(Component::Oil, Direction::Negative), -0.1f, 1.0e-6);
            BOOST_CHECK_CLOSE(rate.flow(Component::Gas, Direction::Negative), -0.2f, 1.0e-6);
            BOOST_CHECK_CLOSE(rate.flow(Component::Water, Direction::Negative), -0.3f, 1.0e-6);
            BOOST_CHECK_CLOSE(rate.flow(Component::Disgas, Direction::Negative), -0.4f, 1.0e-6);
            BOOST_CHECK_CLOSE(rate.flow(Component::Vapoil, Direction::Negative), -0.5f, 1.0e-6);
        }
    }
}

BOOST_AUTO_TEST_CASE(Two_Processes_Inconsistent_ReadWrite)
{
    using Map = Opm::InterRegFlowMap;

    const auto fipnum = all_same_region();
    const auto fipspl = left_right_split_region();
    const auto fipchk = checker_board_region();
    const auto fipsep = all_separate_region();

    auto rank = std::vector<Map>{
        Map {
            fipnum.size(), {
                { "FIPNUM", std::cref(fipnum) },
                { "FIPSPL", std::cref(fipspl) },
                { "FIPCHK", std::cref(fipchk) },
                { "FIPSEP", std::cref(fipsep) },
            }
        },

        Map {
            fipnum.size(), {
                { "FIPSPL", std::cref(fipspl) },
                { "FIPSEP", std::cref(fipsep) },
            }
        },
    };

    addConnections(TwoProc::P1::isInterior, rank[0]);
    rank[0].assignGlobalMaxRegionID({3, 2, 2, 4});

    addConnections(TwoProc::P2::isInterior, rank[1]);
    rank[1].assignGlobalMaxRegionID({2, 4});
    for (auto& r : rank) {
        r.compress();
    }

    auto buffer = makeMessageBuffer();

    rank[1].write(buffer);
    rank[0].read(buffer);

    BOOST_CHECK_MESSAGE(! rank[0].readIsConsistent(),
                        "Inconsistent write() must yield inconsistent read()");
}

BOOST_AUTO_TEST_CASE(Four_Processes)
{
    using Map = Opm::InterRegFlowMap;

    const auto fipnum = all_same_region();
    const auto fipspl = left_right_split_region();
    const auto fipchk = checker_board_region();
    const auto fipsep = all_separate_region();

    auto rank = std::vector<Map>(4, Map {
            fipnum.size(), {
                { "FIPNUM", std::cref(fipnum) },
                { "FIPSPL", std::cref(fipspl) },
                { "FIPCHK", std::cref(fipchk) },
                { "FIPSEP", std::cref(fipsep) },
            }});
    {
        const auto expect = std::vector<std::size_t>{ 1, 2, 2, 4 };

        BOOST_CHECK_MESSAGE(rank[0].getLocalMaxRegionID() == expect,
                            "Local maximum region IDs must match expected on rank 0");

        BOOST_CHECK_MESSAGE(rank[1].getLocalMaxRegionID() == expect,
                            "Local maximum region IDs must match expected on rank 1");
    }

    addConnections(FourProc::P1::isInterior, rank[0]);
    rank[0].assignGlobalMaxRegionID({3, 2, 2, 4});

    addConnections(FourProc::P2::isInterior, rank[1]);
    addConnections(FourProc::P3::isInterior, rank[2]);

    addConnections(FourProc::P4::isInterior, rank[3]);
    rank[3].assignGlobalMaxRegionID({5, 2, 2, 4});

    for (auto& r : rank) {
        r.compress();
    }

    auto buffer = makeMessageBuffer();

    rank[1].write(buffer);
    rank[0].write(buffer);
    rank[3].write(buffer);
    rank[2].read(buffer);
    rank[2].read(buffer);
    rank[2].read(buffer);

    BOOST_CHECK_MESSAGE(rank[2].readIsConsistent(),
                        "Consistent write() must yield consistent read()");

    BOOST_CHECK_THROW(rank[2].addConnection({ 0, TwoProc::P1::isInterior(0) },
                                            { 1, TwoProc::P1::isInterior(1) }, conn_01()),
                      std::logic_error);

    rank[2].compress();

    {
        const auto names = rank[0].names();
        const auto expect = std::vector<std::string> {
            "FIPNUM", "FIPSPL", "FIPCHK", "FIPSEP",
        };

        BOOST_CHECK_MESSAGE(names == expect,
                            "Region array names don't match expected");
    }

    const auto iregFlows = rank[2].getInterRegFlows();
    BOOST_REQUIRE_EQUAL(iregFlows.size(), std::size_t{4});

    // FIPNUM
    {
        const auto& map = iregFlows[0];
        BOOST_CHECK_EQUAL(map.numRegions(), 5);

        for (auto r1 = 0*map.numRegions(); r1 < map.numRegions(); ++r1) {
            for (auto r2 = r1 + 1; r2 < map.numRegions(); ++r2) {
                auto flow = map.getInterRegFlows(r1, r2);
                BOOST_CHECK_MESSAGE(! flow.has_value(),
                                    "There must not be inter-regional flow "
                                    "between regions " << r1 << " and " << r2);
            }
        }
    }

    // FIPSPL
    {
        const auto& map = iregFlows[1];
        BOOST_CHECK_EQUAL(map.numRegions(), 2);

        {
            const auto q12 = map.getInterRegFlows(0, 1);
            BOOST_REQUIRE_MESSAGE(q12.has_value(),
                                  "There must be inter-region flows "
                                  "between regions 1 and 2");

            const auto& [rate, sign] = q12.value();

            BOOST_CHECK_EQUAL(sign, 1.0);

            using FlowView = std::remove_cv_t<std::remove_reference_t<
                decltype(rate)>>;

            using Component = FlowView::Component;
            using Direction = FlowView::Direction;

            BOOST_CHECK_CLOSE(rate.flow(Component::Oil), 0.88f, 1.0e-6);
            BOOST_CHECK_CLOSE(rate.flow(Component::Gas), 1.76f, 1.0e-6);
            BOOST_CHECK_CLOSE(rate.flow(Component::Water), 2.64f, 1.0e-5);
            BOOST_CHECK_CLOSE(rate.flow(Component::Disgas), 3.52f, 1.0e-6);
            BOOST_CHECK_CLOSE(rate.flow(Component::Vapoil), 3.75f, 1.0e-6);

            BOOST_CHECK_CLOSE(rate.flow(Component::Oil, Direction::Positive), 1.0f, 1.0e-6);
            BOOST_CHECK_CLOSE(rate.flow(Component::Gas, Direction::Positive), 2.0f, 1.0e-6);
            BOOST_CHECK_CLOSE(rate.flow(Component::Water, Direction::Positive), 3.0f, 1.0e-6);
            BOOST_CHECK_CLOSE(rate.flow(Component::Disgas, Direction::Positive), 4.0f, 1.0e-6);
            BOOST_CHECK_CLOSE(rate.flow(Component::Vapoil, Direction::Positive), 5.0f, 1.0e-6);

            BOOST_CHECK_CLOSE(rate.flow(Component::Oil, Direction::Negative), -0.12f, 1.0e-6);
            BOOST_CHECK_CLOSE(rate.flow(Component::Gas, Direction::Negative), -0.24f, 1.0e-6);
            BOOST_CHECK_CLOSE(rate.flow(Component::Water, Direction::Negative), -0.36f, 1.0e-6);
            BOOST_CHECK_CLOSE(rate.flow(Component::Disgas, Direction::Negative), -0.48f, 1.0e-6);
            BOOST_CHECK_CLOSE(rate.flow(Component::Vapoil, Direction::Negative), -1.25f, 1.0e-6);
        }

        {
            const auto q21 = map.getInterRegFlows(1, 0);
            BOOST_REQUIRE_MESSAGE(q21.has_value(),
                                  "There must be inter-region flows "
                                  "between regions 2 and 1");

            const auto& [rate, sign] = q21.value();

            BOOST_CHECK_EQUAL(sign, -1.0);

            using FlowView = std::remove_cv_t<std::remove_reference_t<
                decltype(rate)>>;

            using Component = FlowView::Component;
            using Direction = FlowView::Direction;

            BOOST_CHECK_CLOSE(rate.flow(Component::Oil), 0.88f, 1.0e-6);
            BOOST_CHECK_CLOSE(rate.flow(Component::Gas), 1.76f, 1.0e-6);
            BOOST_CHECK_CLOSE(rate.flow(Component::Water), 2.64f, 1.0e-5);
            BOOST_CHECK_CLOSE(rate.flow(Component::Disgas), 3.52f, 1.0e-6);
            BOOST_CHECK_CLOSE(rate.flow(Component::Vapoil), 3.75f, 1.0e-6);

            BOOST_CHECK_CLOSE(rate.flow(Component::Oil, Direction::Positive), 1.0f, 1.0e-6);
            BOOST_CHECK_CLOSE(rate.flow(Component::Gas, Direction::Positive), 2.0f, 1.0e-6);
            BOOST_CHECK_CLOSE(rate.flow(Component::Water, Direction::Positive), 3.0f, 1.0e-6);
            BOOST_CHECK_CLOSE(rate.flow(Component::Disgas, Direction::Positive), 4.0f, 1.0e-6);
            BOOST_CHECK_CLOSE(rate.flow(Component::Vapoil, Direction::Positive), 5.0f, 1.0e-6);

            BOOST_CHECK_CLOSE(rate.flow(Component::Oil, Direction::Negative), -0.12f, 1.0e-6);
            BOOST_CHECK_CLOSE(rate.flow(Component::Gas, Direction::Negative), -0.24f, 1.0e-6);
            BOOST_CHECK_CLOSE(rate.flow(Component::Water, Direction::Negative), -0.36f, 1.0e-6);
            BOOST_CHECK_CLOSE(rate.flow(Component::Disgas, Direction::Negative), -0.48f, 1.0e-6);
            BOOST_CHECK_CLOSE(rate.flow(Component::Vapoil, Direction::Negative), -1.25f, 1.0e-6);
        }
    }

    // FIPCHK
    {
        const auto& map = iregFlows[2];
        BOOST_CHECK_EQUAL(map.numRegions(), 2);

        const auto q12 = map.getInterRegFlows(0, 1);
        BOOST_REQUIRE_MESSAGE(q12.has_value(),
                              "There must be inter-region flows "
                              "between regions 1 and 2");

        const auto& [rate, sign] = q12.value();

        BOOST_CHECK_EQUAL(sign, 1.0);

        using FlowView = std::remove_cv_t<std::remove_reference_t<
            decltype(rate)>>;

        using Component = FlowView::Component;
        using Direction = FlowView::Direction;

        BOOST_CHECK_CLOSE(rate.flow(Component::Oil), 1.42f, 1.0e-5);
        BOOST_CHECK_CLOSE(rate.flow(Component::Gas), 2.84f, 1.0e-5);
        BOOST_CHECK_CLOSE(rate.flow(Component::Water), 4.26f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Disgas), 5.68f, 1.0e-5);
        BOOST_CHECK_CLOSE(rate.flow(Component::Vapoil), 7.75f, 1.0e-6);

        BOOST_CHECK_CLOSE(rate.flow(Component::Oil, Direction::Positive), 1.42f, 1.0e-5);
        BOOST_CHECK_CLOSE(rate.flow(Component::Gas, Direction::Positive), 2.84f, 1.0e-5);
        BOOST_CHECK_CLOSE(rate.flow(Component::Water, Direction::Positive), 4.26f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Disgas, Direction::Positive), 5.68f, 1.0e-5);
        BOOST_CHECK_CLOSE(rate.flow(Component::Vapoil, Direction::Positive), 7.75f, 1.0e-6);

        BOOST_CHECK_CLOSE(rate.flow(Component::Oil, Direction::Negative), 0.0f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Gas, Direction::Negative), 0.0f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Water, Direction::Negative), 0.0f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Disgas, Direction::Negative), 0.0f, 1.0e-6);
        BOOST_CHECK_CLOSE(rate.flow(Component::Vapoil, Direction::Negative), 0.0f, 1.0e-6);
    }

    // FIPSEP
    {
        const auto& map = iregFlows[3];
        BOOST_CHECK_EQUAL(map.numRegions(), 4);

        {
            const auto q12 = map.getInterRegFlows(0, 1);
            BOOST_REQUIRE_MESSAGE(q12.has_value(),
                                  "There must be inter-region flows "
                                  "between regions 1 and 2");

            const auto& [rate, sign] = q12.value();

            BOOST_CHECK_EQUAL(sign, 1.0);

            using FlowView = std::remove_cv_t<std::remove_reference_t<
                decltype(rate)>>;

            using Component = FlowView::Component;
            using Direction = FlowView::Direction;

            BOOST_CHECK_CLOSE(rate.flow(Component::Oil), 0.12f, 1.0e-5);
            BOOST_CHECK_CLOSE(rate.flow(Component::Gas), 0.24f, 1.0e-5);
            BOOST_CHECK_CLOSE(rate.flow(Component::Water), 0.36f, 1.0e-6);
            BOOST_CHECK_CLOSE(rate.flow(Component::Disgas), 0.48f, 1.0e-5);
            BOOST_CHECK_CLOSE(rate.flow(Component::Vapoil), 1.25f, 1.0e-6);

            BOOST_CHECK_CLOSE(rate.flow(Component::Oil, Direction::Positive), 0.12f, 1.0e-5);
            BOOST_CHECK_CLOSE(rate.flow(Component::Gas, Direction::Positive), 0.24f, 1.0e-5);
            BOOST_CHECK_CLOSE(rate.flow(Component::Water, Direction::Positive), 0.36f, 1.0e-6);
            BOOST_CHECK_CLOSE(rate.flow(Component::Disgas, Direction::Positive), 0.48f, 1.0e-5);
            BOOST_CHECK_CLOSE(rate.flow(Component::Vapoil, Direction::Positive), 1.25f, 1.0e-6);

            BOOST_CHECK_CLOSE(rate.flow(Component::Oil, Direction::Negative), 0.0f, 1.0e-6);
            BOOST_CHECK_CLOSE(rate.flow(Component::Gas, Direction::Negative), 0.0f, 1.0e-6);
            BOOST_CHECK_CLOSE(rate.flow(Component::Water, Direction::Negative), 0.0f, 1.0e-6);
            BOOST_CHECK_CLOSE(rate.flow(Component::Disgas, Direction::Negative), 0.0f, 1.0e-6);
            BOOST_CHECK_CLOSE(rate.flow(Component::Vapoil, Direction::Negative), 0.0f, 1.0e-6);
        }

        {
            const auto q13 = map.getInterRegFlows(0, 2);
            BOOST_REQUIRE_MESSAGE(q13.has_value(),
                                  "There must be inter-region flows "
                                  "between regions 1 and 3");

            const auto& [rate, sign] = q13.value();

            BOOST_CHECK_EQUAL(sign, 1.0);

            using FlowView = std::remove_cv_t<std::remove_reference_t<
                decltype(rate)>>;

            using Component = FlowView::Component;
            using Direction = FlowView::Direction;

            BOOST_CHECK_CLOSE(rate.flow(Component::Oil), 0.2f, 1.0e-5);
            BOOST_CHECK_CLOSE(rate.flow(Component::Gas), 0.4f, 1.0e-5);
            BOOST_CHECK_CLOSE(rate.flow(Component::Water), 0.6f, 1.0e-6);
            BOOST_CHECK_CLOSE(rate.flow(Component::Disgas), 0.8f, 1.0e-5);
            BOOST_CHECK_CLOSE(rate.flow(Component::Vapoil), 1.0f, 1.0e-6);

            BOOST_CHECK_CLOSE(rate.flow(Component::Oil, Direction::Positive), 0.2f, 1.0e-5);
            BOOST_CHECK_CLOSE(rate.flow(Component::Gas, Direction::Positive), 0.4f, 1.0e-5);
            BOOST_CHECK_CLOSE(rate.flow(Component::Water, Direction::Positive), 0.6f, 1.0e-6);
            BOOST_CHECK_CLOSE(rate.flow(Component::Disgas, Direction::Positive), 0.8f, 1.0e-5);
            BOOST_CHECK_CLOSE(rate.flow(Component::Vapoil, Direction::Positive), 1.0f, 1.0e-6);

            BOOST_CHECK_CLOSE(rate.flow(Component::Oil, Direction::Negative), 0.0f, 1.0e-6);
            BOOST_CHECK_CLOSE(rate.flow(Component::Gas, Direction::Negative), 0.0f, 1.0e-6);
            BOOST_CHECK_CLOSE(rate.flow(Component::Water, Direction::Negative), 0.0f, 1.0e-6);
            BOOST_CHECK_CLOSE(rate.flow(Component::Disgas, Direction::Negative), 0.0f, 1.0e-6);
            BOOST_CHECK_CLOSE(rate.flow(Component::Vapoil, Direction::Negative), 0.0f, 1.0e-6);
        }

        {
            const auto q14 = map.getInterRegFlows(0, 3);
            BOOST_CHECK_MESSAGE(! q14.has_value(),
                                "There must NOT be inter-region "
                                "flows between regions 1 and 4");
        }

        {
            const auto q23 = map.getInterRegFlows(1, 2);
            BOOST_CHECK_MESSAGE(! q23.has_value(),
                                "There must NOT be inter-region "
                                "flows between regions 2 and 3");
        }

        {
            const auto q24 = map.getInterRegFlows(1, 3);
            BOOST_REQUIRE_MESSAGE(q24.has_value(),
                                  "There must be inter-region flows "
                                  "between regions 2 and 4");

            const auto& [rate, sign] = q24.value();

            BOOST_CHECK_EQUAL(sign, 1.0);

            using FlowView = std::remove_cv_t<std::remove_reference_t<
                decltype(rate)>>;

            using Component = FlowView::Component;
            using Direction = FlowView::Direction;

            BOOST_CHECK_CLOSE(rate.flow(Component::Oil), -0.1f, 1.0e-5);
            BOOST_CHECK_CLOSE(rate.flow(Component::Gas), -0.2f, 1.0e-5);
            BOOST_CHECK_CLOSE(rate.flow(Component::Water), -0.3f, 1.0e-6);
            BOOST_CHECK_CLOSE(rate.flow(Component::Disgas), -0.4f, 1.0e-5);
            BOOST_CHECK_CLOSE(rate.flow(Component::Vapoil), -0.5f, 1.0e-6);

            BOOST_CHECK_CLOSE(rate.flow(Component::Oil, Direction::Positive), 0.0f, 1.0e-5);
            BOOST_CHECK_CLOSE(rate.flow(Component::Gas, Direction::Positive), 0.0f, 1.0e-5);
            BOOST_CHECK_CLOSE(rate.flow(Component::Water, Direction::Positive), 0.0f, 1.0e-6);
            BOOST_CHECK_CLOSE(rate.flow(Component::Disgas, Direction::Positive), 0.0f, 1.0e-5);
            BOOST_CHECK_CLOSE(rate.flow(Component::Vapoil, Direction::Positive), 0.0f, 1.0e-6);

            BOOST_CHECK_CLOSE(rate.flow(Component::Oil, Direction::Negative), -0.1f, 1.0e-6);
            BOOST_CHECK_CLOSE(rate.flow(Component::Gas, Direction::Negative), -0.2f, 1.0e-6);
            BOOST_CHECK_CLOSE(rate.flow(Component::Water, Direction::Negative), -0.3f, 1.0e-6);
            BOOST_CHECK_CLOSE(rate.flow(Component::Disgas, Direction::Negative), -0.4f, 1.0e-6);
            BOOST_CHECK_CLOSE(rate.flow(Component::Vapoil, Direction::Negative), -0.5f, 1.0e-6);
        }
    }
}

BOOST_AUTO_TEST_CASE(Four_Processes_Declared_MaxID)
{
    using Map = Opm::InterRegFlowMap;

    const auto fipnum = all_same_region();
    const auto fipspl = left_right_split_region();
    const auto fipchk = checker_board_region();
    const auto fipsep = all_separate_region();

    auto rank = std::vector<Map>(4, Map {
            fipnum.size(), {
                { "FIPNUM", std::cref(fipnum) },
                { "FIPSPL", std::cref(fipspl) },
                { "FIPCHK", std::cref(fipchk) },
                { "FIPSEP", std::cref(fipsep) },
            }, 42});

    addConnections(FourProc::P1::isInterior, rank[0]);
    rank[0].assignGlobalMaxRegionID({3, 2, 50, 4});

    addConnections(FourProc::P2::isInterior, rank[1]);
    addConnections(FourProc::P3::isInterior, rank[2]);

    addConnections(FourProc::P4::isInterior, rank[3]);
    rank[3].assignGlobalMaxRegionID({5, 2, 2, 4});

    for (auto& r : rank) {
        r.compress();
    }

    auto buffer = makeMessageBuffer();

    rank[1].write(buffer);
    rank[0].write(buffer);
    rank[3].write(buffer);
    rank[2].read(buffer);
    rank[2].read(buffer);
    rank[2].read(buffer);

    BOOST_CHECK_MESSAGE(rank[2].readIsConsistent(),
                        "Consistent write() must yield consistent read()");

    BOOST_CHECK_THROW(rank[2].addConnection({ 0, TwoProc::P1::isInterior(0) },
                                            { 1, TwoProc::P1::isInterior(1) }, conn_01()),
                      std::logic_error);

    rank[2].compress();

    const auto iregFlows = rank[2].getInterRegFlows();
    BOOST_REQUIRE_EQUAL(iregFlows.size(), std::size_t{4});

    // FIPNUM
    {
        const auto& map = iregFlows[0];
        BOOST_CHECK_EQUAL(map.numRegions(), 42);

        for (auto r1 = 0*map.numRegions(); r1 < map.numRegions(); ++r1) {
            for (auto r2 = r1 + 1; r2 < map.numRegions(); ++r2) {
                auto flow = map.getInterRegFlows(r1, r2);
                BOOST_CHECK_MESSAGE(! flow.has_value(),
                                    "There must not be inter-regional flow "
                                    "between regions " << r1 << " and " << r2);
            }
        }
    }

    // FIPSPL
    BOOST_CHECK_EQUAL(iregFlows[1].numRegions(), 42);

    // FIPCHK
    BOOST_CHECK_EQUAL(iregFlows[2].numRegions(), 50);

    // FIPSEP
    BOOST_CHECK_EQUAL(iregFlows[3].numRegions(), 42);
}

BOOST_AUTO_TEST_SUITE_END() // MultiArray_Wrapper
