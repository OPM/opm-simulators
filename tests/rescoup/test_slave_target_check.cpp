/*
  Copyright 2026 Equinor ASA

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

#define BOOST_TEST_MODULE ResCoup_SlaveTargetCheck
#include <boost/test/unit_test.hpp>

#include <opm/simulators/flow/rescoup/ReservoirCouplingSlaveTargetCheck.hpp>

#include <opm/input/eclipse/Deck/Deck.hpp>
#include <opm/input/eclipse/EclipseState/Aquifer/AquiferConfig.hpp>
#include <opm/input/eclipse/EclipseState/Aquifer/NumericalAquifer/NumericalAquifers.hpp>
#include <opm/input/eclipse/EclipseState/Grid/EclipseGrid.hpp>
#include <opm/input/eclipse/EclipseState/Grid/FieldPropsManager.hpp>
#include <opm/input/eclipse/EclipseState/Runspec.hpp>
#include <opm/input/eclipse/EclipseState/SummaryConfig/SummaryConfig.hpp>
#include <opm/input/eclipse/EclipseState/Tables/TableManager.hpp>
#include <opm/input/eclipse/Parser/Parser.hpp>
#include <opm/input/eclipse/Python/Python.hpp>
#include <opm/input/eclipse/Schedule/Schedule.hpp>

#include <map>
#include <memory>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

using Opm::ReservoirCoupling::UnsupportedTargetUse;

// A slave run's schedule and summary configuration, parsed from a deck
// string with SUMMARY and SCHEDULE sections.
struct SlaveInput
{
    explicit SlaveInput(const std::string& deck_string, const bool slave_mode = true)
        : deck { Opm::Parser{}.parseString("RUNSPEC\nDIMENS\n 10 10 10 /\n" + deck_string) }
        , grid { 10, 10, 10 }
        , tables { deck }
        , field_props { deck, Opm::Phases{true, true, true}, grid, tables }
        , schedule { deck, grid, field_props, Opm::NumericalAquifers{},
                     Opm::Runspec{deck}, std::make_shared<Opm::Python>(),
                     /*lowActionParsingStrictness=*/false, slave_mode }
        , summary_config { deck, schedule, field_props, Opm::AquiferConfig{} }
    {}

    Opm::Deck deck;
    Opm::EclipseGrid grid;
    Opm::TableManager tables;
    Opm::FieldPropsManager field_props;
    Opm::Schedule schedule;
    Opm::SummaryConfig summary_config;
};

// Group tree shared by the test decks: slave groups MANI-B and MANI-C
// under PLAT-A, which is not a slave group.  GRUPSLAV follows per test.
const std::string group_tree = R"(
SCHEDULE

GRUPTREE
 'PLAT-A'  'FIELD' /
 'MANI-B'  'PLAT-A' /
 'MANI-C'  'PLAT-A' /
/
)";

const std::string both_slave_groups = R"(
GRUPSLAV
 'MANI-B'  'B1-M' /
 'MANI-C'  'C1-M' /
/
)";

std::vector<UnsupportedTargetUse> findUses(const std::string& deck_string)
{
    return Opm::ReservoirCoupling::findUnsupportedSlaveGroupTargetUses
        (SlaveInput{deck_string}.schedule);
}

// (vector, group) pairs of the uses, to compare without the source text.
std::set<std::pair<std::string, std::string>>
vectorsAndGroups(const std::vector<UnsupportedTargetUse>& uses)
{
    auto result = std::set<std::pair<std::string, std::string>>{};
    for (const auto& use : uses) {
        result.emplace(use.vector, use.group);
    }
    return result;
}

} // Anonymous namespace

BOOST_AUTO_TEST_CASE(UdqDefineNamedSlaveGroup)
{
    const auto uses = findUses(group_tree + both_slave_groups + R"(
UDQ
 DEFINE FUOPT GOPRT 'MANI-C' * 0.9 /
/
)");

    BOOST_REQUIRE_EQUAL(uses.size(), 1u);
    BOOST_CHECK_EQUAL(uses[0].vector, "GOPRT");
    BOOST_CHECK_EQUAL(uses[0].group, "MANI-C");
    BOOST_CHECK_MESSAGE(uses[0].source.find("UDQ FUOPT (DEFINE in") == 0,
                        "Unexpected source: " << uses[0].source);
}

BOOST_AUTO_TEST_CASE(SupportedAndOtherVectorsIgnored)
{
    // GGIRT and GWIRT are reported for slave groups; GOPR is a rate, not a
    // target; PLAT-A is not a slave group.
    const auto uses = findUses(group_tree + both_slave_groups + R"(
UDQ
 DEFINE FUGIT GGIRT 'MANI-C' + GWIRT 'MANI-B' /
 DEFINE FUOPR GOPR 'MANI-C' /
 DEFINE FUOPT GOPRT 'PLAT-A' /
/
)");

    BOOST_CHECK(uses.empty());
}

BOOST_AUTO_TEST_CASE(NameRootAndNoNameMatchSlaveGroups)
{
    const auto uses = findUses(group_tree + both_slave_groups + R"(
UDQ
 DEFINE GUGPT GGPRT 'MANI*' /
 DEFINE GUWPT GWPRT * 0.5 /
/
)");

    const auto expected = std::set<std::pair<std::string, std::string>> {
        {"GGPRT", "MANI-B"}, {"GGPRT", "MANI-C"},
        {"GWPRT", "MANI-B"}, {"GWPRT", "MANI-C"},
    };
    BOOST_CHECK(vectorsAndGroups(uses) == expected);
}

BOOST_AUTO_TEST_CASE(SlaveFlagExemptsPhase)
{
    // MANI-C: oil production and water/liquid production follow the slave's
    // own limits, and so does gas injection -- but not water injection, so
    // GVIRT (water plus gas) still counts.  MANI-B: both injection phases
    // follow the slave, so its GVIRT is fine.
    const auto uses = findUses(group_tree + R"(
GRUPSLAV
 'MANI-B'  'B1-M'  4*  1*  'SLAV'  'SLAV' /
 'MANI-C'  'C1-M'  'SLAV'  'SLAV'  4*  'SLAV' /
/

UDQ
 DEFINE FUC1 GOPRT 'MANI-C' + GWPRT 'MANI-C' + GLPRT 'MANI-C' /
 DEFINE FUC2 GGPRT 'MANI-C' + GVIRT 'MANI-C' /
 DEFINE FUB1 GVIRT 'MANI-B' /
/
)");

    const auto expected = std::set<std::pair<std::string, std::string>> {
        {"GGPRT", "MANI-C"}, {"GVIRT", "MANI-C"},
    };
    BOOST_CHECK(vectorsAndGroups(uses) == expected);
}

BOOST_AUTO_TEST_CASE(ActionxConditions)
{
    const auto uses = findUses(group_tree + both_slave_groups + R"(
ACTIONX
 'A1' 10 /
 GOPRT 'MANI-C' < 500 AND /
 GOPR 'MANI-B' > GLPRT 'MANI-B' /
/
ENDACTIO
)");

    const auto expected = std::set<std::pair<std::string, std::string>> {
        {"GOPRT", "MANI-C"}, {"GLPRT", "MANI-B"},
    };
    BOOST_CHECK(vectorsAndGroups(uses) == expected);
    for (const auto& use : uses) {
        BOOST_CHECK_EQUAL(use.source, "ACTIONX A1");
    }
}

BOOST_AUTO_TEST_CASE(FieldVectorOnlyForFieldSlaveGroup)
{
    const auto udq = std::string { R"(
UDQ
 DEFINE FUOPT FOPRT /
/
)" };

    BOOST_CHECK(findUses(group_tree + both_slave_groups + udq).empty());

    const auto uses = findUses(group_tree + R"(
GRUPSLAV
 'FIELD'  'F-M' /
/
)" + udq);
    const auto expected = std::set<std::pair<std::string, std::string>> {
        {"FOPRT", "FIELD"},
    };
    BOOST_CHECK(vectorsAndGroups(uses) == expected);
}

BOOST_AUTO_TEST_CASE(UseFoundOnceAcrossReportSteps)
{
    const auto uses = findUses(group_tree + both_slave_groups + R"(
UDQ
 DEFINE FUOPT GOPRT 'MANI-C' /
/

TSTEP
 1 /

UDQ
 DEFINE FUX 1.0 /
/

TSTEP
 1 /
)");

    BOOST_CHECK_EQUAL(uses.size(), 1u);
}

BOOST_AUTO_TEST_CASE(NoSlaveGroupsNoUses)
{
    const auto input = SlaveInput { group_tree + R"(
UDQ
 DEFINE FUOPT GOPRT 'MANI-C' /
/
)", /*slave_mode=*/false };

    BOOST_CHECK(Opm::ReservoirCoupling::findUnsupportedSlaveGroupTargetUses
                (input.schedule).empty());
    BOOST_CHECK_NO_THROW(Opm::ReservoirCoupling::checkSlaveGroupTargetVectors
                         (input.schedule, input.summary_config));
}

BOOST_AUTO_TEST_CASE(SummaryRequestsOnlyWarn)
{
    const auto input = SlaveInput { R"(
SUMMARY
GOPRT
/
GGIRT
/
)" + group_tree + R"(
GRUPSLAV
 'MANI-B'  'B1-M'  'SLAV' /
 'MANI-C'  'C1-M' /
/
)" };

    const auto expected = std::map<std::string, std::set<std::string>> {
        {"GOPRT", {"MANI-C"}},
    };
    BOOST_CHECK(Opm::ReservoirCoupling::findUnsupportedSlaveGroupTargetsInSummary
                (input.schedule, input.summary_config) == expected);
    BOOST_CHECK_NO_THROW(Opm::ReservoirCoupling::checkSlaveGroupTargetVectors
                         (input.schedule, input.summary_config));
}

BOOST_AUTO_TEST_CASE(CheckThrowsForUdqUse)
{
    const auto input = SlaveInput { group_tree + both_slave_groups + R"(
UDQ
 DEFINE FUOPT GOPRT 'MANI-C' /
/
)" };

    BOOST_CHECK_THROW(Opm::ReservoirCoupling::checkSlaveGroupTargetVectors
                      (input.schedule, input.summary_config), std::runtime_error);
}
