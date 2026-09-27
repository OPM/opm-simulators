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
#include <opm/simulators/flow/rescoup/ReservoirCouplingSlaveTargetCheck.hpp>

#include <opm/common/OpmLog/OpmLog.hpp>
#include <opm/common/utility/shmatch.hpp>

#include <opm/input/eclipse/EclipseState/SummaryConfig/SummaryConfig.hpp>
#include <opm/input/eclipse/Schedule/Action/ActionX.hpp>
#include <opm/input/eclipse/Schedule/Action/Actions.hpp>
#include <opm/input/eclipse/Schedule/Action/Condition.hpp>
#include <opm/input/eclipse/Schedule/ResCoup/GrupSlav.hpp>
#include <opm/input/eclipse/Schedule/ResCoup/ReservoirCouplingInfo.hpp>
#include <opm/input/eclipse/Schedule/Schedule.hpp>
#include <opm/input/eclipse/Schedule/UDQ/UDQConfig.hpp>
#include <opm/input/eclipse/Schedule/UDQ/UDQDefine.hpp>
#include <opm/input/eclipse/Schedule/UDQ/UDQEnums.hpp>
#include <opm/input/eclipse/Schedule/UDQ/UDQToken.hpp>

#include <algorithm>
#include <array>
#include <cstddef>
#include <optional>
#include <stdexcept>
#include <string_view>
#include <variant>

#include <fmt/format.h>
#include <fmt/ranges.h>

namespace {

using FilterFlag = Opm::ReservoirCoupling::GrupSlav::FilterFlag;

// Group target vectors the slave does not report for a slave group whose
// limits come from the master run.  GGIRT and GWIRT are not here: the slave
// reports the master's injection targets through them.  The master sends
// surface rates only, so GVIRT has no master value to report.
constexpr auto unsupported_target_vectors = std::array {
    std::string_view{"GGPRT"},
    std::string_view{"GLPRT"},
    std::string_view{"GOPRT"},
    std::string_view{"GVIRT"},
    std::string_view{"GVPRT"},
    std::string_view{"GWPRT"},
};

// Group level name of an unsupported target vector: GOPRT for both GOPRT
// and FOPRT.  Nullopt for any other keyword.
std::optional<std::string> groupTargetVector(const std::string& keyword)
{
    if ((keyword.size() != 5) || ((keyword[0] != 'G') && (keyword[0] != 'F'))) {
        return std::nullopt;
    }
    auto vector = "G" + keyword.substr(1);
    if (std::ranges::find(unsupported_target_vectors, vector) == unsupported_target_vectors.end()) {
        return std::nullopt;
    }
    return vector;
}

// Whether the slave group's own limit applies to the vector's phase, so
// that the schedule value is the right target.  The production flags are
// mapped as in GroupStateHelper::getProductionFilterFlag_(): the second
// production flag of GRUPSLAV covers both water and liquid.
bool slaveLimitApplies(const Opm::ReservoirCoupling::GrupSlav& grup_slav,
                       const std::string& vector)
{
    if (vector == "GOPRT") {
        return grup_slav.oilProdFlag() == FilterFlag::SLAV;
    }
    if ((vector == "GWPRT") || (vector == "GLPRT")) {
        return grup_slav.liquidProdFlag() == FilterFlag::SLAV;
    }
    if (vector == "GGPRT") {
        return grup_slav.gasProdFlag() == FilterFlag::SLAV;
    }
    if (vector == "GVPRT") {
        return grup_slav.fluidVolumeProdFlag() == FilterFlag::SLAV;
    }
    // GVIRT: the summary evaluator adds the water and gas injection targets.
    return (grup_slav.waterInjFlag() == FilterFlag::SLAV)
        && (grup_slav.gasInjFlag() == FilterFlag::SLAV);
}

// Add the slave groups for which a summary keyword, followed by the group
// names or name roots in 'names', reads an unsupported target vector.
void addUses(std::set<Opm::ReservoirCoupling::UnsupportedTargetUse>& uses,
             const std::string& source,
             const std::string& keyword,
             const std::vector<std::string>& names,
             const Opm::ReservoirCoupling::CouplingInfo& rescoup)
{
    const auto vector = groupTargetVector(keyword);
    if (!vector.has_value()) {
        return;
    }
    const bool field_vector = (keyword[0] == 'F');
    for (const auto& [group, grup_slav] : rescoup.grupSlavs()) {
        const bool named = field_vector
            ? (group == "FIELD")
            : (names.empty() ||
               std::ranges::any_of(names, [&group](const auto& name)
                                   { return Opm::shmatch(name, group); }));
        if (named && !slaveLimitApplies(grup_slav, *vector)) {
            uses.insert({source, keyword, group});
        }
    }
}

} // Anonymous namespace

namespace Opm::ReservoirCoupling {

void checkSlaveGroupTargetVectors(const Schedule& schedule,
                                  const SummaryConfig& summary_config)
{
    const auto uses = findUnsupportedSlaveGroupTargetUses(schedule);
    if (!uses.empty()) {
        auto lines = std::vector<std::string>{};
        lines.reserve(uses.size());
        for (const auto& use : uses) {
            lines.push_back(fmt::format("  {}: {} for slave group {}",
                                        use.source, use.vector, use.group));
        }
        throw std::runtime_error {
            fmt::format("Reservoir coupling: the master run sets the targets of "
                        "the slave groups, and the slave reports them only through "
                        "GGIRT and GWIRT.  The following would read the slave's own "
                        "schedule value instead:\n{}\n"
                        "Use the slave's own limit for the phase by setting its "
                        "GRUPSLAV flag to SLAV, or remove these uses.",
                        fmt::join(lines, "\n"))
        };
    }

    for (const auto& [vector, groups] :
             findUnsupportedSlaveGroupTargetsInSummary(schedule, summary_config))
    {
        OpmLog::warning(fmt::format("Reservoir coupling: summary vector {} for "
                                    "slave group(s) {} reports the slave's own "
                                    "schedule value, not the target set by the "
                                    "master run.",
                                    vector, fmt::join(groups, ", ")));
    }
}

std::vector<UnsupportedTargetUse>
findUnsupportedSlaveGroupTargetUses(const Schedule& schedule)
{
    auto uses = std::set<UnsupportedTargetUse>{};

    // Consecutive report steps share the UDQ, ACTIONX and GRUPSLAV objects
    // until a keyword changes them, so only the steps that change one of
    // them need a look.
    const UDQConfig* prev_udq = nullptr;
    const Action::Actions* prev_actions = nullptr;
    const CouplingInfo* prev_rescoup = nullptr;
    for (std::size_t step = 0; step < schedule.size(); ++step) {
        const auto& state = schedule[step];
        const auto& udq = state.udq();
        const auto& actions = state.actions();
        const auto& rescoup = state.rescoup();
        if ((&udq == prev_udq) && (&actions == prev_actions) && (&rescoup == prev_rescoup)) {
            continue;
        }
        prev_udq = &udq;
        prev_actions = &actions;
        prev_rescoup = &rescoup;
        if (rescoup.grupSlavCount() == 0) {
            continue;
        }

        for (const auto& define : udq.definitions()) {
            const auto& location = define.location();
            const auto source = fmt::format("UDQ {} (DEFINE in {}, line {})",
                                            define.keyword(), location.filename,
                                            location.lineno);
            for (const auto& token : define.tokens()) {
                const auto* keyword = std::get_if<std::string>(&token.value());
                if ((token.type() == UDQTokenType::ecl_expr) && (keyword != nullptr)) {
                    addUses(uses, source, *keyword, token.selector(), rescoup);
                }
            }
        }

        for (const auto& action : actions) {
            const auto source = fmt::format("ACTIONX {}", action.name());
            for (const auto& condition : action.conditions()) {
                addUses(uses, source, condition.lhs.quantity(), condition.lhs.args(), rescoup);
                addUses(uses, source, condition.rhs.quantity(), condition.rhs.args(), rescoup);
            }
        }
    }

    return { uses.begin(), uses.end() };
}

std::map<std::string, std::set<std::string>>
findUnsupportedSlaveGroupTargetsInSummary(const Schedule& schedule,
                                          const SummaryConfig& summary_config)
{
    using Category = SummaryConfigNode::Category;

    // GRUPSLAV records are only ever added, never changed or removed, so the
    // last report step knows every slave group.
    const auto& rescoup = schedule.back().rescoup();
    auto result = std::map<std::string, std::set<std::string>>{};
    for (const auto& node : summary_config) {
        const auto category = node.category();
        if ((category != Category::Group) && (category != Category::Field)) {
            continue;
        }
        const auto vector = groupTargetVector(node.keyword());
        if (!vector.has_value()) {
            continue;
        }
        const auto group = (category == Category::Field)
            ? std::string{"FIELD"} : node.namedEntity();
        if (rescoup.hasGrupSlav(group) &&
            !slaveLimitApplies(rescoup.grupSlav(group), *vector))
        {
            result[node.keyword()].insert(group);
        }
    }
    return result;
}

} // namespace Opm::ReservoirCoupling
