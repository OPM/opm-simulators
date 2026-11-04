// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*
  This file is part of the Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 2 of the License, or
  (at your option) any later version.

  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.

  Consult the COPYING file in the top-level source directory of this
  module for the precise wording of the license and the list of
  copyright holders.
*/

#include <config.h>

#include <opm/simulators/flow/ActionHandler.hpp>

#include <opm/common/OpmLog/OpmLog.hpp>
#include <opm/common/utility/TimeService.hpp>
#include <opm/common/TimingMacros.hpp>

#include <opm/input/eclipse/EclipseState/EclipseState.hpp>

#include <opm/input/eclipse/Schedule/Action/ActionContext.hpp>
#include <opm/input/eclipse/Schedule/Action/Actions.hpp>
#include <opm/input/eclipse/Schedule/Action/ActionX.hpp>
#include <opm/input/eclipse/Schedule/Action/SimulatorUpdate.hpp>
#include <opm/input/eclipse/Schedule/Action/State.hpp>
#include <opm/input/eclipse/Schedule/Schedule.hpp>
#include <opm/input/eclipse/EclipseState/InitConfig/InitConfig.hpp>
#include <opm/input/eclipse/Schedule/Group/Group.hpp>
#include <opm/input/eclipse/Schedule/UDQ/UDQConfig.hpp>
#include <opm/input/eclipse/Schedule/UDQ/UDQDefine.hpp>
#include <opm/input/eclipse/Schedule/UDQ/UDQEnums.hpp>
#include <opm/input/eclipse/Schedule/UDQ/UDQToken.hpp>
#include <opm/input/eclipse/Schedule/Well/Well.hpp>
#include <opm/input/eclipse/Schedule/Well/WellConnections.hpp>
#include <opm/input/eclipse/Schedule/Well/WellMatcher.hpp>

#include <opm/material/fluidsystems/BlackOilDefaultFluidSystemIndices.hpp>

#include <opm/simulators/wells/BlackoilWellModelGeneric.hpp>

#include <opm/simulators/utils/ParallelCommunication.hpp>
#include <opm/simulators/utils/ParallelSerialization.hpp>

#include <algorithm>
#include <cctype>
#include <chrono>
#include <cstdlib>
#include <map>
#include <sstream>
#include <unordered_set>
#include <variant>
#include <cstddef>
#include <ctime>
#include <string>
#include <unordered_map>
#include <vector>

#include <fmt/chrono.h>
#include <fmt/format.h>
#include <fmt/ranges.h>

namespace {
    std::string formatActionDate(const Opm::TimeStampUTC& timePoint,
                                 const int                reportStep)
    {
        auto time_point = std::tm{};

        time_point.tm_year = timePoint.year()  - 1900;
        time_point.tm_mon  = timePoint.month() -    1;
        time_point.tm_mday = timePoint.day();

        time_point.tm_hour = timePoint.hour();
        time_point.tm_min  = timePoint.minutes();
        time_point.tm_sec  = timePoint.seconds();

        return fmt::format("{:%d-%b-%Y %H:%M:%S} (report interval {} to {})",
                           time_point, reportStep, reportStep + 1);
    }

    void logActivePyAction(const std::string& actionName,
                           const std::string& timeString)
    {
        const auto message =
            fmt::format("Action {} (Python) triggered at {}",
                        actionName, timeString);

        Opm::OpmLog::info("ACTION_TRIGGERED", message);
    }

    template <typename WellNameRange>
    void logActiveAction(const std::string&   actionName,
                         const WellNameRange& matchingWells,
                         const std::string&   timeString)
    {
        const auto wellString = matchingWells.empty()
            ? std::string{}
            : fmt::format(" Well{}: {}",
                          matchingWells.size() != 1 ? "s" : "",
                          fmt::join(matchingWells, ", "));

        const auto message =
            fmt::format("Action {} triggered at {}.{}",
                        actionName, timeString, wellString);

        Opm::OpmLog::info("ACTION_TRIGGERED", message);
    }

    template <typename WellNameRange>
    void logReplayedAction(const std::string&   actionName,
                           const WellNameRange& matchingWells,
                           const std::string&   timeString,
                           const std::string&   parentTimeString)
    {
        const auto wellString = matchingWells.empty()
            ? std::string{}
            : fmt::format(" Well{}: {}",
                          matchingWells.size() != 1 ? "s" : "",
                          fmt::join(matchingWells, ", "));

        const auto message =
            fmt::format("Action {} triggered at {}, replayed: the parent run triggered it at {}.{}",
                        actionName, timeString, parentTimeString, wellString);

        Opm::OpmLog::info("ACTION_TRIGGERED", message);
    }

    void logInactivePyAction(const std::string& actionName,
                             const std::string& timeString)
    {
        const auto message =
            fmt::format("Action {} (Python) NOT triggered at {}.",
                        actionName, timeString);

        Opm::OpmLog::debug("NAMED_ACTION_NOT_TRIGGERED", message);
    }

    void logInactiveAction(const std::string& actionName,
                           const std::string& timeString)
    {
        const auto message =
            fmt::format("Action {} NOT triggered at {}.",
                        actionName, timeString);

        Opm::OpmLog::debug("NAMED_ACTION_NOT_TRIGGERED", message);
    }

    void logInactiveActions(const int          numInactive,
                            const std::string& timeString)
    {
        const auto message =
            fmt::format("{} action{} NOT triggered at {}.",
                        numInactive,
                        (numInactive != 1) ? "s" : "",
                        timeString);

        Opm::OpmLog::debug("ACTION_NOT_TRIGGERED", message);
    }

    template <typename Scalar, class WellModel>
    std::unordered_map<std::string, Scalar>
    fetchWellPI(const WellModel&                             wellModel,
                const std::vector<std::string>&              wellpi_wells,
                const Opm::Parallel::Communication           comm)
    {
        auto wellpi = std::unordered_map<std::string, Scalar> {};

        if (wellpi_wells.empty()) {
            return wellpi;
        }

        auto wellPI = std::vector<Scalar>(wellpi_wells.size());
        for (auto i = 0*wellpi_wells.size(); i < wellpi_wells.size(); ++i) {
            if (wellModel.hasLocalWell(wellpi_wells[i])) {
                wellPI[i] = wellModel.wellPI(wellpi_wells[i]);
            }
        }

        comm.max(wellPI.data(), wellPI.size());

        for (auto i = 0*wellpi_wells.size(); i < wellpi_wells.size(); ++i) {
            wellpi.emplace(wellpi_wells[i], wellPI[i]);
        }

        return wellpi;
    }
} // Anonymous namespace

namespace Opm {

template<typename Scalar, typename IndexTraits>
ActionHandler<Scalar, IndexTraits>::
ActionHandler(EclipseState& ecl_state,
              Schedule& schedule,
              Action::State& actionState,
              SummaryState& summaryState,
              BlackoilWellModelGeneric<Scalar, IndexTraits>& wellModel,
              Parallel::Communication comm)
    : ecl_state_(ecl_state)
    , schedule_(schedule)
    , actionState_(actionState)
    , summaryState_(summaryState)
    , wellModel_(wellModel)
    , comm_(comm)
{}

template<typename Scalar, typename IndexTraits>
void ActionHandler<Scalar, IndexTraits>::
applyActions(const int reportStep,
             const double sim_time,
             const TransFunc& transUp)
{
    OPM_TIMEBLOCK(applyActions);
    const auto& actions = schedule_[reportStep].actions();
    if (actions.empty()) {
        return;
    }

    const Action::Context context{ summaryState_, schedule_[reportStep].wlist_manager() };

    const auto now = TimeStampUTC{ schedule_.getStartTime() } + std::chrono::duration<double>(sim_time);
    const auto ts  = formatActionDate(now, reportStep);


    SimulatorUpdate sim_update;
    for (const auto& pyaction : actions.pending_python(actionState_)) {
        // The std::unordered_map<std::string, Scalar> wellpi contains the well production indices from the last
        // timestep. This map is needed for the keyword WELPI. For a PyAction, we do not know which wells are
        // affected by the PyAction, thus we get the production indices for all wells.
        const std::vector<std::string> wellpi_wells = schedule_[reportStep].well_order().names();
        const auto wellpi = fetchWellPI<Scalar>(this->wellModel_, wellpi_wells, this->comm_);
        auto sim_update_current = schedule_.runPyAction(reportStep, *pyaction, actionState_,
                                                ecl_state_, summaryState_, wellpi);

        if (const auto pyRes = this->actionState_.python_result(pyaction->name());
            !pyRes.has_value() || !*pyRes)
        {
            logInactivePyAction(pyaction->name(), ts);
            continue;
        }
        else {
            logActivePyAction(pyaction->name(), ts);
        }

        sim_update.append(sim_update_current);
    }

    auto non_triggered = 0;
    // Not from 'now': that is rounded down to whole seconds for display,
    // while an action's start time keeps its fraction of a second.
    const auto simTime = TimeService::advance(schedule_.getStartTime(), sim_time);
    for (const auto& action : actions.pending(this->actionState_, simTime)) {
        if (this->replayed_.count(action->name()) > 0) {
            continue;
        }
        const auto actionResult = action->eval(context);
        if (! actionResult.conditionSatisfied()) {
            ++non_triggered;
            logInactiveAction(action->name(), ts);
            continue;
        }

        const auto& matches = actionResult.matches();

        logActiveAction(action->name(), matches.wells(), ts);

        // The std::unordered_map<std::string, Scalar> wellpi contains the well production indices from the last
        // timestep. This map is needed for the keyword WELPI. For an ActionX, we know which wells are affected by
        // the ActionX, thus fetch the well production indices only for the affected wells.
        const std::vector<std::string> wellpi_wells = action->wellpi_wells(schedule_.wellMatcher(reportStep), matches);
        const auto wellpi = fetchWellPI<Scalar>(this->wellModel_, wellpi_wells, this->comm_);

        const auto sim_update_current = this->schedule_
            .applyAction(reportStep, *action, matches, wellpi, true);

        sim_update.append(sim_update_current);
        this->actionState_.add_run(*action, simTime, actionResult);
    }

    // Sector run: the actions the parent triggered up to the end of this
    // time step, in the order they are defined, with the wells it matched
    std::vector<const Action::ActionX*> replayCandidates;
    if (!this->replayTriggers_.empty()) {
        for (const auto& action : actions) {
            replayCandidates.push_back(&action);
        }
    }
    for (const auto* actionPtr : replayCandidates) {
        const auto& action = *actionPtr;
        const auto& name = action.name();
        const auto replayIt = this->replayTriggers_.find(name);
        if (replayIt == this->replayTriggers_.end()) {
            continue;
        }
        const auto& triggers = replayIt->second;
        auto& next = this->replayNext_[name];
        std::vector<std::string> wells;
        std::optional<double> parentTime;
        while ((next < triggers.size()) && (triggers[next].time <= sim_time + 1.0e-3)) {
            if (!parentTime.has_value()) {
                parentTime = triggers[next].time;
            }
            for (const auto& well : triggers[next].wells) {
                if (this->schedule_.hasWell(well, reportStep)
                    && (std::find(wells.begin(), wells.end(), well) == wells.end()))
                {
                    wells.push_back(well);
                }
            }
            ++next;
        }
        if (!parentTime.has_value()) {
            continue;
        }

        auto actionResult = Action::Result { true };
        if (!wells.empty()) {
            actionResult.wells(wells);
        }
        const auto& matches = actionResult.matches();

        const auto parentNow = TimeStampUTC{ schedule_.getStartTime() } + std::chrono::duration<double>(*parentTime);
        auto parentTs = formatActionDate(parentNow, reportStep);
        parentTs = parentTs.substr(0, parentTs.find(" (report"));
        logReplayedAction(name, matches.wells(), ts, parentTs);

        const std::vector<std::string> wellpi_wells = action.wellpi_wells(schedule_.wellMatcher(reportStep), matches);
        const auto wellpi = fetchWellPI<Scalar>(this->wellModel_, wellpi_wells, this->comm_);

        const auto sim_update_current = this->schedule_
            .applyAction(reportStep, action, matches, wellpi, true);

        sim_update.append(sim_update_current);
        this->actionState_.add_run(action, simTime, actionResult);
    }

    if (non_triggered > 0) {
        logInactiveActions(non_triggered, ts);
    }

    bool commit_wellstate = false;
    this->applySimulatorUpdate(reportStep, sim_update, transUp, commit_wellstate);
    // The well state has been stored in a previous object when the time
    // step has completed successfully, the action process might have
    // modified the well state, and to be certain that is not overwritten
    // when starting the next timestep we must commit it.
    if (commit_wellstate) {
        this->wellModel_.commitWGState();
    }
}

template<typename Scalar, typename IndexTraits>
void ActionHandler<Scalar, IndexTraits>::
applySimulatorUpdate(const int report_step,
                     const SimulatorUpdate& sim_update,
                     const TransFunc& updateTrans,
                     bool& commit_wellstate)
{
    OPM_TIMEBLOCK(applySimulatorUpdate);

    this->wellModel_.updateEclWells(report_step, sim_update, this->summaryState_);

    if (!sim_update.affected_wells.empty()) {
        commit_wellstate = true;
    }

    if (sim_update.tran_update) {
        const auto& keywords = schedule_[report_step].geo_keywords();
        ecl_state_.apply_schedule_keywords( keywords );
        eclBroadcast(comm_, ecl_state_.getTransMult() );

        // re-compute transmissibility
        updateTrans(true);
    }
}

template<typename Scalar, typename IndexTraits>
std::optional<std::string> ActionHandler<Scalar, IndexTraits>::
parentQuantityReason_(const std::string& quantity,
                      const std::vector<std::string>& args,
                      const std::unordered_set<std::string>& sectorWells,
                      std::unordered_map<std::string, std::optional<std::string>>& udqMemo) const
{
    if (quantity.empty()) {
        return std::nullopt;
    }

    // Numbers and the date
    char* end = nullptr;
    std::strtod(quantity.c_str(), &end);
    if ((end != nullptr) && (*end == '\0')) {
        return std::nullopt;
    }
    auto name = quantity;
    std::transform(name.begin(), name.end(), name.begin(),
                   [](unsigned char c) { return std::toupper(c); });
    static const std::unordered_set<std::string> dateWords {
        "DAY", "MNTH", "MONTH", "YEAR", "TIME", "TIMESTEP",
        "JAN", "FEB", "MAR", "APR", "MAY", "JUN", "JUL", "AUG", "SEP", "OCT", "NOV", "DEC",
    };
    if (dateWords.count(name) > 0) {
        return std::nullopt;
    }

    const auto lastStep = this->schedule_.size() - 1;
    const auto describe = [&name, &args](const std::string& what)
    {
        return args.empty()
            ? fmt::format("{} ({})", name, what)
            : fmt::format("{} {} ({})", name, fmt::join(args, " "), what);
    };

    // A UDQ is what it is defined from
    if ((name.size() > 1) && (name[1] == 'U')) {
        // The memo holds the non-UDQ quantity a UDQ depends on, if any
        const auto udqReason = [&name](const std::optional<std::string>& leaf)
            -> std::optional<std::string>
        {
            if (!leaf.has_value()) {
                return std::nullopt;
            }
            return fmt::format("{}, a UDQ depending on {}", name, *leaf);
        };
        if (const auto memo = udqMemo.find(name); memo != udqMemo.end()) {
            return udqReason(memo->second);
        }
        udqMemo[name] = std::nullopt; // against cycles

        std::unordered_set<std::string> definitions;
        for (std::size_t step = 0; step < this->schedule_.size(); ++step) {
            const auto& udq = this->schedule_[step].udq();
            if (!udq.has_keyword(name)) {
                continue;
            }
            const UDQDefine* define = nullptr;
            try {
                define = &udq.define(name);
            }
            catch (...) {
                continue; // assigned only
            }
            if (!definitions.insert(define->input_string()).second) {
                continue;
            }
            for (const auto& token : define->tokens()) {
                if ((token.type() != UDQTokenType::ecl_expr)
                    || !std::holds_alternative<std::string>(token.value()))
                {
                    continue;
                }
                auto tokenName = std::get<std::string>(token.value());
                if (auto reason = this->parentQuantityReason_(tokenName, token.selector(),
                                                              sectorWells, udqMemo))
                {
                    std::transform(tokenName.begin(), tokenName.end(), tokenName.begin(),
                                   [](unsigned char c) { return std::toupper(c); });
                    const auto inner = udqMemo.find(tokenName);
                    udqMemo[name] = (inner != udqMemo.end()) && inner->second.has_value()
                        ? inner->second
                        : reason;
                    return udqReason(udqMemo[name]);
                }
            }
        }
        return std::nullopt;
    }

    const auto outsideWell = [&](const std::vector<std::string>& wells) -> std::optional<std::string>
    {
        if (wells.empty()) {
            return std::string{};
        }
        for (const auto& well : wells) {
            if (sectorWells.count(well) == 0) {
                return well;
            }
        }
        return std::nullopt;
    };

    switch (name.front()) {
    case 'W':
    case 'C':
    case 'S': {
        std::vector<std::string> wells;
        if (args.empty()) {
            wells = this->schedule_.wellNames(lastStep);
        }
        else {
            try {
                wells = this->schedule_.wellNames(args.front(), lastStep);
            }
            catch (...) {
            }
        }
        if (const auto well = outsideWell(wells); well.has_value()) {
            return describe(well->empty()
                            ? "no well inside the sector matches"
                            : fmt::format("reads well {}, outside the sector", *well));
        }
        return std::nullopt;
    }

    case 'G': {
        if (args.empty()) {
            return describe("all groups, with wells outside the sector");
        }
        std::vector<std::string> wells;
        std::vector<std::string> pending { args.front() };
        while (!pending.empty()) {
            const auto group = pending.back();
            pending.pop_back();
            if (!this->schedule_.hasGroup(group, lastStep)) {
                return describe("unknown group");
            }
            const auto& g = this->schedule_.getGroup(group, lastStep);
            wells.insert(wells.end(), g.wells().begin(), g.wells().end());
            pending.insert(pending.end(), g.groups().begin(), g.groups().end());
        }
        if (const auto well = outsideWell(wells); well.has_value()) {
            return describe(well->empty()
                            ? "a group without wells"
                            : fmt::format("group with well {}, outside the sector", *well));
        }
        return std::nullopt;
    }

    case 'B': {
        std::string ijk;
        for (const auto& arg : args) {
            ijk += arg + " ";
        }
        std::replace(ijk.begin(), ijk.end(), ',', ' ');
        std::istringstream stream(ijk);
        int i = 0, j = 0, k = 0;
        const auto& grid = this->ecl_state_.getInputGrid();
        if ((stream >> i >> j >> k) && (i >= 1) && (j >= 1) && (k >= 1)
            && (static_cast<std::size_t>(i) <= grid.getNX())
            && (static_cast<std::size_t>(j) <= grid.getNY())
            && (static_cast<std::size_t>(k) <= grid.getNZ())
            && grid.cellActive(i - 1, j - 1, k - 1))
        {
            return std::nullopt;
        }
        return describe("a cell outside the sector");
    }

    case 'F':
        return describe("field quantity");

    case 'R':
        return describe("region quantity");

    default:
        return describe("not a quantity of the sector");
    }
}

template<typename Scalar, typename IndexTraits>
void ActionHandler<Scalar, IndexTraits>::
setupFluxActionReplay(const std::vector<std::string>& triggerNames,
                      const std::vector<double>& triggerTimes,
                      const std::vector<int>& triggerWellCounts,
                      const std::vector<std::string>& triggerWells,
                      const bool replay)
{
    this->replayed_.clear();
    this->replayTriggers_.clear();
    this->replayNext_.clear();

    // Every ACTIONX of the schedule, in the order they are first defined
    std::vector<std::string> order;
    std::unordered_set<std::string> hasActions;
    for (std::size_t step = 0; step < this->schedule_.size(); ++step) {
        for (const auto& action : this->schedule_[step].actions()) {
            if (hasActions.insert(action.name()).second) {
                order.push_back(action.name());
            }
        }
    }
    if (order.empty()) {
        return;
    }

    const auto optionsNote = std::string {
        "Use --flux-replay-actions=false to evaluate every ACTIONX condition in this run, "
        "and --flux-fpr-as-fprflux=true to evaluate FPR as FPRFLUX, the field pressure "
        "averaged over the whole parent model, in UDQs and ACTIONX conditions."
    };

    if (!replay) {
        OpmLog::info(fmt::format("USEFLUX: every ACTIONX condition is evaluated in this run "
                                 "(--flux-replay-actions=false); conditions on field, group or "
                                 "region quantities, or on wells and cells outside the sector, "
                                 "read the parent's values at the parent's sampling."));
        return;
    }
    if (triggerNames.empty()) {
        OpmLog::warning(fmt::format("USEFLUX: the FLUX file does not record which ACTIONX the parent "
                                    "triggered, so every ACTIONX condition is evaluated in this run. "
                                    "Re-create the FLUX file with make_flux next to the parent's PRT "
                                    "file to replay the parent's actions. {}", optionsNote));
        return;
    }

    // Classified on the I/O rank, which has the input grid, and broadcast
    std::string classification;
    if (this->comm_.rank() == 0) {
        std::unordered_set<std::string> sectorWells;
        for (std::size_t step = 0; step < this->schedule_.size(); ++step) {
            for (const auto& wellName : this->schedule_.wellNames(step)) {
                if (!this->schedule_.getWell(wellName, step).getConnections().empty()) {
                    sectorWells.insert(wellName);
                }
            }
        }

        std::unordered_map<std::string, std::optional<std::string>> reasons;
        std::unordered_map<std::string, std::optional<std::string>> udqMemo;
        std::unordered_set<std::string> checked;
        for (std::size_t step = 0; step < this->schedule_.size(); ++step) {
            for (const auto& action : this->schedule_[step].actions()) {
                auto& reason = reasons[action.name()];
                if (reason.has_value()) {
                    continue;
                }
                std::string key = action.name();
                for (const auto& condition : action.conditions()) {
                    key += "|" + condition.lhs.quantity() + fmt::format("{}", fmt::join(condition.lhs.args(), ","))
                         + "|" + condition.rhs.quantity() + fmt::format("{}", fmt::join(condition.rhs.args(), ","));
                }
                if (!checked.insert(key).second) {
                    continue;
                }
                for (const auto& condition : action.conditions()) {
                    for (const auto* quantity : { &condition.lhs, &condition.rhs }) {
                        if (!reason.has_value()) {
                            reason = this->parentQuantityReason_(quantity->quantity(), quantity->args(),
                                                                 sectorWells, udqMemo);
                        }
                    }
                }
            }
        }

        for (const auto& name : order) {
            const auto& reason = reasons[name];
            classification += fmt::format("{}\t{}\n", name, reason.has_value() ? *reason : std::string{"-"});
        }
    }

    {
        int size = static_cast<int>(classification.size());
        this->comm_.broadcast(&size, 1, 0);
        classification.resize(size);
        if (size > 0) {
            this->comm_.broadcast(classification.data(), size, 0);
        }
    }

    std::map<std::string, int> triggerCount;
    for (const auto& name : triggerNames) {
        ++triggerCount[name];
    }

    std::vector<std::string> replayedLines;
    std::vector<std::string> evaluatedLines;
    std::istringstream lines(classification);
    std::string line;
    while (std::getline(lines, line)) {
        const auto tab = line.find('\t');
        const auto name = line.substr(0, tab);
        const auto reason = line.substr(tab + 1);
        if (reason == "-") {
            evaluatedLines.push_back(fmt::format("  {}: reads only quantities of the sector", name));
            continue;
        }
        this->replayed_.insert(name);
        const auto count = triggerCount.count(name) ? triggerCount[name] : 0;
        replayedLines.push_back(fmt::format("  {}: reads {}; the parent triggered it {} time{}",
                                            name, reason, count, (count != 1) ? "s" : ""));
    }

    std::size_t wellIndex = 0;
    for (std::size_t t = 0; t < triggerNames.size(); ++t) {
        const auto count = static_cast<std::size_t>(std::max(triggerWellCounts[t], 0));
        if (this->replayed_.count(triggerNames[t]) > 0) {
            this->replayTriggers_[triggerNames[t]].push_back(
                { triggerTimes[t],
                  { triggerWells.begin() + wellIndex, triggerWells.begin() + wellIndex + count } });
        }
        wellIndex += count;
    }
    const auto& initConfig = this->ecl_state_.getInitConfig();
    const double restartTime = initConfig.restartRequested()
        ? this->schedule_.seconds(initConfig.getRestartStep())
        : 0.0;
    for (auto& [name, triggers] : this->replayTriggers_) {
        std::stable_sort(triggers.begin(), triggers.end(),
                         [](const auto& a, const auto& b) { return a.time < b.time; });
        // A restarted run has replayed the triggers before the restart
        auto& next = this->replayNext_[name];
        next = 0;
        if (restartTime > 0.0) {
            while ((next < triggers.size()) && (triggers[next].time <= restartTime)) {
                ++next;
            }
        }
    }

    OpmLog::info(fmt::format("USEFLUX: ACTIONX replay (--flux-replay-actions=true, the default).\n"
                             "{} action(s) are triggered when the parent run triggered them, with the "
                             "wells it matched, instead of being evaluated, because their conditions "
                             "read quantities this run does not compute itself, but takes from the "
                             "parent at the parent's sampling:\n{}\n"
                             "{} action(s) are evaluated in this run, because their conditions read "
                             "only quantities of the sector:\n{}\n{}",
                             replayedLines.size(), replayedLines.empty() ? "  (none)" : fmt::format("{}", fmt::join(replayedLines, "\n")),
                             evaluatedLines.size(), evaluatedLines.empty() ? "  (none)" : fmt::format("{}", fmt::join(evaluatedLines, "\n")),
                             optionsNote));
}

template<typename Scalar, typename IndexTraits>
void ActionHandler<Scalar, IndexTraits>::
evalUDQAssignments(const unsigned episodeIdx,
                   UDQState& udq_state)
{
    this->schedule_[episodeIdx].udq()
        .eval_assign(this->schedule_.wellMatcher(episodeIdx),
                     this->schedule_[episodeIdx].group_order(),
                     this->schedule_.segmentMatcherFactory(episodeIdx),
                     this->summaryState_,
                     udq_state);
}

template class ActionHandler<double, BlackOilDefaultFluidSystemIndices>;

#if FLOW_INSTANTIATE_FLOAT
template class ActionHandler<float, BlackOilDefaultFluidSystemIndices>;
#endif

} // namespace Opm
