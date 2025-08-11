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
#include <opm/input/eclipse/Schedule/UDQ/UDQConfig.hpp>
#include <opm/input/eclipse/Schedule/Well/Well.hpp>
#include <opm/input/eclipse/Schedule/Well/WellMatcher.hpp>

#include <opm/material/fluidsystems/BlackOilFluidSystem.hpp>

#include <opm/models/blackoil/blackoilindices.hh>
#include <opm/models/blackoil/blackoilonephaseindices.hh>
#include <opm/models/blackoil/blackoiltwophaseindices.hh>



#include <opm/simulators/wells/BlackoilWellModelGeneric.hpp>

#include <opm/simulators/utils/ParallelCommunication.hpp>
#include <opm/simulators/utils/ParallelSerialization.hpp>

#include <chrono>
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

template<typename FluidSystem, typename Indices>
ActionHandler<FluidSystem, Indices>::
ActionHandler(EclipseState& ecl_state,
              Schedule& schedule,
              Action::State& actionState,
              SummaryState& summaryState,
              BlackoilWellModelGeneric<FluidSystem, Indices>& wellModel,
              Parallel::Communication comm)
    : ecl_state_(ecl_state)
    , schedule_(schedule)
    , actionState_(actionState)
    , summaryState_(summaryState)
    , wellModel_(wellModel)
    , comm_(comm)
{}

template<typename FluidSystem, typename Indices>
void ActionHandler<FluidSystem, Indices>::
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
    const auto simTime = asTimeT(now);
    for (const auto& action : actions.pending(this->actionState_, simTime)) {
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

template<typename FluidSystem, typename Indices>
void ActionHandler<FluidSystem, Indices>::
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

template<typename FluidSystem, typename Indices>
void ActionHandler<FluidSystem, Indices>::
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

    template<class Scalar>
    using FS = BlackOilFluidSystem<Scalar, BlackOilDefaultIndexTraits>;

#define INSTANTIATE(T,...) \
    template class ActionHandler<FS<T>, __VA_ARGS__>;

#define INSTANTIATE_TYPE(T)                                                  \
    INSTANTIATE(T,BlackOilOnePhaseIndices<0u,0u,0u,0u,false,false,0u,1u,0u>) \
    INSTANTIATE(T,BlackOilOnePhaseIndices<0u,0u,0u,1u,false,false,0u,1u,0u>) \
    INSTANTIATE(T,BlackOilOnePhaseIndices<0u,0u,0u,0u,false,false,0u,1u,5u>) \
    INSTANTIATE(T,BlackOilTwoPhaseIndices<0u,0u,0u,0u,false,false,0u,0u,0u>) \
    INSTANTIATE(T,BlackOilTwoPhaseIndices<0u,0u,0u,0u,false,false,0u,1u,0u>) \
    INSTANTIATE(T,BlackOilTwoPhaseIndices<0u,0u,0u,0u,false,false,0u,2u,0u>) \
    INSTANTIATE(T,BlackOilTwoPhaseIndices<0u,0u,1u,0u,false,false,0u,2u,0u>) \
    INSTANTIATE(T,BlackOilTwoPhaseIndices<0u,0u,2u,0u,false,false,0u,2u,0u>) \
    INSTANTIATE(T,BlackOilTwoPhaseIndices<0u,0u,0u,0u,false,true,0u,2u,0u>)  \
    INSTANTIATE(T,BlackOilTwoPhaseIndices<0u,0u,0u,1u,false,false,0u,1u,0u>) \
    INSTANTIATE(T,BlackOilTwoPhaseIndices<0u,0u,0u,0u,false,true,0u,0u,0u>)  \
    INSTANTIATE(T,BlackOilTwoPhaseIndices<0u,0u,0u,1u,false,false,0u,0u,0u>) \
    INSTANTIATE(T,BlackOilTwoPhaseIndices<0u,0u,0u,1u,false,true,0u,0u,0u>)  \
    INSTANTIATE(T,BlackOilTwoPhaseIndices<1u,0u,0u,0u,false,false,0u,0u,0u>) \
    INSTANTIATE(T,BlackOilIndices<0u,0u,0u,0u,false,false,0u,0u>)            \
    INSTANTIATE(T,BlackOilIndices<0u,0u,0u,0u,true,false,0u,0u>)             \
    INSTANTIATE(T,BlackOilIndices<0u,0u,0u,0u,false,true,0u,0u>)             \
    INSTANTIATE(T,BlackOilIndices<1u,0u,0u,0u,false,false,0u,0u>)            \
    INSTANTIATE(T,BlackOilIndices<0u,1u,0u,0u,false,false,0u,0u>)            \
    INSTANTIATE(T,BlackOilIndices<0u,0u,1u,0u,false,false,0u,0u>)            \
    INSTANTIATE(T,BlackOilIndices<0u,0u,0u,1u,false,false,0u,0u>)            \
    INSTANTIATE(T,BlackOilIndices<0u,0u,0u,1u,false,true,0u,0u>)             \
    INSTANTIATE(T,BlackOilIndices<1u,0u,0u,0u,true,false,0u,0u>)

    INSTANTIATE_TYPE(double)

#if FLOW_INSTANTIATE_FLOAT
    INSTANTIATE_TYPE(float)
#endif

} // namespace Opm
