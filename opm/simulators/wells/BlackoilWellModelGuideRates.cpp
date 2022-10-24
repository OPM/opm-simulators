/*
  Copyright 2016 SINTEF ICT, Applied Mathematics.
  Copyright 2016 - 2017 Statoil ASA.
  Copyright 2017 Dr. Blatt - HPC-Simulation-Software & Services
  Copyright 2016 - 2018 IRIS AS

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
#include <opm/simulators/wells/BlackoilWellModelGuideRates.hpp>

#include <opm/input/eclipse/Schedule/Schedule.hpp>
#include <opm/input/eclipse/Schedule/Group/Group.hpp>

#include <opm/output/data/Groups.hpp>
#include <opm/output/data/GuideRateValue.hpp>

#include <opm/simulators/wells/BlackoilWellModelGeneric.hpp>
#include <opm/simulators/wells/WellGroupHelpers.hpp>
#include <opm/simulators/wells/WellInterfaceGeneric.hpp>

#include <fmt/format.h>

#include <cstddef>
#include <functional>
#include <stack>
#include <stdexcept>
#include <string_view>
#include <tuple>
#include <unordered_map>
#include <unordered_set>
#include <utility>

namespace {

struct RetrieveWellGuideRate
{
    RetrieveWellGuideRate() = default;

    explicit RetrieveWellGuideRate(const Opm::GuideRate& guideRate,
                                   const std::string&    wgname);

    explicit RetrieveWellGuideRate(const Opm::GuideRate& guideRate,
                                   const Opm::Group&     group);

    bool prod      { false };
    bool inj_water { false };
    bool inj_gas   { false };
};

RetrieveWellGuideRate
operator||(RetrieveWellGuideRate lhs, const RetrieveWellGuideRate& rhs)
{
    lhs.prod      = lhs.prod      || rhs.prod;
    lhs.inj_water = lhs.inj_water || rhs.inj_water;
    lhs.inj_gas   = lhs.inj_gas   || rhs.inj_gas;

    return lhs;
}

RetrieveWellGuideRate::RetrieveWellGuideRate(const Opm::GuideRate& guideRate,
                                             const std::string&    wgname)
    : prod      { guideRate.has(wgname) }
    , inj_water { guideRate.has(wgname, Opm::Phase::WATER) }
    , inj_gas   { guideRate.has(wgname, Opm::Phase::GAS)   }
{}

RetrieveWellGuideRate::RetrieveWellGuideRate(const Opm::GuideRate& guideRate,
                                             const Opm::Group&     group)
    : RetrieveWellGuideRate{ guideRate, group.name() }
{
    if (group.isProductionGroup()) {
        this->prod = true;
    }

    if (group.isInjectionGroup()) {
        this->inj_water = this->inj_water || group.hasInjectionControl(Opm::Phase::WATER);
        this->inj_gas   = this->inj_gas   || group.hasInjectionControl(Opm::Phase::GAS);
    }
}

class GroupTreeWalker
{
public:
    using GroupOp = std::function<void(const Opm::Group&)>;
    using WellOp = std::function<void(const Opm::Well&)>;

    explicit GroupTreeWalker(const Opm::Schedule& sched,
                             const int            reportStepIdx)
        : sched_        (sched)
        , reportStepIdx_(reportStepIdx)
    {}

    GroupTreeWalker& groupOp(GroupOp visit)
    {
        this->visitGroup_ = std::move(visit);
        return *this;
    }

    GroupTreeWalker& wellOp(WellOp visit)
    {
        this->visitWell_ = std::move(visit);
        return *this;
    }

    void clear()
    {
        this->visitGroup_ = GroupOp{};
        this->visitWell_ = WellOp{};
    }

    void traversePreOrder();
    void traversePostOrder();

private:
    using NodeOp = void (GroupTreeWalker::*)(std::string_view) const;

    std::reference_wrapper<const Opm::Schedule> sched_;
    int reportStepIdx_;

    GroupOp visitGroup_{};
    WellOp visitWell_{};

    std::stack<std::string_view, std::vector<std::string_view>> dfsGroupStack_{};
    std::unordered_set<std::size_t> dfsGroupDiscovered_{};

    NodeOp postDiscover_{nullptr};
    NodeOp preFinish_{nullptr};

    void traverse();

    void startWalk();
    void discover(std::string_view group);
    void finish(std::string_view group);

    bool isSeen(std::string_view group) const;
    std::size_t insertIndex(std::string_view group) const;

    void visitGroup(std::string_view group) const;
    void visitWell(std::string_view well) const;

    const Opm::Group& getGroup(std::string_view group) const;
    const Opm::Well& getWell(std::string_view well) const;
};

void GroupTreeWalker::traversePreOrder()
{
    this->preFinish_ = nullptr;
    this->postDiscover_ = &GroupTreeWalker::visitGroup;

    this->traverse();
}

void GroupTreeWalker::traversePostOrder()
{
    this->preFinish_ = &GroupTreeWalker::visitGroup;
    this->postDiscover_ = nullptr;

    this->traverse();
}

void GroupTreeWalker::traverse()
{
    this->startWalk();

    while (! this->dfsGroupStack_.empty()) {
        const auto gname = this->dfsGroupStack_.top();

        if (this->isSeen(gname)) {
            if (this->preFinish_ != nullptr) {
                (this->*preFinish_)(gname);
            }

            this->finish(gname);
            continue;
        }

        this->discover(gname);

        if (this->postDiscover_ != nullptr) {
            (this->*postDiscover_)(gname);
        }

        const auto& group = this->getGroup(gname);

        if (! group.wellgroup()) { // Node group.  Register child groups.
            for (const auto& child : group.groups()) {
                if (! this->isSeen(child)) {
                    this->dfsGroupStack_.push(child);
                }
            }
        }
        else { // Group is a well group--visit its wells.
            for (const auto& well : group.wells()) {
                this->visitWell(well);
            }
        }
    }
}

void GroupTreeWalker::startWalk()
{
    this->dfsGroupDiscovered_.clear();

    while (! this->dfsGroupStack_.empty()) {
        this->dfsGroupStack_.pop();
    }

    this->dfsGroupStack_.push("FIELD");
}

void GroupTreeWalker::discover(std::string_view group)
{
    this->dfsGroupDiscovered_.insert(this->insertIndex(group));
}

void GroupTreeWalker::finish(std::string_view group)
{
    if (this->dfsGroupStack_.top() != group) {
        throw std::invalid_argument {
            fmt::format("Internal Error: Expected group '{}', but got '{}'",
                        group, this->dfsGroupStack_.top())
        };
    }

    this->dfsGroupStack_.pop();
}

bool GroupTreeWalker::isSeen(std::string_view group) const
{
    return this->dfsGroupDiscovered_.find(this->insertIndex(group))
        != this->dfsGroupDiscovered_.end();
}

std::size_t GroupTreeWalker::insertIndex(std::string_view group) const
{
    return this->getGroup(group).insert_index();
}

void GroupTreeWalker::visitGroup(std::string_view group) const
{
    if (! this->visitGroup_) {
        return;
    }

    this->visitGroup_(this->getGroup(group));
}

void GroupTreeWalker::visitWell(std::string_view well) const
{
    if (! this->visitWell_) {
        return;
    }

    this->visitWell_(this->getWell(well));
}

const Opm::Group& GroupTreeWalker::getGroup(std::string_view group) const
{
    return this->sched_.get().getGroup({group.data(), group.size()}, this->reportStepIdx_);
}

const Opm::Well& GroupTreeWalker::getWell(std::string_view well) const
{
    return this->sched_.get().getWell({well.data(), well.size()}, this->reportStepIdx_);
}

} // anonymous namespace

namespace Opm {

void
BlackoilWellModelGuideRates::
getGuideRateValues(const GuideRate::RateVector& qs,
                   const bool                   is_inj,
                   const std::string&           wgname,
                   data::GuideRateValue&        grval) const
{
    auto getGR = [this, &wgname, &qs](const GuideRateModel::Target t)
    {
        return wellModel_.guideRate().getSI(wgname, t, qs);
    };

    // Note: GuideRate does currently (2020-07-20) not support Target::RES.
    grval.set(data::GuideRateValue::Item::Gas,
              getGR(GuideRateModel::Target::GAS));

    grval.set(data::GuideRateValue::Item::Water,
              getGR(GuideRateModel::Target::WAT));

    if (!is_inj) {
        // Producer.  Extract "all" guiderate values.
        grval.set(data::GuideRateValue::Item::Oil,
                  getGR(GuideRateModel::Target::OIL));
    }
}

data::GuideRateValue
BlackoilWellModelGuideRates::
getGuideRateValues(const Well& well) const
{
    auto grval = data::GuideRateValue{};

    const auto& wname = well.name();
    if (!wellModel_.wellState().has(wname)) {
        // No flow rates for 'wname' -- might be before well comes
        // online (e.g., for the initial condition before simulation
        // starts).
        return grval;
    }

    if (!wellModel_.guideRate().has(wname)) {
        // No guiderates exist for 'wname'.
        return grval;
    }

    const auto qs = WellGroupHelpers::
        getWellRateVector(wellModel_.wellState(), wellModel_.phaseUsage(), wname);

    this->getGuideRateValues(qs, well.isInjector(), wname, grval);

    return grval;
}

data::GuideRateValue
BlackoilWellModelGuideRates::
getGuideRateValues(const Group& group) const
{
    auto grval = data::GuideRateValue{};
    const auto& gname = group.name();

    if (!wellModel_.groupState().has_production_rates(gname)) {
        // No flow rates for production group 'gname' -- might be before
        // group comes online (e.g., for the initial condition before
        // simulation starts).
        return grval;
    }

    if (!wellModel_.guideRate().has(gname)) {
        // No guiderates exist for 'gname'.
        return grval;
    }

    const auto qs = WellGroupHelpers::
        getProductionGroupRateVector(wellModel_.groupState(), wellModel_.phaseUsage(), gname);

    const auto is_inj = false; // This procedure only applies to G*PGR.
    this->getGuideRateValues(qs, is_inj, gname, grval);

    return grval;
}

data::GuideRateValue
BlackoilWellModelGuideRates::
getGuideRateInjectionGroupValues(const Group& group) const
{
    auto grval = data::GuideRateValue{};

    const auto& gname = group.name();
    if (wellModel_.guideRate().has(gname, Phase::GAS)) {
        grval.set(data::GuideRateValue::Item::Gas,
                  wellModel_.guideRate().getSI(gname, Phase::GAS));
    }

    if (wellModel_.guideRate().has(gname, Phase::WATER)) {
        grval.set(data::GuideRateValue::Item::Water,
                  wellModel_.guideRate().getSI(gname, Phase::WATER));
    }

    return grval;
}

void BlackoilWellModelGuideRates::
assignWellGuideRates(data::Wells& wsrpt,
                     const int    reportStepIdx) const
{
    auto all = std::unordered_map<std::string, data::GuideRateValue>{};
    auto retrieve = std::unordered_map<std::string, RetrieveWellGuideRate>{};

    auto walker = GroupTreeWalker{wellModel_.schedule(), reportStepIdx};

    // Populates 'retrieve'.
    walker.groupOp([this, &retrieve](const Group& group)
    {
        const auto& gname = group.name();

        const auto parent = (gname == "FIELD")
            ? RetrieveWellGuideRate{}
            : retrieve[group.parent()];

        auto [elm, inserted] =
            retrieve.emplace(std::piecewise_construct,
                             std::forward_as_tuple(gname),
                             std::forward_as_tuple(wellModel_.guideRate(), group));

        if (inserted) {
            elm->second = elm->second || parent;
        }
    });

    // Populates 'all'.
    walker.wellOp([this, &retrieve, &all](const Well& well)
    {
        const auto& wname = well.name();

        const auto is_nontrivial =
            wellModel_.guideRate().has(wname) || wellModel_.guideRate().hasPotentials(wname);

        if (! (is_nontrivial && wellModel_.wellState().has(wname))) {
            all[wname].clear();
            return;
        }

        auto parent_pos = retrieve.find(well.groupName());
        const auto parent = (parent_pos == retrieve.end())
            ? RetrieveWellGuideRate{} // No entry for 'parent'--unexpected.
            : parent_pos->second;

        const auto get_gr = parent
            || RetrieveWellGuideRate{wellModel_.guideRate(), wname};

        const auto qs = WellGroupHelpers::
            getWellRateVector(wellModel_.wellState(), wellModel_.phaseUsage(), wname);

        auto getGR = [this, &wname, &qs](const GuideRateModel::Target t)
        {
            return wellModel_.guideRate().getSI(wname, t, qs);
        };

        auto& grval = all[wname];

        if (well.isInjector()) {
            if (get_gr.inj_gas) { // Well supports WGIGR
                grval.set(data::GuideRateValue::Item::Gas,
                          getGR(GuideRateModel::Target::GAS));
            }
            if (get_gr.inj_water) { // Well supports WWIGR
                grval.set(data::GuideRateValue::Item::Water,
                          getGR(GuideRateModel::Target::WAT));
            }
        }
        else if (get_gr.prod) { // Well is producer AND we want/support WxPGR
            grval
                .set(data::GuideRateValue::Item::Oil  , getGR(GuideRateModel::Target::OIL))
                .set(data::GuideRateValue::Item::Gas  , getGR(GuideRateModel::Target::GAS))
                .set(data::GuideRateValue::Item::Water, getGR(GuideRateModel::Target::WAT));
        }
    });

    // Visit groups before their children, meaning no well is visited until
    // all of its upline parent groups--up to FIELD--have been visited.
    // Upon completion, 'all' contains guide rate values for all wells
    // reachable from 'FIELD' at this time/report step.
    walker.traversePreOrder();

    for (const auto* well : wellModel_.genericWells()) {
        auto xwPos = wsrpt.find(well->name());
        if (xwPos == wsrpt.end()) { // No well results.  Unexpected.
            continue;
        }

        auto grPos = all.find(well->name());
        if (grPos == all.end()) {
            continue;
        }

        xwPos->second.guide_rates = grPos->second;
    }
}

std::unordered_map<std::string, data::GroupGuideRates>
BlackoilWellModelGuideRates::
calculateAllGroupGuideRates(const int reportStepIdx) const
{
    auto gr = std::unordered_map<std::string, data::GroupGuideRates>{};

    auto walker = GroupTreeWalker{wellModel_.schedule(), reportStepIdx};

    // Populates 'gr'.
    walker.groupOp([this, &gr](const Group& group)
    {
        const auto& gname = group.name();

        if (gname == "FIELD") { return; }

        if (wellModel_.guideRate().has(gname)) {
            gr[gname].production = this->getGuideRateValues(group);
        }

        if (wellModel_.guideRate().has(gname, Phase::WATER) ||
            wellModel_.guideRate().has(gname, Phase::GAS))
        {
            gr[gname].injection = this->getGuideRateInjectionGroupValues(group);
        }

        const auto parent = group.parent();
        if (parent == "FIELD") { return; }

        gr[parent].injection  += gr[gname].injection;
        gr[parent].production += gr[gname].production;
    });

    // Populates 'gr'.
    walker.wellOp([this, &gr](const Well& well)
    {
        if (! (wellModel_.guideRate().has(well.name()) ||
               wellModel_.guideRate().hasPotentials(well.name())))
        {
            return;
        }

        const auto& gname = well.groupName();

        auto& grval = well.isInjector()
            ? gr[gname].injection
            : gr[gname].production;

        grval += this->getGuideRateValues(well);
    });

    // Visit wells and groups before their parents, meaning no group is
    // visited until all of its children down to the leaves of the group
    // tree have been visited.  Upon completion, 'gr' contains guide rate
    // values for all groups reachable from 'FIELD' at this time/report
    // step.
    walker.traversePostOrder();

    return gr;
}

void BlackoilWellModelGuideRates::
assignGroupGuideRates(const Group& group,
                      const std::unordered_map<std::string, data::GroupGuideRates>& groupGuideRates,
                      data::GroupData& gdata) const
{
    auto& prod = gdata.guideRates.production;  prod.clear();
    auto& inj  = gdata.guideRates.injection;   inj .clear();

    auto xgrPos = groupGuideRates.find(group.name());
    if (xgrPos == groupGuideRates.end()) {
        // No guiderates defined for this group.
        return;
    }

    const auto& xgr = xgrPos->second;

    if (wellModel_.guideRate().has(group.name())) {
        prod = xgr.production;
    }

    if (wellModel_.guideRate().has(group.name(), Phase::WATER) ||
        wellModel_.guideRate().has(group.name(), Phase::GAS))
    {
        inj = xgr.injection;
    }
}

bool BlackoilWellModelGuideRates::
guideRateUpdateIsNeeded(const int reportStepIdx) const
{
    const auto& genWells = wellModel_.genericWells();
    auto need_update =
    std::any_of(genWells.begin(), genWells.end(),
    [](const WellInterfaceGeneric* well)
    {
        return well->changedToOpenThisStep();
    });
    if (!need_update && wellModel_.reportStepStarts()) {
        const auto& events = wellModel_.schedule()[reportStepIdx].wellgroup_events();
        constexpr auto effective_events_mask = ScheduleEvents::WELL_STATUS_CHANGE
            + ScheduleEvents::INJECTION_TYPE_CHANGED
            + ScheduleEvents::WELL_SWITCHED_INJECTOR_PRODUCER
            + ScheduleEvents::NEW_WELL;

        need_update = std::any_of(genWells.begin(), genWells.end(),
            [&events](const WellInterfaceGeneric* well)
        {
            return events.hasEvent(well->name(), effective_events_mask);
        });
    }
    return wellModel_.comm().max(static_cast<int>(need_update));
}

} // namespace Opm
