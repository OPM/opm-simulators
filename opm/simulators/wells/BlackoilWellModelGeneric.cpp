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

#include <opm/simulators/wells/BlackoilWellModelGeneric.hpp>

#include <opm/output/data/GuideRateValue.hpp>
#include <opm/output/data/Groups.hpp>
#include <opm/output/data/Wells.hpp>
#include <opm/output/eclipse/RestartValue.hpp>

#include <opm/input/eclipse/EclipseState/EclipseState.hpp>
#include <opm/input/eclipse/Schedule/Group/GuideRateConfig.hpp>
#include <opm/input/eclipse/Schedule/Group/GuideRate.hpp>
#include <opm/input/eclipse/Schedule/Schedule.hpp>
#include <opm/input/eclipse/EclipseState/SummaryConfig/SummaryConfig.hpp>

#include <opm/simulators/utils/DeferredLogger.hpp>
#include <opm/simulators/wells/BlackoilWellModelConstraints.hpp>
#include <opm/simulators/wells/BlackoilWellModelRestart.hpp>
#include <opm/simulators/wells/GasLiftStage2.hpp>
#include <opm/simulators/wells/VFPProperties.hpp>
#include <opm/simulators/wells/WellGroupHelpers.hpp>
#include <opm/simulators/wells/WellInterfaceGeneric.hpp>
#include <opm/simulators/wells/WellState.hpp>

#include <algorithm>
#include <cassert>
#include <functional>
#include <stack>
#include <stdexcept>
#include <string_view>
#include <tuple>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include <fmt/format.h>

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
} // Anonymous

namespace Opm {

BlackoilWellModelGeneric::
BlackoilWellModelGeneric(Schedule& schedule,
                         const SummaryState& summaryState,
                         const EclipseState& eclState,
                         const PhaseUsage& phase_usage,
                         const Parallel::Communication& comm)
    : schedule_(schedule)
    , summaryState_(summaryState)
    , eclState_(eclState)
    , comm_(comm)
    , phase_usage_(phase_usage)
    , guideRate_(schedule)
    , active_wgstate_(phase_usage)
    , last_valid_wgstate_(phase_usage)
    , nupcol_wgstate_(phase_usage)
{

  const auto numProcs = comm_.size();
  this->not_on_process_ = [this, numProcs](const Well& well) {
      if (numProcs == decltype(numProcs){1})
          return false;

      // Recall: false indicates NOT active!
      const auto value = std::make_pair(well.name(), true);
      auto candidate = std::lower_bound(this->parallel_well_info_.begin(),
                                        this->parallel_well_info_.end(),
                                        value);

      return (candidate == this->parallel_well_info_.end())
              || (*candidate != value);
  };
}

int
BlackoilWellModelGeneric::
numLocalWells() const
{
    return wells_ecl_.size();
}

int
BlackoilWellModelGeneric::
numPhases() const
{
    return phase_usage_.num_phases;
}

bool
BlackoilWellModelGeneric::
hasWell(const std::string& wname) const
{
    return std::any_of(this->wells_ecl_.begin(), this->wells_ecl_.end(),
        [&wname](const Well& well)
    {
        return well.name() == wname;
    });
}

bool
BlackoilWellModelGeneric::
wellsActive() const
{
    return wells_active_;
}

bool
BlackoilWellModelGeneric::
anyMSWellOpenLocal() const
{
    for (const auto& well : wells_ecl_) {
        if (well.isMultiSegment()) {
            return true;
        }
    }
    return false;
}

const Well&
BlackoilWellModelGeneric::
getWellEcl(const std::string& well_name) const
{
    // finding the iterator of the well in wells_ecl
    auto well_ecl = std::find_if(wells_ecl_.begin(),
                                 wells_ecl_.end(),
                                 [&well_name](const Well& elem)->bool {
                                     return elem.name() == well_name;
                                 });

    assert(well_ecl != wells_ecl_.end());

    return *well_ecl;
}

void
BlackoilWellModelGeneric::
initFromRestartFile(const RestartValue& restartValues,
                    WellTestState wtestState,
                    const size_t numCells,
                    bool handle_ms_well)
{
    // The restart step value is used to identify wells present at the given
    // time step. Wells that are added at the same time step as RESTART is initiated
    // will not be present in a restart file. Use the previous time step to retrieve
    // wells that have information written to the restart file.
    const int report_step = std::max(eclState_.getInitConfig().getRestartStep() - 1, 0);

    const auto& config = this->schedule()[report_step].guide_rate();

    // wells_ecl_ should only contain wells on this processor.
    wells_ecl_ = getLocalWells(report_step);
    this->local_parallel_well_info_ = createLocalParallelWellInfo(wells_ecl_);

    this->initializeWellProdIndCalculators();
    initializeWellPerfData();

    if (! this->wells_ecl_.empty()) {
        handle_ms_well &= anyMSWellOpenLocal();
        // Resize for restart step
        this->wellState().resize(this->wells_ecl_, this->local_parallel_well_info_,
                                 this->schedule(), handle_ms_well, numCells,
                                 this->well_perf_data_, this->summaryState_);

        BlackoilWellModelRestart(*this).loadRestartData(restartValues.wells,
                                                        restartValues.grp_nwrk,
                                                        handle_ms_well,
                                                        this->wellState(),
                                                        this->groupState());

        if (config.has_model()) {
            BlackoilWellModelRestart(*this).loadRestartGuideRates(report_step,
                                                                  config.model().target(),
                                                                  restartValues.wells,
                                                                  this->guideRate_);
        }
    }

    if (config.has_model()) {
        BlackoilWellModelRestart(*this).loadRestartGuideRates(report_step,
                                                              config,
                                                              restartValues.grp_nwrk.groupData,
                                                              this->guideRate_);

        this->guideRate_.updateGuideRateExpiration(this->schedule().seconds(report_step), report_step);
    }

    this->active_wgstate_.wtest_state(std::move(wtestState));
    this->commitWGState();
    initial_step_ = false;
}

std::vector<Well>
BlackoilWellModelGeneric::
getLocalWells(const int timeStepIdx) const
{
    auto w = schedule().getWells(timeStepIdx);
    w.erase(std::remove_if(w.begin(), w.end(), not_on_process_), w.end());
    return w;
}

std::vector<std::reference_wrapper<ParallelWellInfo>>
BlackoilWellModelGeneric::
createLocalParallelWellInfo(const std::vector<Well>& wells)
{
    std::vector<std::reference_wrapper<ParallelWellInfo>> local_parallel_well_info;
    local_parallel_well_info.reserve(wells.size());
    for (const auto& well : wells)
    {
        auto wellPair = std::make_pair(well.name(), true);
        auto pwell = std::lower_bound(parallel_well_info_.begin(),
                                      parallel_well_info_.end(),
                                      wellPair);
        assert(pwell != parallel_well_info_.end() &&
               *pwell == wellPair);
        local_parallel_well_info.push_back(std::ref(*pwell));
    }
    return local_parallel_well_info;
}

void
BlackoilWellModelGeneric::
initializeWellProdIndCalculators()
{
    this->prod_index_calc_.clear();
    this->prod_index_calc_.reserve(this->wells_ecl_.size());
    for (const auto& well : this->wells_ecl_) {
        this->prod_index_calc_.emplace_back(well);
    }
}

void
BlackoilWellModelGeneric::
initializeWellPerfData()
{
    well_perf_data_.resize(wells_ecl_.size());
    int well_index = 0;
    for (const auto& well : wells_ecl_) {
        int connection_index = 0;
        // INVALID_ECL_INDEX marks no above perf available
        int connection_index_above = ParallelWellInfo::INVALID_ECL_INDEX;
        well_perf_data_[well_index].clear();
        well_perf_data_[well_index].reserve(well.getConnections().size());
        CheckDistributedWellConnections checker(well, local_parallel_well_info_[well_index].get());
        bool hasFirstConnection = false;
        bool firstOpenConnection = true;
        auto& parallelWellInfo = this->local_parallel_well_info_[well_index].get();
        parallelWellInfo.beginReset();

        for (const auto& connection : well.getConnections()) {
            const int active_index = compressedIndexForInterior(connection.global_index());
            if (connection.state() == Connection::State::OPEN) {
                if (active_index >= 0) {
                    if (firstOpenConnection)
                    {
                        hasFirstConnection = true;
                    }
                    checker.connectionFound(connection_index);
                    PerforationData pd;
                    pd.cell_index = active_index;
                    pd.connection_transmissibility_factor = connection.CF();
                    pd.satnum_id = connection.satTableId();
                    pd.ecl_index = connection_index;
                    well_perf_data_[well_index].push_back(pd);
                    parallelWellInfo.pushBackEclIndex(connection_index_above,
                                                      connection_index);
                }
                firstOpenConnection = false;
                // Next time this index is the one above as each open connection is
                // is stored somehwere.
                connection_index_above = connection_index;
            } else {
                checker.connectionFound(connection_index);
                if (connection.state() != Connection::State::SHUT) {
                    OPM_THROW(std::runtime_error,
                              "Connection state: " << Connection::State2String(connection.state()) << " not handled");
                }
            }
            // Note: we rely on the connections being filtered! I.e. there are only connections
            // to active cells in the global grid.
            ++connection_index;
        }
        parallelWellInfo.endReset();
        checker.checkAllConnectionsFound();
        parallelWellInfo.communicateFirstPerforation(hasFirstConnection);
        ++well_index;
    }
}

void
BlackoilWellModelGeneric::
checkGconsaleLimits(const Group& group,
                    WellState& well_state,
                    const int reportStepIdx,
                    DeferredLogger& deferred_logger)
{
     // call recursively down the group hiearchy
    for (const std::string& groupName : group.groups()) {
        checkGconsaleLimits( schedule().getGroup(groupName, reportStepIdx), well_state, reportStepIdx, deferred_logger);
    }

    // only for groups with gas injection controls
    if (!group.hasInjectionControl(Phase::GAS)) {
        return;
    }

    // check if gconsale is used for this group
    if (!schedule()[reportStepIdx].gconsale().has(group.name()))
        return;

    std::ostringstream ss;

    const auto& gconsale = schedule()[reportStepIdx].gconsale().get(group.name(), summaryState_);
    const Group::ProductionCMode& oldProductionControl = this->groupState().production_control(group.name());

    int gasPos = phase_usage_.phase_pos[BlackoilPhases::Vapour];
    double production_rate = WellGroupHelpers::sumWellSurfaceRates(group, schedule(), well_state, reportStepIdx, gasPos, /*isInjector*/false);
    double injection_rate = WellGroupHelpers::sumWellSurfaceRates(group, schedule(), well_state, reportStepIdx, gasPos, /*isInjector*/true);

    // sum over all nodes
    injection_rate = comm_.sum(injection_rate);
    production_rate = comm_.sum(production_rate);

    double sales_rate = production_rate - injection_rate;
    double production_target = gconsale.sales_target + injection_rate;

    // add import rate and subtract consumption rate for group for gas
    if (schedule()[reportStepIdx].gconsump().has(group.name())) {
        const auto& gconsump = schedule()[reportStepIdx].gconsump().get(group.name(), summaryState_);
        if (phase_usage_.phase_used[BlackoilPhases::Vapour]) {
            sales_rate += gconsump.import_rate;
            sales_rate -= gconsump.consumption_rate;
            production_target -= gconsump.import_rate;
            production_target += gconsump.consumption_rate;
        }
    }

    if (sales_rate > gconsale.max_sales_rate) {
        switch(gconsale.max_proc) {
        case GConSale::MaxProcedure::NONE: {
            if (oldProductionControl != Group::ProductionCMode::GRAT && oldProductionControl != Group::ProductionCMode::NONE) {
                ss << "Group sales exceed maximum limit, but the action is NONE for " + group.name() + ". Nothing happens";
            }
            break;
            }
        case GConSale::MaxProcedure::CON: {
            OPM_DEFLOG_THROW(std::runtime_error, "Group " + group.name() + "GCONSALE exceed limit CON not implemented", deferred_logger);
            break;
        }
        case GConSale::MaxProcedure::CON_P: {
            OPM_DEFLOG_THROW(std::runtime_error, "Group " + group.name() + "GCONSALE exceed limit CON_P not implemented", deferred_logger);
            break;
        }
        case GConSale::MaxProcedure::WELL: {
            OPM_DEFLOG_THROW(std::runtime_error, "Group " + group.name() + "GCONSALE exceed limit WELL not implemented", deferred_logger);
            break;
        }
        case GConSale::MaxProcedure::PLUG: {
            OPM_DEFLOG_THROW(std::runtime_error, "Group " + group.name() + "GCONSALE exceed limit PLUG not implemented", deferred_logger);
            break;
        }
        case GConSale::MaxProcedure::MAXR: {
            OPM_DEFLOG_THROW(std::runtime_error, "Group " + group.name() + "GCONSALE exceed limit MAXR not implemented", deferred_logger);
            break;
        }
        case GConSale::MaxProcedure::END: {
            OPM_DEFLOG_THROW(std::runtime_error, "Group " + group.name() + "GCONSALE exceed limit END not implemented", deferred_logger);
            break;
        }
        case GConSale::MaxProcedure::RATE: {
            this->groupState().production_control(group.name(), Group::ProductionCMode::GRAT);
            ss << "Maximum GCONSALE limit violated for " << group.name() << ". The group is switched from ";
            ss << Group::ProductionCMode2String(oldProductionControl) << " to " << Group::ProductionCMode2String(Group::ProductionCMode::GRAT);
            ss << " and limited by the maximum sales rate after consumption and import are considered" ;
            this->groupState().update_grat_sales_target(group.name(), production_target);
            break;
        }
        default:
            throw("Invalid procedure for maximum rate limit selected for group" + group.name());
        }
    }
    if (sales_rate < gconsale.min_sales_rate) {
        const Group::ProductionCMode& currentProductionControl = this->groupState().production_control(group.name());
        if ( currentProductionControl == Group::ProductionCMode::GRAT ) {
            ss << "Group " + group.name() + " has sale rate less then minimum permitted value and is under GRAT control. \n";
            ss << "The GRAT is increased to meet the sales minimum rate. \n";
            this->groupState().update_grat_sales_target(group.name(), production_target);
        //} else if () {//TODO add action for WGASPROD
        //} else if () {//TODO add action for drilling queue
        } else {
            ss << "Group " + group.name() + " has sale rate less then minimum permitted value but cannot increase the group production rate \n";
            ss << "or adjust gas production using WGASPROD or drill new wells to meet the sales target. \n";
            ss << "Note that WGASPROD and drilling queues are not implemented in Flow. No action is taken. \n ";
        }
    }
    if (gconsale.sales_target < 0.0) {
        OPM_DEFLOG_THROW(std::runtime_error, "Group " + group.name() + " has sale rate target less then zero. Not implemented in Flow" , deferred_logger);
    }

    if (!ss.str().empty() && comm_.rank() == 0)
        deferred_logger.info(ss.str());
}

bool
BlackoilWellModelGeneric::
checkGroupHigherConstraints(const Group& group,
                            DeferredLogger& deferred_logger,
                            const int reportStepIdx)
{
    // Set up coefficients for RESV <-> surface rate conversion.
    // Use the pvtRegionIdx from the top cell of the first well.
    // TODO fix this!
    // This is only used for converting RESV rates.
    // What is the proper approach?
    const int fipnum = 0;
    int pvtreg = well_perf_data_.empty() || well_perf_data_[0].empty()
        ? pvt_region_idx_[0]
        : pvt_region_idx_[well_perf_data_[0][0].cell_index];

    bool changed = false;
    if ( comm_.size() > 1)
    {
        // Just like in the sequential case the pvtregion is determined
        // by the first cell of the first well. What is the first well
        // is decided by the order in the Schedule using Well::seqIndex()
        int firstWellIndex = well_perf_data_.empty() ?
            std::numeric_limits<int>::max() : wells_ecl_[0].seqIndex();
        auto regIndexPair = std::make_pair(pvtreg, firstWellIndex);
        std::vector<decltype(regIndexPair)> pairs(comm_.size());
        comm_.allgather(&regIndexPair, 1, pairs.data());
        pvtreg = std::min_element(pairs.begin(), pairs.end(),
                                  [](const auto& p1, const auto& p2){ return p1.second < p2.second;})
            ->first;
    }

    std::vector<double> rates(phase_usage_.num_phases, 0.0);

    bool isField = group.name() == "FIELD";
    if (!isField && group.isInjectionGroup()) {
        // Obtain rates for group.
        std::vector<double> resv_coeff_inj(phase_usage_.num_phases, 0.0);
        calcInjRates(fipnum, pvtreg, resv_coeff_inj);

        for (int phasePos = 0; phasePos < phase_usage_.num_phases; ++phasePos) {
            const double local_current_rate = WellGroupHelpers::sumWellSurfaceRates(group, schedule(), this->wellState(), reportStepIdx, phasePos, /* isInjector */ true);
            // Sum over all processes
            rates[phasePos] = comm_.sum(local_current_rate);
        }
        const Phase all[] = { Phase::WATER, Phase::OIL, Phase::GAS };
        for (Phase phase : all) {
            // Check higher up only if under individual (not FLD) control.
            auto currentControl = this->groupState().injection_control(group.name(), phase);
            if (currentControl != Group::InjectionCMode::FLD && group.injectionGroupControlAvailable(phase)) {
                const Group& parentGroup = schedule().getGroup(group.parent(), reportStepIdx);
                const auto [is_changed, scaling_factor] = WellGroupHelpers::checkGroupConstraintsInj(
                group.name(),
                group.parent(),
                parentGroup,
                this->wellState(),
                this->groupState(),
                reportStepIdx,
                &guideRate_,
                rates.data(),
                phase,
                phase_usage_,
                group.getGroupEfficiencyFactor(),
                schedule(),
                summaryState_,
                resv_coeff_inj,
                deferred_logger);
                if (is_changed) {
                    switched_inj_groups_.insert({ {group.name(), phase}, Group::InjectionCMode2String(Group::InjectionCMode::FLD)});
                    BlackoilWellModelConstraints(*this).
                        actionOnBrokenConstraints(group, Group::InjectionCMode::FLD,
                                                  phase, this->groupState(),
                                                  deferred_logger);
                    WellGroupHelpers::updateWellRatesFromGroupTargetScale(scaling_factor, group, schedule(), reportStepIdx, /* isInjector */ true, this->groupState(), this->wellState());
                    changed = true;
                }
            }
        }
    }

    if (!isField && group.isProductionGroup()) {
        // Obtain rates for group.
        for (int phasePos = 0; phasePos < phase_usage_.num_phases; ++phasePos) {
            const double local_current_rate = WellGroupHelpers::sumWellSurfaceRates(group, schedule(), this->wellState(), reportStepIdx, phasePos, /* isInjector */ false);
            // Sum over all processes
            rates[phasePos] = -comm_.sum(local_current_rate);
        }
        std::vector<double> resv_coeff(phase_usage_.num_phases, 0.0);
        calcRates(fipnum, pvtreg, resv_coeff);
        // Check higher up only if under individual (not FLD) control.
        const Group::ProductionCMode& currentControl = this->groupState().production_control(group.name());
        if (currentControl != Group::ProductionCMode::FLD && group.productionGroupControlAvailable()) {
            const Group& parentGroup = schedule().getGroup(group.parent(), reportStepIdx);
            const auto [is_changed, scaling_factor] = WellGroupHelpers::checkGroupConstraintsProd(
                group.name(),
                group.parent(),
                parentGroup,
                this->wellState(),
                this->groupState(),
                reportStepIdx,
                &guideRate_,
                rates.data(),
                phase_usage_,
                group.getGroupEfficiencyFactor(),
                schedule(),
                summaryState_,
                resv_coeff,
                deferred_logger);
            if (is_changed) {
                switched_prod_groups_.insert({group.name(), Group::ProductionCMode2String(Group::ProductionCMode::FLD)});
                const auto exceed_action = group.productionControls(summaryState_).exceed_action;
                BlackoilWellModelConstraints(*this).
                        actionOnBrokenConstraints(group, exceed_action,
                                                  Group::ProductionCMode::FLD,
                                                  this->groupState(),
                                                  deferred_logger);
                WellGroupHelpers::updateWellRatesFromGroupTargetScale(scaling_factor, group, schedule(), reportStepIdx, /* isInjector */ false, this->groupState(), this->wellState());
                changed = true;
            }
        }
    }

    return changed;
}

void
BlackoilWellModelGeneric::
updateEclWells(const int timeStepIdx,
               const std::unordered_set<std::string>& wells,
               const SummaryState& st)
{
    for (const auto& wname : wells) {
        auto well_iter = std::find_if(this->wells_ecl_.begin(), this->wells_ecl_.end(),
            [&wname] (const auto& well) -> bool
        {
            return well.name() == wname;
        });

        if (well_iter == this->wells_ecl_.end()) {
            continue;
        }

        auto well_index = std::distance(this->wells_ecl_.begin(), well_iter);
        this->wells_ecl_[well_index] = schedule_.getWell(wname, timeStepIdx);

        const auto& well = this->wells_ecl_[well_index];
        auto& pd     = this->well_perf_data_[well_index];
        auto  pdIter = pd.begin();
        for (const auto& conn : well.getConnections()) {
            if (conn.state() != Connection::State::SHUT) {
                pdIter->connection_transmissibility_factor = conn.CF();
                ++pdIter;
            }
        }
        auto& ws = this->wellState().well(well_index);

        ws.updateStatus( well.getStatus() );
        ws.reset_connection_factors(pd);
        ws.update_targets(well, st);
        this->prod_index_calc_[well_index].reInit(well);
    }
}

double
BlackoilWellModelGeneric::
wellPI(const int well_index) const
{
    const auto& pu = this->phase_usage_;
    const auto& pi = this->wellState().well(well_index).productivity_index;

    const auto preferred = this->wells_ecl_[well_index].getPreferredPhase();
    switch (preferred) { // Should really have LIQUID = OIL + WATER here too...
    case Phase::WATER:
        return pu.phase_used[BlackoilPhases::PhaseIndex::Aqua]
            ? pi[pu.phase_pos[BlackoilPhases::PhaseIndex::Aqua]]
            : 0.0;

    case Phase::OIL:
        return pu.phase_used[BlackoilPhases::PhaseIndex::Liquid]
            ? pi[pu.phase_pos[BlackoilPhases::PhaseIndex::Liquid]]
            : 0.0;

    case Phase::GAS:
        return pu.phase_used[BlackoilPhases::PhaseIndex::Vapour]
            ? pi[pu.phase_pos[BlackoilPhases::PhaseIndex::Vapour]]
            : 0.0;

    default:
        throw std::invalid_argument {
            "Unsupported preferred phase " +
            std::to_string(static_cast<int>(preferred))
        };
    }
}

double
BlackoilWellModelGeneric::
wellPI(const std::string& well_name) const
{
    auto well_iter = std::find_if(this->wells_ecl_.begin(), this->wells_ecl_.end(),
        [&well_name](const Well& well)
    {
        return well.name() == well_name;
    });

    if (well_iter == this->wells_ecl_.end()) {
        throw std::logic_error { "Could not find well: " + well_name };
    }

    auto well_index = std::distance(this->wells_ecl_.begin(), well_iter);
    return this->wellPI(well_index);
}

bool
BlackoilWellModelGeneric::
wasDynamicallyShutThisTimeStep(const int well_index) const
{
    return this->closed_this_step_.find(this->wells_ecl_[well_index].name()) !=
           this->closed_this_step_.end();
}

void
BlackoilWellModelGeneric::
updateWsolvent(const Group& group,
               const int reportStepIdx,
               const WellState& wellState)
{
    for (const std::string& groupName : group.groups()) {
        const Group& groupTmp = schedule_.getGroup(groupName, reportStepIdx);
        updateWsolvent(groupTmp, reportStepIdx, wellState);
    }

    if (group.isProductionGroup())
        return;

    auto currentGroupControl = this->groupState().injection_control(group.name(), Phase::GAS);
    if( currentGroupControl == Group::InjectionCMode::REIN ) {
        int gasPos = phase_usage_.phase_pos[BlackoilPhases::Vapour];
        const auto& controls = group.injectionControls(Phase::GAS, summaryState_);
        const Group& groupRein = schedule_.getGroup(controls.reinj_group, reportStepIdx);
        double gasProductionRate = WellGroupHelpers::sumWellSurfaceRates(groupRein, schedule_, wellState, reportStepIdx, gasPos, /*isInjector*/false);
        double solventProductionRate = WellGroupHelpers::sumSolventRates(groupRein, schedule_, wellState, reportStepIdx, /*isInjector*/false);

        solventProductionRate = comm_.sum(solventProductionRate);
        gasProductionRate = comm_.sum(gasProductionRate);

        double wsolvent = 0.0;
        if (std::abs(gasProductionRate) > 1e-6)
            wsolvent = solventProductionRate / gasProductionRate;

        setWsolvent(group, reportStepIdx, wsolvent);
    }
}

void
BlackoilWellModelGeneric::
setWsolvent(const Group& group,
            const int reportStepIdx,
            double wsolvent)
{
    for (const std::string& groupName : group.groups()) {
        const Group& groupTmp = schedule_.getGroup(groupName, reportStepIdx);
        setWsolvent(groupTmp, reportStepIdx, wsolvent);
    }

    for (const std::string& wellName : group.wells()) {
        const auto& wellTmp = schedule_.getWell(wellName, reportStepIdx);
        if (wellTmp.getStatus() == Well::Status::SHUT)
            continue;

        getGenWell(wellName)->setWsolvent(wsolvent);
    }
}

data::GuideRateValue
BlackoilWellModelGeneric::
getGuideRateValues(const Well& well) const
{
    auto grval = data::GuideRateValue{};

    const auto& wname = well.name();
    if (!this->wellState().has(wname)) {
        // No flow rates for 'wname' -- might be before well comes
        // online (e.g., for the initial condition before simulation
        // starts).
        return grval;
    }

    if (!this->guideRate_.has(wname)) {
        // No guiderates exist for 'wname'.
        return grval;
    }

    const auto qs = WellGroupHelpers::
        getWellRateVector(this->wellState(), this->phase_usage_, wname);

    this->getGuideRateValues(qs, well.isInjector(), wname, grval);

    return grval;
}

void
BlackoilWellModelGeneric::
getGuideRateValues(const GuideRate::RateVector& qs,
                   const bool                   is_inj,
                   const std::string&           wgname,
                   data::GuideRateValue&        grval) const
{
    auto getGR = [this, &wgname, &qs](const GuideRateModel::Target t)
    {
        return this->guideRate_.getSI(wgname, t, qs);
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
BlackoilWellModelGeneric::
getGuideRateValues(const Group& group) const
{
    auto grval = data::GuideRateValue{};
    const auto& gname = group.name();

    if (!this->groupState().has_production_rates(gname)) {
        // No flow rates for production group 'gname' -- might be before
        // group comes online (e.g., for the initial condition before
        // simulation starts).
        return grval;
    }

    if (!this->guideRate_.has(gname)) {
        // No guiderates exist for 'gname'.
        return grval;
    }

    const auto qs = WellGroupHelpers::
        getProductionGroupRateVector(this->groupState(), this->phase_usage_, gname);

    const auto is_inj = false; // This procedure only applies to G*PGR.
    this->getGuideRateValues(qs, is_inj, gname, grval);

    return grval;
}

data::GuideRateValue
BlackoilWellModelGeneric::
getGuideRateInjectionGroupValues(const Group& group) const
{
    auto grval = data::GuideRateValue{};

    const auto& gname = group.name();
    if (this->guideRate_.has(gname, Phase::GAS)) {
        grval.set(data::GuideRateValue::Item::Gas,
                  this->guideRate_.getSI(gname, Phase::GAS));
    }

    if (this->guideRate_.has(gname, Phase::WATER)) {
        grval.set(data::GuideRateValue::Item::Water,
                  this->guideRate_.getSI(gname, Phase::WATER));
    }

    return grval;
}

void
BlackoilWellModelGeneric::
assignWellGuideRates(data::Wells& wsrpt,
                     const int    reportStepIdx) const
{
    auto all = std::unordered_map<std::string, data::GuideRateValue>{};
    auto retrieve = std::unordered_map<std::string, RetrieveWellGuideRate>{};

    auto walker = GroupTreeWalker{ this->schedule(), reportStepIdx };

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
                             std::forward_as_tuple(this->guideRate_, group));

        if (inserted) {
            elm->second = elm->second || parent;
        }
    });

    // Populates 'all'.
    walker.wellOp([this, &retrieve, &all](const Well& well)
    {
        const auto& wname = well.name();

        const auto is_nontrivial =
            this->guideRate_.has(wname) || this->guideRate_.hasPotentials(wname);

        if (! (is_nontrivial && this->wellState().has(wname))) {
            all[wname].clear();
            return;
        }

        auto parent_pos = retrieve.find(well.groupName());
        const auto parent = (parent_pos == retrieve.end())
            ? RetrieveWellGuideRate{} // No entry for 'parent'--unexpected.
            : parent_pos->second;

        const auto get_gr = parent
            || RetrieveWellGuideRate{ this->guideRate_, wname };

        const auto qs = WellGroupHelpers::
            getWellRateVector(this->wellState(), this->phase_usage_, wname);

        auto getGR = [this, &wname, &qs](const GuideRateModel::Target t)
        {
            return this->guideRate_.getSI(wname, t, qs);
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

    for (const auto& well : this->wells_ecl_) {
        auto xwPos = wsrpt.find(well.name());
        if (xwPos == wsrpt.end()) { // No well results.  Unexpected.
            continue;
        }

        auto grPos = all.find(well.name());
        if (grPos == all.end()) {
            continue;
        }

        xwPos->second.guide_rates = grPos->second;
    }
}

void
BlackoilWellModelGeneric::
assignShutConnections(data::Wells& wsrpt,
                      const int reportStepIndex) const
{
    auto wellID = 0;

    for (const auto& well : this->wells_ecl_) {
        auto& xwel = wsrpt[well.name()]; // data::Wells is a std::map<>

        xwel.dynamicStatus = this->schedule()
                             .getWell(well.name(), reportStepIndex).getStatus();

        const auto wellIsOpen = xwel.dynamicStatus == Well::Status::OPEN;
        auto skip = [wellIsOpen](const Connection& conn)
        {
            return wellIsOpen && (conn.state() != Connection::State::SHUT);
        };

        if (this->wellTestState().well_is_closed(well.name()) &&
            !this->wasDynamicallyShutThisTimeStep(wellID))
        {
            xwel.dynamicStatus = well.getAutomaticShutIn()
                                 ? Well::Status::SHUT : Well::Status::STOP;
        }

        auto& xcon = xwel.connections;
        for (const auto& conn : well.getConnections()) {
            if (skip(conn)) {
                continue;
            }

            auto& xc = xcon.emplace_back();
            xc.index = conn.global_index();
            xc.pressure = xc.reservoir_rate = 0.0;

            xc.effective_Kh = conn.Kh();
            xc.trans_factor = conn.CF();
        }

        ++wellID;
    }
}

std::unordered_map<std::string, data::GroupGuideRates>
BlackoilWellModelGeneric::
calculateAllGroupGuiderates(const int reportStepIdx) const
{
    auto gr = std::unordered_map<std::string, data::GroupGuideRates>{};

    auto walker = GroupTreeWalker{ this->schedule(), reportStepIdx };

    // Populates 'gr'.
    walker.groupOp([this, &gr](const Group& group)
    {
        const auto& gname = group.name();

        if (gname == "FIELD") { return; }

        if (this->guideRate_.has(gname)) {
            gr[gname].production = this->getGuideRateValues(group);
        }

        if (this->guideRate_.has(gname, Phase::WATER) ||
            this->guideRate_.has(gname, Phase::GAS))
        {
            gr[gname].injection =
                this->getGuideRateInjectionGroupValues(group);
        }

        const auto parent = group.parent();
        if (parent == "FIELD") { return; }

        gr[parent].injection  += gr[gname].injection;
        gr[parent].production += gr[gname].production;
    });

    // Populates 'gr'.
    walker.wellOp([this, &gr](const Well& well)
    {
        if (! (this->guideRate_.has(well.name()) ||
               this->guideRate_.hasPotentials(well.name())))
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

void
BlackoilWellModelGeneric::
assignGroupControl(const Group& group,
                   data::GroupData& gdata) const
{
    const auto& gname     = group.name();
    const auto  grup_type = group.getGroupType();
    auto&       cgc       = gdata.currentControl;

    cgc.currentProdConstraint = Group::ProductionCMode::NONE;

    cgc.currentGasInjectionConstraint =
    cgc.currentWaterInjectionConstraint = Group::InjectionCMode::NONE;

    if (this->groupState().has_production_control(gname)) {
        cgc.currentProdConstraint = this->groupState().production_control(gname);
    }

    if ((grup_type == ::Opm::Group::GroupType::INJECTION) ||
        (grup_type == ::Opm::Group::GroupType::MIXED))
    {
        if (this->groupState().has_injection_control(gname, Phase::WATER)) {
            cgc.currentWaterInjectionConstraint = this->groupState().injection_control(gname, Phase::WATER);
        }

        if (this->groupState().has_injection_control(gname, Phase::GAS)) {
            cgc.currentGasInjectionConstraint = this->groupState().injection_control(gname, Phase::GAS);
        }
    }
}

void
BlackoilWellModelGeneric::
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

    if (this->guideRate_.has(group.name())) {
        prod = xgr.production;
    }

    if (this->guideRate_.has(group.name(), Phase::WATER) ||
        this->guideRate_.has(group.name(), Phase::GAS))
    {
        inj = xgr.injection;
    }
}

void
BlackoilWellModelGeneric::
assignGroupValues(const int                               reportStepIdx,
                  std::map<std::string, data::GroupData>& gvalues) const
{
    const auto groupGuideRates =
        this->calculateAllGroupGuiderates(reportStepIdx);

    for (const auto& gname : schedule_.groupNames(reportStepIdx)) {
        const auto& grup = schedule_.getGroup(gname, reportStepIdx);

        auto& gdata = gvalues[gname];
        this->assignGroupControl(grup, gdata);
        this->assignGroupGuideRates(grup, groupGuideRates, gdata);
    }
}

void
BlackoilWellModelGeneric::
assignNodeValues(std::map<std::string, data::NodeData>& nodevalues) const
{
    nodevalues.clear();
    for (const auto& [node, pressure] : node_pressures_) {
        nodevalues.emplace(node, data::NodeData{pressure});
    }
}

data::GroupAndNetworkValues
BlackoilWellModelGeneric::
groupAndNetworkData(const int reportStepIdx) const
{
    auto grp_nwrk_values = data::GroupAndNetworkValues{};

    this->assignGroupValues(reportStepIdx, grp_nwrk_values.groupData);
    this->assignNodeValues(grp_nwrk_values.nodeData);

    return grp_nwrk_values;
}

void
BlackoilWellModelGeneric::
updateAndCommunicateGroupData(const int reportStepIdx,
                              const int iterationIdx)
{
    const Group& fieldGroup = schedule().getGroup("FIELD", reportStepIdx);
    const int nupcol = schedule()[reportStepIdx].nupcol();

    // This builds some necessary lookup structures, so it must be called
    // before we copy to well_state_nupcol_.
    this->wellState().updateGlobalIsGrup(comm_);

    if (iterationIdx < nupcol) {
        this->updateNupcolWGState();
    }

    auto& well_state = this->wellState();
    const auto& well_state_nupcol = this->nupcolWellState();
    // the group target reduction rates needs to be update since wells may have switched to/from GRUP control
    // The group target reduction does not honor NUPCOL.
    std::vector<double> groupTargetReduction(numPhases(), 0.0);
    WellGroupHelpers::updateGroupTargetReduction(fieldGroup, schedule(), reportStepIdx, /*isInjector*/ false, phase_usage_, guideRate_, well_state, this->groupState(), groupTargetReduction);
    std::vector<double> groupTargetReductionInj(numPhases(), 0.0);
    WellGroupHelpers::updateGroupTargetReduction(fieldGroup, schedule(), reportStepIdx, /*isInjector*/ true, phase_usage_, guideRate_, well_state, this->groupState(), groupTargetReductionInj);

    WellGroupHelpers::updateREINForGroups(fieldGroup, schedule(), reportStepIdx, phase_usage_, summaryState_, well_state_nupcol, this->groupState());
    WellGroupHelpers::updateVREPForGroups(fieldGroup, schedule(), reportStepIdx, well_state_nupcol, this->groupState());

    WellGroupHelpers::updateReservoirRatesInjectionGroups(fieldGroup, schedule(), reportStepIdx, well_state_nupcol, this->groupState());
    WellGroupHelpers::updateSurfaceRatesInjectionGroups(fieldGroup, schedule(), reportStepIdx, well_state_nupcol, this->groupState());

    WellGroupHelpers::updateGroupProductionRates(fieldGroup, schedule(), reportStepIdx, well_state_nupcol, this->groupState());

    // We use the rates from the previous time-step to reduce oscillations
    WellGroupHelpers::updateWellRates(fieldGroup, schedule(), reportStepIdx, this->prevWellState(), well_state);

    // Set ALQ for off-process wells to zero
    for (const auto& wname : schedule().wellNames(reportStepIdx)) {
        const bool is_producer = schedule().getWell(wname, reportStepIdx).isProducer();
        const bool not_on_this_process = !well_state.has(wname);
        if (is_producer && not_on_this_process) {
            well_state.setALQ(wname, 0.0);
        }
    }

    well_state.communicateGroupRates(comm_);
    this->groupState().communicate_rates(comm_);
}

bool
BlackoilWellModelGeneric::
hasTHPConstraints() const
{
    return BlackoilWellModelConstraints(*this).hasTHPConstraints();
}

bool
BlackoilWellModelGeneric::
forceShutWellByName(const std::string& wellname,
                    const double simulation_time)
{
    // Only add the well to the closed list on the
    // process that owns it.
    int well_was_shut = 0;
    for (const auto& well : well_container_generic_) {
        if (well->name() == wellname) {
            wellTestState().close_well(wellname, WellTestConfig::Reason::PHYSICAL, simulation_time);
            well_was_shut = 1;
            break;
        }
    }

    // Communicate across processes if a well was shut.
    well_was_shut = comm_.max(well_was_shut);

    // the wellTesteState is updated between timesteps and we also need to update the privous WGstate
    if(well_was_shut)
        this->commitWGState();

    // Only log a message on the output rank.
    if (terminal_output_ && well_was_shut) {
        const std::string msg = "Well " + wellname
            + " will be shut because it cannot get converged.";
        OpmLog::info(msg);
    }

    return (well_was_shut == 1);
}

void
BlackoilWellModelGeneric::
inferLocalShutWells()
{
    this->local_shut_wells_.clear();

    const auto nw = this->numLocalWells();

    auto used = std::vector<bool>(nw, false);
    for (const auto& wellPtr : this->well_container_generic_) {
        used[wellPtr->indexOfWell()] = true;
    }

    for (auto wellID = 0; wellID < nw; ++wellID) {
        if (! used[wellID]) {
            this->local_shut_wells_.push_back(wellID);
        }
    }
}

std::pair<bool, double>
BlackoilWellModelGeneric::
updateNetworkPressures(const int reportStepIdx)
{
    // Get the network and return if inactive.
    const auto& network = schedule()[reportStepIdx].network();
    if (!network.active()) {
        return { false, 0.0 };
    }
    node_pressures_ = WellGroupHelpers::computeNetworkPressures(network,
                                                                this->wellState(),
                                                                this->groupState(),
                                                                *(vfp_properties_->getProd()),
                                                                schedule(),
                                                                reportStepIdx);

    // Set the thp limits of wells
    bool active_limit_change = false;
    double network_imbalance = 0.0;
    for (auto& well : well_container_generic_) {
        // Producers only, since we so far only support the
        // "extended" network model (properties defined by
        // BRANPROP and NODEPROP) which only applies to producers.
        if (well->isProducer()) {
            const auto it = node_pressures_.find(well->wellEcl().groupName());
            if (it != node_pressures_.end()) {
                // The well belongs to a group with has a network pressure constraint,
                // set the dynamic THP constraint of the well accordingly.
                const double new_limit = it->second;
                well->setDynamicThpLimit(new_limit);
                const SingleWellState& ws = this->wellState()[well->indexOfWell()];
                const bool thp_is_limit = ws.production_cmode == Well::ProducerCMode::THP;
                const bool will_switch_to_thp = ws.thp < new_limit;
                if (thp_is_limit || will_switch_to_thp) {
                    active_limit_change = true;
                    network_imbalance = std::max(network_imbalance, std::fabs(new_limit - ws.thp));
                }
            }
        }
    }
    return { active_limit_change, network_imbalance };
}

void
BlackoilWellModelGeneric::
calculateEfficiencyFactors(const int reportStepIdx)
{
    for (auto& well : well_container_generic_) {
        const Well& wellEcl = well->wellEcl();
        double well_efficiency_factor = wellEcl.getEfficiencyFactor();
        WellGroupHelpers::accumulateGroupEfficiencyFactor(schedule().getGroup(wellEcl.groupName(), reportStepIdx), schedule(), reportStepIdx, well_efficiency_factor);
        well->setWellEfficiencyFactor(well_efficiency_factor);
    }
}

WellInterfaceGeneric*
BlackoilWellModelGeneric::
getGenWell(const std::string& well_name)
{
    // finding the iterator of the well in wells_ecl
    auto well = std::find_if(well_container_generic_.begin(),
                             well_container_generic_.end(),
                                 [&well_name](const WellInterfaceGeneric* elem)->bool {
                                     return elem->name() == well_name;
                                 });

    assert(well != well_container_generic_.end());

    return *well;
}

void
BlackoilWellModelGeneric::
setRepRadiusPerfLength()
{
    for (const auto& well : well_container_generic_) {
        well->setRepRadiusPerfLength();
    }
}

void
BlackoilWellModelGeneric::
gliftDebug(const std::string& msg,
           DeferredLogger& deferred_logger) const
{
    if (this->glift_debug && this->terminal_output_) {
        const std::string message = fmt::format(
            "  GLIFT (DEBUG) : BlackoilWellModel : {}", msg);
        deferred_logger.info(message);
    }
}


void
BlackoilWellModelGeneric::
gliftDebugShowALQ(DeferredLogger& deferred_logger)
{
    for (auto& well : this->well_container_generic_) {
        if (well->isProducer()) {
            auto alq = this->wellState().getALQ(well->name());
            const std::string msg = fmt::format("ALQ_REPORT : {} : {}",
                                                well->name(), alq);
            gliftDebug(msg, deferred_logger);
        }
    }
}

// If a group has any production rate constraints, and/or a limit
// on its total rate of lift gas supply,  allocate lift gas
// preferentially to the wells that gain the most benefit from
// it. Lift gas increments are allocated in turn to the well that
// currently has the largest weighted incremental gradient. The
// procedure takes account of any limits on the group production
// rate or lift gas supply applied to any level of group.
void
BlackoilWellModelGeneric::
gasLiftOptimizationStage2(DeferredLogger& deferred_logger,
                          GLiftProdWells& prod_wells,
                          GLiftOptWells& glift_wells,
                          GLiftWellStateMap& glift_well_state_map,
                          const int episodeIndex)
{
    GasLiftStage2 glift {episodeIndex,
                         comm_,
                         schedule_,
                         summaryState_,
                         deferred_logger,
                         this->wellState(),
                         prod_wells,
                         glift_wells,
                         glift_well_state_map,
                         this->glift_debug
    };
    glift.runOptimize();
}

void
BlackoilWellModelGeneric::
updateWellPotentials(const int reportStepIdx,
                     const bool onlyAfterEvent,
                     const SummaryConfig& summaryConfig,
                     DeferredLogger& deferred_logger)
{
    auto well_state_copy = this->wellState();

    const bool write_restart_file = schedule().write_rst_file(reportStepIdx);
    auto exc_type = ExceptionType::NONE;
    std::string exc_msg;
    size_t widx = 0;
    for (const auto& well : well_container_generic_) {
        const bool needed_for_summary =
                ((summaryConfig.hasSummaryKey( "WWPI:" + well->name()) ||
                  summaryConfig.hasSummaryKey( "WOPI:" + well->name()) ||
                  summaryConfig.hasSummaryKey( "WGPI:" + well->name())) && well->isInjector()) ||
                ((summaryConfig.hasKeyword( "GWPI") ||
                  summaryConfig.hasKeyword( "GOPI") ||
                  summaryConfig.hasKeyword( "GGPI")) && well->isInjector()) ||
                ((summaryConfig.hasKeyword( "FWPI") ||
                  summaryConfig.hasKeyword( "FOPI") ||
                  summaryConfig.hasKeyword( "FGPI")) && well->isInjector()) ||
                ((summaryConfig.hasSummaryKey( "WWPP:" + well->name()) ||
                  summaryConfig.hasSummaryKey( "WOPP:" + well->name()) ||
                  summaryConfig.hasSummaryKey( "WGPP:" + well->name())) && well->isProducer()) ||
                ((summaryConfig.hasKeyword( "GWPP") ||
                  summaryConfig.hasKeyword( "GOPP") ||
                  summaryConfig.hasKeyword( "GGPP")) && well->isProducer()) ||
                ((summaryConfig.hasKeyword( "FWPP") ||
                  summaryConfig.hasKeyword( "FOPP") ||
                  summaryConfig.hasKeyword( "FGPP")) && well->isProducer());

        // At the moment, the following events are considered
        // for potentials update
        const uint64_t effective_events_mask = ScheduleEvents::WELL_STATUS_CHANGE
                                             + ScheduleEvents::COMPLETION_CHANGE
                                             + ScheduleEvents::WELL_PRODUCTIVITY_INDEX
                                             + ScheduleEvents::WELL_WELSPECS_UPDATE
                                             + ScheduleEvents::WELLGROUP_EFFICIENCY_UPDATE
                                             + ScheduleEvents::NEW_WELL
                                             + ScheduleEvents::PRODUCTION_UPDATE
                                             + ScheduleEvents::INJECTION_UPDATE;
        const auto& events = schedule()[reportStepIdx].wellgroup_events();
        const bool event = events.hasEvent(well->name(), ScheduleEvents::ACTIONX_WELL_EVENT) ||
                           (report_step_starts_ && events.hasEvent(well->name(), effective_events_mask));
        const bool needPotentialsForGuideRates = well->underPredictionMode() && (!onlyAfterEvent || event);
        const bool needPotentialsForOutput = !onlyAfterEvent && (needed_for_summary || write_restart_file);
        const bool compute_potential = needPotentialsForOutput || needPotentialsForGuideRates;
        if (compute_potential)
        {
            this->computePotentials(widx, well_state_copy, exc_msg, exc_type, deferred_logger);
        }
        ++widx;
    }
    logAndCheckForExceptionsAndThrow(deferred_logger, exc_type,
                                     "computeWellPotentials() failed: " + exc_msg,
                                     terminal_output_, comm_);

}

void
BlackoilWellModelGeneric::
runWellPIScaling(const int timeStepIdx,
                 DeferredLogger& local_deferredLogger)
{
    if (this->last_run_wellpi_.has_value() && (*this->last_run_wellpi_ == timeStepIdx)) {
        // We've already run WELPI scaling for this report step.  Most
        // common for the very first report step.  Don't redo WELPI scaling.
        return;
    }

    auto hasWellPIEvent = [this, timeStepIdx](const int well_index) -> bool
    {
        return this->schedule()[timeStepIdx].wellgroup_events()
            .hasEvent(this->wells_ecl_[well_index].name(),
                      ScheduleEvents::Events::WELL_PRODUCTIVITY_INDEX);
    };

    auto updateEclWell = [this, timeStepIdx](const int well_index) -> void
    {
        const auto& schedule = this->schedule();
        const auto& wname = this->wells_ecl_[well_index].name();
        this->wells_ecl_[well_index] = schedule.getWell(wname, timeStepIdx);

        const auto& well = this->wells_ecl_[well_index];
        auto& pd     = this->well_perf_data_[well_index];
        auto  pdIter = pd.begin();
        for (const auto& conn : well.getConnections()) {
            if (conn.state() != Connection::State::SHUT) {
                pdIter->connection_transmissibility_factor = conn.CF();
                ++pdIter;
            }
        }
        auto& ws = this->wellState().well(well_index);
        ws.reset_connection_factors(pd);
        this->prod_index_calc_[well_index].reInit(well);
    };


    auto rescaleWellPI =
        [this, timeStepIdx](const int    well_index,
                            const double newWellPI) -> void
    {
        const auto& wname = this->wells_ecl_[well_index].name();

        schedule_.applyWellProdIndexScaling(wname, timeStepIdx, newWellPI);
    };

    // Minimal well setup to compute PI/II values
    {
        auto saved_previous_wgstate = this->prevWGState();
        this->commitWGState();

        this->createWellContainer(timeStepIdx);
        this->inferLocalShutWells();

        this->initWellContainer(timeStepIdx);

        this->calculateProductivityIndexValues(local_deferredLogger);
        this->calculateProductivityIndexValuesShutWells(timeStepIdx, local_deferredLogger);

        this->commitWGState(std::move(saved_previous_wgstate));
    }

    const auto nw = this->numLocalWells();
    for (auto wellID = 0*nw; wellID < nw; ++wellID) {
        if (hasWellPIEvent(wellID)) {
            rescaleWellPI(wellID, this->wellPI(wellID));
            updateEclWell(wellID);
        }
    }

    this->last_run_wellpi_ = timeStepIdx;
}


bool
BlackoilWellModelGeneric::
guideRateUpdateIsNeeded(const int reportStepIdx) const {
    auto need_update =
    std::any_of(this->well_container_generic_.begin(),
                this->well_container_generic_.end(),
    [](const WellInterfaceGeneric* well)
    {
        return well->changedToOpenThisStep();
    });
    if (!need_update && this->report_step_starts_) {
        const auto& events = this->schedule()[reportStepIdx].wellgroup_events();
        constexpr auto effective_events_mask = ScheduleEvents::WELL_STATUS_CHANGE
            + ScheduleEvents::INJECTION_TYPE_CHANGED
            + ScheduleEvents::WELL_SWITCHED_INJECTOR_PRODUCER
            + ScheduleEvents::NEW_WELL;

        need_update = std::any_of(this->well_container_generic_.begin(),
                              this->well_container_generic_.end(),
            [&events](const WellInterfaceGeneric* well)
        {
            return events.hasEvent(well->name(), effective_events_mask);
        });
    }
    return this->comm_.max(static_cast<int>(need_update));
}


}
