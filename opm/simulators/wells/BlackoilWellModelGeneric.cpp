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

#include <opm/common/ErrorMacros.hpp>
#include <opm/common/TimingMacros.hpp>
#include <opm/output/data/GuideRateValue.hpp>
#include <opm/output/data/Groups.hpp>
#include <opm/output/data/Wells.hpp>
#include <opm/output/eclipse/RestartValue.hpp>

#include <opm/input/eclipse/EclipseState/EclipseState.hpp>
#include <opm/input/eclipse/EclipseState/SummaryConfig/SummaryConfig.hpp>

#include <opm/input/eclipse/Schedule/Action/SimulatorUpdate.hpp>
#include <opm/input/eclipse/Schedule/Group/GConSale.hpp>
#include <opm/input/eclipse/Schedule/Group/GroupEconProductionLimits.hpp>
#include <opm/input/eclipse/Schedule/Group/GConSump.hpp>
#include <opm/input/eclipse/Schedule/Group/GuideRateConfig.hpp>
#include <opm/input/eclipse/Schedule/Group/GuideRate.hpp>
#include <opm/input/eclipse/Schedule/Network/Balance.hpp>
#include <opm/input/eclipse/Schedule/Network/ExtNetwork.hpp>
#include <opm/input/eclipse/Schedule/Schedule.hpp>
#include <opm/input/eclipse/Schedule/ScheduleTypes.hpp>
#include <opm/input/eclipse/Schedule/Well/NameOrder.hpp>
#include <opm/input/eclipse/Schedule/Well/WellConnections.hpp>
#include <opm/input/eclipse/Schedule/Well/WellTestConfig.hpp>

#include <opm/input/eclipse/Units/Units.hpp>

#include <opm/models/utils/parametersystem.hpp>

#include <opm/material/fluidsystems/BlackOilDefaultFluidSystemIndices.hpp>

#include <opm/simulators/utils/DeferredLogger.hpp>
#include <opm/simulators/wells/BlackoilWellModelConstraints.hpp>
#include <opm/simulators/wells/BlackoilWellModelGasLift.hpp>
#include <opm/simulators/wells/BlackoilWellModelGuideRates.hpp>
#include <opm/simulators/wells/BlackoilWellModelRestart.hpp>
#include <opm/simulators/wells/GasLiftStage2.hpp>
#include <opm/simulators/wells/GroupEconomicLimitsChecker.hpp>
#include <opm/simulators/wells/ParallelWBPCalculation.hpp>
#include <opm/simulators/wells/VFPProperties.hpp>
#include <opm/simulators/wells/WellFilterCake.hpp>
#include <opm/simulators/wells/WellGroupHelper.hpp>
#include <opm/simulators/wells/WellInterfaceGeneric.hpp>
#include <opm/simulators/wells/WellState.hpp>

#if HAVE_MPI
#include <opm/simulators/utils/MPISerializer.hpp>
#endif

#include <algorithm>
#include <cassert>
#include <functional>
#include <iterator>
#include <sstream>
#include <stdexcept>
#include <tuple>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#include <fmt/format.h>

namespace Opm {

template<typename Scalar, typename IndexTraits>
BlackoilWellModelGeneric<Scalar, IndexTraits>::
BlackoilWellModelGeneric(Schedule& schedule,
                         BlackoilWellModelGasLiftGeneric<Scalar, IndexTraits>& gaslift,
                         const SummaryState& summaryState,
                         const EclipseState& eclState,
                         const PhaseUsageInfo<IndexTraits>& pu,
                         const Parallel::Communication& comm)
    : schedule_(schedule)
    , summaryState_(summaryState)
    , eclState_(eclState)
    , comm_(comm)
    , gen_gaslift_(gaslift)
    , wbp_(*this)
    , phase_usage_info_(pu)
    , terminal_output_(comm_.rank() == 0 &&
                       Parameters::Get<Parameters::EnableTerminalOutput>())
    , guideRate_(schedule)
    , active_wgstate_(pu)
    , last_valid_wgstate_(pu)
    , nupcol_wgstate_(pu)
    , wg_helper_(this->wellState(), this->groupState(), this->schedule(), summaryState, guideRate_, pu)
{

    const auto numProcs = comm_.size();
    this->not_on_process_ = [this, numProcs](const std::string& well) {
        if (numProcs == decltype(numProcs){1}) {
            return false;
        }

        // Recall: false indicates NOT active!
        const auto value = std::make_pair(well, true);
        auto candidate = std::lower_bound(this->parallel_well_info_.begin(),
                                          this->parallel_well_info_.end(),
                                          value);

        return (candidate == this->parallel_well_info_.end())
            || (*candidate != value);
    };

    const auto& node_pressures = eclState.getRestartNetworkPressures();
    if (node_pressures.has_value()) {
        if constexpr (std::is_same_v<Scalar,double>) {
            this->node_pressures_ = node_pressures.value();
        } else {
            for (const auto& it : node_pressures.value()) {
                this->node_pressures_[it.first] = it.second;
            }
        }
    }
}

template<typename Scalar, typename IndexTraits>
int BlackoilWellModelGeneric<Scalar, IndexTraits>::
numLocalWells() const
{
    return wells_ecl_.size();
}

template<typename Scalar, typename IndexTraits>
int BlackoilWellModelGeneric<Scalar, IndexTraits>::
numPhases() const
{
    return phase_usage_info_.numActivePhases();
}

template<typename Scalar, typename IndexTraits>
bool BlackoilWellModelGeneric<Scalar, IndexTraits>::
hasLocalWell(const std::string& wname) const
{
    return std::any_of(this->wells_ecl_.begin(),
                       this->wells_ecl_.end(),
        [&wname](const Well& well)
    {
        return well.name() == wname;
    });
}

template<typename Scalar, typename IndexTraits>
bool
BlackoilWellModelGeneric<Scalar, IndexTraits>::
hasOpenLocalWell(const std::string& wname) const
{
    return std::any_of(well_container_generic_.begin(),
                       well_container_generic_.end(),
        [&wname](const auto* elem) -> bool
    {
        return elem->name() == wname;
    });
}

template<typename Scalar, typename IndexTraits>
bool BlackoilWellModelGeneric<Scalar, IndexTraits>::
wellsActive() const
{
    return wells_active_;
}

template<typename Scalar, typename IndexTraits>
bool BlackoilWellModelGeneric<Scalar, IndexTraits>::
networkActive() const
{
    return network_active_;
}

template<typename Scalar, typename IndexTraits>
bool BlackoilWellModelGeneric<Scalar, IndexTraits>::
anyMSWellOpenLocal() const
{
    return std::any_of(wells_ecl_.begin(), wells_ecl_.end(),
                       [](const auto& well) { return well.isMultiSegment(); });
}

template<typename Scalar, typename IndexTraits>
const Well& BlackoilWellModelGeneric<Scalar, IndexTraits>::
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

template<typename Scalar, typename IndexTraits>
void BlackoilWellModelGeneric<Scalar, IndexTraits>::
initFromRestartFile(const RestartValue& restartValues,
                    std::unique_ptr<WellTestState> wtestState,
                    const std::size_t numCells,
                    bool handle_ms_well,
                    bool enable_distributed_wells)
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

    handle_ms_well &= anyMSWellOpenLocal();
    // Resize for restart step
    this->wellState().resize(this->wells_ecl_, this->local_parallel_well_info_,
                             this->schedule(), handle_ms_well, numCells,
                             this->well_perf_data_, this->summaryState_, enable_distributed_wells);

    BlackoilWellModelRestart(*this).
        loadRestartData(restartValues.wells,
                        restartValues.grp_nwrk,
                        handle_ms_well,
                        this->wellState(),
                        this->groupState());

    if (config.has_model()) {
        BlackoilWellModelRestart(*this).
            loadRestartGuideRates(report_step,
                                  config.model().target(),
                                  restartValues.wells,
                                  this->guideRate_);

        BlackoilWellModelRestart(*this).
            loadRestartGuideRates(report_step,
                                  config,
                                  restartValues.grp_nwrk.groupData,
                                  this->guideRate_);

        this->guideRate_.updateGuideRateExpiration(this->schedule().seconds(report_step), report_step);
    }

    this->active_wgstate_.wtest_state(std::move(wtestState));
    this->commitWGState();
    initial_step_ = false;
}

template<typename Scalar, typename IndexTraits>
void BlackoilWellModelGeneric<Scalar, IndexTraits>::
prepareDeserialize(int report_step, const std::size_t numCells, bool handle_ms_well, bool enable_distributed_wells)
{
    // wells_ecl_ should only contain wells on this processor.
    wells_ecl_ = getLocalWells(report_step);
    this->local_parallel_well_info_ = createLocalParallelWellInfo(wells_ecl_);

    this->initializeWellProdIndCalculators();
    initializeWellPerfData();

    if (! this->wells_ecl_.empty()) {
        handle_ms_well &= anyMSWellOpenLocal();
        this->wellState().resize(this->wells_ecl_, this->local_parallel_well_info_,
                                 this->schedule(), handle_ms_well, numCells,
                                 this->well_perf_data_, this->summaryState_, enable_distributed_wells);

    }
    this->wellState().clearWellRates();
    this->commitWGState();
    this->updateNupcolWGState();
}

template<typename Scalar, typename IndexTraits>
std::vector<Well> BlackoilWellModelGeneric<Scalar, IndexTraits>::
getLocalWells(const int timeStepIdx) const
{
    auto w = std::vector<Well>{};

    auto wnames = this->schedule().wellNames(timeStepIdx);
    wnames.erase(std::remove_if(wnames.begin(), wnames.end(),
                                this->not_on_process_),
                 wnames.end());

    w.reserve(wnames.size());

    std::transform(wnames.begin(), wnames.end(),
                   std::back_inserter(w),
                   [&st = this->schedule()[timeStepIdx]]
                   (const std::string& wname)
                   { return st.wells(wname); });

    return w;
}

template<typename Scalar, typename IndexTraits>
std::vector<std::reference_wrapper<ParallelWellInfo<Scalar>>>
BlackoilWellModelGeneric<Scalar, IndexTraits>::
createLocalParallelWellInfo(const std::vector<Well>& wells)
{
    std::vector<std::reference_wrapper<ParallelWellInfo<Scalar>>> local_parallel_well_info;
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

template<typename Scalar, typename IndexTraits>
void BlackoilWellModelGeneric<Scalar, IndexTraits>::
initializeWellProdIndCalculators()
{
    this->prod_index_calc_.clear();
    this->prod_index_calc_.reserve(this->wells_ecl_.size());
    for (const auto& well : this->wells_ecl_) {
        this->prod_index_calc_.emplace_back(well);
    }
}

template<typename Scalar, typename IndexTraits>
void BlackoilWellModelGeneric<Scalar, IndexTraits>::
initializeWellPerfData()
{
    well_perf_data_.resize(wells_ecl_.size());

    this->conn_idx_map_.clear();
    this->conn_idx_map_.reserve(wells_ecl_.size());

    int well_index = 0;

    for (const auto& well : wells_ecl_) {
        int connection_index = 0;

        // INVALID_ECL_INDEX marks no above perf available
        int connection_index_above = ParallelWellInfo<Scalar>::INVALID_ECL_INDEX;

        well_perf_data_[well_index].clear();
        well_perf_data_[well_index].reserve(well.getConnections().size());

        auto& connIdxMap = this->conn_idx_map_
            .emplace_back(well.getConnections().size());

        CheckDistributedWellConnections checker {
            well, this->local_parallel_well_info_[well_index].get()
        };

        bool hasFirstConnection = false;
        bool firstOpenConnection = true;

        auto& parallelWellInfo = this->local_parallel_well_info_[well_index].get();
        parallelWellInfo.beginReset();


        for (const auto& connection : well.getConnections()) {

            const int active_index = well.is_lgr_well()
                ? compressedIndexForInteriorLGR(well.get_lgr_well_tag().value(), connection)
                : this->compressedIndexForInterior(connection.global_index());
            const auto connIsOpen =
                connection.state() == Connection::State::OPEN;

            if (active_index >= 0) {
                connIdxMap.addActiveConnection(connection_index, connIsOpen);
            }

            if ((connIsOpen && (active_index >= 0)) || !connIsOpen) {
                checker.connectionFound(connection_index);
            }

            if (connIsOpen) {
                if (active_index >= 0) {
                    if (firstOpenConnection) {
                        hasFirstConnection = true;
                    }

                    auto& pd = well_perf_data_[well_index].emplace_back();

                    pd.cell_index = active_index;
                    pd.connection_transmissibility_factor = connection.CF();
                    pd.satnum_id = connection.satTableId();
                    pd.ecl_index = connection_index;

                    parallelWellInfo.pushBackEclIndex(connection_index_above,
                                                      connection_index);
                }

                firstOpenConnection = false;

                // Next time this index is the one above as each open
                // connection is stored somewhere.
                connection_index_above = connection_index;
            }
            else if (connection.state() != Connection::State::SHUT) {
                OPM_THROW(std::runtime_error,
                          fmt::format("Connection state '{}' not handled",
                                      Connection::State2String(connection.state())));
            }

            // Note: we rely on the connections being filtered!  I.e., there
            // are only connections to active cells in the global grid.
            ++connection_index;
        }

        parallelWellInfo.endReset();

        checker.checkAllConnectionsFound();

        parallelWellInfo.communicateFirstPerforation(hasFirstConnection);

        ++well_index;
    }
}

template<typename Scalar, typename IndexTraits>
void BlackoilWellModelGeneric<Scalar, IndexTraits>::
checkGEconLimits(
        const Group& group,
        const double simulation_time,
        const int report_step_idx,
        DeferredLogger& deferred_logger)
{
     // call recursively down the group hiearchy
    for (const std::string& group_name : group.groups()) {
        checkGEconLimits( schedule().getGroup(group_name, report_step_idx),
                          simulation_time, report_step_idx, deferred_logger);
    }

    // check if gecon is used for this group
    if (!schedule()[report_step_idx].gecon().has_group(group.name())) {
        return;
    }

    GroupEconomicLimitsChecker<Scalar, IndexTraits> checker {
        *this, wellTestState(), group, simulation_time, report_step_idx, deferred_logger
    };
    if (checker.minOilRate() || checker.minGasRate()) {
        checker.closeWells();
    }
    else if (checker.waterCut() || checker.GOR() || checker.WGR()) {
        checker.doWorkOver();
    }
    if (checker.endRun() && (checker.numProducersOpenInitially() >= 1)
                             && (checker.numProducersOpen() == 0))
    {
        checker.activateEndRun();
    }
}

template<typename Scalar, typename IndexTraits>
void BlackoilWellModelGeneric<Scalar, IndexTraits>::
checkGconsaleLimits(const Group& group,
                    WellState<Scalar, IndexTraits>& well_state,
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

    std::string ss;

    const auto& gconsale = schedule()[reportStepIdx].gconsale().get(group.name(), summaryState_);
    const Group::ProductionCMode& oldProductionControl = this->groupState().production_control(group.name());

    const int gasPos = phase_usage_info_.canonicalToActivePhaseIdx(IndexTraits::gasPhaseIdx);
    Scalar production_rate = this->wgHelper().sumWellSurfaceRates(
        group, gasPos, /*is_injector=*/false
    );
    Scalar injection_rate = this->wgHelper().sumWellSurfaceRates(
        group, gasPos, /*is_injector=*/true
    );
    // sum over all nodes
    injection_rate = comm_.sum(injection_rate);
    production_rate = comm_.sum(production_rate);

    Scalar sales_rate = production_rate - injection_rate;
    Scalar production_target = gconsale.sales_target + injection_rate;

    // add import rate and subtract consumption rate for group for gas
    if (phase_usage_info_.phaseIsActive(IndexTraits::gasPhaseIdx)) {
        const auto& [consumption_rate, import_rate] = this->groupState().gconsump_rates(group.name());
        sales_rate += import_rate;
        sales_rate -= consumption_rate;
        production_target -= import_rate;
        production_target += consumption_rate;
    }

    if (sales_rate > gconsale.max_sales_rate) {
        switch (gconsale.max_proc) {
        case GConSale::MaxProcedure::NONE: {
            if (oldProductionControl != Group::ProductionCMode::GRAT && oldProductionControl != Group::ProductionCMode::NONE) {
                ss = fmt::format("Group sales exceed maximum limit, but the action is NONE for {}. Nothing happens",
                                 group.name());
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
            ss = fmt::format("Maximum GCONSALE limit violated for {}. "
                             "The group is switched from {} to {} "
                             "and limited by the maximum sales rate after "
                             "consumption and import are considered",
                             group.name(),
                             Group::ProductionCMode2String(oldProductionControl),
                             Group::ProductionCMode2String(Group::ProductionCMode::GRAT));
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
            ss = fmt::format("Group {} has sale rate less then minimum permitted value and is under GRAT control.\n"
                             "The GRAT is increased to meet the sales minimum rate.",
                             group.name());
            this->groupState().update_grat_sales_target(group.name(), production_target);
        //} else if () {//TODO add action for WGASPROD
        //} else if () {//TODO add action for drilling queue
        } else {
            ss = fmt::format("Group {} has sale rate less then minimum permitted value but cannot increase the group production rate \n"
                             "or adjust gas production using WGASPROD or drill new wells to meet the sales target. \n"
                             "Note that WGASPROD and drilling queues are not implemented in Flow. No action is taken.",
                             group.name());
        }
    }
    if (gconsale.sales_target < 0.0) {
        OPM_DEFLOG_THROW(std::runtime_error, "Group " + group.name() + " has sale rate target less then zero. Not implemented in Flow" , deferred_logger);
    }

    if (!ss.empty() && comm_.rank() == 0)
        deferred_logger.info(ss);
}

template<typename Scalar, typename IndexTraits>
bool BlackoilWellModelGeneric<Scalar, IndexTraits>::
checkGroupHigherConstraints(const Group& group,
                            DeferredLogger& deferred_logger,
                            const int reportStepIdx,
                            const int max_number_of_group_switch,
                            const bool update_group_switching_log)
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

    std::vector<Scalar> rates(this->numPhases(), 0.0);

    bool isField = group.name() == "FIELD";
    if (!isField && group.isInjectionGroup()) {
        // Obtain rates for group.
        std::vector<Scalar> resv_coeff_inj(this->numPhases(), 0.0);
        calcInjResvCoeff(fipnum, pvtreg, resv_coeff_inj);

        // checkGroupConstraintsInj considers 'available' rates (e.g., group rates minus reduction rates).
        // So when checking constraints, current groups rate must also be subtracted it's reduction rate
        const std::vector<Scalar> reduction_rates = this->groupState().injection_reduction_rates(group.name());

        for (int phasePos = 0; phasePos < this->numPhases(); ++phasePos) {
            const Scalar local_current_rate = this->wgHelper().sumWellSurfaceRates(
                group, phasePos, /*is_injector=*/true
            );
            // Sum over all processes
            rates[phasePos] = comm_.sum(local_current_rate) - reduction_rates[phasePos];
        }
        const Phase all[] = { Phase::WATER, Phase::OIL, Phase::GAS };
        for (Phase phase : all) {
            const auto currentControl = this->groupState().injection_control(group.name(), phase);
            bool group_is_oscillating = false;
            if (auto groupPos = switched_inj_groups_.find(group.name()); groupPos != switched_inj_groups_.end()) {
                auto& ctrls = groupPos->second[static_cast<std::underlying_type_t<Phase>>(phase)];
                const int number_of_switches = std::count(ctrls.begin(), ctrls.end(), currentControl);
                group_is_oscillating = (number_of_switches >= max_number_of_group_switch);
                if (group_is_oscillating) {
                    const bool output_first_time = (number_of_switches == max_number_of_group_switch);
                    if (output_first_time) {
                        if (comm_.rank() == 0) {
                            std::ostringstream os;
                            os << phase;
                            const std::string msg =
                                fmt::format("Group control for {} injector group {} is oscillating. Group control kept at {}.",
                                            std::move(os).str(),
                                            group.name(),
                                            Group::InjectionCMode2String(currentControl));
                            deferred_logger.info(msg);
                        }
                        ctrls.push_back(currentControl);
                    }
                }
            }

            if (group_is_oscillating) {
                continue;
            }

            // Check higher up only if under individual (not FLD) control.
            if (currentControl != Group::InjectionCMode::FLD && group.injectionGroupControlAvailable(phase)) {
                const Group& parentGroup = schedule().getGroup(group.parent(), reportStepIdx);
                const auto [is_changed, scaling_factor] = this->wgHelper().checkGroupConstraintsInj(
                    group.name(),
                    group.parent(),
                    parentGroup,
                    rates.data(),
                    phase,
                    group.getGroupEfficiencyFactor(),
                    resv_coeff_inj,
                    /*check_guide_rate*/true,
                    deferred_logger
                );
                if (is_changed) {
                    auto& group_log = switched_inj_groups_[group.name()][static_cast<std::underlying_type_t<Phase>>(phase)];
                    if (update_group_switching_log || group_log.empty()) {
                        group_log.push_back(currentControl);
                    }
                    BlackoilWellModelConstraints(*this).
                        actionOnBrokenConstraints(group, Group::InjectionCMode::FLD,
                                                  phase, this->groupState(),
                                                  deferred_logger);
                    this->wgHelper().updateWellRatesFromGroupTargetScale(
                        scaling_factor, group, /*is_injector=*/true, this->wellState()
                    );
                    changed = true;
                }
            }
        }
    }

    if (!isField && group.isProductionGroup()) {
        // Obtain rates for group.
        // checkGroupConstraintsProd considers 'available' rates (e.g., group rates minus reduction rates).
        // So when checking constraints, current groups rate must also be subtracted it's reduction rate
        const std::vector<Scalar> reduction_rates = this->groupState().production_reduction_rates(group.name());
        const Group::ProductionCMode currentControl = this->groupState().production_control(group.name());
        if (auto groupPos = switched_prod_groups_.find(group.name()); groupPos != switched_prod_groups_.end()) {
            auto& ctrls = groupPos->second;
            const int number_of_switches = std::count(ctrls.begin(), ctrls.end(), currentControl);
            const bool group_is_oscillating = (number_of_switches >= max_number_of_group_switch);
            if (group_is_oscillating) {
                const bool output_first_time = (number_of_switches== max_number_of_group_switch);
                if (output_first_time) {
                    if (comm_.rank() == 0) {
                        const std::string msg =
                        fmt::format("Group control for production group {} is oscillating. Group control kept at {}.",
                                    group.name(),
                                    Group::ProductionCMode2String(currentControl));
                        deferred_logger.info(msg);
                    }
                    ctrls.push_back(currentControl);
                }
                return false;
            }
        }
        for (int phasePos = 0; phasePos < this->numPhases(); ++phasePos) {
            const Scalar local_current_rate = this->wgHelper().sumWellSurfaceRates(
                group, phasePos, /*is_injector=*/false
            );
            // Sum over all processes
            rates[phasePos] = -comm_.sum(local_current_rate) - reduction_rates[phasePos];
        }
        std::vector<Scalar> resv_coeff(this->numPhases(), 0.0);
        calcResvCoeff(fipnum, pvtreg, this->groupState().production_rates(group.name()), resv_coeff);
        // Check higher up only if under individual (not FLD) control.
        if (currentControl != Group::ProductionCMode::FLD && group.productionGroupControlAvailable()) {
            const Group& parentGroup = schedule().getGroup(group.parent(), reportStepIdx);
            const auto [is_changed, scaling_factor] = this->wgHelper().checkGroupConstraintsProd(
                group.name(),
                group.parent(),
                parentGroup,
                rates.data(),
                group.getGroupEfficiencyFactor(),
                resv_coeff,
                /*check_guide_rate*/true,
                deferred_logger
            );
            if (is_changed) {
                const auto group_limit_action = group.productionControls(summaryState_).group_limit_action;
                std::optional<std::string> worst_offending_well = std::nullopt;
                changed = BlackoilWellModelConstraints(*this).
                        actionOnBrokenConstraints(group,
                                                  group_limit_action,
                                                  Group::ProductionCMode::FLD,
                                                  worst_offending_well,
                                                  this->groupState(),
                                                  deferred_logger);

                if (changed) {
                    if (update_group_switching_log || switched_prod_groups_[group.name()].empty()) {
                        switched_prod_groups_[group.name()].push_back(currentControl);
                    }
                    this->wgHelper().updateWellRatesFromGroupTargetScale(
                        scaling_factor, group, /*is_injector=*/false, this->wellState()
                    );
                }
            }
        }
    }

    return changed;
}

template<typename Scalar, typename IndexTraits>
void BlackoilWellModelGeneric<Scalar, IndexTraits>::
updateEclWells(const int timeStepIdx,
               const SimulatorUpdate& sim_update,
               const SummaryState& st)
{
    this->updateEclWellsConstraints(timeStepIdx, sim_update, st);

    if (! sim_update.well_structure_changed &&
        ! this->wellStructureChangedDynamically_)
    {
        this->updateEclWellsCTFFromAction(timeStepIdx, sim_update);
    }

    if (sim_update.well_structure_changed) {
        // Note: Record change if this update triggered a well structure
        // change.  Otherwise, we risk setting ChangedDynamically_ to false
        // if a subsequent action at the same time step does *not* change
        // the well topology.
        this->wellStructureChangedDynamically_ = true;
    }
}

template<typename Scalar, typename IndexTraits>
template<typename Iter, typename Body>
void BlackoilWellModelGeneric<Scalar, IndexTraits>::
wellUpdateLoop(Iter first, Iter last, const int timeStepIdx, Body&& body)
{
    std::for_each(first, last,
                  [this, timeStepIdx,
                   loopBody = std::forward<Body>(body)]
                  (const auto& wname)
    {
        auto well_iter = std::find_if(this->wells_ecl_.begin(),
                                      this->wells_ecl_.end(),
                                      [&wname](const auto& well)
                                      { return well.name() == wname; });

        if (well_iter == this->wells_ecl_.end()) {
            return;
        }

        const auto wellIdx =
            std::distance(this->wells_ecl_.begin(), well_iter);

        const auto& well = this->wells_ecl_[wellIdx] =
            this->schedule_.getWell(wname, timeStepIdx);

        loopBody(wellIdx, well);
    });
}

template<typename Scalar, typename IndexTraits>
void BlackoilWellModelGeneric<Scalar, IndexTraits>::
updateEclWellsConstraints(const int              timeStepIdx,
                          const SimulatorUpdate& sim_update,
                          const SummaryState&    st)
{
    this->wellUpdateLoop(sim_update.affected_wells.begin(),
                         sim_update.affected_wells.end(),
                         timeStepIdx,
                         [this, &st]
                         (const auto wellIdx, const auto& well)
    {
        auto& ws = this->wellState().well(wellIdx);
        ws.updateStatus(well.getStatus());
        ws.update_type_and_targets(well, st);
    });
}

template<typename Scalar, typename IndexTraits>
void BlackoilWellModelGeneric<Scalar, IndexTraits>::
updateEclWellsCTFFromAction(const int              timeStepIdx,
                            const SimulatorUpdate& sim_update)
{
    this->wellUpdateLoop(sim_update.welpi_wells.begin(),
                         sim_update.welpi_wells.end(),
                         timeStepIdx,
                         [this](const auto wellIdx, const auto& well)
    {
        auto& pd = this->well_perf_data_[wellIdx];

        {
            auto pdIter = pd.begin();

            for (const auto& conn : well.getConnections()) {
                if (conn.state() != Connection::State::SHUT) {
                    pdIter->connection_transmissibility_factor = conn.CF();
                    ++pdIter;
                }
            }
        }

        this->wellState().well(wellIdx).reset_connection_factors(pd);
        this->prod_index_calc_[wellIdx].reInit(well);
    });
}

template<typename Scalar, typename IndexTraits>
Scalar
BlackoilWellModelGeneric<Scalar, IndexTraits>::
wellPI(const int well_index) const
{
    const auto& pi = this->wellState().well(well_index).productivity_index;

    const auto preferred = this->wells_ecl_[well_index].getPreferredPhase();
    switch (preferred) { // Should really have LIQUID = OIL + WATER here too...
    case Phase::WATER:
        return phase_usage_info_.phaseIsActive(IndexTraits::waterPhaseIdx)
            ? pi[phase_usage_info_.canonicalToActivePhaseIdx(IndexTraits::waterPhaseIdx)]
            : 0.0;

    case Phase::OIL:
        return phase_usage_info_.phaseIsActive(IndexTraits::oilPhaseIdx)
            ? pi[phase_usage_info_.canonicalToActivePhaseIdx(IndexTraits::oilPhaseIdx)]
            : 0.0;

    case Phase::GAS:
        return phase_usage_info_.phaseIsActive(IndexTraits::gasPhaseIdx)
            ? pi[phase_usage_info_.canonicalToActivePhaseIdx(IndexTraits::gasPhaseIdx)]
            : 0.0;

    default:
        throw std::invalid_argument {
            "Unsupported preferred phase " +
            std::to_string(static_cast<int>(preferred))
        };
    }
}

template<typename Scalar, typename IndexTraits>
Scalar
BlackoilWellModelGeneric<Scalar, IndexTraits>::
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

template<typename Scalar, typename IndexTraits>
bool BlackoilWellModelGeneric<Scalar, IndexTraits>::
wasDynamicallyShutThisTimeStep(const int well_index) const
{
    return wasDynamicallyShutThisTimeStep(this->wells_ecl_[well_index].name());
}

template<typename Scalar, typename IndexTraits>
bool
BlackoilWellModelGeneric<Scalar, IndexTraits>::
wasDynamicallyShutThisTimeStep(const std::string& well_name) const
{
    return this->closed_this_step_.find(well_name) !=
           this->closed_this_step_.end();
}

template<typename Scalar, typename IndexTraits>
void BlackoilWellModelGeneric<Scalar, IndexTraits>::
updateWsolvent(const Group& group,
               const int reportStepIdx,
               const WellState<Scalar, IndexTraits>& wellState)
{
    for (const std::string& groupName : group.groups()) {
        const Group& groupTmp = schedule_.getGroup(groupName, reportStepIdx);
        updateWsolvent(groupTmp, reportStepIdx, wellState);
    }

    if (group.isProductionGroup())
        return;

    auto currentGroupControl = this->groupState().injection_control(group.name(), Phase::GAS);
    if( currentGroupControl == Group::InjectionCMode::REIN ) {
        const int gasPos = phase_usage_info_.canonicalToActivePhaseIdx(IndexTraits::gasPhaseIdx);
        const auto& controls = group.injectionControls(Phase::GAS, summaryState_);
        const Group& groupRein = schedule_.getGroup(controls.reinj_group, reportStepIdx);
        Scalar gasProductionRate = this->wgHelper().sumWellSurfaceRates(
            groupRein, gasPos, /*is_injector=*/false
        );
        Scalar solventProductionRate = this->wgHelper().sumSolventRates(
            groupRein, /*is_injector=*/false
        );

        solventProductionRate = comm_.sum(solventProductionRate);
        gasProductionRate = comm_.sum(gasProductionRate);

        Scalar wsolvent = 0.0;
        if (std::abs(gasProductionRate) > 1e-6)
            wsolvent = solventProductionRate / gasProductionRate;

        setWsolvent(group, reportStepIdx, wsolvent);
    }
}

template<typename Scalar, typename IndexTraits>
void BlackoilWellModelGeneric<Scalar, IndexTraits>::
setWsolvent(const Group& group,
            const int reportStepIdx,
            Scalar wsolvent)
{
    for (const std::string& groupName : group.groups()) {
        const Group& groupTmp = schedule_.getGroup(groupName, reportStepIdx);
        setWsolvent(groupTmp, reportStepIdx, wsolvent);
    }

    for (const std::string& wellName : group.wells()) {
        if (! hasOpenLocalWell(wellName))
            continue;

        getGenWell(wellName)->setWsolvent(wsolvent);
    }
}

template<typename Scalar, typename IndexTraits>
template <typename LoopBody>
void BlackoilWellModelGeneric<Scalar, IndexTraits>::
loopOwnedWells(LoopBody&& loopBody) const
{
    auto wellIndex = 0 * this->wells_ecl_.size();

    for (const auto& pwInfo : this->local_parallel_well_info_) {
        if (pwInfo.get().isOwner()) {
            loopBody(wellIndex, this->wells_ecl_[wellIndex]);
        }

        ++wellIndex;
    }
}

template<typename Scalar, typename IndexTraits>
void BlackoilWellModelGeneric<Scalar, IndexTraits>::
assignWellTargets(data::Wells& wsrpt) const
{
    this->loopOwnedWells([this, &wsrpt]
                         ([[maybe_unused]] const auto wellIndex,
                          const Well&                 well)
    {
        // data::Wells (i.e., 'wsrpt') is a std::map<>
        auto& limits = wsrpt[well.name()].limits;

        if (well.isProducer()) {
            this->assignProductionWellTargets(well, limits);
        }
        else {
            this->assignInjectionWellTargets(well, limits);
        }
    });
}

template<typename Scalar, typename IndexTraits>
void BlackoilWellModelGeneric<Scalar, IndexTraits>::
assignProductionWellTargets(const Well& well, data::WellControlLimits& limits) const
{
    using Item = data::WellControlLimits::Item;

    const auto ctrl = well.productionControls(this->summaryState());

    limits
        .set(Item::Bhp, ctrl.bhp_limit)
        .set(Item::OilRate, ctrl.oil_rate)
        .set(Item::WaterRate, ctrl.water_rate)
        .set(Item::GasRate, ctrl.gas_rate)
        .set(Item::ResVRate, ctrl.resv_rate)
        .set(Item::LiquidRate, ctrl.liquid_rate);
}

template<typename Scalar, typename IndexTraits>
void BlackoilWellModelGeneric<Scalar, IndexTraits>::
assignInjectionWellTargets(const Well& well, data::WellControlLimits& limits) const
{
    using Item = data::WellControlLimits::Item;

    const auto ctrl = well.injectionControls(this->summaryState());

    limits
        .set(Item::Bhp, ctrl.bhp_limit)
        .set(Item::ResVRate, ctrl.reservoir_rate);

    if (ctrl.injector_type == InjectorType::MULTI) {
        // Not supported
        return;
    }

    auto rateItem = Item::WaterRate;
    if (ctrl.injector_type == InjectorType::GAS) {
        rateItem = Item::GasRate;
    }
    else if (ctrl.injector_type == InjectorType::OIL) {
        rateItem = Item::OilRate;
    }

    limits.set(rateItem, ctrl.surface_rate);
}

template<typename Scalar, typename IndexTraits>
void BlackoilWellModelGeneric<Scalar, IndexTraits>::
assignDynamicWellStatus(data::Wells& wsrpt) const
{
    this->loopOwnedWells([this, &wsrpt]
                         (const auto wellID,
                          const Well& well)
    {
        auto& xwel = wsrpt[well.name()]; // data::Wells is a std::map<>
        xwel.dynamicStatus = this->wellState().well(well.name()).status;

        // Well was shut this time step we keep it open one last time
        // for output
        if (wasDynamicallyShutThisTimeStep(wellID)) {
            xwel.dynamicStatus = Well::Status::OPEN;
        }
    });
}

template<typename Scalar, typename IndexTraits>
void BlackoilWellModelGeneric<Scalar, IndexTraits>::
assignShutConnections(data::Wells& wsrpt,
                      const int reportStepIndex) const
{
    this->loopOwnedWells([this, reportStepIndex, &wsrpt]
                         ([[maybe_unused]] const auto wellID,
                          const Well&                 well)
    {
        const auto wellIsOpen =
            this->schedule()[reportStepIndex].wells(well.name())
            .getStatus() == Well::Status::OPEN;

        auto skip = [wellIsOpen](const Connection& conn)
        {
            return wellIsOpen && (conn.state() != Connection::State::SHUT);
        };

        auto& xwel = wsrpt[well.name()]; // data::Wells is a std::map<>

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
            xc.d_factor = conn.dFactor();
        }
    });
}

template<typename Scalar, typename IndexTraits>
void BlackoilWellModelGeneric<Scalar, IndexTraits>::
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

template<typename Scalar, typename IndexTraits>
void BlackoilWellModelGeneric<Scalar, IndexTraits>::
assignGroupValues(const int                               reportStepIdx,
                  std::map<std::string, data::GroupData>& gvalues) const
{
    const auto groupGuideRates =
        BlackoilWellModelGuideRates(*this).calculateAllGroupGuideRates(reportStepIdx);

    for (const auto& gname : schedule_.groupNames(reportStepIdx)) {
        const auto& grup = schedule_.getGroup(gname, reportStepIdx);

        auto& gdata = gvalues[gname];
        this->assignGroupControl(grup, gdata);
        BlackoilWellModelGuideRates(*this).assignGroupGuideRates(grup, groupGuideRates, gdata);
    }
}

template<typename Scalar, typename IndexTraits>
void BlackoilWellModelGeneric<Scalar, IndexTraits>::
assignNodeValues(std::map<std::string, data::NodeData>& nodevalues,
                 const int reportStepIdx) const
{
    nodevalues.clear();
    if (reportStepIdx < 0) return;

    for (const auto& [node, pressure] : node_pressures_) {
        nodevalues.emplace(node, data::NodeData{pressure});
        // Assign node values of well groups to GPR:WELLNAME
        const auto& sched = schedule();
        if (!sched.hasGroup(node, reportStepIdx)) {
            continue;
        }
        const auto& group = sched.getGroup(node, reportStepIdx);
        for (const std::string& wellname : group.wells()) {
            nodevalues.emplace(wellname, data::NodeData{pressure});
        }
    }

    const auto& network = schedule()[reportStepIdx].network();
    if (!network.active()) {
        return;
    }

    auto converged_pressures = this->wgHelper().computeNetworkPressures(
        network,
        *(this->vfp_properties_->getProd()),
        this->comm_
    );
    for (const auto& [node, converged_pressure] : converged_pressures) {
        auto it = nodevalues.find(node);
        assert(it != nodevalues.end() );
        it->second.converged_pressure = converged_pressure;
        // Assign node values of group to GPR:WELLNAME
        const auto& sched = schedule();
        if (!sched.hasGroup(node, reportStepIdx)) {
            continue;
        }
        const auto& group = sched.getGroup(node, reportStepIdx);
        for (const std::string& wellname : group.wells()) {
            auto it2 = nodevalues.find(wellname);
            assert(it2 != nodevalues.end());
            it2->second.converged_pressure = converged_pressure;
        }
    }
}

template<typename Scalar, typename IndexTraits>
data::GroupAndNetworkValues
BlackoilWellModelGeneric<Scalar, IndexTraits>::
groupAndNetworkData(const int reportStepIdx) const
{
    auto grp_nwrk_values = data::GroupAndNetworkValues{};

    this->assignGroupValues(reportStepIdx, grp_nwrk_values.groupData);
    this->assignNodeValues(grp_nwrk_values.nodeData, reportStepIdx - 1); // Schedule state info at previous step

    return grp_nwrk_values;
}

template<typename Scalar, typename IndexTraits>
void BlackoilWellModelGeneric<Scalar, IndexTraits>::
updateAndCommunicateGroupData(const int reportStepIdx,
                              const int iterationIdx,
                              const Scalar tol_nupcol,
                              const bool update_wellgrouptarget,
                              DeferredLogger& deferred_logger)
{
    OPM_TIMEFUNCTION();
    const Group& fieldGroup = schedule().getGroup("FIELD", reportStepIdx);
    const int nupcol = schedule()[reportStepIdx].nupcol();

    // Update accumulated group consumption/import rates for current report step
    if(iterationIdx == 0) {
        this->groupState().update_gconsump(schedule(), reportStepIdx, this->summaryState_);
    }

    // This builds some necessary lookup structures, so it must be called
    // before we copy to well_state_nupcol_.
    // Wells may be temporarily stopped due to convergence or operability issues
    // during local well solves. Temporarily stopped wells do not contribute
    // to groups and are therefore not considered in the group target calculations.
    std::vector<WellStatus> well_status(this->numLocalWells(), WellStatus::SHUT);
    for (const auto& well : well_container_generic_) {
        well_status[well->indexOfWell()] = well->wellStatus();
    }

    this->wellState().updateGlobalIsGrup(comm_, well_status);

    WellGroupHelperType &wg_helper = this->wgHelper();
    if (iterationIdx < nupcol) {
        OPM_TIMEBLOCK(updateNupcol);
        this->updateNupcolWGState();
    } else {
        for (const auto& gr_name : schedule().groupNames(reportStepIdx)) {
            const Phase all[] = { Phase::WATER, Phase::OIL, Phase::GAS };
            for (Phase phase : all) {
                if (this->groupState().has_injection_control(gr_name, phase)) {
                    if (this->groupState().injection_control(gr_name, phase) == Group::InjectionCMode::VREP ||
                        this->groupState().injection_control(gr_name, phase) == Group::InjectionCMode::REIN) {
		        OPM_TIMEBLOCK(extraIterationsAfterNupcol);
                        const bool is_vrep = this->groupState().injection_control(gr_name, phase) == Group::InjectionCMode::VREP;
                        const Group& group = schedule().getGroup(gr_name, reportStepIdx);
                        const int np = this->wellState().numPhases();
                        Scalar gr_rate_nupcol = 0.0;
                        {
                            // Temporarily use the nupcol well state for all helper functions
                            // At the end of this scope, the well state will be restored to its original value
                            auto guard = wg_helper.pushWellState(this->nupcolWellState());
                            for (int phaseIdx = 0; phaseIdx < np; ++phaseIdx) {
                                gr_rate_nupcol += wg_helper.sumWellPhaseRates(
                                    /*res_rates=*/is_vrep,
                                    group,
                                    phaseIdx,
                                    /*is_injector=*/false
                                );
                            }
                        }
                        Scalar gr_rate = 0.0;
                        for (int phaseIdx = 0; phaseIdx < np; ++phaseIdx) {
                            gr_rate += wg_helper.sumWellPhaseRates(
                                /*res_rates=*/is_vrep,
                                group,
                                phaseIdx,
                                /*is_injector=*/false
                            );
                        }
                        // sum contributions from owned wells to everybody
                        gr_rate_nupcol = comm_.sum(gr_rate_nupcol);
                        gr_rate = comm_.sum(gr_rate);

                        Scalar small_rate = 1e-12; // m3/s
                        Scalar denominator = (0.5*gr_rate_nupcol + 0.5*gr_rate);
                        Scalar rel_change = denominator > small_rate ? std::abs( (gr_rate_nupcol - gr_rate) / denominator) : 0.0;
                        if ( rel_change > tol_nupcol) {
                            this->updateNupcolWGState();
                            if (comm_.rank() == 0) {
                                const std::string control_str = is_vrep? "VREP" : "REIN";
                                const std::string msg = fmt::format("Group prodution relative change {} larger than tolerance {} "
                                                        "at iteration {}. Update {} for Group {} even if iteration is larger than {} given by NUPCOL." ,
                                                        rel_change, tol_nupcol, iterationIdx, control_str, gr_name, nupcol);
                                deferred_logger.debug(msg);
                            }
                        }
                    }
                }
            }
        }
    }
    {
        constexpr int num_configs = 4;
        constexpr std::array<bool, num_configs> is_production_group = {true, false, false, false};
        constexpr std::array<Phase, num_configs> phases = { Phase::OIL, Phase::WATER, Phase::OIL, Phase::GAS };
        for (int i = 0; i < num_configs; i++) {
            wg_helper.updateGroupControlledWells(
                is_production_group[i], phases[i], this->groupState(), deferred_logger
            );
        }
    }
    // the group target reduction rates needs to be update since wells may have switched to/from GRUP control
    // The group target reduction does not honor NUPCOL.
    wg_helper.updateGroupTargetReduction(fieldGroup, /*is_injector=*/false, this->groupState());
    wg_helper.updateGroupTargetReduction(fieldGroup, /*is_injector=*/true, this->groupState());
    {
        // Temporarily use the nupcol well state for all helper functions
        // At the end of this scope, the well state will be restored to its original value
        auto guard = wg_helper.pushWellState(this->nupcolWellState());
        wg_helper.updateREINForGroups(fieldGroup, /*sum_rank=*/comm_.rank() == 0, this->groupState());
        wg_helper.updateVREPForGroups(fieldGroup, this->groupState());
        wg_helper.updateReservoirRatesInjectionGroups(fieldGroup, this->groupState());
        wg_helper.updateSurfaceRatesInjectionGroups(fieldGroup, this->groupState());
        wg_helper.updateNetworkLeafNodeProductionRates(this->groupState());
        wg_helper.updateGroupProductionRates(fieldGroup, this->groupState());
    }
    wg_helper.updateWellRates(fieldGroup, this->nupcolWellState(), this->wellState());
    this->wellState().communicateGroupRates(comm_);
    this->groupState().communicate_rates(comm_);

    if (update_wellgrouptarget) {
        for (const auto& well : well_container_generic_) {
            const auto& ws = this->wellState().well(well->indexOfWell());
            const auto& group = this->schedule().getGroup(well->wellEcl().groupName(), well->currentStep());
            std::vector<Scalar> resv_coeff(this->numPhases(), 0.0);
            const int fipnum = 0;
            int pvtreg = well->pvtRegionIdx();
            calcResvCoeff(fipnum, pvtreg, this->groupState().production_rates(group.name()), resv_coeff);
            const Scalar efficiencyFactor = well->wellEcl().getEfficiencyFactor() *
                                    ws.efficiency_scaling_factor;
            if (well->isProducer()) {
                std::pair<WellProducerCMode, Scalar> group_target;
                group_target.second = std::numeric_limits<Scalar>::max();
                group_target = wg_helper.getWellGroupTargetProducer(
                    well->name(),
                    well->wellEcl().groupName(),
                    group,
                    ws.surface_rates.data(),
                    efficiencyFactor,
                    resv_coeff,
                    deferred_logger
                );
                auto& ws_update = this->wellState().well(well->indexOfWell());
                ws_update.production_cmode_group_translated = group_target.first;
                ws_update.group_target = group_target.second;
            } else {
                const auto& well_controls = well->wellEcl().injectionControls(summaryState_);
                auto injectorType = well_controls.injector_type;
                Phase injectionPhase;
                switch (injectorType) {
                case InjectorType::WATER:
                {
                    injectionPhase = Phase::WATER;
                    break;
                }
                case InjectorType::OIL:
                {
                    injectionPhase = Phase::OIL;
                    break;
                }
                case InjectorType::GAS:
                {
                    injectionPhase = Phase::GAS;
                    break;
                }
                default:
                    throw std::logic_error("MULTI-phase injection is not supported, but was requested for well " + well->name());
                }
                std::pair<WellInjectorCMode, Scalar> group_target;
                group_target = wg_helper.getWellGroupTargetInjector(
                    well->name(),
                    well->wellEcl().groupName(),
                    group,
                    ws.surface_rates.data(),
                    injectionPhase,
                    efficiencyFactor,
                    resv_coeff,
                    deferred_logger
                );
                auto& ws_update = this->wellState().well(well->indexOfWell());
                ws_update.injection_cmode_group_translated = group_target.first;
                ws_update.group_target = group_target.second;
            }
        }
    }
}

template<typename Scalar, typename IndexTraits>
void BlackoilWellModelGeneric<Scalar, IndexTraits>::
updateNetworkActiveState(const int report_step) {
    const auto& network = schedule()[report_step].network();
    if (!network.active()) {
        this->network_active_ = false;
        return;
    }

    bool network_active = false;
    for (const auto& well : well_container_generic_) {
        const bool is_partof_network = network.has_node(well->wellEcl().groupName());
        const bool prediction_mode = well->wellEcl().predictionMode();
        if (is_partof_network && prediction_mode) {
            network_active = true;
            break;
        }
    }
    this->network_active_ = comm_.max(network_active);
}

template<typename Scalar, typename IndexTraits>
bool BlackoilWellModelGeneric<Scalar, IndexTraits>::
needPreStepNetworkRebalance(const int report_step) const
{
    const auto& network = schedule()[report_step].network();
    bool network_rebalance_necessary = false;
    for (const auto& well : well_container_generic_) {
        const bool is_partof_network = network.has_node(well->wellEcl().groupName());
        // TODO: we might find more relevant events to be included here (including network change events?)
        const auto& events = this->wellState().well(well->indexOfWell()).events;
        if (is_partof_network && events.hasEvent(ScheduleEvents::WELL_STATUS_CHANGE)) {
            network_rebalance_necessary = true;
            break;
        }
    }
    network_rebalance_necessary = comm_.max(network_rebalance_necessary);
    return network_rebalance_necessary;
}

template<typename Scalar, typename IndexTraits>
bool BlackoilWellModelGeneric<Scalar, IndexTraits>::
forceShutWellByName(const std::string& wellname,
                    const double simulation_time,
                    const bool dont_shut_grup_wells)
{
    // Only add the well to the closed list on the
    // process that owns it.
    int well_was_shut = 0;
    const auto it = std::find_if(well_container_generic_.begin(),
                                 well_container_generic_.end(),
                                 [&wellname, &wState = this->wellState(), dont_shut_grup_wells](const auto& well)
                                 {
                                     if (well->name() == wellname) {
                                         // if one well on individual control (typical thp/bhp)
                                         // in a group struggles to converge
                                         // it may lead to problems for the other wells in the group
                                         // we dont want to shut all the wells in a group only the one
                                         // creating the problems.
                                          if (dont_shut_grup_wells) {
                                             const auto& ws = wState.well(well->indexOfWell());
                                             if (well->isInjector()) {
                                                 return ws.injection_cmode != Well::InjectorCMode::GRUP;
                                             } else {
                                                 return ws.production_cmode != Well::ProducerCMode::GRUP;
                                             }
                                         }
                                         return true;
                                     }
                                     else {
                                         return false;
                                     }
                                 });
    if (it != well_container_generic_.end()) {
        wellTestState().close_well(wellname, WellTestConfig::Reason::PHYSICAL, simulation_time);
        well_was_shut = 1;
    }

    // Communicate across processes if a well was shut.
    well_was_shut = comm_.max(well_was_shut);

    // the wellTesteState is updated between timesteps and we also need to update the privous WGstate
    if(well_was_shut)
        this->commitWGState();

    // Only log a message on the output rank.
    if (terminal_output_ && well_was_shut) {
        const std::string msg = "Well " + wellname
            + " will be shut because it fails to converge.";
        OpmLog::info(msg);
    }

    return (well_was_shut == 1);
}

template<typename Scalar, typename IndexTraits>
void BlackoilWellModelGeneric<Scalar, IndexTraits>::
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

template<typename Scalar, typename IndexTraits>
Scalar
BlackoilWellModelGeneric<Scalar, IndexTraits>::
updateNetworkPressures(const int reportStepIdx, const Scalar damping_factor, const Scalar upper_update_bound)
{
    OPM_TIMEFUNCTION();
    // Get the network and return if inactive (no wells in network at this time)
    const auto& network = schedule()[reportStepIdx].network();
    if (!network.active()) {
        return 0.0;
    }

    const auto previous_node_pressures = node_pressures_;

    node_pressures_ = this->wgHelper().computeNetworkPressures(
        network, *(vfp_properties_->getProd()), comm_);

    // here, the network imbalance is the difference between the previous nodal pressure and the new nodal pressure
    Scalar network_imbalance = 0.;
    if (!this->networkActive())
        return network_imbalance;

    if (!previous_node_pressures.empty()) {
        for (const auto& [name, new_pressure]: node_pressures_) {
            if (previous_node_pressures.count(name) <= 0) {
                if (std::abs(new_pressure) > network_imbalance) {
                    network_imbalance = std::abs(new_pressure);
                }
                continue;
            }
            const auto pressure = previous_node_pressures.at(name);
            const Scalar change = (new_pressure - pressure);
            if (std::abs(change) > network_imbalance) {
                network_imbalance = std::abs(change);
            }
            // We dampen the nodal pressure change during one iteration since our nodal pressure calculation
            // is somewhat explicit. There is a relative dampening factor applied to the update value, and also
            // the maximum update is limited (to 5 bar by default, can be changed with --network-max-pressure-update-in-bars).
            const Scalar damped_change = std::min(damping_factor * std::abs(change), upper_update_bound);
            const Scalar sign = change > 0 ? 1. : -1.;
            node_pressures_[name] = pressure + sign * damped_change;
        }
    } else {
        for (const auto& [name, pressure]: node_pressures_) {
            if (std::abs(pressure) > network_imbalance) {
                network_imbalance = std::abs(pressure);
            }
        }
    }

    for (auto& well : well_container_generic_) {

        // Producers only, since we so far only support the
        // "extended" network model (properties defined by
        // BRANPROP and NODEPROP) which only applies to producers.
        if (well->isProducer() && well->wellEcl().predictionMode()) {
            const auto it = node_pressures_.find(well->wellEcl().groupName());
            if (it != node_pressures_.end()) {
                // The well belongs to a group with has a network pressure constraint,
                // set the dynamic THP constraint of the well accordingly.
                const Scalar new_limit = it->second;
                well->setDynamicThpLimit(new_limit);
                SingleWellState<Scalar, IndexTraits>& ws = this->wellState()[well->indexOfWell()];
                const bool thp_is_limit = ws.production_cmode == Well::ProducerCMode::THP;
                // TODO: not sure why the thp is NOT updated properly elsewhere
                if (thp_is_limit) {
                    ws.thp = well->getTHPConstraint(summaryState_);
                }
            }
        }
    }
    return network_imbalance;
}

template<typename Scalar, typename IndexTraits>
void BlackoilWellModelGeneric<Scalar, IndexTraits>::
calculateEfficiencyFactors(const int reportStepIdx)
{
    for (auto& well : well_container_generic_) {
        const Well& wellEcl = well->wellEcl();
        Scalar well_efficiency_factor = wellEcl.getEfficiencyFactor() *
                                        wellState().getGlobalEfficiencyScalingFactor(well->name());
        this->wgHelper().accumulateGroupEfficiencyFactor(
            schedule().getGroup(wellEcl.groupName(), reportStepIdx),
            well_efficiency_factor
        );
        well->setWellEfficiencyFactor(well_efficiency_factor);
    }
}

template<typename Scalar, typename IndexTraits>
WellInterfaceGeneric<Scalar, IndexTraits>*
BlackoilWellModelGeneric<Scalar, IndexTraits>::
getGenWell(const std::string& well_name)
{
    // finding the iterator of the well in wells_ecl
    auto well = std::find_if(well_container_generic_.begin(),
                             well_container_generic_.end(),
                                [&well_name](const WellInterfaceGeneric<Scalar, IndexTraits>* elem)->bool {
                                     return elem->name() == well_name;
                                 });

    assert(well != well_container_generic_.end());

    return *well;
}

template<typename Scalar, typename IndexTraits>
void BlackoilWellModelGeneric<Scalar, IndexTraits>::
setRepRadiusPerfLength()
{
    for (const auto& well : well_container_generic_) {
        well->setRepRadiusPerfLength();
    }
}

template<typename Scalar, typename IndexTraits>
void BlackoilWellModelGeneric<Scalar, IndexTraits>::
updateWellPotentials(const int reportStepIdx,
                     const bool onlyAfterEvent,
                     const SummaryConfig& summaryConfig,
                     DeferredLogger& deferred_logger)
{
    auto well_state_copy = this->wellState();

    const bool write_restart_file = schedule().write_rst_file(reportStepIdx);
    auto exc_type = ExceptionType::NONE;
    std::string exc_msg;
    std::size_t widx = 0;
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
    logAndCheckForProblemsAndThrow(deferred_logger, exc_type,
                                   "updateWellPotentials() failed: " + exc_msg,
                                   terminal_output_, comm_);
}

template<typename Scalar, typename IndexTraits>
void BlackoilWellModelGeneric<Scalar, IndexTraits>::
runWellPIScaling(const int reportStepIdx,
                 DeferredLogger& local_deferredLogger)
{
    if (this->last_run_wellpi_.has_value() && (*this->last_run_wellpi_ == reportStepIdx)) {
        // We've already run WELPI scaling for this report step.  Most
        // common for the very first report step.  Don't redo WELPI scaling.
        return;
    }

    auto hasWellPIEvent = [this, reportStepIdx](const int well_index) -> bool
    {
        return this->schedule()[reportStepIdx].wellgroup_events()
            .hasEvent(this->wells_ecl_[well_index].name(),
                      ScheduleEvents::Events::WELL_PRODUCTIVITY_INDEX);
    };

    auto updateEclWell = [this, reportStepIdx](const int well_index) -> void
    {
        const auto& schedule = this->schedule();
        const auto& wname = this->wells_ecl_[well_index].name();
        this->wells_ecl_[well_index] = schedule.getWell(wname, reportStepIdx);

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
        [this, reportStepIdx](const int    well_index,
                              const Scalar newWellPI) -> void
    {
        const auto& wname = this->wells_ecl_[well_index].name();

        schedule_.applyWellProdIndexScaling(wname, reportStepIdx, newWellPI);
    };

    // Minimal well setup to compute PI/II values
    {
        auto saved_previous_wgstate = this->prevWGState();
        this->commitWGState();

        this->createWellContainer(reportStepIdx);
        this->inferLocalShutWells();

        this->initWellContainer(reportStepIdx);

        this->calculateProductivityIndexValues(local_deferredLogger);
        this->calculateProductivityIndexValuesShutWells(reportStepIdx, local_deferredLogger);

        this->commitWGState(std::move(saved_previous_wgstate));
    }

    const auto nw = this->numLocalWells();
    for (auto wellID = 0*nw; wellID < nw; ++wellID) {
        if (hasWellPIEvent(wellID)) {
            rescaleWellPI(wellID, this->wellPI(wellID));
            updateEclWell(wellID);
        }
    }

    this->last_run_wellpi_ = reportStepIdx;
}

template<typename Scalar, typename IndexTraits>
bool BlackoilWellModelGeneric<Scalar, IndexTraits>::
shouldBalanceNetwork(const int reportStepIdx, const int iterationIdx) const
{
    // if network is not active, we do not need to balance the network
    const auto& network = schedule()[reportStepIdx].network();
    if (!network.active()) {
        return false;
    }

    const auto& balance = schedule()[reportStepIdx].network_balance();
    if (balance.mode() == Network::Balance::CalcMode::TimeStepStart) {
        return iterationIdx == 0;
    } else if (balance.mode() == Network::Balance::CalcMode::NUPCOL) {
        const int nupcol = schedule()[reportStepIdx].nupcol();
        return iterationIdx < nupcol;
    } else {
        // We do not support any other rebalancing modes,
        // i.e. TimeInterval based rebalancing is not available.
        // This should be warned about elsewhere, so we choose to
        // avoid spamming with a warning here.
        return false;
    }
}

template<typename Scalar, typename IndexTraits>
std::vector<int> BlackoilWellModelGeneric<Scalar, IndexTraits>::
getCellsForConnections(const Well& well) const
{
    std::vector<int> wellCells;
    // All possible connections of the well
    const auto& connectionSet = well.getConnections();
    wellCells.reserve(connectionSet.size());

    for (const auto& connection : connectionSet)
    {
        int compressed_idx = well.is_lgr_well()
            ? compressedIndexForInteriorLGR(well.get_lgr_well_tag().value(), connection)
            : this->compressedIndexForInterior(connection.global_index());

        if (compressed_idx >= 0) { // Ignore connections in inactive/remote cells.
            wellCells.push_back(compressed_idx);
        }
    }

    return wellCells;
}

template<typename Scalar, typename IndexTraits>
std::vector<std::string> BlackoilWellModelGeneric<Scalar, IndexTraits>::
getWellsForTesting(const int timeStepIdx,
                   const double simulationTime)
{
  const auto& wtest_config = schedule()[timeStepIdx].wtest_config();
  if (!wtest_config.empty()) { // there is a WTEST request
      return wellTestState().test_wells(wtest_config, simulationTime);
  } else
      return {};
}

template<typename Scalar, typename IndexTraits>
void BlackoilWellModelGeneric<Scalar, IndexTraits>::
assignMassGasRate(data::Wells& wsrpt,
                  const Scalar gasDensity) const
{
    using rt = data::Rates::opt;

    for (auto& wrpt : wsrpt) {
        auto& well_rates = wrpt.second.rates;
        const auto w_mass_rate = well_rates.get(rt::gas, 0.0) * gasDensity;

        well_rates.set(rt::mass_gas, w_mass_rate);
    }
}

template<typename Scalar, typename IndexTraits>
void BlackoilWellModelGeneric<Scalar, IndexTraits>::
assignWellTracerRates(data::Wells& wsrpt,
                      const WellTracerRates& wellTracerRates,
                      const unsigned reportStep) const
{
    if (wellTracerRates.empty()) {
        return; // no tracers
    }

    for (const auto& wTR : wellTracerRates) {
        const auto& eclWell = schedule_.getWell(wTR.first, reportStep);
        auto xwPos = wsrpt.find(eclWell.name());
        if (xwPos == wsrpt.end()) { // No well results.
            continue;
        }
        for (const auto& tr : wTR.second) {
            xwPos->second.rates.set(data::Rates::opt::tracer, tr.rate, tr.name);
        }
    }
}

template<typename Scalar, typename IndexTraits>
void BlackoilWellModelGeneric<Scalar, IndexTraits>::
assignMswTracerRates(data::Wells& wsrpt,
                     const MswTracerRates& mswTracerRates,
                     const unsigned reportStep) const
{
    if (mswTracerRates.empty()) {
        return;
    }

    for (const auto& mswTR : mswTracerRates) {
        const auto& eclWell = schedule_.getWell(mswTR.first, reportStep);
        auto xwPos = wsrpt.find(eclWell.name());
        if (xwPos == wsrpt.end()) { // No well results.
            continue;
        }

        auto& wData = xwPos->second;
        for (const auto& tr : mswTR.second) {
            for (const auto& [seg, rate] : tr.rate) {
                auto segPos = wData.segments.find(seg);
                if (segPos != wData.segments.end()) {
                    segPos->second.rates.set(data::Rates::opt::tracer, rate, tr.name);
                }
            }
        }
    }
}

template<typename Scalar, typename IndexTraits>
std::vector<std::vector<int>> BlackoilWellModelGeneric<Scalar, IndexTraits>::
getMaxWellConnections() const
{
    auto wellConnections = std::vector<std::vector<int>>{};

    auto schedule_wells = this->schedule().wellNames();
    schedule_wells.erase(std::remove_if(schedule_wells.begin(),
                                        schedule_wells.end(),
                                        this->not_on_process_),
                         schedule_wells.end());

    wellConnections.reserve(schedule_wells.size());

    const auto possibleFutureConnections = schedule().getPossibleFutureConnections();

#if HAVE_MPI
    // Communicate Map to other processes, since it is only available on rank 0
    Parallel::MpiSerializer ser(comm_);
    ser.broadcast(Parallel::RootRank{0}, possibleFutureConnections);
#endif

    // initialize the additional cell connections introduced by wells.
    for (const auto& well : schedule_wells) {
        auto& compressed_well_perforations = wellConnections.emplace_back
            (this->getCellsForConnections(this->schedule().back().wells(well)));

        const auto possibleFutureConnectionSetIt = possibleFutureConnections.find(well);
        if (possibleFutureConnectionSetIt != possibleFutureConnections.end()) {
            for (const auto& global_index : possibleFutureConnectionSetIt->second) {
                const int compressed_idx = compressedIndexForInterior(global_index);
                if (compressed_idx >= 0) { // Ignore connections in inactive/remote cells.
                    compressed_well_perforations.push_back(compressed_idx);
                }
            }
        }

        // also include wells with no perforations in case
        std::sort(compressed_well_perforations.begin(),
                  compressed_well_perforations.end());
    }

    return wellConnections;
}

template<typename Scalar, typename IndexTraits>
int BlackoilWellModelGeneric<Scalar, IndexTraits>::numLocalWellsEnd() const
{
    const auto& wnames = schedule().back().well_order().names();

    return std::count_if(wnames.begin(), wnames.end(),
                         [this](const std::string& wname)
                         { return ! this->not_on_process_(wname); });
}

template<typename Scalar, typename IndexTraits>
int BlackoilWellModelGeneric<Scalar, IndexTraits>::numLocalNonshutWells() const
{
    return well_container_generic_.size();
}

template<typename Scalar, typename IndexTraits>
void BlackoilWellModelGeneric<Scalar, IndexTraits>::initInjMult()
{
    for (auto& well : this->well_container_generic_) {
        if (well->isInjector() && well->wellEcl().getInjMultMode() != Well::InjMultMode::NONE) {
            const auto& ws = this->wellState().well(well->indexOfWell());
            const auto& perf_data = ws.perf_data;

            auto& values = this->prev_inj_multipliers_[well->name()];
            if (values.empty()) {
                values.assign(perf_data.size(), 1.0);
            }
            well->initInjMult(values);
        }
    }
}

template<typename Scalar, typename IndexTraits>
void BlackoilWellModelGeneric<Scalar, IndexTraits>::
updateFiltrationModelsPostStep(const double dt,
                               const std::size_t water_index,
                               DeferredLogger& deferred_logger)
{
    for (auto& well : this->well_container_generic_) {
        if (well->isInjector()) {
            const Scalar conc = well->wellEcl().evalFilterConc(this->summaryState_);
            if (conc > 0.) {
                // Update filter cake build-ups (external to the wellbore)
                auto retval = this->filter_cake_
                        .emplace(std::piecewise_construct,
                                 std::forward_as_tuple(well->name()),
                                 std::tuple{});
                auto& fc = retval.first->second;
                fc.updatePostStep(*well, this->wellState(), dt, conc, water_index, deferred_logger);
            }
        }
    }
}

template<typename Scalar, typename IndexTraits>
void BlackoilWellModelGeneric<Scalar, IndexTraits>::
updateInjMult(DeferredLogger& deferred_logger)
{
    for (const auto& well : this->well_container_generic_) {
        if (well->isInjector() && well->wellEcl().getInjMultMode() != Well::InjMultMode::NONE) {
            well->updateInjMult(this->prev_inj_multipliers_[well->name()], deferred_logger);
        }
    }
}

template<typename Scalar, typename IndexTraits>
void BlackoilWellModelGeneric<Scalar, IndexTraits>::
updateFiltrationModelsPreStep(DeferredLogger& deferred_logger)
{
    for (auto& well : this->well_container_generic_) {
        if (well->isInjector()) {
            const auto it = filter_cake_.find(well->name());
            if (it != filter_cake_.end()) {
                it->second.updatePreStep(*well, deferred_logger);
                well->updateFilterCakeMultipliers(it->second.multipliers());
            }
        }
    }
}

template<typename Scalar, typename IndexTraits>
void BlackoilWellModelGeneric<Scalar, IndexTraits>::
logPrimaryVars() const
{
    std::ostringstream os;
    for (const auto& w : this->well_container_generic_) {
        os << w->name() << ":";
        auto pv = w->getPrimaryVars();
        for (const Scalar v : pv) {
            os << ' ' << v;
        }
        os << '\n';
    }
    OpmLog::debug(os.str());
}

template<typename Scalar, typename IndexTraits>
void BlackoilWellModelGeneric<Scalar, IndexTraits>::
reportGroupSwitching(DeferredLogger& local_deferredLogger) const
{
    for (const auto& [name, ctrls] : this->switched_prod_groups_) {
        if (ctrls.empty()) {
            continue;
        }
        const Group::ProductionCMode& currentControl = this->groupState().production_control(name);
        if (ctrls[0] != currentControl) {
            const std::string msg =
                fmt::format("    Production Group {} control model changed from {} to {}",
                            name,
                            Group::ProductionCMode2String(ctrls[0]),
                            Group::ProductionCMode2String(currentControl));
            local_deferredLogger.info(msg);
        }
    }
    for (const auto& [grname, grdata] : this->switched_inj_groups_) {
        const Phase all[] = {Phase::WATER, Phase::OIL, Phase::GAS};
        for (Phase phase : all) {
            if (!this->groupState().has_injection_control(grname, phase)) {
                continue;
            }
            const auto& ctrls = grdata[static_cast<std::underlying_type_t<Phase>>(phase)];
            if (ctrls.empty()) {
                continue;
            }
            const Group::InjectionCMode currentControl =
                this->groupState().injection_control(grname, phase);
            if (ctrls[0] != currentControl) {
                std::ostringstream ss;
                ss << phase;
                const std::string msg =
                    fmt::format("    Injection Group {} (phase = {}) control model changed from {} to {}",
                                grname,
                                ss.str(),
                                Group::InjectionCMode2String(ctrls[0]),
                                Group::InjectionCMode2String(currentControl));
                local_deferredLogger.info(msg);
            }
        }
    }
}

template<typename Scalar, typename IndexTraits>
bool BlackoilWellModelGeneric<Scalar, IndexTraits>::
operator==(const BlackoilWellModelGeneric& rhs) const
{
    return this->initial_step_ == rhs.initial_step_
        && this->report_step_starts_ == rhs.report_step_starts_
        && this->last_run_wellpi_ == rhs.last_run_wellpi_
        && this->local_shut_wells_ == rhs.local_shut_wells_
        && this->closed_this_step_ == rhs.closed_this_step_
        && this->node_pressures_ == rhs.node_pressures_
        && this->last_valid_node_pressures_ == rhs.last_valid_node_pressures_
        && this->prev_inj_multipliers_ == rhs.prev_inj_multipliers_
        && this->active_wgstate_ == rhs.active_wgstate_
        && this->last_valid_wgstate_ == rhs.last_valid_wgstate_
        && this->nupcol_wgstate_ == rhs.nupcol_wgstate_
        && this->switched_prod_groups_ == rhs.switched_prod_groups_
        && this->switched_inj_groups_ == rhs.switched_inj_groups_
        && this->closed_offending_wells_ == rhs.closed_offending_wells_
        && this->gen_gaslift_ == rhs.gen_gaslift_;
}

template class BlackoilWellModelGeneric<double, BlackOilDefaultFluidSystemIndices>;

#if FLOW_INSTANTIATE_FLOAT
template class BlackoilWellModelGeneric<float, BlackOilDefaultFluidSystemIndices>;
#endif

}
