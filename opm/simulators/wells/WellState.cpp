/*
  Copyright 2014 SINTEF ICT, Applied Mathematics.
  Copyright 2017 IRIS AS

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
#include <opm/simulators/wells/WellState.hpp>

#include <opm/common/ErrorMacros.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/Schedule.hpp>
#include <opm/simulators/wells/ParallelWellInfo.hpp>

#include <algorithm>
#include <cassert>
#include <numeric>

namespace Opm
{

void WellState::base_init(const std::vector<double>& cellPressures,
                                               const std::vector<Well>& wells_ecl,
                                               const std::vector<ParallelWellInfo*>& parallel_well_info,
                                               const std::vector<std::vector<PerforationData>>& well_perf_data,
                                               const SummaryState& summary_state)
{
    // clear old name mapping
    this->wellMap_.clear();
    this->perfpress_.clear();
    this->perf_skin_pressure_.clear();
    this->perf_water_throughput_.clear();
    this->perf_water_velocity_.clear();
    this->perfphaserates_.clear();
    this->perfrates_.clear();
    this->perfRateBrine_.clear();
    this->perfRateSolvent_.clear();
    this->perfRatePolymer_.clear();
    this->status_.clear();
    this->well_perf_data_.clear();
    this->parallel_well_info_.clear();
    this->wellrates_.clear();
    this->bhp_.clear();
    this->thp_.clear();
    this->temperature_.clear();
    this->segment_state.clear();
    this->well_potentials_.clear();
    this->productivity_index_.clear();
    this->conn_productivity_index_.clear();
    {
        // const int nw = wells->number_of_wells;
        const int nw = wells_ecl.size();
        // const int np = wells->number_of_phases;
        int connpos = 0;
        for (int w = 0; w < nw; ++w) {
            const Well& well = wells_ecl[w];

            // Initialize bhp(), thp(), wellRates(), temperature().
            initSingleWell(cellPressures, w, well, well_perf_data[w], parallel_well_info[w], summary_state);

            // Setup wellname -> well index mapping.
            const int num_perf_this_well = well_perf_data[w].size();
            std::string name = well.name();
            assert( name.size() > 0 );
            mapentry_t& wellMapEntry = wellMap_[name];
            wellMapEntry[ 0 ] = w;
            wellMapEntry[ 1 ] = connpos;
            wellMapEntry[ 2 ] = num_perf_this_well;
            connpos += num_perf_this_well;
        }
    }
}





void WellState::initSingleWell(const std::vector<double>& cellPressures,
                                                    const int w,
                                                    const Well& well,
                                                    const std::vector<PerforationData>& well_perf_data,
                                                    const ParallelWellInfo* well_info,
                                                    const SummaryState& summary_state)
{
    assert(well.isInjector() || well.isProducer());

    // Set default zero initial well rates.
    // May be overwritten below.
    const auto& pu = this->phase_usage_;
    const int np = pu.num_phases;

    this->status_.add(well.name(), Well::Status::OPEN);
    this->well_perf_data_.add(well.name(), well_perf_data);
    this->parallel_well_info_.add(well.name(), well_info);
    this->wellrates_.add(well.name(), std::vector<double>(np, 0));
    this->well_potentials_.add(well.name(), std::vector<double>(np, 0));
    const int num_perf_this_well = well_info->communication().sum(well_perf_data_[w].size());
    this->segment_state.add(well.name(), SegmentState{});
    this->perfpress_.add(well.name(), std::vector<double>(num_perf_this_well, -1e100));
    this->perfrates_.add(well.name(), std::vector<double>(num_perf_this_well, 0));
    this->perfphaserates_.add(well.name(), std::vector<double>(np*num_perf_this_well, 0));
    this->perf_skin_pressure_.add(well.name(), std::vector<double>(num_perf_this_well, 0));
    this->perf_water_velocity_.add(well.name(), std::vector<double>(num_perf_this_well, 0));
    this->perf_water_throughput_.add(well.name(), std::vector<double>(num_perf_this_well, 0));
    this->perfRatePolymer_.add(well.name(), std::vector<double>(num_perf_this_well, 0));
    this->perfRateSolvent_.add(well.name(), std::vector<double>(num_perf_this_well, 0));
    this->perfRateBrine_.add(well.name(), std::vector<double>(num_perf_this_well, 0));
    this->bhp_.add(well.name(), 0.0);
    this->thp_.add(well.name(), 0.0);
    this->productivity_index_.add(well.name(), std::vector<double>(np, 0));
    this->conn_productivity_index_.add(well.name(), std::vector<double>(num_perf_this_well * np, 0));
    if ( well.isInjector() )
        this->temperature_.add(well.name(), well.injectionControls(summary_state).temperature);
    else
        this->temperature_.add(well.name(), 273.15 + 15.56); // standard condition temperature

    if ( num_perf_this_well == 0 )
        return;

    const auto inj_controls = well.isInjector() ? well.injectionControls(summary_state) : Well::InjectionControls(0);
    const auto prod_controls = well.isProducer() ? well.productionControls(summary_state) : Well::ProductionControls(0);

    const bool is_bhp = well.isInjector() ? (inj_controls.cmode == Well::InjectorCMode::BHP)
                                          : (prod_controls.cmode == Well::ProducerCMode::BHP);
    const double bhp_limit = well.isInjector() ? inj_controls.bhp_limit : prod_controls.bhp_limit;
    const bool is_grup = well.isInjector() ? (inj_controls.cmode == Well::InjectorCMode::GRUP)
                                           : (prod_controls.cmode == Well::ProducerCMode::GRUP);

    const double inj_surf_rate = well.isInjector() ? inj_controls.surface_rate : 0.0; // To avoid a "maybe-uninitialized" warning.

    const double local_pressure = well_perf_data_[w].empty() ?
                                      0 : cellPressures[well_perf_data_[w][0].cell_index];
    const double global_pressure = well_info->broadcastFirstPerforationValue(local_pressure);

    if (well.getStatus() == Well::Status::OPEN) {
        this->status_[w] = Well::Status::OPEN;
    }

    if (well.getStatus() == Well::Status::STOP) {
        // Stopped well:
        // 1. Rates: zero well rates.
        // 2. Bhp: assign bhp equal to bhp control, if
        //    applicable, otherwise assign equal to
        //    first perforation cell pressure.
        if (is_bhp) {
            bhp_[w] = bhp_limit;
        } else {
            bhp_[w] = global_pressure;
        }
    } else if (is_grup) {
        // Well under group control.
        // 1. Rates: zero well rates.
        // 2. Bhp: initialize bhp to be a
        //    little above or below (depending on if
        //    the well is an injector or producer)
        //    pressure in first perforation cell.
        const double safety_factor = well.isInjector() ? 1.01 : 0.99;
        bhp_[w] = safety_factor * global_pressure;
    } else {
        // Open well, under own control:
        // 1. Rates: initialize well rates to match
        //    controls if type is ORAT/GRAT/WRAT
        //    (producer) or RATE (injector).
        //    Otherwise, we cannot set the correct
        //    value here and initialize to zero rate.
        auto& rates = this->wellrates_[w];
        if (well.isInjector()) {
            if (inj_controls.cmode == Well::InjectorCMode::RATE) {
                switch (inj_controls.injector_type) {
                case InjectorType::WATER:
                    assert(pu.phase_used[BlackoilPhases::Aqua]);
                    rates[pu.phase_pos[BlackoilPhases::Aqua]] = inj_surf_rate;
                    break;
                case InjectorType::GAS:
                    assert(pu.phase_used[BlackoilPhases::Vapour]);
                    rates[pu.phase_pos[BlackoilPhases::Vapour]] = inj_surf_rate;
                    break;
                case InjectorType::OIL:
                    assert(pu.phase_used[BlackoilPhases::Liquid]);
                    rates[pu.phase_pos[BlackoilPhases::Liquid]] = inj_surf_rate;
                    break;
                case InjectorType::MULTI:
                    // Not currently handled, keep zero init.
                    break;
                }
            } else {
                // Keep zero init.
            }
        } else {
            assert(well.isProducer());
            // Note negative rates for producing wells.
            switch (prod_controls.cmode) {
            case Well::ProducerCMode::ORAT:
                assert(pu.phase_used[BlackoilPhases::Liquid]);
                rates[pu.phase_pos[BlackoilPhases::Liquid]] = -prod_controls.oil_rate;
                break;
            case Well::ProducerCMode::WRAT:
                assert(pu.phase_used[BlackoilPhases::Aqua]);
                rates[pu.phase_pos[BlackoilPhases::Aqua]] = -prod_controls.water_rate;
                break;
            case Well::ProducerCMode::GRAT:
                assert(pu.phase_used[BlackoilPhases::Vapour]);
                rates[pu.phase_pos[BlackoilPhases::Vapour]] = -prod_controls.gas_rate;
                break;
            default:
                // Keep zero init.
                break;
            }
        }
        // 2. Bhp: initialize bhp to be target pressure if
        //    bhp-controlled well, otherwise set to a
        //    little above or below (depending on if
        //    the well is an injector or producer)
        //    pressure in first perforation cell.
        if (is_bhp) {
            bhp_[w] = bhp_limit;
        } else {
            const double safety_factor = well.isInjector() ? 1.01 : 0.99;
            bhp_[w] = safety_factor * global_pressure;
        }
    }

    // 3. Thp: assign thp equal to thp target/limit, if such a limit exists,
    //    otherwise keep it zero.
    const bool has_thp = well.isInjector() ? inj_controls.hasControl(Well::InjectorCMode::THP)
                                           : prod_controls.hasControl(Well::ProducerCMode::THP);
    const double thp_limit = well.isInjector() ? inj_controls.thp_limit : prod_controls.thp_limit;
    if (has_thp) {
        thp_[w] = thp_limit;
    }
}



void WellState::init(const std::vector<double>& cellPressures,
                                          const Schedule& schedule,
                                          const std::vector<Well>& wells_ecl,
                                          const std::vector<ParallelWellInfo*>& parallel_well_info,
                                          const int report_step,
                                          const WellState* prevState,
                                          const std::vector<std::vector<PerforationData>>& well_perf_data,
                                          const SummaryState& summary_state)
{
    // call init on base class
    this->base_init(cellPressures, wells_ecl, parallel_well_info, well_perf_data, summary_state);
    this->global_well_info = std::make_optional<GlobalWellInfo>( schedule, report_step, wells_ecl );
    for (const auto& wname : schedule.wellNames(report_step))
    {
        well_rates.insert({wname, std::make_pair(false, std::vector<double>(this->numPhases()))});
    }
    for (const auto& winfo: parallel_well_info)
    {
        well_rates[winfo->name()].first = winfo->isOwner();
    }

    const int nw = wells_ecl.size();

    if( nw == 0 ) return ;

    // Initialize perfphaserates_, which must be done here.
    const auto& pu = this->phaseUsage();
    const int np = pu.num_phases;

    int nperf = 0;
    for (const auto& wpd : well_perf_data) {
        nperf += wpd.size();
    }

    well_reservoir_rates_.clear();
    well_dissolved_gas_rates_.clear();
    well_vaporized_oil_rates_.clear();

    this->events_.clear();
    {
        const auto& wg_events = schedule[report_step].wellgroup_events();
        for (const auto& ecl_well : wells_ecl) {
            const auto& wname = ecl_well.name();
            if (wg_events.has(wname))
                this->events_.add( wname, wg_events.at(wname) );
            else
                this->events_.add( wname, Events() );
        }
    }

    first_perf_index_.resize(nw, 0);
    first_perf_index_[0] = 0;
    for (int w = 0; w < nw; ++w) {
        // Initialize perfphaserates_ to well
        // rates divided by the number of perforations.
        const auto& wname = wells_ecl[w].name();
        const auto& well_info = this->wellMap().at(wname);
        const int connpos = well_info[1];
        const int num_perf_this_well = well_info[2];
        const int global_num_perf_this_well = parallel_well_info[w]->communication().sum(num_perf_this_well);
        auto& perf_press = this->perfPress(w);

        first_perf_index_[w] = connpos;
        auto& phase_rates = this->perfPhaseRates(w);

        for (int perf = 0; perf < num_perf_this_well; ++perf) {
            if (wells_ecl[w].getStatus() == Well::Status::OPEN) {
                for (int p = 0; p < this->numPhases(); ++p) {
                    phase_rates[this->numPhases()*perf + p] = wellRates(w)[p] / double(global_num_perf_this_well);
                }
            }
            perf_press[perf] = cellPressures[well_perf_data[w][perf].cell_index];
        }

        this->well_reservoir_rates_.add(wname, std::vector<double>(np, 0));
        this->well_dissolved_gas_rates_.add(wname, 0);
        this->well_vaporized_oil_rates_.add(wname, 0);
    }

    is_producer_.clear();
    for (int w = 0; w < nw; ++w) {
        const auto& ecl_well = wells_ecl[w];
        this->is_producer_.add( ecl_well.name(), ecl_well.isProducer());
    }

    current_injection_controls_.clear();
    current_production_controls_.clear();
    for (int w = 0; w < nw; ++w) {
        const auto& wname = wells_ecl[w].name();
        current_production_controls_.add(wname, Well::ProducerCMode::CMODE_UNDEFINED);
        current_injection_controls_.add(wname, Well::InjectorCMode::CMODE_UNDEFINED);
        if (wells_ecl[w].isProducer()) {
            const auto controls = wells_ecl[w].productionControls(summary_state);
            currentProductionControl(w, controls.cmode);
        }
        else {
            const auto controls = wells_ecl[w].injectionControls(summary_state);
            currentInjectionControl(w, controls.cmode);
        }
    }


    for (int w = 0; w < nw; ++w) {
        switch (wells_ecl[w].getStatus()) {
        case Well::Status::SHUT:
            this->shutWell(w);
            break;

        case Well::Status::STOP:
            this->stopWell(w);
            break;

        default:
            this->openWell(w);
            break;
        }
    }

    // intialize wells that have been there before
    // order may change so the mapping is based on the well name
    if (prevState && !prevState->wellMap().empty()) {
        auto end = prevState->wellMap().end();
        for (int w = 0; w < nw; ++w) {
            const Well& well = wells_ecl[w];
            if (well.getStatus() == Well::Status::SHUT) {
                continue;
            }
            const auto& wname = well.name();
            auto it = prevState->wellMap().find(well.name());
            if (it != end)
            {
                const int newIndex = w;
                const int oldIndex = it->second[ 0 ];
                if (prevState->status_[oldIndex] == Well::Status::SHUT) {
                    // Well was shut in previous state, do not use its values.
                    continue;
                }

                if (is_producer_[newIndex] != prevState->is_producer_[oldIndex]) {
                    // Well changed to/from injector from/to producer, do not use its privious values.
                    continue;
                }

                // bhp
                this->update_bhp( newIndex, prevState->bhp( oldIndex ));

                // thp
                this->update_thp( newIndex, prevState->thp( oldIndex ));

                // If new target is set using WCONPROD, WCONINJE etc. we use the new control
                if (!this->events_[w].hasEvent(WellState::event_mask)) {
                    current_injection_controls_[ newIndex ] = prevState->currentInjectionControl(oldIndex);
                    current_production_controls_[ newIndex ] = prevState->currentProductionControl(oldIndex);
                }

                wellRates(w) = prevState->wellRates(oldIndex);
                wellReservoirRates(w) = prevState->wellReservoirRates(oldIndex);

                // Well potentials
                for (int p=0; p < np; p++) {
                    this->wellPotentials(newIndex)[p] = prevState->wellPotentials(oldIndex)[p];
                }

                // perfPhaseRates
                const int num_perf_old_well = (*it).second[ 2 ];
                const auto new_iter = this->wellMap().find(well.name());
                if (new_iter == this->wellMap().end()) {
                    throw std::logic_error {
                        well.name() + " is not in internal well map - "
                        "Bug in WellState"
                    };
                }

                const int num_perf_this_well = new_iter->second[2];

                const int num_perf_changed = parallel_well_info[w]->communication()
                    .sum(static_cast<int>(num_perf_old_well != num_perf_this_well));
                const bool global_num_perf_same = num_perf_changed == 0;


                // copy perforation rates when the number of
                // perforations is equal, otherwise initialize
                // perfphaserates to well rates divided by the
                // number of perforations.
                if (global_num_perf_same)
                {
                    const auto& src_rates = prevState->perfPhaseRates(oldIndex);
                    auto& target_rates = this->perfPhaseRates(newIndex);
                    for (int perf_index = 0; perf_index < num_perf_this_well; perf_index++) {
                        for (int p = 0; p < np; p++) {
                            target_rates[perf_index*np + p] = src_rates[perf_index*np + p];
                        }
                    }
                } else {
                    const int global_num_perf_this_well = parallel_well_info[w]->communication().sum(num_perf_this_well);
                    auto& target_rates = this->perfPhaseRates(newIndex);
                    for (int perf_index = 0; perf_index < num_perf_this_well; perf_index++) {
                        for (int p = 0; p < np; ++p) {
                            target_rates[perf_index*np + p] = wellRates(w)[p] / double(global_num_perf_this_well);
                        }
                    }
                }

                // perfPressures
                if (global_num_perf_same)
                {
                    auto& target_press = perfPress(w);
                    const auto& src_press = prevState->perfPress(well.name());
                    for (int perf = 0; perf < num_perf_this_well; ++perf)
                    {
                        target_press[perf] = src_press[perf];
                    }
                }

                // perfSolventRates
                if (pu.has_solvent) {
                    if (global_num_perf_same)
                    {
                        this->perfRateSolvent_.copy_welldata(prevState->perfRateSolvent_, wname);
                    }
                }

                // polymer injectivity related
                //
                // here we did not consider the case that we close
                // some perforation during the running and also,
                // wells can be shut and re-opened
                if (pu.has_polymermw) {
                    if (global_num_perf_same)
                    {
                        this->perf_water_velocity_.copy_welldata(prevState->perf_water_velocity_, wname);
                        this->perf_skin_pressure_.copy_welldata(prevState->perf_skin_pressure_, wname);
                        this->perf_water_throughput_.copy_welldata(prevState->perf_water_throughput_, wname);
                    }
                }

                    // Productivity index.
                this->productivity_index_.copy_welldata( prevState->productivity_index_, wname );
            }

            // If in the new step, there is no THP related
            // target/limit anymore, its thp value should be set to
            // zero.
            const bool has_thp = well.isInjector()
                ? well.injectionControls (summary_state).hasControl(Well::InjectorCMode::THP)
                : well.productionControls(summary_state).hasControl(Well::ProducerCMode::THP);

            if (!has_thp) {
                this->update_thp(w, 0.0);
            }
        }
    }


    updateWellsDefaultALQ(wells_ecl);
    do_glift_optimization_ = true;
}

void WellState::resize(const std::vector<Well>& wells_ecl,
                                            const std::vector<ParallelWellInfo*>& parallel_well_info,
                                            const Schedule& schedule,
                                            const bool handle_ms_well,
                                            const size_t numCells,
                                            const std::vector<std::vector<PerforationData>>& well_perf_data,
                                            const SummaryState& summary_state)
{
    const std::vector<double> tmp(numCells, 0.0); // <- UGLY HACK to pass the size
    init(tmp, schedule, wells_ecl, parallel_well_info, 0, nullptr, well_perf_data, summary_state);

    if (handle_ms_well) {
        initWellStateMSWell(wells_ecl, nullptr);
    }
}

const std::vector<double>&
WellState::currentWellRates(const std::string& wellName) const
{
    auto it = well_rates.find(wellName);

    if (it == well_rates.end())
        OPM_THROW(std::logic_error, "Could not find any rates for well  " << wellName);

    return it->second.second;
}

template<class Communication>
void WellState::gatherVectorsOnRoot(const std::vector<data::Connection>& from_connections,
                                                         std::vector<data::Connection>& to_connections,
                                                         const Communication& comm) const
{
    int size = from_connections.size();
    std::vector<int> sizes;
    std::vector<int> displ;
    if (comm.rank()==0){
        sizes.resize(comm.size());
    }
    comm.gather(&size, sizes.data(), 1, 0);

    if (comm.rank()==0){
        displ.resize(comm.size()+1, 0);
        std::partial_sum(sizes.begin(), sizes.end(), displ.begin()+1);
        to_connections.resize(displ.back());
    }
    comm.gatherv(from_connections.data(), size, to_connections.data(),
                 sizes.data(), displ.data(), 0);
}

data::Wells
WellState::report(const int* globalCellIdxMap,
                                       const std::function<bool(const int)>& wasDynamicallyClosed) const
{
    if (this->numWells() == 0)
        return {};

    using rt = data::Rates::opt;
    const auto& pu = this->phaseUsage();
    const int np = pu.num_phases;

    data::Wells res;
    for( const auto& itr : this->wellMap() ) {
        const auto well_index = itr.second[ 0 ];
        if ((this->status_[well_index] == Well::Status::SHUT) &&
            ! wasDynamicallyClosed(well_index))
        {
            continue;
        }

        const auto& pwinfo = *this->parallel_well_info_[well_index];
        using WellT = std::remove_reference_t<decltype(res[ itr.first ])>;
        WellT dummyWell; // dummy if we are not owner
        auto& well = pwinfo.isOwner() ? res[ itr.first ] : dummyWell;
        well.bhp = this->bhp(well_index);
        well.thp = this->thp( well_index );
        well.temperature = this->temperature( well_index );

        const auto& wv = this->wellRates(well_index);
        if( pu.phase_used[BlackoilPhases::Aqua] ) {
            well.rates.set( rt::wat, wv[ pu.phase_pos[BlackoilPhases::Aqua] ] );
        }

        if( pu.phase_used[BlackoilPhases::Liquid] ) {
            well.rates.set( rt::oil, wv[ pu.phase_pos[BlackoilPhases::Liquid] ] );
        }

        if( pu.phase_used[BlackoilPhases::Vapour] ) {
            well.rates.set( rt::gas, wv[ pu.phase_pos[BlackoilPhases::Vapour] ] );
        }

        if (pwinfo.communication().size()==1)
        {
            reportConnections(well, pu, itr, globalCellIdxMap);
        }
        else
        {
            assert(pwinfo.communication().rank() != 0 || &dummyWell != &well);
            // report the local connections
            reportConnections(dummyWell, pu, itr, globalCellIdxMap);
            // gather them to well on root.
            gatherVectorsOnRoot(dummyWell.connections, well.connections,
                                pwinfo.communication());
        }
    }

    std::vector<rt> phs(np);
    if (pu.phase_used[Water]) {
        phs.at( pu.phase_pos[Water] ) = rt::wat;
    }

    if (pu.phase_used[Oil]) {
        phs.at( pu.phase_pos[Oil] ) = rt::oil;
    }

    if (pu.phase_used[Gas]) {
        phs.at( pu.phase_pos[Gas] ) = rt::gas;
    }

    // This is a reference or example on **how** to convert from
    // WellState to something understood by opm-common's output
    // layer.  It is intended to be properly implemented and
    // maintained as a part of simulators, as it relies on simulator
    // internals, details and representations.

    for (const auto& wt : this->wellMap()) {
        const auto w = wt.second[ 0 ];
        if (((this->status_[w] == Well::Status::SHUT) &&
             ! wasDynamicallyClosed(w)) ||
            ! this->parallel_well_info_[w]->isOwner())
        {
            continue;
        }

        auto& well = res.at(wt.first);
        const auto& reservoir_rates = this->well_reservoir_rates_[w];
        const auto& well_potentials = this->well_potentials_[w];
        const auto& wpi = this->productivity_index_[w];

        if (pu.phase_used[Water]) {
            well.rates.set(rt::reservoir_water, reservoir_rates[pu.phase_pos[Water]]);
            well.rates.set(rt::productivity_index_water, wpi[pu.phase_pos[Water]]);
            well.rates.set(rt::well_potential_water, well_potentials[pu.phase_pos[Water]]);
        }

        if (pu.phase_used[Oil]) {
            well.rates.set(rt::reservoir_oil, reservoir_rates[pu.phase_pos[Oil]]);
            well.rates.set(rt::productivity_index_oil, wpi[pu.phase_pos[Oil]]);
            well.rates.set(rt::well_potential_oil, well_potentials[pu.phase_pos[Oil]]);
        }

        if (pu.phase_used[Gas]) {
            well.rates.set(rt::reservoir_gas, reservoir_rates[pu.phase_pos[Gas]]);
            well.rates.set(rt::productivity_index_gas, wpi[pu.phase_pos[Gas]]);
            well.rates.set(rt::well_potential_gas, well_potentials[pu.phase_pos[Gas]]);
        }

        if (pu.has_solvent || pu.has_zFraction) {
            well.rates.set(rt::solvent, solventWellRate(w));
        }

        if (pu.has_polymer) {
            well.rates.set(rt::polymer, polymerWellRate(w));
        }

        if (pu.has_brine) {
            well.rates.set(rt::brine, brineWellRate(w));
        }

        if (is_producer_[w]) {
            well.rates.set(rt::alq, getALQ(/*wellName=*/wt.first));
        }
        else {
            well.rates.set(rt::alq, 0.0);
        }

        well.rates.set(rt::dissolved_gas, this->well_dissolved_gas_rates_[w]);
        well.rates.set(rt::vaporized_oil, this->well_vaporized_oil_rates_[w]);

        {
            auto& curr = well.current_control;

            curr.isProducer = this->is_producer_[w];
            curr.prod = this->currentProductionControl(w);
            curr.inj  = this->currentInjectionControl(w);
        }

        const auto nseg = this->numSegments(w);
        for (auto seg_ix = 0*nseg; seg_ix < nseg; ++seg_ix) {
            const auto seg_no = this->segmentNumber(w, seg_ix);
            well.segments[seg_no] =
                    this->reportSegmentResults(pu, w, seg_ix, seg_no);
        }
    }

    return res;
}

void WellState::reportConnections(data::Well& well,
                                                       const PhaseUsage &pu,
                                                       const WellMapType::value_type& wt,
                                                       const int* globalCellIdxMap) const
{
    using rt = data::Rates::opt;
    const auto well_index = wt.second[ 0 ];
    const auto& pd = this->well_perf_data_[well_index];
    const int num_perf_well = pd.size();
    well.connections.resize(num_perf_well);

    const auto& perf_rates = this->perfRates(well_index);
    const auto& perf_pressure = this->perfPress(well_index);
    for( int i = 0; i < num_perf_well; ++i ) {
        const auto active_index = this->well_perf_data_[well_index][i].cell_index;
        auto& connection = well.connections[ i ];
        connection.index = globalCellIdxMap[active_index];
        connection.pressure = perf_pressure[i];
        connection.reservoir_rate = perf_rates[i];
        connection.trans_factor = pd[i].connection_transmissibility_factor;
    }
    assert(num_perf_well == int(well.connections.size()));


    const int np = pu.num_phases;
    size_t local_comp_index = 0;
    std::vector< rt > phs( np );
    std::vector<rt> pi(np);
    if( pu.phase_used[Water] ) {
        phs.at( pu.phase_pos[Water] ) = rt::wat;
        pi .at( pu.phase_pos[Water] ) = rt::productivity_index_water;
    }

    if( pu.phase_used[Oil] ) {
        phs.at( pu.phase_pos[Oil] ) = rt::oil;
        pi .at( pu.phase_pos[Oil] ) = rt::productivity_index_oil;
    }

    if( pu.phase_used[Gas] ) {
        phs.at( pu.phase_pos[Gas] ) = rt::gas;
        pi .at( pu.phase_pos[Gas] ) = rt::productivity_index_gas;
    }
    for( auto& comp : well.connections) {
        const auto * rates = &this->perfPhaseRates(well_index)[np*local_comp_index];
        const auto& connPI  = this->connectionProductivityIndex(well_index);

        for( int i = 0; i < np; ++i ) {
            comp.rates.set( phs[ i ], rates[i] );
            comp.rates.set( pi [ i ], connPI[i] );
        }
        if ( pu.has_polymer ) {
            const auto& perf_polymer_rate = this->perfRatePolymer(well_index);
            comp.rates.set( rt::polymer, perf_polymer_rate[local_comp_index]);
        }
        if ( pu.has_brine ) {
            const auto& perf_brine_rate = this->perfRateBrine(well_index);
            comp.rates.set( rt::brine, perf_brine_rate[local_comp_index]);
        }
        if ( pu.has_solvent ) {
            const auto& perf_solvent_rate = this->perfRateSolvent(well_index);
            comp.rates.set( rt::solvent, perf_solvent_rate[local_comp_index] );
        }

        ++local_comp_index;
    }
    assert(local_comp_index == this->well_perf_data_[wt.second[0]].size());
}

void WellState::initWellStateMSWell(const std::vector<Well>& wells_ecl,
                                                         const WellState* prev_well_state)
{
    // still using the order in wells
    const int nw = wells_ecl.size();
    if (nw == 0) {
        return;
    }
    const auto& pu = this->phaseUsage();
    const int np = pu.num_phases;

    // in the init function, the well rates and perforation rates have been initialized or copied from prevState
    // what we do here, is to set the segment rates and perforation rates
    for (int w = 0; w < nw; ++w) {
        const auto& well_ecl = wells_ecl[w];

        if ( well_ecl.isMultiSegment() ) {
            const WellSegments& segment_set = well_ecl.getSegments();
            // assuming the order of the perforations in well_ecl is the same with Wells
            const WellConnections& completion_set = well_ecl.getConnections();
            // number of segment for this single well
            this->segment_state.update(w, SegmentState{np, segment_set});
            const int well_nseg = segment_set.size();
            int n_activeperf = 0;

            // we need to know for each segment, how many perforation it has and how many segments using it as outlet_segment
            // that is why I think we should use a well model to initialize the WellState here
            std::vector<std::vector<int>> segment_perforations(well_nseg);
            for (size_t perf = 0; perf < completion_set.size(); ++perf) {
                const Connection& connection = completion_set.get(perf);
                if (connection.state() == Connection::State::OPEN) {
                    const int segment_index = segment_set.segmentNumberToIndex(connection.segment());
                    segment_perforations[segment_index].push_back(n_activeperf);
                    n_activeperf++;
                }
            }

            std::vector<std::vector<int>> segment_inlets(well_nseg);
            for (int seg = 0; seg < well_nseg; ++seg) {
                const Segment& segment = segment_set[seg];
                const int segment_number = segment.segmentNumber();
                const int outlet_segment_number = segment.outletSegment();
                if (outlet_segment_number > 0) {
                    const int segment_index = segment_set.segmentNumberToIndex(segment_number);
                    const int outlet_segment_index = segment_set.segmentNumberToIndex(outlet_segment_number);
                    segment_inlets[outlet_segment_index].push_back(segment_index);
                }
            }


            // for the seg_rates_, now it becomes a recursive solution procedure.
            {
                // make sure the information from wells_ecl consistent with wells
                assert((n_activeperf == num_perf_this_well) &&
                       "Inconsistent number of reservoir connections in well");

                if (pu.phase_used[Gas]) {
                    auto& perf_rates = this->perfPhaseRates(w);
                    const int gaspos = pu.phase_pos[Gas];
                    // scale the phase rates for Gas to avoid too bad initial guess for gas fraction
                    // it will probably benefit the standard well too, while it needs to be justified
                    // TODO: to see if this strategy can benefit StandardWell too
                    // TODO: it might cause big problem for gas rate control or if there is a gas rate limit
                    // maybe the best way is to initialize the fractions first then get the rates
                    for (int perf = 0; perf < n_activeperf; perf++)
                        perf_rates[perf*np + gaspos] *= 100;
                }

                const auto& perf_rates = this->perfPhaseRates(w);
                std::vector<double> perforation_rates(perf_rates.begin(), perf_rates.end());

                auto& segments = this->segments(w);
                calculateSegmentRates(segment_inlets, segment_perforations, perforation_rates, np, 0 /* top segment */, segments.rates);
            }
            // for the segment pressure, the segment pressure is the same with the first perforation belongs to the segment
            // if there is no perforation associated with this segment, it uses the pressure from the outlet segment
            // which requres the ordering is successful
            // Not sure what is the best way to handle the initialization, hopefully, the bad initialization can be
            // improved during the solveWellEq process
            {
                // top segment is always the first one, and its pressure is the well bhp
                auto& segment_pressure = this->segments(w).pressure;
                segment_pressure[0] = bhp(w);
                const auto& perf_press = this->perfPress(w);
                for (int seg = 1; seg < well_nseg; ++seg) {
                    if ( !segment_perforations[seg].empty() ) {
                        const int first_perf = segment_perforations[seg][0];
                        segment_pressure[seg] = perf_press[first_perf];
                    } else {
                        // seg_press_.push_back(bhp); // may not be a good decision
                        // using the outlet segment pressure // it needs the ordering is correct
                        const int outlet_seg = segment_set[seg].outletSegment();
                        segment_pressure[seg] = segment_pressure[segment_set.segmentNumberToIndex(outlet_seg)];
                    }
                }
            }
        }
    }


    if (prev_well_state) {
        for (int w = 0; w < nw; ++w) {
            const Well& well = wells_ecl[w];
            if (well.getStatus() == Well::Status::SHUT)
                continue;

            if ( !well.isMultiSegment() )
                continue;

            const auto& wname = well.name();
            if (prev_well_state->segment_state.has(wname)) {
                if (prev_well_state->status_[wname] == Well::Status::SHUT) {
                    continue;
                }

                // TODO: the well with same name can change a lot, like they might not have same number of segments
                // we need to handle that later.
                // for now, we just copy them.
                this->segment_state.copy_welldata(prev_well_state->segment_state, wname);
            }
        }
    }
}

void
WellState::calculateSegmentRates(const std::vector<std::vector<int>>& segment_inlets,
                                                      const std::vector<std::vector<int>>&segment_perforations,
                                                      const std::vector<double>& perforation_rates,
                                                      const int np, const int segment,
                                                      std::vector<double>& segment_rates)
{
    // the rate of the segment equals to the sum of the contribution from the perforations and inlet segment rates.
    // the first segment is always the top segment, its rates should be equal to the well rates.
    assert(segment_inlets.size() == segment_perforations.size());
    const int well_nseg = segment_inlets.size();
    if (segment == 0) { // beginning the calculation
        segment_rates.resize(np * well_nseg, 0.0);
    }
    // contributions from the perforations belong to this segment
    for (const int& perf : segment_perforations[segment]) {
        for (int p = 0; p < np; ++p) {
            segment_rates[np * segment + p] += perforation_rates[np * perf + p];
        }
    }
    for (const int& inlet_seg : segment_inlets[segment]) {
        calculateSegmentRates(segment_inlets, segment_perforations, perforation_rates, np, inlet_seg, segment_rates);
        for (int p = 0; p < np; ++p) {
            segment_rates[np * segment + p] += segment_rates[np * inlet_seg + p];
        }
    }
}

double WellState::solventWellRate(const int w) const
{
    const auto& perf_rates_solvent = this->perfRateSolvent(w);
    return parallel_well_info_[w]->sumPerfValues(perf_rates_solvent.begin(), perf_rates_solvent.end());
}

double WellState::polymerWellRate(const int w) const
{
    const auto& perf_rates_polymer = this->perfRatePolymer(w);
    return parallel_well_info_[w]->sumPerfValues(perf_rates_polymer.begin(), perf_rates_polymer.end());
}

double WellState::brineWellRate(const int w) const
{
    const auto& perf_rates_brine = this->perfRateBrine(w);
    return parallel_well_info_[w]->sumPerfValues(perf_rates_brine.begin(), perf_rates_brine.end());
}


void WellState::stopWell(int well_index)
{
    this->status_[well_index] = Well::Status::STOP;
    this->thp_[well_index] = 0;
}

void WellState::shutWell(int well_index)
{
    this->status_[well_index] = Well::Status::SHUT;
    this->thp_[well_index] = 0;
    this->bhp_[well_index] = 0;
    const int np = numPhases();
    this->wellrates_[well_index].assign(np, 0);

    auto& resv = this->well_reservoir_rates_[well_index];
    auto& wpi  = this->productivity_index_[well_index];

    for (int p = 0; p < np; ++p) {
        resv[p] = 0.0;
        wpi[p]  = 0.0;
    }

    auto& connpi = this->conn_productivity_index_[well_index];
    connpi.assign(connpi.size(), 0);
}

void WellState::updateStatus(int well_index, Well::Status status)
{
    switch (status) {
    case Well::Status::OPEN:
        this->openWell(well_index);
        break;
    case Well::Status::SHUT:
        this->shutWell(well_index);
        break;
    case Well::Status::STOP:
        this->stopWell(well_index);
        break;
    default:
        throw std::logic_error("Invalid well status");
    }
}



template<class Comm>
void WellState::communicateGroupRates(const Comm& comm)
{
    // Compute the size of the data.
    std::size_t sz = 0;
    for (const auto& [_, owner_rates] : this->well_rates) {
        (void)_;
        const auto& [__, rates] = owner_rates;
        (void)__;
        sz += rates.size();
    }
    sz += this->alq_state.pack_size();


    // Make a vector and collect all data into it.
    std::vector<double> data(sz);
    std::size_t pos = 0;
    for (const auto& [_, owner_rates] : this->well_rates) {
        (void)_;
        const auto& [owner, rates] = owner_rates;
        for (const auto& value : rates) {
            if (owner)
                data[pos++] = value;
            else
                data[pos++] = 0;
        }
    }
    pos += this->alq_state.pack_data(&data[pos]);
    assert(pos == sz);

    // Communicate it with a single sum() call.
    comm.sum(data.data(), data.size());

    pos = 0;
    for (auto& [_, owner_rates] : this->well_rates) {
        (void)_;
        auto& [__, rates] = owner_rates;
        (void)__;
        for (auto& value : rates)
            value = data[pos++];
    }
    pos += this->alq_state.unpack_data(&data[pos]);
    assert(pos == sz);
}


template<class Comm>
void WellState::updateGlobalIsGrup(const Comm& comm)
{
    this->global_well_info.value().update_group(this->status_.data(), this->current_injection_controls_.data(), this->current_production_controls_.data());
    this->global_well_info.value().communicate(comm);
}

data::Segment
WellState::reportSegmentResults(const PhaseUsage& pu,
                                                     const int         well_id,
                                                     const int         seg_ix,
                                                     const int         seg_no) const
{
    const auto& segments = this->segments(well_id);
    if (segments.empty())
        return {};

    auto seg_res = data::Segment{};
    {
        using Value = data::SegmentPressures::Value;
        auto& segpress = seg_res.pressures;
        segpress[Value::Pressure] = segments.pressure[seg_ix];
        segpress[Value::PDrop] = segments.pressure_drop(seg_ix);
        segpress[Value::PDropHydrostatic] = segments.pressure_drop_hydrostatic[seg_ix];
        segpress[Value::PDropFriction] = segments.pressure_drop_friction[seg_ix];
        segpress[Value::PDropAccel] = segments.pressure_drop_accel[seg_ix];
    }

    const auto rate = &segments.rates[seg_ix * pu.num_phases];
    if (pu.phase_used[Water]) {
        seg_res.rates.set(data::Rates::opt::wat,
                          rate[pu.phase_pos[Water]]);
    }

    if (pu.phase_used[Oil]) {
        seg_res.rates.set(data::Rates::opt::oil,
                          rate[pu.phase_pos[Oil]]);
    }

    if (pu.phase_used[Gas]) {
        seg_res.rates.set(data::Rates::opt::gas,
                          rate[pu.phase_pos[Gas]]);
    }

    seg_res.segNumber = seg_no;

    return seg_res;
}

bool WellState::wellIsOwned(std::size_t well_index,
                                                 [[maybe_unused]] const std::string& wellName) const
{
    const auto& well_info = parallelWellInfo(well_index);
    assert(well_info.name() == wellName);

    return well_info.isOwner();
}

bool WellState::wellIsOwned(const std::string& wellName) const
{
    const auto& it = this->wellMap_.find( wellName );
    if (it == this->wellMap_.end()) {
        OPM_THROW(std::logic_error, "Could not find well " << wellName << " in well map");
    }
    const int well_index = it->second[0];
    return wellIsOwned(well_index, wellName);
}

int WellState::numSegments(const int well_id) const
{
    const auto& segments = this->segments(well_id);
    return segments.size();
}

int WellState::segmentNumber(const int well_id, const int seg_id) const
{
    const auto& segments = this->segment_state[well_id];
    return segments.segment_number()[seg_id];
}

void WellState::updateWellsDefaultALQ( const std::vector<Well>& wells_ecl )
{
    const int nw = wells_ecl.size();
    for (int i = 0; i<nw; i++) {
        const Well &well = wells_ecl[i];
        if (well.isProducer()) {
            // NOTE: This is the value set in item 12 of WCONPROD, or with WELTARG
            auto alq = well.alq_value();
            this->alq_state.update_default(well.name(), alq);
        }
    }
}

void WellState::resetConnectionTransFactors(const int well_index,
                                                                 const std::vector<PerforationData>& well_perf_data)
{
    if (this->well_perf_data_[well_index].size() != well_perf_data.size()) {
        throw std::invalid_argument {
            "Size mismatch for perforation data in well "
            + std::to_string(well_index)
        };
    }

    auto connID = std::size_t{0};
    auto dst = this->well_perf_data_[well_index].begin();
    for (const auto& src : well_perf_data) {
        if (dst->cell_index != src.cell_index) {
            throw std::invalid_argument {
                "Cell index mismatch in connection "
                + std::to_string(connID)
                        + " of well "
                        + std::to_string(well_index)
            };
        }

        if (dst->satnum_id != src.satnum_id) {
            throw std::invalid_argument {
                "Saturation function table mismatch in connection "
                + std::to_string(connID)
                        + " of well "
                        + std::to_string(well_index)
            };
        }

        dst->connection_transmissibility_factor =
                src.connection_transmissibility_factor;

        ++dst;
        ++connID;
    }
}

const ParallelWellInfo&
WellState::parallelWellInfo(std::size_t well_index) const
{
    return *parallel_well_info_[well_index];
}

template void WellState::updateGlobalIsGrup<ParallelWellInfo::Communication>(const ParallelWellInfo::Communication& comm);
template void WellState::communicateGroupRates<ParallelWellInfo::Communication>(const ParallelWellInfo::Communication& comm);
} // namespace Opm
