/*
  Copyright 2012 SINTEF ICT, Applied Mathematics.

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

#include <opm/simulators/wells/ParallelWellInfo.hpp>

#include <cassert>
#include <stdexcept>

namespace Opm
{

void WellState::init(const std::vector<double>& cellPressures,
                     const std::vector<Well>& wells_ecl,
                     const std::vector<ParallelWellInfo*>& parallel_well_info,
                     const std::vector<std::vector<PerforationData>>& well_perf_data,
                     const SummaryState& summary_state)
{
    // clear old name mapping
    this->wellMap_.clear();
    this->perfpress_.clear();
    this->perfrates_.clear();
    this->status_.clear();
    this->well_perf_data_.clear();
    this->parallel_well_info_.clear();
    this->wellrates_.clear();
    {
        // const int nw = wells->number_of_wells;
        const int nw = wells_ecl.size();
        // const int np = wells->number_of_phases;
        bhp_.resize(nw, 0.0);
        thp_.resize(nw, 0.0);
        temperature_.resize(nw, 273.15 + 15.56); // standard condition temperature
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


const ParallelWellInfo&
WellState::parallelWellInfo(std::size_t well_index) const
{
    return *parallel_well_info_[well_index];
}


void WellState::shutWell(int well_index)
{
    this->status_[well_index] = Well::Status::SHUT;
    this->thp_[well_index] = 0;
    this->bhp_[well_index] = 0;
    const int np = numPhases();
    this->wellrates_[well_index].assign(np, 0);
}

void WellState::stopWell(int well_index)
{
    this->status_[well_index] = Well::Status::STOP;
    this->thp_[well_index] = 0;
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

    if ( well.isInjector() ) {
        temperature_[w] = well.injectionControls(summary_state).temperature;
    }
    this->status_.add(well.name(), Well::Status::OPEN);
    this->well_perf_data_.add(well.name(), well_perf_data);
    this->parallel_well_info_.add(well.name(), well_info);
    this->wellrates_.add(well.name(), std::vector<double>(np, 0));

    const int num_perf_this_well = well_info->communication().sum(well_perf_data_[w].size());
    this->perfpress_.add(well.name(), std::vector<double>(num_perf_this_well, -1e100));
    this->perfrates_.add(well.name(), std::vector<double>(num_perf_this_well, 0));
    if ( num_perf_this_well == 0 ) {
        // No perforations of the well. Initialize to zero.
        bhp_[w] = 0.;
        thp_[w] = 0.;
        return;
    }

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
        this->openWell(w);
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
        auto & well_rates = this->wellrates_[w];
        if (well.isInjector()) {
            if (inj_controls.cmode == Well::InjectorCMode::RATE) {
                switch (inj_controls.injector_type) {
                case InjectorType::WATER:
                    assert(pu.phase_used[BlackoilPhases::Aqua]);
                    well_rates[pu.phase_pos[BlackoilPhases::Aqua]] = inj_surf_rate;
                    break;
                case InjectorType::GAS:
                    assert(pu.phase_used[BlackoilPhases::Vapour]);
                    well_rates[pu.phase_pos[BlackoilPhases::Vapour]] = inj_surf_rate;
                    break;
                case InjectorType::OIL:
                    assert(pu.phase_used[BlackoilPhases::Liquid]);
                    well_rates[pu.phase_pos[BlackoilPhases::Liquid]] = inj_surf_rate;
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
                well_rates[pu.phase_pos[BlackoilPhases::Liquid]] = -prod_controls.oil_rate;
                break;
            case Well::ProducerCMode::WRAT:
                assert(pu.phase_used[BlackoilPhases::Aqua]);
                well_rates[pu.phase_pos[BlackoilPhases::Aqua]] = -prod_controls.water_rate;
                break;
            case Well::ProducerCMode::GRAT:
                assert(pu.phase_used[BlackoilPhases::Vapour]);
                well_rates[pu.phase_pos[BlackoilPhases::Vapour]] = -prod_controls.gas_rate;
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

} // namespace Opm
