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
    wellMap_.clear();

    well_perf_data_ = well_perf_data;
    parallel_well_info_ = parallel_well_info;

    {
        // const int nw = wells->number_of_wells;
        const int nw = wells_ecl.size();
        const int np = this->phase_usage_.num_phases;
        // const int np = wells->number_of_phases;
        status_.assign(nw, Well::Status::OPEN);
        bhp_.resize(nw, 0.0);
        thp_.resize(nw, 0.0);
        temperature_.resize(nw, 273.15 + 15.56); // standard condition temperature
        wellrates_.resize(nw * np, 0.0);
        int connpos = 0;
        for (int w = 0; w < nw; ++w) {
            const Well& well = wells_ecl[w];

            // Initialize bhp(), thp(), wellRates(), temperature().
            initSingleWell(cellPressures, w, well, *parallel_well_info[w], summary_state);

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

        // The perforation rates and perforation pressures are
        // not expected to be consistent with bhp_ and wellrates_
        // after init().
        perfrates_.resize(connpos, 0.0);
        perfpress_.resize(connpos, -1e100);
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

bool WellState::wellIsOwned(std::size_t well_index,
                            [[maybe_unused]] const std::string& wellName) const
{
    const auto& well_info = parallelWellInfo(well_index);
    assert(well_info.name() == wellName);

    return well_info.isOwner();
}

bool WellState::wellIsOwned(const std::string& wellName) const
{
    const auto& it = wellMap().find( wellName );
    if (it == wellMap().end()) {
        OPM_THROW(std::logic_error, "Could not find well " << wellName << " in well map");
    }
    const int well_index = it->second[0];
    return wellIsOwned(well_index, wellName);
}

void WellState::shutWell(int well_index)
{
    this->status_[well_index] = Well::Status::SHUT;
    this->thp_[well_index] = 0;
    this->bhp_[well_index] = 0;
    const int np = numPhases();
    for (int p = 0; p < np; ++p)
        this->wellrates_[np * well_index + p] = 0;
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

data::Wells WellState::report(const int* globalCellIdxMap,
                              const std::function<bool(const int)>& wasDynamicallyClosed) const
{
    using rt = data::Rates::opt;

    const auto& pu = this->phaseUsage();
    data::Wells dw;
    for( const auto& itr : this->wellMap_ ) {
        const auto well_index = itr.second[ 0 ];
        if ((this->status_[well_index] == Well::Status::SHUT) &&
            ! wasDynamicallyClosed(well_index))
        {
            continue;
        }

        const auto& pwinfo = *parallel_well_info_[well_index];
        using WellT = std::remove_reference_t<decltype(dw[ itr.first ])>;
        WellT dummyWell; // dummy if we are not owner
        auto& well = pwinfo.isOwner() ? dw[ itr.first ] : dummyWell;
        well.bhp = this->bhp(well_index);
        well.thp = this->thp( well_index );
        well.temperature = this->temperature().at( well_index );

        const auto wellrate_index = well_index * pu.num_phases;
        const auto& wv = this->wellRates();
        if( pu.phase_used[BlackoilPhases::Aqua] ) {
            well.rates.set( rt::wat, wv[ wellrate_index + pu.phase_pos[BlackoilPhases::Aqua] ] );
        }

        if( pu.phase_used[BlackoilPhases::Liquid] ) {
            well.rates.set( rt::oil, wv[ wellrate_index + pu.phase_pos[BlackoilPhases::Liquid] ] );
        }

        if( pu.phase_used[BlackoilPhases::Vapour] ) {
            well.rates.set( rt::gas, wv[ wellrate_index + pu.phase_pos[BlackoilPhases::Vapour] ] );
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

    return dw;

}

void WellState::reportConnections(data::Well& well,
                                  const PhaseUsage&,
                                  const WellMapType::value_type& itr,
                                  const int* globalCellIdxMap) const
{
    const auto well_index = itr.second[ 0 ];
    const auto& pd = this->well_perf_data_[well_index];
    const int num_perf_well = pd.size();
    well.connections.resize(num_perf_well);

    const auto * perf_rates = &this->perfRates()[itr.second[1]];
    const auto * perf_pressure = &this->perfPress()[itr.second[1]];
    for( int i = 0; i < num_perf_well; ++i ) {
        const auto active_index = this->well_perf_data_[well_index][i].cell_index;
        auto& connection = well.connections[ i ];
        connection.index = globalCellIdxMap[active_index];
        connection.pressure = perf_pressure[i];
        connection.reservoir_rate = perf_rates[i];
        connection.trans_factor = pd[i].connection_transmissibility_factor;
    }
    assert(num_perf_well == int(well.connections.size()));
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

void WellState::initSingleWell(const std::vector<double>& cellPressures,
                               const int w,
                               const Well& well,
                               const ParallelWellInfo& well_info,
                               const SummaryState& summary_state)
{
    assert(well.isInjector() || well.isProducer());

    // Set default zero initial well rates.
    // May be overwritten below.
    const auto& pu = this->phase_usage_;
    const int np = pu.num_phases;
    for (int p = 0; p < np; ++p) {
        wellrates_[np*w + p] = 0.0;
    }

    if ( well.isInjector() ) {
        temperature_[w] = well.injectionControls(summary_state).temperature;
    }

    const int num_perf_this_well = well_info.communication().sum(well_perf_data_[w].size());
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
    const double global_pressure = well_info.broadcastFirstPerforationValue(local_pressure);

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
        if (well.isInjector()) {
            if (inj_controls.cmode == Well::InjectorCMode::RATE) {
                switch (inj_controls.injector_type) {
                case InjectorType::WATER:
                    assert(pu.phase_used[BlackoilPhases::Aqua]);
                    wellrates_[np*w + pu.phase_pos[BlackoilPhases::Aqua]] = inj_surf_rate;
                    break;
                case InjectorType::GAS:
                    assert(pu.phase_used[BlackoilPhases::Vapour]);
                    wellrates_[np*w + pu.phase_pos[BlackoilPhases::Vapour]] = inj_surf_rate;
                    break;
                case InjectorType::OIL:
                    assert(pu.phase_used[BlackoilPhases::Liquid]);
                    wellrates_[np*w + pu.phase_pos[BlackoilPhases::Liquid]] = inj_surf_rate;
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
                wellrates_[np*w + pu.phase_pos[BlackoilPhases::Liquid]] = -prod_controls.oil_rate;
                break;
            case Well::ProducerCMode::WRAT:
                assert(pu.phase_used[BlackoilPhases::Aqua]);
                wellrates_[np*w + pu.phase_pos[BlackoilPhases::Aqua]] = -prod_controls.water_rate;
                break;
            case Well::ProducerCMode::GRAT:
                assert(pu.phase_used[BlackoilPhases::Vapour]);
                wellrates_[np*w + pu.phase_pos[BlackoilPhases::Vapour]] = -prod_controls.gas_rate;
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
