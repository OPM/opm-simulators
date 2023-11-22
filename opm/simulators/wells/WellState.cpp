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
#include <opm/input/eclipse/Schedule/MSW/WellSegments.hpp>
#include <opm/input/eclipse/Schedule/Schedule.hpp>
#include <opm/input/eclipse/Schedule/Well/Well.hpp>
#include <opm/input/eclipse/Schedule/Well/WellConnections.hpp>

#include <opm/simulators/utils/ParallelCommunication.hpp>
#include <opm/simulators/wells/ParallelWellInfo.hpp>
#include <opm/grid/common/p2pcommunicator.hh>
#include <opm/output/data/Wells.hpp>

#include <algorithm>
#include <cassert>
#include <numeric>
#include <set>
#include <stdexcept>
#include <vector>

#include <fmt/format.h>

namespace {

using P2PCommunicatorType = Dune::Point2PointCommunicator<Dune::SimpleMessageBuffer>;
using MessageBufferType = P2PCommunicatorType::MessageBufferType;

class PackUnpackXConn : public P2PCommunicatorType::DataHandleInterface
{
public:
    using XConn = std::vector<Opm::data::Connection>;

    explicit PackUnpackXConn(const bool   isOwner,
                             const XConn& local,
                             XConn&       global);

    // Pack all data associated with link.
    void pack(const int link, MessageBufferType& buffer);

    // Unpack all data associated with link.
    void unpack([[maybe_unused]] const int link,
                MessageBufferType&         buffer);

private:
    const XConn& local_;
    XConn& global_;
};

PackUnpackXConn::PackUnpackXConn(const bool   isOwner,
                                 const XConn& local,
                                 XConn&       global)
    : local_ (local)
    , global_(global)
{
    if (! isOwner) {
        return;
    }

    this->global_.insert(this->global_.end(),
                         this->local_.begin(),
                         this->local_.end());
}

void PackUnpackXConn::pack(const int          link,
                           MessageBufferType& buffer)
{
    // We should only get one link
    if (link != 0) {
        throw std::logic_error {
            "link in pack() does not match expected value 0"
        };
    }

    // Write all local connection results
    {
        const auto nconn = this->local_.size();
        buffer.write(nconn);
    }

    for (const auto& conn : this->local_) {
        conn.write(buffer);
    }
}

void PackUnpackXConn::unpack([[maybe_unused]] const int link,
                             MessageBufferType&         buffer)
{
    const auto nconn = [this, &buffer]()
    {
        auto nc = 0 * this->local_.size();
        buffer.read(nc);

        return nc;
    }();

    this->global_.reserve(this->global_.size() + nconn);

    for (auto conn = 0*nconn; conn < nconn; ++conn) {
        this->global_.emplace_back().read(buffer);
    }
}

} // Anonymous namespace

namespace Opm {

WellState::WellState(const ParallelWellInfo& pinfo)
    : phase_usage_{}
{
    wells_.add("test4",
               SingleWellState{"dummy", pinfo, false, 0.0, {}, phase_usage_, 0.0});
}

WellState WellState::serializationTestObject(const ParallelWellInfo& pinfo)
{
    WellState result(PhaseUsage{});
    result.alq_state = ALQState::serializationTestObject();
    result.well_rates = {{"test2", {true, {1.0}}}, {"test3", {false, {2.0}}}};
    result.wells_.add("test4", SingleWellState::serializationTestObject(pinfo));

    return result;
}

void WellState::base_init(const std::vector<double>& cellPressures,
                          const std::vector<Well>& wells_ecl,
                          const std::vector<std::reference_wrapper<ParallelWellInfo>>& parallel_well_info,
                          const std::vector<std::vector<PerforationData>>& well_perf_data,
                          const SummaryState& summary_state)
{
    // clear old name mapping
    this->wells_.clear();
    {
        // const int nw = wells->number_of_wells;
        const int nw = wells_ecl.size();
        // const int np = wells->number_of_phases;
        for (int w = 0; w < nw; ++w) {
            const Well& well = wells_ecl[w];

            // Initialize bhp(), thp(), wellRates(), temperature().
            initSingleWell(cellPressures, well, well_perf_data[w], parallel_well_info[w], summary_state);
        }
    }
}

void WellState::initSingleProducer(const Well& well,
                                   const ParallelWellInfo& well_info,
                                   double pressure_first_connection,
                                   const std::vector<PerforationData>& well_perf_data,
                                   const SummaryState& summary_state) {
    const auto& pu = this->phase_usage_;
    const double temp = 273.15 + 15.56;

    auto& ws = this->wells_.add(well.name(),
                                SingleWellState{well.name(),
                                                well_info,
                                                true,
                                                pressure_first_connection,
                                                well_perf_data,
                                                pu,
                                                temp});

    // the rest of the code needs to executed even if ws.perf_data is empty
    // as this does not say anything for the whole well if it is distributed.
    // Hence never ever return here!
    if (well.getStatus() == Well::Status::OPEN) {
        ws.status = Well::Status::OPEN;
    }

    ws.update_producer_targets(well, summary_state);
}

void WellState::initSingleInjector(const Well& well,
                                   const ParallelWellInfo& well_info,
                                   double pressure_first_connection,
                                   const std::vector<PerforationData>& well_perf_data,
                                   const SummaryState& summary_state) {

    const auto& pu = this->phase_usage_;
    const double temp = well.temperature();

    auto& ws = this->wells_.add(well.name(), SingleWellState{well.name(),
                                                             well_info,
                                                             false,
                                                             pressure_first_connection,
                                                             well_perf_data,
                                                             pu,
                                                             temp});

    // the rest of the code needs to executed even if ws.perf_data is empty
    // as this does not say anything for the whole well if it is distributed.
    // Hence never ever return here!
    if (well.getStatus() == Well::Status::OPEN) {
        ws.status = Well::Status::OPEN;
    }

    ws.update_injector_targets(well, summary_state);
}

void WellState::initSingleWell(const std::vector<double>& cellPressures,
                               const Well& well,
                               const std::vector<PerforationData>& well_perf_data,
                               const ParallelWellInfo& well_info,
                               const SummaryState& summary_state)
{
    double pressure_first_connection = -1;
    if (!well_perf_data.empty())
        pressure_first_connection = cellPressures[well_perf_data[0].cell_index];
    pressure_first_connection = well_info.broadcastFirstPerforationValue(pressure_first_connection);

    if (well.isInjector()) {
        this->initSingleInjector(well, well_info, pressure_first_connection,
                                 well_perf_data, summary_state);
    } else {
        this->initSingleProducer(well, well_info, pressure_first_connection,
                                 well_perf_data, summary_state);
    }
}

void WellState::init(const std::vector<double>& cellPressures,
                     const Schedule& schedule,
                     const std::vector<Well>& wells_ecl,
                     const std::vector<std::reference_wrapper<ParallelWellInfo>>& parallel_well_info,
                     const int report_step,
                     const WellState* prevState,
                     const std::vector<std::vector<PerforationData>>& well_perf_data,
                     const SummaryState& summary_state)
{
    // call init on base class
    this->base_init(cellPressures, wells_ecl, parallel_well_info,
                    well_perf_data, summary_state);
    this->global_well_info = std::make_optional<GlobalWellInfo>(schedule,
                                                                report_step,
                                                                wells_ecl);
    for (const auto& wname : schedule.wellNames(report_step))
    {
        well_rates.insert({wname, std::make_pair(false, std::vector<double>(this->numPhases()))});
    }
    for (const auto& winfo: parallel_well_info)
    {
        well_rates[winfo.get().name()].first = winfo.get().isOwner();
    }

    const int nw = wells_ecl.size();

    if( nw == 0 ) return ;

    // Initialize perfphaserates_, which must be done here.
    const auto& pu = this->phaseUsage();
    const int np = pu.num_phases;

    {
        const auto& wg_events = schedule[report_step].wellgroup_events();
        for (const auto& ecl_well : wells_ecl) {
            const auto& wname = ecl_well.name();
            if (wg_events.has(wname))
                this->well(wname).events = wg_events.at(wname);
        }
    }

    for (int w = 0; w < nw; ++w) {
        // Initialize perfphaserates_ to well
        // rates divided by the number of perforations.
        const auto& ecl_well = wells_ecl[w];
        auto& ws = this->well(w);
        auto& perf_data = ws.perf_data;
        const int num_perf_this_well = perf_data.size();
        const int global_num_perf_this_well = ecl_well.getConnections().num_open();

        for (int perf = 0; perf < num_perf_this_well; ++perf) {
            if (wells_ecl[w].getStatus() == Well::Status::OPEN) {
                for (int p = 0; p < this->numPhases(); ++p) {
                    perf_data.phase_rates[this->numPhases()*perf + p] = ws.surface_rates[p] / double(global_num_perf_this_well);
                }
            }
            perf_data.pressure[perf] = cellPressures[well_perf_data[w][perf].cell_index];
        }
    }

    for (int w = 0; w < nw; ++w) {
        auto& ws = this->well(w);
        if (wells_ecl[w].isProducer()) {
            const auto controls = wells_ecl[w].productionControls(summary_state);
            if (controls.cmode == Well::ProducerCMode::GRUP && !wells_ecl[w].isAvailableForGroupControl()) {
                ws.production_cmode = Well::ProducerCMode::BHP; // wells always has a BHP control
            } else {
                ws.production_cmode = controls.cmode;
            }
        }
        else {
            const auto controls = wells_ecl[w].injectionControls(summary_state);
            if (controls.cmode == Well::InjectorCMode::GRUP && !wells_ecl[w].isAvailableForGroupControl()) {
                ws.injection_cmode = Well::InjectorCMode::BHP; // wells always has a BHP control
            } else {
                ws.injection_cmode = controls.cmode;
            }
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
    if (prevState && prevState->size() > 0) {
        for (int w = 0; w < nw; ++w) {
            const Well& well = wells_ecl[w];
            if (well.getStatus() == Well::Status::SHUT) {
                continue;
            }
            auto& new_well = this->well(w);
            const auto& old_index = prevState->index(well.name());
            if (old_index.has_value()) {
                const auto& prev_well = prevState->well(old_index.value());
                new_well.init_timestep(prev_well);


                if (prev_well.status == Well::Status::SHUT) {
                    // Well was shut in previous state, do not use its values.
                    continue;
                }

                if (new_well.producer != prev_well.producer)
                    // Well changed to/from injector from/to producer, do not use its privious values.
                    continue;

                // If new target is set using WCONPROD, WCONINJE etc. we use the new control
                if (!new_well.events.hasEvent(WellState::event_mask)) {
                    new_well.injection_cmode = prev_well.injection_cmode;
                    new_well.production_cmode = prev_well.production_cmode;
                }

                new_well.surface_rates = prev_well.surface_rates;
                new_well.reservoir_rates = prev_well.reservoir_rates;
                new_well.well_potentials = prev_well.well_potentials;

                // perfPhaseRates
                const int num_perf_old_well = prev_well.perf_data.size();
                const int num_perf_this_well = new_well.perf_data.size();
                const bool global_num_perf_same = (num_perf_this_well == num_perf_old_well);

                // copy perforation rates when the number of
                // perforations is equal, otherwise initialize
                // perfphaserates to well rates divided by the
                // number of perforations.
                // TODO: we might still need the values from the prev_well if the connection structure changes
                if (global_num_perf_same)
                {
                    auto& perf_data = new_well.perf_data;
                    const auto& prev_perf_data = prev_well.perf_data;
                    perf_data.try_assign( prev_perf_data );
                } else {
                    const int global_num_perf_this_well = well.getConnections().num_open();
                    auto& perf_data = new_well.perf_data;
                    auto& target_rates = perf_data.phase_rates;
                    for (int perf_index = 0; perf_index < num_perf_this_well; perf_index++) {
                        for (int p = 0; p < np; ++p) {
                            target_rates[perf_index*np + p] = new_well.surface_rates[p] / double(global_num_perf_this_well);
                        }
                    }
                }

                // Productivity index.
                new_well.productivity_index = prev_well.productivity_index;

                // if there is no valid VFP table associated, we set the THP value to be 0.
                if (well.vfp_table_number() == 0) {
                    new_well.thp = 0.;
                }
            }
        }
    }


    updateWellsDefaultALQ(wells_ecl);
}

void WellState::resize(const std::vector<Well>& wells_ecl,
                       const std::vector<std::reference_wrapper<ParallelWellInfo>>& parallel_well_info,
                       const Schedule& schedule,
                       const bool handle_ms_well,
                       const std::size_t numCells,
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
        OPM_THROW(std::logic_error,
                  "Could not find any rates for well " + wellName);

    return it->second.second;
}

template<class Communication>
void WellState::gatherVectorsOnRoot(const std::vector<data::Connection>& from_connections,
                                                         std::vector<data::Connection>& to_connections,
                                                         const Communication& comm) const
{
    auto send = std::set<int>{};
    auto recv = std::set<int>{};

    const auto isOwner = comm.rank() == 0;
    if (isOwner) {
        for (auto other = 0*comm.size() + 1; other < comm.size(); ++other) {
            recv.insert(other);
        }
    }
    else {
        send.insert(0);
    }

    auto toOwnerComm = P2PCommunicatorType{ comm };
    toOwnerComm.insertRequest(send, recv);

    PackUnpackXConn lineariser { isOwner, from_connections, to_connections };
    toOwnerComm.exchange(lineariser);
}

data::Wells
WellState::report(const int* globalCellIdxMap,
                  const std::function<bool(const int)>& wasDynamicallyClosed) const
{
    if (this->numWells() == 0) {
        return {};
    }

    using rt = data::Rates::opt;
    const auto& pu = this->phaseUsage();

    data::Wells res;
    for (std::size_t well_index = 0; well_index < this->size(); ++well_index) {
        const auto& ws = this->well(well_index);
        if ((ws.status == Well::Status::SHUT) && !wasDynamicallyClosed(well_index))
        {
            continue;
        }

        const auto& reservoir_rates = ws.reservoir_rates;
        const auto& well_potentials = ws.well_potentials;
        const auto& wpi = ws.productivity_index;
        const auto& wv = ws.surface_rates;
        const auto& wname = this->name(well_index);

        auto dummyWell = data::Well{};
        auto& well = ws.parallel_info.get().isOwner() ? res[wname] : dummyWell;

        well.bhp = ws.bhp;
        well.thp = ws.thp;
        well.temperature = ws.temperature;
        well.filtrate.rate = ws.sum_filtrate_rate();
        well.filtrate.total = ws.sum_filtrate_total();
        well.filtrate.concentration = ws.filtrate_conc;

        if (pu.phase_used[BlackoilPhases::Aqua]) {
            well.rates.set(rt::wat, wv[ pu.phase_pos[BlackoilPhases::Aqua] ] );
            well.rates.set(rt::reservoir_water, reservoir_rates[pu.phase_pos[BlackoilPhases::Aqua]]);
            well.rates.set(rt::productivity_index_water, wpi[pu.phase_pos[BlackoilPhases::Aqua]]);
            well.rates.set(rt::well_potential_water, well_potentials[pu.phase_pos[BlackoilPhases::Aqua]]);
        }

        if (pu.phase_used[BlackoilPhases::Liquid]) {
            well.rates.set(rt::oil, wv[ pu.phase_pos[BlackoilPhases::Liquid] ] );
            well.rates.set(rt::reservoir_oil, reservoir_rates[pu.phase_pos[BlackoilPhases::Liquid]]);
            well.rates.set(rt::productivity_index_oil, wpi[pu.phase_pos[BlackoilPhases::Liquid]]);
            well.rates.set(rt::well_potential_oil, well_potentials[pu.phase_pos[BlackoilPhases::Liquid]]);
        }

        if( pu.phase_used[BlackoilPhases::Vapour] ) {
            well.rates.set(rt::gas, wv[ pu.phase_pos[BlackoilPhases::Vapour] ] );
            well.rates.set(rt::reservoir_gas, reservoir_rates[pu.phase_pos[BlackoilPhases::Vapour]]);
            well.rates.set(rt::productivity_index_gas, wpi[pu.phase_pos[BlackoilPhases::Vapour]]);
            well.rates.set(rt::well_potential_gas, well_potentials[pu.phase_pos[BlackoilPhases::Vapour]]);
        }

        if (pu.has_solvent || pu.has_zFraction) {
            well.rates.set(rt::solvent, ws.sum_solvent_rates());
        }

        if (pu.has_polymer) {
            well.rates.set(rt::polymer, ws.sum_polymer_rates());
        }

        if (pu.has_brine) {
            well.rates.set(rt::brine, ws.sum_brine_rates());
        }

        if (ws.producer) {
            well.rates.set(rt::alq, getALQ(wname));
        }
        else {
            well.rates.set(rt::alq, 0.0);
        }

        well.rates.set(rt::dissolved_gas,
                       ws.phase_mixing_rates[ws.dissolved_gas] +
                       ws.phase_mixing_rates[ws.dissolved_gas_in_water]);
        well.rates.set(rt::vaporized_oil, ws.phase_mixing_rates[ws.vaporized_oil]);
        well.rates.set(rt::vaporized_water, ws.phase_mixing_rates[ws.vaporized_water]);

        {
            auto& curr = well.current_control;

            curr.isProducer = ws.producer;
            curr.prod = ws.production_cmode;
            curr.inj  = ws.injection_cmode;
        }

        const auto& pwinfo = ws.parallel_info.get();
        if (pwinfo.communication().size() == 1) {
            reportConnections(well.connections, pu, well_index, globalCellIdxMap);
        }
        else {
            std::vector<data::Connection> connections;

            reportConnections(connections, pu, well_index, globalCellIdxMap);
            gatherVectorsOnRoot(connections, well.connections, pwinfo.communication());
        }

        const auto nseg = ws.segments.size();
        for (auto seg_ix = 0*nseg; seg_ix < nseg; ++seg_ix) {
            const auto seg_no = ws.segments.segment_number()[seg_ix];
            well.segments[seg_no] = this->reportSegmentResults(well_index, seg_ix, seg_no);
        }
    }

    return res;
}

void WellState::reportConnections(std::vector<data::Connection>& connections,
                                  const PhaseUsage &pu,
                                  std::size_t well_index,
                                  const int* globalCellIdxMap) const
{
    using rt = data::Rates::opt;
    const auto& ws = this->well(well_index);
    const auto& perf_data = ws.perf_data;
    const int num_perf_well = perf_data.size();
    connections.resize(num_perf_well);
    const auto& perf_rates = perf_data.rates;
    const auto& perf_pressure = perf_data.pressure;
    const auto& perf_mixing_rates = perf_data.phase_mixing_rates;
    for (int i = 0; i < num_perf_well; ++i) {
      const auto active_index = perf_data.cell_index[i];
        auto& connection = connections[ i ];
        connection.index = globalCellIdxMap[active_index];
        connection.pressure = perf_pressure[i];
        connection.reservoir_rate = perf_rates[i];
        connection.trans_factor = perf_data.connection_transmissibility_factor[i];
        connection.d_factor = perf_data.connection_d_factor[i];
        connection.rates.set(rt::dissolved_gas, perf_mixing_rates[i][ws.dissolved_gas]);
        connection.rates.set(rt::vaporized_oil, perf_mixing_rates[i][ws.vaporized_oil]);
        if (!ws.producer) {
            const auto& filtrate_data = perf_data.filtrate_data;
            auto& filtrate = connection.filtrate;
            filtrate.rate = filtrate_data.rates[i];
            filtrate.total = filtrate_data.total[i];
            filtrate.skin_factor = filtrate_data.skin_factor[i];
            filtrate.thickness = filtrate_data.thickness[i];
            filtrate.poro = filtrate_data.poro[i];
            filtrate.perm = filtrate_data.perm[i];
            filtrate.radius = filtrate_data.radius[i];
            filtrate.area_of_flow = filtrate_data.area_of_flow[i];
        }
    }

    const int np = pu.num_phases;
    std::vector< rt > phs( np );
    std::vector<rt> pi(np);
    if (pu.phase_used[Water]) {
        phs.at( pu.phase_pos[Water] ) = rt::wat;
        pi .at( pu.phase_pos[Water] ) = rt::productivity_index_water;
    }

    if (pu.phase_used[Oil]) {
        phs.at( pu.phase_pos[Oil] ) = rt::oil;
        pi .at( pu.phase_pos[Oil] ) = rt::productivity_index_oil;
    }

    if (pu.phase_used[Gas]) {
        phs.at( pu.phase_pos[Gas] ) = rt::gas;
        pi .at( pu.phase_pos[Gas] ) = rt::productivity_index_gas;
    }

    std::size_t local_conn_index = 0;
    for (auto& comp : connections) {
        const auto * rates = &perf_data.phase_rates[np * local_conn_index];
        const auto * connPI = &perf_data.prod_index[np * local_conn_index];


        for (int i = 0; i < np; ++i) {
            comp.rates.set( phs[ i ], rates[i] );
            comp.rates.set( pi [ i ], connPI[i] );
        }
        if (pu.has_polymer) {
            const auto& perf_polymer_rate = perf_data.polymer_rates;
            comp.rates.set( rt::polymer, perf_polymer_rate[local_conn_index]);
        }
        if (pu.has_brine) {
            const auto& perf_brine_rate = perf_data.brine_rates;
            comp.rates.set( rt::brine, perf_brine_rate[local_conn_index]);
        }
        if (pu.has_solvent) {
            const auto& perf_solvent_rate = perf_data.solvent_rates;
            comp.rates.set( rt::solvent, perf_solvent_rate[local_conn_index] );
        }
        ++local_conn_index;
    }
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
        auto& ws = this->well(w);

        if (well_ecl.isMultiSegment()) {
            const WellSegments& segment_set = well_ecl.getSegments();
            // assuming the order of the perforations in well_ecl is the same with Wells
            const WellConnections& completion_set = well_ecl.getConnections();
            // number of segment for this single well
            ws.segments = SegmentState{np, segment_set};
            const int well_nseg = segment_set.size();
            int n_activeperf = 0;

            // we need to know for each segment, how many perforation it has and how many segments using it as outlet_segment
            // that is why I think we should use a well model to initialize the WellState here
            std::vector<std::vector<int>> segment_perforations(well_nseg);
            for (std::size_t perf = 0; perf < completion_set.size(); ++perf) {
                const Connection& connection = completion_set.get(perf);
                if (connection.state() == Connection::State::OPEN) {
                    const int segment_index = segment_set.segmentNumberToIndex(connection.segment());
                    if (segment_index == -1) {
                        OPM_THROW(std::logic_error,
                                  fmt::format("COMPSEGS: Well {} has connection in cell {}, {}, {} "
                                              "without associated segment.", well_ecl.name(),
                                              connection.getI() + 1 , connection.getJ() + 1,
                                              connection.getK() + 1 ));
                    }

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

            auto& perf_data = ws.perf_data;
            // for the seg_rates_, now it becomes a recursive solution procedure.
            if (pu.phase_used[Gas]) {
                auto& perf_rates = perf_data.phase_rates;
                const int gaspos = pu.phase_pos[Gas];
                // scale the phase rates for Gas to avoid too bad initial guess for gas fraction
                // it will probably benefit the standard well too, while it needs to be justified
                // TODO: to see if this strategy can benefit StandardWell too
                // TODO: it might cause big problem for gas rate control or if there is a gas rate limit
                // maybe the best way is to initialize the fractions first then get the rates
                for (int perf = 0; perf < n_activeperf; perf++)
                    perf_rates[perf*np + gaspos] *= 100;
            }

            const auto& perf_rates = perf_data.phase_rates;
            std::vector<double> perforation_rates(perf_rates.begin(), perf_rates.end());

            calculateSegmentRates(segment_inlets, segment_perforations, perforation_rates, np, 0 /* top segment */, ws.segments.rates);
            // for the segment pressure, the segment pressure is the same with the first perforation belongs to the segment
            // if there is no perforation associated with this segment, it uses the pressure from the outlet segment
            // which requres the ordering is successful
            // Not sure what is the best way to handle the initialization, hopefully, the bad initialization can be
            // improved during the solveWellEq process
            {
                // top segment is always the first one, and its pressure is the well bhp
                auto& segment_pressure = ws.segments.pressure;
                segment_pressure[0] = ws.bhp;
                const auto& perf_press = perf_data.pressure;
                for (int seg = 1; seg < well_nseg; ++seg) {
                    if (!segment_perforations[seg].empty()) {
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

            if (!well.isMultiSegment())
                continue;

            const auto& wname = well.name();
            if (prev_well_state->has(wname)) {
                auto& ws = this->well(w);
                const auto& prev_ws = prev_well_state->well(wname);
                if (prev_ws.status == Well::Status::SHUT) {
                    continue;
                }

                // we do not copy the segment information if the number of segments have changed
                // safer way should be to check the relevant Event
                if (ws.segments.size() == prev_ws.segments.size()) {
                    ws.segments = prev_ws.segments;
                }
            }
        }
    }
}

void
WellState::calculateSegmentRates(const std::vector<std::vector<int>>& segment_inlets,
                                 const std::vector<std::vector<int>>& segment_perforations,
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

void WellState::stopWell(int well_index)
{
    auto& ws = this->well(well_index);
    ws.stop();
}

void WellState::openWell(int well_index)
{
    auto& ws = this->well(well_index);
    ws.open();
}

void WellState::shutWell(int well_index)
{
    auto& ws = this->well(well_index);
    ws.shut();
}

void WellState::updateStatus(int well_index, WellStatus status)
{
    auto& ws = this->well(well_index);
    ws.updateStatus(status);
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
    this->global_well_info.value().clear();
    for (std::size_t well_index = 0; well_index < this->size(); well_index++) {
        const auto& ws = this->well(well_index);
        if (ws.producer)
            this->global_well_info.value().update_producer(well_index, ws.status, ws.production_cmode);
        else
            this->global_well_info.value().update_injector(well_index, ws.status, ws.injection_cmode);
    }
    this->global_well_info.value().communicate(comm);
}

data::Segment
WellState::reportSegmentResults(const int well_id,
                                const int seg_ix,
                                const int seg_no) const
{
    using PhaseQuant = data::SegmentPhaseQuantity::Item;
    using PhaseDensity = data::SegmentPhaseDensity::Item;

    const auto& segments = this->well(well_id).segments;
    if (segments.empty()) {
        return {};
    }

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

    const auto& pu = this->phaseUsage();
    const auto* rate = &segments.rates[seg_ix * pu.num_phases];
    const auto* resv = &segments.phase_resv_rates[seg_ix * pu.num_phases];
    const auto* velocity = &segments.phase_velocity[seg_ix * pu.num_phases];
    const auto* holdup = &segments.phase_holdup[seg_ix * pu.num_phases];
    const auto* viscosity = &segments.phase_viscosity[seg_ix * pu.num_phases];
    const auto* density = &segments.phase_density[seg_ix * (pu.num_phases + 2)]; // +2 for mixture densities

    if (pu.phase_used[Water]) {
        const auto iw = pu.phase_pos[Water];

        seg_res.rates.set(data::Rates::opt::wat, rate[iw]);
        seg_res.rates.set(data::Rates::opt::reservoir_water, resv[iw]);
        seg_res.velocity.set(PhaseQuant::Water, velocity[iw]);
        seg_res.holdup.set(PhaseQuant::Water, holdup[iw]);
        seg_res.viscosity.set(PhaseQuant::Water, viscosity[iw]);
        seg_res.density.set(PhaseDensity::Water, density[iw]);
    }

    if (pu.phase_used[Oil]) {
        const auto io = pu.phase_pos[Oil];

        seg_res.rates.set(data::Rates::opt::oil, rate[io]);
        seg_res.rates.set(data::Rates::opt::vaporized_oil, segments.vaporized_oil_rate[seg_ix]);
        seg_res.rates.set(data::Rates::opt::reservoir_oil, resv[io]);
        seg_res.velocity.set(PhaseQuant::Oil, velocity[io]);
        seg_res.holdup.set(PhaseQuant::Oil, holdup[io]);
        seg_res.viscosity.set(PhaseQuant::Oil, viscosity[io]);
        seg_res.density.set(PhaseDensity::Oil, density[io]);
    }

    if (pu.phase_used[Gas]) {
        const auto ig = pu.phase_pos[Gas];

        seg_res.rates.set(data::Rates::opt::gas, rate[ig]);
        seg_res.rates.set(data::Rates::opt::dissolved_gas, segments.dissolved_gas_rate[seg_ix]);
        seg_res.rates.set(data::Rates::opt::reservoir_gas, resv[ig]);
        seg_res.velocity.set(PhaseQuant::Gas, velocity[ig]);
        seg_res.holdup.set(PhaseQuant::Gas, holdup[ig]);
        seg_res.viscosity.set(PhaseQuant::Gas, viscosity[ig]);
        seg_res.density.set(PhaseDensity::Gas, density[ig]);
    }

    seg_res.segNumber = seg_no;

    // Recall: The mixture density *without* exponents is stored at offset
    // 'num_phases' and the mixture density *with* exponents is stored at
    // offset 'num_phases + 1'.
    seg_res.density.set(PhaseDensity::Mixture, density[pu.num_phases]);
    seg_res.density.set(PhaseDensity::MixtureWithExponents, density[pu.num_phases + 1]);

    return seg_res;
}

bool WellState::wellIsOwned(std::size_t well_index,
                            [[maybe_unused]] const std::string& wellName) const
{
    const auto& well_info = this->parallelWellInfo(well_index);
    assert(well_info.name() == wellName);

    return well_info.isOwner();
}

bool WellState::wellIsOwned(const std::string& wellName) const
{
    const auto& well_index = this->index(wellName);
    if (!well_index.has_value()) {
        OPM_THROW(std::logic_error,
                  fmt::format("Could not find well {} in well map", wellName));
    }

    return wellIsOwned(well_index.value(), wellName);
}

void WellState::updateWellsDefaultALQ(const std::vector<Well>& wells_ecl)
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

bool WellState::operator==(const WellState& rhs) const
{
    return this->alq_state == rhs.alq_state &&
           this->well_rates == rhs.well_rates &&
           this->wells_ == rhs.wells_;
}

const ParallelWellInfo&
WellState::parallelWellInfo(std::size_t well_index) const
{
    const auto& ws = this->well(well_index);
    return ws.parallel_info;
}

template void WellState::updateGlobalIsGrup<Parallel::Communication>(const Parallel::Communication& comm);
template void WellState::communicateGroupRates<Parallel::Communication>(const Parallel::Communication& comm);

} // namespace Opm
