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

#include <opm/output/data/Wells.hpp>

#include <opm/simulators/wells/ConnFracStatistics.hpp>
#include <opm/simulators/wells/ParallelWellInfo.hpp>
#include <opm/simulators/wells/PerforationData.hpp>
#include <opm/simulators/wells/RunningStatistics.hpp>

#include <opm/simulators/utils/ParallelCommunication.hpp>

#include <opm/grid/common/p2pcommunicator.hh>

#include <algorithm>
#include <cassert>
#include <initializer_list>
#include <numeric>
#include <optional>
#include <set>
#include <stdexcept>
#include <type_traits>
#include <utility>
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

template<class Scalar>
WellState<Scalar>::WellState(const ParallelWellInfo<Scalar>& pinfo)
    : phase_usage_{{BlackoilPhases::Aqua, BlackoilPhases::Liquid}}
{
    wells_.add("test4",
               SingleWellState<Scalar>{"dummy", pinfo, false, 0.0, {}, phase_usage_, 0.0});
}

template<class Scalar>
WellState<Scalar> WellState<Scalar>::
serializationTestObject(const ParallelWellInfo<Scalar>& pinfo)
{
    WellState result(PhaseUsage{});
    result.well_rates = {{"test2", {true, {1.0}}}, {"test3", {false, {2.0}}}};
    result.wells_.add("test4", SingleWellState<Scalar>::serializationTestObject(pinfo));

    return result;
}

template<class Scalar>
void WellState<Scalar>::base_init(const std::vector<Scalar>& cellPressures,
                                  const std::vector<Scalar>& cellTemperatures,
                                  const std::vector<Well>& wells_ecl,
                                  const std::vector<std::reference_wrapper<ParallelWellInfo<Scalar>>>& parallel_well_info,
                                  const std::vector<std::vector<PerforationData<Scalar>>>& well_perf_data,
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
            initSingleWell(cellPressures, cellTemperatures, well, well_perf_data[w], parallel_well_info[w], summary_state);
        }
    }
}

template<class Scalar>
void WellState<Scalar>::initSingleProducer(const Well& well,
                                           const ParallelWellInfo<Scalar>& well_info,
                                           Scalar pressure_first_connection,
                                           const std::vector<PerforationData<Scalar>>& well_perf_data,
                                           const SummaryState& summary_state)
{
    const auto& pu = this->phase_usage_;
    const Scalar temp = 273.15 + 15.56;

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

template<class Scalar>
void WellState<Scalar>::initSingleInjector(const Well& well,
                                           const ParallelWellInfo<Scalar>& well_info,
                                           Scalar pressure_first_connection,
                                           Scalar temperature_first_connection,
                                           const std::vector<PerforationData<Scalar>>& well_perf_data,
                                           const SummaryState& summary_state)
{
    const auto& pu = this->phase_usage_;
    const Scalar temp = well.hasInjTemperature() ? well.inj_temperature() : temperature_first_connection;
    auto& ws = this->wells_.add(well.name(), SingleWellState<Scalar>{well.name(),
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

template<class Scalar>
void WellState<Scalar>::initSingleWell(const std::vector<Scalar>& cellPressures,
                                       const std::vector<Scalar>& cellTemperatures,
                                       const Well& well,
                                       const std::vector<PerforationData<Scalar>>& well_perf_data,
                                       const ParallelWellInfo<Scalar>& well_info,
                                       const SummaryState& summary_state)
{
    Scalar pressure_first_connection = -1;
    if (!well_perf_data.empty()) {
        pressure_first_connection = cellPressures[well_perf_data[0].cell_index];
    }
    // The following call is necessary to ensure that processes that do not contain the first perforation get the correct value
    pressure_first_connection = well_info.broadcastFirstPerforationValue(pressure_first_connection);

    if (well.isInjector()) {
        Scalar temperature_first_connection = -1;
        if (!well_perf_data.empty()) {
            temperature_first_connection = cellTemperatures[well_perf_data[0].cell_index];
        }
        // The following call is necessary to ensure that processes that do not contain the first perforation get the correct value
        temperature_first_connection = well_info.broadcastFirstPerforationValue(temperature_first_connection);
        this->initSingleInjector(well, well_info, pressure_first_connection, temperature_first_connection,
                                 well_perf_data, summary_state);
    } else {
        this->initSingleProducer(well, well_info, pressure_first_connection,
                                 well_perf_data, summary_state);
    }
}

template<class Scalar>
void WellState<Scalar>::init(const std::vector<Scalar>& cellPressures,
                             const std::vector<Scalar>& cellTemperatures,
                             const Schedule& schedule,
                             const std::vector<Well>& wells_ecl,
                             const std::vector<std::reference_wrapper<ParallelWellInfo<Scalar>>>& parallel_well_info,
                             const int report_step,
                             const WellState* prevState,
                             const std::vector<std::vector<PerforationData<Scalar>>>& well_perf_data,
                             const SummaryState& summary_state,
                             const bool enableDistributedWells)
{
    // Call init on base class.
    this->base_init(cellPressures, cellTemperatures,
                    wells_ecl, parallel_well_info,
                    well_perf_data, summary_state);

    this->enableDistributedWells_ = enableDistributedWells;

    this->global_well_info.emplace(schedule, report_step, wells_ecl);

    well_rates.clear();

    this->permanently_inactive_well_names_ = schedule.getInactiveWellNamesAtEnd();

    for (const auto& wname : schedule.wellNames(report_step)) {
        well_rates.insert({wname, std::make_pair(false, std::vector<Scalar>(this->numPhases()))});
    }

    for (const auto& winfo : parallel_well_info) {
        well_rates[winfo.get().name()].first = winfo.get().isOwner();
    }

    if (wells_ecl.empty()) {
        return;
    }

    const int nw = wells_ecl.size();

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
                    perf_data.phase_rates[this->numPhases()*perf + p] = ws.surface_rates[p] / Scalar(global_num_perf_this_well);
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
    if ((prevState != nullptr) && (prevState->size() > 0)) {
        for (int w = 0; w < nw; ++w) {
            if (wells_ecl[w].getStatus() == Well::Status::SHUT) {
                continue;
            }

            const auto old_index = prevState->index(wells_ecl[w].name());
            if (! old_index.has_value()) {
                continue;
            }

            const auto& prev_well = prevState->well(*old_index);

            auto& new_well = this->well(w);
            new_well.init_timestep(prev_well);

            if (prev_well.status == Well::Status::SHUT) {
                // Well was shut in previous state, do not use its values.
                continue;
            }

            if (new_well.producer != prev_well.producer) {
                // Well changed to/from injector from/to producer, do not
                // use its previous values.
                continue;
            }

            // If new target is set using WCONPROD, WCONINJE etc. we use the new control
            if (!new_well.events.hasEvent(WellState::event_mask)) {
                new_well.injection_cmode = prev_well.injection_cmode;
                new_well.production_cmode = prev_well.production_cmode;
            }

            new_well.surface_rates = prev_well.surface_rates;
            new_well.prev_surface_rates = prev_well.prev_surface_rates;
            new_well.reservoir_rates = prev_well.reservoir_rates;
            new_well.well_potentials = prev_well.well_potentials;
            new_well.group_target = prev_well.group_target;

            // perfPhaseRates
            //
            // Copy perforation rates when the number of perforations is
            // equal, otherwise initialize perfphaserates to well rates
            // divided by the number of perforations.
            //
            // TODO: we might still need the values from the prev_well if
            // the connection structure changes.
            if (const auto num_perf_this_well = new_well.perf_data.size();
                num_perf_this_well == prev_well.perf_data.size())
            {
                new_well.perf_data.try_assign(prev_well.perf_data);
            }
            else {
                const auto global_num_perf_this_well =
                    static_cast<Scalar>(wells_ecl[w].getConnections().num_open());

                auto target_rate = new_well.perf_data.phase_rates.begin();
                for (auto perf_index = 0*num_perf_this_well; perf_index < num_perf_this_well; ++perf_index) {
                    for (int p = 0; p < np; ++p, ++target_rate) {
                        *target_rate = new_well.surface_rates[p] / global_num_perf_this_well;
                    }
                }
            }

            // Productivity index.
            new_well.productivity_index = prev_well.productivity_index;

            // If there is no valid VFP table associated, we set the THP
            // value to zero.
            if (wells_ecl[w].vfp_table_number() == 0) {
                new_well.thp = Scalar{};
            }
        }
    }

    updateWellsDefaultALQ(schedule, report_step, summary_state);
}

template<class Scalar>
void WellState<Scalar>::resize(const std::vector<Well>& wells_ecl,
                               const std::vector<std::reference_wrapper<ParallelWellInfo<Scalar>>>& parallel_well_info,
                               const Schedule& schedule,
                               const bool handle_ms_well,
                               const std::size_t numCells,
                               const std::vector<std::vector<PerforationData<Scalar>>>& well_perf_data,
                               const SummaryState& summary_state,
                               const bool enable_distributed_wells)
{
    this->enableDistributedWells_ = enable_distributed_wells;
    const std::vector<Scalar> tmp(numCells, 0.0); // <- UGLY HACK to pass the size
    init(tmp, tmp, schedule, wells_ecl, parallel_well_info, 0, nullptr, well_perf_data, summary_state, this->enableDistributedWells_);

    if (handle_ms_well) {
        initWellStateMSWell(wells_ecl, nullptr);
    }
}

template<class Scalar>
const std::vector<Scalar>&
WellState<Scalar>::currentWellRates(const std::string& wellName) const
{
    auto it = well_rates.find(wellName);

    if (it == well_rates.end())
        OPM_THROW(std::logic_error,
                  "Could not find any rates for well " + wellName);

    return it->second.second;
}

template<class Scalar>
void WellState<Scalar>::
gatherVectorsOnRoot(const std::vector<data::Connection>& from_connections,
                    std::vector<data::Connection>& to_connections,
                    const Parallel::Communication& comm) const
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

template<class Scalar>
data::Wells
WellState<Scalar>::report(const int* globalCellIdxMap,
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
        well.efficiency_scaling_factor = ws.efficiency_scaling_factor;
        well.filtrate.rate = ws.sum_filtrate_rate();
        well.filtrate.total = ws.sum_filtrate_total();
        well.filtrate.concentration = ws.filtrate_conc;

        if (pu.phase_used[BlackoilPhases::Aqua]) {
            well.rates.set(rt::wat, wv[ pu.phase_pos[BlackoilPhases::Aqua] ] );
            well.rates.set(rt::reservoir_water, reservoir_rates[pu.phase_pos[BlackoilPhases::Aqua]]);
            well.rates.set(rt::productivity_index_water, wpi[pu.phase_pos[BlackoilPhases::Aqua]]);
            well.rates.set(rt::well_potential_water, well_potentials[pu.phase_pos[BlackoilPhases::Aqua]]);
            well.rates.set(rt::mass_wat, ws.sum_wat_mass_rates());
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

        if (pu.has_micp) {
            well.rates.set(rt::microbial, ws.sum_microbial_rates());
            well.rates.set(rt::oxygen, ws.sum_oxygen_rates());
            well.rates.set(rt::urea, ws.sum_urea_rates());
        }

        if (ws.producer) {
            well.rates.set(rt::alq, ws.alq_state.get());
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

        if (const auto& pwinfo = ws.parallel_info.get();
            pwinfo.communication().size() == 1)
        {
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

template<class Scalar>
void WellState<Scalar>::reportConnections(std::vector<data::Connection>& connections,
                                          const PhaseUsage& pu,
                                          const std::size_t well_index,
                                          const int* globalCellIdxMap) const
{
    const auto& ws = this->well(well_index);

    const auto& perf_data = ws.perf_data;
    const auto  num_perf_well = perf_data.size();

    connections.resize(num_perf_well);

    for (auto i = 0*num_perf_well; i < num_perf_well; ++i) {
        connections[i].index = globalCellIdxMap[perf_data.cell_index[i]];
    }

    this->reportConnectionFactors(well_index, connections);
    this->reportConnectionPressuresAndRates(well_index, pu, connections);

    if (! ws.producer) {
        this->reportConnectionFilterCake(well_index, connections);
    }

    if (! perf_data.connFracStatistics.empty()) {
        this->reportFractureStatistics(perf_data.connFracStatistics,
                                       connections);
    }
}

template<class Scalar>
void WellState<Scalar>::initWellStateMSWell(const std::vector<Well>& wells_ecl,
                                            const WellState* prev_well_state)
{
    // still using the order in wells
    const int nw = wells_ecl.size();
    if (nw == 0) {
        return;
    }
    const auto& pu = this->phaseUsage();
    const int np = pu.num_phases;

    std::optional<std::string> distributedMSWellLocalErrorMessage;
    std::optional<std::string> connectionWithoutAssociatedSegmentsLocalErrorMessage;
    std::optional<std::string> initializePerforationMismatchLocalErrorMessage;

    // in the init function, the well rates and perforation rates have been initialized or copied from prevState
    // what we do here, is to set the segment rates and perforation rates
    for (int w = 0; w < nw; ++w) {
        const auto& well_ecl = wells_ecl[w];
        if (this->is_permanently_inactive_well(well_ecl.name()))
            continue;

        auto& ws = this->well(w);

        if (well_ecl.isMultiSegment()) {
            const WellSegments& segment_set = well_ecl.getSegments();
            // assuming the order of the perforations in well_ecl is the same with Wells
            const WellConnections& completion_set = well_ecl.getConnections();
            // number of segment for this single well
            ws.segments = SegmentState<Scalar>{np, segment_set};
            const int well_nseg = segment_set.size();
            int n_activeperf = 0;
            int n_activeperf_local = 0;

            // we need to know for each segment, how many perforation it has and how many segments using it as outlet_segment
            // that is why I think we should use a well model to initialize the WellState here
            std::vector<std::vector<int>> segment_perforations(well_nseg);
            std::unordered_map<int,int> active_perf_index_local_to_global = {};
            std::unordered_map<int,int> active_to_local = {};
            for (std::size_t perf = 0; perf < completion_set.size(); ++perf) {
                const Connection& connection = completion_set.get(perf);
                if (connection.state() == Connection::State::OPEN) {
                    const int segment_index = segment_set.segmentNumberToIndex(connection.segment());
                    if (segment_index == -1) {
                        if (!connectionWithoutAssociatedSegmentsLocalErrorMessage) {
                            connectionWithoutAssociatedSegmentsLocalErrorMessage = ""; // First well with error: initialize the string
                        }

                        // Append to the existing error message for errors in further wells
                        *connectionWithoutAssociatedSegmentsLocalErrorMessage +=
                            fmt::format("COMPSEGS: Well {} has connection in cell {}, {}, {} "
                            "without associated segment.\n", well_ecl.name(),
                            connection.getI() + 1 , connection.getJ() + 1,
                            connection.getK() + 1 );
                    }

                    segment_perforations[segment_index].push_back(n_activeperf);
                    if (ws.parallel_info.get().globalPerfToLocalPerf(perf) > -1) {
                        active_perf_index_local_to_global.insert({n_activeperf_local, n_activeperf});
                        active_to_local.insert({n_activeperf, ws.parallel_info.get().globalPerfToLocalPerf(perf)});
                        n_activeperf_local++;
                    }
                    n_activeperf++;
                }
            }
            ws.parallel_info.get().setActivePerfToLocalPerfMap(active_to_local);

            // Check if the multi-segment well is distributed across several processes by comparing the local number
            // of active perforations (ws.perf_data.size()) with the total number of active perforations (n_activeperf).
            // If they differ on any rank, the well is considered distributed.

            if (static_cast<int>(ws.perf_data.size()) != n_activeperf) {
                if (!distributedMSWellLocalErrorMessage) {
                    distributedMSWellLocalErrorMessage = ""; // First well with error: initialize the string
                }

                // Append to the existing error message for errors in further wells
                *distributedMSWellLocalErrorMessage += fmt::format(
                    "Distributed multi-segment well {} detected on rank {}.\n"
                    "The number of active perforations on rank {} is {}, but the total number of active perforations is {}.\n\n",
                    well_ecl.name(),
                    ws.parallel_info.get().communication().rank(),
                    ws.parallel_info.get().communication().rank(),
                    ws.perf_data.size(),
                    n_activeperf
                );
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
                for (std::size_t perf = 0; perf < perf_data.size(); perf++)
                    perf_rates[perf*np + gaspos] *= 100;
            }

            const auto& perf_rates = perf_data.phase_rates;
            const auto& perf_press = perf_data.pressure;
            // The function calculateSegmentRates as well as the loop filling the segment_pressure work
            // with *global* containers. Now we create global vectors containing the phase_rates and
            // pressures of all processes.
            size_t number_of_global_perfs = 0;

            if (ws.parallel_info.get().communication().size() > 1) {
                number_of_global_perfs = ws.parallel_info.get().communication().sum(perf_data.size());
            } else {
                number_of_global_perfs = perf_data.size();
            }

            std::vector<Scalar> perforation_rates(number_of_global_perfs * np, 0.0);
            std::vector<Scalar> perforation_pressures(number_of_global_perfs, 0.0);

            assert(perf_data.size() == perf_press.size());
            assert(perf_data.size() * np == perf_rates.size());
            for (size_t perf = 0; perf < perf_data.size(); ++perf) {
                if (auto candidate = active_perf_index_local_to_global.find(perf); candidate != active_perf_index_local_to_global.end()) {
                    const int global_active_perf_index = candidate->second;
                    perforation_pressures[global_active_perf_index] = perf_press[perf];
                    for (int i = 0; i < np; i++) {
                        perforation_rates[global_active_perf_index * np + i] = perf_rates[perf * np + i];
                    }
                } else {
                    if (!initializePerforationMismatchLocalErrorMessage) {
                        initializePerforationMismatchLocalErrorMessage = ""; // First well with error: initialize the string
                    }

                    // Append to the existing error message for errors in further wells
                    *initializePerforationMismatchLocalErrorMessage += fmt::format("Error when initializing MS Well state of well {}, there is no active perforation index for the local index {}.\n", well_ecl.name(), perf);
                }
            }
            if (ws.parallel_info.get().communication().size() > 1) {
                ws.parallel_info.get().communication().sum(perforation_rates.data(), perforation_rates.size());
                ws.parallel_info.get().communication().sum(perforation_pressures.data(), perforation_pressures.size());
            }

            calculateSegmentRates(ws.parallel_info, segment_inlets, segment_perforations, perforation_rates, np, 0 /* top segment */, ws.segments.rates);

            // for the segment pressure, the segment pressure is the same with the first perforation belongs to the segment
            // if there is no perforation associated with this segment, it uses the pressure from the outlet segment
            // which requres the ordering is successful
            // Not sure what is the best way to handle the initialization, hopefully, the bad initialization can be
            // improved during the solveWellEq process
            {
                // top segment is always the first one, and its pressure is the well bhp
                auto& segment_pressure = ws.segments.pressure;
                segment_pressure[0] = ws.bhp;
                // The segment_indices contain the indices of the segments, that are only available on one process.
                std::vector<int> segment_indices;
                for (int seg = 1; seg < well_nseg; ++seg) {
                    if (!segment_perforations[seg].empty()) {
                        const int first_perf_global_index = segment_perforations[seg][0];
                        segment_pressure[seg] = perforation_pressures[first_perf_global_index];
                        segment_indices.push_back(seg);
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

    {
        // Instead of throwing inside the loop over all wells, collect the findings and throw at the end.
        // Throwing inside of the loop over all wells leads to a deadlock in some cases.

        std::string fullErrorMessage = {};

        if (!this->enableDistributedWells_ && distributedMSWellLocalErrorMessage) {
            *distributedMSWellLocalErrorMessage += "Either adjust the partitioning or set the flag --allow-distributed-wells to true.\n";
            fullErrorMessage += *distributedMSWellLocalErrorMessage;
        }

        if (connectionWithoutAssociatedSegmentsLocalErrorMessage) {
            fullErrorMessage += *connectionWithoutAssociatedSegmentsLocalErrorMessage;
        }

        if (initializePerforationMismatchLocalErrorMessage) {
            fullErrorMessage += *initializePerforationMismatchLocalErrorMessage;
        }

        if (!fullErrorMessage.empty()) {
            throw std::logic_error(fullErrorMessage);
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

template<class Scalar>
void WellState<Scalar>::
calculateSegmentRatesBeforeSum(const ParallelWellInfo<Scalar>& pw_info,
                               const std::vector<std::vector<int>>& segment_inlets,
                               const std::vector<std::vector<int>>& segment_perforations,
                               const std::vector<Scalar>& perforation_rates,
                               const int np, const int segment,
                               std::vector<Scalar>& segment_rates)
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
        auto local_perf = pw_info.activePerfToLocalPerf(perf);
        // If local_perf == -1, then the perforation is not on this process.
        // The perforation of the other processes are added in calculateSegmentRates.
        if (local_perf > -1) {
            for (int p = 0; p < np; ++p) {
                segment_rates[np * segment + p] += perforation_rates[np * local_perf + p];
            }
        }
    }
    for (const int& inlet_seg : segment_inlets[segment]) {
        calculateSegmentRatesBeforeSum(pw_info, segment_inlets, segment_perforations, perforation_rates, np, inlet_seg, segment_rates);
        for (int p = 0; p < np; ++p) {
            segment_rates[np * segment + p] += segment_rates[np * inlet_seg + p];
        }
    }
}
template<class Scalar>
void WellState<Scalar>::
calculateSegmentRates(const ParallelWellInfo<Scalar>& pw_info,
                      const std::vector<std::vector<int>>& segment_inlets,
                      const std::vector<std::vector<int>>& segment_perforations,
                      const std::vector<Scalar>& perforation_rates,
                      const int np, const int segment,
                      std::vector<Scalar>& segment_rates)
{
    calculateSegmentRatesBeforeSum(pw_info, segment_inlets, segment_perforations, perforation_rates, np, segment, segment_rates);
    pw_info.communication().sum(segment_rates.data(), segment_rates.size());
}

template<class Scalar>
void WellState<Scalar>::stopWell(int well_index)
{
    auto& ws = this->well(well_index);
    ws.stop();
}

template<class Scalar>
void WellState<Scalar>::openWell(int well_index)
{
    auto& ws = this->well(well_index);
    ws.open();
}

template<class Scalar>
void WellState<Scalar>::shutWell(int well_index)
{
    auto& ws = this->well(well_index);
    ws.shut();
}

template<class Scalar>
void WellState<Scalar>::updateStatus(int well_index, WellStatus status)
{
    auto& ws = this->well(well_index);
    ws.updateStatus(status);
}

template<class Scalar>
void WellState<Scalar>::communicateGroupRates(const Parallel::Communication& comm)
{
    // Compute the size of the data.
    std::size_t sz = std::accumulate(this->well_rates.begin(), this->well_rates.end(), std::size_t{0},
                     [](const std::size_t acc, const auto& rates)
                     { return acc + rates.second.second.size(); });

    // Make a vector and collect all data into it.
    std::vector<Scalar> data(sz);
    auto pos = data.begin();
    std::for_each(this->well_rates.begin(), this->well_rates.end(),
                  [&pos](const auto& input)
                  {
                      const auto& [owner, rates] = input.second;
                      if (owner) {
                          std::copy(rates.begin(), rates.end(), pos);
                      }
                      pos += rates.size();
                  });

    assert(pos == data.end());

    // Communicate it with a single sum() call.
    comm.sum(data.data(), data.size());

    pos = data.begin();
    std::for_each(this->well_rates.begin(), this->well_rates.end(),
                  [&pos](auto& input)
                  {
                      auto& rates = input.second.second;
                      std::copy(pos, pos + rates.size(), rates.begin());
                      pos += rates.size();
                  });
    assert(pos == data.end());
}

template<class Scalar>
void WellState<Scalar>::updateGlobalIsGrup(const Parallel::Communication& comm)
{
    this->global_well_info.value().clear();
    for (std::size_t well_index = 0; well_index < this->size(); well_index++) {
        const auto& ws = this->well(well_index);
        this->global_well_info.value().update_efficiency_scaling_factor(well_index, ws.efficiency_scaling_factor);
        if (ws.producer)
            this->global_well_info.value().update_producer(well_index, ws.status, ws.production_cmode);
        else
            this->global_well_info.value().update_injector(well_index, ws.status, ws.injection_cmode);
    }
    this->global_well_info.value().communicate(comm);
}

template<class Scalar>
void WellState<Scalar>::
updateEfficiencyScalingFactor(const std::string& wellName,
                              const Scalar value)
{
    const auto idx = this->index(wellName);
    this->global_well_info.value().update_efficiency_scaling_factor(*idx, value);
    this->well(*idx).efficiency_scaling_factor = value;
}

template<class Scalar>
data::Segment
WellState<Scalar>::reportSegmentResults(const int well_id,
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

template <class Scalar>
void WellState<Scalar>::
reportConnectionFactors(const std::size_t well_index,
                        std::vector<data::Connection>& connections) const
{
    const auto& perf_data = this->well(well_index).perf_data;
    const auto num_perf_well = perf_data.size();

    for (auto i = 0*num_perf_well; i < num_perf_well; ++i) {
        auto& connection = connections[i];

        connection.trans_factor = perf_data.connection_transmissibility_factor[i];
        connection.d_factor = perf_data.connection_d_factor[i];
        connection.compact_mult = perf_data.connection_compaction_tmult[i];
    }
}

template <class Scalar>
void WellState<Scalar>::
reportConnectionPressuresAndRates(const std::size_t well_index,
                                  const PhaseUsage& pu,
                                  std::vector<data::Connection>& connections) const
{
    using rt = data::Rates::opt;

    const int np = pu.num_phases;
    std::vector<rt> phs(np);
    std::vector<rt> pi(np);

    if (pu.phase_used[Water]) {
        phs.at(pu.phase_pos[Water]) = rt::wat;
        pi .at(pu.phase_pos[Water]) = rt::productivity_index_water;
    }

    if (pu.phase_used[Oil]) {
        phs.at(pu.phase_pos[Oil]) = rt::oil;
        pi .at(pu.phase_pos[Oil]) = rt::productivity_index_oil;
    }

    if (pu.phase_used[Gas]) {
        phs.at(pu.phase_pos[Gas]) = rt::gas;
        pi .at(pu.phase_pos[Gas]) = rt::productivity_index_gas;
    }

    const auto& ws = this->well(well_index);
    const auto& perf_data = ws.perf_data;
    const auto  num_perf_well = perf_data.size();

    for (auto i = 0*num_perf_well; i < num_perf_well; ++i) {
        auto& connection = connections[i];

        {
            const auto* rates = &perf_data.phase_rates[np * i];
            const auto* connPI = &perf_data.prod_index[np * i];

            for (int p = 0; p < np; ++p) {
                connection.rates.set(phs[p], rates [p]);
                connection.rates.set(pi [p], connPI[p]);
            }
        }

        connection.pressure = perf_data.pressure[i];
        connection.reservoir_rate = perf_data.rates[i];

        connection.rates.set(rt::dissolved_gas, perf_data.phase_mixing_rates[i][ws.dissolved_gas]);
        connection.rates.set(rt::vaporized_oil, perf_data.phase_mixing_rates[i][ws.vaporized_oil]);
    }

    if (pu.has_polymer) {
        for (auto i = 0*num_perf_well; i < num_perf_well; ++i) {
            connections[i].rates.set(rt::polymer, perf_data.polymer_rates[i]);
        }
    }

    if (pu.has_brine) {
        for (auto i = 0*num_perf_well; i < num_perf_well; ++i) {
            connections[i].rates.set(rt::brine, perf_data.brine_rates[i]);
        }
    }

    if (pu.has_solvent) {
        for (auto i = 0*num_perf_well; i < num_perf_well; ++i) {
            connections[i].rates.set(rt::solvent, perf_data.solvent_rates[i]);
        }
    }

    if (pu.has_micp) {
        for (auto i = 0*num_perf_well; i < num_perf_well; ++i) {
            connections[i].rates.set(rt::microbial, perf_data.microbial_rates[i]);
            connections[i].rates.set(rt::oxygen, perf_data.oxygen_rates[i]);
            connections[i].rates.set(rt::urea, perf_data.urea_rates[i]);
        }
    }
    
    if (pu.has_co2_or_h2store) {
        for (auto i = 0*num_perf_well; i < num_perf_well; ++i) {
            connections[i].rates.set(rt::mass_gas, perf_data.gas_mass_rates[i]);
        }
    }

    if (pu.phase_used[Water]) {
        for (auto i = 0*num_perf_well; i < num_perf_well; ++i) {
            connections[i].rates.set(rt::mass_wat, perf_data.wat_mass_rates[i]);
        }
    }
}

template <class Scalar>
void WellState<Scalar>::
reportConnectionFilterCake(const std::size_t well_index,
                           std::vector<data::Connection>& connections) const
{
    const auto& perf_data = this->well(well_index).perf_data;
    const auto num_perf_well = perf_data.size();

    const auto& filtrate_data = perf_data.filtrate_data;

    for (auto i = 0*num_perf_well; i < num_perf_well; ++i) {
        auto& filtrate = connections[i].filtrate;

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

template <class Scalar>
void WellState<Scalar>::
reportFractureStatistics(const std::vector<ConnFracStatistics<Scalar>>& stats,
                         std::vector<data::Connection>& connections) const
{
    using Quantity = typename ConnFracStatistics<Scalar>::Quantity;
    using StatResult = data::ConnectionFracturing;

    auto connIx = 0*connections.size();
    for (auto& connection : connections) {
        for (const auto& [q, result] : {
                std::pair { Quantity::Pressure, &StatResult::press },
                std::pair { Quantity::FlowRate, &StatResult::rate  },
                std::pair { Quantity::Width   , &StatResult::width },
            })
        {
            const auto& stat = stats[connIx].statistics(q);

            if (stat.sampleSize() > 0) {
                auto& x = connection.fract.*result;

                x.avg = stat.mean();
                x.min = stat.min();
                x.max = stat.max();

                if (const auto stdev = stat.stdev(); stdev.has_value()) {
                    x.stdev = *stdev;
                }

                connection.fract.numCells = stat.sampleSize();
            }
        }

        ++connIx;
    }
}

template<class Scalar>
bool WellState<Scalar>::wellIsOwned(std::size_t well_index,
                                    [[maybe_unused]] const std::string& wellName) const
{
    const auto& well_info = this->parallelWellInfo(well_index);
    assert(well_info.name() == wellName);

    return well_info.isOwner();
}

template<class Scalar>
bool WellState<Scalar>::wellIsOwned(const std::string& wellName) const
{
    const auto& well_index = this->index(wellName);
    if (!well_index.has_value()) {
        OPM_THROW(std::logic_error,
                  fmt::format("Could not find well {} in well map", wellName));
    }

    return wellIsOwned(well_index.value(), wellName);
}

template <typename Scalar>
void WellState<Scalar>::updateWellsDefaultALQ(const Schedule& schedule,
                                              const int report_step,
                                              const SummaryState& summary_state)
{
    const auto wells = schedule.wellNames(report_step);
    for (const auto& wname : wells) {
        const auto& well = schedule.getWell(wname, report_step);
        const auto& well_index = this->index(wname);
        if (! well.isProducer() || !well_index) {
            continue;
        }
        const auto alq = well.alq_value(summary_state);
        this->well(wname).alq_state.update_default(alq);
    }
}

template<class Scalar>
bool WellState<Scalar>::operator==(const WellState& rhs) const
{
    return this->well_rates == rhs.well_rates &&
           this->wells_ == rhs.wells_ &&
           this->permanently_inactive_well_names_ == rhs.permanently_inactive_well_names_;
}

template<class Scalar>
const ParallelWellInfo<Scalar>&
WellState<Scalar>::parallelWellInfo(std::size_t well_index) const
{
    const auto& ws = this->well(well_index);
    return ws.parallel_info;
}


template class WellState<double>;

#if FLOW_INSTANTIATE_FLOAT
template class WellState<float>;
#endif

} // namespace Opm
