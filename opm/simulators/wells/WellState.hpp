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

#ifndef OPM_WELLSTATEFULLYIMPLICITBLACKOIL_HEADER_INCLUDED
#define OPM_WELLSTATEFULLYIMPLICITBLACKOIL_HEADER_INCLUDED

#include <opm/common/ErrorMacros.hpp>

#include <opm/core/props/BlackoilPhases.hpp>

#include <opm/simulators/wells/ALQState.hpp>
#include <opm/simulators/wells/GlobalWellInfo.hpp>
#include <opm/simulators/wells/PerfData.hpp>
#include <opm/simulators/wells/PerforationData.hpp>
#include <opm/simulators/wells/SegmentState.hpp>
#include <opm/simulators/wells/SingleWellState.hpp>
#include <opm/simulators/wells/WellContainer.hpp>

#include <opm/output/data/Wells.hpp>

#include <opm/input/eclipse/Schedule/Events.hpp>

#include <dune/common/version.hh>
#include <dune/common/parallel/mpihelper.hh>

#include <functional>
#include <map>
#include <optional>
#include <string>
#include <utility>
#include <vector>

namespace Opm
{

class ParallelWellInfo;
class Schedule;
enum class WellStatus;

/// The state of a set of wells, tailored for use by the fully
/// implicit blackoil simulator.
class WellState
{
public:
    static const uint64_t event_mask = ScheduleEvents::WELL_STATUS_CHANGE + ScheduleEvents::PRODUCTION_UPDATE + ScheduleEvents::INJECTION_UPDATE;
    // TODO: same definition with WellInterface, eventually they should go to a common header file.
    static const int Water = BlackoilPhases::Aqua;
    static const int Oil = BlackoilPhases::Liquid;
    static const int Gas = BlackoilPhases::Vapour;

    // Only usable for testing purposes
    explicit WellState(const ParallelWellInfo& pinfo);

    explicit WellState(const PhaseUsage& pu)
        : phase_usage_(pu)
    {}

    static WellState serializationTestObject(const ParallelWellInfo& pinfo);

    std::size_t size() const {
        return this->wells_.size();
    }

    std::vector<std::string> wells() const {
        return this->wells_.wells();
    }


    int numWells() const
    {
        return this->size();
    }

    const ParallelWellInfo& parallelWellInfo(std::size_t well_index) const;



    /// Allocate and initialize if wells is non-null.  Also tries
    /// to give useful initial values to the bhp(), wellRates()
    /// and perfPhaseRatesORG() fields, depending on controls
    void init(const std::vector<double>& cellPressures,
              const Schedule& schedule,
              const std::vector<Well>& wells_ecl,
              const std::vector<std::reference_wrapper<ParallelWellInfo>>& parallel_well_info,
              const int report_step,
              const WellState* prevState,
              const std::vector<std::vector<PerforationData>>& well_perf_data,
              const SummaryState& summary_state);

    void resize(const std::vector<Well>& wells_ecl,
                const std::vector<std::reference_wrapper<ParallelWellInfo>>& parallel_well_info,
                const Schedule& schedule,
                const bool handle_ms_well,
                const std::size_t numCells,
                const std::vector<std::vector<PerforationData>>& well_perf_data,
                const SummaryState& summary_state);

    void setCurrentWellRates(const std::string& wellName, const std::vector<double>& new_rates ) {
        auto& [owner, rates] = this->well_rates.at(wellName);
        if (owner)
            rates = new_rates;
    }

    const std::vector<double>& currentWellRates(const std::string& wellName) const;

    bool hasWellRates(const std::string& wellName) const {
        return this->well_rates.find(wellName) != this->well_rates.end();
    }

    void clearWellRates()
    {
        this->well_rates.clear();
    }

    template<class Communication>
    void gatherVectorsOnRoot(const std::vector< data::Connection >& from_connections,
                             std::vector< data::Connection >& to_connections,
                             const Communication& comm) const;

    data::Wells
    report(const int* globalCellIdxMap,
           const std::function<bool(const int)>& wasDynamicallyClosed) const;

    void reportConnections(std::vector<data::Connection>& connections, const PhaseUsage &pu,
                           std::size_t well_index,
                           const int* globalCellIdxMap) const;

    /// init the MS well related.
    void initWellStateMSWell(const std::vector<Well>& wells_ecl,
                             const WellState* prev_well_state);

    static void calculateSegmentRates(const std::vector<std::vector<int>>& segment_inlets, const std::vector<std::vector<int>>&segment_perforations,
                                      const std::vector<double>& perforation_rates, const int np, const int segment, std::vector<double>& segment_rates);


    template<class Comm>
    void communicateGroupRates(const Comm& comm);

    template<class Comm>
    void updateGlobalIsGrup(const Comm& comm);

    bool isInjectionGrup(const std::string& name) const {
        return this->global_well_info.value().in_injecting_group(name);
    }

    bool isProductionGrup(const std::string& name) const {
        return this->global_well_info.value().in_producing_group(name);
    }

    double getALQ( const std::string& name) const
    {
        return this->alq_state.get(name);
    }

    void setALQ( const std::string& name, double value)
    {
        this->alq_state.set(name, value);
    }

    int gliftGetDebugCounter() {
        return this->alq_state.get_debug_counter();
    }

    void gliftSetDebugCounter(int value) {
        return this->alq_state.set_debug_counter(value);
    }

    int gliftUpdateDebugCounter() {
        return this->alq_state.update_debug_counter();
    }

    bool gliftCheckAlqOscillation(const std::string &name) const {
        return this->alq_state.oscillation(name);
    }

    int gliftGetAlqDecreaseCount(const std::string &name) {
        return this->alq_state.get_decrement_count(name);
    }

    int gliftGetAlqIncreaseCount(const std::string &name) {
        return this->alq_state.get_increment_count(name);
    }

    void gliftUpdateAlqIncreaseCount(const std::string &name, bool increase) {
        this->alq_state.update_count(name, increase);
    }

    void gliftTimeStepInit() {
        this->alq_state.reset_count();
    }

    int wellNameToGlobalIdx(const std::string &name) {
        return this->global_well_info.value().well_index(name);
    }

    std::string globalIdxToWellName(const int index) {
        return this->global_well_info.value().well_name(index);
    }

    bool wellIsOwned(std::size_t well_index,
                     const std::string& wellName) const;

    bool wellIsOwned(const std::string& wellName) const;

    void updateStatus(int well_index, WellStatus status);

    void openWell(int well_index);
    void shutWell(int well_index);
    void stopWell(int well_index);

    /// The number of phases present.
    int numPhases() const
    {
        return this->phase_usage_.num_phases;
    }

    const PhaseUsage& phaseUsage() const {
        return this->phase_usage_;
    }

    /// One rate per well and phase.
    std::vector<double>& wellRates(std::size_t well_index) { return this->wells_[well_index].surface_rates; }
    const std::vector<double>& wellRates(std::size_t well_index) const { return this->wells_[well_index].surface_rates; }

    const std::string& name(std::size_t well_index) const {
        return this->wells_.well_name(well_index);
    }

    std::optional<std::size_t> index(const std::string& well_name) const {
        return this->wells_.well_index(well_name);
    }

    const SingleWellState& operator[](std::size_t well_index) const {
        return this->wells_[well_index];
    }

    const SingleWellState& operator[](const std::string& well_name) const {
        return this->wells_[well_name];
    }

    SingleWellState& operator[](std::size_t well_index) {
        return this->wells_[well_index];
    }

    SingleWellState& operator[](const std::string& well_name) {
        return this->wells_[well_name];
    }

    const SingleWellState& well(std::size_t well_index) const {
        return this->operator[](well_index);
    }

    const SingleWellState& well(const std::string& well_name) const {
        return this->operator[](well_name);
    }

    SingleWellState& well(std::size_t well_index) {
        return this->operator[](well_index);
    }

    SingleWellState& well(const std::string& well_name) {
        return this->operator[](well_name);
    }

    bool has(const std::string& well_name) const {
        return this->wells_.has(well_name);
    }

    bool operator==(const WellState&) const;

    template<class Serializer>
    void serializeOp(Serializer& serializer)
    {
        serializer(alq_state);
        serializer(well_rates);
        if (serializer.isSerializing()) {
            serializer(wells_.size());
        } else {
            std::size_t size = 0;
            serializer(size);
            if (size != wells_.size()) {
                OPM_THROW(std::runtime_error, "Error deserializing WellState: size mismatch");
            }
        }
        for (auto& w : wells_) {
            serializer(w);
        }
    }

private:
    PhaseUsage phase_usage_;

    // The wells_ variable is essentially a map of all the wells on the current
    // process. Observe that since a well can be split over several processes a
    // well might appear in the WellContainer on different processes.
    WellContainer<SingleWellState> wells_;

    // The members alq_state, global_well_info and well_rates are map like
    // structures which will have entries for *all* the wells in the system.

    // Use of std::optional<> here is a technical crutch, the
    // WellStateFullyImplicitBlackoil class should be default constructible,
    // whereas the GlobalWellInfo is not.
    std::optional<GlobalWellInfo> global_well_info;
    ALQState alq_state;

    // The well_rates variable is defined for all wells on all processors. The
    // bool in the value pair is whether the current process owns the well or
    // not.
    std::map<std::string, std::pair<bool, std::vector<double>>> well_rates;

    data::Segment
    reportSegmentResults(const int         well_id,
                         const int         seg_ix,
                         const int         seg_no) const;

    // If the ALQ has changed since the previous report step,
    // reset current_alq and update default_alq. ALQ is used for
    // constant lift gas injection and for gas lift optimization
    // (THP controlled wells).

    void updateWellsDefaultALQ(const std::vector<Well>& wells_ecl, const SummaryState& summary_state);


    /// Allocate and initialize if wells is non-null.
    /// Also tries to give useful initial values to the bhp() and
    /// wellRates() fields, depending on controls.  The
    /// perfRates() field is filled with zero, and perfPress()
    /// with -1e100.
    void base_init(const std::vector<double>& cellPressures,
                   const std::vector<Well>& wells_ecl,
                   const std::vector<std::reference_wrapper<ParallelWellInfo>>& parallel_well_info,
                   const std::vector<std::vector<PerforationData>>& well_perf_data,
                   const SummaryState& summary_state);

    void initSingleWell(const std::vector<double>& cellPressures,
                        const Well& well,
                        const std::vector<PerforationData>& well_perf_data,
                        const ParallelWellInfo& well_info,
                        const SummaryState& summary_state);

    void initSingleProducer(const Well& well,
                            const ParallelWellInfo& well_info,
                            double pressure_first_connection,
                            const std::vector<PerforationData>& well_perf_data,
                            const SummaryState& summary_state);

    void initSingleInjector(const Well& well,
                            const ParallelWellInfo& well_info,
                            double pressure_first_connection,
                            const std::vector<PerforationData>& well_perf_data,
                            const SummaryState& summary_state);

};

} // namespace Opm


#endif // OPM_WELLSTATEFULLYIMPLICITBLACKOIL_HEADER_INCLUDED
