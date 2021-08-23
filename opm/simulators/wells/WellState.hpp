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

#include <opm/simulators/wells/ALQState.hpp>
#include <opm/simulators/wells/SingleWellState.hpp>
#include <opm/simulators/wells/GlobalWellInfo.hpp>
#include <opm/simulators/wells/SegmentState.hpp>
#include <opm/simulators/wells/WellContainer.hpp>
#include <opm/core/props/BlackoilPhases.hpp>
#include <opm/simulators/wells/PerforationData.hpp>
#include <opm/simulators/wells/PerfData.hpp>
#include <opm/output/data/Wells.hpp>

#include <opm/parser/eclipse/EclipseState/Schedule/Events.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/Well/Well.hpp>

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

/// The state of a set of wells, tailored for use by the fully
/// implicit blackoil simulator.
class WellState
{
public:
    using mapentry_t = std::array<int, 3>;
    using WellMapType = std::map<std::string, mapentry_t>;

    static const uint64_t event_mask = ScheduleEvents::WELL_STATUS_CHANGE + ScheduleEvents::PRODUCTION_UPDATE + ScheduleEvents::INJECTION_UPDATE;

    virtual ~WellState() = default;

    // TODO: same definition with WellInterface, eventually they should go to a common header file.
    static const int Water = BlackoilPhases::Aqua;
    static const int Oil = BlackoilPhases::Liquid;
    static const int Gas = BlackoilPhases::Vapour;

    explicit WellState(const PhaseUsage& pu)
    {
        this->phase_usage_ = pu;
    }

    const WellMapType& wellMap() const { return wellMap_; }
    WellMapType& wellMap() { return wellMap_; }

    std::size_t size() const {
        return this->wellMap_.size();
    }


    int numWells() const
    {
        return this->size();
    }

    int wellIndex(const std::string& wellName) const;

    const ParallelWellInfo& parallelWellInfo(std::size_t well_index) const;



    /// Allocate and initialize if wells is non-null.  Also tries
    /// to give useful initial values to the bhp(), wellRates()
    /// and perfPhaseRatesORG() fields, depending on controls
    void init(const std::vector<double>& cellPressures,
              const Schedule& schedule,
              const std::vector<Well>& wells_ecl,
              const std::vector<ParallelWellInfo*>& parallel_well_info,
              const int report_step,
              const WellState* prevState,
              const std::vector<std::vector<PerforationData>>& well_perf_data,
              const SummaryState& summary_state);

    void resize(const std::vector<Well>& wells_ecl,
                const std::vector<ParallelWellInfo*>& parallel_well_info,
                const Schedule& schedule,
                const bool handle_ms_well,
                const size_t numCells,
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

    /// One rate pr well
    double solventWellRate(const int w) const;

    /// One rate pr well
    double polymerWellRate(const int w) const;

    /// One rate pr well
    double brineWellRate(const int w) const;

    const WellContainer<std::vector<double>>& wellReservoirRates() const { return well_reservoir_rates_; }

    std::vector<double>& wellReservoirRates(std::size_t well_index)
    {
        return well_reservoir_rates_[well_index];
    }

    const std::vector<double>& wellReservoirRates(std::size_t well_index) const
    {
        return well_reservoir_rates_[well_index];
    }


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

    bool gliftOptimizationEnabled() const {
        return do_glift_optimization_;
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

    /// Special purpose method to support dynamically rescaling a well's
    /// CTFs through WELPI.
    ///
    /// \param[in] well_index Process-local linear index of single well.
    ///    Must be in the range 0..numWells()-1.
    ///
    /// \param[in] well_perf_data New perforation data.  Only
    ///    PerforationData::connection_transmissibility_factor actually
    ///    used (overwrites existing internal values).
    void resetConnectionTransFactors(const int well_index,
                                     const std::vector<PerforationData>& well_perf_data);

    void updateStatus(int well_index, Well::Status status);

    void openWell(int well_index) {
        this->wells_[well_index].status = Well::Status::OPEN;
    }

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
    const WellContainer<std::vector<double>>& wellRates() const { return wellrates_; }
    std::vector<double>& wellRates(std::size_t well_index) { return wellrates_[well_index]; }
    const std::vector<double>& wellRates(std::size_t well_index) const { return wellrates_[well_index]; }

    std::size_t numPerf(std::size_t well_index) const { return this->perfdata[well_index].size(); }

    PerfData& perfData(const std::string& wname) {
        return this->perfdata[wname];
    }

    const PerfData& perfData(const std::string& wname) const {
        return this->perfdata[wname];
    }

    PerfData& perfData(std::size_t well_index) {
        return this->perfdata[well_index];
    }

    const PerfData& perfData(std::size_t well_index) const {
        return this->perfdata[well_index];
    }

    const std::string& name(std::size_t well_index) const {
        return this->wells_.well_name(well_index);
    }

    const SingleWellState& well(std::size_t well_index) const {
        return this->wells_[well_index];
    }

    const SingleWellState& well(const std::string& well_name) const {
        return this->wells_[well_name];
    }

    SingleWellState& well(std::size_t well_index) {
        return this->wells_[well_index];
    }

    SingleWellState& well(const std::string& well_name) {
        return this->wells_[well_name];
    }

    bool has(const std::string& well_name) const {
        return this->wells_.has(well_name);
    }

private:
    WellMapType wellMap_;
    // Use of std::optional<> here is a technical crutch, the
    // WellStateFullyImplicitBlackoil class should be default constructible,
    // whereas the GlobalWellInfo is not.
    std::optional<GlobalWellInfo> global_well_info;
    ALQState alq_state;
    bool do_glift_optimization_;
    PhaseUsage phase_usage_;

    WellContainer<SingleWellState> wells_;
    WellContainer<const ParallelWellInfo*> parallel_well_info_;
    WellContainer<std::vector<double>> wellrates_;
    WellContainer<PerfData> perfdata;

    // The well_rates variable is defined for all wells on all processors. The
    // bool in the value pair is whether the current process owns the well or
    // not.
    std::map<std::string, std::pair<bool, std::vector<double>>> well_rates;

    // phase rates under reservoir condition for wells
    // or voidage phase rates
    WellContainer<std::vector<double>> well_reservoir_rates_;

    data::Segment
    reportSegmentResults(const PhaseUsage& pu,
                         const int         well_id,
                         const int         seg_ix,
                         const int         seg_no) const;

    int numSegments(const int well_id) const;

    int segmentNumber(const int well_id, const int seg_id) const;

    // If the ALQ has changed since the previous report step,
    // reset current_alq and update default_alq. ALQ is used for
    // constant lift gas injection and for gas lift optimization
    // (THP controlled wells).
    //
    // NOTE: If a well is no longer used (e.g. it is shut down)
    // it is still kept in the maps "default_alq_" and "current_alq_". Since the
    // number of unused entries should be small (negligible memory
    // overhead) this is simpler than writing code to delete it.
    //
    void updateWellsDefaultALQ(const std::vector<Well>& wells_ecl);

    /// Allocate and initialize if wells is non-null.
    /// Also tries to give useful initial values to the bhp() and
    /// wellRates() fields, depending on controls.  The
    /// perfRates() field is filled with zero, and perfPress()
    /// with -1e100.
    void base_init(const std::vector<double>& cellPressures,
                   const std::vector<Well>& wells_ecl,
                   const std::vector<ParallelWellInfo*>& parallel_well_info,
                   const std::vector<std::vector<PerforationData>>& well_perf_data,
                   const SummaryState& summary_state);

    void initSingleWell(const std::vector<double>& cellPressures,
                        const int w,
                        const Well& well,
                        const std::vector<PerforationData>& well_perf_data,
                        const ParallelWellInfo* well_info,
                        const SummaryState& summary_state);



};

} // namespace Opm


#endif // OPM_WELLSTATEFULLYIMPLICITBLACKOIL_HEADER_INCLUDED
