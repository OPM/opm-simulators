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

#include <opm/simulators/wells/WellState.hpp>
#include <opm/simulators/wells/ALQState.hpp>
#include <opm/simulators/wells/GlobalWellInfo.hpp>
#include <opm/simulators/wells/WellContainer.hpp>
#include <opm/core/props/BlackoilPhases.hpp>

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
class WellStateFullyImplicitBlackoil
    : public WellState
{
    typedef WellState  BaseType;
public:
    static const uint64_t event_mask = ScheduleEvents::WELL_STATUS_CHANGE + ScheduleEvents::PRODUCTION_UPDATE + ScheduleEvents::INJECTION_UPDATE;
    typedef BaseType :: WellMapType WellMapType;

    virtual ~WellStateFullyImplicitBlackoil() = default;

    // TODO: same definition with WellInterface, eventually they should go to a common header file.
    static const int Water = BlackoilPhases::Aqua;
    static const int Oil = BlackoilPhases::Liquid;
    static const int Gas = BlackoilPhases::Vapour;

    using BaseType :: wellRates;
    using BaseType :: bhp;
    using BaseType :: perfPress;
    using BaseType :: wellMap;
    using BaseType :: numWells;
    using BaseType :: numPhases;
    using BaseType :: resetConnectionTransFactors;
    using BaseType :: updateStatus;

    explicit WellStateFullyImplicitBlackoil(const PhaseUsage& pu) :
        WellState(pu)
    {
    }

    /// Allocate and initialize if wells is non-null.  Also tries
    /// to give useful initial values to the bhp(), wellRates()
    /// and perfPhaseRates() fields, depending on controls
    void init(const std::vector<double>& cellPressures,
              const Schedule& schedule,
              const std::vector<Well>& wells_ecl,
              const std::vector<ParallelWellInfo*>& parallel_well_info,
              const int report_step,
              const WellStateFullyImplicitBlackoil* prevState,
              const std::vector<std::vector<PerforationData>>& well_perf_data,
              const SummaryState& summary_state);

    void resize(const std::vector<Well>& wells_ecl,
                const std::vector<ParallelWellInfo*>& parallel_well_info,
                const Schedule& schedule,
                const bool handle_ms_well,
                const size_t numCells,
                const std::vector<std::vector<PerforationData>>& well_perf_data,
                const SummaryState& summary_state);

    /// One rate per phase and well connection.
    std::vector<double>& mutable_perfPhaseRates() { return perfphaserates_; }
    const std::vector<double>& perfPhaseRates() const { return perfphaserates_; }

    /// One current control per injecting well.
    Well::InjectorCMode currentInjectionControl(std::size_t well_index) const { return current_injection_controls_[well_index]; }
    void currentInjectionControl(std::size_t well_index, Well::InjectorCMode cmode) { current_injection_controls_[well_index] = cmode; }

    /// One current control per producing well.
    Well::ProducerCMode currentProductionControl(std::size_t well_index) const { return current_production_controls_[well_index]; }
    void currentProductionControl(std::size_t well_index, Well::ProducerCMode cmode) { current_production_controls_[well_index] = cmode; }

    void setCurrentWellRates(const std::string& wellName, const std::vector<double>& rates ) {
        well_rates[wellName].second = rates;
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

    void reportConnections(data::Well& well, const PhaseUsage &pu,
                           const WellMapType::value_type& wt,
                           const int* globalCellIdxMap) const;

    /// init the MS well related.
    void initWellStateMSWell(const std::vector<Well>& wells_ecl,
                             const WellStateFullyImplicitBlackoil* prev_well_state);

    static void calculateSegmentRates(const std::vector<std::vector<int>>& segment_inlets, const std::vector<std::vector<int>>&segment_perforations,
                                      const std::vector<double>& perforation_rates, const int np, const int segment, std::vector<double>& segment_rates);

    Events& events(std::size_t well_index) {
        return this->events_[well_index];
    }

    const std::vector<int>& firstPerfIndex() const
    {
        return first_perf_index_;
    }

    /// One rate pr well connection.
    std::vector<double>& perfRateSolvent() { return perfRateSolvent_; }
    const std::vector<double>& perfRateSolvent() const { return perfRateSolvent_; }

    /// One rate pr well
    double solventWellRate(const int w) const;

    /// One rate pr well connection.
    std::vector<double>& perfRatePolymer() { return perfRatePolymer_; }
    const std::vector<double>& perfRatePolymer() const { return perfRatePolymer_; }

    /// One rate pr well
    double polymerWellRate(const int w) const;

    /// One rate pr well connection.
    std::vector<double>& perfRateBrine() { return perfRateBrine_; }
    const std::vector<double>& perfRateBrine() const { return perfRateBrine_; }

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

    std::vector<double>& wellDissolvedGasRates()
    {
        return well_dissolved_gas_rates_;
    }

    std::vector<double>& wellVaporizedOilRates()
    {
        return well_vaporized_oil_rates_;
    }

    const std::vector<double>& segRates() const
    {
        return seg_rates_;
    }

    std::vector<double>& segRates()
    {
        return seg_rates_;
    }

    const std::vector<double>& segPress() const
    {
        return seg_press_;
    }

    std::vector<double>& segPressDrop()
    {
        return seg_pressdrop_;
    }

    const std::vector<double>& segPressDrop() const
    {
        return seg_pressdrop_;
    }

    std::vector<double>& segPressDropFriction()
    {
        return seg_pressdrop_friction_;
    }

    const std::vector<double>& segPressDropFriction() const
    {
        return seg_pressdrop_friction_;
    }

    std::vector<double>& segPressDropHydroStatic()
    {
        return seg_pressdrop_hydorstatic_;
    }

    const std::vector<double>& segPressDropHydroStatic() const
    {
        return seg_pressdrop_hydorstatic_;
    }

    std::vector<double>& segPressDropAcceleration()
    {
        return seg_pressdrop_acceleration_;
    }

    const std::vector<double>& segPressDropAcceleration() const
    {
        return seg_pressdrop_acceleration_;
    }

    std::vector<double>& segPress()
    {
        return seg_press_;
    }

    int numSegment() const
    {
        return nseg_;
    }

    int topSegmentIndex(const int w) const;

    std::vector<double>& productivityIndex() {
        return productivity_index_;
    }

    const std::vector<double>& productivityIndex() const {
        return productivity_index_;
    }

    std::vector<double>& connectionProductivityIndex() {
        return this->conn_productivity_index_;
    }

    const std::vector<double>& connectionProductivityIndex() const {
        return this->conn_productivity_index_;
    }

    std::vector<double>& wellPotentials() {
        return well_potentials_;
    }

    const std::vector<double>& wellPotentials() const {
        return well_potentials_;
    }

    std::vector<double>& perfThroughput() {
        return perf_water_throughput_;
    }

    const std::vector<double>& perfThroughput() const {
        return perf_water_throughput_;
    }

    std::vector<double>& perfSkinPressure() {
        return perf_skin_pressure_;
    }

    const std::vector<double>& perfSkinPressure() const {
        return perf_skin_pressure_;
    }

    std::vector<double>& perfWaterVelocity() {
        return perf_water_velocity_;
    }

    const std::vector<double>& perfWaterVelocity() const {
        return perf_water_velocity_;
    }

    void shutWell(int well_index) override;

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
        disableGliftOptimization();
    }

    void disableGliftOptimization() {
        do_glift_optimization_ = false;
    }

    void enableGliftOptimization() {
        do_glift_optimization_ = true;
    }

    int wellNameToGlobalIdx(const std::string &name) {
        return this->global_well_info.value().well_index(name);
    }

    std::string globalIdxToWellName(const int index) {
        return this->global_well_info.value().well_name(index);
    }

private:
    std::vector<double> perfphaserates_;
    WellContainer<int> is_producer_; // Size equal to number of local wells.

    // vector with size number of wells +1.
    // iterate over all perforations of a given well
    // for (int perf = first_perf_index_[well_index]; perf < first_perf_index_[well_index] + num_perf_[well_index]; ++perf)
    std::vector<int> first_perf_index_;
    std::vector<int> num_perf_;
    WellContainer<Opm::Well::InjectorCMode> current_injection_controls_;
    WellContainer<Well::ProducerCMode> current_production_controls_;

    // Use of std::optional<> here is a technical crutch, the
    // WellStateFullyImplicitBlackoil class should be default constructible,
    // whereas the GlobalWellInfo is not.
    std::optional<GlobalWellInfo> global_well_info;
    std::map<std::string, std::pair<bool, std::vector<double>>> well_rates;

    ALQState alq_state;
    bool do_glift_optimization_;

    std::vector<double> perfRateSolvent_;

    // only for output
    std::vector<double> perfRatePolymer_;
    std::vector<double> perfRateBrine_;

    // it is the throughput of water flow through the perforations
    // it is used as a measure of formation damage around well-bore due to particle deposition
    // it will only be used for injectors to check the injectivity
    std::vector<double> perf_water_throughput_;

    // skin pressure of peforation
    // it will only be used for injectors to check the injectivity
    std::vector<double> perf_skin_pressure_;

    // it will only be used for injectors to check the injectivity
    // water velocity of perforation
    std::vector<double> perf_water_velocity_;

    // phase rates under reservoir condition for wells
    // or voidage phase rates
    WellContainer<std::vector<double>> well_reservoir_rates_;

    // dissolved gas rates or solution gas production rates
    // should be zero for injection wells
    std::vector<double> well_dissolved_gas_rates_;

    // vaporized oil rates or solution oil producation rates
    // should be zero for injection wells
    std::vector<double> well_vaporized_oil_rates_;

    // some events happens to the well, like this well is a new well
    // or new well control keywords happens
    // \Note: for now, only WCON* keywords, and well status change is considered
    WellContainer<Events> events_;

    // MS well related
    // for StandardWell, the number of segments will be one
    std::vector<double> seg_rates_;
    std::vector<double> seg_press_;
    // the index of the top segments, which is used to locate the
    // multisegment well related information in WellState
    std::vector<int> top_segment_index_;
    int nseg_; // total number of the segments

    // The following data are only recorded for output
    // pressure drop
    std::vector<double> seg_pressdrop_;
    // frictional pressure drop
    std::vector<double> seg_pressdrop_friction_;
    // hydrostatic pressure drop
    std::vector<double> seg_pressdrop_hydorstatic_;
    // accelerational pressure drop
    std::vector<double> seg_pressdrop_acceleration_;

    // Productivity Index
    std::vector<double> productivity_index_;

    // Connection-level Productivity Index
    std::vector<double> conn_productivity_index_;

    // Well potentials
    std::vector<double> well_potentials_;

    /// Map segment index to segment number, mostly for MS wells.
    ///
    /// Segment number (one-based) of j-th segment in i-th well is
    /// \code
    ///    const auto top    = topSegmentIndex(i);
    ///    const auto seg_No = seg_number_[top + j];
    /// \end
    std::vector<int> seg_number_;

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
};

} // namespace Opm


#endif // OPM_WELLSTATEFULLYIMPLICITBLACKOIL_HEADER_INCLUDED
