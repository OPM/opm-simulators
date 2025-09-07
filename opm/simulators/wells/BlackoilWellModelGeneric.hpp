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

#ifndef OPM_BLACKOILWELLMODEL_GENERIC_HEADER_INCLUDED
#define OPM_BLACKOILWELLMODEL_GENERIC_HEADER_INCLUDED

#include <opm/output/data/GuideRateValue.hpp>

#include <opm/input/eclipse/Schedule/Group/GuideRate.hpp>
#include <opm/input/eclipse/Schedule/Well/PAvg.hpp>
#include <opm/input/eclipse/Schedule/Well/PAvgCalculator.hpp>
#include <opm/input/eclipse/Schedule/Well/PAvgCalculatorCollection.hpp>
#include <opm/input/eclipse/Schedule/Well/WellTestState.hpp>

#include <opm/simulators/utils/DeferredLoggingErrorHelpers.hpp>
#include <opm/simulators/wells/BlackoilWellModelWBP.hpp>
#include <opm/simulators/wells/ConnectionIndexMap.hpp>
#include <opm/simulators/wells/ParallelPAvgDynamicSourceData.hpp>
#include <opm/simulators/wells/ParallelWBPCalculation.hpp>
#include <opm/simulators/wells/ParallelWellInfo.hpp>
#include <opm/simulators/wells/PerforationData.hpp>
#include <opm/simulators/wells/WellFilterCake.hpp>
#include <opm/simulators/wells/WellProdIndexCalculator.hpp>
#include <opm/simulators/wells/WellTracerRate.hpp>
#include <opm/simulators/wells/WGState.hpp>

#include <algorithm>
#include <cstddef>
#include <functional>
#include <map>
#include <memory>
#include <optional>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>


namespace Opm {
    class DeferredLogger;
    class EclipseState;
    template<typename Scalar, typename IndexTraits> class BlackoilWellModelGasLiftGeneric;
    template<typename Scalar, typename IndexTraits> class GasLiftGroupInfo;
    template<typename Scalar, typename IndexTraits> class GasLiftSingleWellGeneric;
    template<class Scalar> class GasLiftWellState;
    class Group;
    class GuideRateConfig;
    class RestartValue;
    class Schedule;
    struct SimulatorUpdate;
    class SummaryConfig;
    template<typename Scalar, typename IndexTraits> class VFPProperties;
    template<typename Scalar, typename IndexTraits> class WellGroupHelpers;
    template<typename Scalar, typename IndexTraits> class WellInterfaceGeneric;
    template<typename Scalar, typename IndexTraits> class WellState;
} // namespace Opm

namespace Opm { namespace data {
    struct GroupData;
    struct GroupGuideRates;
    class GroupAndNetworkValues;
    struct NodeData;
}} // namespace Opm::data

namespace Opm::Parameters {

struct EnableTerminalOutput { static constexpr bool value = true; };

} // namespace Opm::Parameters

namespace Opm {

/// Class for handling the blackoil well model.
template<typename Scalar, typename IndexTraits>
class BlackoilWellModelGeneric
{
    using WellGroupHelpersType =  WellGroupHelpers<Scalar, IndexTraits>;
public:
    BlackoilWellModelGeneric(Schedule& schedule,
                             BlackoilWellModelGasLiftGeneric<Scalar, IndexTraits>& gaslift,
                             const SummaryState& summaryState,
                             const EclipseState& eclState,
                             const PhaseUsageInfo<IndexTraits>& phase_usage,
                             const Parallel::Communication& comm);

    virtual ~BlackoilWellModelGeneric() = default;
    virtual int compressedIndexForInteriorLGR([[maybe_unused]] const std::string& lgr_tag,
                                              [[maybe_unused]] const Connection& conn) const
    {
        throw std::runtime_error("compressedIndexForInteriorLGR not implemented");
    }

    int numLocalWells() const;
    int numLocalWellsEnd() const;
    int numLocalNonshutWells() const;
    int numPhases() const;

    /// return true if wells are available in the reservoir
    bool wellsActive() const;

    //! \brief Returns true if well is defined and has connections on current rank.
    bool hasLocalWell(const std::string& wname) const;

    //! \brief Returns true if well is defined, open and has connections on current rank.
    bool hasOpenLocalWell(const std::string& well_name) const;

    /// return true if network is active (at least one network well in prediction mode)
    bool networkActive() const;

    // whether there exists any multisegment well open on this process
    bool anyMSWellOpenLocal() const;

    const std::vector<Well>& eclWells() const
    { return wells_ecl_; }

    bool terminalOutput() const
    { return terminal_output_; }

    const Well& getWellEcl(const std::string& well_name) const;
    std::vector<Well> getLocalWells(const int timeStepIdx) const;
    const Schedule& schedule() const { return schedule_; }
    const PhaseUsageInfo<IndexTraits>& phaseUsage() const { return phase_usage_info_; }
    const GroupState<Scalar>& groupState() const { return this->active_wgstate_.group_state; }
    std::vector<const WellInterfaceGeneric<Scalar, IndexTraits>*> genericWells() const
    { return {well_container_generic_.begin(), well_container_generic_.end()}; }

    std::vector<WellInterfaceGeneric<Scalar, IndexTraits>*> genericWells()
    { return well_container_generic_; }

    /*
      Immutable version of the currently active wellstate.
    */
    const WellState<Scalar, IndexTraits>& wellState() const
    {
        return this->active_wgstate_.well_state;
    }

    /*
      Mutable version of the currently active wellstate.
    */
    WellState<Scalar, IndexTraits>& wellState()
    {
        return this->active_wgstate_.well_state;
    }

    /*
      Will return the currently active nupcolWellState; must update
      the internal nupcol wellstate with updateNupcolWGState() first.
    */
    const WellState<Scalar, IndexTraits>& nupcolWellState() const
    {
        return this->nupcol_wgstate_.well_state;
    }
    GroupState<Scalar>& groupState() { return this->active_wgstate_.group_state; }

    WellTestState& wellTestState() { return this->active_wgstate_.well_test_state; }

    const WellTestState& wellTestState() const { return this->active_wgstate_.well_test_state; }

    Scalar wellPI(const int well_index) const;
    Scalar wellPI(const std::string& well_name) const;

    void updateEclWells(const int timeStepIdx,
                        const SimulatorUpdate& sim_update,
                        const SummaryState& st);

    void initFromRestartFile(const RestartValue& restartValues,
                             std::unique_ptr<WellTestState> wtestState,
                             const std::size_t numCells,
                             bool handle_ms_well,
                             bool enable_distributed_wells);

    void prepareDeserialize(int report_step,
                            const std::size_t numCells,
                            bool handle_ms_well,
                            bool enable_distributed_wells);

    /*
      Will assign the internal member last_valid_well_state_ to the
      current value of the this->active_well_state_. The state stored
      with storeWellState() can then subsequently be recovered with the
      resetWellState() method.
    */
    void commitWGState()
    {
        this->last_valid_wgstate_ = this->active_wgstate_;
        this->last_valid_node_pressures_ = this->node_pressures_;
    }

    data::GroupAndNetworkValues groupAndNetworkData(const int reportStepIdx) const;

    /// Checks if network is active (at least one network well on prediction).
    void updateNetworkActiveState(const int report_step);

    /// Checks if there are reasons to perform a pre-step network re-balance.
    /// (Currently, the only reasons are network well status changes.)
    /// (TODO: Consider if adding network change events would be helpful.)
    bool needPreStepNetworkRebalance(const int report_step) const;

    /// Shut down any single well
    /// Returns true if the well was actually found and shut.
    bool forceShutWellByName(const std::string& wellname,
                             const double simulation_time,
                             const bool dont_shut_grup_wells);

    const std::vector<PerforationData<Scalar>>& perfData(const int well_idx) const
    { return well_perf_data_[well_idx]; }

    const Parallel::Communication& comm() const { return comm_; }

    const EclipseState& eclipseState() const { return eclState_; }

    const SummaryState& summaryState() const { return summaryState_; }

    const GuideRate& guideRate() const { return guideRate_; }
    GuideRate& guideRate() { return guideRate_; }

    const std::map<std::string, double>& wellOpenTimes() const { return well_open_times_; }
    const std::map<std::string, double>& wellCloseTimes() const { return well_close_times_; }
    const WellGroupEvents& reportStepStartEvents() const { return report_step_start_events_; }

    std::vector<int> getCellsForConnections(const Well& well) const;

    bool reportStepStarts() const { return report_step_starts_; }

    bool shouldBalanceNetwork(const int reportStepIndex,
                              const int iterationIdx) const;

    void updateClosedWellsThisStep(const std::string& well_name) const
    {
        this->closed_this_step_.insert(well_name);
    }
    bool wasDynamicallyShutThisTimeStep(const std::string& well_name) const;

    void logPrimaryVars() const;

    template<class Serializer>
    void serializeOp(Serializer& serializer)
    {
        serializer(initial_step_);
        serializer(report_step_starts_);
        serializer(last_run_wellpi_);
        serializer(local_shut_wells_);
        serializer(closed_this_step_);
        serializer(guideRate_);
        serializer(node_pressures_);
        serializer(last_valid_node_pressures_);
        serializer(prev_inj_multipliers_);
        serializer(active_wgstate_);
        serializer(last_valid_wgstate_);
        serializer(nupcol_wgstate_);
        serializer(switched_prod_groups_);
        serializer(switched_inj_groups_);
        serializer(closed_offending_wells_);
        serializer(gen_gaslift_);
    }

    bool operator==(const BlackoilWellModelGeneric& rhs) const;

    const ParallelWellInfo<Scalar>&
    parallelWellInfo(const std::size_t idx) const
    { return local_parallel_well_info_[idx].get(); }

    bool isOwner(const std::string& wname) const
    {
        auto pwInfoPos = std::find_if(this->parallel_well_info_.begin(),
                                      this->parallel_well_info_.end(),
                                      [&wname](const auto& pwInfo)
                                      { return pwInfo.name() == wname; });

        return (pwInfoPos != this->parallel_well_info_.end())
            && pwInfoPos->isOwner();
    }

    const ConnectionIndexMap& connectionIndexMap(const std::size_t idx)
    { return conn_idx_map_[idx]; }

protected:
    /*
      The dynamic state of the well model is maintained with an instance
      of the WellState class. Currently we have
      three different wellstate instances:

       1. The currently active wellstate is in the active_well_state_
          member. That is the state which is mutated by the simulator.

       2. In the case timestep fails to converge and we must go back and
          try again with a smaller timestep we need to recover the last
          valid wellstate. This is maintained with the
          last_valid_well_state_ member and the functions
          commitWGState() and resetWellState().

        3. For the NUPCOL functionality we should either use the
           currently active wellstate or a wellstate frozen at max
           nupcol iterations. This is handled with the member
           nupcol_well_state_ and the updateNupcolWGState() function.
    */

    /*
      Will return the last good wellstate. This is typcially used when
      initializing a new report step where the Schedule object might
      have introduced new wells. The wellstate returned by
      prevWellState() must have been stored with the commitWGState()
      function first.
    */
    const WellState<Scalar, IndexTraits>& prevWellState() const
    {
        return this->last_valid_wgstate_.well_state;
    }

    const WGState<Scalar, IndexTraits>& prevWGState() const
    {
        return this->last_valid_wgstate_;
    }

    /*
      Will store a copy of the input argument well_state in the
      last_valid_well_state_ member, that state can then be recovered
      with a subsequent call to resetWellState().
    */
    void commitWGState(WGState<Scalar, IndexTraits> wgstate)
    {
        this->last_valid_wgstate_ = std::move(wgstate);
    }

    /*
      Will update the internal variable active_well_state_ to whatever
      was stored in the last_valid_well_state_ member. This function
      works in pair with commitWellState() which should be called first.
    */
    void resetWGState()
    {
        this->active_wgstate_ = this->last_valid_wgstate_;
        this->node_pressures_ = this->last_valid_node_pressures_;
    }

    /*
      Will store the current active wellstate in the nupcol_well_state_
      member. This can then be subsequently retrieved with accessor
      nupcolWellState().
    */
    void updateNupcolWGState()
    {
        this->nupcol_wgstate_ = this->active_wgstate_;
    }

    void reportGroupSwitching(DeferredLogger& local_deferredLogger) const;

    /// \brief Create the parallel well information
    /// \param localWells The local wells from ECL schedule
    std::vector<std::reference_wrapper<ParallelWellInfo<Scalar>>>
    createLocalParallelWellInfo(const std::vector<Well>& wells);

    void initializeWellProdIndCalculators();
    void initializeWellPerfData();

    bool wasDynamicallyShutThisTimeStep(const int well_index) const;

    Scalar updateNetworkPressures(const int reportStepIdx,
                                  const Scalar damping_factor,
                                  const Scalar update_upper_bound);

    void updateWsolvent(const Group& group,
                        const int reportStepIdx,
                        const WellState<Scalar, IndexTraits>& wellState);
    void setWsolvent(const Group& group,
                     const int reportStepIdx,
                     Scalar wsolvent);
    virtual void calcResvCoeff(const int fipnum,
                               const int pvtreg,
                               const std::vector<Scalar>& production_rates,
                               std::vector<Scalar>& resv_coeff) = 0;
    virtual void calcInjResvCoeff(const int fipnum,
                                  const int pvtreg,
                                  std::vector<Scalar>& resv_coeff) = 0;

    /// Assign dynamic well status for each well owned by current rank
    ///
    /// \param[in,out] wsrpt Well solution object.  On exit, holds current
    /// values for \code data::Well::dynamicStatus \endcode.
    void assignDynamicWellStatus(data::Wells& wsrpt) const;

    /// Assign basic result quantities for shut connections of wells owned
    /// by current rank.
    ///
    /// Mostly provided for summary file output purposes.  Applies to fully
    /// shut/stopped wells and shut connections of open/flowing wells.
    ///
    /// \param[in,out] wsrpt Well solution object.  On exit, also contains a
    /// few quantities, like the D factor, the Kh product and the CTF, for
    /// shut connections.
    ///
    /// \param[in] reportStepIdx Zero-based index of current report step.
    void assignShutConnections(data::Wells& wsrpt,
                               const int reportStepIndex) const;

    void assignWellTargets(data::Wells& wsrpt) const;
    void assignProductionWellTargets(const Well& well, data::WellControlLimits& limits) const;
    void assignInjectionWellTargets(const Well& well, data::WellControlLimits& limits) const;

    void assignGroupControl(const Group& group,
                            data::GroupData& gdata) const;
    void assignGroupValues(const int reportStepIdx,
                           std::map<std::string, data::GroupData>& gvalues) const;
    void assignNodeValues(std::map<std::string, data::NodeData>& nodevalues,
                          const int reportStepIdx) const;

    void calculateEfficiencyFactors(const int reportStepIdx);

    void checkGconsaleLimits(const Group& group,
                             WellState<Scalar, IndexTraits>& well_state,
                             const int reportStepIdx,
                             DeferredLogger& deferred_logger);

    void checkGEconLimits(const Group& group,
                          const double simulation_time,
                          const int report_step_idx,
                          DeferredLogger& deferred_logger);

    bool checkGroupHigherConstraints(const Group& group,
                                     DeferredLogger& deferred_logger,
                                     const int reportStepIdx,
                                     const int max_number_of_group_switch,
                                     const bool update_group_switching_log);

    void updateAndCommunicateGroupData(const int reportStepIdx,
                                       const int iterationIdx,
                                       const Scalar tol_nupcol,
                                       const bool update_wellgrouptarget, // we only want to update the wellgrouptarget after the groups have found their controls
                                       DeferredLogger& deferred_logger);

    void inferLocalShutWells();

    void setRepRadiusPerfLength();

    virtual void computePotentials(const std::size_t widx,
                                   const WellState<Scalar, IndexTraits>& well_state_copy,
                                   std::string& exc_msg,
                                   ExceptionType::ExcEnum& exc_type,
                                   DeferredLogger& deferred_logger) = 0;

    // Calculating well potentials for each well
    void updateWellPotentials(const int reportStepIdx,
                              const bool onlyAfterEvent,
                              const SummaryConfig& summaryConfig,
                              DeferredLogger& deferred_logger);

    void initInjMult();

    void updateInjMult(DeferredLogger& deferred_logger);
    void updateInjFCMult(DeferredLogger& deferred_logger);

    void updateFiltrationModelsPostStep(const double dt,
                                        const std::size_t water_index,
                                        DeferredLogger& deferred_logger);

    void updateFiltrationModelsPreStep(DeferredLogger& deferred_logger);

    // create the well container
    virtual void createWellContainer(const int time_step) = 0;
    virtual void initWellContainer(const int reportStepIdx) = 0;

    virtual void calculateProductivityIndexValuesShutWells(const int reportStepIdx,
                                                           DeferredLogger& deferred_logger) = 0;
    virtual void calculateProductivityIndexValues(DeferredLogger& deferred_logger) = 0;

    void runWellPIScaling(const int reportStepIdx,
                          DeferredLogger& local_deferredLogger);

    /// \brief get compressed index for interior cells (-1, otherwise
    virtual int compressedIndexForInterior(int cartesian_cell_idx) const = 0;

    std::vector<std::vector<int>> getMaxWellConnections() const;

    std::vector<std::string> getWellsForTesting(const int timeStepIdx,
                                                const double simulationTime);

    using WellTracerRates = std::unordered_map<int, std::vector<WellTracerRate<Scalar>>>;
    void assignWellTracerRates(data::Wells& wsrpt,
                               const WellTracerRates& wellTracerRates,
                               const unsigned reportStep) const;

    using MswTracerRates = std::unordered_map<int, std::vector<MSWellTracerRate<Scalar>>>;
    void assignMswTracerRates(data::Wells& wsrpt,
                              const MswTracerRates& mswTracerRates,
                              const unsigned reportStep) const;

    void assignMassGasRate(data::Wells& wsrpt,
                           const Scalar gasDensity) const;

    Schedule& schedule_;
    
    const SummaryState& summaryState_;
    const EclipseState& eclState_;
    const Parallel::Communication& comm_;
    BlackoilWellModelGasLiftGeneric<Scalar, IndexTraits>& gen_gaslift_;
    BlackoilWellModelWBP<Scalar, IndexTraits> wbp_;


    const PhaseUsageInfo<IndexTraits>& phase_usage_info_;
    bool terminal_output_{false};
    bool wells_active_{false};
    bool network_active_{false};
    bool initial_step_{};
    bool report_step_starts_{};

    std::optional<int> last_run_wellpi_{};

    std::vector<Well> wells_ecl_;
    std::vector<std::vector<PerforationData<Scalar>>> well_perf_data_;

    // Times at which wells were opened (for WCYCLE)
    std::map<std::string, double> well_open_times_;

    // Times at which wells were shut (for WCYCLE)
    std::map<std::string, double> well_close_times_;

    std::vector<ConnectionIndexMap> conn_idx_map_{};
    std::function<bool(const std::string&)> not_on_process_{};

    // a vector of all the wells.
    std::vector<WellInterfaceGeneric<Scalar, IndexTraits>*> well_container_generic_{};

    std::vector<int> local_shut_wells_{};

    std::vector<ParallelWellInfo<Scalar>> parallel_well_info_;
    std::vector<std::reference_wrapper<ParallelWellInfo<Scalar>>> local_parallel_well_info_;

    std::vector<WellProdIndexCalculator<Scalar>> prod_index_calc_;

    std::vector<int> pvt_region_idx_;

    mutable std::unordered_set<std::string> closed_this_step_;

    GuideRate guideRate_;
    std::unique_ptr<VFPProperties<Scalar, IndexTraits>> vfp_properties_{};

    // Network pressures for output and initialization
    std::map<std::string, Scalar> node_pressures_;
    // Valid network pressures for output and initialization for safe restart after failed iterations
    std::map<std::string, Scalar> last_valid_node_pressures_;

    // previous injection multiplier, it is used in the injection multiplier calculation for WINJMULT keyword
    std::unordered_map<std::string, std::vector<Scalar>> prev_inj_multipliers_;

    // Handling for filter cake injection multipliers
    std::unordered_map<std::string, WellFilterCake<Scalar, IndexTraits>> filter_cake_;

    /*
      The various wellState members should be accessed and modified
      through the accessor functions wellState(), prevWellState(),
      commitWellState(), resetWellState(), nupcolWellState() and
      updateNupcolWGState().
    */
    WGState<Scalar, IndexTraits> active_wgstate_;
    WGState<Scalar, IndexTraits> last_valid_wgstate_;
    WGState<Scalar, IndexTraits> nupcol_wgstate_;
    WellGroupEvents report_step_start_events_; //!< Well group events at start of report step

    bool wellStructureChangedDynamically_{false};

    // Store maps of group name and new group controls for output
    std::map<std::string, std::vector<Group::ProductionCMode>> switched_prod_groups_;
    std::map<std::string, std::array<std::vector<Group::InjectionCMode>, 3>> switched_inj_groups_;
    // Store map of group name and close offending well for output
    std::map<std::string, std::pair<std::string, std::string>> closed_offending_wells_;

private:
    WellInterfaceGeneric<Scalar, IndexTraits>* getGenWell(const std::string& well_name);

    template <typename Iter, typename Body>
    void wellUpdateLoop(Iter first, Iter last, const int timeStepIdx, Body&& body);

    void updateEclWellsConstraints(const int timeStepIdx,
                                   const SimulatorUpdate& sim_update,
                                   const SummaryState& st);

    void updateEclWellsCTFFromAction(const int timeStepIdx,
                                     const SimulatorUpdate& sim_update);

    /// Run caller-defined code for each well owned by current rank
    ///
    /// 'const' version.
    ///
    /// \tparam LoopBody Call-back type for user-defined code.  Expected to
    /// support a function call operator of the form
    /// \code
    ///   void operator()(const std::size_t i, const Well& well) const
    /// \endcode
    /// which will be invoked for each well owned by the local process.  The
    /// parameters are the index into \c wells_ecl_ and the object at that
    /// index, respectively.
    ///
    /// \param[in] loopBody Call-back function.  Typically defined as a
    /// lambda expression.
    template <typename LoopBody>
    void loopOwnedWells(LoopBody&& loopBody) const;
};

} // namespace Opm

#endif
