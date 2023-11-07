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

#include <opm/simulators/wells/ParallelPAvgDynamicSourceData.hpp>
#include <opm/simulators/wells/ParallelWBPCalculation.hpp>
#include <opm/simulators/wells/PerforationData.hpp>
#include <opm/simulators/wells/WellFilterCake.hpp>
#include <opm/simulators/wells/WellProdIndexCalculator.hpp>
#include <opm/simulators/wells/WGState.hpp>

#include <cstddef>
#include <functional>
#include <map>
#include <memory>
#include <optional>
#include <string>
#include <type_traits>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace Opm {
    class DeferredLogger;
    class EclipseState;
    class GasLiftGroupInfo;
    class GasLiftSingleWellGeneric;
    class GasLiftWellState;
    class Group;
    class GuideRateConfig;
    class ParallelWellInfo;
    class RestartValue;
    class Schedule;
    struct SimulatorUpdate;
    class SummaryConfig;
    class VFPProperties;
    class WellInterfaceGeneric;
    class WellState;
} // namespace Opm

namespace Opm { namespace data {
    struct GroupData;
    struct GroupGuideRates;
    class GroupAndNetworkValues;
    struct NodeData;
}} // namespace Opm::data

namespace Opm {

/// Class for handling the blackoil well model.
class BlackoilWellModelGeneric
{
public:
    // ---------      Types      ---------
    using GLiftOptWells = std::map<std::string, std::unique_ptr<GasLiftSingleWellGeneric>>;
    using GLiftProdWells = std::map<std::string, const WellInterfaceGeneric*>;
    using GLiftWellStateMap = std::map<std::string, std::unique_ptr<GasLiftWellState>>;

    BlackoilWellModelGeneric(Schedule& schedule,
                             const SummaryState& summaryState,
                             const EclipseState& eclState,
                             const PhaseUsage& phase_usage,
                             const Parallel::Communication& comm);

    virtual ~BlackoilWellModelGeneric() = default;

    int numLocalWells() const;
    int numLocalWellsEnd() const;
    int numLocalNonshutWells() const;
    int numPhases() const;

    /// return true if wells are available in the reservoir
    bool wellsActive() const;
    bool hasWell(const std::string& wname) const;

    /// return true if network is active (at least one network well in prediction mode)
    bool networkActive() const;

    // whether there exists any multisegment well open on this process
    bool anyMSWellOpenLocal() const;

    const Well& getWellEcl(const std::string& well_name) const;
    std::vector<Well> getLocalWells(const int timeStepIdx) const;
    const Schedule& schedule() const { return schedule_; }
    const PhaseUsage& phaseUsage() const { return phase_usage_; }
    const GroupState& groupState() const { return this->active_wgstate_.group_state; }
    std::vector<const WellInterfaceGeneric*> genericWells() const
    { return {well_container_generic_.begin(), well_container_generic_.end()}; }

    /*
      Immutable version of the currently active wellstate.
    */
    const WellState& wellState() const
    {
        return this->active_wgstate_.well_state;
    }

    /*
      Mutable version of the currently active wellstate.
    */
    WellState& wellState()
    {
        return this->active_wgstate_.well_state;
    }

    /*
      Will return the currently active nupcolWellState; must initialize
      the internal nupcol wellstate with initNupcolWellState() first.
    */
    const WellState& nupcolWellState() const
    {
        return this->nupcol_wgstate_.well_state;
    }
    GroupState& groupState() { return this->active_wgstate_.group_state; }

    WellTestState& wellTestState() { return this->active_wgstate_.well_test_state; }

    const WellTestState& wellTestState() const { return this->active_wgstate_.well_test_state; }


    double wellPI(const int well_index) const;
    double wellPI(const std::string& well_name) const;

    void updateEclWells(const int timeStepIdx,
                        const SimulatorUpdate& sim_update,
                        const SummaryState& st);

    void initFromRestartFile(const RestartValue& restartValues,
                             WellTestState wtestState,
                             const std::size_t numCells,
                             bool handle_ms_well);

    void prepareDeserialize(int report_step,
                            const std::size_t numCells,
                            bool handle_ms_well);

    /*
      Will assign the internal member last_valid_well_state_ to the
      current value of the this->active_well_state_. The state stored
      with storeWellState() can then subsequently be recovered with the
      resetWellState() method.
    */
    void commitWGState()
    {
        this->last_valid_wgstate_ = this->active_wgstate_;
    }

    data::GroupAndNetworkValues groupAndNetworkData(const int reportStepIdx) const;

    /// Return true if any well has a THP constraint.
    bool hasTHPConstraints() const;

    /// Checks if network is active (at least one network well on prediction).
    void updateNetworkActiveState(const int report_step);

    /// Checks if there are reasons to perform a pre-step network re-balance.
    /// (Currently, the only reasons are network well status changes.)
    /// (TODO: Consider if adding network change events would be helpful.)
    bool needPreStepNetworkRebalance(const int report_step) const;

    /// Shut down any single well
    /// Returns true if the well was actually found and shut.
    bool forceShutWellByName(const std::string& wellname,
                             const double simulation_time);

    const std::vector<PerforationData>& perfData(const int well_idx) const
    { return well_perf_data_[well_idx]; }

    const Parallel::Communication& comm() const { return comm_; }

    const EclipseState& eclipseState() const { return eclState_; }

    const SummaryState& summaryState() const { return summaryState_; }

    const GuideRate& guideRate() const { return guideRate_; }

    bool reportStepStarts() const { return report_step_starts_; }

    bool shouldBalanceNetwork(const int reportStepIndex,
                              const int iterationIdx) const;

    void updateClosedWellsThisStep(const std::string& well_name) const {
        this->closed_this_step_.insert(well_name);
    }

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
        serializer(prev_inj_multipliers_);
        serializer(active_wgstate_);
        serializer(last_valid_wgstate_);
        serializer(nupcol_wgstate_);
        serializer(last_glift_opt_time_);
        serializer(switched_prod_groups_);
        serializer(switched_inj_groups_);
    }

    bool operator==(const BlackoilWellModelGeneric& rhs) const
    {
        return this->initial_step_ == rhs.initial_step_ &&
               this->report_step_starts_ == rhs.report_step_starts_ &&
               this->last_run_wellpi_ == rhs.last_run_wellpi_ &&
               this->local_shut_wells_ == rhs.local_shut_wells_ &&
               this->closed_this_step_ == rhs.closed_this_step_ &&
               this->node_pressures_ == rhs.node_pressures_ &&
               this->prev_inj_multipliers_ == rhs.prev_inj_multipliers_ &&
               this->active_wgstate_ == rhs.active_wgstate_ &&
               this->last_valid_wgstate_ == rhs.last_valid_wgstate_ &&
               this->nupcol_wgstate_ == rhs.nupcol_wgstate_ &&
               this->last_glift_opt_time_ == rhs.last_glift_opt_time_ &&
               this->switched_prod_groups_ == rhs.switched_prod_groups_ &&
               this->switched_inj_groups_ == rhs.switched_inj_groups_;
    }

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
          commitWellState() and resetWellState().

        3. For the NUPCOL functionality we should either use the
           currently active wellstate or a wellstate frozen at max
           nupcol iterations. This is handled with the member
           nupcol_well_state_ and the initNupcolWellState() function.
    */


    /*
      Will return the last good wellstate. This is typcially used when
      initializing a new report step where the Schedule object might
      have introduced new wells. The wellstate returned by
      prevWellState() must have been stored with the commitWellState()
      function first.
    */
    const WellState& prevWellState() const
    {
        return this->last_valid_wgstate_.well_state;
    }


    const WGState& prevWGState() const
    {
        return this->last_valid_wgstate_;
    }



    /*
      Will store a copy of the input argument well_state in the
      last_valid_well_state_ member, that state can then be recovered
      with a subsequent call to resetWellState().
    */
    void commitWGState(WGState wgstate)
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

    /// \brief Create the parallel well information
    /// \param localWells The local wells from ECL schedule
    std::vector<std::reference_wrapper<ParallelWellInfo>> createLocalParallelWellInfo(const std::vector<Well>& wells);

    void initializeWellProdIndCalculators();
    void initializeWellPerfData();

    bool wasDynamicallyShutThisTimeStep(const int well_index) const;

    double updateNetworkPressures(const int reportStepIdx);

    void updateWsolvent(const Group& group,
                        const int reportStepIdx,
                        const WellState& wellState);
    void setWsolvent(const Group& group,
                     const int reportStepIdx,
                     double wsolvent);
    virtual void calcRates(const int fipnum,
                           const int pvtreg,
                           const std::vector<double>& production_rates,
                           std::vector<double>& resv_coeff) = 0;
    virtual void calcInjRates(const int fipnum,
                           const int pvtreg,
                           std::vector<double>& resv_coeff) = 0;

    void assignShutConnections(data::Wells& wsrpt,
                               const int reportStepIndex) const;
    void assignGroupControl(const Group& group,
                            data::GroupData& gdata) const;
    void assignGroupValues(const int reportStepIdx,
                           std::map<std::string, data::GroupData>& gvalues) const;
    void assignNodeValues(std::map<std::string, data::NodeData>& nodevalues) const;

    void calculateEfficiencyFactors(const int reportStepIdx);

    void checkGconsaleLimits(const Group& group,
                             WellState& well_state,
                             const int reportStepIdx,
                             DeferredLogger& deferred_logger);

    void checkGEconLimits(const Group& group,
                          const double simulation_time,
                          const int report_step_idx,
                          DeferredLogger& deferred_logger);

    bool checkGroupHigherConstraints(const Group& group,
                                     DeferredLogger& deferred_logger,
                                     const int reportStepIdx);

    void updateAndCommunicateGroupData(const int reportStepIdx,
                                       const int iterationIdx);

    void inferLocalShutWells();

    void setRepRadiusPerfLength();

    void gliftDebug(const std::string& msg,
                    DeferredLogger& deferred_logger) const;

    void gliftDebugShowALQ(DeferredLogger& deferred_logger);

    void gasLiftOptimizationStage2(DeferredLogger& deferred_logger,
                                   GLiftProdWells& prod_wells,
                                   GLiftOptWells& glift_wells,
                                   GasLiftGroupInfo& group_info,
                                   GLiftWellStateMap& map,
                                   const int episodeIndex);

    virtual void computePotentials(const std::size_t widx,
                                   const WellState& well_state_copy,
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

    void updateFiltrationParticleVolume(const double dt, const std::size_t water_index);

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

    std::vector<int> getCellsForConnections(const Well& well) const;
    std::vector<std::vector<int>> getMaxWellConnections() const;

    std::vector<std::string> getWellsForTesting(const int timeStepIdx,
                                                const double simulationTime);

    using WellTracerRates = std::map<std::pair<std::string, std::string>, double>;
    void assignWellTracerRates(data::Wells& wsrpt,
                               const WellTracerRates& wellTracerRates) const;

    Schedule& schedule_;
    const SummaryState& summaryState_;
    const EclipseState& eclState_;
    const Parallel::Communication& comm_;

    PhaseUsage phase_usage_;
    bool terminal_output_{false};
    bool wells_active_{false};
    bool network_active_{false};
    bool initial_step_{};
    bool report_step_starts_{};

    std::optional<int> last_run_wellpi_{};

    std::vector<Well> wells_ecl_;
    std::vector<std::vector<PerforationData>> well_perf_data_;

    /// Connection index mappings
    class ConnectionIndexMap
    {
    public:
        /// Constructor.
        ///
        /// \param[in] numConns Total number of well connections, both open
        ///   and closed/shut.  Typically \code WellConnections::size() \endcode.
        explicit ConnectionIndexMap(const std::size_t numConns)
            : local_(numConns, -1)
        {
            this->global_.reserve(numConns);
            this->open_.reserve(numConns);
        }

        /// Enumerate/map new active connection.
        ///
        /// \param[in] connIdx Global well connection index.  Must be an
        ///   integer in the range 0..numConns-1.
        ///
        /// \param[in] connIsOpen Whether or not the connection is
        ///   open/flowing.
        void addActiveConnection(const int  connIdx,
                                 const bool connIsOpen)
        {
            this->local_[connIdx] =
                static_cast<int>(this->global_.size());

            this->global_.push_back(connIdx);

            const auto open_conn_idx = connIsOpen
                ? this->num_open_conns_++
                : -1;

            this->open_.push_back(open_conn_idx);
        }

        /// Get local connection IDs/indices of every existing well
        /// connection.
        ///
        /// Negative value (-1) for connections that don't intersect the
        /// current rank.
        const std::vector<int>& local() const
        {
            return this->local_;
        }

        /// Get global connection ID of local (on-rank) connection.
        ///
        /// \param[in] connIdx Local connection index.
        ///
        /// \return Global connection ID of \p connIdx.
        int global(const int connIdx) const
        {
            return this->global_[connIdx];
        }

        /// Get open connection ID of local (on-rank) connection.
        ///
        /// \param[in] connIdx Local connection index.
        ///
        /// \return Open connection ID of \p connIdx.  Integer in the range
        ///   0..#open connections - 1 if the connection is open or negative
        ///   value (-1) otherwise.
        int open(const int connIdx) const
        {
            return this->open_[connIdx];
        }

    private:
        /// Local connection IDs/indices of every existing well connection.
        /// Negative value (-1) for connections that don't intersect the
        /// current rank.
        std::vector<int> local_{};

        /// Global connection index of each on-rank reservoir connection.
        /// Reverse/transpose mapping of \c local_.
        std::vector<int> global_{};

        /// Open connection index of each on-rank reservoir connection.
        std::vector<int> open_{};

        /// Number of open connections on this rank.
        int num_open_conns_{0};
    };

    std::vector<ConnectionIndexMap> conn_idx_map_{};
    std::function<bool(const Well&)> not_on_process_{};

    // a vector of all the wells.
    std::vector<WellInterfaceGeneric*> well_container_generic_{};

    std::vector<int> local_shut_wells_{};

    std::vector<ParallelWellInfo> parallel_well_info_;
    std::vector<std::reference_wrapper<ParallelWellInfo>> local_parallel_well_info_;

    std::vector<WellProdIndexCalculator> prod_index_calc_;
    mutable ParallelWBPCalculation wbpCalculationService_;

    std::vector<int> pvt_region_idx_;

    mutable std::unordered_set<std::string> closed_this_step_;

    GuideRate guideRate_;
    std::unique_ptr<VFPProperties> vfp_properties_{};
    std::map<std::string, double> node_pressures_; // Storing network pressures for output.

    // previous injection multiplier, it is used in the injection multiplier calculation for WINJMULT keyword
    std::unordered_map<std::string, std::vector<double>> prev_inj_multipliers_;

    // Handling for filter cake injection multipliers
    std::unordered_map<std::string, WellFilterCake> filter_cake_;

    /*
      The various wellState members should be accessed and modified
      through the accessor functions wellState(), prevWellState(),
      commitWellState(), resetWellState(), nupcolWellState() and
      updateNupcolWellState().
    */
    WGState active_wgstate_;
    WGState last_valid_wgstate_;
    WGState nupcol_wgstate_;

    bool glift_debug = false;

    double last_glift_opt_time_ = -1.0;

    bool wellStructureChangedDynamically_{false};

    std::map<std::string, std::string> switched_prod_groups_;
    std::map<std::pair<std::string, Opm::Phase>, std::string> switched_inj_groups_;

private:
    WellInterfaceGeneric* getGenWell(const std::string& well_name);
};


} // namespace Opm

#endif
