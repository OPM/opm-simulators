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

#include <functional>
#include <map>
#include <memory>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include <opm/output/data/GuideRateValue.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/Well/WellTestState.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/Group/GuideRate.hpp>

#include <opm/simulators/utils/DeferredLoggingErrorHelpers.hpp>
#include <opm/simulators/wells/ParallelWellInfo.hpp>
#include <opm/simulators/wells/PerforationData.hpp>
#include <opm/simulators/wells/WellInterfaceGeneric.hpp>
#include <opm/simulators/wells/WellProdIndexCalculator.hpp>
#include <opm/simulators/wells/WGState.hpp>

namespace Opm {

namespace data {
struct GroupData;
struct GroupGuideRates;
class GroupAndNetworkValues;
struct NodeData;
}

class DeferredLogger;
class EclipseState;
class GasLiftSingleWellGeneric;
class GasLiftWellState;
class Group;
class RestartValue;
class Schedule;
class SummaryConfig;
class VFPProperties;
class WellState;

/// Class for handling the blackoil well model.
class BlackoilWellModelGeneric
{
public:
    // ---------      Types      ---------
    using Comm = Dune::CollectiveCommunication<Dune::MPIHelper::MPICommunicator>;
    using GLiftOptWells = std::map<std::string,std::unique_ptr<GasLiftSingleWellGeneric>>;
    using GLiftProdWells = std::map<std::string,const WellInterfaceGeneric*>;
    using GLiftWellStateMap = std::map<std::string,std::unique_ptr<GasLiftWellState>>;

    BlackoilWellModelGeneric(Schedule& schedule,
                             const SummaryState& summaryState,
                             const EclipseState& eclState,
                             const PhaseUsage& phase_usage,
                             const Comm& comm);

    virtual ~BlackoilWellModelGeneric() = default;

    int numLocalWells() const;
    int numPhases() const;

    /// return true if wells are available in the reservoir
    bool wellsActive() const;
    bool hasWell(const std::string& wname);
    /// return true if wells are available on this process
    bool localWellsActive() const;
    // whether there exists any multisegment well open on this process
    bool anyMSWellOpenLocal() const;

    const Well& getWellEcl(const std::string& well_name) const;
    std::vector<Well> getLocalWells(const int timeStepIdx) const;
    const Schedule& schedule() const { return schedule_; }
    const PhaseUsage& phaseUsage() const { return phase_usage_; }
    const GroupState& groupState() const { return this->active_wgstate_.group_state; }

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

    GroupState& groupState() { return this->active_wgstate_.group_state; }


    double wellPI(const int well_index) const;
    double wellPI(const std::string& well_name) const;

    void updateEclWells(const int timeStepIdx,
                        const std::unordered_set<std::string>& wells);


    void loadRestartData(std::size_t report_step,
                         const data::Wells& rst_wells,
                         const data::GroupAndNetworkValues& grpNwrkValues,
                         const PhaseUsage& phases,
                         const bool handle_ms_well,
                         WellState& well_state);

    void initFromRestartFile(const RestartValue& restartValues,
                             const size_t numCells,
                             bool handle_ms_well);

    void setWellsActive(const bool wells_active);

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

    /// Shut down any single well, but only if it is in prediction mode.
    /// Returns true if the well was actually found and shut.
    bool forceShutWellByNameIfPredictionMode(const std::string& wellname,
                                             const double simulation_time);

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
      Will return the currently active nupcolWellState; must initialize
      the internal nupcol wellstate with initNupcolWellState() first.
    */
    const WellState& nupcolWellState() const
    {
        return this->nupcol_wgstate_.well_state;
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
    std::vector<ParallelWellInfo*> createLocalParallelWellInfo(const std::vector<Well>& wells);

    void initializeWellProdIndCalculators();
    void initializeWellPerfData();

    bool wasDynamicallyShutThisTimeStep(const int well_index) const;

    void updateNetworkPressures(const int reportStepIdx);

    void updateWsolvent(const Group& group,
                        const int reportStepIdx,
                        const WellState& wellState);
    void setWsolvent(const Group& group,
                     const int reportStepIdx,
                     double wsolvent);
    virtual void calcRates(const int fipnum,
                           const int pvtreg,
                           std::vector<double>& resv_coeff) = 0;
    virtual void calcInjRates(const int fipnum,
                           const int pvtreg,
                           std::vector<double>& resv_coeff) = 0;

    data::GuideRateValue getGuideRateValues(const Group& group) const;
    data::GuideRateValue getGuideRateValues(const Well& well) const;
    data::GuideRateValue getGuideRateInjectionGroupValues(const Group& group) const;
    void getGuideRateValues(const GuideRate::RateVector& qs,
                            const bool                   is_inj,
                            const std::string&           wgname,
                            data::GuideRateValue&        grval) const;

    void assignWellGuideRates(data::Wells& wsrpt) const;
    void assignShutConnections(data::Wells& wsrpt,
                               const int reportStepIndex) const;
    void assignGroupControl(const Group& group,
                            data::GroupData& gdata) const;
    void assignGroupGuideRates(const Group& group,
                               const std::unordered_map<std::string, data::GroupGuideRates>& groupGuideRates,
                               data::GroupData& gdata) const;
    void assignGroupValues(const int reportStepIdx,
                           std::map<std::string, data::GroupData>& gvalues) const;
    void assignNodeValues(std::map<std::string, data::NodeData>& nodevalues) const;

    std::unordered_map<std::string, data::GroupGuideRates>
    calculateAllGroupGuiderates(const int reportStepIdx) const;

    void calculateEfficiencyFactors(const int reportStepIdx);

    bool checkGroupConstraints(const Group& group,
                               const int reportStepIdx,
                               DeferredLogger& deferred_logger) const;

    std::pair<Group::InjectionCMode, double> checkGroupInjectionConstraints(const Group& group,
                                                         const int reportStepIdx,
                                                         const Phase& phase) const;
    std::pair<Group::ProductionCMode, double> checkGroupProductionConstraints(const Group& group,
                                                           const int reportStepIdx,
                                                           DeferredLogger& deferred_logger) const;

    void checkGconsaleLimits(const Group& group,
                             WellState& well_state,
                             const int reportStepIdx,
                             DeferredLogger& deferred_logger);

    void checkGroupHigherConstraints(const Group& group,
                                     DeferredLogger& deferred_logger,
                                     const int reportStepIdx,
                                     std::set<std::string>& switched_groups);

    void updateGroupIndividualControl(const Group& group,
                                      DeferredLogger& deferred_logger,
                                      const int reportStepIdx,
                                      std::set<std::string>& switched_groups);

    void updateGroupIndividualControls(DeferredLogger& deferred_logger,
                                       std::set<std::string>& switched_groups,
                                       const int reportStepIdx,
                                       const int iterationIdx);

    void updateGroupHigherControls(DeferredLogger& deferred_logger,
                                   const int reportStepIdx,
                                   std::set<std::string>& switched_groups);

    void actionOnBrokenConstraints(const Group& group,
                                   const Group::ExceedAction& exceed_action,
                                   const Group::ProductionCMode& newControl,
                                   DeferredLogger& deferred_logger);
    void actionOnBrokenConstraints(const Group& group,
                                   const Group::InjectionCMode& newControl,
                                   const Phase& controlPhase,
                                   DeferredLogger& deferred_logger);

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

    // create the well container
    virtual void createWellContainer(const int time_step) = 0;
    virtual void initWellContainer() = 0;

    virtual void calculateProductivityIndexValuesShutWells(const int reportStepIdx,
                                                           DeferredLogger& deferred_logger) = 0;
    virtual void calculateProductivityIndexValues(DeferredLogger& deferred_logger) = 0;

    void runWellPIScaling(const int timeStepIdx,
                          DeferredLogger& local_deferredLogger);

    Schedule& schedule_;
    const SummaryState& summaryState_;
    const EclipseState& eclState_;
    const Comm& comm_;

    PhaseUsage phase_usage_;
    bool terminal_output_{false};
    bool wells_active_{false};
    bool initial_step_{};
    bool report_step_starts_{};

    std::optional<int> last_run_wellpi_{};

    std::vector<Well> wells_ecl_;
    std::vector<std::vector<PerforationData>> well_perf_data_;
    std::function<bool(const Well&)> not_on_process_{};

    // a vector of all the wells.
    std::vector<WellInterfaceGeneric*> well_container_generic_{};

    std::vector<int> local_shut_wells_{};

    std::vector<ParallelWellInfo> parallel_well_info_;
    std::vector<ParallelWellInfo*> local_parallel_well_info_;

    std::vector<WellProdIndexCalculator> prod_index_calc_;

    // Map from logically cartesian cell indices to compressed ones.
    // Cells not in the interior are not mapped. This deactivates
    // these for distributed wells and makes the distribution non-overlapping.
    std::vector<int> cartesian_to_compressed_;

    std::vector<int> pvt_region_idx_;

    mutable std::unordered_set<std::string> closed_this_step_;

    WellTestState wellTestState_{};
    GuideRate guideRate_;
    std::unique_ptr<VFPProperties> vfp_properties_{};
    std::map<std::string, double> node_pressures_; // Storing network pressures for output.

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

  private:
    WellInterfaceGeneric* getGenWell(const std::string& well_name);
};


} // namespace Opm

#endif
