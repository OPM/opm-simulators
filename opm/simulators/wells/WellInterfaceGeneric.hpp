/*
  Copyright 2017 SINTEF Digital, Mathematics and Cybernetics.
  Copyright 2017 Statoil ASA.
  Copyright 2017 IRIS
  Copyright 2019 Norce

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


#ifndef OPM_WELLINTERFACE_GENERIC_HEADER_INCLUDED
#define OPM_WELLINTERFACE_GENERIC_HEADER_INCLUDED

#include <opm/input/eclipse/Schedule/Well/Well.hpp>
#include <opm/simulators/flow/BlackoilModelParameters.hpp>
#include <opm/simulators/wells/RuntimePerforation.hpp>

#include <map>
#include <optional>
#include <string>
#include <vector>

namespace Opm
{

class DeferredLogger;
class GuideRate;
template<class Scalar> class ParallelWellInfo;
template<class Scalar> struct PerforationData;
struct PhaseUsage;
class SummaryState;
template<typename FluidSystem, typename Indices> class VFPProperties;
class WellTestState;
template<typename FluidSystem, typename Indices> class WellState;
template<typename FluidSystem, typename Indices> class SingleWellState;
class Group;
class Schedule;

template<typename FluidSystem, typename Indices>
class WellInterfaceGeneric {
public:
    using Scalar = typename FluidSystem::Scalar;
    using ModelParameters = BlackoilModelParameters<Scalar>;

    using WellStateType = WellState<FluidSystem, Indices>;

    WellInterfaceGeneric(const Well& well,
                         const ParallelWellInfo<Scalar>& parallel_well_info,
                         const int time_step,
                         const ModelParameters& param,
                         const int pvtRegionIdx,
                         const int num_components,
                         const int num_phases,
                         const int index_of_well,
                         const std::vector<PerforationData<Scalar>>& perf_data);

    /// \brief Get the perforations of the well
    const std::vector<PerforationData<Scalar>>& perforationData() const;

    /// Well name.
    const std::string& name() const;

    /// True if the well is an injector.
    bool isInjector() const;

    /// True if the well is a producer.
    bool isProducer() const;

    /// Well cells.
    const std::vector<int>& cells() const { return well_cells_; }

    /// Index of well in the wells struct and wellState
    int indexOfWell() const;

    void adaptRatesForVFP(std::vector<Scalar>& rates) const;

    const Well& wellEcl() const;
    Well& wellEcl();
    const PhaseUsage& phaseUsage() const;

    /// Returns true if the well is currently in prediction mode (i.e. not history mode).
    bool underPredictionMode() const;

    // whether the well is operable
    bool isOperableAndSolvable() const;
    bool useVfpExplicit () const;
    bool thpLimitViolatedButNotSwitched() const;

    void initCompletions();
    void closeCompletions(const WellTestState& wellTestState);

    void setVFPProperties(const VFPProperties<FluidSystem, Indices>* vfp_properties_arg);
    void setPrevSurfaceRates(WellStateType& well_state,
                             const WellStateType& prev_well_state) const;
    void setGuideRate(const GuideRate* guide_rate_arg);
    void setWellEfficiencyFactor(const Scalar efficiency_factor);
    void setRepRadiusPerfLength();
    void setWsolvent(const Scalar wsolvent);
    void setDynamicThpLimit(const Scalar thp_limit);
    std::optional<Scalar> getDynamicThpLimit() const;
    void setDynamicThpLimit(const std::optional<Scalar> thp_limit);
    void updatePerforatedCell(std::vector<bool>& is_cell_perforated);

    /// Returns true if the well has one or more THP limits/constraints.
    bool wellHasTHPConstraints(const SummaryState& summaryState) const;

    void stopWell() { this->wellStatus_ = Well::Status::STOP; }
    void openWell() { this->wellStatus_ = Well::Status::OPEN; }

    bool wellIsStopped() const { return this->wellStatus_ == Well::Status::STOP; }

    int currentStep() const { return this->current_step_; }

    int pvtRegionIdx() const { return pvtRegionIdx_; }

    const GuideRate* guideRate() const { return guide_rate_; }

    int numComponents() const { return num_components_; }

    int numPhases() const { return number_of_phases_; }

    int numLocalPerfs() const { return number_of_local_perforations_; }

    Scalar refDepth() const { return ref_depth_; }

    Scalar gravity() const { return gravity_; }

    const VFPProperties<FluidSystem, Indices>* vfpProperties() const { return vfp_properties_; }

    const ParallelWellInfo<Scalar>& parallelWellInfo() const { return parallel_well_info_; }

    const std::vector<Scalar>& perfDepth() const { return perf_depth_; }

    std::vector<Scalar>& perfDepth() { return perf_depth_; }

    const std::vector<Scalar>& wellIndex() const { return well_index_; }

    const std::map<int,std::vector<int>>& getCompletions() const { return completions_; }

    Scalar getTHPConstraint(const SummaryState& summaryState) const;
    Scalar getALQ(const WellStateType& well_state) const;
    Scalar wsolvent() const;
    Scalar rsRvInj() const;

    // at the beginning of the time step, we check what inj_multiplier from the previous running
    void initInjMult(const std::vector<Scalar>& max_inj_mult);

    // update the InjMult information at the end of the time step, so it can be used for later.
    void updateInjMult(std::vector<Scalar>& inj_multipliers,
                       DeferredLogger& deferred_logger) const;

    // Note:: for multisegment wells, bhp is actually segment pressure in practice based on observation
    // it might change in the future
    Scalar getInjMult(const int local_perf_index, const Scalar bhp, const Scalar perf_pres, DeferredLogger& dlogger) const;

    // whether a well is specified with a non-zero and valid VFP table number
    bool isVFPActive(DeferredLogger& deferred_logger) const;

    void reportWellSwitching(const SingleWellState<FluidSystem, Indices>& ws,
                             DeferredLogger& deferred_logger) const;

    bool changedToOpenThisStep() const { return this->changed_to_open_this_step_; }

    void updateWellTestState(const SingleWellState<FluidSystem, Indices>& ws,
                             const double& simulationTime,
                             const bool& writeMessageToOPMLog,
                             const bool zero_group_target,
                             WellTestState& wellTestState,
                             DeferredLogger& deferred_logger) const;

    bool isPressureControlled(const WellStateType& well_state) const;

    Scalar wellEfficiencyFactor() const { return well_efficiency_factor_; }

    //! \brief Update filter cake multipliers.
    void updateFilterCakeMultipliers(const std::vector<Scalar>& inj_fc_multiplier)
    {
        inj_fc_multiplier_ = inj_fc_multiplier;
    }

    void resetWellOperability();

    virtual std::vector<Scalar> getPrimaryVars() const
    {
        return {};
    }

    virtual int setPrimaryVars(typename std::vector<Scalar>::const_iterator)
    {
        return 0;
    }

    virtual Scalar connectionDensity(const int globalConnIdx,
                                     const int openConnIdx) const = 0;

    void addPerforations(const std::vector<RuntimePerforation>& perfs);

protected:
    bool getAllowCrossFlow() const;

    Scalar wmicrobes_() const;
    Scalar wfoam_() const;
    Scalar woxygen_() const;
    Scalar wpolymer_() const;
    Scalar wsalt_() const;
    Scalar wurea_() const;

    int polymerTable_() const;
    int polymerInjTable_() const;
    int polymerWaterTable_() const;

    bool wellUnderZeroRateTargetIndividual(const SummaryState& summary_state,
                                           const WellState<FluidSystem, Indices>& well_state) const;

    bool wellUnderGroupControl(const SingleWellState<FluidSystem, Indices>& ws) const;

    std::pair<bool,bool>
    computeWellPotentials(std::vector<Scalar>& well_potentials,
                          const WellStateType& well_state);

    void checkNegativeWellPotentials(std::vector<Scalar>& well_potentials,
                                     const bool checkOperability,
                                     DeferredLogger& deferred_logger);

    void prepareForPotentialCalculations(const SummaryState& summary_state,
                                         WellStateType& well_state,
                                         Well::InjectionControls& inj_controls,
                                         Well::ProductionControls& prod_controls) const;

    void resetDampening() {
        std::fill(this->inj_multiplier_damp_factor_.begin(), this->inj_multiplier_damp_factor_.end(), 1.0);
    }

    // definition of the struct OperabilityStatus
    struct OperabilityStatus
    {
        bool isOperableAndSolvable() const
        {
            if (!operable_under_only_bhp_limit || !solvable || has_negative_potentials) {
                return false;
            } else {
                return ( (isOperableUnderBHPLimit() || isOperableUnderTHPLimit()) );
            }
        }

        bool isOperableUnderBHPLimit() const
        {
            return operable_under_only_bhp_limit && obey_thp_limit_under_bhp_limit;
        }

        bool isOperableUnderTHPLimit() const
        {
            return can_obtain_bhp_with_thp_limit && obey_bhp_limit_with_thp_limit;
        }

        void resetOperability()
        {
            operable_under_only_bhp_limit = true;
            obey_thp_limit_under_bhp_limit = true;
            can_obtain_bhp_with_thp_limit = true;
            obey_bhp_limit_with_thp_limit = true;
        }

        // whether the well can be operated under bhp limit
        // without considering other limits.
        // if it is false, then the well is not operable for sure.
        bool operable_under_only_bhp_limit = true;
        // if the well can be operated under bhp limit, will it obey(not violate)
        // the thp limit when operated under bhp limit
        bool obey_thp_limit_under_bhp_limit = true;
        // whether the well operate under the thp limit only
        bool can_obtain_bhp_with_thp_limit = true;
        // whether the well obey bhp limit when operated under thp limit
        bool obey_bhp_limit_with_thp_limit = true;
        // the well is solveable
        bool solvable = true;
        // the well have non positive potentials
        bool has_negative_potentials = false;
        //thp limit violated but not switched
        mutable bool thp_limit_violated_but_not_switched = false;

        bool use_vfpexplicit = false;
    };

    OperabilityStatus operability_status_;

    Well well_ecl_;

    const ParallelWellInfo<Scalar>& parallel_well_info_;
    const int current_step_;
    const ModelParameters& param_;

    // The pvt region of the well. We assume
    // We assume a well to not penetrate more than one pvt region.
    const int pvtRegionIdx_;

    const int num_components_;

    // number of phases
    int number_of_phases_;

    // the index of well in Wells struct
    int index_of_well_;

    const std::vector<PerforationData<Scalar>>* perf_data_;

    // the vectors used to describe the inflow performance relationship (IPR)
    // Q = IPR_A - BHP * IPR_B
    // TODO: it minght need to go to WellInterface, let us implement it in StandardWell first
    // it is only updated and used for producers for now
    mutable std::vector<Scalar> ipr_a_;
    mutable std::vector<Scalar> ipr_b_;

    // cell index for each well perforation
    std::vector<int> well_cells_;

    // well index for each perforation
    std::vector<Scalar> well_index_;

    // number of the perforations for this well on this process
    int number_of_local_perforations_;

    // depth for each perforation
    std::vector<Scalar> perf_depth_;

    // representative radius of the perforations, used in shear calculation
    std::vector<Scalar> perf_rep_radius_;

    // length of the perforations, use in shear calculation
    std::vector<Scalar> perf_length_;

    // well bore diameter
    std::vector<Scalar> bore_diameters_;

    /*
     *  completions_ contains the mapping from completion id to connection indices
     *  {
     *      2 : [ConnectionIndex, ConnectionIndex],
     *      1 : [ConnectionIndex, ConnectionIndex, ConnectionIndex],
     *      5 : [ConnectionIndex],
     *      7 : [ConnectionIndex]
     *      ...
     *   }
     *   The integer IDs correspond to the COMPLETION id given by the COMPLUMP keyword.
     *   When there is no COMPLUMP keyword used, a default completion number will be assigned
     *   based on the order of the declaration of the connections.
     *   Since the connections not OPEN is not included in the Wells, so they will not be considered
     *   in this mapping relation.
     */
    std::map<int, std::vector<int>> completions_;

    // reference depth for the BHP
    Scalar ref_depth_;

    // saturation table nubmer for each well perforation
    std::vector<int> saturation_table_number_;

    Well::Status wellStatus_;

    const PhaseUsage* phase_usage_;

    Scalar gravity_;
    Scalar wsolvent_;
    std::optional<Scalar> dynamic_thp_limit_;

    // recording the multiplier calculate from the keyword WINJMULT during the time step
    mutable std::vector<Scalar> inj_multiplier_;

    // the injection multiplier from the previous running, it is mostly used for CIRR mode
    // which intends to keep the fracturing open
    std::vector<Scalar> prev_inj_multiplier_;

    // WINJMULT multipliers for previous iteration (used for oscillation detection)
    mutable std::vector<Scalar> inj_multiplier_previter_;
    // WINJMULT dampening factors (used in case of oscillations)
    mutable std::vector<Scalar> inj_multiplier_damp_factor_;

    // the multiplier due to injection filtration cake
    std::vector<Scalar> inj_fc_multiplier_;

    Scalar well_efficiency_factor_;
    const VFPProperties<FluidSystem, Indices>* vfp_properties_;
    const GuideRate* guide_rate_;

    std::vector<std::string> well_control_log_;

    bool changed_to_open_this_step_ = true;
};

}

#endif // OPM_WELLINTERFACE_GENERIC_HEADER_INCLUDED
