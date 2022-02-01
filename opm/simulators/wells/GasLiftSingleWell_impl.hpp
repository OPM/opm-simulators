/*
  Copyright 2020 Equinor ASA.

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

namespace Opm {

template<typename TypeTag>
GasLiftSingleWell<TypeTag>::
GasLiftSingleWell(const StdWell &std_well,
                  const Simulator &ebos_simulator,
                  const SummaryState &summary_state,
                  DeferredLogger &deferred_logger,
                  WellState &well_state,
                  const GroupState &group_state,
                  GasLiftGroupInfo &group_info,
                  GLiftSyncGroups &sync_groups,
                  bool glift_debug
                 )
    // The parent class GasLiftSingleWellGeneric contains all stuff
    //   that is not dependent on TypeTag
    : GasLiftSingleWellGeneric(
        deferred_logger,
        well_state,
        group_state,
        std_well.wellEcl(),
        summary_state,
        group_info,
        ebos_simulator.vanguard().schedule(),
        ebos_simulator.episodeIndex(),
        sync_groups,
        glift_debug
    )
   , ebos_simulator_{ebos_simulator}
   , std_well_{std_well}
{
    const auto& gl_well = *gl_well_;
    if(useFixedAlq_(gl_well)) {
        updateWellStateAlqFixedValue_(gl_well);
        this->optimize_ = false; // lift gas supply is fixed
    }
    else {
        setAlqMaxRate_(gl_well);
        this->optimize_ = true;
    }

    const auto& pu = std_well_.phaseUsage();
    this->oil_pos_ = pu.phase_pos[Oil];
    this->gas_pos_ = pu.phase_pos[Gas];
    this->water_pos_ = pu.phase_pos[Water];
    // get the alq value used for this well for the previous iteration (a
    //   nonlinear iteration in assemble() in BlackoilWellModel).
    //   If gas lift optimization has not been applied to this well yet, the
    //   default value is used.
    this->orig_alq_ = this->well_state_.getALQ(this->well_name_);
    if(this->optimize_) {
        setAlqMinRate_(gl_well);
        // NOTE: According to item 4 in WLIFTOPT, this value does not
        //    have to be positive.
        // TODO: Does it make sense to have a negative value?
        this->alpha_w_ = gl_well.weight_factor();
        if (this->alpha_w_ <= 0 ) {
            displayWarning_("Nonpositive value for alpha_w ignored");
            this->alpha_w_ = 1.0;
        }

        // NOTE: According to item 6 in WLIFTOPT:
        //   "If this value is greater than zero, the incremental gas rate will influence
        //    the calculation of the incremental gradient and may be used
        //    to discourage the allocation of lift gas to wells which produce more gas."
        // TODO: Does this mean that we should ignore this value if it
        //   is negative?
        this->alpha_g_ = gl_well.inc_weight_factor();

        // TODO: adhoc value.. Should we keep max_iterations_ as a safety measure
        //   or does it not make sense to have it?
        this->max_iterations_ = 1000;
    }
}

/****************************************
 * Private methods in alphabetical order
 ****************************************/

template<typename TypeTag>
GasLiftSingleWellGeneric::BasicRates
GasLiftSingleWell<TypeTag>::
computeWellRates_( double bhp, bool bhp_is_limited, bool debug_output ) const
{
    std::vector<double> potentials(this->num_phases_, 0.0);
    this->std_well_.computeWellRatesWithBhp(
        this->ebos_simulator_, bhp, potentials, this->deferred_logger_);
    if (debug_output) {
        const std::string msg = fmt::format("computed well potentials given bhp {}, "
            "oil: {}, gas: {}, water: {}", bhp,
            -potentials[this->oil_pos_], -potentials[this->gas_pos_],
            -potentials[this->water_pos_]);
        displayDebugMessage_(msg);
    }

    for (auto& potential : potentials) {
        potential = std::min(0.0, potential);
    }
    return {-potentials[this->oil_pos_],
            -potentials[this->gas_pos_],
            -potentials[this->water_pos_],
            bhp_is_limited
    };
}

template<typename TypeTag>
std::optional<double>
GasLiftSingleWell<TypeTag>::
computeBhpAtThpLimit_(double alq) const
{
    auto bhp_at_thp_limit = this->std_well_.computeBhpAtThpLimitProdWithAlq(
        this->ebos_simulator_,
        this->summary_state_,
        this->deferred_logger_,
        alq);
    if (bhp_at_thp_limit) {
        if (*bhp_at_thp_limit < this->controls_.bhp_limit) {
            const std::string msg = fmt::format(
                "Computed bhp ({}) from thp limit is below bhp limit ({}), (ALQ = {})."
                " Using bhp limit instead",
                *bhp_at_thp_limit, this->controls_.bhp_limit, alq);
            displayDebugMessage_(msg);
            bhp_at_thp_limit = this->controls_.bhp_limit;
        }
        //bhp_at_thp_limit = std::max(*bhp_at_thp_limit, this->controls_.bhp_limit);
    }
    else {
        const std::string msg = fmt::format(
            "Failed in getting converged bhp potential from thp limit (ALQ = {})", alq);
        displayDebugMessage_(msg);
    }
    return bhp_at_thp_limit;
}

template<typename TypeTag>
void
GasLiftSingleWell<TypeTag>::
setAlqMaxRate_(const GasLiftOpt::Well &well)
{
    auto& max_alq_optional = well.max_rate();
    if (max_alq_optional) {
        // NOTE: To prevent extrapolation of the VFP tables, any value
        // entered here must not exceed the largest ALQ value in the well's VFP table.
        this->max_alq_ = *max_alq_optional;
    }
    else { // i.e. WLIFTOPT, item 3 has been defaulted
        // According to the manual for WLIFTOPT, item 3:
        //   The default value should be set to the largest ALQ
        //   value in the well's VFP table
        const auto& table = std_well_.vfp_properties_->getProd()->getTable(
                this->controls_.vfp_table_number);
        const auto& alq_values = table.getALQAxis();
        // Assume the alq_values are sorted in ascending order, so
        // the last item should be the largest value:
        this->max_alq_ = alq_values.back();
    }
}

}
