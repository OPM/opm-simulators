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

#ifndef OPM_GASLIFT_SINGLE_WELL_IMPL_HEADER_INCLUDED
#define OPM_GASLIFT_SINGLE_WELL_IMPL_HEADER_INCLUDED

// Improve IDE experience
#ifndef OPM_GASLIFT_SINGLE_WELL_HEADER_INCLUDED
#include <config.h>
#include <opm/simulators/wells/GasLiftSingleWell.hpp>
#endif

#include <opm/common/TimingMacros.hpp>

#include <opm/input/eclipse/Schedule/GasLiftOpt.hpp>
#include <opm/input/eclipse/Schedule/Well/Well.hpp>

#include <string>
#include <vector>

#include <fmt/format.h>

namespace Opm {

template<typename TypeTag>
GasLiftSingleWell<TypeTag>::
GasLiftSingleWell(const WellInterface<TypeTag>& well,
                  const Simulator& simulator,
                  const SummaryState& summary_state,
                  DeferredLogger& deferred_logger,
                  WellState<Scalar, IndexTraits>& well_state,
                  const GroupState<Scalar>& group_state,
                  GasLiftGroupInfo<Scalar, IndexTraits>& group_info,
                  GLiftSyncGroups &sync_groups,
                  const Parallel::Communication& comm,
                  bool glift_debug)
    // The parent class GasLiftSingleWellGeneric contains all stuff
    //   that is not dependent on TypeTag
    : GasLiftSingleWellGeneric<Scalar, IndexTraits>(deferred_logger,
                                       well_state,
                                       group_state,
                                       well.wellEcl(),
                                       summary_state,
                                       group_info,
                                       simulator.vanguard().schedule(),
                                       simulator.episodeIndex(),
                                       sync_groups,
                                       comm,
                                       glift_debug)
   , simulator_{simulator}
   , well_{well}
{
    const auto& gl_well = *this->gl_well_;
    if (this->useFixedAlq_(gl_well)) {
        this->updateWellStateAlqFixedValue_(gl_well);
        this->optimize_ = false; // lift gas supply is fixed
    }
    else {
        setAlqMaxRate_(gl_well);
        this->optimize_ = true;
    }

    setupPhaseVariables_();
    // get the alq value used for this well for the previous iteration (a
    //   nonlinear iteration in assemble() in BlackoilWellModel).
    //   If gas lift optimization has not been applied to this well yet, the
    //   default value is used.
    this->orig_alq_ = this->well_state_.well(this->well_name_).alq_state.get();
    if (this->optimize_) {
        this->setAlqMinRate_(gl_well);
        // NOTE: According to item 4 in WLIFTOPT, this value does not
        //    have to be positive.
        // TODO: Does it make sense to have a negative value?
        this->alpha_w_ = gl_well.weight_factor();
        if (this->alpha_w_ <= 0 ) {
            this->displayWarning_("Nonpositive value for alpha_w ignored");
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
typename GasLiftSingleWell<TypeTag>::BasicRates
GasLiftSingleWell<TypeTag>::
computeWellRates_(Scalar bhp, bool bhp_is_limited, bool debug_output ) const
{
    std::vector<Scalar> potentials(this->NUM_PHASES, 0.0);
    this->well_.computeWellRatesWithBhp(this->simulator_,
                                        bhp,
                                        potentials,
                                        this->deferred_logger_);
    if (debug_output) {
        const std::string msg = fmt::format("computed well potentials given bhp {}, "
            "oil: {}, gas: {}, water: {}", bhp,
            -potentials[this->oil_pos_], -potentials[this->gas_pos_],
            -potentials[this->water_pos_]);
        this->displayDebugMessage_(msg);
    }

    std::transform(potentials.begin(), potentials.end(),
                   potentials.begin(),
                   [](const auto& potential)
                   { return std::min(Scalar{0}, potential); });
    return {-potentials[this->oil_pos_],
            -potentials[this->gas_pos_],
            -potentials[this->water_pos_],
            bhp_is_limited
    };
}

template<typename TypeTag>
std::optional<typename GasLiftSingleWell<TypeTag>::Scalar>
GasLiftSingleWell<TypeTag>::
computeBhpAtThpLimit_(Scalar alq, bool debug_output) const
{
    OPM_TIMEFUNCTION();
    auto bhp_at_thp_limit = this->well_.computeBhpAtThpLimitProdWithAlq(
        this->simulator_,
        this->summary_state_,
        alq,
        this->deferred_logger_,
        /*iterate_if_no_solution */ false);
    if (bhp_at_thp_limit) {
        if (*bhp_at_thp_limit < this->controls_.bhp_limit) {
            if (debug_output && this->debug) {
                const std::string msg = fmt::format(
                    "Computed bhp ({}) from thp limit is below bhp limit ({}), (ALQ = {})."
                    " Using bhp limit instead",
                    *bhp_at_thp_limit, this->controls_.bhp_limit, alq
                );
                this->displayDebugMessage_(msg);
            }
            bhp_at_thp_limit = this->controls_.bhp_limit;
        }
        //bhp_at_thp_limit = std::max(*bhp_at_thp_limit, this->controls_.bhp_limit);
    }
    else {
        const std::string msg = fmt::format(
            "Failed in getting converged bhp potential from thp limit (ALQ = {})", alq);
        this->displayDebugMessage_(msg);
    }
    return bhp_at_thp_limit;
}

template<typename TypeTag>
void
GasLiftSingleWell<TypeTag>::
setupPhaseVariables_()
{
#ifndef NDEBUG
    bool num_phases_ok = (FluidSystem::numActivePhases()== 3);
#endif
    if (FluidSystem::numActivePhases()== 2) {
        // NOTE: We support two-phase oil-water flow, by setting the gas flow rate
        //   to zero. This is done by initializing the potential vector to zero:
        //
        //     std::vector<double> potentials(NUM_PHASES, 0.0);
        //
        //  see e.g. runOptimizeLoop_() in GasLiftSingleWellGeneric.cpp
        //  In addition the VFP calculations, e.g. to calculate BHP from THP
        //  has been adapted to the two-phase oil-water case, see the comment
        //  in WellInterfaceGeneric.cpp for the method adaptRatesForVFP() for
        //  more information.
        if (    FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx)
             && FluidSystem::phaseIsActive(FluidSystem::oilPhaseIdx)
             && !FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx) )
        {
#ifndef NDEBUG
            num_phases_ok = true;  // two-phase oil-water is also supported
#endif
        }
        else {
            throw std::logic_error("Two-phase gas lift optimization only supported"
                                  " for oil and water");
        }
    }
    assert(num_phases_ok);
    this->oil_pos_ = FluidSystem::canonicalToActivePhaseIdx(FluidSystem::oilPhaseIdx);
    this->gas_pos_ = FluidSystem::canonicalToActivePhaseIdx(FluidSystem::gasPhaseIdx);
    this->water_pos_ = FluidSystem::canonicalToActivePhaseIdx(FluidSystem::waterPhaseIdx);
}

template<typename TypeTag>
void
GasLiftSingleWell<TypeTag>::
setAlqMaxRate_(const GasLiftWell& well)
{
    const auto& max_alq_optional = well.max_rate();
    if (max_alq_optional) {
        // NOTE: To prevent extrapolation of the VFP tables, any value
        // entered here must not exceed the largest ALQ value in the well's VFP table.
        this->max_alq_ = *max_alq_optional;
    }
    else { // i.e. WLIFTOPT, item 3 has been defaulted
        // According to the manual for WLIFTOPT, item 3:
        //   The default value should be set to the largest ALQ
        //   value in the well's VFP table
        const auto& table = well_.vfpProperties()->getProd()->getTable(
                this->controls_.vfp_table_number);
        const auto& alq_values = table.getALQAxis();
        // Assume the alq_values are sorted in ascending order, so
        // the last item should be the largest value:
        this->max_alq_ = alq_values.back();
    }
}

template<typename TypeTag>
bool
GasLiftSingleWell<TypeTag>::
checkThpControl_() const
{
    const int well_index = this->well_state_.index(this->well_name_).value();
    const Well::ProducerCMode& control_mode =
                         this->well_state_.well(well_index).production_cmode;
    bool thp_control = control_mode == Well::ProducerCMode::THP;
    const auto& well = getWell();
    thp_control = thp_control || well.thpLimitViolatedButNotSwitched();
    if (this->debug) {
        if (!thp_control) {
            this->displayDebugMessage_("Well is not under THP control, skipping iteration..");
        }
    }
    return thp_control;
}

}

#endif
