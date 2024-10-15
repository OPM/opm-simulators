/*
  Copyright 2017 SINTEF Digital, Mathematics and Cybernetics.
  Copyright 2017 Statoil ASA.
  Copyright 2018 IRIS

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
#include <opm/simulators/wells/WellInterfaceGeneric.hpp>

#include <opm/common/ErrorMacros.hpp>

#include <opm/input/eclipse/Schedule/Schedule.hpp>
#include <opm/input/eclipse/Schedule/Well/FilterCake.hpp>
#include <opm/input/eclipse/Schedule/Well/WellBrineProperties.hpp>
#include <opm/input/eclipse/Schedule/Well/WellConnections.hpp>
#include <opm/input/eclipse/Schedule/Well/WellFoamProperties.hpp>
#include <opm/input/eclipse/Schedule/Well/WellMICPProperties.hpp>
#include <opm/input/eclipse/Schedule/Well/WellPolymerProperties.hpp>
#include <opm/input/eclipse/Schedule/Well/WellTestState.hpp>
#include <opm/input/eclipse/Schedule/Well/WVFPEXP.hpp>
#include <opm/simulators/utils/DeferredLoggingErrorHelpers.hpp>
#include <opm/simulators/wells/PerforationData.hpp>
#include <opm/simulators/wells/ParallelWellInfo.hpp>
#include <opm/simulators/wells/VFPHelpers.hpp>
#include <opm/simulators/wells/VFPProperties.hpp>
#include <opm/simulators/wells/WellBhpThpCalculator.hpp>
#include <opm/simulators/wells/WellHelpers.hpp>
#include <opm/simulators/wells/WellState.hpp>
#include <opm/simulators/wells/WellTest.hpp>

#include <fmt/format.h>

#include <cassert>
#include <cmath>
#include <cstddef>
#include <stdexcept>

namespace Opm {

template<class Scalar>
WellInterfaceGeneric<Scalar>::
WellInterfaceGeneric(const Well& well,
                     const ParallelWellInfo<Scalar>& pw_info,
                     const int time_step,
                     const ModelParameters& param,
                     const int pvtRegionIdx,
                     const int num_components,
                     const int num_phases,
                     const int index_of_well,
                     const std::vector<PerforationData<Scalar>>& perf_data)
      : well_ecl_(well)
      , parallel_well_info_(pw_info)
      , current_step_(time_step)
      , param_(param)
      , pvtRegionIdx_(pvtRegionIdx)
      , num_components_(num_components)
      , number_of_phases_(num_phases)
      , index_of_well_(index_of_well)
      , perf_data_(&perf_data)
      , ipr_a_(num_components)
      , ipr_b_(num_components)
{
    assert(well.name()==pw_info.name());
    assert(std::is_sorted(perf_data.begin(), perf_data.end(),
                          [](const auto& perf1, const auto& perf2){
        return perf1.ecl_index < perf2.ecl_index;
    }));
    if (time_step < 0) {
        OPM_THROW(std::invalid_argument, "Negative time step is used to construct WellInterface");
    }

    ref_depth_ = well.getRefDepth();

    // We do not want to count SHUT perforations here, so
    // it would be wrong to use wells.getConnections().size().
    // This is the number_of_perforations_ on this process only!
    number_of_perforations_ = perf_data.size();

    // perforations related
    {
        well_cells_.resize(number_of_perforations_);
        well_index_.resize(number_of_perforations_);
        saturation_table_number_.resize(number_of_perforations_);
        int perf = 0;
        for (const auto& pd : perf_data) {
            well_cells_[perf] = pd.cell_index;
            well_index_[perf] = pd.connection_transmissibility_factor;
            saturation_table_number_[perf] = pd.satnum_id;
            ++perf;
        }
    }

    // initialization of the completions mapping
    initCompletions();

    well_efficiency_factor_ = 1.0;

    this->wellStatus_ = Well::Status::OPEN;
    if (well.getStatus() == Well::Status::STOP) {
        this->wellStatus_ = Well::Status::STOP;
    }

    wsolvent_ = 0.0;

    well_control_log_.clear();
}

// Currently the VFP calculations requires three-phase input data, see
//  the documentation for keyword VFPPROD and its implementation in
//  VFPProdProperties.cpp. However, by setting the gas flow rate to a dummy
//  value in VFPPROD record 5 (GFR values) and supplying a dummy input value
//  for the gas rate to the methods in VFPProdProperties.cpp, we can extend
//  the VFP calculations to the two-phase oil-water case.
template<class Scalar>
void WellInterfaceGeneric<Scalar>::
adaptRatesForVFP(std::vector<Scalar>& rates) const
{
    const auto& pu = this->phaseUsage();
    if (pu.num_phases == 2) {
        if (    pu.phase_used[BlackoilPhases::Aqua] == 1
             && pu.phase_used[BlackoilPhases::Liquid] == 1
             && pu.phase_used[BlackoilPhases::Vapour] == 0)
        {
            assert(rates.size() == 2);
            rates.push_back(0.0);  // set gas rate to zero
        }
        else {
            throw std::logic_error("Two-phase VFP calculation only "
                                   "supported for oil and water");
        }
    }
}

template<class Scalar>
const std::vector<PerforationData<Scalar>>&
WellInterfaceGeneric<Scalar>::perforationData() const
{
    return *perf_data_;
}

template<class Scalar>
const std::string&
WellInterfaceGeneric<Scalar>::name() const
{
    return well_ecl_.name();
}

template<class Scalar>
bool WellInterfaceGeneric<Scalar>::isInjector() const
{
    return well_ecl_.isInjector();
}

template<class Scalar>
bool WellInterfaceGeneric<Scalar>::isProducer() const
{
    return well_ecl_.isProducer();
}

template<class Scalar>
int WellInterfaceGeneric<Scalar>::indexOfWell() const
{
    return index_of_well_;
}

template<class Scalar>
bool WellInterfaceGeneric<Scalar>::getAllowCrossFlow() const
{
    return well_ecl_.getAllowCrossFlow();
}

template<class Scalar>
const Well& WellInterfaceGeneric<Scalar>::wellEcl() const
{
    return well_ecl_;
}

template<class Scalar>
Well& WellInterfaceGeneric<Scalar>::wellEcl()
{
    return well_ecl_;
}

template<class Scalar>
const PhaseUsage& WellInterfaceGeneric<Scalar>::phaseUsage() const
{
    assert(phase_usage_ != nullptr);

    return *phase_usage_;
}

template<class Scalar>
Scalar WellInterfaceGeneric<Scalar>::wsolvent() const
{
    return wsolvent_;
}

template<class Scalar>
Scalar WellInterfaceGeneric<Scalar>::rsRvInj() const
{
    return well_ecl_.getInjectionProperties().rsRvInj;
}

template<class Scalar>
void WellInterfaceGeneric<Scalar>::
initInjMult(const std::vector<Scalar>& max_inj_mult)
{
    // prev_inj_multiplier_ will stay unchanged during the time step
    // while inj_multiplier_ might be updated during the time step
    this->prev_inj_multiplier_ = max_inj_mult;
    const auto nperf = max_inj_mult.size();
    // initializing the (report?) step specific multipliers and damping factors to be 1.0
    this->inj_multiplier_ = std::vector<Scalar>(nperf, 1.);
    this->inj_multiplier_previter_ = std::vector<Scalar>(nperf, 1.);
    this->inj_multiplier_damp_factor_ = std::vector<Scalar>(nperf, 1.);
}

template<class Scalar>
void WellInterfaceGeneric<Scalar>::
updateInjMult(std::vector<Scalar>& inj_multipliers,
              DeferredLogger& deferred_logger) const
{
    if (inj_multipliers.size() != this->inj_multiplier_.size()) {
        OPM_DEFLOG_THROW(std::runtime_error,
                         fmt::format("We do not support changing connection numbers during simulation with WINJMULT "
                                     "for well {}", name()),
                         deferred_logger);
    }

    inj_multipliers = this->inj_multiplier_;
}

template<class Scalar>
Scalar WellInterfaceGeneric<Scalar>::
getInjMult(const int perf,
           const Scalar bhp,
           const Scalar perf_pres,
           DeferredLogger& dlogger) const
{
    assert(!this->isProducer());

    Scalar multiplier = 1.;

    const auto perf_ecl_index = this->perforationData()[perf].ecl_index;
    const bool is_wrev = this->well_ecl_.getInjMultMode() == Well::InjMultMode::WREV;

    const bool active_injmult = (is_wrev && this->well_ecl_.aciveWellInjMult()) ||
                              this->well_ecl_.getConnections()[perf_ecl_index].activeInjMult();

    if (active_injmult) {
        const auto& injmult= is_wrev ? this->well_ecl_.getWellInjMult()
                                                  : this->well_ecl_.getConnections()[perf_ecl_index].injmult();
        const Scalar pres = is_wrev ? bhp : perf_pres;

        const auto frac_press = injmult.fracture_pressure;
        const auto gradient = injmult.multiplier_gradient;
        if (pres > frac_press) {
            multiplier = 1. + (pres - frac_press) * gradient;
        }
        const Scalar osc_threshold = param_.inj_mult_osc_threshold_;
        const Scalar prev_multiplier = this->inj_multiplier_[perf_ecl_index] > 0 ? this->inj_multiplier_[perf_ecl_index] : 1.0;
        if (std::abs(multiplier - prev_multiplier) > osc_threshold) {
            const Scalar prev2_multiplier = this->inj_multiplier_previter_[perf_ecl_index] > 0 ? this->inj_multiplier_previter_[perf_ecl_index] : 1.0;
            const bool oscillating = (multiplier - prev_multiplier) * (prev_multiplier - prev2_multiplier) < 0;

            const Scalar min_damp_factor = this->param_.inj_mult_min_damp_factor_;
            Scalar damp_factor = this->inj_multiplier_damp_factor_[perf_ecl_index];
            if (oscillating) {
                if (damp_factor > min_damp_factor) {
                    const Scalar damp_multiplier = this->param_.inj_mult_damp_mult_;
                    damp_factor = std::max(min_damp_factor, damp_factor * damp_multiplier);
                    this->inj_multiplier_damp_factor_[perf_ecl_index] = damp_factor;
                } else {
                    const auto msg = fmt::format("Well {}, perforation {}: Injectivity multiplier dampening factor reached minimum (= {})",
                                                 this->name(),
                                                 perf_ecl_index,
                                                 min_damp_factor);
                    dlogger.warning(msg);
                }
            }
            multiplier = multiplier * damp_factor + (1.0 - damp_factor) * prev_multiplier;
        }
        // for CIRR mode, if there is no active WINJMULT setup, we will use the previous injection multiplier,
        // to mimic keeping the existing fracturing open
        if (this->well_ecl_.getInjMultMode() == Well::InjMultMode::CIRR) {
            multiplier = std::max(multiplier, this->prev_inj_multiplier_[perf_ecl_index]);
        }

        this->inj_multiplier_[perf_ecl_index] = multiplier;
        this->inj_multiplier_previter_[perf_ecl_index] = prev_multiplier;
    }

    return multiplier;
}

template<class Scalar>
bool WellInterfaceGeneric<Scalar>::
wellHasTHPConstraints(const SummaryState& summaryState) const
{
    // only wells under prediction mode can have THP constraint
    if (!this->wellEcl().predictionMode()) {
        return false;
    }

    if (dynamic_thp_limit_) {
        return true;
    }

    return WellBhpThpCalculator(*this).wellHasTHPConstraints(summaryState);
}

template<class Scalar>
void WellInterfaceGeneric<Scalar>::
updateWellTestState(const SingleWellState<Scalar>& ws,
                    const double& simulationTime,
                    const bool& writeMessageToOPMLog,
                    const bool zero_group_target,
                    WellTestState& wellTestState,
                    DeferredLogger& deferred_logger) const
{
    // updating well test state based on Economic limits for operable wells
    if (this->isOperableAndSolvable()) {
        WellTest(*this).updateWellTestStateEconomic(ws, simulationTime, writeMessageToOPMLog, wellTestState,
                                                    zero_group_target, deferred_logger);
    } else {
        // updating well test state based on physical (THP/BHP) limits.
        WellTest(*this).updateWellTestStatePhysical(simulationTime, writeMessageToOPMLog, wellTestState, deferred_logger);
    }

    // TODO: well can be shut/closed due to other reasons
}

template<class Scalar>
Scalar WellInterfaceGeneric<Scalar>::
getTHPConstraint(const SummaryState& summaryState) const
{
    if (dynamic_thp_limit_) {
        return *dynamic_thp_limit_;
    }

    return WellBhpThpCalculator(*this).getTHPConstraint(summaryState);
}

template<class Scalar>
bool WellInterfaceGeneric<Scalar>::underPredictionMode() const
{
    return well_ecl_.predictionMode();
}

template<class Scalar>
void WellInterfaceGeneric<Scalar>::initCompletions()
{
    assert(completions_.empty() );

    const WellConnections& connections = well_ecl_.getConnections();
    const std::size_t num_conns = connections.size();

    int num_active_connections = 0;
    auto my_next_perf = perf_data_->begin();
    for (std::size_t c = 0; c < num_conns; ++c) {
        if (my_next_perf == perf_data_->end())
        {
            break;
        }
        if (my_next_perf->ecl_index > c)
        {
            continue;
        }
        assert(my_next_perf->ecl_index == c);
        if (connections[c].state() == Connection::State::OPEN) {
            completions_[connections[c].complnum()].push_back(num_active_connections++);
        }
        ++my_next_perf;
    }
    assert(my_next_perf == perf_data_->end());
}

template<class Scalar>
void WellInterfaceGeneric<Scalar>::
closeCompletions(const WellTestState& wellTestState)
{
    const auto& connections = well_ecl_.getConnections();
    int perfIdx = 0;
    for (const auto& connection : connections) {
        if (connection.state() == Connection::State::OPEN) {
            if (wellTestState.completion_is_closed(name(), connection.complnum())) {
                this->well_index_[perfIdx] = 0.0;
            }
            perfIdx++;
        }
    }
}

template<class Scalar>
void WellInterfaceGeneric<Scalar>::
setVFPProperties(const VFPProperties<Scalar>* vfp_properties_arg)
{
    vfp_properties_ = vfp_properties_arg;
}

template<class Scalar>
void WellInterfaceGeneric<Scalar>::
setPrevSurfaceRates(WellState<Scalar>& well_state,
                    const WellState<Scalar>& prev_well_state) const
{
    auto& ws = well_state.well(this->index_of_well_);
    auto& ws_prev = prev_well_state.well(this->index_of_well_);
    // The logic here is a bit fragile:
    // We need non-zero prev_surface_rates for the purpose of providing explicit fractions
    // (if needed) for vfp interpolation.
    // We assume that current surface rates either are initialized from previous step
    // or (if newly opened) from updateWellStateRates. This is fine unless well was
    // stopped in previous step in which case it's rates will be zero. In this case,
    // we select the previous rates of the previous well state (and hope for the best).
    const bool zero_rates = std::all_of(ws.surface_rates.begin(), ws.surface_rates.end(),
            [](Scalar rate) {
                return rate == 0.; // TODO: should we use a threshhold for comparison?
            } );

    if (zero_rates) {
        ws.prev_surface_rates = ws_prev.prev_surface_rates;
    } else {
        ws.prev_surface_rates = ws.surface_rates;
    }
}

template<class Scalar>
void WellInterfaceGeneric<Scalar>::
setGuideRate(const GuideRate* guide_rate_arg)
{
    guide_rate_ = guide_rate_arg;
}

template<class Scalar>
void WellInterfaceGeneric<Scalar>::
setWellEfficiencyFactor(const Scalar efficiency_factor)
{
    well_efficiency_factor_ = efficiency_factor;
}

template<class Scalar>
void WellInterfaceGeneric<Scalar>::setRepRadiusPerfLength()
{
    const int nperf = number_of_perforations_;

    perf_rep_radius_.clear();
    perf_length_.clear();
    bore_diameters_.clear();

    perf_rep_radius_.reserve(nperf);
    perf_length_.reserve(nperf);
    bore_diameters_.reserve(nperf);

    const WellConnections& connections = well_ecl_.getConnections();
    const std::size_t num_conns = connections.size();
    int num_active_connections = 0;
    auto my_next_perf = perf_data_->begin();
    for (std::size_t c = 0; c < num_conns; ++c) {
        if (my_next_perf == perf_data_->end())
        {
            break;
        }
        if (my_next_perf->ecl_index > c)
        {
            continue;
        }
        assert(my_next_perf->ecl_index == c);
        const auto& connection = connections[c];
        if (connection.state() == Connection::State::OPEN) {
            Scalar radius = connection.rw();
            Scalar re = connection.re(); // area equivalent radius of the grid block
            Scalar perf_length = connection.connectionLength(); // the length of the well perforation
            const Scalar repR = std::sqrt(re * radius);
            perf_rep_radius_.push_back(repR);
            perf_length_.push_back(perf_length);
            bore_diameters_.push_back(2. * radius);
            num_active_connections++;
        }
        ++my_next_perf;
    }
    assert(my_next_perf == perf_data_->end());
    assert(num_active_connections == nperf);
}

template<class Scalar>
void WellInterfaceGeneric<Scalar>::
setWsolvent(const Scalar wsolvent)
{
    wsolvent_ = wsolvent;
}

template<class Scalar>
void WellInterfaceGeneric<Scalar>::
setDynamicThpLimit(const Scalar thp_limit)
{
    dynamic_thp_limit_ = thp_limit;
}

template<class Scalar>
std::optional<Scalar>
WellInterfaceGeneric<Scalar>::getDynamicThpLimit() const
{
    return dynamic_thp_limit_;
}

template<class Scalar>
void WellInterfaceGeneric<Scalar>::
updatePerforatedCell(std::vector<bool>& is_cell_perforated)
{
    for (int perf_idx = 0; perf_idx < number_of_perforations_; ++perf_idx) {
        is_cell_perforated[well_cells_[perf_idx]] = true;
    }
}

template<class Scalar>
bool WellInterfaceGeneric<Scalar>::
isVFPActive(DeferredLogger& deferred_logger) const
{
    // since the well_controls only handles the VFP number when THP constraint/target is there.
    // we need to get the table number through the parser, in case THP constraint/target is not there.
    // When THP control/limit is not active, if available VFP table is provided, we will still need to
    // update THP value. However, it will only used for output purpose.
    if (isProducer()) { // producer
        const int table_id = well_ecl_.vfp_table_number();
        if (table_id <= 0) {
            return false;
        } else {
            if (vfp_properties_->getProd()->hasTable(table_id)) {
                return true;
            } else {
                OPM_DEFLOG_THROW(std::runtime_error,
                                 fmt::format("VFPPROD table {} is specified "
                                             "for well {}, while we could not access it during simulation",
                                             table_id, name()),
                                 deferred_logger);
            }
        }

    } else { // injector
        const int table_id = well_ecl_.vfp_table_number();
        if (table_id <= 0) {
            return false;
        } else {
            if (vfp_properties_->getInj()->hasTable(table_id)) {
                return true;
            } else {
                OPM_DEFLOG_THROW(std::runtime_error,
                                 fmt::format("VFPINJ table {} is specified "
                                             "for well {}, while we could not access it during simulation",
                                             table_id, name()),
                                 deferred_logger);
            }
        }
    }
}

template<class Scalar>
bool WellInterfaceGeneric<Scalar>::isOperableAndSolvable() const
{
    return operability_status_.isOperableAndSolvable();
}

template<class Scalar>
bool WellInterfaceGeneric<Scalar>::useVfpExplicit() const
{
    const auto& wvfpexp = well_ecl_.getWVFPEXP();
    return (wvfpexp.explicit_lookup() || operability_status_.use_vfpexplicit);
}

template<class Scalar>
bool WellInterfaceGeneric<Scalar>::thpLimitViolatedButNotSwitched() const
{
    return operability_status_.thp_limit_violated_but_not_switched;
}

template<class Scalar>
Scalar WellInterfaceGeneric<Scalar>::
getALQ(const WellState<Scalar>& well_state) const
{
    // no alq for injectors.
    if (isInjector())
        return 0.0;

    return well_state.getALQ(name());
}

template<class Scalar>
void WellInterfaceGeneric<Scalar>::
reportWellSwitching(const SingleWellState<Scalar> &ws,
                    DeferredLogger& deferred_logger) const
{
    if (well_control_log_.empty())
        return;

    std::string from = well_control_log_[0];
    std::string to;
    if (isInjector()) {
        to = WellInjectorCMode2String(ws.injection_cmode);
    } else {
        to = WellProducerCMode2String(ws.production_cmode);
    }
    // only report the final switching
    if (from != to) {
        deferred_logger.info(fmt::format("    Well {} control mode changed from {} to {}",
                                         name(), from, to));
    }
}

template<class Scalar>
bool WellInterfaceGeneric<Scalar>::
isPressureControlled(const WellState<Scalar>& well_state) const
{
    const auto& ws = well_state.well(this->index_of_well_);
    if (this->isInjector()) {
        const Well::InjectorCMode& current = ws.injection_cmode;
        return current == Well::InjectorCMode::THP ||
               current == Well::InjectorCMode::BHP;
    } else {
        const Well::ProducerCMode& current = ws.production_cmode;
        return current == Well::ProducerCMode::THP ||
               current == Well::ProducerCMode::BHP;
    }
}

template<class Scalar>
bool WellInterfaceGeneric<Scalar>::
wellUnderZeroRateTargetIndividual(const SummaryState& summary_state,
                                  const WellState<Scalar>& well_state) const
{
    if (this->isProducer()) { // producers
        const auto prod_controls = this->well_ecl_.productionControls(summary_state);
        const auto prod_mode = well_state.well(this->indexOfWell()).production_cmode;
        return wellhelpers::rateControlWithZeroProdTarget(prod_controls, prod_mode);
    } else { // injectors
        const auto inj_controls = this->well_ecl_.injectionControls(summary_state);
        const auto inj_mode = well_state.well(this->indexOfWell()).injection_cmode;
        return wellhelpers::rateControlWithZeroInjTarget(inj_controls, inj_mode);
    }
}

template<class Scalar>
bool WellInterfaceGeneric<Scalar>::
wellUnderGroupControl(const SingleWellState<Scalar>& ws) const
{
    // Check if well is under group control
    const bool isGroupControlled = (this->isInjector() && ws.injection_cmode == Well::InjectorCMode::GRUP) ||
        (this->isProducer() && ws.production_cmode == Well::ProducerCMode::GRUP);
    return isGroupControlled;
}

template<class Scalar>
void WellInterfaceGeneric<Scalar>::resetWellOperability()
{
    this->operability_status_.resetOperability();
}

template<class Scalar>
Scalar WellInterfaceGeneric<Scalar>::wmicrobes_() const
{
    auto injectorType = this->well_ecl_.injectorType();

    if (injectorType == InjectorType::WATER) {
        WellMICPProperties microbes = this->well_ecl_.getMICPProperties();
        const Scalar microbial_injection_concentration = microbes.m_microbialConcentration;
        return microbial_injection_concentration;
    } else {
        // Not a water injection well => no microbes.
        return 0.0;
    }
}

template<class Scalar>
Scalar WellInterfaceGeneric<Scalar>::wfoam_() const
{
    auto injectorType = this->well_ecl_.injectorType();

    if (injectorType == InjectorType::GAS) {
        WellFoamProperties fprop = this->well_ecl_.getFoamProperties();
        return fprop.m_foamConcentration;
    } else {
        // Not a gas injection well => no foam.
        return 0.0;
    }
}

template<class Scalar>
Scalar WellInterfaceGeneric<Scalar>::wsalt_() const
{
    auto injectorType = this->well_ecl_.injectorType();

    if (injectorType == InjectorType::WATER) {
        WellBrineProperties fprop = this->well_ecl_.getBrineProperties();
        return fprop.m_saltConcentration;
    } else {
        // Not a water injection well => no salt (?).
        return 0.0;
    }
}

template<class Scalar>
Scalar WellInterfaceGeneric<Scalar>::woxygen_() const
{
    auto injectorType = this->well_ecl_.injectorType();

    if (injectorType == InjectorType::WATER) {
        WellMICPProperties oxygen = this->well_ecl_.getMICPProperties();
        const Scalar oxygen_injection_concentration = oxygen.m_oxygenConcentration;
        return oxygen_injection_concentration;
    } else {
        // Not a water injection well => no oxygen.
        return 0.0;
    }
}

template<class Scalar>
Scalar WellInterfaceGeneric<Scalar>::wpolymer_() const
{
    auto injectorType = this->well_ecl_.injectorType();

    if (injectorType == InjectorType::WATER) {
        WellPolymerProperties polymer = this->well_ecl_.getPolymerProperties();
        const Scalar polymer_injection_concentration = polymer.m_polymerConcentration;
        return polymer_injection_concentration;
    } else {
        // Not a water injection well => no polymer.
        return 0.0;
    }
}

template<class Scalar>
Scalar WellInterfaceGeneric<Scalar>::wurea_() const
{
    auto injectorType = this->well_ecl_.injectorType();

    if (injectorType == InjectorType::WATER) {
        WellMICPProperties urea = this->well_ecl_.getMICPProperties();
        const Scalar urea_injection_concentration = urea.m_ureaConcentration / 10.; //Dividing by scaling factor 10
        return urea_injection_concentration;
    } else {
        // Not a water injection well => no urea.
        return 0.0;
    }
}

template<class Scalar>
int WellInterfaceGeneric<Scalar>::polymerTable_() const
{
    return this->well_ecl_.getPolymerProperties().m_skprpolytable;
}

template<class Scalar>
int WellInterfaceGeneric<Scalar>::polymerWaterTable_() const
{
    return this->well_ecl_.getPolymerProperties().m_skprwattable;
}

template<class Scalar>
int WellInterfaceGeneric<Scalar>::polymerInjTable_() const
{
    return this->well_ecl_.getPolymerProperties().m_plymwinjtable;
}

template<class Scalar>
std::pair<bool,bool> WellInterfaceGeneric<Scalar>::
computeWellPotentials(std::vector<Scalar>& well_potentials,
                      const WellState<Scalar>& well_state)
{
    const int np = this->number_of_phases_;
    well_potentials.resize(np, 0.0);

    // Stopped wells have zero potential.
    if (this->wellIsStopped()) {
        return {false, false};
    }
    this->operability_status_.has_negative_potentials = false;

    // If the well is pressure controlled the potential equals the rate.
    bool thp_controlled_well = false;
    bool bhp_controlled_well = false;
    bool compute_potential = true;
    const auto& ws = well_state.well(this->index_of_well_);
    if (this->isInjector()) {
        const Well::InjectorCMode& current = ws.injection_cmode;
        if (current == Well::InjectorCMode::THP) {
            thp_controlled_well = true;
        }
        if (current == Well::InjectorCMode::BHP) {
            bhp_controlled_well = true;
        }
    } else {
        const Well::ProducerCMode& current = ws.production_cmode;
        if (current == Well::ProducerCMode::THP) {
            thp_controlled_well = true;
        }
        if (current == Well::ProducerCMode::BHP) {
            bhp_controlled_well = true;
        }
    }

    if (!this->changed_to_open_this_step_ &&
        (thp_controlled_well || bhp_controlled_well)) {
        Scalar total_rate = 0.0;
        const Scalar sign = this->isInjector() ? 1.0 : -1.0;
        for (int phase = 0; phase < np; ++phase){
            total_rate += sign * ws.surface_rates[phase];
        }
        // for pressure controlled wells the well rates are the potentials
        // if the rates are trivial we are most probably looking at the newly
        // opened well, and we therefore make the effort of computing the potentials anyway.
        if (total_rate > 0) {
            for (int phase = 0; phase < np; ++phase){
                well_potentials[phase] = sign * ws.surface_rates[phase];
            }
            compute_potential = false;
        }
    }

    return {compute_potential, bhp_controlled_well};
}

template<class Scalar>
void WellInterfaceGeneric<Scalar>::
checkNegativeWellPotentials(std::vector<Scalar>& well_potentials,
                            const bool checkOperability,
                            DeferredLogger& deferred_logger)
{
    const Scalar sign = this->isInjector() ? 1.0 : -1.0;
    Scalar total_potential = 0.0;
    for (int phase = 0; phase < this->number_of_phases_; ++phase) {
        well_potentials[phase] *= sign;
        total_potential += well_potentials[phase];
    }
    if (total_potential < 0.0 && checkOperability) {
        // wells with negative potentials are not operable
        this->operability_status_.has_negative_potentials = true;
        const std::string msg = std::string("well ") + this->name() +
                                ": has negative potentials and is not operable";
        deferred_logger.warning("NEGATIVE_POTENTIALS_INOPERABLE", msg);
    }
}

template<class Scalar>
void WellInterfaceGeneric<Scalar>::
prepareForPotentialCalculations(const SummaryState& summary_state,
                                WellState<Scalar>& well_state,
                                Well::InjectionControls& inj_controls,
                                Well::ProductionControls& prod_controls) const
{
    const bool has_thp = this->wellHasTHPConstraints(summary_state);
    auto& ws = well_state.well(this->index_of_well_);
    // Modify control (only pressure constraints) and set new target if needed.
    // Also set value for current target in state
    if (this->isInjector()) {
        inj_controls.clearControls();
        inj_controls.addControl(Well::InjectorCMode::BHP);
        if (has_thp){
            inj_controls.addControl(Well::InjectorCMode::THP);
        }
        if (!(ws.injection_cmode == Well::InjectorCMode::BHP)){
            if (has_thp){
                ws.injection_cmode = Well::InjectorCMode::THP;
                ws.thp = this->getTHPConstraint(summary_state);
            } else {
                ws.injection_cmode = Well::InjectorCMode::BHP;
                ws.bhp = inj_controls.bhp_limit;
            }
        } 
    } else { // producer
        prod_controls.clearControls();
        prod_controls.addControl(Well::ProducerCMode::BHP);
        if (has_thp){
            prod_controls.addControl(Well::ProducerCMode::THP);
        }
        if (!(ws.production_cmode == Well::ProducerCMode::BHP)){
            if (has_thp){
                ws.production_cmode = Well::ProducerCMode::THP;
                ws.thp = this->getTHPConstraint(summary_state);
            } else {
                ws.production_cmode = Well::ProducerCMode::BHP;
                ws.bhp = prod_controls.bhp_limit;
            }
        } 
    }    
}

template class WellInterfaceGeneric<double>;

#if FLOW_INSTANTIATE_FLOAT
template class WellInterfaceGeneric<float>;
#endif

} // namespace Opm
