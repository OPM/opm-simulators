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

namespace Opm
{

WellInterfaceGeneric::WellInterfaceGeneric(const Well& well,
                                           const ParallelWellInfo& pw_info,
                                           const int time_step,
                                           const int pvtRegionIdx,
                                           const int num_components,
                                           const int num_phases,
                                           const int index_of_well,
                                           const std::vector<PerforationData>& perf_data)
      : well_ecl_(well)
      , parallel_well_info_(pw_info)
      , current_step_(time_step)
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
void WellInterfaceGeneric::adaptRatesForVFP(std::vector<double>& rates) const
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

const std::vector<PerforationData>& WellInterfaceGeneric::perforationData() const
{
    return *perf_data_;
}

const std::string& WellInterfaceGeneric::name() const
{
    return well_ecl_.name();
}

bool WellInterfaceGeneric::isInjector() const
{
    return well_ecl_.isInjector();
}

bool WellInterfaceGeneric::isProducer() const
{
    return well_ecl_.isProducer();
}

int WellInterfaceGeneric::indexOfWell() const
{
    return index_of_well_;
}

bool WellInterfaceGeneric::getAllowCrossFlow() const
{
    return well_ecl_.getAllowCrossFlow();
}

const Well& WellInterfaceGeneric::wellEcl() const
{
    return well_ecl_;
}

Well& WellInterfaceGeneric::wellEcl()
{
    return well_ecl_;
}

const PhaseUsage& WellInterfaceGeneric::phaseUsage() const
{
    assert(phase_usage_ != nullptr);

    return *phase_usage_;
}

double WellInterfaceGeneric::wsolvent() const
{
    return wsolvent_;
}

double WellInterfaceGeneric::rsRvInj() const
{
    return well_ecl_.getInjectionProperties().rsRvInj;
}

void WellInterfaceGeneric::initInjMult(const std::vector<double>& max_inj_mult)
{
    // prev_inj_multiplier_ will stay unchanged during the time step
    // while inj_multiplier_ might be updated during the time step
    this->prev_inj_multiplier_ = max_inj_mult;
    // initializing the inj_multipler_ to be 1.0
    this->inj_multiplier_ = std::vector<double>(max_inj_mult.size(), 1.);
}

void WellInterfaceGeneric::updateInjMult(std::vector<double>& inj_multipliers, DeferredLogger& deferred_logger) const
{
    if (inj_multipliers.size() != this->inj_multiplier_.size()) {
        OPM_DEFLOG_THROW(std::runtime_error,
                         fmt::format("We do not support changing connection numbers during simulation with WINJMULT "
                                     "for well {}", name()),
                         deferred_logger);
    }

    inj_multipliers = this->inj_multiplier_;
}



double WellInterfaceGeneric::getInjMult(const int perf,
                                        const double bhp,
                                        const double perf_pres) const
{
    assert(!this->isProducer());

    double multiplier = 1.;

    const auto perf_ecl_index = this->perforationData()[perf].ecl_index;
    const bool is_wrev = this->well_ecl_.getInjMultMode() == Well::InjMultMode::WREV;

    const bool active_injmult = (is_wrev && this->well_ecl_.aciveWellInjMult()) ||
                              this->well_ecl_.getConnections()[perf_ecl_index].activeInjMult();

    if (active_injmult) {
        const auto& injmult= is_wrev ? this->well_ecl_.getWellInjMult()
                                                  : this->well_ecl_.getConnections()[perf_ecl_index].injmult();
        const double pres = is_wrev ? bhp : perf_pres;

        const auto frac_press = injmult.fracture_pressure;
        const auto gradient = injmult.multiplier_gradient;
        if (pres > frac_press) {
            multiplier = 1. + (pres - frac_press) * gradient;
        }
    }

    // for CIRR mode, if there is no active WINJMULT setup, we will use the previous injection multiplier,
    // to mimic keeping the existing fracturing open
    if (this->well_ecl_.getInjMultMode() == Well::InjMultMode::CIRR) {
        multiplier = std::max(multiplier, this->prev_inj_multiplier_[perf_ecl_index]);
    }

    this->inj_multiplier_[perf_ecl_index] = multiplier;
    return multiplier;
}




bool WellInterfaceGeneric::wellHasTHPConstraints(const SummaryState& summaryState) const
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

void WellInterfaceGeneric::updateWellTestState(const SingleWellState& ws,
                                               const double& simulationTime,
                                               const bool& writeMessageToOPMLog,
                                               WellTestState& wellTestState,
                                               DeferredLogger& deferred_logger) const
{
    // updating well test state based on Economic limits for operable wells
    if (this->isOperableAndSolvable()) {
        WellTest(*this).updateWellTestStateEconomic(ws, simulationTime, writeMessageToOPMLog, wellTestState, deferred_logger);
    } else {
        // updating well test state based on physical (THP/BHP) limits.
        WellTest(*this).updateWellTestStatePhysical(simulationTime, writeMessageToOPMLog, wellTestState, deferred_logger);
    }

    // TODO: well can be shut/closed due to other reasons
}

double WellInterfaceGeneric::getTHPConstraint(const SummaryState& summaryState) const
{
    if (dynamic_thp_limit_) {
        return *dynamic_thp_limit_;
    }

    return WellBhpThpCalculator(*this).getTHPConstraint(summaryState);
}

bool WellInterfaceGeneric::underPredictionMode() const
{
    return well_ecl_.predictionMode();
}

void WellInterfaceGeneric::initCompletions()
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

void WellInterfaceGeneric::closeCompletions(const WellTestState& wellTestState)
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

void WellInterfaceGeneric::setVFPProperties(const VFPProperties* vfp_properties_arg)
{
    vfp_properties_ = vfp_properties_arg;
}

void WellInterfaceGeneric::setPrevSurfaceRates(WellState& well_state,
                                               const WellState& prev_well_state) const
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
                [](double rate) {
                    return rate == 0.; // TODO: should we use a threshhold for comparison?
                } );

        if (zero_rates) {
            ws.prev_surface_rates = ws_prev.prev_surface_rates;
        } else {
            ws.prev_surface_rates = ws.surface_rates;
        }
    }

void WellInterfaceGeneric::setGuideRate(const GuideRate* guide_rate_arg)
{
    guide_rate_ = guide_rate_arg;
}

void WellInterfaceGeneric::setWellEfficiencyFactor(const double efficiency_factor)
{
    well_efficiency_factor_ = efficiency_factor;
}

void WellInterfaceGeneric::setRepRadiusPerfLength()
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
            double radius = connection.rw();
            double re = connection.re(); // area equivalent radius of the grid block
            double perf_length = connection.connectionLength(); // the length of the well perforation
            const double repR = std::sqrt(re * radius);
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

void WellInterfaceGeneric::setWsolvent(const double wsolvent)
{
    wsolvent_ = wsolvent;
}

void WellInterfaceGeneric::setDynamicThpLimit(const double thp_limit)
{
    dynamic_thp_limit_ = thp_limit;
}

std::optional<double> WellInterfaceGeneric::getDynamicThpLimit() const
{
    return dynamic_thp_limit_;
}

void WellInterfaceGeneric::updatePerforatedCell(std::vector<bool>& is_cell_perforated)
{

    for (int perf_idx = 0; perf_idx<number_of_perforations_; ++perf_idx) {
        is_cell_perforated[well_cells_[perf_idx]] = true;
    }
}

bool WellInterfaceGeneric::isVFPActive(DeferredLogger& deferred_logger) const
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

bool WellInterfaceGeneric::isOperableAndSolvable() const
{
    return operability_status_.isOperableAndSolvable();
}

bool WellInterfaceGeneric::useVfpExplicit() const
{
    const auto& wvfpexp = well_ecl_.getWVFPEXP();
    return (wvfpexp.explicit_lookup() || operability_status_.use_vfpexplicit);
}

bool WellInterfaceGeneric::thpLimitViolatedButNotSwitched() const
{
    return operability_status_.thp_limit_violated_but_not_switched;
}

double WellInterfaceGeneric::getALQ(const WellState& well_state) const
{
    // no alq for injectors.
    if (isInjector())
        return 0.0;

    return well_state.getALQ(name());
}

void WellInterfaceGeneric::reportWellSwitching(const SingleWellState& ws, DeferredLogger& deferred_logger) const
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

bool WellInterfaceGeneric::isPressureControlled(const WellState& well_state) const
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



bool WellInterfaceGeneric::wellUnderZeroRateTarget(const SummaryState& summary_state,
                                                   const WellState& well_state) const
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

bool WellInterfaceGeneric::stopppedOrZeroRateTarget(const SummaryState& summary_state,
                                                    const WellState& well_state) const
{
    return (this->wellIsStopped() || this->wellUnderZeroRateTarget(summary_state, well_state));

}

void WellInterfaceGeneric::resetWellOperability()
{
    this->operability_status_.resetOperability();
}

double WellInterfaceGeneric::wmicrobes_() const
{
    auto injectorType = this->well_ecl_.injectorType();

    if (injectorType == InjectorType::WATER) {
        WellMICPProperties microbes = this->well_ecl_.getMICPProperties();
        const double microbial_injection_concentration = microbes.m_microbialConcentration;
        return microbial_injection_concentration;
    } else {
        // Not a water injection well => no microbes.
        return 0.0;
    }
}

double WellInterfaceGeneric::wfoam_() const
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

double WellInterfaceGeneric::wsalt_() const
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

double WellInterfaceGeneric::woxygen_() const
{
    auto injectorType = this->well_ecl_.injectorType();

    if (injectorType == InjectorType::WATER) {
        WellMICPProperties oxygen = this->well_ecl_.getMICPProperties();
        const double oxygen_injection_concentration = oxygen.m_oxygenConcentration;
        return oxygen_injection_concentration;
    } else {
        // Not a water injection well => no oxygen.
        return 0.0;
    }
}

double WellInterfaceGeneric::wpolymer_() const
{
    auto injectorType = this->well_ecl_.injectorType();

    if (injectorType == InjectorType::WATER) {
        WellPolymerProperties polymer = this->well_ecl_.getPolymerProperties();
        const double polymer_injection_concentration = polymer.m_polymerConcentration;
        return polymer_injection_concentration;
    } else {
        // Not a water injection well => no polymer.
        return 0.0;
    }
}

double WellInterfaceGeneric::wurea_() const
{
    auto injectorType = this->well_ecl_.injectorType();

    if (injectorType == InjectorType::WATER) {
        WellMICPProperties urea = this->well_ecl_.getMICPProperties();
        const double urea_injection_concentration = urea.m_ureaConcentration / 10.; //Dividing by scaling factor 10
        return urea_injection_concentration;
    } else {
        // Not a water injection well => no urea.
        return 0.0;
    }
}

int WellInterfaceGeneric::polymerTable_() const
{
    return this->well_ecl_.getPolymerProperties().m_skprpolytable;
}

int WellInterfaceGeneric::polymerWaterTable_() const
{
    return this->well_ecl_.getPolymerProperties().m_skprwattable;
}

int WellInterfaceGeneric::polymerInjTable_() const
{
    return this->well_ecl_.getPolymerProperties().m_plymwinjtable;
}

std::pair<bool,bool> WellInterfaceGeneric::
computeWellPotentials(std::vector<double>& well_potentials,
                      const WellState& well_state)
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
        double total_rate = 0.0;
        const double sign = this->isInjector() ? 1.0 : -1.0;
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

void WellInterfaceGeneric::
checkNegativeWellPotentials(std::vector<double>& well_potentials,
                            const bool checkOperability,
                            DeferredLogger& deferred_logger)
{
    const double sign = this->isInjector() ? 1.0 : -1.0;
    double total_potential = 0.0;
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

void WellInterfaceGeneric::
prepareForPotentialCalculations(const SummaryState& summary_state,
                                WellState& well_state, 
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

} // namespace Opm
