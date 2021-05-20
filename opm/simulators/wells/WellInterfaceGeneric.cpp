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

#include <opm/parser/eclipse/EclipseState/Schedule/Well/WellTestState.hpp>

#include <opm/simulators/utils/DeferredLoggingErrorHelpers.hpp>
#include <opm/simulators/wells/PerforationData.hpp>
#include <opm/simulators/wells/ParallelWellInfo.hpp>
#include <opm/simulators/wells/VFPProperties.hpp>
#include <opm/simulators/wells/WellState.hpp>

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
                                           const int first_perf_index,
                                           const std::vector<PerforationData>& perf_data)
      : well_ecl_(well)
      , parallel_well_info_(pw_info)
      , current_step_(time_step)
      , pvtRegionIdx_(pvtRegionIdx)
      , num_components_(num_components)
      , number_of_phases_(num_phases)
      , index_of_well_(index_of_well)
      , first_perf_(first_perf_index)
      , perf_data_(&perf_data)
      , ipr_a_(number_of_phases_)
      , ipr_b_(number_of_phases_)
{
    assert(well.name()==pw_info.name());
    assert(std::is_sorted(perf_data.begin(), perf_data.end(),
                          [](const auto& perf1, const auto& perf2){
        return perf1.ecl_index < perf2.ecl_index;
    }));
    if (time_step < 0) {
        OPM_THROW(std::invalid_argument, "Negtive time step is used to construct WellInterface");
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

const PhaseUsage& WellInterfaceGeneric::phaseUsage() const
{
    assert(phase_usage_ != nullptr);

    return *phase_usage_;
}

double WellInterfaceGeneric::wsolvent() const
{
    return wsolvent_;
}

bool WellInterfaceGeneric::wellHasTHPConstraints(const SummaryState& summaryState) const
{
    if (dynamic_thp_limit_) {
        return true;
    }

    if (well_ecl_.isInjector()) {
        const auto controls = well_ecl_.injectionControls(summaryState);
        if (controls.hasControl(Well::InjectorCMode::THP))
            return true;
    }

    if (well_ecl_.isProducer( )) {
        const auto controls = well_ecl_.productionControls(summaryState);
        if (controls.hasControl(Well::ProducerCMode::THP))
            return true;
    }

    return false;

}

double WellInterfaceGeneric::mostStrictBhpFromBhpLimits(const SummaryState& summaryState) const
{
    if (well_ecl_.isInjector()) {
        const auto& controls = well_ecl_.injectionControls(summaryState);
        return controls.bhp_limit;
    }

    if (well_ecl_.isProducer( )) {
        const auto& controls = well_ecl_.productionControls(summaryState);
        return controls.bhp_limit;
    }

    return 0.0;
}

double WellInterfaceGeneric::getTHPConstraint(const SummaryState& summaryState) const
{
    if (dynamic_thp_limit_) {
        return *dynamic_thp_limit_;
    }
    if (well_ecl_.isInjector()) {
        const auto& controls = well_ecl_.injectionControls(summaryState);
        return controls.thp_limit;
    }

    if (well_ecl_.isProducer( )) {
        const auto& controls = well_ecl_.productionControls(summaryState);
        return controls.thp_limit;
    }

    return 0.0;
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

void WellInterfaceGeneric::closeCompletions(WellTestState& wellTestState)
{
    const auto& connections = well_ecl_.getConnections();
    int perfIdx = 0;
    for (const auto& connection : connections) {
        if (connection.state() == Connection::State::OPEN) {
            if (wellTestState.hasCompletion(name(), connection.complnum())) {
                well_index_[perfIdx] = 0.0;
            }
            perfIdx++;
        }
    }
}

void WellInterfaceGeneric::setVFPProperties(const VFPProperties* vfp_properties_arg)
{
    vfp_properties_ = vfp_properties_arg;
}

void WellInterfaceGeneric::setGuideRate(const GuideRate* guide_rate_arg)
{
    guide_rate_ = guide_rate_arg;
}

void WellInterfaceGeneric::setWellEfficiencyFactor(const double efficiency_factor)
{
    well_efficiency_factor_ = efficiency_factor;
}

void WellInterfaceGeneric::setRepRadiusPerfLength(const std::vector<int>& cartesian_to_compressed)
{
    const int nperf = number_of_perforations_;

    perf_rep_radius_.clear();
    perf_length_.clear();
    bore_diameters_.clear();

    perf_rep_radius_.reserve(nperf);
    perf_length_.reserve(nperf);
    bore_diameters_.reserve(nperf);

    // COMPDAT handling
    const auto& connectionSet = well_ecl_.getConnections();
    CheckDistributedWellConnections checker(well_ecl_, parallel_well_info_);
    for (size_t c=0; c<connectionSet.size(); c++) {
        const auto& connection = connectionSet.get(c);
        const int cell =
                cartesian_to_compressed[connection.global_index()];
        if (connection.state() != Connection::State::OPEN || cell >= 0)
        {
            checker.connectionFound(c);
        }
        if (connection.state() == Connection::State::OPEN) {

            if (cell >= 0) {
                double radius = connection.rw();
                double re = connection.re(); // area equivalent radius of the grid block
                double perf_length = connection.connectionLength(); // the length of the well perforation
                const double repR = std::sqrt(re * radius);
                perf_rep_radius_.push_back(repR);
                perf_length_.push_back(perf_length);
                bore_diameters_.push_back(2. * radius);
            }
        }
    }
    checker.checkAllConnectionsFound();
}

void WellInterfaceGeneric::setWsolvent(const double wsolvent)
{
    wsolvent_ = wsolvent;
}

void WellInterfaceGeneric::setDynamicThpLimit(const double thp_limit)
{
    dynamic_thp_limit_ = thp_limit;
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
                OPM_DEFLOG_THROW(std::runtime_error, "VFPPROD table " << std::to_string(table_id) << " is specfied,"
                              << " for well " << name() << ", while we could not access it during simulation", deferred_logger);
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
                OPM_DEFLOG_THROW(std::runtime_error, "VFPINJ table " << std::to_string(table_id) << " is specfied,"
                              << " for well " << name() << ", while we could not access it during simulation", deferred_logger);
            }
        }
    }
}

void WellInterfaceGeneric::updateWellTestStatePhysical(const WellState& /* well_state */,
                                                       const double simulation_time,
                                                       const bool write_message_to_opmlog,
                                                       WellTestState& well_test_state,
                                                       DeferredLogger& deferred_logger) const
{
    if (!isOperable()) {
        if (well_test_state.hasWellClosed(name(), WellTestConfig::Reason::ECONOMIC) ||
            well_test_state.hasWellClosed(name(), WellTestConfig::Reason::PHYSICAL) ) {
            // Already closed, do nothing.
        } else {
            well_test_state.closeWell(name(), WellTestConfig::Reason::PHYSICAL, simulation_time);
            if (write_message_to_opmlog) {
                const std::string action = well_ecl_.getAutomaticShutIn() ? "shut" : "stopped";
                const std::string msg = "Well " + name()
                    + " will be " + action + " as it can not operate under current reservoir conditions.";
                deferred_logger.info(msg);
            }
        }
    }
}

bool WellInterfaceGeneric::isOperable() const
{
    return operability_status_.isOperable();
}

double WellInterfaceGeneric::getALQ(const WellState& well_state) const
{
    return well_state.getALQ(name());
}


} // namespace Opm
