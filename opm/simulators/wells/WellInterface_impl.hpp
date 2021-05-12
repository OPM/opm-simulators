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

#include <opm/parser/eclipse/EclipseState/Schedule/ScheduleTypes.hpp>
#include <opm/simulators/utils/DeferredLoggingErrorHelpers.hpp>
#include <opm/simulators/wells/TargetCalculator.hpp>

namespace Opm
{


    template<typename TypeTag>
    WellInterface<TypeTag>::
    WellInterface(const Well& well,
                  const ParallelWellInfo& pw_info,
                  const int time_step,
                  const ModelParameters& param,
                  const RateConverterType& rate_converter,
                  const int pvtRegionIdx,
                  const int num_components,
                  const int num_phases,
                  const int index_of_well,
                  const int first_perf_index,
                  const std::vector<PerforationData>& perf_data)
      : well_ecl_(well)
      , parallel_well_info_(pw_info)
      , current_step_(time_step)
      , param_(param)
      , rateConverter_(rate_converter)
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

        connectionRates_.resize(number_of_perforations_);

        this->wellStatus_ = Well::Status::OPEN;
        if (well.getStatus() == Well::Status::STOP) {
            this->wellStatus_ = Well::Status::STOP;
        }

        wsolvent_ = 0.0;

        if constexpr (has_solvent || has_zFraction) {
            if (well.isInjector()) {
                auto injectorType = well_ecl_.injectorType();
                if (injectorType == InjectorType::GAS) {
                    wsolvent_ = well_ecl_.getSolventFraction();
                }
            }
        }
    }

    template<typename TypeTag>
    void
    WellInterface<TypeTag>::
    updatePerforatedCell(std::vector<bool>& is_cell_perforated)
    {

        for (int perf_idx = 0; perf_idx<number_of_perforations_; ++perf_idx) {
            is_cell_perforated[well_cells_[perf_idx]] = true;
        }
    }



    template<typename TypeTag>
    void
    WellInterface<TypeTag>::
    init(const PhaseUsage* phase_usage_arg,
         const std::vector<double>& /* depth_arg */,
         const double gravity_arg,
         const int /* num_cells */,
         const std::vector< Scalar >& B_avg)
    {
        phase_usage_ = phase_usage_arg;
        gravity_ = gravity_arg;
        B_avg_ = B_avg;
    }





    template<typename TypeTag>
    void
    WellInterface<TypeTag>::
    initCompletions()
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






    template<typename TypeTag>
    void
    WellInterface<TypeTag>::
    setVFPProperties(const VFPProperties* vfp_properties_arg)
    {
        vfp_properties_ = vfp_properties_arg;
    }

    template<typename TypeTag>
    void
    WellInterface<TypeTag>::
    setGuideRate(const GuideRate* guide_rate_arg)
    {
        guide_rate_ = guide_rate_arg;
    }


    template<typename TypeTag>
    const std::string&
    WellInterface<TypeTag>::
    name() const
    {
        return well_ecl_.name();
    }





    template<typename TypeTag>
    bool
    WellInterface<TypeTag>::
    isInjector() const
    {
        return well_ecl_.isInjector();
    }




    template<typename TypeTag>
    bool
    WellInterface<TypeTag>::
    isProducer() const
    {
        return well_ecl_.isProducer();
    }




    template<typename TypeTag>
    int
    WellInterface<TypeTag>::
    indexOfWell() const
    {
        return index_of_well_;
    }






    template<typename TypeTag>
    bool
    WellInterface<TypeTag>::
    getAllowCrossFlow() const
    {
        return well_ecl_.getAllowCrossFlow();
    }




    template<typename TypeTag>
    void
    WellInterface<TypeTag>::
    setWellEfficiencyFactor(const double efficiency_factor)
    {
        well_efficiency_factor_ = efficiency_factor;
    }



    template<typename TypeTag>
    const Well&
    WellInterface<TypeTag>::
    wellEcl() const
    {
      return well_ecl_;
    }



    template<typename TypeTag>
    const PhaseUsage&
    WellInterface<TypeTag>::
    phaseUsage() const
    {
        assert(phase_usage_ != nullptr);

        return *phase_usage_;
    }





    template<typename TypeTag>
    int
    WellInterface<TypeTag>::
    flowPhaseToEbosCompIdx( const int phaseIdx ) const
    {
        const auto& pu = phaseUsage();
        if (FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx) && pu.phase_pos[Water] == phaseIdx)
            return Indices::canonicalToActiveComponentIndex(FluidSystem::waterCompIdx);
        if (FluidSystem::phaseIsActive(FluidSystem::oilPhaseIdx) && pu.phase_pos[Oil] == phaseIdx)
            return Indices::canonicalToActiveComponentIndex(FluidSystem::oilCompIdx);
        if (FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx) && pu.phase_pos[Gas] == phaseIdx)
            return Indices::canonicalToActiveComponentIndex(FluidSystem::gasCompIdx);

        // for other phases return the index
        return phaseIdx;
    }

    template<typename TypeTag>
    int
    WellInterface<TypeTag>::
    flowPhaseToEbosPhaseIdx( const int phaseIdx ) const
    {
        const auto& pu = phaseUsage();
        if (FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx) && pu.phase_pos[Water] == phaseIdx)
            return FluidSystem::waterPhaseIdx;
        if (FluidSystem::phaseIsActive(FluidSystem::oilPhaseIdx) && pu.phase_pos[Oil] == phaseIdx)
            return FluidSystem::oilPhaseIdx;
        if (FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx) && pu.phase_pos[Gas] == phaseIdx)
            return FluidSystem::gasPhaseIdx;

        // for other phases return the index
        return phaseIdx;
    }

    template<typename TypeTag>
    int
    WellInterface<TypeTag>::
    ebosCompIdxToFlowCompIdx( const unsigned compIdx ) const
    {
        const auto& pu = phaseUsage();
        if (FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx) && Indices::canonicalToActiveComponentIndex(FluidSystem::waterCompIdx) == compIdx)
            return pu.phase_pos[Water];
        if (FluidSystem::phaseIsActive(FluidSystem::oilPhaseIdx) && Indices::canonicalToActiveComponentIndex(FluidSystem::oilCompIdx) == compIdx)
            return pu.phase_pos[Oil];
        if (FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx) && Indices::canonicalToActiveComponentIndex(FluidSystem::gasCompIdx) == compIdx)
            return pu.phase_pos[Gas];

        // for other phases return the index
        return compIdx;
    }




    template<typename TypeTag>
    double
    WellInterface<TypeTag>::
    wsolvent() const
    {
        return wsolvent_;
    }



    template<typename TypeTag>
    void
    WellInterface<TypeTag>::
    setWsolvent(const double wsolvent)
    {
       wsolvent_ = wsolvent;
    }





    template<typename TypeTag>
    void
    WellInterface<TypeTag>::
    setDynamicThpLimit(const double thp_limit)
    {
       dynamic_thp_limit_ = thp_limit;
    }





    template<typename TypeTag>
    double
    WellInterface<TypeTag>::
    wpolymer() const
    {
        if constexpr (has_polymer) {
            auto injectorType = well_ecl_.injectorType();

            if (injectorType == InjectorType::WATER) {
                WellPolymerProperties polymer = well_ecl_.getPolymerProperties();
                const double polymer_injection_concentration = polymer.m_polymerConcentration;
                return polymer_injection_concentration;
            } else {
                // Not a water injection well => no polymer.
                return 0.0;
            }
        }

        return 0.0;
    }





    template<typename TypeTag>
    double
    WellInterface<TypeTag>::
    wfoam() const
    {
        if constexpr (has_foam) {
            auto injectorType = well_ecl_.injectorType();

            if (injectorType == InjectorType::GAS) {
                WellFoamProperties fprop = well_ecl_.getFoamProperties();
                return fprop.m_foamConcentration;
            } else {
                // Not a gas injection well => no foam.
                return 0.0;
            }
        }

        return 0.0;
    }



    template<typename TypeTag>
    double
    WellInterface<TypeTag>::
    wsalt() const
    {
        if constexpr (has_brine) {
            auto injectorType = well_ecl_.injectorType();

            if (injectorType == InjectorType::WATER) {
                WellBrineProperties fprop = well_ecl_.getBrineProperties();
                return fprop.m_saltConcentration;
            } else {
                // Not a water injection well => no salt (?).
                return 0.0;
            }
        }

        return 0.0;
    }




    template<typename TypeTag>
    bool
    WellInterface<TypeTag>::
    wellHasTHPConstraints(const SummaryState& summaryState) const
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

    template<typename TypeTag>
    double
    WellInterface<TypeTag>::
    mostStrictBhpFromBhpLimits(const SummaryState& summaryState) const
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

    template<typename TypeTag>
    double
    WellInterface<TypeTag>::
    getTHPConstraint(const SummaryState& summaryState) const
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






    template<typename TypeTag>
    bool
    WellInterface<TypeTag>::
    updateWellControl(const Simulator& ebos_simulator,
                      const IndividualOrGroup iog,
                      WellState& well_state,
                      const GroupState& group_state,
                      Opm::DeferredLogger& deferred_logger) /* const */
    {
        if (this->wellIsStopped()) {
            return false;
        }

        const auto& summaryState = ebos_simulator.vanguard().summaryState();
        const auto& schedule = ebos_simulator.vanguard().schedule();
        const auto& well = well_ecl_;
        std::string from;
        if (well.isInjector()) {
            from = Well::InjectorCMode2String(well_state.currentInjectionControl(index_of_well_));
        } else {
            from = Well::ProducerCMode2String(well_state.currentProductionControl(index_of_well_));
        }

        bool changed = false;
        if (iog == IndividualOrGroup::Individual) {
            changed = checkIndividualConstraints(well_state, summaryState);
        } else if (iog == IndividualOrGroup::Group) {
            changed = checkGroupConstraints(well_state, group_state, schedule, summaryState, deferred_logger);
        } else {
            assert(iog == IndividualOrGroup::Both);
            changed = checkConstraints(well_state, group_state, schedule, summaryState, deferred_logger);
        }

        auto cc = Dune::MPIHelper::getCollectiveCommunication();

        // checking whether control changed
        if (changed) {
            std::string to;
            if (well.isInjector()) {
                to = Well::InjectorCMode2String(well_state.currentInjectionControl(index_of_well_));
            } else {
                to = Well::ProducerCMode2String(well_state.currentProductionControl(index_of_well_));
            }
            std::ostringstream ss;
            ss << "    Switching control mode for well " << name()
               << " from " << from
               << " to " <<  to;
            if (cc.size() > 1) {
               ss << " on rank " << cc.rank();
            }
            deferred_logger.info(ss.str());
            updateWellStateWithTarget(ebos_simulator, well_state, deferred_logger);
            updatePrimaryVariables(well_state, deferred_logger);
        }

        return changed;
    }





    template<typename TypeTag>
    bool
    WellInterface<TypeTag>::
    underPredictionMode() const
    {
        return well_ecl_.predictionMode();
    }





    template<typename TypeTag>
    bool
    WellInterface<TypeTag>::
    checkRateEconLimits(const WellEconProductionLimits& econ_production_limits,
                        const std::vector<double>& well_rates,
                        Opm::DeferredLogger& deferred_logger) const
    {
        const Opm::PhaseUsage& pu = phaseUsage();
        const int np = number_of_phases_;

        if (econ_production_limits.onMinOilRate()) {
            assert(FluidSystem::phaseIsActive(FluidSystem::oilPhaseIdx));
            const double oil_rate = well_rates[index_of_well_ * np + pu.phase_pos[ Oil ] ];
            const double min_oil_rate = econ_production_limits.minOilRate();
            if (std::abs(oil_rate) < min_oil_rate) {
                return true;
            }
        }

        if (econ_production_limits.onMinGasRate() ) {
            assert(FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx));
            const double gas_rate = well_rates[index_of_well_ * np + pu.phase_pos[ Gas ] ];
            const double min_gas_rate = econ_production_limits.minGasRate();
            if (std::abs(gas_rate) < min_gas_rate) {
                return true;
            }
        }

        if (econ_production_limits.onMinLiquidRate() ) {
            assert(FluidSystem::phaseIsActive(FluidSystem::oilPhaseIdx));
            assert(FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx));
            const double oil_rate = well_rates[index_of_well_ * np + pu.phase_pos[ Oil ] ];
            const double water_rate = well_rates[index_of_well_ * np + pu.phase_pos[ Water ] ];
            const double liquid_rate = oil_rate + water_rate;
            const double min_liquid_rate = econ_production_limits.minLiquidRate();
            if (std::abs(liquid_rate) < min_liquid_rate) {
                return true;
            }
        }

        if (econ_production_limits.onMinReservoirFluidRate()) {
            deferred_logger.warning("NOT_SUPPORTING_MIN_RESERVOIR_FLUID_RATE", "Minimum reservoir fluid production rate limit is not supported yet");
        }

        return false;
    }



    template<typename TypeTag>
    template<class ValueType>
    ValueType
    WellInterface<TypeTag>::
    calculateBhpFromThp(const WellState &well_state,
                        const std::vector<ValueType>& rates,
                        const Well& well,
                        const SummaryState& summaryState,
                        Opm::DeferredLogger& deferred_logger) const
    {
        // TODO: when well is under THP control, the BHP is dependent on the rates,
        // the well rates is also dependent on the BHP, so it might need to do some iteration.
        // However, when group control is involved, change of the rates might impacts other wells
        // so iterations on a higher level will be required. Some investigation might be needed when
        // we face problems under THP control.

        assert(int(rates.size()) == 3); // the vfp related only supports three phases now.

        const ValueType aqua = rates[Water];
        const ValueType liquid = rates[Oil];
        const ValueType vapour = rates[Gas];

        // pick the reference density
        // typically the reference in the top layer
        const double rho = getRefDensity();

        if (this->isInjector() )
        {
            const auto& controls = well.injectionControls(summaryState);
            const double vfp_ref_depth = vfp_properties_->getInj()->getTable(controls.vfp_table_number).getDatumDepth();
            const double dp = wellhelpers::computeHydrostaticCorrection(ref_depth_, vfp_ref_depth, rho, gravity_);
            return vfp_properties_->getInj()->bhp(controls.vfp_table_number, aqua, liquid, vapour, this->getTHPConstraint(summaryState)) - dp;
         }
         else if (this->isProducer()) {
             const auto& controls = well.productionControls(summaryState);
             const double vfp_ref_depth = vfp_properties_->getProd()->getTable(controls.vfp_table_number).getDatumDepth();
             const double dp = wellhelpers::computeHydrostaticCorrection(ref_depth_, vfp_ref_depth, rho, gravity_);
             return vfp_properties_->getProd()->bhp(controls.vfp_table_number, aqua, liquid, vapour, this->getTHPConstraint(summaryState), getALQ(well_state)) - dp;
         }
         else {
             OPM_DEFLOG_THROW(std::logic_error, "Expected INJECTOR or PRODUCER for well " + name(), deferred_logger);
         }



    }


    template<typename TypeTag>
    double
    WellInterface<TypeTag>::
    getALQ(const WellState &well_state) const
    {
        return well_state.getALQ(name());
    }



    template<typename TypeTag>
    void
    WellInterface<TypeTag>::
    checkMaxWaterCutLimit(const WellEconProductionLimits& econ_production_limits,
                          const WellState& well_state,
                          RatioLimitCheckReport& report) const
    {

        assert(FluidSystem::phaseIsActive(FluidSystem::oilPhaseIdx));
        assert(FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx));

        // function to calculate water cut based on rates
        auto waterCut = [](const std::vector<double>& rates,
                           const PhaseUsage& pu) {

            const double oil_rate = rates[pu.phase_pos[Oil]];
            const double water_rate = rates[pu.phase_pos[Water]];

            // both rate should be in the same direction
            assert(oil_rate * water_rate >= 0.);

            const double liquid_rate = oil_rate + water_rate;
            if (liquid_rate != 0.) {
                return (water_rate / liquid_rate);
            } else {
                return 0.;
            }
        };

        const double max_water_cut_limit = econ_production_limits.maxWaterCut();
        assert(max_water_cut_limit > 0.);

        const bool watercut_limit_violated = checkMaxRatioLimitWell(well_state, max_water_cut_limit, waterCut);

        if (watercut_limit_violated) {
            report.ratio_limit_violated = true;
            checkMaxRatioLimitCompletions(well_state, max_water_cut_limit, waterCut, report);
        }
    }





    template<typename TypeTag>
    void
    WellInterface<TypeTag>::
    checkMaxGORLimit(const WellEconProductionLimits& econ_production_limits,
                     const WellState& well_state,
                     RatioLimitCheckReport& report) const
    {

        assert(FluidSystem::phaseIsActive(FluidSystem::oilPhaseIdx));
        assert(FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx));

        // function to calculate gor based on rates
        auto gor = [](const std::vector<double>& rates,
                      const PhaseUsage& pu) {

            const double oil_rate = rates[pu.phase_pos[Oil]];
            const double gas_rate = rates[pu.phase_pos[Gas]];

            // both rate should be in the same direction
            assert(oil_rate * gas_rate >= 0.);

            double gas_oil_ratio = 0.;

            if (oil_rate != 0.) {
                gas_oil_ratio = gas_rate / oil_rate;
            } else {
                if (gas_rate != 0.) {
                    gas_oil_ratio = 1.e100; // big value to mark it as violated
                } else {
                    gas_oil_ratio = 0.0;
                }
            }

            return gas_oil_ratio;
        };

        const double max_gor_limit = econ_production_limits.maxGasOilRatio();
        assert(max_gor_limit > 0.);

        const bool gor_limit_violated = checkMaxRatioLimitWell(well_state, max_gor_limit, gor);

        if (gor_limit_violated) {
            report.ratio_limit_violated = true;
            checkMaxRatioLimitCompletions(well_state, max_gor_limit, gor, report);
        }
    }





    template<typename TypeTag>
    void
    WellInterface<TypeTag>::
    checkMaxWGRLimit(const WellEconProductionLimits& econ_production_limits,
                     const WellState& well_state,
                     RatioLimitCheckReport& report) const
    {

        assert(FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx));
        assert(FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx));

        // function to calculate wgr based on rates
        auto wgr = [](const std::vector<double>& rates,
                      const PhaseUsage& pu) {

            const double water_rate = rates[pu.phase_pos[Water]];
            const double gas_rate = rates[pu.phase_pos[Gas]];

            // both rate should be in the same direction
            assert(water_rate * gas_rate >= 0.);

            double water_gas_ratio = 0.;

            if (gas_rate != 0.) {
                water_gas_ratio = water_rate / gas_rate;
            } else {
                if (water_rate != 0.) {
                    water_gas_ratio = 1.e100; // big value to mark it as violated
                } else {
                    water_gas_ratio = 0.0;
                }
            }

            return water_gas_ratio;
        };

        const double max_wgr_limit = econ_production_limits.maxWaterGasRatio();
        assert(max_wgr_limit > 0.);

        const bool wgr_limit_violated = checkMaxRatioLimitWell(well_state, max_wgr_limit, wgr);

        if (wgr_limit_violated) {
            report.ratio_limit_violated = true;
            checkMaxRatioLimitCompletions(well_state, max_wgr_limit, wgr, report);
        }
    }





    template<typename TypeTag>
    void
    WellInterface<TypeTag>::
    checkRatioEconLimits(const WellEconProductionLimits& econ_production_limits,
                         const WellState& well_state,
                         RatioLimitCheckReport& report,
                         Opm::DeferredLogger& deferred_logger) const
    {
        // TODO: not sure how to define the worst-offending completion when more than one
        //       ratio related limit is violated.
        //       The defintion used here is that we define the violation extent based on the
        //       ratio between the value and the corresponding limit.
        //       For each violated limit, we decide the worst-offending completion separately.
        //       Among the worst-offending completions, we use the one has the biggest violation
        //       extent.

        if (econ_production_limits.onMaxWaterCut()) {
            checkMaxWaterCutLimit(econ_production_limits, well_state, report);
        }

        if (econ_production_limits.onMaxGasOilRatio()) {
            checkMaxGORLimit(econ_production_limits, well_state, report);
        }

        if (econ_production_limits.onMaxWaterGasRatio()) {
            checkMaxWGRLimit(econ_production_limits, well_state, report);
        }

        if (econ_production_limits.onMaxGasLiquidRatio()) {
            deferred_logger.warning("NOT_SUPPORTING_MAX_GLR", "the support for max Gas-Liquid ratio is not implemented yet!");
        }

        if (report.ratio_limit_violated) {
            assert(report.worst_offending_completion != INVALIDCOMPLETION);
            assert(report.violation_extent > 1.);
        }
    }




    template<typename TypeTag>
    template<typename RatioFunc>
    bool
    WellInterface<TypeTag>::
    checkMaxRatioLimitWell(const WellState& well_state,
                           const double max_ratio_limit,
                           const RatioFunc& ratioFunc) const
    {
        const int np = number_of_phases_;

        std::vector<double> well_rates(np, 0.0);

        for (int p = 0; p < np; ++p) {
            well_rates[p] = well_state.wellRates()[index_of_well_ * np + p];
        }

        const double well_ratio = ratioFunc(well_rates, phaseUsage());

        return (well_ratio > max_ratio_limit);
    }




    template<typename TypeTag>
    template<typename RatioFunc>
    void
    WellInterface<TypeTag>::
    checkMaxRatioLimitCompletions(const WellState& well_state,
                                  const double max_ratio_limit,
                                  const RatioFunc& ratioFunc,
                                  RatioLimitCheckReport& report) const
    {
        int worst_offending_completion = INVALIDCOMPLETION;

        // the maximum water cut value of the completions
        // it is used to identify the most offending completion
        double max_ratio_completion = 0;

        // look for the worst_offending_completion
        for (const auto& completion : completions_) {

            const int np = number_of_phases_;
            std::vector<double> completion_rates(np, 0.0);

            // looping through the connections associated with the completion
            const std::vector<int>& conns = completion.second;
            for (const int c : conns) {
                const int index_con = c + first_perf_;

                for (int p = 0; p < np; ++p) {
                    const double connection_rate = well_state.perfPhaseRates()[index_con * np + p];
                    completion_rates[p] += connection_rate;
                }
            } // end of for (const int c : conns)

            parallel_well_info_.communication().sum(completion_rates.data(), completion_rates.size());
            const double ratio_completion = ratioFunc(completion_rates, phaseUsage());

            if (ratio_completion > max_ratio_completion) {
                worst_offending_completion = completion.first;
                max_ratio_completion = ratio_completion;
            }
        } // end of for (const auto& completion : completions_)

        assert(max_ratio_completion > max_ratio_limit);
        assert(worst_offending_completion != INVALIDCOMPLETION);
        const double violation_extent = max_ratio_completion / max_ratio_limit;
        assert(violation_extent > 1.0);

        if (violation_extent > report.violation_extent) {
            report.worst_offending_completion = worst_offending_completion;
            report.violation_extent = violation_extent;
        }
    }





    template<typename TypeTag>
    void
    WellInterface<TypeTag>::
    updateWellTestState(const WellState& well_state,
                        const double& simulationTime,
                        const bool& writeMessageToOPMLog,
                        WellTestState& wellTestState,
                        Opm::DeferredLogger& deferred_logger) const
    {

        // currently, we only updateWellTestState for producers
        if (this->isInjector()) {
            return;
        }

        // Based on current understanding, only under prediction mode, we need to shut well due to various
        // reasons or limits. With more knowlage or testing cases later, this might need to be corrected.
        if (!underPredictionMode() ) {
            return;
        }

        // updating well test state based on physical (THP/BHP) limits.
        updateWellTestStatePhysical(well_state, simulationTime, writeMessageToOPMLog, wellTestState, deferred_logger);

        // updating well test state based on Economic limits.
        updateWellTestStateEconomic(well_state, simulationTime, writeMessageToOPMLog, wellTestState, deferred_logger);

        // TODO: well can be shut/closed due to other reasons
    }





    template<typename TypeTag>
    void
    WellInterface<TypeTag>::
    updateWellTestStatePhysical(const WellState& /* well_state */,
                                const double simulation_time,
                                const bool write_message_to_opmlog,
                                WellTestState& well_test_state,
                                Opm::DeferredLogger& deferred_logger) const
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





    template<typename TypeTag>
    void
    WellInterface<TypeTag>::
    updateWellTestStateEconomic(const WellState& well_state,
                                const double simulation_time,
                                const bool write_message_to_opmlog,
                                WellTestState& well_test_state,
                                Opm::DeferredLogger& deferred_logger) const
    {
        if (this->wellIsStopped())
            return;

        const WellEconProductionLimits& econ_production_limits = well_ecl_.getEconLimits();

        // if no limit is effective here, then continue to the next well
        if ( !econ_production_limits.onAnyEffectiveLimit() ) {
            return;
        }

        // flag to check if the mim oil/gas rate limit is violated
        bool rate_limit_violated = false;

        const auto& quantity_limit = econ_production_limits.quantityLimit();
        if (econ_production_limits.onAnyRateLimit()) {
            if (quantity_limit == WellEconProductionLimits::QuantityLimit::POTN)
                rate_limit_violated = checkRateEconLimits(econ_production_limits, well_state.wellPotentials(), deferred_logger);
            else {
                rate_limit_violated = checkRateEconLimits(econ_production_limits, well_state.wellRates(), deferred_logger);
            }
        }

        if (rate_limit_violated) {
            if (econ_production_limits.endRun()) {
                const std::string warning_message = std::string("ending run after well closed due to economic limits")
                                                  + std::string("is not supported yet \n")
                                                  + std::string("the program will keep running after ") + name()
                                                  + std::string(" is closed");
                deferred_logger.warning("NOT_SUPPORTING_ENDRUN", warning_message);
            }

            if (econ_production_limits.validFollowonWell()) {
                deferred_logger.warning("NOT_SUPPORTING_FOLLOWONWELL", "opening following on well after well closed is not supported yet");
            }

            well_test_state.closeWell(name(), WellTestConfig::Reason::ECONOMIC, simulation_time);
            if (write_message_to_opmlog) {
                if (well_ecl_.getAutomaticShutIn()) {
                    const std::string msg = std::string("well ") + name() + std::string(" will be shut due to rate economic limit");
                    deferred_logger.info(msg);
                } else {
                    const std::string msg = std::string("well ") + name() + std::string(" will be stopped due to rate economic limit");
                    deferred_logger.info(msg);
                }
            }
            // the well is closed, not need to check other limits
            return;
        }


        if ( !econ_production_limits.onAnyRatioLimit() ) {
            // there is no need to check the ratio limits
            return;
        }

        // checking for ratio related limits, mostly all kinds of ratio.
        RatioLimitCheckReport ratio_report;

        checkRatioEconLimits(econ_production_limits, well_state, ratio_report, deferred_logger);

        if (ratio_report.ratio_limit_violated) {
            const auto workover = econ_production_limits.workover();
            switch (workover) {
            case WellEconProductionLimits::EconWorkover::CON:
                {
                    const int worst_offending_completion = ratio_report.worst_offending_completion;

                    well_test_state.addClosedCompletion(name(), worst_offending_completion, simulation_time);
                    if (write_message_to_opmlog) {
                        if (worst_offending_completion < 0) {
                            const std::string msg = std::string("Connection ") + std::to_string(- worst_offending_completion)
                                    + std::string(" for well ") + name() + std::string(" will be closed due to economic limit");
                            deferred_logger.info(msg);
                        } else {
                            const std::string msg = std::string("Completion ") + std::to_string(worst_offending_completion)
                                    + std::string(" for well ") + name() + std::string(" will be closed due to economic limit");
                            deferred_logger.info(msg);
                        }
                    }

                    bool allCompletionsClosed = true;
                    const auto& connections = well_ecl_.getConnections();
                    for (const auto& connection : connections) {
                        if (connection.state() == Connection::State::OPEN
                            && !well_test_state.hasCompletion(name(), connection.complnum())) {
                            allCompletionsClosed = false;
                        }
                    }

                    if (allCompletionsClosed) {
                        well_test_state.closeWell(name(), WellTestConfig::Reason::ECONOMIC, simulation_time);
                        if (write_message_to_opmlog) {
                            if (well_ecl_.getAutomaticShutIn()) {
                                const std::string msg = name() + std::string(" will be shut due to last completion closed");
                                deferred_logger.info(msg);
                            } else {
                                const std::string msg = name() + std::string(" will be stopped due to last completion closed");
                                deferred_logger.info(msg);
                            }
                        }
                    }
                    break;
                }
            case WellEconProductionLimits::EconWorkover::WELL:
                {
                well_test_state.closeWell(name(), WellTestConfig::Reason::ECONOMIC, simulation_time);
                if (write_message_to_opmlog) {
                    if (well_ecl_.getAutomaticShutIn()) {
                        // tell the control that the well is closed
                        const std::string msg = name() + std::string(" will be shut due to ratio economic limit");
                        deferred_logger.info(msg);
                    } else {
                        const std::string msg = name() + std::string(" will be stopped due to ratio economic limit");
                        deferred_logger.info(msg);
                    }
                }
                    break;
                }
            case WellEconProductionLimits::EconWorkover::NONE:
                break;
                default:
                {
                    deferred_logger.warning("NOT_SUPPORTED_WORKOVER_TYPE",
                                            "not supporting workover type " + WellEconProductionLimits::EconWorkover2String(workover) );
                }
            }
        }
    }





    template<typename TypeTag>
    void
    WellInterface<TypeTag>::
    wellTesting(const Simulator& simulator,
                const double simulation_time, const int report_step,
                const WellTestConfig::Reason testing_reason,
                /* const */ WellState& well_state,
                const GroupState& group_state,
                WellTestState& well_test_state,
                Opm::DeferredLogger& deferred_logger)
    {
        if (testing_reason == WellTestConfig::Reason::PHYSICAL) {
            wellTestingPhysical(simulator, simulation_time, report_step,
                                well_state, group_state, well_test_state, deferred_logger);
        }

        if (testing_reason == WellTestConfig::Reason::ECONOMIC) {
            wellTestingEconomic(simulator, simulation_time,
                                well_state, group_state, well_test_state, deferred_logger);
        }
    }





    template<typename TypeTag>
    void
    WellInterface<TypeTag>::
    wellTestingEconomic(const Simulator& simulator,
                        const double simulation_time, const WellState& well_state, const GroupState& group_state,
                        WellTestState& welltest_state, Opm::DeferredLogger& deferred_logger)
    {
        deferred_logger.info(" well " + name() + " is being tested for economic limits");

        WellState well_state_copy = well_state;

        updateWellStateWithTarget(simulator, well_state_copy, deferred_logger);
        calculateExplicitQuantities(simulator, well_state_copy, deferred_logger);
        updatePrimaryVariables(well_state_copy, deferred_logger);
        initPrimaryVariablesEvaluation();

        WellTestState welltest_state_temp;

        bool testWell = true;
        // if a well is closed because all completions are closed, we need to check each completion
        // individually. We first open all completions, then we close one by one by calling updateWellTestState
        // untill the number of closed completions do not increase anymore.
        while (testWell) {
            const size_t original_number_closed_completions = welltest_state_temp.sizeCompletions();
            solveWellForTesting(simulator, well_state_copy, group_state, deferred_logger);
            updateWellTestState(well_state_copy, simulation_time, /*writeMessageToOPMLog=*/ false, welltest_state_temp, deferred_logger);
            closeCompletions(welltest_state_temp);

            // Stop testing if the well is closed or shut due to all completions shut
            // Also check if number of completions has increased. If the number of closed completions do not increased
            // we stop the testing.
            // TODO: it can be tricky here, if the well is shut/closed due to other reasons
            if ( welltest_state_temp.sizeWells() > 0 ||
                (original_number_closed_completions == welltest_state_temp.sizeCompletions()) ) {
                    testWell = false; // this terminates the while loop
            }
        }

        // update wellTestState if the well test succeeds
        if (!welltest_state_temp.hasWellClosed(name(), WellTestConfig::Reason::ECONOMIC)) {
            welltest_state.openWell(name(), WellTestConfig::Reason::ECONOMIC);
            const std::string msg = std::string("well ") + name() + std::string(" is re-opened through ECONOMIC testing");
            deferred_logger.info(msg);

            // also reopen completions
            for (auto& completion : well_ecl_.getCompletions()) {
                if (!welltest_state_temp.hasCompletion(name(), completion.first)) {
                    welltest_state.dropCompletion(name(), completion.first);
                }
            }
        }
    }





    template<typename TypeTag>
    void
    WellInterface<TypeTag>::
    setRepRadiusPerfLength(const std::vector<int>& cartesian_to_compressed)
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

    template<typename TypeTag>
    double
    WellInterface<TypeTag>::scalingFactor(const int phaseIdx) const
    {
        const auto& pu = phaseUsage();
        if (FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx) && pu.phase_pos[Water] == phaseIdx)
            return 1.0;
        if (FluidSystem::phaseIsActive(FluidSystem::oilPhaseIdx) && pu.phase_pos[Oil] == phaseIdx)
            return 1.0;
        if (FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx) && pu.phase_pos[Gas] == phaseIdx)
            return 0.01;
        if (has_solvent && phaseIdx == contiSolventEqIdx )
            return 0.01;

        // we should not come this far
        assert(false);
        return 1.0;
    }





    template<typename TypeTag>
    bool
    WellInterface<TypeTag>::isVFPActive(Opm::DeferredLogger& deferred_logger) const
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





    template<typename TypeTag>
    bool
    WellInterface<TypeTag>::
    iterateWellEquations(const Simulator& ebosSimulator,
                         const double dt,
                         WellState& well_state,
                         const GroupState& group_state,
                         Opm::DeferredLogger& deferred_logger)
    {
        const auto& summary_state = ebosSimulator.vanguard().summaryState();
        const auto inj_controls = well_ecl_.isInjector() ? well_ecl_.injectionControls(summary_state) : Well::InjectionControls(0);
        const auto prod_controls = well_ecl_.isProducer() ? well_ecl_.productionControls(summary_state) : Well::ProductionControls(0);

        return this->iterateWellEqWithControl(ebosSimulator, dt, inj_controls, prod_controls, well_state, group_state, deferred_logger);
    }





    template<typename TypeTag>
    void
    WellInterface<TypeTag>::calculateReservoirRates(WellState& well_state) const
    {
        const int fipreg = 0; // not considering the region for now
        const int np = number_of_phases_;

        std::vector<double> surface_rates(np, 0.0);
        const int well_rate_index = np * index_of_well_;
        for (int p = 0; p < np; ++p) {
            surface_rates[p] = well_state.wellRates()[well_rate_index + p];
        }

        std::vector<double> voidage_rates(np, 0.0);
        rateConverter_.calcReservoirVoidageRates(fipreg, pvtRegionIdx_, surface_rates, voidage_rates);

        for (int p = 0; p < np; ++p) {
            well_state.wellReservoirRates()[well_rate_index + p] = voidage_rates[p];
        }
    }

    template<typename TypeTag>
    void
    WellInterface<TypeTag>::closeCompletions(WellTestState& wellTestState)
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

    template<typename TypeTag>
    void
    WellInterface<TypeTag>::
    solveWellForTesting(const Simulator& ebosSimulator, WellState& well_state, const GroupState& group_state,
                        Opm::DeferredLogger& deferred_logger)
    {
        // keep a copy of the original well state
        const WellState well_state0 = well_state;
        const double dt = ebosSimulator.timeStepSize();
        const bool converged = iterateWellEquations(ebosSimulator, dt, well_state, group_state, deferred_logger);
        if (converged) {
            deferred_logger.debug("WellTest: Well equation for well " + name() +  " converged");
        } else {
            const int max_iter = param_.max_welleq_iter_;
            deferred_logger.debug("WellTest: Well equation for well " +name() + " failed converging in "
                          + std::to_string(max_iter) + " iterations");
            well_state = well_state0;
        }
    }

    template<typename TypeTag>
    void
    WellInterface<TypeTag>::
    solveWellEquation(const Simulator& ebosSimulator,
                      WellState& well_state,
                      const GroupState& group_state,
                      Opm::DeferredLogger& deferred_logger)
    {
        if (!this->isOperable())
            return;

        // keep a copy of the original well state
        const WellState well_state0 = well_state;
        const double dt = ebosSimulator.timeStepSize();
        const bool converged = iterateWellEquations(ebosSimulator, dt, well_state, group_state, deferred_logger);
        if (converged) {
            deferred_logger.debug("Compute initial well solution for well " + name() +  ". Converged");
        } else {
            const int max_iter = param_.max_welleq_iter_;
            deferred_logger.debug("Compute initial well solution for well " +name() + ". Failed to converge in "
                                  + std::to_string(max_iter) + " iterations");
            well_state = well_state0;
        }
    }



    template <typename TypeTag>
    void
    WellInterface<TypeTag>::
    assembleWellEq(const Simulator& ebosSimulator,
                   const double dt,
                   WellState& well_state,
                   const GroupState& group_state,
                   DeferredLogger& deferred_logger)
    {

        checkWellOperability(ebosSimulator, well_state, deferred_logger);

        if (this->useInnerIterations()) {
            this->iterateWellEquations(ebosSimulator, dt, well_state, group_state, deferred_logger);
        }

        const auto& summary_state = ebosSimulator.vanguard().summaryState();
        const auto inj_controls = well_ecl_.isInjector() ? well_ecl_.injectionControls(summary_state) : Well::InjectionControls(0);
        const auto prod_controls = well_ecl_.isProducer() ? well_ecl_.productionControls(summary_state) : Well::ProductionControls(0);
        assembleWellEqWithoutIteration(ebosSimulator, dt, inj_controls, prod_controls, well_state, group_state, deferred_logger);
    }



    template<typename TypeTag>
    void
    WellInterface<TypeTag>::addCellRates(RateVector& rates, int cellIdx) const
    {
        for (int perfIdx = 0; perfIdx < number_of_perforations_; ++perfIdx) {
            if (cells()[perfIdx] == cellIdx) {
                for (int i = 0; i < RateVector::dimension; ++i) {
                    rates[i] += connectionRates_[perfIdx][i];
                }
            }
        }
    }

    template<typename TypeTag>
    typename WellInterface<TypeTag>::Scalar
    WellInterface<TypeTag>::volumetricSurfaceRateForConnection(int cellIdx, int phaseIdx) const {
        for (int perfIdx = 0; perfIdx < number_of_perforations_; ++perfIdx) {
            if (cells()[perfIdx] == cellIdx) {
                const unsigned activeCompIdx = Indices::canonicalToActiveComponentIndex(FluidSystem::solventComponentIndex(phaseIdx));
                return connectionRates_[perfIdx][activeCompIdx].value();
            }
        }
        // this is not thread safe
        OPM_THROW(std::invalid_argument, "The well with name " + name()
                  + " does not perforate cell " + std::to_string(cellIdx));
        return 0.0;
    }


    template<typename TypeTag>
    void
    WellInterface<TypeTag>::
    wellTestingPhysical(const Simulator& ebos_simulator,
                        const double /* simulation_time */, const int /* report_step */,
                        WellState& well_state,
                        const GroupState& group_state,
                        WellTestState& welltest_state,
                        Opm::DeferredLogger& deferred_logger)
    {
        deferred_logger.info(" well " + name() + " is being tested for physical limits");

        // some most difficult things are the explicit quantities, since there is no information
        // in the WellState to do a decent initialization

        // TODO: Let us assume that the simulator is updated

        // Let us try to do a normal simualtion running, to keep checking the operability status
        // If the well is not operable during any of the time. It means it does not pass the physical
        // limit test.

        // create a copy of the well_state to use. If the operability checking is sucessful, we use this one
        // to replace the original one
        WellState well_state_copy = well_state;

        // TODO: well state for this well is kind of all zero status
        // we should be able to provide a better initialization
        calculateExplicitQuantities(ebos_simulator, well_state_copy, deferred_logger);

        updateWellOperability(ebos_simulator, well_state_copy, deferred_logger);

        if ( !this->isOperable() ) {
            const std::string msg = " well " + name() + " is not operable during well testing for physical reason";
            deferred_logger.debug(msg);
            return;
        }

        updateWellStateWithTarget(ebos_simulator, well_state_copy, deferred_logger);

        calculateExplicitQuantities(ebos_simulator, well_state_copy, deferred_logger);

        const double dt = ebos_simulator.timeStepSize();
        const bool converged = this->iterateWellEquations(ebos_simulator, dt, well_state_copy, group_state, deferred_logger);

        if (!converged) {
            const std::string msg = " well " + name() + " did not get converged during well testing for physical reason";
            deferred_logger.debug(msg);
            return;
        }

        if (this->isOperable() ) {
            welltest_state.openWell(name(), WellTestConfig::PHYSICAL );
            const std::string msg = " well " + name() + " is re-opened through well testing for physical reason";
            deferred_logger.info(msg);
            well_state = well_state_copy;
        } else {
            const std::string msg = " well " + name() + " is not operable during well testing for physical reason";
            deferred_logger.debug(msg);
        }
    }



    template<typename TypeTag>
    void
    WellInterface<TypeTag>::
    checkWellOperability(const Simulator& ebos_simulator,
                         const WellState& well_state,
                         Opm::DeferredLogger& deferred_logger)
    {

        const bool checkOperability = EWOMS_GET_PARAM(TypeTag, bool, EnableWellOperabilityCheck);
        if (!checkOperability) {
            return;
        }

        // focusing on PRODUCER for now
        if (this->isInjector()) {
            return;
        }

        if (!this->underPredictionMode() ) {
            return;
        }

        if (this->wellIsStopped() && !changed_to_stopped_this_step_) {
            return;
        }

        const bool old_well_operable = this->operability_status_.isOperable();

        updateWellOperability(ebos_simulator, well_state, deferred_logger);

        const bool well_operable = this->operability_status_.isOperable();

        if (!well_operable && old_well_operable) {
            if (well_ecl_.getAutomaticShutIn()) {
                deferred_logger.info(" well " + name() + " gets SHUT during iteration ");
            } else {
                if (!this->wellIsStopped()) {
                    deferred_logger.info(" well " + name() + " gets STOPPED during iteration ");
                    this->stopWell();
                    changed_to_stopped_this_step_ = true;
                }
            }
        } else if (well_operable && !old_well_operable) {
            deferred_logger.info(" well " + name() + " gets REVIVED during iteration ");
            this->openWell();
            changed_to_stopped_this_step_ = false;
        }
    }





    template<typename TypeTag>
    void
    WellInterface<TypeTag>::
    updateWellOperability(const Simulator& ebos_simulator,
                          const WellState& well_state,
                          Opm::DeferredLogger& deferred_logger)
    {
        this->operability_status_.reset();

        auto current_control = well_state.currentProductionControl(this->index_of_well_);
        // Operability checking is not free
        // Only check wells under BHP and THP control
        if(current_control == Well::ProducerCMode::BHP || current_control == Well::ProducerCMode::THP) {
            updateIPR(ebos_simulator, deferred_logger);
            checkOperabilityUnderBHPLimitProducer(well_state, ebos_simulator, deferred_logger);
        }
        // we do some extra checking for wells under THP control.
        if (current_control == Well::ProducerCMode::THP) {
            checkOperabilityUnderTHPLimitProducer(ebos_simulator, well_state, deferred_logger);
        }
    }


    template<typename TypeTag>
    void
    WellInterface<TypeTag>::
    updateWellStateWithTarget(const Simulator& ebos_simulator,
                              WellState& well_state,
                              Opm::DeferredLogger& deferred_logger) const
    {

        // only bhp and wellRates are used to initilize the primaryvariables for standard wells
        const auto& well = well_ecl_;
        const int well_index = index_of_well_;
        const auto& pu = phaseUsage();
        const int np = well_state.numPhases();
        const auto& summaryState = ebos_simulator.vanguard().summaryState();

        if (this->wellIsStopped()) {
            for (int p = 0; p<np; ++p) {
                well_state.wellRates()[well_index*np + p] = 0.0;
            }
            well_state.thp()[well_index] = 0.0;
            return;
        }

        if (this->isInjector() )
        {
            const auto& controls = well.injectionControls(summaryState);

            InjectorType injectorType = controls.injector_type;
            int phasePos;
            switch (injectorType) {
            case InjectorType::WATER:
            {
                phasePos = pu.phase_pos[BlackoilPhases::Aqua];
                break;
            }
            case InjectorType::OIL:
            {
                phasePos = pu.phase_pos[BlackoilPhases::Liquid];
                break;
            }
            case InjectorType::GAS:
            {
                phasePos = pu.phase_pos[BlackoilPhases::Vapour];
                break;
            }
            default:
                OPM_DEFLOG_THROW(std::runtime_error, "Expected WATER, OIL or GAS as type for injectors "  + name(), deferred_logger );
            }

            auto current = well_state.currentInjectionControl(well_index);

            switch(current) {
            case Well::InjectorCMode::RATE:
            {
                well_state.wellRates()[well_index*np + phasePos] = controls.surface_rate;
                break;
            }

            case Well::InjectorCMode::RESV:
            {
                std::vector<double> convert_coeff(number_of_phases_, 1.0);
                rateConverter_.calcCoeff(/*fipreg*/ 0, pvtRegionIdx_, convert_coeff);
                const double coeff = convert_coeff[phasePos];
                well_state.wellRates()[well_index*np + phasePos] = controls.reservoir_rate/coeff;
                break;
            }

            case Well::InjectorCMode::THP:
            {
                std::vector<double> rates(3, 0.0);
                for (int p = 0; p<np; ++p) {
                    rates[p] = well_state.wellRates()[well_index*np + p];
                }
                double bhp = calculateBhpFromThp(well_state, rates, well, summaryState, deferred_logger);
                well_state.bhp()[well_index] = bhp;

                // if the total rates are negative or zero
                // we try to provide a better intial well rate
                // using the well potentials
                double total_rate = std::accumulate(rates.begin(), rates.end(), 0.0);
                if (total_rate <= 0.0){
                    for (int p = 0; p<np; ++p) {
                        well_state.wellRates()[well_index*np + p] = well_state.wellPotentials()[well_index*np + p];
                    }
                }
                break;
            }
            case Well::InjectorCMode::BHP:
            {
                well_state.bhp()[well_index] = controls.bhp_limit;
                double total_rate = 0.0;
                for (int p = 0; p<np; ++p) {
                    total_rate += well_state.wellRates()[well_index*np + p];
                }
                // if the total rates are negative or zero
                // we try to provide a better intial well rate
                // using the well potentials
                if (total_rate <= 0.0){
                    for (int p = 0; p<np; ++p) {
                        well_state.wellRates()[well_index*np + p] = well_state.wellPotentials()[well_index*np + p];
                    }
                }
                break;
            }
            case Well::InjectorCMode::GRUP:
            {
                //do nothing at the moment
                break;
            }
            case Well::InjectorCMode::CMODE_UNDEFINED:
            {
                OPM_DEFLOG_THROW(std::runtime_error, "Well control must be specified for well "  + name(), deferred_logger );
            }

            }
        }
        //Producer
        else
        {
            auto current = well_state.currentProductionControl(well_index);
            const auto& controls = well.productionControls(summaryState);
            switch (current) {
            case Well::ProducerCMode::ORAT:
            {
                double current_rate = -well_state.wellRates()[ well_index*np + pu.phase_pos[Oil] ];
                // for trivial rates or opposite direction we don't just scale the rates
                // but use either the potentials or the mobility ratio to initial the well rates
                if (current_rate > 0.0) {
                    for (int p = 0; p<np; ++p) {
                        well_state.wellRates()[well_index*np + p] *= controls.oil_rate/current_rate;
                    }
                } else {
                    const std::vector<double> fractions = initialWellRateFrations(ebos_simulator, well_state.wellPotentials());
                    double control_fraction = fractions[pu.phase_pos[Oil]];
                    if (control_fraction != 0.0) {
                        for (int p = 0; p<np; ++p) {
                            well_state.wellRates()[well_index*np + p] = - fractions[p] * controls.oil_rate/control_fraction;
                        }
                    }
                }
                break;
            }
            case Well::ProducerCMode::WRAT:
            {
                double current_rate = -well_state.wellRates()[ well_index*np + pu.phase_pos[Water] ];
                // for trivial rates or opposite direction we don't just scale the rates
                // but use either the potentials or the mobility ratio to initial the well rates
                if (current_rate > 0.0) {
                    for (int p = 0; p<np; ++p) {
                        well_state.wellRates()[well_index*np + p] *= controls.water_rate/current_rate;
                    }
                } else {
                    const std::vector<double> fractions = initialWellRateFrations(ebos_simulator, well_state.wellPotentials());
                    double control_fraction = fractions[pu.phase_pos[Water]];
                    if (control_fraction != 0.0) {
                        for (int p = 0; p<np; ++p) {
                            well_state.wellRates()[well_index*np + p] = - fractions[p] * controls.water_rate/control_fraction;
                        }
                    }
                }
                break;
            }
            case Well::ProducerCMode::GRAT:
            {
                double current_rate = -well_state.wellRates()[ well_index*np + pu.phase_pos[Gas] ];
                // or trivial rates or opposite direction we don't just scale the rates
                // but use either the potentials or the mobility ratio to initial the well rates
                if (current_rate > 0.0) {
                    for (int p = 0; p<np; ++p) {
                        well_state.wellRates()[well_index*np + p] *= controls.gas_rate/current_rate;
                    }
                } else {
                    const std::vector<double> fractions = initialWellRateFrations(ebos_simulator, well_state.wellPotentials());
                    double control_fraction = fractions[pu.phase_pos[Gas]];
                    if (control_fraction != 0.0) {
                        for (int p = 0; p<np; ++p) {
                            well_state.wellRates()[well_index*np + p] = - fractions[p] * controls.gas_rate/control_fraction;
                        }
                    }
                }

                break;

            }
            case Well::ProducerCMode::LRAT:
            {
                double current_rate = -well_state.wellRates()[ well_index*np + pu.phase_pos[Water] ]
                        - well_state.wellRates()[ well_index*np + pu.phase_pos[Oil] ];
                // or trivial rates or opposite direction we don't just scale the rates
                // but use either the potentials or the mobility ratio to initial the well rates
                if (current_rate > 0.0) {
                    for (int p = 0; p<np; ++p) {
                        well_state.wellRates()[well_index*np + p] *= controls.liquid_rate/current_rate;
                    }
                } else {
                    const std::vector<double> fractions = initialWellRateFrations(ebos_simulator, well_state.wellPotentials());
                    double control_fraction = fractions[pu.phase_pos[Water]] + fractions[pu.phase_pos[Oil]];
                    if (control_fraction != 0.0) {
                        for (int p = 0; p<np; ++p) {
                            well_state.wellRates()[well_index*np + p] = - fractions[p] * controls.liquid_rate / control_fraction;
                        }
                    }
                }
                break;
            }
            case Well::ProducerCMode::CRAT:
            {
                OPM_DEFLOG_THROW(std::runtime_error, "CRAT control not supported " << name(), deferred_logger);
            }
            case Well::ProducerCMode::RESV:
            {
                std::vector<double> convert_coeff(number_of_phases_, 1.0);
                rateConverter_.calcCoeff(/*fipreg*/ 0, pvtRegionIdx_, convert_coeff);
                double total_res_rate = 0.0;
                for (int p = 0; p<np; ++p) {
                    total_res_rate -= well_state.wellRates()[well_index*np + p] * convert_coeff[p];
                }
                if (controls.prediction_mode) {
                    // or trivial rates or opposite direction we don't just scale the rates
                    // but use either the potentials or the mobility ratio to initial the well rates
                    if (total_res_rate > 0.0) {
                        for (int p = 0; p<np; ++p) {
                            well_state.wellRates()[well_index*np + p] *= controls.resv_rate/total_res_rate;
                        }
                    } else {
                        const std::vector<double> fractions = initialWellRateFrations(ebos_simulator, well_state.wellPotentials());
                        for (int p = 0; p<np; ++p) {
                            well_state.wellRates()[well_index*np + p] = - fractions[p] * controls.resv_rate / convert_coeff[p];
                        }
                    }
                } else {
                    std::vector<double> hrates(number_of_phases_,0.);
                    if (FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx)) {
                        hrates[pu.phase_pos[Water]] = controls.water_rate;
                    }
                    if (FluidSystem::phaseIsActive(FluidSystem::oilPhaseIdx)) {
                        hrates[pu.phase_pos[Oil]] = controls.oil_rate;
                    }
                    if (FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx)) {
                        hrates[pu.phase_pos[Gas]] = controls.gas_rate;
                    }
                    std::vector<double> hrates_resv(number_of_phases_,0.);
                    rateConverter_.calcReservoirVoidageRates(/*fipreg*/ 0, pvtRegionIdx_, hrates, hrates_resv);
                    double target = std::accumulate(hrates_resv.begin(), hrates_resv.end(), 0.0);
                    // or trivial rates or opposite direction we don't just scale the rates
                    // but use either the potentials or the mobility ratio to initial the well rates
                    if (total_res_rate > 0.0) {
                        for (int p = 0; p<np; ++p) {
                            well_state.wellRates()[well_index*np + p] *= target/total_res_rate;
                        }
                    } else {
                        const std::vector<double> fractions = initialWellRateFrations(ebos_simulator, well_state.wellPotentials());
                        for (int p = 0; p<np; ++p) {
                            well_state.wellRates()[well_index*np + p] = - fractions[p] * target / convert_coeff[p];
                        }
                    }

                }
                break;
            }
            case Well::ProducerCMode::BHP:
            {
                well_state.bhp()[well_index] = controls.bhp_limit;
                double total_rate = 0.0;
                for (int p = 0; p<np; ++p) {
                    total_rate -= well_state.wellRates()[well_index*np + p];
                }
                // if the total rates are negative or zero
                // we try to provide a better intial well rate
                // using the well potentials
                if (total_rate <= 0.0){
                    for (int p = 0; p<np; ++p) {
                        well_state.wellRates()[well_index*np + p] = -well_state.wellPotentials()[well_index*np + p];
                    }
                }
                break;
            }
            case Well::ProducerCMode::THP:
            {
                std::vector<double> rates(3, 0.0);
                for (int p = 0; p<np; ++p) {
                    rates[p] = well_state.wellRates()[well_index*np + p];
                }
                double bhp = calculateBhpFromThp(well_state, rates, well, summaryState, deferred_logger);
                well_state.bhp()[well_index] = bhp;

                // if the total rates are negative or zero
                // we try to provide a better intial well rate
                // using the well potentials
                double total_rate = -std::accumulate(rates.begin(), rates.end(), 0.0);
                if (total_rate <= 0.0){
                    for (int p = 0; p<np; ++p) {
                        well_state.wellRates()[well_index*np + p] = -well_state.wellPotentials()[well_index*np + p];
                    }
                }
                break;
            }
            case Well::ProducerCMode::GRUP:
            {
                //do nothing at the moment
                break;
            }
            case Well::ProducerCMode::CMODE_UNDEFINED:
            case Well::ProducerCMode::NONE:
            {
                OPM_DEFLOG_THROW(std::runtime_error, "Well control must be specified for well "  + name() , deferred_logger);
            }

                break;
            } // end of switch
        }
    }

    template<typename TypeTag>
    std::vector<double>
    WellInterface<TypeTag>::
    initialWellRateFrations(const Simulator& ebosSimulator, const std::vector<double>& potentials) const
    {
        const int np = number_of_phases_;
        std::vector<double> scaling_factor(np);

        double total_potentials = 0.0;
        for (int p = 0; p<np; ++p) {
            total_potentials += potentials[index_of_well_*np + p];
        }
        if (total_potentials > 0) {
            for (int p = 0; p<np; ++p) {
                scaling_factor[p] = potentials[index_of_well_*np + p] / total_potentials;
            }
            return scaling_factor;
        }
        // if we don't have any potentials we weight it using the mobilites
        // We only need approximation so we don't bother with the vapporized oil and dissolved gas
        double total_tw = 0;
        const int nperf = number_of_perforations_;
        for (int perf = 0; perf < nperf; ++perf) {
            total_tw += well_index_[perf];
        }
        for (int perf = 0; perf < nperf; ++perf) {
            const int cell_idx = well_cells_[perf];
            const auto& intQuants = *(ebosSimulator.model().cachedIntensiveQuantities(cell_idx, /*timeIdx=*/0));
            const auto& fs = intQuants.fluidState();
            const double well_tw_fraction = well_index_[perf] / total_tw;
            double total_mobility = 0.0;
            for (int p = 0; p < np; ++p) {
                int ebosPhaseIdx = flowPhaseToEbosPhaseIdx(p);
                total_mobility += fs.invB(ebosPhaseIdx).value() * intQuants.mobility(ebosPhaseIdx).value();
            }
            for (int p = 0; p < np; ++p) {
                int ebosPhaseIdx = flowPhaseToEbosPhaseIdx(p);
                scaling_factor[p] += well_tw_fraction * fs.invB(ebosPhaseIdx).value() * intQuants.mobility(ebosPhaseIdx).value() / total_mobility;
            }
        }
        return scaling_factor;
    }

    template<typename TypeTag>
    bool
    WellInterface<TypeTag>::
    isOperable() const {
        return operability_status_.isOperable();
    }





    template <typename TypeTag>
    bool
    WellInterface<TypeTag>::checkConstraints(WellState& well_state,
                                             const GroupState& group_state,
                                             const Schedule& schedule,
                                             const SummaryState& summaryState,
                                             DeferredLogger& deferred_logger) const
    {
        const bool ind_broken = checkIndividualConstraints(well_state, summaryState);
        if (ind_broken) {
            return true;
        } else {
            return checkGroupConstraints(well_state, group_state, schedule, summaryState, deferred_logger);
        }
    }





    template <typename TypeTag>
    bool
    WellInterface<TypeTag>::checkIndividualConstraints(WellState& well_state,
                                                       const SummaryState& summaryState) const
    {
        const auto& well = well_ecl_;
        const PhaseUsage& pu = phaseUsage();
        const int well_index = index_of_well_;
        const auto wellrate_index = well_index * pu.num_phases;

        if (well.isInjector()) {
            const auto controls = well.injectionControls(summaryState);
            auto currentControl = well_state.currentInjectionControl(well_index);

            if (controls.hasControl(Well::InjectorCMode::BHP) && currentControl != Well::InjectorCMode::BHP)
            {
                const auto& bhp = controls.bhp_limit;
                double current_bhp = well_state.bhp()[well_index];
                if (bhp < current_bhp) {
                    well_state.currentInjectionControl(well_index, Well::InjectorCMode::BHP);
                    return true;
                }
            }

            if (controls.hasControl(Well::InjectorCMode::RATE) && currentControl != Well::InjectorCMode::RATE)
            {
                InjectorType injectorType = controls.injector_type;
                double current_rate = 0.0;

                switch (injectorType) {
                case InjectorType::WATER:
                {
                    current_rate = well_state.wellRates()[ wellrate_index + pu.phase_pos[BlackoilPhases::Aqua] ];
                    break;
                }
                case InjectorType::OIL:
                {
                    current_rate = well_state.wellRates()[ wellrate_index + pu.phase_pos[BlackoilPhases::Liquid] ];
                    break;
                }
                case InjectorType::GAS:
                {
                    current_rate = well_state.wellRates()[ wellrate_index + pu.phase_pos[BlackoilPhases::Vapour] ];
                    break;
                }
                default:
                    throw("Expected WATER, OIL or GAS as type for injectors " + well.name());
                }

                if (controls.surface_rate < current_rate) {
                    well_state.currentInjectionControl(well_index, Well::InjectorCMode::RATE);
                    return true;
                }

            }

            if (controls.hasControl(Well::InjectorCMode::RESV) && currentControl != Well::InjectorCMode::RESV)
            {
                double current_rate = 0.0;
                if( pu.phase_used[BlackoilPhases::Aqua] )
                    current_rate += well_state.wellReservoirRates()[ wellrate_index + pu.phase_pos[BlackoilPhases::Aqua] ];

                if( pu.phase_used[BlackoilPhases::Liquid] )
                    current_rate += well_state.wellReservoirRates()[ wellrate_index + pu.phase_pos[BlackoilPhases::Liquid] ];

                if( pu.phase_used[BlackoilPhases::Vapour] )
                    current_rate += well_state.wellReservoirRates()[ wellrate_index + pu.phase_pos[BlackoilPhases::Vapour] ];

                if (controls.reservoir_rate < current_rate) {
                    currentControl = Well::InjectorCMode::RESV;
                    return true;
                }
            }

            if (controls.hasControl(Well::InjectorCMode::THP) && currentControl != Well::InjectorCMode::THP)
            {
                const auto& thp = this->getTHPConstraint(summaryState);
                double current_thp = well_state.thp()[well_index];
                if (thp < current_thp) {
                    currentControl = Well::InjectorCMode::THP;
                    return true;
                }
            }

        }

        if (well.isProducer( )) {
            const auto controls = well.productionControls(summaryState);
            auto currentControl = well_state.currentProductionControl(well_index);

            if (controls.hasControl(Well::ProducerCMode::BHP) && currentControl != Well::ProducerCMode::BHP )
            {
                const double bhp = controls.bhp_limit;
                double current_bhp = well_state.bhp()[well_index];
                if (bhp > current_bhp) {
                    well_state.currentProductionControl(well_index, Well::ProducerCMode::BHP);
                    return true;
                }
            }

            if (controls.hasControl(Well::ProducerCMode::ORAT) && currentControl != Well::ProducerCMode::ORAT) {
                double current_rate = -well_state.wellRates()[ wellrate_index + pu.phase_pos[BlackoilPhases::Liquid] ];
                if (controls.oil_rate < current_rate  ) {
                    well_state.currentProductionControl(well_index, Well::ProducerCMode::ORAT);
                    return true;
                }
            }

            if (controls.hasControl(Well::ProducerCMode::WRAT) && currentControl != Well::ProducerCMode::WRAT ) {
                double current_rate = -well_state.wellRates()[ wellrate_index + pu.phase_pos[BlackoilPhases::Aqua] ];
                if (controls.water_rate < current_rate  ) {
                    well_state.currentProductionControl(well_index, Well::ProducerCMode::WRAT);
                    return true;
                }
            }

            if (controls.hasControl(Well::ProducerCMode::GRAT) && currentControl != Well::ProducerCMode::GRAT ) {
                double current_rate = -well_state.wellRates()[ wellrate_index + pu.phase_pos[BlackoilPhases::Vapour] ];
                if (controls.gas_rate < current_rate  ) {
                    well_state.currentProductionControl(well_index, Well::ProducerCMode::GRAT);
                    return true;
                }
            }

            if (controls.hasControl(Well::ProducerCMode::LRAT) && currentControl != Well::ProducerCMode::LRAT) {
                double current_rate = -well_state.wellRates()[ wellrate_index + pu.phase_pos[BlackoilPhases::Liquid] ];
                current_rate -= well_state.wellRates()[ wellrate_index + pu.phase_pos[BlackoilPhases::Aqua] ];
                if (controls.liquid_rate < current_rate  ) {
                    well_state.currentProductionControl(well_index, Well::ProducerCMode::LRAT);
                    return true;
                }
            }

            if (controls.hasControl(Well::ProducerCMode::RESV) && currentControl != Well::ProducerCMode::RESV ) {
                double current_rate = 0.0;
                if( pu.phase_used[BlackoilPhases::Aqua] )
                    current_rate -= well_state.wellReservoirRates()[ wellrate_index + pu.phase_pos[BlackoilPhases::Aqua] ];

                if( pu.phase_used[BlackoilPhases::Liquid] )
                    current_rate -= well_state.wellReservoirRates()[ wellrate_index + pu.phase_pos[BlackoilPhases::Liquid] ];

                if( pu.phase_used[BlackoilPhases::Vapour] )
                    current_rate -= well_state.wellReservoirRates()[ wellrate_index + pu.phase_pos[BlackoilPhases::Vapour] ];

                if (controls.prediction_mode && controls.resv_rate > current_rate) {
                    well_state.currentProductionControl(well_index, Well::ProducerCMode::RESV);
                    return true;
                }

                if (!controls.prediction_mode) {
                    const int fipreg = 0; // not considering the region for now
                    const int np = number_of_phases_;

                    std::vector<double> surface_rates(np, 0.0);
                    if( pu.phase_used[BlackoilPhases::Aqua] )
                        surface_rates[pu.phase_pos[BlackoilPhases::Aqua]] = controls.water_rate;
                    if( pu.phase_used[BlackoilPhases::Liquid] )
                        surface_rates[pu.phase_pos[BlackoilPhases::Liquid]] = controls.oil_rate;
                    if( pu.phase_used[BlackoilPhases::Vapour] )
                        surface_rates[pu.phase_pos[BlackoilPhases::Vapour]] = controls.gas_rate;

                    std::vector<double> voidage_rates(np, 0.0);
                    rateConverter_.calcReservoirVoidageRates(fipreg, pvtRegionIdx_, surface_rates, voidage_rates);

                    double resv_rate = 0.0;
                    for (int p = 0; p < np; ++p) {
                        resv_rate += voidage_rates[p];
                    }

                    if (resv_rate < current_rate) {
                        well_state.currentProductionControl(well_index, Well::ProducerCMode::RESV);
                        return true;
                    }
                }
            }

            if (controls.hasControl(Well::ProducerCMode::THP) && currentControl != Well::ProducerCMode::THP)
            {
                const auto& thp = this->getTHPConstraint(summaryState);
                double current_thp =  well_state.thp()[well_index];
                if (thp > current_thp) {
                    well_state.currentProductionControl(well_index, Well::ProducerCMode::THP);
                    return true;
                }
            }

        }

        return false;
    }





    template <typename TypeTag>
    bool
    WellInterface<TypeTag>::checkGroupConstraints(WellState& well_state,
                                                  const GroupState& group_state,
                                                  const Schedule& schedule,
                                                  const SummaryState& summaryState,
                                                  DeferredLogger& deferred_logger) const
    {
        const auto& well = well_ecl_;
        const int well_index = index_of_well_;

        if (well.isInjector()) {
            auto currentControl = well_state.currentInjectionControl(well_index);

            if (currentControl != Well::InjectorCMode::GRUP) {
                // This checks only the first encountered group limit,
                // in theory there could be several, and then we should
                // test all but the one currently applied. At that point,
                // this if-statement should be removed and we should always
                // check, skipping over only the single group parent whose
                // control is the active one for the well (if any).
                const auto& group = schedule.getGroup( well.groupName(), current_step_ );
                const double efficiencyFactor = well.getEfficiencyFactor();
                const std::pair<bool, double> group_constraint = checkGroupConstraintsInj(
                                                                                          group, well_state, group_state, efficiencyFactor, schedule, summaryState, deferred_logger);
                // If a group constraint was broken, we set the current well control to
                // be GRUP.
                if (group_constraint.first) {
                    well_state.currentInjectionControl(index_of_well_, Well::InjectorCMode::GRUP);
                    const int np = well_state.numPhases();
                    for (int p = 0; p<np; ++p) {
                        well_state.wellRates()[index_of_well_*np + p] *= group_constraint.second;
                    }
                }
                return group_constraint.first;
            }
        }

        if (well.isProducer( )) {
            auto currentControl = well_state.currentProductionControl(well_index);

            if (currentControl != Well::ProducerCMode::GRUP) {
                // This checks only the first encountered group limit,
                // in theory there could be several, and then we should
                // test all but the one currently applied. At that point,
                // this if-statement should be removed and we should always
                // check, skipping over only the single group parent whose
                // control is the active one for the well (if any).
                const auto& group = schedule.getGroup( well.groupName(), current_step_ );
                const double efficiencyFactor = well.getEfficiencyFactor();
                const std::pair<bool, double> group_constraint = checkGroupConstraintsProd(
                                                                                           group, well_state, group_state, efficiencyFactor, schedule, summaryState, deferred_logger);
                // If a group constraint was broken, we set the current well control to
                // be GRUP.
                if (group_constraint.first) {
                    well_state.currentProductionControl(index_of_well_, Well::ProducerCMode::GRUP);
                    const int np = well_state.numPhases();
                    for (int p = 0; p<np; ++p) {
                        well_state.wellRates()[index_of_well_*np + p] *= group_constraint.second;
                    }
                }
                return group_constraint.first;
            }
        }

        return false;
    }





    template <typename TypeTag>
    std::pair<bool, double>
    WellInterface<TypeTag>::checkGroupConstraintsInj(const Group& group,
                                                     const WellState& well_state,
                                                     const GroupState& group_state,
                                                     const double efficiencyFactor,
                                                     const Schedule& schedule,
                                                     const SummaryState& summaryState,
                                                     DeferredLogger& deferred_logger) const
    {
        // Translate injector type from control to Phase.
        const auto& well_controls = well_ecl_.injectionControls(summaryState);
        auto injectorType = well_controls.injector_type;
        Phase injectionPhase;
        switch (injectorType) {
        case InjectorType::WATER:
        {
            injectionPhase = Phase::WATER;
            break;
        }
        case InjectorType::OIL:
        {
            injectionPhase = Phase::OIL;
            break;
        }
        case InjectorType::GAS:
        {
            injectionPhase = Phase::GAS;
            break;
        }
        default:
            throw("Expected WATER, OIL or GAS as type for injector " + name());
        }

        // Make conversion factors for RESV <-> surface rates.
        std::vector<double> resv_coeff(phaseUsage().num_phases, 1.0);
        rateConverter_.calcCoeff(0, pvtRegionIdx_, resv_coeff); // FIPNUM region 0 here, should use FIPNUM from WELSPECS.

        // Call check for the well's injection phase.
        return WellGroupHelpers::checkGroupConstraintsInj(name(),
                                                          well_ecl_.groupName(),
                                                          group,
                                                          well_state,
                                                          group_state,
                                                          current_step_,
                                                          guide_rate_,
                                                          well_state.wellRates().data() + index_of_well_ * phaseUsage().num_phases,
                                                          injectionPhase,
                                                          phaseUsage(),
                                                          efficiencyFactor,
                                                          schedule,
                                                          summaryState,
                                                          resv_coeff,
                                                          deferred_logger);
    }





    template <typename TypeTag>
    std::pair<bool, double>
    WellInterface<TypeTag>::checkGroupConstraintsProd(const Group& group,
                                                      const WellState& well_state,
                                                      const GroupState& group_state,
                                                      const double efficiencyFactor,
                                                      const Schedule& schedule,
                                                      const SummaryState& summaryState,
                                                      DeferredLogger& deferred_logger) const
    {
        // Make conversion factors for RESV <-> surface rates.
        std::vector<double> resv_coeff(phaseUsage().num_phases, 1.0);
        rateConverter_.calcCoeff(0, pvtRegionIdx_, resv_coeff); // FIPNUM region 0 here, should use FIPNUM from WELSPECS.

        return WellGroupHelpers::checkGroupConstraintsProd(name(),
                                                           well_ecl_.groupName(),
                                                           group,
                                                           well_state,
                                                           group_state,
                                                           current_step_,
                                                           guide_rate_,
                                                           well_state.wellRates().data() + index_of_well_ * phaseUsage().num_phases,
                                                           phaseUsage(),
                                                           efficiencyFactor,
                                                           schedule,
                                                           summaryState,
                                                           resv_coeff,
                                                           deferred_logger);
    }





    template <typename TypeTag>
    template <class EvalWell, class BhpFromThpFunc>
    void
    WellInterface<TypeTag>::assembleControlEqInj(const WellState& well_state,
                                                 const GroupState& group_state,
                                                 const Opm::Schedule& schedule,
                                                 const SummaryState& summaryState,
                                                 const Well::InjectionControls& controls,
                                                 const EvalWell& bhp,
                                                 const EvalWell& injection_rate,
                                                 BhpFromThpFunc bhp_from_thp,
                                                 EvalWell& control_eq,
                                                 Opm::DeferredLogger& deferred_logger)
    {
        auto current = well_state.currentInjectionControl(index_of_well_);
        const InjectorType injectorType = controls.injector_type;
        const auto& pu = phaseUsage();
        const double efficiencyFactor = well_ecl_.getEfficiencyFactor();

        switch (current) {
        case Well::InjectorCMode::RATE: {
            control_eq = injection_rate - controls.surface_rate;
            break;
        }
        case Well::InjectorCMode::RESV: {
            std::vector<double> convert_coeff(number_of_phases_, 1.0);
            rateConverter_.calcCoeff(/*fipreg*/ 0, pvtRegionIdx_, convert_coeff);

            double coeff = 1.0;

            switch (injectorType) {
            case InjectorType::WATER: {
                coeff = convert_coeff[pu.phase_pos[BlackoilPhases::Aqua]];
                break;
            }
            case InjectorType::OIL: {
                coeff = convert_coeff[pu.phase_pos[BlackoilPhases::Liquid]];
                break;
            }
            case InjectorType::GAS: {
                coeff = convert_coeff[pu.phase_pos[BlackoilPhases::Vapour]];
                break;
            }
            default:
                throw("Expected WATER, OIL or GAS as type for injectors " + well_ecl_.name());
            }

            control_eq = coeff * injection_rate - controls.reservoir_rate;
            break;
        }
        case Well::InjectorCMode::THP: {
            control_eq = bhp - bhp_from_thp();
            break;
        }
        case Well::InjectorCMode::BHP: {
            control_eq = bhp - controls.bhp_limit;
            break;
        }
        case Well::InjectorCMode::GRUP: {
            assert(well_ecl_.isAvailableForGroupControl());
            const auto& group = schedule.getGroup(well_ecl_.groupName(), current_step_);
            getGroupInjectionControl(group,
                                     well_state,
                                     group_state,
                                     schedule,
                                     summaryState,
                                     injectorType,
                                     bhp,
                                     injection_rate,
                                     control_eq,
                                     efficiencyFactor,
                                     deferred_logger);
            break;
        }
        case Well::InjectorCMode::CMODE_UNDEFINED: {
            OPM_DEFLOG_THROW(std::runtime_error, "Well control must be specified for well " + name(), deferred_logger);
        }
        }
    }




    template <typename TypeTag>
    template <class EvalWell, class BhpFromThpFunc>
    void
    WellInterface<TypeTag>::assembleControlEqProd(const WellState& well_state,
                                                  const GroupState& group_state,
                                                  const Opm::Schedule& schedule,
                                                  const SummaryState& summaryState,
                                                  const Well::ProductionControls& controls,
                                                  const EvalWell& bhp,
                                                  const std::vector<EvalWell>& rates, // Always 3 canonical rates.
                                                  BhpFromThpFunc bhp_from_thp,
                                                  EvalWell& control_eq,
                                                  Opm::DeferredLogger& deferred_logger)
    {
        auto current = well_state.currentProductionControl(index_of_well_);
        const auto& pu = phaseUsage();
        const double efficiencyFactor = well_ecl_.getEfficiencyFactor();

        switch (current) {
        case Well::ProducerCMode::ORAT: {
            assert(FluidSystem::phaseIsActive(FluidSystem::oilPhaseIdx));
            const EvalWell rate = -rates[BlackoilPhases::Liquid];
            control_eq = rate - controls.oil_rate;
            break;
        }
        case Well::ProducerCMode::WRAT: {
            assert(FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx));
            const EvalWell rate = -rates[BlackoilPhases::Aqua];
            control_eq = rate - controls.water_rate;
            break;
        }
        case Well::ProducerCMode::GRAT: {
            assert(FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx));
            const EvalWell rate = -rates[BlackoilPhases::Vapour];
            control_eq = rate - controls.gas_rate;
            break;
        }
        case Well::ProducerCMode::LRAT: {
            assert(FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx));
            assert(FluidSystem::phaseIsActive(FluidSystem::oilPhaseIdx));
            EvalWell rate = -rates[BlackoilPhases::Aqua] - rates[BlackoilPhases::Liquid];
            control_eq = rate - controls.liquid_rate;
            break;
        }
        case Well::ProducerCMode::CRAT: {
            OPM_DEFLOG_THROW(std::runtime_error, "CRAT control not supported " << name(), deferred_logger);
        }
        case Well::ProducerCMode::RESV: {
            auto total_rate = rates[0]; // To get the correct type only.
            total_rate = 0.0;
            std::vector<double> convert_coeff(number_of_phases_, 1.0);
            rateConverter_.calcCoeff(/*fipreg*/ 0, pvtRegionIdx_, convert_coeff);
            for (int phase = 0; phase < 3; ++phase) {
                if (pu.phase_used[phase]) {
                    const int pos = pu.phase_pos[phase];
                    total_rate -= rates[phase] * convert_coeff[pos]; // Note different indices.
                }
            }
            if (controls.prediction_mode) {
                control_eq = total_rate - controls.resv_rate;
            } else {
                std::vector<double> hrates(number_of_phases_, 0.);
                if (FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx)) {
                    hrates[pu.phase_pos[Water]] = controls.water_rate;
                }
                if (FluidSystem::phaseIsActive(FluidSystem::oilPhaseIdx)) {
                    hrates[pu.phase_pos[Oil]] = controls.oil_rate;
                }
                if (FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx)) {
                    hrates[pu.phase_pos[Gas]] = controls.gas_rate;
                }
                std::vector<double> hrates_resv(number_of_phases_, 0.);
                rateConverter_.calcReservoirVoidageRates(/*fipreg*/ 0, pvtRegionIdx_, hrates, hrates_resv);
                double target = std::accumulate(hrates_resv.begin(), hrates_resv.end(), 0.0);
                control_eq = total_rate - target;
            }
            break;
        }
        case Well::ProducerCMode::BHP: {
            control_eq = bhp - controls.bhp_limit;
            break;
        }
        case Well::ProducerCMode::THP: {
            control_eq = bhp - bhp_from_thp();
            break;
        }
        case Well::ProducerCMode::GRUP: {
            assert(well_ecl_.isAvailableForGroupControl());
            const auto& group = schedule.getGroup(well_ecl_.groupName(), current_step_);
            // Annoying thing: the rates passed to this function are
            // always of size 3 and in canonical (for PhaseUsage)
            // order. This is what is needed for VFP calculations if
            // they are required (THP controlled well). But for the
            // group production control things we must pass only the
            // active phases' rates.
            std::vector<EvalWell> active_rates(pu.num_phases);
            for (int canonical_phase = 0; canonical_phase < 3; ++canonical_phase) {
                if (pu.phase_used[canonical_phase]) {
                    active_rates[pu.phase_pos[canonical_phase]] = rates[canonical_phase];
                }
            }
            getGroupProductionControl(group, well_state, group_state, schedule, summaryState, bhp, active_rates, control_eq, efficiencyFactor);
            break;
        }
        case Well::ProducerCMode::CMODE_UNDEFINED: {
            OPM_DEFLOG_THROW(std::runtime_error, "Well control must be specified for well " + name(), deferred_logger);
        }
        case Well::ProducerCMode::NONE: {
            OPM_DEFLOG_THROW(std::runtime_error, "Well control must be specified for well " + name(), deferred_logger);
        }
        }
    }



    template <typename TypeTag>
    template <class EvalWell>
    void
    WellInterface<TypeTag>::getGroupInjectionControl(const Group& group,
                                                      const WellState& well_state,
                                                     const GroupState& group_state,
                                                      const Opm::Schedule& schedule,
                                                      const SummaryState& summaryState,
                                                      const InjectorType& injectorType,
                                                      const EvalWell& bhp,
                                                      const EvalWell& injection_rate,
                                                      EvalWell& control_eq,
                                                      double efficiencyFactor,
                                                      Opm::DeferredLogger& deferred_logger)
    {
        // Setting some defaults to silence warnings below.
        // Will be overwritten in the switch statement.
        Phase injectionPhase = Phase::WATER;
        switch (injectorType) {
        case InjectorType::WATER:
        {
            injectionPhase = Phase::WATER;
            break;
        }
        case InjectorType::OIL:
        {
            injectionPhase = Phase::OIL;
            break;
        }
        case InjectorType::GAS:
        {
            injectionPhase = Phase::GAS;
            break;
        }
        default:
            // Should not be here.
            assert(false);
        }

        auto currentGroupControl = group_state.injection_control(group.name(), injectionPhase);
        if (currentGroupControl == Group::InjectionCMode::FLD ||
            currentGroupControl == Group::InjectionCMode::NONE) {
            if (!group.injectionGroupControlAvailable(injectionPhase)) {
                // We cannot go any further up the hierarchy. This could
                // be the FIELD group, or any group for which this has
                // been set in GCONINJE or GCONPROD. If we are here
                // anyway, it is likely that the deck set inconsistent
                // requirements, such as GRUP control mode on a well with
                // no appropriate controls defined on any of its
                // containing groups. We will therefore use the wells' bhp
                // limit equation as a fallback.
                const auto& controls = well_ecl_.injectionControls(summaryState);
                control_eq = bhp - controls.bhp_limit;
                return;
            } else {
                // Inject share of parents control
                const auto& parent = schedule.getGroup( group.parent(), current_step_ );
                efficiencyFactor *= group.getGroupEfficiencyFactor();
                getGroupInjectionControl(parent, well_state, group_state, schedule, summaryState, injectorType, bhp, injection_rate, control_eq, efficiencyFactor, deferred_logger);
                return;
            }
        }

        efficiencyFactor *= group.getGroupEfficiencyFactor();
        const auto& well = well_ecl_;
        const auto pu = phaseUsage();

        if (!group.isInjectionGroup()) {
            // use bhp as control eq and let the updateControl code find a valid control
            const auto& controls = well.injectionControls(summaryState);
            control_eq = bhp - controls.bhp_limit;
            return;
        }

        // If we are here, we are at the topmost group to be visited in the recursion.
        // This is the group containing the control we will check against.

        // Make conversion factors for RESV <-> surface rates.
        std::vector<double> resv_coeff(phaseUsage().num_phases, 1.0);
        rateConverter_.calcCoeff(0, pvtRegionIdx_, resv_coeff); // FIPNUM region 0 here, should use FIPNUM from WELSPECS.

        double sales_target = 0;
        if (schedule[current_step_].gconsale().has(group.name())) {
            const auto& gconsale = schedule[current_step_].gconsale().get(group.name(), summaryState);
            sales_target = gconsale.sales_target;
        }
        WellGroupHelpers::InjectionTargetCalculator tcalc(currentGroupControl, pu, resv_coeff, group.name(), sales_target, group_state, injectionPhase, deferred_logger);
        WellGroupHelpers::FractionCalculator fcalc(schedule, well_state, group_state, current_step_, guide_rate_, tcalc.guideTargetMode(), pu, false, injectionPhase);

        auto localFraction = [&](const std::string& child) {
            return fcalc.localFraction(child, "");
        };

        auto localReduction = [&](const std::string& group_name) {
            const std::vector<double>& groupTargetReductions = group_state.injection_reduction_rates(group_name);
            return tcalc.calcModeRateFromRates(groupTargetReductions);
        };

        const double orig_target = tcalc.groupTarget(group.injectionControls(injectionPhase, summaryState), deferred_logger);
        const auto chain = WellGroupHelpers::groupChainTopBot(name(), group.name(), schedule, current_step_);
        // Because 'name' is the last of the elements, and not an ancestor, we subtract one below.
        const size_t num_ancestors = chain.size() - 1;
        double target = orig_target;
        for (size_t ii = 0; ii < num_ancestors; ++ii) {
            if ((ii == 0) || guide_rate_->has(chain[ii], injectionPhase)) {
                // Apply local reductions only at the control level
                // (top) and for levels where we have a specified
                // group guide rate.
                target -= localReduction(chain[ii]);
            }
            target *= localFraction(chain[ii+1]);
        }
        // Avoid negative target rates coming from too large local reductions.
        const double target_rate = std::max(0.0, target / efficiencyFactor);
        const auto current_rate = injection_rate; // Switch sign since 'rates' are negative for producers.
        control_eq = current_rate - target_rate;
    }



    template <typename TypeTag>
    template <class EvalWell>
    void
    WellInterface<TypeTag>::getGroupProductionControl(const Group& group,
                                                      const WellState& well_state,
                                                      const GroupState& group_state,
                                                      const Opm::Schedule& schedule,
                                                      const SummaryState& summaryState,
                                                      const EvalWell& bhp,
                                                      const std::vector<EvalWell>& rates,
                                                      EvalWell& control_eq,
                                                      double efficiencyFactor)
    {
        const Group::ProductionCMode& currentGroupControl = group_state.production_control(group.name());
        if (currentGroupControl == Group::ProductionCMode::FLD ||
            currentGroupControl == Group::ProductionCMode::NONE) {
            if (!group.productionGroupControlAvailable()) {
                // We cannot go any further up the hierarchy. This could
                // be the FIELD group, or any group for which this has
                // been set in GCONINJE or GCONPROD. If we are here
                // anyway, it is likely that the deck set inconsistent
                // requirements, such as GRUP control mode on a well with
                // no appropriate controls defined on any of its
                // containing groups. We will therefore use the wells' bhp
                // limit equation as a fallback.
                const auto& controls = well_ecl_.productionControls(summaryState);
                control_eq = bhp - controls.bhp_limit;
                return;
            } else {
                // Produce share of parents control
                const auto& parent = schedule.getGroup( group.parent(), current_step_ );
                efficiencyFactor *= group.getGroupEfficiencyFactor();
                getGroupProductionControl(parent, well_state, group_state, schedule, summaryState, bhp, rates, control_eq, efficiencyFactor);
                return;
            }
        }

        efficiencyFactor *= group.getGroupEfficiencyFactor();
        const auto& well = well_ecl_;
        const auto pu = phaseUsage();

        if (!group.isProductionGroup()) {
            // use bhp as control eq and let the updateControl code find a valid control
            const auto& controls = well.productionControls(summaryState);
            control_eq = bhp - controls.bhp_limit;
            return;
        }

        // If we are here, we are at the topmost group to be visited in the recursion.
        // This is the group containing the control we will check against.

        // Make conversion factors for RESV <-> surface rates.
        std::vector<double> resv_coeff(phaseUsage().num_phases, 1.0);
        rateConverter_.calcCoeff(0, pvtRegionIdx_, resv_coeff); // FIPNUM region 0 here, should use FIPNUM from WELSPECS.

        // gconsale may adjust the grat target.
        // the adjusted rates is send to the targetCalculator
        double gratTargetFromSales = 0.0;
        if (group_state.has_grat_sales_target(group.name()))
            gratTargetFromSales = group_state.grat_sales_target(group.name());

        WellGroupHelpers::TargetCalculator tcalc(currentGroupControl, pu, resv_coeff, gratTargetFromSales);
        WellGroupHelpers::FractionCalculator fcalc(schedule, well_state, group_state, current_step_, guide_rate_, tcalc.guideTargetMode(), pu, true, Phase::OIL);

        auto localFraction = [&](const std::string& child) {
            return fcalc.localFraction(child, "");
        };

        auto localReduction = [&](const std::string& group_name) {
            const std::vector<double>& groupTargetReductions = group_state.production_reduction_rates(group_name);
            return tcalc.calcModeRateFromRates(groupTargetReductions);
        };

        const double orig_target = tcalc.groupTarget(group.productionControls(summaryState));
        const auto chain = WellGroupHelpers::groupChainTopBot(name(), group.name(), schedule, current_step_);
        // Because 'name' is the last of the elements, and not an ancestor, we subtract one below.
        const size_t num_ancestors = chain.size() - 1;
        double target = orig_target;
        for (size_t ii = 0; ii < num_ancestors; ++ii) {
            if ((ii == 0) || guide_rate_->has(chain[ii])) {
                // Apply local reductions only at the control level
                // (top) and for levels where we have a specified
                // group guide rate.
                target -= localReduction(chain[ii]);
            }
            target *= localFraction(chain[ii+1]);
        }
        // Avoid negative target rates coming from too large local reductions.
        const double target_rate = std::max(0.0, target / efficiencyFactor);
        const auto current_rate = -tcalc.calcModeRateFromRates(rates); // Switch sign since 'rates' are negative for producers.
        control_eq = current_rate - target_rate;
    }




    template <typename TypeTag>
    void
    WellInterface<TypeTag>::
    updateWellStateRates(const Simulator& ebosSimulator,
                         WellState& well_state,
                         DeferredLogger& deferred_logger) const
    {
        // Check if the rates of this well only are single-phase, do nothing
        // if more than one nonzero rate.
        int nonzero_rate_index = -1;
        for (int p = 0; p < number_of_phases_; ++p) {
            if (well_state.wellRates()[index_of_well_ * number_of_phases_ + p] != 0.0) {
                if (nonzero_rate_index == -1) {
                    nonzero_rate_index = p;
                } else {
                    // More than one nonzero rate.
                    return;
                }
            }
        }
        if (nonzero_rate_index == -1) {
            // No nonzero rates.
            return;
        }

        // Calculate the rates that follow from the current primary variables.
        std::vector<double> well_q_s = computeCurrentWellRates(ebosSimulator, deferred_logger);

        // Set the currently-zero phase flows to be nonzero in proportion to well_q_s.
        const double initial_nonzero_rate = well_state.wellRates()[index_of_well_ * number_of_phases_ + nonzero_rate_index];
        const int comp_idx_nz = flowPhaseToEbosCompIdx(nonzero_rate_index);
        for (int p = 0; p < number_of_phases_; ++p) {
            if (p != nonzero_rate_index) {
                const int comp_idx = flowPhaseToEbosCompIdx(p);
                double& rate = well_state.wellRates()[index_of_well_ * number_of_phases_ + p];
                rate = (initial_nonzero_rate/well_q_s[comp_idx_nz]) * (well_q_s[comp_idx]);
            }
        }
    }


} // namespace Opm
