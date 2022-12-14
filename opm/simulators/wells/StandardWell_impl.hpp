/*
  Copyright 2017 SINTEF Digital, Mathematics and Cybernetics.
  Copyright 2017 Statoil ASA.
  Copyright 2016 - 2017 IRIS AS.

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

#include <opm/common/Exceptions.hpp>

#include <opm/simulators/utils/DeferredLoggingErrorHelpers.hpp>
#include <opm/simulators/wells/StandardWellAssemble.hpp>
#include <opm/simulators/wells/VFPHelpers.hpp>
#include <opm/simulators/wells/WellBhpThpCalculator.hpp>
#include <opm/simulators/wells/WellConvergence.hpp>

#include <algorithm>
#include <functional>
#include <numeric>

namespace Opm
{

    template<typename TypeTag>
    StandardWell<TypeTag>::
    StandardWell(const Well& well,
                 const ParallelWellInfo& pw_info,
                 const int time_step,
                 const ModelParameters& param,
                 const RateConverterType& rate_converter,
                 const int pvtRegionIdx,
                 const int num_components,
                 const int num_phases,
                 const int index_of_well,
                 const std::vector<PerforationData>& perf_data)
    : Base(well, pw_info, time_step, param, rate_converter, pvtRegionIdx, num_components, num_phases, index_of_well, perf_data)
    , StdWellEval(static_cast<const WellInterfaceIndices<FluidSystem,Indices,Scalar>&>(*this))
    , regularize_(false)
    {
        assert(this->num_components_ == numWellConservationEq);
    }





    template<typename TypeTag>
    void
    StandardWell<TypeTag>::
    init(const PhaseUsage* phase_usage_arg,
         const std::vector<double>& depth_arg,
         const double gravity_arg,
         const int num_cells,
         const std::vector< Scalar >& B_avg,
         const bool changed_to_open_this_step)
    {
        Base::init(phase_usage_arg, depth_arg, gravity_arg, num_cells, B_avg, changed_to_open_this_step);
        this->StdWellEval::init(this->perf_depth_, depth_arg, num_cells, Base::has_polymermw);
    }





    template<typename TypeTag>
    void StandardWell<TypeTag>::
    initPrimaryVariablesEvaluation()
    {
        this->primary_variables_.init();
    }





    template<typename TypeTag>
    void
    StandardWell<TypeTag>::
    computePerfRateEval(const IntensiveQuantities& intQuants,
                        const std::vector<EvalWell>& mob,
                        const EvalWell& bhp,
                        const double Tw,
                        const int perf,
                        const bool allow_cf,
                        std::vector<EvalWell>& cq_s,
                        double& perf_dis_gas_rate,
                        double& perf_vap_oil_rate,
                        double& perf_vap_wat_rate,
                        DeferredLogger& deferred_logger) const
    {
        const auto& fs = intQuants.fluidState();
        const EvalWell pressure = this->extendEval(this->getPerfCellPressure(fs));
        const EvalWell rs = this->extendEval(fs.Rs());
        const EvalWell rv = this->extendEval(fs.Rv());
        const EvalWell rvw = this->extendEval(fs.Rvw());

        std::vector<EvalWell> b_perfcells_dense(this->num_components_, EvalWell{this->primary_variables_.numWellEq() + Indices::numEq, 0.0});
        for (unsigned phaseIdx = 0; phaseIdx < FluidSystem::numPhases; ++phaseIdx) {
            if (!FluidSystem::phaseIsActive(phaseIdx)) {
                continue;
            }
            const unsigned compIdx = Indices::canonicalToActiveComponentIndex(FluidSystem::solventComponentIndex(phaseIdx));
            b_perfcells_dense[compIdx] =  this->extendEval(fs.invB(phaseIdx));
        }
        if constexpr (has_solvent) {
            b_perfcells_dense[Indices::contiSolventEqIdx] = this->extendEval(intQuants.solventInverseFormationVolumeFactor());
        }

        if constexpr (has_zFraction) {
            if (this->isInjector()) {
                const unsigned gasCompIdx = Indices::canonicalToActiveComponentIndex(FluidSystem::gasCompIdx);
                b_perfcells_dense[gasCompIdx] *= (1.0 - this->wsolvent());
                b_perfcells_dense[gasCompIdx] += this->wsolvent()*intQuants.zPureInvFormationVolumeFactor().value();
            }
        }

        EvalWell skin_pressure = EvalWell{this->primary_variables_.numWellEq() + Indices::numEq, 0.0};
        if (has_polymermw) {
            if (this->isInjector()) {
                const int pskin_index = Bhp + 1 + this->numPerfs() + perf;
                skin_pressure = this->primary_variables_.eval(pskin_index);
            }
        }

        // surface volume fraction of fluids within wellbore
        std::vector<EvalWell> cmix_s(this->numComponents(), EvalWell{this->primary_variables_.numWellEq() + Indices::numEq});
        for (int componentIdx = 0; componentIdx < this->numComponents(); ++componentIdx) {
            cmix_s[componentIdx] = this->primary_variables_.surfaceVolumeFraction(componentIdx);
        }

        computePerfRate(mob,
                        pressure,
                        bhp,
                        rs,
                        rv,
                        rvw,
                        b_perfcells_dense,
                        Tw,
                        perf,
                        allow_cf,
                        skin_pressure,
                        cmix_s,
                        cq_s,
                        perf_dis_gas_rate,
                        perf_vap_oil_rate,
                        perf_vap_wat_rate,
                        deferred_logger);
    }

    template<typename TypeTag>
    void
    StandardWell<TypeTag>::
    computePerfRateScalar(const IntensiveQuantities& intQuants,
                          const std::vector<Scalar>& mob,
                          const Scalar& bhp,
                          const double Tw,
                          const int perf,
                          const bool allow_cf,
                          std::vector<Scalar>& cq_s,
                          DeferredLogger& deferred_logger) const
    {
        const auto& fs = intQuants.fluidState();
        const Scalar pressure = this->getPerfCellPressure(fs).value();
        const Scalar rs = fs.Rs().value();
        const Scalar rv = fs.Rv().value();
        const Scalar rvw = fs.Rvw().value();
        std::vector<Scalar> b_perfcells_dense(this->num_components_, 0.0);
        for (unsigned phaseIdx = 0; phaseIdx < FluidSystem::numPhases; ++phaseIdx) {
            if (!FluidSystem::phaseIsActive(phaseIdx)) {
                continue;
            }
            const unsigned compIdx = Indices::canonicalToActiveComponentIndex(FluidSystem::solventComponentIndex(phaseIdx));
            b_perfcells_dense[compIdx] =  fs.invB(phaseIdx).value();
        }
        if constexpr (has_solvent) {
            b_perfcells_dense[Indices::contiSolventEqIdx] = intQuants.solventInverseFormationVolumeFactor().value();
        }

        if constexpr (has_zFraction) {
            if (this->isInjector()) {
                const unsigned gasCompIdx = Indices::canonicalToActiveComponentIndex(FluidSystem::gasCompIdx);
                b_perfcells_dense[gasCompIdx] *= (1.0 - this->wsolvent());
                b_perfcells_dense[gasCompIdx] += this->wsolvent()*intQuants.zPureInvFormationVolumeFactor().value();
            }
        }

        Scalar skin_pressure =0.0;
        if (has_polymermw) {
            if (this->isInjector()) {
                const int pskin_index = Bhp + 1 + this->numPerfs() + perf;
                skin_pressure = getValue(this->primary_variables_.eval(pskin_index));
            }
        }

        Scalar perf_dis_gas_rate = 0.0;
        Scalar perf_vap_oil_rate = 0.0;
        Scalar perf_vap_wat_rate = 0.0;

        // surface volume fraction of fluids within wellbore
        std::vector<Scalar> cmix_s(this->numComponents(), 0.0);
        for (int componentIdx = 0; componentIdx < this->numComponents(); ++componentIdx) {
            cmix_s[componentIdx] = getValue(this->primary_variables_.surfaceVolumeFraction(componentIdx));
        }

        computePerfRate(mob,
                        pressure,
                        bhp,
                        rs,
                        rv,
                        rvw,
                        b_perfcells_dense,
                        Tw,
                        perf,
                        allow_cf,
                        skin_pressure,
                        cmix_s,
                        cq_s,
                        perf_dis_gas_rate,
                        perf_vap_oil_rate,
                        perf_vap_wat_rate,
                        deferred_logger);
    }

    template<typename TypeTag>
    template<class Value>
    void
    StandardWell<TypeTag>::
    computePerfRate(const std::vector<Value>& mob,
                    const Value& pressure,
                    const Value& bhp,
                    const Value& rs,
                    const Value& rv,
                    const Value& rvw,
                    std::vector<Value>& b_perfcells_dense,
                    const double Tw,
                    const int perf,
                    const bool allow_cf,
                    const Value& skin_pressure,
                    const std::vector<Value>& cmix_s,
                    std::vector<Value>& cq_s,
                    double& perf_dis_gas_rate,
                    double& perf_vap_oil_rate,
                    double& perf_vap_wat_rate,
                    DeferredLogger& deferred_logger) const
    {
        // Pressure drawdown (also used to determine direction of flow)
        const Value well_pressure = bhp + this->connections_.pressure_diff(perf);
        Value drawdown = pressure - well_pressure;
        if (this->isInjector()) {
            drawdown += skin_pressure;
        }

        // producing perforations
        if ( drawdown > 0 )  {
            //Do nothing if crossflow is not allowed
            if (!allow_cf && this->isInjector()) {
                return;
            }

            // compute component volumetric rates at standard conditions
            for (int componentIdx = 0; componentIdx < this->numComponents(); ++componentIdx) {
                const Value cq_p = - Tw * (mob[componentIdx] * drawdown);
                cq_s[componentIdx] = b_perfcells_dense[componentIdx] * cq_p;
            }

            if (FluidSystem::phaseIsActive(FluidSystem::oilPhaseIdx) && FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx)) {
                const unsigned oilCompIdx = Indices::canonicalToActiveComponentIndex(FluidSystem::oilCompIdx);
                const unsigned gasCompIdx = Indices::canonicalToActiveComponentIndex(FluidSystem::gasCompIdx);
                const Value cq_sOil = cq_s[oilCompIdx];
                const Value cq_sGas = cq_s[gasCompIdx];
                const Value dis_gas = rs * cq_sOil;
                const Value vap_oil = rv * cq_sGas;

                cq_s[gasCompIdx] += dis_gas;
                cq_s[oilCompIdx] += vap_oil;

                // recording the perforation solution gas rate and solution oil rates
                if (this->isProducer()) {
                    perf_dis_gas_rate = getValue(dis_gas);
                    perf_vap_oil_rate = getValue(vap_oil);
                }
                if (FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx)) {
                    const unsigned waterCompIdx = Indices::canonicalToActiveComponentIndex(FluidSystem::waterCompIdx);
                    const Value vap_wat = rvw * cq_sGas;
                    cq_s[waterCompIdx] += vap_wat;
                    if (this->isProducer())
                        perf_vap_wat_rate = getValue(vap_wat);
                }
            }

        } else {
            //Do nothing if crossflow is not allowed
            if (!allow_cf && this->isProducer()) {
                return;
            }

            // Using total mobilities
            Value total_mob_dense = mob[0];
            for (int componentIdx = 1; componentIdx < this->numComponents(); ++componentIdx) {
                total_mob_dense += mob[componentIdx];
            }

            // injection perforations total volume rates
            const Value cqt_i = - Tw * (total_mob_dense * drawdown);

            // compute volume ratio between connection at standard conditions
            Value volumeRatio = bhp * 0.0; // initialize it with the correct type
;
            if (FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx)) {
                const unsigned waterCompIdx = Indices::canonicalToActiveComponentIndex(FluidSystem::waterCompIdx);
                volumeRatio += cmix_s[waterCompIdx] / b_perfcells_dense[waterCompIdx];
            }

            if constexpr (Indices::enableSolvent) {
                volumeRatio += cmix_s[Indices::contiSolventEqIdx] / b_perfcells_dense[Indices::contiSolventEqIdx];
            }

            if (FluidSystem::phaseIsActive(FluidSystem::oilPhaseIdx) && FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx)) {
                const unsigned oilCompIdx = Indices::canonicalToActiveComponentIndex(FluidSystem::oilCompIdx);
                const unsigned gasCompIdx = Indices::canonicalToActiveComponentIndex(FluidSystem::gasCompIdx);
                // Incorporate RS/RV factors if both oil and gas active
                const Value d = 1.0 - rv * rs;

                if (d <= 0.0) {
                    std::ostringstream sstr;
                    sstr << "Problematic d value " << d << " obtained for well " << this->name()
                         << " during computePerfRate calculations with rs " << rs
                         << ", rv " << rv << " and pressure " << pressure
                         << " obtaining d " << d
                         << " Continue as if no dissolution (rs = 0) and vaporization (rv = 0) "
                         << " for this connection.";
                    deferred_logger.debug(sstr.str());
                }
                const Value tmp_oil = d > 0.0? (cmix_s[oilCompIdx] - rv * cmix_s[gasCompIdx]) / d : cmix_s[oilCompIdx];
                volumeRatio += tmp_oil / b_perfcells_dense[oilCompIdx];

                const Value tmp_gas =  d > 0.0? (cmix_s[gasCompIdx] - rs * cmix_s[oilCompIdx]) / d : cmix_s[gasCompIdx];
                volumeRatio += tmp_gas / b_perfcells_dense[gasCompIdx];
            }
            else {
                if (FluidSystem::phaseIsActive(FluidSystem::oilPhaseIdx)) {
                    const unsigned oilCompIdx = Indices::canonicalToActiveComponentIndex(FluidSystem::oilCompIdx);
                    volumeRatio += cmix_s[oilCompIdx] / b_perfcells_dense[oilCompIdx];
                }
                if (FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx)) {
                    const unsigned gasCompIdx = Indices::canonicalToActiveComponentIndex(FluidSystem::gasCompIdx);
                    volumeRatio += cmix_s[gasCompIdx] / b_perfcells_dense[gasCompIdx];
                }
            }

            // injecting connections total volumerates at standard conditions
            Value cqt_is = cqt_i/volumeRatio;
            for (int componentIdx = 0; componentIdx < this->numComponents(); ++componentIdx) {
                cq_s[componentIdx] = cmix_s[componentIdx] * cqt_is;
            }

            // calculating the perforation solution gas rate and solution oil rates
            if (this->isProducer()) {
                if (FluidSystem::phaseIsActive(FluidSystem::oilPhaseIdx) && FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx)) {
                    const unsigned oilCompIdx = Indices::canonicalToActiveComponentIndex(FluidSystem::oilCompIdx);
                    const unsigned gasCompIdx = Indices::canonicalToActiveComponentIndex(FluidSystem::gasCompIdx);
                    // TODO: the formulations here remain to be tested with cases with strong crossflow through production wells
                    // s means standard condition, r means reservoir condition
                    // q_os = q_or * b_o + rv * q_gr * b_g
                    // q_gs = q_gr * b_g + rs * q_or * b_o
                    // d = 1.0 - rs * rv
                    // q_or = 1 / (b_o * d) * (q_os - rv * q_gs)
                    // q_gr = 1 / (b_g * d) * (q_gs - rs * q_os)

                    const double d = 1.0 - getValue(rv) * getValue(rs);

                    if (d <= 0.0) {
                        std::ostringstream sstr;
                        sstr << "Problematic d value " << d << " obtained for well " << this->name()
                             << " during computePerfRate calculations with rs " << rs
                            << ", rv " << rv << " and pressure " << pressure
                            << " obtaining d " << d
                            << " Continue as if no dissolution (rs = 0) and vaporization (rv = 0) "
                            << " for this connection.";
                        deferred_logger.debug(sstr.str());
                    } else {
                        // vaporized oil into gas
                        // rv * q_gr * b_g = rv * (q_gs - rs * q_os) / d
                        perf_vap_oil_rate = getValue(rv) * (getValue(cq_s[gasCompIdx]) - getValue(rs) * getValue(cq_s[oilCompIdx])) / d;
                        // dissolved of gas in oil
                        // rs * q_or * b_o = rs * (q_os - rv * q_gs) / d
                        perf_dis_gas_rate = getValue(rs) * (getValue(cq_s[oilCompIdx]) - getValue(rv) * getValue(cq_s[gasCompIdx])) / d;
                    }
                    if (FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx)) {
                        // q_ws = q_wr * b_w + rvw * q_gr * b_g
                        // q_wr = 1 / b_w * (q_ws - rvw * q_gr * b_g) = 1 / b_w * (q_ws - rvw * 1 / d  (q_gs - rs * q_os))
                        // vaporized water in gas
                        // rvw * q_gr * b_g = q_ws -q_wr *b_w = rvw * (q_gs -rs *q_os) / d
                        perf_vap_wat_rate = getValue(rvw) * (getValue(cq_s[gasCompIdx]) - getValue(rs) * getValue(cq_s[oilCompIdx])) / d;
                    }
                }
                else if (FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx) && FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx)) {
                    //no oil
                    const unsigned gasCompIdx = Indices::canonicalToActiveComponentIndex(FluidSystem::gasCompIdx);
                    perf_vap_wat_rate = getValue(rvw) * getValue(cq_s[gasCompIdx]);
                }
            }
        }
    }


    template<typename TypeTag>
    void
    StandardWell<TypeTag>::
    assembleWellEqWithoutIteration(const Simulator& ebosSimulator,
                                   const double dt,
                                   const Well::InjectionControls& /*inj_controls*/,
                                   const Well::ProductionControls& /*prod_controls*/,
                                   WellState& well_state,
                                   const GroupState& group_state,
                                   DeferredLogger& deferred_logger)
    {
        // TODO: only_wells should be put back to save some computation
        // for example, the matrices B C does not need to update if only_wells
        if (!this->isOperableAndSolvable() && !this->wellIsStopped()) return;

        // clear all entries
        this->linSys_.clear();

        assembleWellEqWithoutIterationImpl(ebosSimulator, dt, well_state, group_state, deferred_logger);
    }




    template<typename TypeTag>
    void
    StandardWell<TypeTag>::
    assembleWellEqWithoutIterationImpl(const Simulator& ebosSimulator,
                                       const double dt,
                                       WellState& well_state,
                                       const GroupState& group_state,
                                       DeferredLogger& deferred_logger)
    {
        // try to regularize equation if the well does not converge
        const Scalar regularization_factor =  this->regularize_? this->param_.regularization_factor_wells_ : 1.0;
        const double volume = 0.1 * unit::cubic(unit::feet) * regularization_factor;

        auto& ws = well_state.well(this->index_of_well_);

        ws.vaporized_oil_rate = 0;
        ws.dissolved_gas_rate = 0;
        ws.vaporized_wat_rate = 0;

        const int np = this->number_of_phases_;

        std::vector<RateVector> connectionRates = this->connectionRates_; // Copy to get right size.
        auto& perf_data = ws.perf_data;
        auto& perf_rates = perf_data.phase_rates;
        for (int perf = 0; perf < this->number_of_perforations_; ++perf) {
            // Calculate perforation quantities.
            std::vector<EvalWell> cq_s(this->num_components_, {this->primary_variables_.numWellEq() + Indices::numEq, 0.0});
            EvalWell water_flux_s{this->primary_variables_.numWellEq() + Indices::numEq, 0.0};
            EvalWell cq_s_zfrac_effective{this->primary_variables_.numWellEq() + Indices::numEq, 0.0};
            calculateSinglePerf(ebosSimulator, perf, well_state, connectionRates, cq_s, water_flux_s, cq_s_zfrac_effective, deferred_logger);

            // Equation assembly for this perforation.
            if constexpr (has_polymer && Base::has_polymermw) {
                if (this->isInjector()) {
                    handleInjectivityEquations(ebosSimulator, well_state, perf, water_flux_s, deferred_logger);
                }
            }
            const int cell_idx = this->well_cells_[perf];
            for (int componentIdx = 0; componentIdx < this->num_components_; ++componentIdx) {
                // the cq_s entering mass balance equations need to consider the efficiency factors.
                const EvalWell cq_s_effective = cq_s[componentIdx] * this->well_efficiency_factor_;

                connectionRates[perf][componentIdx] = Base::restrictEval(cq_s_effective);

                StandardWellAssemble<FluidSystem,Indices,Scalar>(*this).
                    assemblePerforationEq(cq_s_effective,
                                          componentIdx,
                                          cell_idx,
                                          this->primary_variables_.numWellEq(),
                                          this->linSys_);

                // Store the perforation phase flux for later usage.
                if (has_solvent && componentIdx == Indices::contiSolventEqIdx) {
                    auto& perf_rate_solvent = perf_data.solvent_rates;
                    perf_rate_solvent[perf] = cq_s[componentIdx].value();
                } else {
                    perf_rates[perf*np + this->ebosCompIdxToFlowCompIdx(componentIdx)] = cq_s[componentIdx].value();
                }
            }

            if constexpr (has_zFraction) {
                StandardWellAssemble<FluidSystem,Indices,Scalar>(*this).
                    assembleZFracEq(cq_s_zfrac_effective,
                                    cell_idx,
                                    this->primary_variables_.numWellEq(),
                                    this->linSys_);
            }
        }
        // Update the connection
        this->connectionRates_ = connectionRates;

        // Accumulate dissolved gas and vaporized oil flow rates across all
        // ranks sharing this well (this->index_of_well_).
        {
            const auto& comm = this->parallel_well_info_.communication();
            ws.dissolved_gas_rate = comm.sum(ws.dissolved_gas_rate);
            ws.vaporized_oil_rate = comm.sum(ws.vaporized_oil_rate);
            ws.vaporized_wat_rate = comm.sum(ws.vaporized_wat_rate);
        }

        // accumulate resWell_ and duneD_ in parallel to get effects of all perforations (might be distributed)
        this->linSys_.sumDistributed(this->parallel_well_info_.communication());

        // add vol * dF/dt + Q to the well equations;
        for (int componentIdx = 0; componentIdx < numWellConservationEq; ++componentIdx) {
            // TODO: following the development in MSW, we need to convert the volume of the wellbore to be surface volume
            // since all the rates are under surface condition
            EvalWell resWell_loc(this->primary_variables_.numWellEq() + Indices::numEq, 0.0);
            if (FluidSystem::numActivePhases() > 1) {
                assert(dt > 0);
                resWell_loc += (this->primary_variables_.surfaceVolumeFraction(componentIdx) -
                                this->F0_[componentIdx]) * volume / dt;
            }
            resWell_loc -= this->primary_variables_.getQs(componentIdx) * this->well_efficiency_factor_;
            StandardWellAssemble<FluidSystem,Indices,Scalar>(*this).
                assembleSourceEq(resWell_loc,
                                 componentIdx,
                                 this->primary_variables_.numWellEq(),
                                 this->linSys_);
        }

        const auto& summaryState = ebosSimulator.vanguard().summaryState();
        const Schedule& schedule = ebosSimulator.vanguard().schedule();
        StandardWellAssemble<FluidSystem,Indices,Scalar>(*this).
            assembleControlEq(well_state, group_state,
                              schedule, summaryState,
                              this->primary_variables_,
                              this->connections_.rho(),
                              this->linSys_,
                              deferred_logger);


        // do the local inversion of D.
        try {
            this->linSys_.invert();
        } catch( ... ) {
            OPM_DEFLOG_THROW(NumericalProblem, "Error when inverting local well equations for well " + name(), deferred_logger);
        }
    }




    template<typename TypeTag>
    void
    StandardWell<TypeTag>::
    calculateSinglePerf(const Simulator& ebosSimulator,
                        const int perf,
                        WellState& well_state,
                        std::vector<RateVector>& connectionRates,
                        std::vector<EvalWell>& cq_s,
                        EvalWell& water_flux_s,
                        EvalWell& cq_s_zfrac_effective,
                        DeferredLogger& deferred_logger) const
    {
        const bool allow_cf = this->getAllowCrossFlow() || openCrossFlowAvoidSingularity(ebosSimulator);
        const EvalWell& bhp = this->primary_variables_.eval(Bhp);
        const int cell_idx = this->well_cells_[perf];
        const auto& intQuants = *(ebosSimulator.model().cachedIntensiveQuantities(cell_idx, /*timeIdx=*/ 0));
        std::vector<EvalWell> mob(this->num_components_, {this->primary_variables_.numWellEq() + Indices::numEq, 0.});
        getMobilityEval(ebosSimulator, perf, mob, deferred_logger);

        double perf_dis_gas_rate = 0.;
        double perf_vap_oil_rate = 0.;
        double perf_vap_wat_rate = 0.;
        double trans_mult = ebosSimulator.problem().template rockCompTransMultiplier<double>(intQuants,  cell_idx);
        const double Tw = this->well_index_[perf] * trans_mult;
        computePerfRateEval(intQuants, mob, bhp, Tw, perf, allow_cf,
                            cq_s, perf_dis_gas_rate, perf_vap_oil_rate, perf_vap_wat_rate, deferred_logger);

        auto& ws = well_state.well(this->index_of_well_);
        auto& perf_data = ws.perf_data;
        if constexpr (has_polymer && Base::has_polymermw) {
            if (this->isInjector()) {
                // Store the original water flux computed from the reservoir quantities.
                // It will be required to assemble the injectivity equations.
                const unsigned water_comp_idx = Indices::canonicalToActiveComponentIndex(FluidSystem::waterCompIdx);
                water_flux_s = cq_s[water_comp_idx];
                // Modify the water flux for the rest of this function to depend directly on the
                // local water velocity primary variable.
                handleInjectivityRate(ebosSimulator, perf, cq_s);
            }
        }

        // updating the solution gas rate and solution oil rate
        if (this->isProducer()) {
            ws.dissolved_gas_rate += perf_dis_gas_rate;
            ws.vaporized_oil_rate += perf_vap_oil_rate;
            ws.vaporized_wat_rate += perf_vap_wat_rate;
        }

        if constexpr (has_energy) {
            connectionRates[perf][Indices::contiEnergyEqIdx] = 0.0;
        }

        if constexpr (has_energy) {

            auto fs = intQuants.fluidState();
            for (unsigned phaseIdx = 0; phaseIdx < FluidSystem::numPhases; ++phaseIdx) {
                if (!FluidSystem::phaseIsActive(phaseIdx)) {
                    continue;
                }

                // convert to reservoir conditions
                EvalWell cq_r_thermal(this->primary_variables_.numWellEq() + Indices::numEq, 0.);
                const unsigned activeCompIdx = Indices::canonicalToActiveComponentIndex(FluidSystem::solventComponentIndex(phaseIdx));
                const bool both_oil_gas = FluidSystem::phaseIsActive(FluidSystem::oilPhaseIdx) && FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx);
                if ( !both_oil_gas || FluidSystem::waterPhaseIdx == phaseIdx ) {
                    cq_r_thermal = cq_s[activeCompIdx] / this->extendEval(fs.invB(phaseIdx));
                } else {
                    // remove dissolved gas and vapporized oil
                    const unsigned oilCompIdx = Indices::canonicalToActiveComponentIndex(FluidSystem::oilCompIdx);
                    const unsigned gasCompIdx = Indices::canonicalToActiveComponentIndex(FluidSystem::gasCompIdx);
                    // q_os = q_or * b_o + rv * q_gr * b_g
                    // q_gs = q_gr * g_g + rs * q_or * b_o
                    // q_gr = 1 / (b_g * d) * (q_gs - rs * q_os)
                    // d = 1.0 - rs * rv
                    const EvalWell d = this->extendEval(1.0 - fs.Rv() * fs.Rs());
                    if (d <= 0.0) {
                        std::ostringstream sstr;
                        sstr << "Problematic d value " << d << " obtained for well " << this->name()
                            << " during calculateSinglePerf with rs " << fs.Rs()
                            << ", rv " << fs.Rv()
                            << " obtaining d " << d
                            << " Continue as if no dissolution (rs = 0) and vaporization (rv = 0) "
                            << " for this connection.";
                        deferred_logger.debug(sstr.str());
                        cq_r_thermal = cq_s[activeCompIdx] / this->extendEval(fs.invB(phaseIdx));
                    } else {
                        if(FluidSystem::gasPhaseIdx == phaseIdx) {
                            cq_r_thermal = (cq_s[gasCompIdx] - this->extendEval(fs.Rs()) * cq_s[oilCompIdx]) / (d * this->extendEval(fs.invB(phaseIdx)) );
                        } else if(FluidSystem::oilPhaseIdx == phaseIdx) {
                            // q_or = 1 / (b_o * d) * (q_os - rv * q_gs)
                            cq_r_thermal = (cq_s[oilCompIdx] - this->extendEval(fs.Rv()) * cq_s[gasCompIdx]) / (d * this->extendEval(fs.invB(phaseIdx)) );
                        }
                    }
                }

                // change temperature for injecting fluids
                if (this->isInjector() && cq_s[activeCompIdx] > 0.0){
                    // only handles single phase injection now
                    assert(this->well_ecl_.injectorType() != InjectorType::MULTI);
                    fs.setTemperature(this->well_ecl_.temperature());
                    typedef typename std::decay<decltype(fs)>::type::Scalar FsScalar;
                    typename FluidSystem::template ParameterCache<FsScalar> paramCache;
                    const unsigned pvtRegionIdx = intQuants.pvtRegionIndex();
                    paramCache.setRegionIndex(pvtRegionIdx);
                    paramCache.setMaxOilSat(ebosSimulator.problem().maxOilSaturation(cell_idx));
                    paramCache.updatePhase(fs, phaseIdx);

                    const auto& rho = FluidSystem::density(fs, paramCache, phaseIdx);
                    fs.setDensity(phaseIdx, rho);
                    const auto& h = FluidSystem::enthalpy(fs, paramCache, phaseIdx);
                    fs.setEnthalpy(phaseIdx, h);
                    cq_r_thermal *= this->extendEval(fs.enthalpy(phaseIdx)) * this->extendEval(fs.density(phaseIdx));
                    connectionRates[perf][Indices::contiEnergyEqIdx] += getValue(cq_r_thermal);
                } else {
                    // compute the thermal flux
                    cq_r_thermal *= this->extendEval(fs.enthalpy(phaseIdx)) * this->extendEval(fs.density(phaseIdx));
                    connectionRates[perf][Indices::contiEnergyEqIdx] += Base::restrictEval(cq_r_thermal);
                }
            }
        }

        if constexpr (has_polymer) {
            // TODO: the application of well efficiency factor has not been tested with an example yet
            const unsigned waterCompIdx = Indices::canonicalToActiveComponentIndex(FluidSystem::waterCompIdx);
            EvalWell cq_s_poly = cq_s[waterCompIdx];
            if (this->isInjector()) {
                cq_s_poly *= this->wpolymer();
            } else {
                cq_s_poly *= this->extendEval(intQuants.polymerConcentration() * intQuants.polymerViscosityCorrection());
            }
            // Note. Efficiency factor is handled in the output layer
            auto& perf_rate_polymer = perf_data.polymer_rates;
            perf_rate_polymer[perf] = cq_s_poly.value();

            cq_s_poly *= this->well_efficiency_factor_;
            connectionRates[perf][Indices::contiPolymerEqIdx] = Base::restrictEval(cq_s_poly);

            if constexpr (Base::has_polymermw) {
                updateConnectionRatePolyMW(cq_s_poly, intQuants, well_state, perf, connectionRates, deferred_logger);
            }
        }

        if constexpr (has_foam) {
            // TODO: the application of well efficiency factor has not been tested with an example yet
            const unsigned gasCompIdx = Indices::canonicalToActiveComponentIndex(FluidSystem::gasCompIdx);
            EvalWell cq_s_foam = cq_s[gasCompIdx] * this->well_efficiency_factor_;
            if (this->isInjector()) {
                cq_s_foam *= this->wfoam();
            } else {
                cq_s_foam *= this->extendEval(intQuants.foamConcentration());
            }
            connectionRates[perf][Indices::contiFoamEqIdx] = Base::restrictEval(cq_s_foam);
        }

        if constexpr (has_zFraction) {
            // TODO: the application of well efficiency factor has not been tested with an example yet
            const unsigned gasCompIdx = Indices::canonicalToActiveComponentIndex(FluidSystem::gasCompIdx);
            cq_s_zfrac_effective = cq_s[gasCompIdx];
            if (this->isInjector()) {
                cq_s_zfrac_effective *= this->wsolvent();
            } else if (cq_s_zfrac_effective.value() != 0.0) {
                const double dis_gas_frac = perf_dis_gas_rate / cq_s_zfrac_effective.value();
                cq_s_zfrac_effective *= this->extendEval(dis_gas_frac*intQuants.xVolume() + (1.0-dis_gas_frac)*intQuants.yVolume());
            }
            auto& perf_rate_solvent = perf_data.solvent_rates;
            perf_rate_solvent[perf] = cq_s_zfrac_effective.value();

            cq_s_zfrac_effective *= this->well_efficiency_factor_;
            connectionRates[perf][Indices::contiZfracEqIdx] = Base::restrictEval(cq_s_zfrac_effective);
        }

        if constexpr (has_brine) {
            // TODO: the application of well efficiency factor has not been tested with an example yet
            const unsigned waterCompIdx = Indices::canonicalToActiveComponentIndex(FluidSystem::waterCompIdx);
            // Correction salt rate; evaporated water does not contain salt 
            EvalWell cq_s_sm = cq_s[waterCompIdx] - perf_vap_wat_rate;
            if (this->isInjector()) {
                cq_s_sm *= this->wsalt();
            } else {
                cq_s_sm *= this->extendEval(intQuants.fluidState().saltConcentration());
            }
            // Note. Efficiency factor is handled in the output layer
            auto& perf_rate_brine = perf_data.brine_rates;
            perf_rate_brine[perf] = cq_s_sm.value();

            cq_s_sm *= this->well_efficiency_factor_;
            connectionRates[perf][Indices::contiBrineEqIdx] = Base::restrictEval(cq_s_sm);
        }

        if constexpr (has_micp) {
            const unsigned waterCompIdx = Indices::canonicalToActiveComponentIndex(FluidSystem::waterCompIdx);
            EvalWell cq_s_microbe = cq_s[waterCompIdx];
            if (this->isInjector()) {
                cq_s_microbe *= this->wmicrobes();
            } else {
                cq_s_microbe *= this->extendEval(intQuants.microbialConcentration());
            }
            connectionRates[perf][Indices::contiMicrobialEqIdx] = Base::restrictEval(cq_s_microbe);
            EvalWell cq_s_oxygen = cq_s[waterCompIdx];
            if (this->isInjector()) {
                cq_s_oxygen *= this->woxygen();
            } else {
                cq_s_oxygen *= this->extendEval(intQuants.oxygenConcentration());
            }
            connectionRates[perf][Indices::contiOxygenEqIdx] = Base::restrictEval(cq_s_oxygen);
            EvalWell cq_s_urea = cq_s[waterCompIdx];
            if (this->isInjector()) {
                cq_s_urea *= this->wurea();
            } else {
                cq_s_urea *= this->extendEval(intQuants.ureaConcentration());
            }
            connectionRates[perf][Indices::contiUreaEqIdx] = Base::restrictEval(cq_s_urea);
        }

        // Store the perforation pressure for later usage.
        perf_data.pressure[perf] = ws.bhp + this->connections_.pressure_diff(perf);
    }




    template<typename TypeTag>
    void
    StandardWell<TypeTag>::
    getMobilityEval(const Simulator& ebosSimulator,
                    const int perf,
                    std::vector<EvalWell>& mob,
                    DeferredLogger& deferred_logger) const
    {
        const int cell_idx = this->well_cells_[perf];
        assert (int(mob.size()) == this->num_components_);
        const auto& intQuants = *(ebosSimulator.model().cachedIntensiveQuantities(cell_idx, /*timeIdx=*/0));
        const auto& materialLawManager = ebosSimulator.problem().materialLawManager();

        // either use mobility of the perforation cell or calcualte its own
        // based on passing the saturation table index
        const int satid = this->saturation_table_number_[perf] - 1;
        const int satid_elem = materialLawManager->satnumRegionIdx(cell_idx);
        if( satid == satid_elem ) { // the same saturation number is used. i.e. just use the mobilty from the cell

            for (unsigned phaseIdx = 0; phaseIdx < FluidSystem::numPhases; ++phaseIdx) {
                if (!FluidSystem::phaseIsActive(phaseIdx)) {
                    continue;
                }

                const unsigned activeCompIdx = Indices::canonicalToActiveComponentIndex(FluidSystem::solventComponentIndex(phaseIdx));
                mob[activeCompIdx] = this->extendEval(intQuants.mobility(phaseIdx));
            }
            if (has_solvent) {
                mob[Indices::contiSolventEqIdx] = this->extendEval(intQuants.solventMobility());
            }
        } else {

            const auto& paramsCell = materialLawManager->connectionMaterialLawParams(satid, cell_idx);
            std::array<Eval,3> relativePerms = { 0.0, 0.0, 0.0 };
            MaterialLaw::relativePermeabilities(relativePerms, paramsCell, intQuants.fluidState());

            // reset the satnumvalue back to original
            materialLawManager->connectionMaterialLawParams(satid_elem, cell_idx);

            // compute the mobility
            for (unsigned phaseIdx = 0; phaseIdx < FluidSystem::numPhases; ++phaseIdx) {
                if (!FluidSystem::phaseIsActive(phaseIdx)) {
                    continue;
                }

                const unsigned activeCompIdx = Indices::canonicalToActiveComponentIndex(FluidSystem::solventComponentIndex(phaseIdx));
                mob[activeCompIdx] = this->extendEval(relativePerms[phaseIdx] / intQuants.fluidState().viscosity(phaseIdx));
            }

            // this may not work if viscosity and relperms has been modified?
            if constexpr (has_solvent) {
                OPM_DEFLOG_THROW(std::runtime_error, "individual mobility for wells does not work in combination with solvent", deferred_logger);
            }
        }

        // modify the water mobility if polymer is present
        if constexpr (has_polymer) {
            if (!FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx)) {
                OPM_DEFLOG_THROW(std::runtime_error, "Water is required when polymer is active", deferred_logger);
            }

            // for the cases related to polymer molecular weight, we assume fully mixing
            // as a result, the polymer and water share the same viscosity
            if constexpr (!Base::has_polymermw) {
                updateWaterMobilityWithPolymer(ebosSimulator, perf, mob, deferred_logger);
            }
        }
    }

    template<typename TypeTag>
    void
    StandardWell<TypeTag>::
    getMobilityScalar(const Simulator& ebosSimulator,
                      const int perf,
                      std::vector<Scalar>& mob,
                      DeferredLogger& deferred_logger) const
    {
        const int cell_idx = this->well_cells_[perf];
        assert (int(mob.size()) == this->num_components_);
        const auto& intQuants = *(ebosSimulator.model().cachedIntensiveQuantities(cell_idx, /*timeIdx=*/0));
        const auto& materialLawManager = ebosSimulator.problem().materialLawManager();

        // either use mobility of the perforation cell or calcualte its own
        // based on passing the saturation table index
        const int satid = this->saturation_table_number_[perf] - 1;
        const int satid_elem = materialLawManager->satnumRegionIdx(cell_idx);
        if( satid == satid_elem ) { // the same saturation number is used. i.e. just use the mobilty from the cell

            for (unsigned phaseIdx = 0; phaseIdx < FluidSystem::numPhases; ++phaseIdx) {
                if (!FluidSystem::phaseIsActive(phaseIdx)) {
                    continue;
                }

                const unsigned activeCompIdx = Indices::canonicalToActiveComponentIndex(FluidSystem::solventComponentIndex(phaseIdx));
                mob[activeCompIdx] = getValue(intQuants.mobility(phaseIdx));
            }
            if (has_solvent) {
                mob[Indices::contiSolventEqIdx] = getValue(intQuants.solventMobility());
            }
        } else {

            const auto& paramsCell = materialLawManager->connectionMaterialLawParams(satid, cell_idx);
            std::array<Eval,3> relativePerms = { 0.0, 0.0, 0.0 };
            MaterialLaw::relativePermeabilities(relativePerms, paramsCell, intQuants.fluidState());

            // reset the satnumvalue back to original
            materialLawManager->connectionMaterialLawParams(satid_elem, cell_idx);

            // compute the mobility
            for (unsigned phaseIdx = 0; phaseIdx < FluidSystem::numPhases; ++phaseIdx) {
                if (!FluidSystem::phaseIsActive(phaseIdx)) {
                    continue;
                }

                const unsigned activeCompIdx = Indices::canonicalToActiveComponentIndex(FluidSystem::solventComponentIndex(phaseIdx));
                mob[activeCompIdx] = getValue(relativePerms[phaseIdx]) / getValue(intQuants.fluidState().viscosity(phaseIdx));
            }

            // this may not work if viscosity and relperms has been modified?
            if constexpr (has_solvent) {
                OPM_DEFLOG_THROW(std::runtime_error, "individual mobility for wells does not work in combination with solvent", deferred_logger);
            }
        }

        // modify the water mobility if polymer is present
        if constexpr (has_polymer) {
            if (!FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx)) {
                OPM_DEFLOG_THROW(std::runtime_error, "Water is required when polymer is active", deferred_logger);
            }

            // for the cases related to polymer molecular weight, we assume fully mixing
            // as a result, the polymer and water share the same viscosity
            if constexpr (!Base::has_polymermw) {
                std::vector<EvalWell> mob_eval(this->num_components_, {this->primary_variables_.numWellEq() + Indices::numEq, 0.});
                updateWaterMobilityWithPolymer(ebosSimulator, perf, mob_eval, deferred_logger);
                for (size_t i = 0; i < mob.size(); ++i) {
                    mob[i] = getValue(mob_eval[i]);
                }
            }
        }
    }



    template<typename TypeTag>
    void
    StandardWell<TypeTag>::
    updateWellState(const BVectorWell& dwells,
                    WellState& well_state,
                    DeferredLogger& deferred_logger)
    {
        if (!this->isOperableAndSolvable() && !this->wellIsStopped()) return;

        updatePrimaryVariablesNewton(dwells, well_state, deferred_logger);

        updateWellStateFromPrimaryVariables(well_state, deferred_logger);
        Base::calculateReservoirRates(well_state.well(this->index_of_well_));
    }





    template<typename TypeTag>
    void
    StandardWell<TypeTag>::
    updatePrimaryVariablesNewton(const BVectorWell& dwells,
                                 const WellState& /* well_state */,
                                 DeferredLogger& deferred_logger)
    {
        const double dFLimit = this->param_.dwell_fraction_max_;
        const double dBHPLimit = this->param_.dbhp_max_rel_;
        this->primary_variables_.updateNewton(dwells, dFLimit, dBHPLimit);

        // for the water velocity and skin pressure
        if constexpr (Base::has_polymermw) {
            this->primary_variables_.updateNewtonPolyMW(dwells);
        }

        this->primary_variables_.checkFinite(deferred_logger);
    }





    template<typename TypeTag>
    void
    StandardWell<TypeTag>::
    updateWellStateFromPrimaryVariables(WellState& well_state, DeferredLogger& deferred_logger) const
    {
        this->StdWellEval::updateWellStateFromPrimaryVariables(well_state, deferred_logger);

        // other primary variables related to polymer injectivity study
        if constexpr (Base::has_polymermw) {
            this->primary_variables_.copyToWellStatePolyMW(well_state);
        }
    }





    template<typename TypeTag>
    void
    StandardWell<TypeTag>::
    updateIPR(const Simulator& ebos_simulator, DeferredLogger& deferred_logger) const
    {
        // TODO: not handling solvent related here for now

        // initialize all the values to be zero to begin with
        std::fill(this->ipr_a_.begin(), this->ipr_a_.end(), 0.);
        std::fill(this->ipr_b_.begin(), this->ipr_b_.end(), 0.);

        for (int perf = 0; perf < this->number_of_perforations_; ++perf) {
            std::vector<Scalar> mob(this->num_components_, 0.0);
            getMobilityScalar(ebos_simulator, perf, mob, deferred_logger);

            const int cell_idx = this->well_cells_[perf];
            const auto& int_quantities = *(ebos_simulator.model().cachedIntensiveQuantities(cell_idx, /*timeIdx=*/ 0));
            const auto& fs = int_quantities.fluidState();
            // the pressure of the reservoir grid block the well connection is in
            double p_r = this->getPerfCellPressure(fs).value();

            // calculating the b for the connection
            std::vector<double> b_perf(this->num_components_);
            for (size_t phase = 0; phase < FluidSystem::numPhases; ++phase) {
                if (!FluidSystem::phaseIsActive(phase)) {
                    continue;
                }
                const unsigned comp_idx = Indices::canonicalToActiveComponentIndex(FluidSystem::solventComponentIndex(phase));
                b_perf[comp_idx] = fs.invB(phase).value();
            }
            if constexpr (has_solvent) {
                b_perf[Indices::contiSolventEqIdx] = int_quantities.solventInverseFormationVolumeFactor().value();
            }

            // the pressure difference between the connection and BHP
            const double h_perf = this->connections_.pressure_diff(perf);
            const double pressure_diff = p_r - h_perf;

            // Let us add a check, since the pressure is calculated based on zero value BHP
            // it should not be negative anyway. If it is negative, we might need to re-formulate
            // to taking into consideration the crossflow here.
            if ( (this->isProducer() && pressure_diff < 0.) || (this->isInjector() && pressure_diff > 0.) ) {
                deferred_logger.debug("CROSSFLOW_IPR",
                                "cross flow found when updateIPR for well " + name()
                                + " . The connection is ignored in IPR calculations");
                // we ignore these connections for now
                continue;
            }

            // the well index associated with the connection
            const double tw_perf = this->well_index_[perf]*ebos_simulator.problem().template rockCompTransMultiplier<double>(int_quantities, cell_idx);

            std::vector<double> ipr_a_perf(this->ipr_a_.size());
            std::vector<double> ipr_b_perf(this->ipr_b_.size());
            for (int comp_idx = 0; comp_idx < this->num_components_; ++comp_idx) {
                const double tw_mob = tw_perf * mob[comp_idx] * b_perf[comp_idx];
                ipr_a_perf[comp_idx] += tw_mob * pressure_diff;
                ipr_b_perf[comp_idx] += tw_mob;
            }

            // we need to handle the rs and rv when both oil and gas are present
            if (FluidSystem::phaseIsActive(FluidSystem::oilPhaseIdx) && FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx)) {
                const unsigned oil_comp_idx = Indices::canonicalToActiveComponentIndex(FluidSystem::oilCompIdx);
                const unsigned gas_comp_idx = Indices::canonicalToActiveComponentIndex(FluidSystem::gasCompIdx);
                const double rs = (fs.Rs()).value();
                const double rv = (fs.Rv()).value();

                const double dis_gas_a = rs * ipr_a_perf[oil_comp_idx];
                const double vap_oil_a = rv * ipr_a_perf[gas_comp_idx];

                ipr_a_perf[gas_comp_idx] += dis_gas_a;
                ipr_a_perf[oil_comp_idx] += vap_oil_a;

                const double dis_gas_b = rs * ipr_b_perf[oil_comp_idx];
                const double vap_oil_b = rv * ipr_b_perf[gas_comp_idx];

                ipr_b_perf[gas_comp_idx] += dis_gas_b;
                ipr_b_perf[oil_comp_idx] += vap_oil_b;
            }

            for (size_t comp_idx = 0; comp_idx < ipr_a_perf.size(); ++comp_idx) {
                this->ipr_a_[comp_idx] += ipr_a_perf[comp_idx];
                this->ipr_b_[comp_idx] += ipr_b_perf[comp_idx];
            }
        }
        this->parallel_well_info_.communication().sum(this->ipr_a_.data(), this->ipr_a_.size());
        this->parallel_well_info_.communication().sum(this->ipr_b_.data(), this->ipr_b_.size());
    }


    template<typename TypeTag>
    void
    StandardWell<TypeTag>::
    checkOperabilityUnderBHPLimit(const WellState& well_state, const Simulator& ebos_simulator, DeferredLogger& deferred_logger)
    {
        const auto& summaryState = ebos_simulator.vanguard().summaryState();
        const double bhp_limit = WellBhpThpCalculator(*this).mostStrictBhpFromBhpLimits(summaryState);
        // Crude but works: default is one atmosphere.
        // TODO: a better way to detect whether the BHP is defaulted or not
        const bool bhp_limit_not_defaulted = bhp_limit > 1.5 * unit::barsa;
        if ( bhp_limit_not_defaulted || !this->wellHasTHPConstraints(summaryState) ) {
            // if the BHP limit is not defaulted or the well does not have a THP limit
            // we need to check the BHP limit
            double total_ipr_mass_rate = 0.0;
            for (unsigned phaseIdx = 0; phaseIdx < FluidSystem::numPhases; ++phaseIdx)
            {
                if (!FluidSystem::phaseIsActive(phaseIdx)) {
                    continue;
                }

                const unsigned compIdx = Indices::canonicalToActiveComponentIndex(FluidSystem::solventComponentIndex(phaseIdx));
                const double ipr_rate = this->ipr_a_[compIdx] - this->ipr_b_[compIdx] * bhp_limit;

                const double rho = FluidSystem::referenceDensity( phaseIdx, Base::pvtRegionIdx() );
                total_ipr_mass_rate += ipr_rate * rho;
            }
            if ( (this->isProducer() && total_ipr_mass_rate < 0.) || (this->isInjector() && total_ipr_mass_rate > 0.) ) {
                this->operability_status_.operable_under_only_bhp_limit = false;
            }

            // checking whether running under BHP limit will violate THP limit
            if (this->operability_status_.operable_under_only_bhp_limit && this->wellHasTHPConstraints(summaryState)) {
                // option 1: calculate well rates based on the BHP limit.
                // option 2: stick with the above IPR curve
                // we use IPR here
                std::vector<double> well_rates_bhp_limit;
                computeWellRatesWithBhp(ebos_simulator, bhp_limit, well_rates_bhp_limit, deferred_logger);

                this->adaptRatesForVFP(well_rates_bhp_limit);
                const double thp = WellBhpThpCalculator(*this).calculateThpFromBhp(well_rates_bhp_limit,
                                                                                   bhp_limit,
                                                                                   this->connections_.rho(),
                                                                                   this->getALQ(well_state),
                                                                                   deferred_logger);
                const double thp_limit = this->getTHPConstraint(summaryState);
                if ( (this->isProducer() && thp < thp_limit) || (this->isInjector() && thp > thp_limit) ) {
                    this->operability_status_.obey_thp_limit_under_bhp_limit = false;
                }
            }
        } else {
            // defaulted BHP and there is a THP constraint
            // default BHP limit is about 1 atm.
            // when applied the hydrostatic pressure correction dp,
            // most likely we get a negative value (bhp + dp)to search in the VFP table,
            // which is not desirable.
            // we assume we can operate under defaulted BHP limit and will violate the THP limit
            // when operating under defaulted BHP limit.
            this->operability_status_.operable_under_only_bhp_limit = true;
            this->operability_status_.obey_thp_limit_under_bhp_limit = false;
        }
    }





    template<typename TypeTag>
    void
    StandardWell<TypeTag>::
    checkOperabilityUnderTHPLimit(const Simulator& ebos_simulator, const WellState& well_state, DeferredLogger& deferred_logger)
    {
        const auto& summaryState = ebos_simulator.vanguard().summaryState();
        const auto obtain_bhp = this->isProducer() ? computeBhpAtThpLimitProd(well_state, ebos_simulator, summaryState, deferred_logger)
        : computeBhpAtThpLimitInj(ebos_simulator, summaryState, deferred_logger);

        if (obtain_bhp) {
            this->operability_status_.can_obtain_bhp_with_thp_limit = true;

            const double  bhp_limit = WellBhpThpCalculator(*this).mostStrictBhpFromBhpLimits(summaryState);
            this->operability_status_.obey_bhp_limit_with_thp_limit = (*obtain_bhp >= bhp_limit);

            const double thp_limit = this->getTHPConstraint(summaryState);
            if (this->isProducer() && *obtain_bhp < thp_limit) {
                const std::string msg = " obtained bhp " + std::to_string(unit::convert::to(*obtain_bhp, unit::barsa))
                                        + " bars is SMALLER than thp limit "
                                        + std::to_string(unit::convert::to(thp_limit, unit::barsa))
                                        + " bars as a producer for well " + name();
                deferred_logger.debug(msg);
            }
            else if (this->isInjector() && *obtain_bhp > thp_limit) {
                const std::string msg = " obtained bhp " + std::to_string(unit::convert::to(*obtain_bhp, unit::barsa))
                                        + " bars is LARGER than thp limit "
                                        + std::to_string(unit::convert::to(thp_limit, unit::barsa))
                                        + " bars as a injector for well " + name();
                deferred_logger.debug(msg);
            }
        } else {
            this->operability_status_.can_obtain_bhp_with_thp_limit = false;
            this->operability_status_.obey_bhp_limit_with_thp_limit = false;
            if (!this->wellIsStopped()) {
                const double thp_limit = this->getTHPConstraint(summaryState);
                deferred_logger.debug(" could not find bhp value at thp limit "
                                      + std::to_string(unit::convert::to(thp_limit, unit::barsa))
                                      + " bar for well " + name() + ", the well might need to be closed ");
            }
        }
    }





    template<typename TypeTag>
    bool
    StandardWell<TypeTag>::
    allDrawDownWrongDirection(const Simulator& ebos_simulator) const
    {
        bool all_drawdown_wrong_direction = true;

        for (int perf = 0; perf < this->number_of_perforations_; ++perf) {
            const int cell_idx = this->well_cells_[perf];
            const auto& intQuants = *(ebos_simulator.model().cachedIntensiveQuantities(cell_idx, /*timeIdx=*/0));
            const auto& fs = intQuants.fluidState();

            const double pressure = this->getPerfCellPressure(fs).value();
            const double bhp = this->primary_variables_.eval(Bhp).value();

            // Pressure drawdown (also used to determine direction of flow)
            const double well_pressure = bhp + this->connections_.pressure_diff(perf);
            const double drawdown = pressure - well_pressure;

            // for now, if there is one perforation can produce/inject in the correct
            // direction, we consider this well can still produce/inject.
            // TODO: it can be more complicated than this to cause wrong-signed rates
            if ( (drawdown < 0. && this->isInjector()) ||
                 (drawdown > 0. && this->isProducer()) )  {
                all_drawdown_wrong_direction = false;
                break;
            }
        }

        const auto& comm = this->parallel_well_info_.communication();
        if (comm.size() > 1)
        {
            all_drawdown_wrong_direction =
                (comm.min(all_drawdown_wrong_direction ? 1 : 0) == 1);
        }

        return all_drawdown_wrong_direction;
    }




    template<typename TypeTag>
    bool
    StandardWell<TypeTag>::
    canProduceInjectWithCurrentBhp(const Simulator& ebos_simulator,
                                   const WellState& well_state,
                                   DeferredLogger& deferred_logger)
    {
        const double bhp = well_state.well(this->index_of_well_).bhp;
        std::vector<double> well_rates;
        computeWellRatesWithBhp(ebos_simulator, bhp, well_rates, deferred_logger);

        const double sign = (this->isProducer()) ? -1. : 1.;
        const double threshold = sign * std::numeric_limits<double>::min();

        bool can_produce_inject = false;
        for (const auto value : well_rates) {
            if (this->isProducer() && value < threshold) {
                can_produce_inject = true;
                break;
            } else if (this->isInjector() && value > threshold) {
                can_produce_inject = true;
                break;
            }
        }

        if (!can_produce_inject) {
            deferred_logger.debug(" well " + name() + " CANNOT produce or inejct ");
        }

        return can_produce_inject;
    }





    template<typename TypeTag>
    bool
    StandardWell<TypeTag>::
    openCrossFlowAvoidSingularity(const Simulator& ebos_simulator) const
    {
        return !this->getAllowCrossFlow() && allDrawDownWrongDirection(ebos_simulator);
    }




    template<typename TypeTag>
    void
    StandardWell<TypeTag>::
    computePropertiesForWellConnectionPressures(const Simulator& ebosSimulator,
                                                const WellState& well_state,
                                                std::vector<double>& b_perf,
                                                std::vector<double>& rsmax_perf,
                                                std::vector<double>& rvmax_perf,
                                                std::vector<double>& rvwmax_perf,
                                                std::vector<double>& surf_dens_perf) const
    {
        std::function<Scalar(int,int)> getTemperature =
        [&ebosSimulator](int cell_idx, int phase_idx)
        {
            return ebosSimulator.model().cachedIntensiveQuantities(cell_idx, 0)->fluidState().temperature(phase_idx).value();
        };
        std::function<Scalar(int)> getSaltConcentration =
        [&ebosSimulator](int cell_idx)
        {
            return ebosSimulator.model().cachedIntensiveQuantities(cell_idx, 0)->fluidState().saltConcentration().value();
        };
        std::function<int(int)> getPvtRegionIdx =
        [&ebosSimulator](int cell_idx)
        {
            return ebosSimulator.model().cachedIntensiveQuantities(cell_idx, 0)->fluidState().pvtRegionIndex();
        };
        std::function<Scalar(int)> getInvFac =
        [&ebosSimulator](int cell_idx)
        {
            return ebosSimulator.model().cachedIntensiveQuantities(cell_idx, 0)->solventInverseFormationVolumeFactor().value();
        };
        std::function<Scalar(int)> getSolventDensity =
        [&ebosSimulator](int cell_idx)
        {
            return ebosSimulator.model().cachedIntensiveQuantities(cell_idx, 0)->solventRefDensity();
        };

        this->connections_.computePropertiesForPressures(well_state,
                                                         getTemperature,
                                                         getSaltConcentration,
                                                         getPvtRegionIdx,
                                                         getInvFac,
                                                         getSolventDensity,
                                                         b_perf,
                                                         rsmax_perf,
                                                         rvmax_perf,
                                                         rvwmax_perf,
                                                         surf_dens_perf);
    }





    template<typename TypeTag>
    ConvergenceReport
    StandardWell<TypeTag>::
    getWellConvergence(const WellState& well_state,
                       const std::vector<double>& B_avg,
                       DeferredLogger& deferred_logger,
                       const bool relax_tolerance) const
    {
        // the following implementation assume that the polymer is always after the w-o-g phases
        // For the polymer, energy and foam cases, there is one more mass balance equations of reservoir than wells
        assert((int(B_avg.size()) == this->num_components_) || has_polymer || has_energy || has_foam || has_brine || has_zFraction || has_micp);

        std::vector<double> res;
        ConvergenceReport report = this->StdWellEval::getWellConvergence(well_state,
                                                                         B_avg,
                                                                         this->param_.max_residual_allowed_,
                                                                         this->param_.tolerance_wells_,
                                                                         this->param_.relaxed_tolerance_flow_well_,
                                                                         relax_tolerance,
                                                                         res,
                                                                         deferred_logger);
        checkConvergenceExtraEqs(res, report);

        return report;
    }





    template<typename TypeTag>
    void
    StandardWell<TypeTag>::
    updateProductivityIndex(const Simulator& ebosSimulator,
                            const WellProdIndexCalculator& wellPICalc,
                            WellState& well_state,
                            DeferredLogger& deferred_logger) const
    {
        auto fluidState = [&ebosSimulator, this](const int perf)
        {
            const auto cell_idx = this->well_cells_[perf];
            return ebosSimulator.model()
               .cachedIntensiveQuantities(cell_idx, /*timeIdx=*/ 0)->fluidState();
        };

        const int np = this->number_of_phases_;
        auto setToZero = [np](double* x) -> void
        {
            std::fill_n(x, np, 0.0);
        };

        auto addVector = [np](const double* src, double* dest) -> void
        {
            std::transform(src, src + np, dest, dest, std::plus<>{});
        };

        auto& ws = well_state.well(this->index_of_well_);
        auto& perf_data = ws.perf_data;
        auto* wellPI = ws.productivity_index.data();
        auto* connPI = perf_data.prod_index.data();

        setToZero(wellPI);

        const auto preferred_phase = this->well_ecl_.getPreferredPhase();
        auto subsetPerfID = 0;

        for (const auto& perf : *this->perf_data_) {
            auto allPerfID = perf.ecl_index;

            auto connPICalc = [&wellPICalc, allPerfID](const double mobility) -> double
            {
                return wellPICalc.connectionProdIndStandard(allPerfID, mobility);
            };

            std::vector<EvalWell> mob(this->num_components_, {this->primary_variables_.numWellEq() + Indices::numEq, 0.0});
            getMobilityEval(ebosSimulator, static_cast<int>(subsetPerfID), mob, deferred_logger);

            const auto& fs = fluidState(subsetPerfID);
            setToZero(connPI);

            if (this->isInjector()) {
                this->computeConnLevelInjInd(fs, preferred_phase, connPICalc,
                                             mob, connPI, deferred_logger);
            }
            else {  // Production or zero flow rate
                this->computeConnLevelProdInd(fs, connPICalc, mob, connPI);
            }

            addVector(connPI, wellPI);

            ++subsetPerfID;
            connPI += np;
        }

        // Sum with communication in case of distributed well.
        const auto& comm = this->parallel_well_info_.communication();
        if (comm.size() > 1) {
            comm.sum(wellPI, np);
        }

        assert ((static_cast<int>(subsetPerfID) == this->number_of_perforations_) &&
                "Internal logic error in processing connections for PI/II");
    }



    template<typename TypeTag>
    void
    StandardWell<TypeTag>::
    computeWellConnectionDensitesPressures(const Simulator& ebosSimulator,
                                           const WellState& well_state,
                                           const std::vector<double>& b_perf,
                                           const std::vector<double>& rsmax_perf,
                                           const std::vector<double>& rvmax_perf,
                                           const std::vector<double>& rvwmax_perf,
                                           const std::vector<double>& surf_dens_perf,
                                           DeferredLogger& deferred_logger)
    {
        std::function<Scalar(int,int)> invB =
        [&ebosSimulator](int cell_idx, int phase_idx)
        {
            return ebosSimulator.model().cachedIntensiveQuantities(cell_idx, 0)->fluidState().invB(phase_idx).value();
        };
        std::function<Scalar(int,int)> mobility =
        [&ebosSimulator](int cell_idx, int phase_idx)
        {
            return ebosSimulator.model().cachedIntensiveQuantities(cell_idx, 0)->mobility(phase_idx).value();
        };
        std::function<Scalar(int)> invFac =
        [&ebosSimulator](int cell_idx)
        {
            return ebosSimulator.model().cachedIntensiveQuantities(cell_idx, 0)->solventInverseFormationVolumeFactor().value();
        };
        std::function<Scalar(int)> solventMobility =
        [&ebosSimulator](int cell_idx)
        {
            return ebosSimulator.model().cachedIntensiveQuantities(cell_idx, 0)->solventMobility().value();
        };

        this->connections_.computeProperties(well_state,
                                             invB,
                                             mobility,
                                             invFac,
                                             solventMobility,
                                             b_perf,
                                             rsmax_perf,
                                             rvmax_perf,
                                             rvwmax_perf,
                                             surf_dens_perf,
                                             deferred_logger);
    }





    template<typename TypeTag>
    void
    StandardWell<TypeTag>::
    computeWellConnectionPressures(const Simulator& ebosSimulator,
                                   const WellState& well_state,
                                   DeferredLogger& deferred_logger)
    {
         // 1. Compute properties required by computePressureDelta().
         //    Note that some of the complexity of this part is due to the function
         //    taking std::vector<double> arguments, and not Eigen objects.
         std::vector<double> b_perf;
         std::vector<double> rsmax_perf;
         std::vector<double> rvmax_perf;
         std::vector<double> rvwmax_perf;
         std::vector<double> surf_dens_perf;
         computePropertiesForWellConnectionPressures(ebosSimulator, well_state, b_perf, rsmax_perf, rvmax_perf, rvwmax_perf, surf_dens_perf);
         computeWellConnectionDensitesPressures(ebosSimulator, well_state, b_perf, rsmax_perf, rvmax_perf, rvwmax_perf, surf_dens_perf, deferred_logger);
    }





    template<typename TypeTag>
    void
    StandardWell<TypeTag>::
    solveEqAndUpdateWellState(WellState& well_state, DeferredLogger& deferred_logger)
    {
        if (!this->isOperableAndSolvable() && !this->wellIsStopped()) return;

        // We assemble the well equations, then we check the convergence,
        // which is why we do not put the assembleWellEq here.
        BVectorWell dx_well(1);
        dx_well[0].resize(this->primary_variables_.numWellEq());
        this->linSys_.solve( dx_well);

        updateWellState(dx_well, well_state, deferred_logger);
    }





    template<typename TypeTag>
    void
    StandardWell<TypeTag>::
    calculateExplicitQuantities(const Simulator& ebosSimulator,
                                const WellState& well_state,
                                DeferredLogger& deferred_logger)
    {
        updatePrimaryVariables(well_state, deferred_logger);
        initPrimaryVariablesEvaluation();
        computeWellConnectionPressures(ebosSimulator, well_state, deferred_logger);
        this->computeAccumWell();
    }



    template<typename TypeTag>
    void
    StandardWell<TypeTag>::
    apply(const BVector& x, BVector& Ax) const
    {
        if (!this->isOperableAndSolvable() && !this->wellIsStopped()) return;

        if (this->param_.matrix_add_well_contributions_)
        {
            // Contributions are already in the matrix itself
            return;
        }

        this->linSys_.apply(x, Ax);
    }




    template<typename TypeTag>
    void
    StandardWell<TypeTag>::
    apply(BVector& r) const
    {
        if (!this->isOperableAndSolvable() && !this->wellIsStopped()) return;

        this->linSys_.apply(r);
    }




    template<typename TypeTag>
    void
    StandardWell<TypeTag>::
    recoverWellSolutionAndUpdateWellState(const BVector& x,
                                          WellState& well_state,
                                          DeferredLogger& deferred_logger)
    {
        if (!this->isOperableAndSolvable() && !this->wellIsStopped()) return;

        BVectorWell xw(1);
        xw[0].resize(this->primary_variables_.numWellEq());

        this->linSys_.recoverSolutionWell(x, xw);
        updateWellState(xw, well_state, deferred_logger);
    }




    template<typename TypeTag>
    void
    StandardWell<TypeTag>::
    computeWellRatesWithBhp(const Simulator& ebosSimulator,
                            const double& bhp,
                            std::vector<double>& well_flux,
                            DeferredLogger& deferred_logger) const
    {

        const int np = this->number_of_phases_;
        well_flux.resize(np, 0.0);

        const bool allow_cf = this->getAllowCrossFlow();

        for (int perf = 0; perf < this->number_of_perforations_; ++perf) {
            const int cell_idx = this->well_cells_[perf];
            const auto& intQuants = *(ebosSimulator.model().cachedIntensiveQuantities(cell_idx, /*timeIdx=*/ 0));
            // flux for each perforation
            std::vector<Scalar> mob(this->num_components_, 0.);
            getMobilityScalar(ebosSimulator, perf, mob, deferred_logger);
            double trans_mult = ebosSimulator.problem().template rockCompTransMultiplier<double>(intQuants, cell_idx);
            const double Tw = this->well_index_[perf] * trans_mult;

            std::vector<Scalar> cq_s(this->num_components_, 0.);
            computePerfRateScalar(intQuants, mob, bhp, Tw, perf, allow_cf,
                            cq_s, deferred_logger);

            for(int p = 0; p < np; ++p) {
                well_flux[this->ebosCompIdxToFlowCompIdx(p)] += cq_s[p];
            }
        }
        this->parallel_well_info_.communication().sum(well_flux.data(), well_flux.size());
    }



    template<typename TypeTag>
    void
    StandardWell<TypeTag>::
    computeWellRatesWithBhpIterations(const Simulator& ebosSimulator,
                                      const double& bhp,
                                      std::vector<double>& well_flux,
                                      DeferredLogger& deferred_logger) const
    {

        // iterate to get a more accurate well density
        // create a copy of the well_state to use. If the operability checking is sucessful, we use this one
        // to replace the original one
        WellState well_state_copy = ebosSimulator.problem().wellModel().wellState();
        const auto& group_state  = ebosSimulator.problem().wellModel().groupState();
        auto& ws = well_state_copy.well(this->index_of_well_);

        //  Set current control to bhp, and bhp value in state, modify bhp limit in control object.
        if (this->well_ecl_.isInjector()) {
            ws.injection_cmode = Well::InjectorCMode::BHP;
        } else {
            ws.production_cmode = Well::ProducerCMode::BHP;
        }
        ws.bhp = bhp;

        // initialized the well rates with the potentials i.e. the well rates based on bhp
        const int np = this->number_of_phases_;
        const double sign = this->well_ecl_.isInjector() ? 1.0 : -1.0;
        for (int phase = 0; phase < np; ++phase){
            well_state_copy.wellRates(this->index_of_well_)[phase]
                    = sign * ws.well_potentials[phase];
        }
        // creating a copy of the well itself, to avoid messing up the explicit informations
        // during this copy, the only information not copied properly is the well controls
        StandardWell<TypeTag> well(*this);
        well.calculateExplicitQuantities(ebosSimulator, well_state_copy, deferred_logger);

        const double dt = ebosSimulator.timeStepSize();
        bool converged = well.iterateWellEquations(ebosSimulator, dt, well_state_copy, group_state, deferred_logger);
        if (!converged) {
            const std::string msg = " well " + name() + " did not get converged during well potential calculations "
                                                        " potentials are computed based on unconverged solution";
            deferred_logger.debug(msg);
        }
        well.updatePrimaryVariables(well_state_copy, deferred_logger);
        well.computeWellConnectionPressures(ebosSimulator, well_state_copy, deferred_logger);
        well.initPrimaryVariablesEvaluation();
        well.computeWellRatesWithBhp(ebosSimulator, bhp, well_flux, deferred_logger);
    }




    template<typename TypeTag>
    std::vector<double>
    StandardWell<TypeTag>::
    computeWellPotentialWithTHP(const Simulator& ebos_simulator,
                               DeferredLogger& deferred_logger,
                               const WellState &well_state) const
    {
        std::vector<double> potentials(this->number_of_phases_, 0.0);
        const auto& summary_state = ebos_simulator.vanguard().summaryState();

        const auto& well = this->well_ecl_;
        if (well.isInjector()){
            const auto& controls = this->well_ecl_.injectionControls(summary_state);
            auto bhp_at_thp_limit = computeBhpAtThpLimitInj(ebos_simulator, summary_state, deferred_logger);
            if (bhp_at_thp_limit) {
                const double bhp = std::min(*bhp_at_thp_limit, controls.bhp_limit);
                computeWellRatesWithBhp(ebos_simulator, bhp, potentials, deferred_logger);
            } else {
                deferred_logger.warning("FAILURE_GETTING_CONVERGED_POTENTIAL",
                                        "Failed in getting converged thp based potential calculation for well "
                                        + name() + ". Instead the bhp based value is used");
                const double bhp = controls.bhp_limit;
                computeWellRatesWithBhp(ebos_simulator, bhp, potentials, deferred_logger);
            }
        } else {
            computeWellRatesWithThpAlqProd(
                ebos_simulator, summary_state,
                deferred_logger, potentials, this->getALQ(well_state)
            );
        }

        return potentials;
    }

    template<typename TypeTag>
    double
    StandardWell<TypeTag>::
    computeWellRatesAndBhpWithThpAlqProd(const Simulator &ebos_simulator,
                               const SummaryState &summary_state,
                               DeferredLogger &deferred_logger,
                               std::vector<double> &potentials,
                               double alq) const
    {
        double bhp;
        auto bhp_at_thp_limit = computeBhpAtThpLimitProdWithAlq(
                              ebos_simulator, summary_state, alq, deferred_logger);
        if (bhp_at_thp_limit) {
            const auto& controls = this->well_ecl_.productionControls(summary_state);
            bhp = std::max(*bhp_at_thp_limit, controls.bhp_limit);
            computeWellRatesWithBhp(ebos_simulator, bhp, potentials, deferred_logger);
        }
        else {
            deferred_logger.warning("FAILURE_GETTING_CONVERGED_POTENTIAL",
                "Failed in getting converged thp based potential calculation for well "
                + name() + ". Instead the bhp based value is used");
            const auto& controls = this->well_ecl_.productionControls(summary_state);
            bhp = controls.bhp_limit;
            computeWellRatesWithBhp(ebos_simulator, bhp, potentials, deferred_logger);
        }
        return bhp;
    }

    template<typename TypeTag>
    void
    StandardWell<TypeTag>::
    computeWellRatesWithThpAlqProd(const Simulator &ebos_simulator,
                               const SummaryState &summary_state,
                               DeferredLogger &deferred_logger,
                               std::vector<double> &potentials,
                               double alq) const
    {
        /*double bhp =*/
        computeWellRatesAndBhpWithThpAlqProd(ebos_simulator,
                                             summary_state,
                                             deferred_logger,
                                             potentials,
                                             alq);
    }

    template<typename TypeTag>
    void
    StandardWell<TypeTag>::
    computeWellPotentials(const Simulator& ebosSimulator,
                          const WellState& well_state,
                          std::vector<double>& well_potentials,
                          DeferredLogger& deferred_logger) // const
    {
        const int np = this->number_of_phases_;
        well_potentials.resize(np, 0.0);

        if (this->wellIsStopped()) {
            return;
        }

        this->operability_status_.has_negative_potentials = false;
        // If the well is pressure controlled the potential equals the rate.
        bool thp_controlled_well = false;
        bool bhp_controlled_well = false;
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
        if (!this->changed_to_open_this_step_ && (thp_controlled_well || bhp_controlled_well)) {

            double total_rate = 0.0;
            const double sign = this->isInjector() ? 1.0:-1.0;
            for (int phase = 0; phase < np; ++phase){
                total_rate += sign * ws.surface_rates[phase];
            }
            // for pressure controlled wells the well rates are the potentials
            // if the rates are trivial we are most probably looking at the newly
            // opened well and we therefore make the affort of computing the potentials anyway.
            if (total_rate > 0) {
                for (int phase = 0; phase < np; ++phase){
                    well_potentials[phase] = sign * ws.surface_rates[phase];
                }
                return;
            }
        }

        // does the well have a THP related constraint?
        const auto& summaryState = ebosSimulator.vanguard().summaryState();
        if (!Base::wellHasTHPConstraints(summaryState) || bhp_controlled_well) {
            // get the bhp value based on the bhp constraints
            double bhp = WellBhpThpCalculator(*this).mostStrictBhpFromBhpLimits(summaryState);

            // In some very special cases the bhp pressure target are
            // temporary violated. This may lead to too small or negative potentials
            // that could lead to premature shutting of wells.
            // As a remedy the bhp that gives the largest potential is used.
            // For converged cases, ws.bhp <=bhp for injectors and ws.bhp >= bhp,
            // and the potentials will be computed using the limit as expected.
            if (this->isInjector())
                bhp = std::max(ws.bhp, bhp);
            else
                bhp = std::min(ws.bhp, bhp);

            assert(std::abs(bhp) != std::numeric_limits<double>::max());
            computeWellRatesWithBhpIterations(ebosSimulator, bhp, well_potentials, deferred_logger);
        } else {
            // the well has a THP related constraint
            well_potentials = computeWellPotentialWithTHP(ebosSimulator, deferred_logger, well_state);
        }

        const double sign = this->isInjector() ? 1.0:-1.0;
        double total_potential = 0.0;
        for (int phase = 0; phase < np; ++phase){
            well_potentials[phase] *= sign;
            total_potential += well_potentials[phase];
        }
        if (total_potential < 0.0 && this->param_.check_well_operability_) {
            // wells with negative potentials are not operable
            this->operability_status_.has_negative_potentials = true;
            const std::string msg = std::string("well ") + this->name() + std::string(": has negative potentials and is not operable");
            deferred_logger.warning("NEGATIVE_POTENTIALS_INOPERABLE", msg);
        }
    }





    template<typename TypeTag>
    void
    StandardWell<TypeTag>::
    updatePrimaryVariables(const WellState& well_state, DeferredLogger& deferred_logger)
    {
        if (!this->isOperableAndSolvable() && !this->wellIsStopped()) return;

        this->primary_variables_.update(well_state, deferred_logger);

        // other primary variables related to polymer injection
        if constexpr (Base::has_polymermw) {
            this->primary_variables_.updatePolyMW(well_state);
        }

        this->primary_variables_.checkFinite(deferred_logger);
    }




    template<typename TypeTag>
    double
    StandardWell<TypeTag>::
    getRefDensity() const
    {
        return this->connections_.rho();
    }




    template<typename TypeTag>
    void
    StandardWell<TypeTag>::
    updateWaterMobilityWithPolymer(const Simulator& ebos_simulator,
                                   const int perf,
                                   std::vector<EvalWell>& mob,
                                   DeferredLogger& deferred_logger) const
    {
        const int cell_idx = this->well_cells_[perf];
        const auto& int_quant = *(ebos_simulator.model().cachedIntensiveQuantities(cell_idx, /*timeIdx=*/ 0));
        const EvalWell polymer_concentration = this->extendEval(int_quant.polymerConcentration());

        // TODO: not sure should based on the well type or injecting/producing peforations
        // it can be different for crossflow
        if (this->isInjector()) {
            // assume fully mixing within injecting wellbore
            const auto& visc_mult_table = PolymerModule::plyviscViscosityMultiplierTable(int_quant.pvtRegionIndex());
            const unsigned waterCompIdx = Indices::canonicalToActiveComponentIndex(FluidSystem::waterCompIdx);
            mob[waterCompIdx] /= (this->extendEval(int_quant.waterViscosityCorrection()) * visc_mult_table.eval(polymer_concentration, /*extrapolate=*/true) );
        }

        if (PolymerModule::hasPlyshlog()) {
            // we do not calculate the shear effects for injection wells when they do not
            // inject polymer.
            if (this->isInjector() && this->wpolymer() == 0.) {
                return;
            }
            // compute the well water velocity with out shear effects.
            // TODO: do we need to turn on crossflow here?
            const bool allow_cf = this->getAllowCrossFlow() || openCrossFlowAvoidSingularity(ebos_simulator);
            const EvalWell& bhp = this->primary_variables_.eval(Bhp);

            std::vector<EvalWell> cq_s(this->num_components_, {this->primary_variables_.numWellEq() + Indices::numEq, 0.});
            double perf_dis_gas_rate = 0.;
            double perf_vap_oil_rate = 0.;
            double perf_vap_wat_rate = 0.;
            double trans_mult = ebos_simulator.problem().template rockCompTransMultiplier<double>(int_quant, cell_idx);
            const double Tw = this->well_index_[perf] * trans_mult;
            computePerfRateEval(int_quant, mob, bhp, Tw, perf, allow_cf,
                                cq_s, perf_dis_gas_rate, perf_vap_oil_rate, perf_vap_wat_rate, deferred_logger);
            // TODO: make area a member
            const double area = 2 * M_PI * this->perf_rep_radius_[perf] * this->perf_length_[perf];
            const auto& material_law_manager = ebos_simulator.problem().materialLawManager();
            const auto& scaled_drainage_info =
                        material_law_manager->oilWaterScaledEpsInfoDrainage(cell_idx);
            const double swcr = scaled_drainage_info.Swcr;
            const EvalWell poro = this->extendEval(int_quant.porosity());
            const EvalWell sw = this->extendEval(int_quant.fluidState().saturation(FluidSystem::waterPhaseIdx));
            // guard against zero porosity and no water
            const EvalWell denom = max( (area * poro * (sw - swcr)), 1e-12);
            const unsigned waterCompIdx = Indices::canonicalToActiveComponentIndex(FluidSystem::waterCompIdx);
            EvalWell water_velocity = cq_s[waterCompIdx] / denom * this->extendEval(int_quant.fluidState().invB(FluidSystem::waterPhaseIdx));

            if (PolymerModule::hasShrate()) {
                // the equation for the water velocity conversion for the wells and reservoir are from different version
                // of implementation. It can be changed to be more consistent when possible.
                water_velocity *= PolymerModule::shrate( int_quant.pvtRegionIndex() ) / this->bore_diameters_[perf];
            }
            const EvalWell shear_factor = PolymerModule::computeShearFactor(polymer_concentration,
                                                                int_quant.pvtRegionIndex(),
                                                                water_velocity);
             // modify the mobility with the shear factor.
            mob[waterCompIdx] /= shear_factor;
        }
    }

    template<typename TypeTag>
    void
    StandardWell<TypeTag>::addWellContributions(SparseMatrixAdapter& jacobian) const
    {
        this->linSys_.extract(jacobian);
    }

    
    template <typename TypeTag>
    void
    StandardWell<TypeTag>::addWellPressureEquations(PressureMatrix& jacobian,
                                                    const BVector& weights,
                                                    const int pressureVarIndex,
                                                    const bool use_well_weights,
                                                    const WellState& well_state) const
    {
        this->linSys_.extractCPRPressureMatrix(jacobian,
                                               weights,
                                               pressureVarIndex,
                                               use_well_weights,
                                               *this,
                                               Bhp,
                                               well_state);
    }

    

    template<typename TypeTag>
    typename StandardWell<TypeTag>::EvalWell
    StandardWell<TypeTag>::
    pskinwater(const double throughput,
               const EvalWell& water_velocity,
              DeferredLogger& deferred_logger) const
    {
        if constexpr (Base::has_polymermw) {
            const int water_table_id = this->well_ecl_.getPolymerProperties().m_skprwattable;
            if (water_table_id <= 0) {
                OPM_DEFLOG_THROW(std::runtime_error, "Unused SKPRWAT table id used for well " << name(), deferred_logger);
            }
            const auto& water_table_func = PolymerModule::getSkprwatTable(water_table_id);
            const EvalWell throughput_eval(this->primary_variables_.numWellEq() + Indices::numEq, throughput);
            // the skin pressure when injecting water, which also means the polymer concentration is zero
            EvalWell pskin_water(this->primary_variables_.numWellEq() + Indices::numEq, 0.0);
            pskin_water = water_table_func.eval(throughput_eval, water_velocity);
            return pskin_water;
        } else {
            OPM_DEFLOG_THROW(std::runtime_error, "Polymermw is not activated, "
                                          "while injecting skin pressure is requested for well " << name(), deferred_logger);
        }
    }





    template<typename TypeTag>
    typename StandardWell<TypeTag>::EvalWell
    StandardWell<TypeTag>::
    pskin(const double throughput,
              const EvalWell& water_velocity,
              const EvalWell& poly_inj_conc,
              DeferredLogger& deferred_logger) const
    {
        if constexpr (Base::has_polymermw) {
            const double sign = water_velocity >= 0. ? 1.0 : -1.0;
            const EvalWell water_velocity_abs = abs(water_velocity);
            if (poly_inj_conc == 0.) {
                return sign * pskinwater(throughput, water_velocity_abs, deferred_logger);
            }
            const int polymer_table_id = this->well_ecl_.getPolymerProperties().m_skprpolytable;
            if (polymer_table_id <= 0) {
                OPM_DEFLOG_THROW(std::runtime_error, "Unavailable SKPRPOLY table id used for well " << name(), deferred_logger);
            }
            const auto& skprpolytable = PolymerModule::getSkprpolyTable(polymer_table_id);
            const double reference_concentration = skprpolytable.refConcentration;
            const EvalWell throughput_eval(this->primary_variables_.numWellEq() + Indices::numEq, throughput);
            // the skin pressure when injecting water, which also means the polymer concentration is zero
            EvalWell pskin_poly(this->primary_variables_.numWellEq() + Indices::numEq, 0.0);
            pskin_poly = skprpolytable.table_func.eval(throughput_eval, water_velocity_abs);
            if (poly_inj_conc == reference_concentration) {
                return sign * pskin_poly;
            }
            // poly_inj_conc != reference concentration of the table, then some interpolation will be required
            const EvalWell pskin_water = pskinwater(throughput, water_velocity_abs, deferred_logger);
            const EvalWell pskin = pskin_water + (pskin_poly - pskin_water) / reference_concentration * poly_inj_conc;
            return sign * pskin;
        } else {
            OPM_DEFLOG_THROW(std::runtime_error, "Polymermw is not activated, "
                                          "while injecting skin pressure is requested for well " << name(), deferred_logger);
        }
    }





    template<typename TypeTag>
    typename StandardWell<TypeTag>::EvalWell
    StandardWell<TypeTag>::
    wpolymermw(const double throughput,
               const EvalWell& water_velocity,
               DeferredLogger& deferred_logger) const
    {
        if constexpr (Base::has_polymermw) {
            const int table_id = this->well_ecl_.getPolymerProperties().m_plymwinjtable;
            const auto& table_func = PolymerModule::getPlymwinjTable(table_id);
            const EvalWell throughput_eval(this->primary_variables_.numWellEq() + Indices::numEq, throughput);
            EvalWell molecular_weight(this->primary_variables_.numWellEq() + Indices::numEq, 0.);
            if (this->wpolymer() == 0.) { // not injecting polymer
                return molecular_weight;
            }
            molecular_weight = table_func.eval(throughput_eval, abs(water_velocity));
            return molecular_weight;
        } else {
            OPM_DEFLOG_THROW(std::runtime_error, "Polymermw is not activated, "
                                          "while injecting polymer molecular weight is requested for well " << name(), deferred_logger);
        }
    }





    template<typename TypeTag>
    void
    StandardWell<TypeTag>::
    updateWaterThroughput(const double dt, WellState &well_state) const
    {
        if constexpr (Base::has_polymermw) {
            if (this->isInjector()) {
                auto& ws = well_state.well(this->index_of_well_);
                auto& perf_water_throughput = ws.perf_data.water_throughput;
                for (int perf = 0; perf < this->number_of_perforations_; ++perf) {
                    const double perf_water_vel = this->primary_variables_.value(Bhp + 1 + perf);
                    // we do not consider the formation damage due to water flowing from reservoir into wellbore
                    if (perf_water_vel > 0.) {
                        perf_water_throughput[perf] += perf_water_vel * dt;
                    }
                }
            }
        }
    }





    template<typename TypeTag>
    void
    StandardWell<TypeTag>::
    handleInjectivityRate(const Simulator& ebosSimulator,
                          const int perf,
                          std::vector<EvalWell>& cq_s) const
    {
        const int cell_idx = this->well_cells_[perf];
        const auto& int_quants = *(ebosSimulator.model().cachedIntensiveQuantities(cell_idx, /*timeIdx=*/ 0));
        const auto& fs = int_quants.fluidState();
        const EvalWell b_w = this->extendEval(fs.invB(FluidSystem::waterPhaseIdx));
        const double area = M_PI * this->bore_diameters_[perf] * this->perf_length_[perf];
        const int wat_vel_index = Bhp + 1 + perf;
        const unsigned water_comp_idx = Indices::canonicalToActiveComponentIndex(FluidSystem::waterCompIdx);

        // water rate is update to use the form from water velocity, since water velocity is
        // a primary variable now
        cq_s[water_comp_idx] = area * this->primary_variables_.eval(wat_vel_index) * b_w;
    }




    template<typename TypeTag>
    void
    StandardWell<TypeTag>::
    handleInjectivityEquations(const Simulator& ebosSimulator,
                               const WellState& well_state,
                               const int perf,
                               const EvalWell& water_flux_s,
                               DeferredLogger& deferred_logger)
    {
        const int cell_idx = this->well_cells_[perf];
        const auto& int_quants = *(ebosSimulator.model().cachedIntensiveQuantities(cell_idx, /*timeIdx=*/ 0));
        const auto& fs = int_quants.fluidState();
        const EvalWell b_w = this->extendEval(fs.invB(FluidSystem::waterPhaseIdx));
        const EvalWell water_flux_r = water_flux_s / b_w;
        const double area = M_PI * this->bore_diameters_[perf] * this->perf_length_[perf];
        const EvalWell water_velocity = water_flux_r / area;
        const int wat_vel_index = Bhp + 1 + perf;

        // equation for the water velocity
        const EvalWell eq_wat_vel = this->primary_variables_.eval(wat_vel_index) - water_velocity;

        const auto& ws = well_state.well(this->index_of_well_);
        const auto& perf_data = ws.perf_data;
        const auto& perf_water_throughput = perf_data.water_throughput;
        const double throughput = perf_water_throughput[perf];
        const int pskin_index = Bhp + 1 + this->number_of_perforations_ + perf;

        EvalWell poly_conc(this->primary_variables_.numWellEq() + Indices::numEq, 0.0);
        poly_conc.setValue(this->wpolymer());

        // equation for the skin pressure
        const EvalWell eq_pskin = this->primary_variables_.eval(pskin_index)
                                  - pskin(throughput, this->primary_variables_.eval(wat_vel_index), poly_conc, deferred_logger);

        StandardWellAssemble<FluidSystem,Indices,Scalar>(*this).
                assembleInjectivityEq(eq_pskin,
                                      eq_wat_vel,
                                      pskin_index,
                                      wat_vel_index,
                                      cell_idx,
                                      this->primary_variables_.numWellEq(),
                                      this->linSys_);
    }





    template<typename TypeTag>
    void
    StandardWell<TypeTag>::
    checkConvergenceExtraEqs(const std::vector<double>& res,
                             ConvergenceReport& report) const
    {
        // if different types of extra equations are involved, this function needs to be refactored further

        // checking the convergence of the extra equations related to polymer injectivity
        if constexpr (Base::has_polymermw) {
            WellConvergence(*this).
                checkConvergencePolyMW(res, Bhp, this->param_.max_residual_allowed_, report);
        }
    }





    template<typename TypeTag>
    void
    StandardWell<TypeTag>::
    updateConnectionRatePolyMW(const EvalWell& cq_s_poly,
                               const IntensiveQuantities& int_quants,
                               const WellState& well_state,
                               const int perf,
                               std::vector<RateVector>& connectionRates,
                               DeferredLogger& deferred_logger) const
    {
        // the source term related to transport of molecular weight
        EvalWell cq_s_polymw = cq_s_poly;
        if (this->isInjector()) {
            const int wat_vel_index = Bhp + 1 + perf;
            const EvalWell water_velocity = this->primary_variables_.eval(wat_vel_index);
            if (water_velocity > 0.) { // injecting
                const auto& ws = well_state.well(this->index_of_well_);
                const auto& perf_water_throughput = ws.perf_data.water_throughput;
                const double throughput = perf_water_throughput[perf];
                const EvalWell molecular_weight = wpolymermw(throughput, water_velocity, deferred_logger);
                cq_s_polymw *= molecular_weight;
            } else {
                // we do not consider the molecular weight from the polymer
                // going-back to the wellbore through injector
                cq_s_polymw *= 0.;
            }
        } else if (this->isProducer()) {
            if (cq_s_polymw < 0.) {
                cq_s_polymw *= this->extendEval(int_quants.polymerMoleWeight() );
            } else {
                // we do not consider the molecular weight from the polymer
                // re-injecting back through producer
                cq_s_polymw *= 0.;
            }
        }
        connectionRates[perf][Indices::contiPolymerMWEqIdx] = Base::restrictEval(cq_s_polymw);
    }






    template<typename TypeTag>
    std::optional<double>
    StandardWell<TypeTag>::
    computeBhpAtThpLimitProd(const WellState& well_state,
                             const Simulator& ebos_simulator,
                             const SummaryState& summary_state,
                             DeferredLogger& deferred_logger) const
    {
        return computeBhpAtThpLimitProdWithAlq(ebos_simulator,
                                               summary_state,
                                               this->getALQ(well_state),
                                               deferred_logger);
    }

    template<typename TypeTag>
    std::optional<double>
    StandardWell<TypeTag>::
    computeBhpAtThpLimitProdWithAlq(const Simulator& ebos_simulator,
                                    const SummaryState& summary_state,
                                    const double alq_value,
                                    DeferredLogger& deferred_logger) const
    {
        // Make the frates() function.
        auto frates = [this, &ebos_simulator, &deferred_logger](const double bhp) {
            // Not solving the well equations here, which means we are
            // calculating at the current Fg/Fw values of the
            // well. This does not matter unless the well is
            // crossflowing, and then it is likely still a good
            // approximation.
            std::vector<double> rates(3);
            computeWellRatesWithBhp(ebos_simulator, bhp, rates, deferred_logger);
            this->adaptRatesForVFP(rates);
            return rates;
        };

        double max_pressure = 0.0;
        for (int perf = 0; perf < this->number_of_perforations_; ++perf) {
            const int cell_idx = this->well_cells_[perf];
            const auto& int_quants = *(ebos_simulator.model().cachedIntensiveQuantities(cell_idx, /*timeIdx=*/ 0));
            const auto& fs = int_quants.fluidState();
            double pressure_cell = this->getPerfCellPressure(fs).value();
            max_pressure = std::max(max_pressure, pressure_cell);
        }
        auto bhpAtLimit = WellBhpThpCalculator(*this).computeBhpAtThpLimitProd(frates,
                                                                               summary_state,
                                                                               max_pressure,
                                                                               this->connections_.rho(),
                                                                               alq_value,
                                                                               this->getTHPConstraint(summary_state),
                                                                               deferred_logger);
        auto v = frates(*bhpAtLimit);
        if (bhpAtLimit && std::all_of(v.cbegin(), v.cend(), [](double i){ return i <= 0; }))
            return bhpAtLimit;

        auto fratesIter = [this, &ebos_simulator, &deferred_logger](const double bhp) {
            // Solver the well iterations to see if we are
            // able to get a solution with an update
            // solution
            std::vector<double> rates(3);
            computeWellRatesWithBhpIterations(ebos_simulator, bhp, rates, deferred_logger);
            this->adaptRatesForVFP(rates);
            return rates;
        };

        bhpAtLimit = WellBhpThpCalculator(*this).computeBhpAtThpLimitProd(fratesIter,
                                                                          summary_state,
                                                                          max_pressure,
                                                                          this->connections_.rho(),
                                                                          alq_value,
                                                                          this->getTHPConstraint(summary_state),
                                                                          deferred_logger);
        v = frates(*bhpAtLimit);
        if(bhpAtLimit && std::all_of(v.cbegin(), v.cend(), [](double i){ return i <= 0; }))
            return bhpAtLimit;

        // we still don't get a valied solution.
        return std::nullopt;
    }



    template<typename TypeTag>
    std::optional<double>
    StandardWell<TypeTag>::
    computeBhpAtThpLimitInj(const Simulator& ebos_simulator,
                            const SummaryState& summary_state,
                            DeferredLogger& deferred_logger) const
    {
        // Make the frates() function.
        auto frates = [this, &ebos_simulator, &deferred_logger](const double bhp) {
            // Not solving the well equations here, which means we are
            // calculating at the current Fg/Fw values of the
            // well. This does not matter unless the well is
            // crossflowing, and then it is likely still a good
            // approximation.
            std::vector<double> rates(3);
            computeWellRatesWithBhp(ebos_simulator, bhp, rates, deferred_logger);
            return rates;
        };

        return WellBhpThpCalculator(*this).computeBhpAtThpLimitInj(frates,
                                                                   summary_state,
                                                                   this->connections_.rho(),
                                                                   1e-6,
                                                                   50,
                                                                   true,
                                                                   deferred_logger);
    }





    template<typename TypeTag>
    bool
    StandardWell<TypeTag>::
    iterateWellEqWithControl(const Simulator& ebosSimulator,
                             const double dt,
                             const Well::InjectionControls& inj_controls,
                             const Well::ProductionControls& prod_controls,
                             WellState& well_state,
                             const GroupState& group_state,
                             DeferredLogger& deferred_logger)
    {
        const int max_iter = this->param_.max_inner_iter_wells_;
        int it = 0;
        bool converged;
        bool relax_convergence = false;
        this->regularize_ = false;
        do {
            assembleWellEqWithoutIteration(ebosSimulator, dt, inj_controls, prod_controls, well_state, group_state, deferred_logger);

            if (it > this->param_.strict_inner_iter_wells_) {
                relax_convergence = true;
                this->regularize_ = true;
            }

            auto report = getWellConvergence(well_state, Base::B_avg_, deferred_logger, relax_convergence);

            converged = report.converged();
            if (converged) {
                break;
            }

            ++it;
            solveEqAndUpdateWellState(well_state, deferred_logger);

            // TODO: when this function is used for well testing purposes, will need to check the controls, so that we will obtain convergence
            // under the most restrictive control. Based on this converged results, we can check whether to re-open the well. Either we refactor
            // this function or we use different functions for the well testing purposes.
            // We don't allow for switching well controls while computing well potentials and testing wells
            // updateWellControl(ebosSimulator, well_state, deferred_logger);
            initPrimaryVariablesEvaluation();
        } while (it < max_iter);

        return converged;
    }


    template<typename TypeTag>
    std::vector<double>
    StandardWell<TypeTag>::
    computeCurrentWellRates(const Simulator& ebosSimulator,
                            DeferredLogger& deferred_logger) const
    {
        // Calculate the rates that follow from the current primary variables.
        std::vector<double> well_q_s(this->num_components_, 0.);
        const EvalWell& bhp = this->primary_variables_.eval(Bhp);
        const bool allow_cf = this->getAllowCrossFlow() || openCrossFlowAvoidSingularity(ebosSimulator);
        for (int perf = 0; perf < this->number_of_perforations_; ++perf) {
            const int cell_idx = this->well_cells_[perf];
            const auto& intQuants = *(ebosSimulator.model().cachedIntensiveQuantities(cell_idx, /*timeIdx=*/ 0));
            std::vector<Scalar> mob(this->num_components_, 0.);
            getMobilityScalar(ebosSimulator, perf, mob, deferred_logger);
            std::vector<Scalar> cq_s(this->num_components_, 0.);
            double trans_mult = ebosSimulator.problem().template rockCompTransMultiplier<double>(intQuants,  cell_idx);
            const double Tw = this->well_index_[perf] * trans_mult;
            computePerfRateScalar(intQuants, mob, bhp.value(), Tw, perf, allow_cf,
                            cq_s, deferred_logger);
            for (int comp = 0; comp < this->num_components_; ++comp) {
                well_q_s[comp] += cq_s[comp];
            }
        }
        const auto& comm = this->parallel_well_info_.communication();
        if (comm.size() > 1)
        {
            comm.sum(well_q_s.data(), well_q_s.size());
        }
        return well_q_s;
    }





    template <typename TypeTag>
    void
    StandardWell<TypeTag>::
    computeConnLevelProdInd(const typename StandardWell<TypeTag>::FluidState& fs,
                            const std::function<double(const double)>& connPICalc,
                            const std::vector<EvalWell>& mobility,
                            double* connPI) const
    {
        const auto& pu = this->phaseUsage();
        const int   np = this->number_of_phases_;
        for (int p = 0; p < np; ++p) {
            // Note: E100's notion of PI value phase mobility includes
            // the reciprocal FVF.
            const auto connMob =
                mobility[ this->flowPhaseToEbosCompIdx(p) ].value()
                    * fs.invB(this->flowPhaseToEbosPhaseIdx(p)).value();

            connPI[p] = connPICalc(connMob);
        }

        if (FluidSystem::phaseIsActive(FluidSystem::oilPhaseIdx) &&
            FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx))
        {
            const auto io = pu.phase_pos[Oil];
            const auto ig = pu.phase_pos[Gas];

            const auto vapoil = connPI[ig] * fs.Rv().value();
            const auto disgas = connPI[io] * fs.Rs().value();

            connPI[io] += vapoil;
            connPI[ig] += disgas;
        }
    }





    template <typename TypeTag>
    void
    StandardWell<TypeTag>::
    computeConnLevelInjInd(const typename StandardWell<TypeTag>::FluidState& fs,
                           const Phase preferred_phase,
                           const std::function<double(const double)>& connIICalc,
                           const std::vector<EvalWell>& mobility,
                           double* connII,
                           DeferredLogger& deferred_logger) const
    {
        // Assumes single phase injection
        const auto& pu = this->phaseUsage();

        auto phase_pos = 0;
        if (preferred_phase == Phase::GAS) {
            phase_pos = pu.phase_pos[Gas];
        }
        else if (preferred_phase == Phase::OIL) {
            phase_pos = pu.phase_pos[Oil];
        }
        else if (preferred_phase == Phase::WATER) {
            phase_pos = pu.phase_pos[Water];
        }
        else {
            OPM_DEFLOG_THROW(NotImplemented,
                             "Unsupported Injector Type ("
                             << static_cast<int>(preferred_phase)
                             << ") for well " << this->name()
                             << " during connection I.I. calculation",
                             deferred_logger);
        }

        const auto zero   = EvalWell{this->primary_variables_.numWellEq() + Indices::numEq, 0.0};
        const auto mt     = std::accumulate(mobility.begin(), mobility.end(), zero);
        connII[phase_pos] = connIICalc(mt.value() * fs.invB(this->flowPhaseToEbosPhaseIdx(phase_pos)).value());
    }
} // namespace Opm
