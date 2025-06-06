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

#ifndef OPM_STANDARDWELL_IMPL_HEADER_INCLUDED
#define OPM_STANDARDWELL_IMPL_HEADER_INCLUDED

// Improve IDE experience
#ifndef OPM_STANDARDWELL_HEADER_INCLUDED
#include <config.h>
#include <opm/simulators/wells/StandardWell.hpp>
#endif

#include <opm/common/Exceptions.hpp>

#include <opm/input/eclipse/Units/Units.hpp>

#include <opm/simulators/utils/DeferredLoggingErrorHelpers.hpp>
#include <opm/simulators/wells/StandardWellAssemble.hpp>
#include <opm/simulators/wells/VFPHelpers.hpp>
#include <opm/simulators/wells/WellBhpThpCalculator.hpp>
#include <opm/simulators/wells/WellConvergence.hpp>

#include <algorithm>
#include <cstddef>
#include <functional>

#include <fmt/format.h>

namespace Opm
{

    template<typename TypeTag>
    StandardWell<TypeTag>::
    StandardWell(const Well& well,
                 const ParallelWellInfo<Scalar>& pw_info,
                 const int time_step,
                 const ModelParameters& param,
                 const RateConverterType& rate_converter,
                 const int pvtRegionIdx,
                 const int num_components,
                 const int num_phases,
                 const int index_of_well,
                 const std::vector<PerforationData<Scalar>>& perf_data)
    : Base(well, pw_info, time_step, param, rate_converter, pvtRegionIdx, num_components, num_phases, index_of_well, perf_data)
    , StdWellEval(static_cast<const WellInterfaceIndices<FluidSystem,Indices>&>(*this))
    , regularize_(false)
    {
        assert(this->num_components_ == numWellConservationEq);
    }





    template<typename TypeTag>
    void
    StandardWell<TypeTag>::
    init(const PhaseUsage* phase_usage_arg,
         const std::vector<Scalar>& depth_arg,
         const Scalar gravity_arg,
         const std::vector< Scalar >& B_avg,
         const bool changed_to_open_this_step)
    {
        Base::init(phase_usage_arg, depth_arg, gravity_arg, B_avg, changed_to_open_this_step);
        this->StdWellEval::init(this->perf_depth_, depth_arg, Base::has_polymermw);
    }





    template<typename TypeTag>
    template<class Value>
    void
    StandardWell<TypeTag>::
    computePerfRate(const IntensiveQuantities& intQuants,
                    const std::vector<Value>& mob,
                    const Value& bhp,
                    const std::vector<Scalar>& Tw,
                    const int perf,
                    const bool allow_cf,
                    std::vector<Value>& cq_s,
                    PerforationRates<Scalar>& perf_rates,
                    DeferredLogger& deferred_logger) const
    {
        auto obtain = [this](const Eval& value)
                      {
                          if constexpr (std::is_same_v<Value, Scalar>) {
                              static_cast<void>(this); // suppress clang warning
                              return getValue(value);
                          } else {
                              return this->extendEval(value);
                          }
                      };
        auto obtainN = [](const auto& value)
        {
            if constexpr (std::is_same_v<Value, Scalar>) {
                return getValue(value);
            } else {
                return value;
            }
        };
        auto zeroElem = [this]()
                        {
                            if constexpr (std::is_same_v<Value, Scalar>) {
                                static_cast<void>(this); // suppress clang warning
                                return 0.0;
                            } else {
                                return Value{this->primary_variables_.numWellEq() + Indices::numEq, 0.0};
                            }
                        };

        const auto& fs = intQuants.fluidState();
        const Value pressure = obtain(this->getPerfCellPressure(fs));
        const Value rs = obtain(fs.Rs());
        const Value rv = obtain(fs.Rv());
        const Value rvw = obtain(fs.Rvw());
        const Value rsw = obtain(fs.Rsw());

        std::vector<Value> b_perfcells_dense(this->numComponents(), zeroElem());
        for (unsigned phaseIdx = 0; phaseIdx < FluidSystem::numPhases; ++phaseIdx) {
            if (!FluidSystem::phaseIsActive(phaseIdx)) {
                continue;
            }
            const unsigned compIdx = Indices::canonicalToActiveComponentIndex(FluidSystem::solventComponentIndex(phaseIdx));
            b_perfcells_dense[compIdx] =  obtain(fs.invB(phaseIdx));
        }
        if constexpr (has_solvent) {
            b_perfcells_dense[Indices::contiSolventEqIdx] = obtain(intQuants.solventInverseFormationVolumeFactor());
        }

        if constexpr (has_zFraction) {
            if (this->isInjector()) {
                const unsigned gasCompIdx = Indices::canonicalToActiveComponentIndex(FluidSystem::gasCompIdx);
                b_perfcells_dense[gasCompIdx] *= (1.0 - this->wsolvent());
                b_perfcells_dense[gasCompIdx] += this->wsolvent()*intQuants.zPureInvFormationVolumeFactor().value();
            }
        }

        Value skin_pressure = zeroElem();
        if (has_polymermw) {
            if (this->isInjector()) {
                const int pskin_index = Bhp + 1 + this->numLocalPerfs() + perf;
                skin_pressure = obtainN(this->primary_variables_.eval(pskin_index));
            }
        }

        // surface volume fraction of fluids within wellbore
        std::vector<Value> cmix_s(this->numComponents(), zeroElem());
        for (int componentIdx = 0; componentIdx < this->numComponents(); ++componentIdx) {
            cmix_s[componentIdx] = obtainN(this->primary_variables_.surfaceVolumeFraction(componentIdx));
        }

        computePerfRate(mob,
                        pressure,
                        bhp,
                        rs,
                        rv,
                        rvw,
                        rsw,
                        b_perfcells_dense,
                        Tw,
                        perf,
                        allow_cf,
                        skin_pressure,
                        cmix_s,
                        cq_s,
                        perf_rates,
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
                    const Value& rsw,
                    std::vector<Value>& b_perfcells_dense,
                    const std::vector<Scalar>& Tw,
                    const int perf,
                    const bool allow_cf,
                    const Value& skin_pressure,
                    const std::vector<Value>& cmix_s,
                    std::vector<Value>& cq_s,
                    PerforationRates<Scalar>& perf_rates,
                    DeferredLogger& deferred_logger) const
    {
        // Pressure drawdown (also used to determine direction of flow)
        const Value well_pressure = bhp + this->connections_.pressure_diff(perf);
        Value drawdown = pressure - well_pressure;
        if (this->isInjector()) {
            drawdown += skin_pressure;
        }

        RatioCalculator<Value> ratioCalc{
            FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx)
                ? Indices::canonicalToActiveComponentIndex(FluidSystem::gasCompIdx)
                : -1,
            FluidSystem::phaseIsActive(FluidSystem::oilPhaseIdx)
                ? Indices::canonicalToActiveComponentIndex(FluidSystem::oilCompIdx)
                : -1,
            FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx)
                ? Indices::canonicalToActiveComponentIndex(FluidSystem::waterCompIdx)
                : -1,
            this->name()
        };

        // producing perforations
        if (drawdown > 0)  {
            // Do nothing if crossflow is not allowed
            if (!allow_cf && this->isInjector()) {
                return;
            }

            // compute component volumetric rates at standard conditions
            for (int componentIdx = 0; componentIdx < this->numComponents(); ++componentIdx) {
                const Value cq_p = - Tw[componentIdx] * (mob[componentIdx] * drawdown);
                cq_s[componentIdx] = b_perfcells_dense[componentIdx] * cq_p;
            }

            if (FluidSystem::phaseIsActive(FluidSystem::oilPhaseIdx) &&
                FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx))
            {
                ratioCalc.gasOilPerfRateProd(cq_s, perf_rates, rv, rs, rvw,
                                             FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx),
                                             this->isProducer());
            } else if (FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx) &&
                       FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx))
            {
                ratioCalc.gasWaterPerfRateProd(cq_s, perf_rates, rvw, rsw, this->isProducer());
            }
        } else {
            // Do nothing if crossflow is not allowed
            if (!allow_cf && this->isProducer()) {
                return;
            }

            // Using total mobilities
            Value total_mob_dense = mob[0];
            for (int componentIdx = 1; componentIdx < this->numComponents(); ++componentIdx) {
                total_mob_dense += mob[componentIdx];
            }

            // compute volume ratio between connection at standard conditions
            Value volumeRatio = bhp * 0.0; // initialize it with the correct type

            if (FluidSystem::enableVaporizedWater() && FluidSystem::enableDissolvedGasInWater()) {
                ratioCalc.disOilVapWatVolumeRatio(volumeRatio, rvw, rsw, pressure,
                                                  cmix_s, b_perfcells_dense, deferred_logger);
                // DISGASW only supported for gas-water CO2STORE/H2STORE case
                // and the simulator will throw long before it reach to this point in the code
                // For blackoil support of DISGASW we need to add the oil component here
                assert(FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx));
                assert(FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx));
                assert(!FluidSystem::phaseIsActive(FluidSystem::oilPhaseIdx));
            } else {

                if (FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx)) {
                    const unsigned waterCompIdx = Indices::canonicalToActiveComponentIndex(FluidSystem::waterCompIdx);
                    volumeRatio += cmix_s[waterCompIdx] / b_perfcells_dense[waterCompIdx];
                }

                if constexpr (Indices::enableSolvent) {
                    volumeRatio += cmix_s[Indices::contiSolventEqIdx] / b_perfcells_dense[Indices::contiSolventEqIdx];
                }

                if (FluidSystem::phaseIsActive(FluidSystem::oilPhaseIdx) &&
                    FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx))
                {
                    ratioCalc.gasOilVolumeRatio(volumeRatio, rv, rs, pressure,
                                                cmix_s, b_perfcells_dense,
                                                deferred_logger);
                } else {
                    if (FluidSystem::phaseIsActive(FluidSystem::oilPhaseIdx)) {
                        const unsigned oilCompIdx = Indices::canonicalToActiveComponentIndex(FluidSystem::oilCompIdx);
                        volumeRatio += cmix_s[oilCompIdx] / b_perfcells_dense[oilCompIdx];
                    }
                    if (FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx)) {
                        const unsigned gasCompIdx = Indices::canonicalToActiveComponentIndex(FluidSystem::gasCompIdx);
                        volumeRatio += cmix_s[gasCompIdx] / b_perfcells_dense[gasCompIdx];
                    }
                }
            }

            // injecting connections total volumerates at standard conditions
            for (int componentIdx = 0; componentIdx < this->numComponents(); ++componentIdx) {
                const Value cqt_i = - Tw[componentIdx] * (total_mob_dense * drawdown);
                Value cqt_is = cqt_i / volumeRatio;
                cq_s[componentIdx] = cmix_s[componentIdx] * cqt_is;
            }

            // calculating the perforation solution gas rate and solution oil rates
            if (this->isProducer()) {
                if (FluidSystem::phaseIsActive(FluidSystem::oilPhaseIdx) &&
                    FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx))
                {
                    ratioCalc.gasOilPerfRateInj(cq_s, perf_rates,
                                                rv, rs, pressure, rvw,
                                                FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx),
                                                deferred_logger);
                }
                if (FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx) &&
                    FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx))
                {
                    //no oil
                    ratioCalc.gasWaterPerfRateInj(cq_s, perf_rates, rvw, rsw,
                                                  pressure, deferred_logger);
                }
            }
        }
    }


    template<typename TypeTag>
    void
    StandardWell<TypeTag>::
    assembleWellEqWithoutIteration(const Simulator& simulator,
                                   const double dt,
                                   const Well::InjectionControls& inj_controls,
                                   const Well::ProductionControls& prod_controls,
                                   WellState<Scalar>& well_state,
                                   const GroupState<Scalar>& group_state,
                                   DeferredLogger& deferred_logger)
    {
        // TODO: only_wells should be put back to save some computation
        // for example, the matrices B C does not need to update if only_wells
        if (!this->isOperableAndSolvable() && !this->wellIsStopped()) return;

        // clear all entries
        this->linSys_.clear();

        assembleWellEqWithoutIterationImpl(simulator, dt, inj_controls,
                                           prod_controls, well_state,
                                           group_state, deferred_logger);
    }




    template<typename TypeTag>
    void
    StandardWell<TypeTag>::
    assembleWellEqWithoutIterationImpl(const Simulator& simulator,
                                       const double dt,
                                       const Well::InjectionControls& inj_controls,
                                       const Well::ProductionControls& prod_controls,
                                       WellState<Scalar>& well_state,
                                       const GroupState<Scalar>& group_state,
                                       DeferredLogger& deferred_logger)
    {
        // try to regularize equation if the well does not converge
        const Scalar regularization_factor =  this->regularize_? this->param_.regularization_factor_wells_ : 1.0;
        const Scalar volume = 0.1 * unit::cubic(unit::feet) * regularization_factor;

        auto& ws = well_state.well(this->index_of_well_);
        ws.phase_mixing_rates.fill(0.0);


        const int np = this->number_of_phases_;

        std::vector<RateVector> connectionRates = this->connectionRates_; // Copy to get right size.

        auto& perf_data = ws.perf_data;
        auto& perf_rates = perf_data.phase_rates;
        for (int perf = 0; perf < this->number_of_local_perforations_; ++perf) {
            // Calculate perforation quantities.
            std::vector<EvalWell> cq_s(this->num_components_, {this->primary_variables_.numWellEq() + Indices::numEq, 0.0});
            EvalWell water_flux_s{this->primary_variables_.numWellEq() + Indices::numEq, 0.0};
            EvalWell cq_s_zfrac_effective{this->primary_variables_.numWellEq() + Indices::numEq, 0.0};
            calculateSinglePerf(simulator, perf, well_state, connectionRates,
                                cq_s, water_flux_s, cq_s_zfrac_effective, deferred_logger);

            // Equation assembly for this perforation.
            if constexpr (has_polymer && Base::has_polymermw) {
                if (this->isInjector()) {
                    handleInjectivityEquations(simulator, well_state, perf,
                                               water_flux_s, deferred_logger);
                }
            }
            for (int componentIdx = 0; componentIdx < this->num_components_; ++componentIdx) {
                // the cq_s entering mass balance equations need to consider the efficiency factors.
                const EvalWell cq_s_effective = cq_s[componentIdx] * this->well_efficiency_factor_;

                connectionRates[perf][componentIdx] = Base::restrictEval(cq_s_effective);

                StandardWellAssemble<FluidSystem,Indices>(*this).
                    assemblePerforationEq(cq_s_effective,
                                          componentIdx,
                                          perf,
                                          this->primary_variables_.numWellEq(),
                                          this->linSys_);

                // Store the perforation phase flux for later usage.
                if (has_solvent && componentIdx == Indices::contiSolventEqIdx) {
                    auto& perf_rate_solvent = perf_data.solvent_rates;
                    perf_rate_solvent[perf] = cq_s[componentIdx].value();
                } else {
                    perf_rates[perf*np + this->modelCompIdxToFlowCompIdx(componentIdx)] = cq_s[componentIdx].value();
                }
            }

            if constexpr (has_zFraction) {
                StandardWellAssemble<FluidSystem,Indices>(*this).
                    assembleZFracEq(cq_s_zfrac_effective,
                                    perf,
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
            comm.sum(ws.phase_mixing_rates.data(), ws.phase_mixing_rates.size());
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
            StandardWellAssemble<FluidSystem,Indices>(*this).
                assembleSourceEq(resWell_loc,
                                 componentIdx,
                                 this->primary_variables_.numWellEq(),
                                 this->linSys_);
        }

        const auto& summaryState = simulator.vanguard().summaryState();
        const Schedule& schedule = simulator.vanguard().schedule();
        const bool stopped_or_zero_target = this->stoppedOrZeroRateTarget(simulator, well_state, deferred_logger);
        StandardWellAssemble<FluidSystem,Indices>(*this).
            assembleControlEq(well_state, group_state,
                              schedule, summaryState,
                              inj_controls, prod_controls,
                              this->primary_variables_,
                              this->connections_.rho(),
                              this->linSys_,
                              stopped_or_zero_target,
                              deferred_logger);


        // do the local inversion of D.
        try {
            this->linSys_.invert();
        } catch( ... ) {
            OPM_DEFLOG_PROBLEM(NumericalProblem, "Error when inverting local well equations for well " + name(), deferred_logger);
        }
    }




    template<typename TypeTag>
    void
    StandardWell<TypeTag>::
    calculateSinglePerf(const Simulator& simulator,
                        const int perf,
                        WellState<Scalar>& well_state,
                        std::vector<RateVector>& connectionRates,
                        std::vector<EvalWell>& cq_s,
                        EvalWell& water_flux_s,
                        EvalWell& cq_s_zfrac_effective,
                        DeferredLogger& deferred_logger) const
    {
        const bool allow_cf = this->getAllowCrossFlow() || openCrossFlowAvoidSingularity(simulator);
        const EvalWell& bhp = this->primary_variables_.eval(Bhp);
        const int cell_idx = this->well_cells_[perf];
        const auto& intQuants = simulator.model().intensiveQuantities(cell_idx, /*timeIdx=*/ 0);
        std::vector<EvalWell> mob(this->num_components_, {this->primary_variables_.numWellEq() + Indices::numEq, 0.});
        getMobility(simulator, perf, mob, deferred_logger);

        PerforationRates<Scalar> perf_rates;
        Scalar trans_mult = simulator.problem().template wellTransMultiplier<Scalar>(intQuants,  cell_idx);
        const auto& wellstate_nupcol = simulator.problem().wellModel().nupcolWellState().well(this->index_of_well_);
        const std::vector<Scalar> Tw = this->wellIndex(perf, intQuants, trans_mult, wellstate_nupcol);
        computePerfRate(intQuants, mob, bhp, Tw, perf, allow_cf,
                        cq_s, perf_rates, deferred_logger);

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
                handleInjectivityRate(simulator, perf, cq_s);
            }
        }

        // updating the solution gas rate and solution oil rate
        if (this->isProducer()) {
            ws.phase_mixing_rates[ws.dissolved_gas] += perf_rates.dis_gas;
            ws.phase_mixing_rates[ws.dissolved_gas_in_water] += perf_rates.dis_gas_in_water;
            ws.phase_mixing_rates[ws.vaporized_oil] += perf_rates.vap_oil;
            ws.phase_mixing_rates[ws.vaporized_water] += perf_rates.vap_wat;
            perf_data.phase_mixing_rates[perf][ws.dissolved_gas] = perf_rates.dis_gas;
            perf_data.phase_mixing_rates[perf][ws.dissolved_gas_in_water] = perf_rates.dis_gas_in_water;
            perf_data.phase_mixing_rates[perf][ws.vaporized_oil] = perf_rates.vap_oil;
            perf_data.phase_mixing_rates[perf][ws.vaporized_water] = perf_rates.vap_wat;
        }

        if constexpr (has_energy) {
            connectionRates[perf][Indices::contiEnergyEqIdx] =
                connectionRateEnergy(simulator.problem().maxOilSaturation(cell_idx),
                                     cq_s, intQuants, deferred_logger);
        }

        if constexpr (has_polymer) {
            std::variant<Scalar,EvalWell> polymerConcentration;
            if (this->isInjector()) {
                polymerConcentration = this->wpolymer();
            } else {
                polymerConcentration = this->extendEval(intQuants.polymerConcentration() *
                                                        intQuants.polymerViscosityCorrection());
            }

            [[maybe_unused]] EvalWell cq_s_poly;
            std::tie(connectionRates[perf][Indices::contiPolymerEqIdx],
                     cq_s_poly) =
                this->connections_.connectionRatePolymer(perf_data.polymer_rates[perf],
                                                         cq_s, polymerConcentration);

            if constexpr (Base::has_polymermw) {
                updateConnectionRatePolyMW(cq_s_poly, intQuants, well_state,
                                           perf, connectionRates, deferred_logger);
            }
        }

        if constexpr (has_foam) {
            std::variant<Scalar,EvalWell> foamConcentration;
            if (this->isInjector()) {
                foamConcentration = this->wfoam();
            } else {
                foamConcentration = this->extendEval(intQuants.foamConcentration());
            }
            connectionRates[perf][Indices::contiFoamEqIdx] =
                this->connections_.connectionRateFoam(cq_s, foamConcentration,
                                                      FoamModule::transportPhase(),
                                                      deferred_logger);
        }

        if constexpr (has_zFraction) {
            std::variant<Scalar,std::array<EvalWell,2>> solventConcentration;
            if (this->isInjector()) {
                solventConcentration = this->wsolvent();
            } else {
                solventConcentration = std::array{this->extendEval(intQuants.xVolume()),
                                                  this->extendEval(intQuants.yVolume())};
            }
            std::tie(connectionRates[perf][Indices::contiZfracEqIdx],
                     cq_s_zfrac_effective) =
                this->connections_.connectionRatezFraction(perf_data.solvent_rates[perf],
                                                           perf_rates.dis_gas, cq_s,
                                                           solventConcentration);
        }

        if constexpr (has_brine) {
            std::variant<Scalar,EvalWell> saltConcentration;
            if (this->isInjector()) {
                saltConcentration = this->wsalt();
            } else {
                saltConcentration = this->extendEval(intQuants.fluidState().saltConcentration());
            }

            connectionRates[perf][Indices::contiBrineEqIdx] =
                this->connections_.connectionRateBrine(perf_data.brine_rates[perf],
                                                       perf_rates.vap_wat, cq_s,
                                                       saltConcentration);
        }

        if constexpr (has_micp) {
            std::variant<Scalar,EvalWell> microbialConcentration;
            std::variant<Scalar,EvalWell> oxygenConcentration;
            std::variant<Scalar,EvalWell> ureaConcentration;
            if (this->isInjector()) {
                microbialConcentration = this->wmicrobes();
                oxygenConcentration = this->woxygen();
                ureaConcentration = this->wurea();
            } else {
                microbialConcentration = this->extendEval(intQuants.microbialConcentration());
                oxygenConcentration = this->extendEval(intQuants.oxygenConcentration());
                ureaConcentration = this->extendEval(intQuants.ureaConcentration());
            }
            std::tie(connectionRates[perf][Indices::contiMicrobialEqIdx],
                     connectionRates[perf][Indices::contiOxygenEqIdx],
                     connectionRates[perf][Indices::contiUreaEqIdx]) =
                this->connections_.connectionRatesMICP(perf_data.microbial_rates[perf],
                                                       perf_data.oxygen_rates[perf],
                                                       perf_data.urea_rates[perf],
                                                       cq_s,
                                                       microbialConcentration,
                                                       oxygenConcentration,
                                                       ureaConcentration);
        }

        // Store the perforation pressure for later usage.
        perf_data.pressure[perf] = ws.bhp + this->connections_.pressure_diff(perf);

        // Store the perforation gass mass rate.
        const auto& pu = well_state.phaseUsage();
        if (pu.has_co2_or_h2store) {
            const unsigned gas_comp_idx = Indices::canonicalToActiveComponentIndex(FluidSystem::gasCompIdx);
            const Scalar rho = FluidSystem::referenceDensity( FluidSystem::gasPhaseIdx, Base::pvtRegionIdx() );
            perf_data.gas_mass_rates[perf] = cq_s[gas_comp_idx].value() * rho;
        }

        // Store the perforation water mass rate.
        if (FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx)) {
            const unsigned wat_comp_idx = Indices::canonicalToActiveComponentIndex(FluidSystem::waterCompIdx);
            const Scalar rho = FluidSystem::referenceDensity( FluidSystem::waterPhaseIdx, Base::pvtRegionIdx() );
            perf_data.wat_mass_rates[perf] = cq_s[wat_comp_idx].value() * rho;
        }
    }



    template<typename TypeTag>
    template<class Value>
    void
    StandardWell<TypeTag>::
    getMobility(const Simulator& simulator,
                const int perf,
                std::vector<Value>& mob,
                DeferredLogger& deferred_logger) const
    {
        auto obtain = [this](const Eval& value)
                      {
                          if constexpr (std::is_same_v<Value, Scalar>) {
                              static_cast<void>(this); // suppress clang warning
                              return getValue(value);
                          } else {
                              return this->extendEval(value);
                          }
                      };
        WellInterface<TypeTag>::getMobility(simulator, perf, mob,
                                            obtain, deferred_logger);

        // modify the water mobility if polymer is present
        if constexpr (has_polymer) {
            if (!FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx)) {
                OPM_DEFLOG_THROW(std::runtime_error, "Water is required when polymer is active", deferred_logger);
            }

            // for the cases related to polymer molecular weight, we assume fully mixing
            // as a result, the polymer and water share the same viscosity
            if constexpr (!Base::has_polymermw) {
                if constexpr (std::is_same_v<Value, Scalar>) {
                    std::vector<EvalWell> mob_eval(this->num_components_, {this->primary_variables_.numWellEq() + Indices::numEq, 0.});
                    for (std::size_t i = 0; i < mob.size(); ++i) {
                        mob_eval[i].setValue(mob[i]);
                    }
                    updateWaterMobilityWithPolymer(simulator, perf, mob_eval, deferred_logger);
                    for (std::size_t i = 0; i < mob.size(); ++i) {
                        mob[i] = getValue(mob_eval[i]);
                    }
                } else {
                    updateWaterMobilityWithPolymer(simulator, perf, mob, deferred_logger);
                }
            }
        }

        // if the injecting well has WINJMULT setup, we update the mobility accordingly
        if (this->isInjector() && this->well_ecl_.getInjMultMode() != Well::InjMultMode::NONE) {
            const Scalar bhp = this->primary_variables_.value(Bhp);
            const Scalar perf_press = bhp +  this->connections_.pressure_diff(perf);
            const Scalar multiplier = this->getInjMult(perf, bhp, perf_press, deferred_logger);
            for (std::size_t i = 0; i < mob.size(); ++i) {
                mob[i] *= multiplier;
            }
        }
    }


    template<typename TypeTag>
    void
    StandardWell<TypeTag>::
    updateWellState(const Simulator& simulator,
                    const BVectorWell& dwells,
                    WellState<Scalar>& well_state,
                    DeferredLogger& deferred_logger)
    {
        if (!this->isOperableAndSolvable() && !this->wellIsStopped()) return;

        const bool stop_or_zero_rate_target = this->stoppedOrZeroRateTarget(simulator, well_state, deferred_logger);
        updatePrimaryVariablesNewton(dwells, stop_or_zero_rate_target, deferred_logger);

        const auto& summary_state = simulator.vanguard().summaryState();
        updateWellStateFromPrimaryVariables(well_state, summary_state, deferred_logger);
        Base::calculateReservoirRates(simulator.vanguard().eclState().runspec().co2Storage(), well_state.well(this->index_of_well_));
    }





    template<typename TypeTag>
    void
    StandardWell<TypeTag>::
    updatePrimaryVariablesNewton(const BVectorWell& dwells,
                                 const bool stop_or_zero_rate_target,
                                 DeferredLogger& deferred_logger)
    {
        const Scalar dFLimit = this->param_.dwell_fraction_max_;
        const Scalar dBHPLimit = this->param_.dbhp_max_rel_;
        this->primary_variables_.updateNewton(dwells, stop_or_zero_rate_target, dFLimit, dBHPLimit, deferred_logger);

        // for the water velocity and skin pressure
        if constexpr (Base::has_polymermw) {
            this->primary_variables_.updateNewtonPolyMW(dwells);
        }

        this->primary_variables_.checkFinite(deferred_logger);
    }





    template<typename TypeTag>
    void
    StandardWell<TypeTag>::
    updateWellStateFromPrimaryVariables(WellState<Scalar>& well_state,
                                        const SummaryState& summary_state,
                                        DeferredLogger& deferred_logger) const
    {
        this->StdWellEval::updateWellStateFromPrimaryVariables(well_state, summary_state, deferred_logger);

        // other primary variables related to polymer injectivity study
        if constexpr (Base::has_polymermw) {
            this->primary_variables_.copyToWellStatePolyMW(well_state);
        }
    }





    template<typename TypeTag>
    void
    StandardWell<TypeTag>::
    updateIPR(const Simulator& simulator, DeferredLogger& deferred_logger) const
    {
        // TODO: not handling solvent related here for now

        // initialize all the values to be zero to begin with
        std::fill(this->ipr_a_.begin(), this->ipr_a_.end(), 0.);
        std::fill(this->ipr_b_.begin(), this->ipr_b_.end(), 0.);

        for (int perf = 0; perf < this->number_of_local_perforations_; ++perf) {
            std::vector<Scalar> mob(this->num_components_, 0.0);
            getMobility(simulator, perf, mob, deferred_logger);

            const int cell_idx = this->well_cells_[perf];
            const auto& int_quantities = simulator.model().intensiveQuantities(cell_idx, /*timeIdx=*/ 0);
            const auto& fs = int_quantities.fluidState();
            // the pressure of the reservoir grid block the well connection is in
            Scalar p_r = this->getPerfCellPressure(fs).value();

            // calculating the b for the connection
            std::vector<Scalar> b_perf(this->num_components_);
            for (std::size_t phase = 0; phase < FluidSystem::numPhases; ++phase) {
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
            const Scalar h_perf = this->connections_.pressure_diff(perf);
            const Scalar pressure_diff = p_r - h_perf;

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
            Scalar trans_mult = simulator.problem().template wellTransMultiplier<Scalar>(int_quantities, cell_idx);
            const auto& wellstate_nupcol = simulator.problem().wellModel().nupcolWellState().well(this->index_of_well_);
            const std::vector<Scalar> tw_perf = this->wellIndex(perf,
                                                                int_quantities,
                                                                trans_mult,
                                                                wellstate_nupcol);
            std::vector<Scalar> ipr_a_perf(this->ipr_a_.size());
            std::vector<Scalar> ipr_b_perf(this->ipr_b_.size());
            for (int comp_idx = 0; comp_idx < this->num_components_; ++comp_idx) {
                const Scalar tw_mob = tw_perf[comp_idx] * mob[comp_idx] * b_perf[comp_idx];
                ipr_a_perf[comp_idx] += tw_mob * pressure_diff;
                ipr_b_perf[comp_idx] += tw_mob;
            }

            // we need to handle the rs and rv when both oil and gas are present
            if (FluidSystem::phaseIsActive(FluidSystem::oilPhaseIdx) && FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx)) {
                const unsigned oil_comp_idx = Indices::canonicalToActiveComponentIndex(FluidSystem::oilCompIdx);
                const unsigned gas_comp_idx = Indices::canonicalToActiveComponentIndex(FluidSystem::gasCompIdx);
                const Scalar rs = (fs.Rs()).value();
                const Scalar rv = (fs.Rv()).value();

                const Scalar dis_gas_a = rs * ipr_a_perf[oil_comp_idx];
                const Scalar vap_oil_a = rv * ipr_a_perf[gas_comp_idx];

                ipr_a_perf[gas_comp_idx] += dis_gas_a;
                ipr_a_perf[oil_comp_idx] += vap_oil_a;

                const Scalar dis_gas_b = rs * ipr_b_perf[oil_comp_idx];
                const Scalar vap_oil_b = rv * ipr_b_perf[gas_comp_idx];

                ipr_b_perf[gas_comp_idx] += dis_gas_b;
                ipr_b_perf[oil_comp_idx] += vap_oil_b;
            }

            for (std::size_t comp_idx = 0; comp_idx < ipr_a_perf.size(); ++comp_idx) {
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
    updateIPRImplicit(const Simulator& simulator,
                      WellState<Scalar>& well_state,
                      DeferredLogger& deferred_logger)
    {
        // Compute IPR based on *converged* well-equation:
        // For a component rate r the derivative dr/dbhp is obtained by
        // dr/dbhp = - (partial r/partial x) * inv(partial Eq/partial x) * (partial Eq/partial bhp_target)
        // where Eq(x)=0 is the well equation setup with bhp control and primary variables x

        // We shouldn't have zero rates at this stage, but check
        bool zero_rates;
        auto rates = well_state.well(this->index_of_well_).surface_rates;
        zero_rates = true;
        for (std::size_t p = 0; p < rates.size(); ++p) {
            zero_rates &= rates[p] == 0.0;
        }
        auto& ws = well_state.well(this->index_of_well_);
        if (zero_rates) {
            const auto msg = fmt::format("updateIPRImplicit: Well {} has zero rate, IPRs might be problematic", this->name());
            deferred_logger.debug(msg);
            /*
            // could revert to standard approach here:
            updateIPR(simulator, deferred_logger);
            for (int comp_idx = 0; comp_idx < this->num_components_; ++comp_idx){
                const int idx = this->modelCompIdxToFlowCompIdx(comp_idx);
                ws.implicit_ipr_a[idx] = this->ipr_a_[comp_idx];
                ws.implicit_ipr_b[idx] = this->ipr_b_[comp_idx];
            }
            return;
            */
        }
        const auto& group_state  = simulator.problem().wellModel().groupState();

        std::fill(ws.implicit_ipr_a.begin(), ws.implicit_ipr_a.end(), 0.);
        std::fill(ws.implicit_ipr_b.begin(), ws.implicit_ipr_b.end(), 0.);

        auto inj_controls = Well::InjectionControls(0);
        auto prod_controls = Well::ProductionControls(0);
        prod_controls.addControl(Well::ProducerCMode::BHP);
        prod_controls.bhp_limit = well_state.well(this->index_of_well_).bhp;

        //  Set current control to bhp, and bhp value in state, modify bhp limit in control object.
        const auto cmode = ws.production_cmode;
        ws.production_cmode = Well::ProducerCMode::BHP;
        const double dt = simulator.timeStepSize();
        assembleWellEqWithoutIteration(simulator, dt, inj_controls, prod_controls, well_state, group_state, deferred_logger);

        const size_t nEq = this->primary_variables_.numWellEq();
        BVectorWell rhs(1);
        rhs[0].resize(nEq);
        // rhs = 0 except -1 for control eq
        for (size_t i=0; i < nEq; ++i){
            rhs[0][i] = 0.0;
        }
        rhs[0][Bhp] = -1.0;

        BVectorWell x_well(1);
        x_well[0].resize(nEq);
        this->linSys_.solve(rhs, x_well);

        for (int comp_idx = 0; comp_idx < this->num_components_; ++comp_idx){
            EvalWell comp_rate = this->primary_variables_.getQs(comp_idx);
            const int idx = this->modelCompIdxToFlowCompIdx(comp_idx);
            for (size_t pvIdx = 0; pvIdx < nEq; ++pvIdx) {
                // well primary variable derivatives in EvalWell start at position Indices::numEq
                ws.implicit_ipr_b[idx] -= x_well[0][pvIdx]*comp_rate.derivative(pvIdx+Indices::numEq);
            }
            ws.implicit_ipr_a[idx] = ws.implicit_ipr_b[idx]*ws.bhp - comp_rate.value();
        }
        // reset cmode
        ws.production_cmode = cmode;
    }

    template<typename TypeTag>
    void
    StandardWell<TypeTag>::
    checkOperabilityUnderBHPLimit(const WellState<Scalar>& well_state,
                                  const Simulator& simulator,
                                  DeferredLogger& deferred_logger)
    {
        const auto& summaryState = simulator.vanguard().summaryState();
        const Scalar bhp_limit = WellBhpThpCalculator(*this).mostStrictBhpFromBhpLimits(summaryState);
        // Crude but works: default is one atmosphere.
        // TODO: a better way to detect whether the BHP is defaulted or not
        const bool bhp_limit_not_defaulted = bhp_limit > 1.5 * unit::barsa;
        if ( bhp_limit_not_defaulted || !this->wellHasTHPConstraints(summaryState) ) {
            // if the BHP limit is not defaulted or the well does not have a THP limit
            // we need to check the BHP limit
            Scalar total_ipr_mass_rate = 0.0;
            for (unsigned phaseIdx = 0; phaseIdx < FluidSystem::numPhases; ++phaseIdx)
            {
                if (!FluidSystem::phaseIsActive(phaseIdx)) {
                    continue;
                }

                const unsigned compIdx = Indices::canonicalToActiveComponentIndex(FluidSystem::solventComponentIndex(phaseIdx));
                const Scalar ipr_rate = this->ipr_a_[compIdx] - this->ipr_b_[compIdx] * bhp_limit;

                const Scalar rho = FluidSystem::referenceDensity( phaseIdx, Base::pvtRegionIdx() );
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
                std::vector<Scalar> well_rates_bhp_limit;
                computeWellRatesWithBhp(simulator, bhp_limit, well_rates_bhp_limit, deferred_logger);

                this->adaptRatesForVFP(well_rates_bhp_limit);
                const Scalar thp_limit = this->getTHPConstraint(summaryState);
                const Scalar thp = WellBhpThpCalculator(*this).calculateThpFromBhp(well_rates_bhp_limit,
                                                                                   bhp_limit,
                                                                                   this->connections_.rho(),
                                                                                   this->getALQ(well_state),
                                                                                   thp_limit,
                                                                                   deferred_logger);
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
    checkOperabilityUnderTHPLimit(const Simulator& simulator,
                                  const WellState<Scalar>& well_state,
                                  DeferredLogger& deferred_logger)
    {
        const auto& summaryState = simulator.vanguard().summaryState();
        const auto obtain_bhp = this->isProducer() ? computeBhpAtThpLimitProd(well_state, simulator, summaryState, deferred_logger)
        : computeBhpAtThpLimitInj(simulator, summaryState, deferred_logger);

        if (obtain_bhp) {
            this->operability_status_.can_obtain_bhp_with_thp_limit = true;

            const Scalar bhp_limit = WellBhpThpCalculator(*this).mostStrictBhpFromBhpLimits(summaryState);
            this->operability_status_.obey_bhp_limit_with_thp_limit = this->isProducer() ?
                                                              *obtain_bhp >= bhp_limit : *obtain_bhp <= bhp_limit ;

            const Scalar thp_limit = this->getTHPConstraint(summaryState);
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
                const Scalar thp_limit = this->getTHPConstraint(summaryState);
                deferred_logger.debug(" could not find bhp value at thp limit "
                                      + std::to_string(unit::convert::to(thp_limit, unit::barsa))
                                      + " bar for well " + name() + ", the well might need to be closed ");
            }
        }
    }





    template<typename TypeTag>
    bool
    StandardWell<TypeTag>::
    allDrawDownWrongDirection(const Simulator& simulator) const
    {
        bool all_drawdown_wrong_direction = true;

        for (int perf = 0; perf < this->number_of_local_perforations_; ++perf) {
            const int cell_idx = this->well_cells_[perf];
            const auto& intQuants = simulator.model().intensiveQuantities(cell_idx, /*timeIdx=*/0);
            const auto& fs = intQuants.fluidState();

            const Scalar pressure = this->getPerfCellPressure(fs).value();
            const Scalar bhp = this->primary_variables_.eval(Bhp).value();

            // Pressure drawdown (also used to determine direction of flow)
            const Scalar well_pressure = bhp + this->connections_.pressure_diff(perf);
            const Scalar drawdown = pressure - well_pressure;

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
    openCrossFlowAvoidSingularity(const Simulator& simulator) const
    {
        return !this->getAllowCrossFlow() && allDrawDownWrongDirection(simulator);
    }




    template<typename TypeTag>
    typename StandardWell<TypeTag>::WellConnectionProps
    StandardWell<TypeTag>::
    computePropertiesForWellConnectionPressures(const Simulator&         simulator,
                                                const WellState<Scalar>& well_state) const
    {
        auto prop_func = typename StdWellEval::StdWellConnections::PressurePropertyFunctions {
            // getTemperature
            [&model = simulator.model()](int cell_idx, int phase_idx)
            {
                return model.intensiveQuantities(cell_idx, /* time_idx = */ 0)
                    .fluidState().temperature(phase_idx).value();
            },

            // getSaltConcentration
            [&model = simulator.model()](int cell_idx)
            {
                return model.intensiveQuantities(cell_idx, /* time_idx = */ 0)
                    .fluidState().saltConcentration().value();
            },

            // getPvtRegionIdx
            [&model = simulator.model()](int cell_idx)
            {
                return model.intensiveQuantities(cell_idx, /* time_idx = */ 0)
                    .fluidState().pvtRegionIndex();
            }
        };

        if constexpr (Indices::enableSolvent) {
            prop_func.solventInverseFormationVolumeFactor =
                [&model = simulator.model()](int cell_idx)
            {
                return model.intensiveQuantities(cell_idx, /* time_idx = */ 0)
                    .solventInverseFormationVolumeFactor().value();
            };

            prop_func.solventRefDensity = [&model = simulator.model()](int cell_idx)
            {
                return model.intensiveQuantities(cell_idx, /* time_idx = */ 0)
                    .solventRefDensity();
            };
        }

        return this->connections_.computePropertiesForPressures(well_state, prop_func);
    }





    template<typename TypeTag>
    ConvergenceReport
    StandardWell<TypeTag>::
    getWellConvergence(const Simulator& simulator,
                       const WellState<Scalar>& well_state,
                       const std::vector<Scalar>& B_avg,
                       DeferredLogger& deferred_logger,
                       const bool relax_tolerance) const
    {
        // the following implementation assume that the polymer is always after the w-o-g phases
        // For the polymer, energy and foam cases, there is one more mass balance equations of reservoir than wells
        assert((int(B_avg.size()) == this->num_components_) || has_polymer || has_energy || has_foam || has_brine || has_zFraction || has_micp);

        Scalar tol_wells = this->param_.tolerance_wells_;
        // use stricter tolerance for stopped wells and wells under zero rate target control.
        constexpr Scalar stopped_factor = 1.e-4;
        // use stricter tolerance for dynamic thp to ameliorate network convergence
        constexpr Scalar dynamic_thp_factor = 1.e-1;
        if (this->stoppedOrZeroRateTarget(simulator, well_state, deferred_logger)) {
            tol_wells = tol_wells*stopped_factor;
        } else if (this->getDynamicThpLimit()) {
            tol_wells = tol_wells*dynamic_thp_factor;
        }

        std::vector<Scalar> res;
        ConvergenceReport report = this->StdWellEval::getWellConvergence(well_state,
                                                                         B_avg,
                                                                         this->param_.max_residual_allowed_,
                                                                         tol_wells,
                                                                         this->param_.relaxed_tolerance_flow_well_,
                                                                         relax_tolerance,
                                                                         this->wellIsStopped(),
                                                                         res,
                                                                         deferred_logger);

        checkConvergenceExtraEqs(res, report);

        return report;
    }





    template<typename TypeTag>
    void
    StandardWell<TypeTag>::
    updateProductivityIndex(const Simulator& simulator,
                            const WellProdIndexCalculator<Scalar>& wellPICalc,
                            WellState<Scalar>& well_state,
                            DeferredLogger& deferred_logger) const
    {
        auto fluidState = [&simulator, this](const int perf)
        {
            const auto cell_idx = this->well_cells_[perf];
            return simulator.model()
               .intensiveQuantities(cell_idx, /*timeIdx=*/ 0).fluidState();
        };

        const int np = this->number_of_phases_;
        auto setToZero = [np](Scalar* x) -> void
        {
            std::fill_n(x, np, 0.0);
        };

        auto addVector = [np](const Scalar* src, Scalar* dest) -> void
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

            auto connPICalc = [&wellPICalc, allPerfID](const Scalar mobility) -> Scalar
            {
                return wellPICalc.connectionProdIndStandard(allPerfID, mobility);
            };

            std::vector<Scalar> mob(this->num_components_, 0.0);
            getMobility(simulator, static_cast<int>(subsetPerfID), mob, deferred_logger);

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

        assert ((static_cast<int>(subsetPerfID) == this->number_of_local_perforations_) &&
                "Internal logic error in processing connections for PI/II");
    }



    template<typename TypeTag>
    void StandardWell<TypeTag>::
    computeWellConnectionDensitesPressures(const Simulator& simulator,
                                           const WellState<Scalar>& well_state,
                                           const WellConnectionProps& props,
                                           DeferredLogger& deferred_logger)
    {
        // Cell level dynamic property call-back functions as fall-back
        // option for calculating connection level mixture densities in
        // stopped or zero-rate producer wells.
        const auto prop_func = typename StdWellEval::StdWellConnections::DensityPropertyFunctions {
            // This becomes slightly more palatable with C++20's designated
            // initialisers.

            // mobility: Phase mobilities in specified cell.
            [&model = simulator.model()](const int               cell,
                                         const std::vector<int>& phases,
                                         std::vector<Scalar>&    mob)
            {
                const auto& iq = model.intensiveQuantities(cell, /* time_idx = */ 0);

                std::transform(phases.begin(), phases.end(), mob.begin(),
                               [&iq](const int phase) { return iq.mobility(phase).value(); });
            },

            // densityInCell: Reservoir condition phase densities in
            // specified cell.
            [&model = simulator.model()](const int               cell,
                                         const std::vector<int>& phases,
                                         std::vector<Scalar>&    rho)
            {
                const auto& fs = model.intensiveQuantities(cell, /* time_idx = */ 0).fluidState();

                std::transform(phases.begin(), phases.end(), rho.begin(),
                               [&fs](const int phase) { return fs.density(phase).value(); });
            }
        };

        const auto stopped_or_zero_rate_target = this->
            stoppedOrZeroRateTarget(simulator, well_state, deferred_logger);

        this->connections_
            .computeProperties(stopped_or_zero_rate_target, well_state,
                               prop_func, props, deferred_logger);
    }





    template<typename TypeTag>
    void
    StandardWell<TypeTag>::
    computeWellConnectionPressures(const Simulator& simulator,
                                   const WellState<Scalar>& well_state,
                                   DeferredLogger& deferred_logger)
    {
         const auto props = computePropertiesForWellConnectionPressures
             (simulator, well_state);

         computeWellConnectionDensitesPressures(simulator, well_state,
                                                props, deferred_logger);
    }





    template<typename TypeTag>
    void
    StandardWell<TypeTag>::
    solveEqAndUpdateWellState(const Simulator& simulator,
                              WellState<Scalar>& well_state,
                              DeferredLogger& deferred_logger)
    {
        if (!this->isOperableAndSolvable() && !this->wellIsStopped()) return;

        // We assemble the well equations, then we check the convergence,
        // which is why we do not put the assembleWellEq here.
        BVectorWell dx_well(1);
        dx_well[0].resize(this->primary_variables_.numWellEq());
        this->linSys_.solve( dx_well);

        updateWellState(simulator, dx_well, well_state, deferred_logger);
    }





    template<typename TypeTag>
    void
    StandardWell<TypeTag>::
    calculateExplicitQuantities(const Simulator& simulator,
                                const WellState<Scalar>& well_state,
                                DeferredLogger& deferred_logger)
    {
        updatePrimaryVariables(simulator, well_state, deferred_logger);
        computeWellConnectionPressures(simulator, well_state, deferred_logger);
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
    recoverWellSolutionAndUpdateWellState(const Simulator& simulator,
                                          const BVector& x,
                                          WellState<Scalar>& well_state,
                                          DeferredLogger& deferred_logger)
    {
        if (!this->isOperableAndSolvable() && !this->wellIsStopped()) return;

        BVectorWell xw(1);
        xw[0].resize(this->primary_variables_.numWellEq());

        this->linSys_.recoverSolutionWell(x, xw);
        updateWellState(simulator, xw, well_state, deferred_logger);
    }




    template<typename TypeTag>
    void
    StandardWell<TypeTag>::
    computeWellRatesWithBhp(const Simulator& simulator,
                            const Scalar& bhp,
                            std::vector<Scalar>& well_flux,
                            DeferredLogger& deferred_logger) const
    {
        OPM_TIMEFUNCTION();
        const int np = this->number_of_phases_;
        well_flux.resize(np, 0.0);

        const bool allow_cf = this->getAllowCrossFlow();

        for (int perf = 0; perf < this->number_of_local_perforations_; ++perf) {
            const int cell_idx = this->well_cells_[perf];
            const auto& intQuants = simulator.model().intensiveQuantities(cell_idx, /*timeIdx=*/ 0);
            // flux for each perforation
            std::vector<Scalar> mob(this->num_components_, 0.);
            getMobility(simulator, perf, mob, deferred_logger);
            Scalar trans_mult = simulator.problem().template wellTransMultiplier<Scalar>(intQuants, cell_idx);
            const auto& wellstate_nupcol = simulator.problem().wellModel().nupcolWellState().well(this->index_of_well_);
            const std::vector<Scalar> Tw = this->wellIndex(perf, intQuants, trans_mult, wellstate_nupcol);

            std::vector<Scalar> cq_s(this->num_components_, 0.);
            PerforationRates<Scalar> perf_rates;
            computePerfRate(intQuants, mob, bhp, Tw, perf, allow_cf,
                            cq_s, perf_rates, deferred_logger);

            for(int p = 0; p < np; ++p) {
                well_flux[this->modelCompIdxToFlowCompIdx(p)] += cq_s[p];
            }

            // the solvent contribution is added to the gas potentials
            if constexpr (has_solvent) {
                const auto& pu = this->phaseUsage();
                assert(pu.phase_used[Gas]);
                const int gas_pos = pu.phase_pos[Gas];
                well_flux[gas_pos] += cq_s[Indices::contiSolventEqIdx];
            }
        }
        this->parallel_well_info_.communication().sum(well_flux.data(), well_flux.size());
    }



    template<typename TypeTag>
    void
    StandardWell<TypeTag>::
    computeWellRatesWithBhpIterations(const Simulator& simulator,
                                      const Scalar& bhp,
                                      std::vector<Scalar>& well_flux,
                                      DeferredLogger& deferred_logger) const
    {
        // creating a copy of the well itself, to avoid messing up the explicit information
        // during this copy, the only information not copied properly is the well controls
        StandardWell<TypeTag> well_copy(*this);
        well_copy.resetDampening();

        // iterate to get a more accurate well density
        // create a copy of the well_state to use. If the operability checking is sucessful, we use this one
        // to replace the original one
        WellState<Scalar> well_state_copy = simulator.problem().wellModel().wellState();
        const auto& group_state  = simulator.problem().wellModel().groupState();

        // Get the current controls.
        const auto& summary_state = simulator.vanguard().summaryState();
        auto inj_controls = well_copy.well_ecl_.isInjector()
                            ? well_copy.well_ecl_.injectionControls(summary_state)
                            : Well::InjectionControls(0);
        auto prod_controls = well_copy.well_ecl_.isProducer()
                             ? well_copy.well_ecl_.productionControls(summary_state) :
                             Well::ProductionControls(0);

        //  Set current control to bhp, and bhp value in state, modify bhp limit in control object.
        auto& ws = well_state_copy.well(this->index_of_well_);
        if (well_copy.well_ecl_.isInjector()) {
            inj_controls.bhp_limit = bhp;
            ws.injection_cmode = Well::InjectorCMode::BHP;
        } else {
            prod_controls.bhp_limit = bhp;
            ws.production_cmode = Well::ProducerCMode::BHP;
        }
        ws.bhp = bhp;

        // initialized the well rates with the potentials i.e. the well rates based on bhp
        const int np = this->number_of_phases_;
        const Scalar sign = this->well_ecl_.isInjector() ? 1.0 : -1.0;
        for (int phase = 0; phase < np; ++phase){
            well_state_copy.wellRates(this->index_of_well_)[phase]
                    = sign * ws.well_potentials[phase];
        }
        well_copy.updatePrimaryVariables(simulator, well_state_copy, deferred_logger);
        well_copy.computeAccumWell();

        const double dt = simulator.timeStepSize();
        const bool converged = well_copy.iterateWellEqWithControl(simulator, dt, inj_controls, prod_controls, well_state_copy, group_state, deferred_logger);
        if (!converged) {
            const std::string msg = " well " + name() + " did not get converged during well potential calculations "
                                                        " potentials are computed based on unconverged solution";
            deferred_logger.debug(msg);
        }
        well_copy.updatePrimaryVariables(simulator, well_state_copy, deferred_logger);
        well_copy.computeWellConnectionPressures(simulator, well_state_copy, deferred_logger);
        well_copy.computeWellRatesWithBhp(simulator, bhp, well_flux, deferred_logger);
    }




    template<typename TypeTag>
    std::vector<typename StandardWell<TypeTag>::Scalar>
    StandardWell<TypeTag>::
    computeWellPotentialWithTHP(const Simulator& simulator,
                               DeferredLogger& deferred_logger,
                               const WellState<Scalar>& well_state) const
    {
        std::vector<Scalar> potentials(this->number_of_phases_, 0.0);
        const auto& summary_state = simulator.vanguard().summaryState();

        const auto& well = this->well_ecl_;
        if (well.isInjector()){
            const auto& controls = this->well_ecl_.injectionControls(summary_state);
            auto bhp_at_thp_limit = computeBhpAtThpLimitInj(simulator, summary_state, deferred_logger);
            if (bhp_at_thp_limit) {
                const Scalar bhp = std::min(*bhp_at_thp_limit,
                                            static_cast<Scalar>(controls.bhp_limit));
                computeWellRatesWithBhp(simulator, bhp, potentials, deferred_logger);
            } else {
                deferred_logger.warning("FAILURE_GETTING_CONVERGED_POTENTIAL",
                                        "Failed in getting converged thp based potential calculation for well "
                                        + name() + ". Instead the bhp based value is used");
                const Scalar bhp = controls.bhp_limit;
                computeWellRatesWithBhp(simulator, bhp, potentials, deferred_logger);
            }
        } else {
            computeWellRatesWithThpAlqProd(
                simulator, summary_state,
                deferred_logger, potentials, this->getALQ(well_state)
            );
        }

        return potentials;
    }

    template<typename TypeTag>
    bool
    StandardWell<TypeTag>::
    computeWellPotentialsImplicit(const Simulator& simulator,
                                  const WellState<Scalar>& well_state,
                                  std::vector<Scalar>& well_potentials,
                                  DeferredLogger& deferred_logger) const
    {
        // Create a copy of the well.
        // TODO: check if we can avoid taking multiple copies. Call from updateWellPotentials
        // is allready a copy, but not from other calls.
        StandardWell<TypeTag> well_copy(*this);

        // store a copy of the well state, we don't want to update the real well state
        WellState<Scalar> well_state_copy = well_state;
        const auto& group_state = simulator.problem().wellModel().groupState();
        auto& ws = well_state_copy.well(this->index_of_well_);

        // get current controls
        const auto& summary_state = simulator.vanguard().summaryState();
        auto inj_controls = well_copy.well_ecl_.isInjector()
            ? well_copy.well_ecl_.injectionControls(summary_state)
            : Well::InjectionControls(0);
        auto prod_controls = well_copy.well_ecl_.isProducer()
            ? well_copy.well_ecl_.productionControls(summary_state) :
            Well::ProductionControls(0);

        // prepare/modify well state and control
        well_copy.prepareForPotentialCalculations(summary_state, well_state_copy, inj_controls, prod_controls);

       // update connection pressures relative to updated bhp to get better estimate of connection dp
        const int num_perf = ws.perf_data.size();
        for (int perf = 0; perf < num_perf; ++perf) {
            ws.perf_data.pressure[perf] = ws.bhp + well_copy.connections_.pressure_diff(perf);
        }
        // initialize rates from previous potentials
        const int np = this->number_of_phases_;
        bool trivial = true;
        for (int phase = 0; phase < np; ++phase){
            trivial = trivial && (ws.well_potentials[phase] == 0.0) ;
        }
        if (!trivial) {
            const Scalar sign = well_copy.well_ecl_.isInjector() ? 1.0 : -1.0;
            for (int phase = 0; phase < np; ++phase) {
                ws.surface_rates[phase] = sign * ws.well_potentials[phase];
            }
        }

        well_copy.calculateExplicitQuantities(simulator, well_state_copy, deferred_logger);
        const double dt = simulator.timeStepSize();
        // iterate to get a solution at the given bhp.
        bool converged = false;
        if (this->well_ecl_.isProducer()) {
            converged = well_copy.solveWellWithOperabilityCheck(simulator, dt, inj_controls, prod_controls, well_state_copy, group_state, deferred_logger);
        } else {
            converged = well_copy.iterateWellEqWithSwitching(simulator, dt, inj_controls, prod_controls, well_state_copy, group_state, deferred_logger);
        }

        // fetch potentials (sign is updated on the outside).
        well_potentials.clear();
        well_potentials.resize(np, 0.0);
        for (int comp_idx = 0; comp_idx < this->num_components_; ++comp_idx) {
            if (has_solvent && comp_idx == Indices::contiSolventEqIdx) continue; // we do not store the solvent in the well_potentials
            const EvalWell rate = well_copy.primary_variables_.getQs(comp_idx);
            well_potentials[this->modelCompIdxToFlowCompIdx(comp_idx)] = rate.value();
        }

        // the solvent contribution is added to the gas potentials
        if constexpr (has_solvent) {
            const auto& pu = this->phaseUsage();
            assert(pu.phase_used[Gas]);
            const int gas_pos = pu.phase_pos[Gas];
            const EvalWell rate = well_copy.primary_variables_.getQs(Indices::contiSolventEqIdx);
            well_potentials[gas_pos] += rate.value();
        }
        return converged;
    }


    template<typename TypeTag>
    typename StandardWell<TypeTag>::Scalar
    StandardWell<TypeTag>::
    computeWellRatesAndBhpWithThpAlqProd(const Simulator &simulator,
                               const SummaryState &summary_state,
                               DeferredLogger& deferred_logger,
                               std::vector<Scalar>& potentials,
                               Scalar alq) const
    {
        Scalar bhp;
        auto bhp_at_thp_limit = computeBhpAtThpLimitProdWithAlq(
                              simulator, summary_state, alq, deferred_logger, /*iterate_if_no_solution */ true);
        if (bhp_at_thp_limit) {
            const auto& controls = this->well_ecl_.productionControls(summary_state);
            bhp = std::max(*bhp_at_thp_limit,
                           static_cast<Scalar>(controls.bhp_limit));
            computeWellRatesWithBhp(simulator, bhp, potentials, deferred_logger);
        }
        else {
            deferred_logger.warning("FAILURE_GETTING_CONVERGED_POTENTIAL",
                "Failed in getting converged thp based potential calculation for well "
                + name() + ". Instead the bhp based value is used");
            const auto& controls = this->well_ecl_.productionControls(summary_state);
            bhp = controls.bhp_limit;
            computeWellRatesWithBhp(simulator, bhp, potentials, deferred_logger);
        }
        return bhp;
    }

    template<typename TypeTag>
    void
    StandardWell<TypeTag>::
    computeWellRatesWithThpAlqProd(const Simulator& simulator,
                                   const SummaryState& summary_state,
                                   DeferredLogger& deferred_logger,
                                   std::vector<Scalar>& potentials,
                                   Scalar alq) const
    {
        /*double bhp =*/
        computeWellRatesAndBhpWithThpAlqProd(simulator,
                                             summary_state,
                                             deferred_logger,
                                             potentials,
                                             alq);
    }

    template<typename TypeTag>
    void
    StandardWell<TypeTag>::
    computeWellPotentials(const Simulator& simulator,
                          const WellState<Scalar>& well_state,
                          std::vector<Scalar>& well_potentials,
                          DeferredLogger& deferred_logger) // const
    {
        const auto [compute_potential, bhp_controlled_well] =
            this->WellInterfaceGeneric<Scalar>::computeWellPotentials(well_potentials, well_state);

        if (!compute_potential) {
            return;
        }

        bool converged_implicit = false;
        // for newly opened wells we dont compute the potentials implicit
        // group controlled wells with defaulted guiderates will have zero targets as
        // the potentials are used to compute the well fractions.
        if (this->param_.local_well_solver_control_switching_ && !(this->changed_to_open_this_step_ && this->wellUnderZeroRateTarget(simulator, well_state, deferred_logger))) {
            converged_implicit = computeWellPotentialsImplicit(simulator, well_state, well_potentials, deferred_logger);
        }
        if (!converged_implicit) {
            // does the well have a THP related constraint?
            const auto& summaryState = simulator.vanguard().summaryState();
            if (!Base::wellHasTHPConstraints(summaryState) || bhp_controlled_well) {
                // get the bhp value based on the bhp constraints
                Scalar bhp = WellBhpThpCalculator(*this).mostStrictBhpFromBhpLimits(summaryState);

                // In some very special cases the bhp pressure target are
                // temporary violated. This may lead to too small or negative potentials
                // that could lead to premature shutting of wells.
                // As a remedy the bhp that gives the largest potential is used.
                // For converged cases, ws.bhp <=bhp for injectors and ws.bhp >= bhp,
                // and the potentials will be computed using the limit as expected.
                const auto& ws = well_state.well(this->index_of_well_);
                if (this->isInjector())
                    bhp = std::max(ws.bhp, bhp);
                else
                    bhp = std::min(ws.bhp, bhp);

                assert(std::abs(bhp) != std::numeric_limits<Scalar>::max());
                computeWellRatesWithBhpIterations(simulator, bhp, well_potentials, deferred_logger);
            } else {
                // the well has a THP related constraint
                well_potentials = computeWellPotentialWithTHP(simulator, deferred_logger, well_state);
            }
        }

        this->checkNegativeWellPotentials(well_potentials,
                                          this->param_.check_well_operability_,
                                          deferred_logger);
    }







    template<typename TypeTag>
    typename StandardWell<TypeTag>::Scalar
    StandardWell<TypeTag>::
    connectionDensity([[maybe_unused]] const int globalConnIdx,
                      const int openConnIdx) const
    {
        return (openConnIdx < 0)
            ? 0.0
            : this->connections_.rho(openConnIdx);
    }





    template<typename TypeTag>
    void
    StandardWell<TypeTag>::
    updatePrimaryVariables(const Simulator& simulator,
                           const WellState<Scalar>& well_state,
                           DeferredLogger& deferred_logger)
    {
        if (!this->isOperableAndSolvable() && !this->wellIsStopped()) return;

        const bool stop_or_zero_rate_target = this->stoppedOrZeroRateTarget(simulator, well_state, deferred_logger);
        this->primary_variables_.update(well_state, stop_or_zero_rate_target, deferred_logger);

        // other primary variables related to polymer injection
        if constexpr (Base::has_polymermw) {
            this->primary_variables_.updatePolyMW(well_state);
        }

        this->primary_variables_.checkFinite(deferred_logger);
    }




    template<typename TypeTag>
    typename StandardWell<TypeTag>::Scalar
    StandardWell<TypeTag>::
    getRefDensity() const
    {
        return this->connections_.rho();
    }




    template<typename TypeTag>
    void
    StandardWell<TypeTag>::
    updateWaterMobilityWithPolymer(const Simulator& simulator,
                                   const int perf,
                                   std::vector<EvalWell>& mob,
                                   DeferredLogger& deferred_logger) const
    {
        const int cell_idx = this->well_cells_[perf];
        const auto& int_quant = simulator.model().intensiveQuantities(cell_idx, /*timeIdx=*/ 0);
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
            const bool allow_cf = this->getAllowCrossFlow() || openCrossFlowAvoidSingularity(simulator);
            const EvalWell& bhp = this->primary_variables_.eval(Bhp);

            std::vector<EvalWell> cq_s(this->num_components_, {this->primary_variables_.numWellEq() + Indices::numEq, 0.});
            PerforationRates<Scalar> perf_rates;
            Scalar trans_mult = simulator.problem().template wellTransMultiplier<Scalar>(int_quant, cell_idx);
            const auto& wellstate_nupcol = simulator.problem().wellModel().nupcolWellState().well(this->index_of_well_);
            const std::vector<Scalar> Tw = this->wellIndex(perf, int_quant, trans_mult, wellstate_nupcol);
            computePerfRate(int_quant, mob, bhp, Tw, perf, allow_cf, cq_s,
                            perf_rates, deferred_logger);
            // TODO: make area a member
            const Scalar area = 2 * M_PI * this->perf_rep_radius_[perf] * this->perf_length_[perf];
            const auto& material_law_manager = simulator.problem().materialLawManager();
            const auto& scaled_drainage_info =
                        material_law_manager->oilWaterScaledEpsInfoDrainage(cell_idx);
            const Scalar swcr = scaled_drainage_info.Swcr;
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
                                                    const WellState<Scalar>& well_state) const
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
    pskinwater(const Scalar throughput,
               const EvalWell& water_velocity,
              DeferredLogger& deferred_logger) const
    {
        if constexpr (Base::has_polymermw) {
            const int water_table_id = this->polymerWaterTable_();
            if (water_table_id <= 0) {
                OPM_DEFLOG_THROW(std::runtime_error,
                                 fmt::format("Unused SKPRWAT table id used for well {}", name()),
                                 deferred_logger);
            }
            const auto& water_table_func = PolymerModule::getSkprwatTable(water_table_id);
            const EvalWell throughput_eval(this->primary_variables_.numWellEq() + Indices::numEq, throughput);
            // the skin pressure when injecting water, which also means the polymer concentration is zero
            EvalWell pskin_water(this->primary_variables_.numWellEq() + Indices::numEq, 0.0);
            pskin_water = water_table_func.eval(throughput_eval, water_velocity);
            return pskin_water;
        } else {
            OPM_DEFLOG_THROW(std::runtime_error,
                             fmt::format("Polymermw is not activated, while injecting "
                                         "skin pressure is requested for well {}", name()),
                             deferred_logger);
        }
    }





    template<typename TypeTag>
    typename StandardWell<TypeTag>::EvalWell
    StandardWell<TypeTag>::
    pskin(const Scalar throughput,
              const EvalWell& water_velocity,
              const EvalWell& poly_inj_conc,
              DeferredLogger& deferred_logger) const
    {
        if constexpr (Base::has_polymermw) {
            const Scalar sign = water_velocity >= 0. ? 1.0 : -1.0;
            const EvalWell water_velocity_abs = abs(water_velocity);
            if (poly_inj_conc == 0.) {
                return sign * pskinwater(throughput, water_velocity_abs, deferred_logger);
            }
            const int polymer_table_id = this->polymerTable_();
            if (polymer_table_id <= 0) {
                OPM_DEFLOG_THROW(std::runtime_error,
                                 fmt::format("Unavailable SKPRPOLY table id used for well {}", name()),
                                 deferred_logger);
            }
            const auto& skprpolytable = PolymerModule::getSkprpolyTable(polymer_table_id);
            const Scalar reference_concentration = skprpolytable.refConcentration;
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
            OPM_DEFLOG_THROW(std::runtime_error,
                             fmt::format("Polymermw is not activated, while injecting "
                                         "skin pressure is requested for well {}", name()),
                             deferred_logger);
        }
    }





    template<typename TypeTag>
    typename StandardWell<TypeTag>::EvalWell
    StandardWell<TypeTag>::
    wpolymermw(const Scalar throughput,
               const EvalWell& water_velocity,
               DeferredLogger& deferred_logger) const
    {
        if constexpr (Base::has_polymermw) {
            const int table_id = this->polymerInjTable_();
            const auto& table_func = PolymerModule::getPlymwinjTable(table_id);
            const EvalWell throughput_eval(this->primary_variables_.numWellEq() + Indices::numEq, throughput);
            EvalWell molecular_weight(this->primary_variables_.numWellEq() + Indices::numEq, 0.);
            if (this->wpolymer() == 0.) { // not injecting polymer
                return molecular_weight;
            }
            molecular_weight = table_func.eval(throughput_eval, abs(water_velocity));
            return molecular_weight;
        } else {
            OPM_DEFLOG_THROW(std::runtime_error,
                             fmt::format("Polymermw is not activated, while injecting "
                                         "polymer molecular weight is requested for well {}", name()),
                             deferred_logger);
        }
    }





    template<typename TypeTag>
    void
    StandardWell<TypeTag>::
    updateWaterThroughput([[maybe_unused]] const double dt,
                          WellState<Scalar>&            well_state) const
    {
        if constexpr (Base::has_polymermw) {
            if (!this->isInjector()) {
                return;
            }

            auto& perf_water_throughput = well_state.well(this->index_of_well_)
                .perf_data.water_throughput;

            for (int perf = 0; perf < this->number_of_local_perforations_; ++perf) {
                const Scalar perf_water_vel =
                    this->primary_variables_.value(Bhp + 1 + perf);

                // we do not consider the formation damage due to water
                // flowing from reservoir into wellbore
                if (perf_water_vel > Scalar{0}) {
                    perf_water_throughput[perf] += perf_water_vel * dt;
                }
            }
        }
    }





    template<typename TypeTag>
    void
    StandardWell<TypeTag>::
    handleInjectivityRate(const Simulator& simulator,
                          const int perf,
                          std::vector<EvalWell>& cq_s) const
    {
        const int cell_idx = this->well_cells_[perf];
        const auto& int_quants = simulator.model().intensiveQuantities(cell_idx, /*timeIdx=*/ 0);
        const auto& fs = int_quants.fluidState();
        const EvalWell b_w = this->extendEval(fs.invB(FluidSystem::waterPhaseIdx));
        const Scalar area = M_PI * this->bore_diameters_[perf] * this->perf_length_[perf];
        const int wat_vel_index = Bhp + 1 + perf;
        const unsigned water_comp_idx = Indices::canonicalToActiveComponentIndex(FluidSystem::waterCompIdx);

        // water rate is update to use the form from water velocity, since water velocity is
        // a primary variable now
        cq_s[water_comp_idx] = area * this->primary_variables_.eval(wat_vel_index) * b_w;
    }




    template<typename TypeTag>
    void
    StandardWell<TypeTag>::
    handleInjectivityEquations(const Simulator& simulator,
                               const WellState<Scalar>& well_state,
                               const int perf,
                               const EvalWell& water_flux_s,
                               DeferredLogger& deferred_logger)
    {
        const int cell_idx = this->well_cells_[perf];
        const auto& int_quants = simulator.model().intensiveQuantities(cell_idx, /*timeIdx=*/ 0);
        const auto& fs = int_quants.fluidState();
        const EvalWell b_w = this->extendEval(fs.invB(FluidSystem::waterPhaseIdx));
        const EvalWell water_flux_r = water_flux_s / b_w;
        const Scalar area = M_PI * this->bore_diameters_[perf] * this->perf_length_[perf];
        const EvalWell water_velocity = water_flux_r / area;
        const int wat_vel_index = Bhp + 1 + perf;

        // equation for the water velocity
        const EvalWell eq_wat_vel = this->primary_variables_.eval(wat_vel_index) - water_velocity;

        const auto& ws = well_state.well(this->index_of_well_);
        const auto& perf_data = ws.perf_data;
        const auto& perf_water_throughput = perf_data.water_throughput;
        const Scalar throughput = perf_water_throughput[perf];
        const int pskin_index = Bhp + 1 + this->number_of_local_perforations_ + perf;

        EvalWell poly_conc(this->primary_variables_.numWellEq() + Indices::numEq, 0.0);
        poly_conc.setValue(this->wpolymer());

        // equation for the skin pressure
        const EvalWell eq_pskin = this->primary_variables_.eval(pskin_index)
                                  - pskin(throughput, this->primary_variables_.eval(wat_vel_index), poly_conc, deferred_logger);

        StandardWellAssemble<FluidSystem,Indices>(*this).
                assembleInjectivityEq(eq_pskin,
                                      eq_wat_vel,
                                      pskin_index,
                                      wat_vel_index,
                                      perf,
                                      this->primary_variables_.numWellEq(),
                                      this->linSys_);
    }





    template<typename TypeTag>
    void
    StandardWell<TypeTag>::
    checkConvergenceExtraEqs(const std::vector<Scalar>& res,
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
                               const WellState<Scalar>& well_state,
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
                const Scalar throughput = perf_water_throughput[perf];
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
    std::optional<typename StandardWell<TypeTag>::Scalar>
    StandardWell<TypeTag>::
    computeBhpAtThpLimitProd(const WellState<Scalar>& well_state,
                             const Simulator& simulator,
                             const SummaryState& summary_state,
                             DeferredLogger& deferred_logger) const
    {
        return computeBhpAtThpLimitProdWithAlq(simulator,
                                               summary_state,
                                               this->getALQ(well_state),
                                               deferred_logger,
                                               /*iterate_if_no_solution */ true);
    }

    template<typename TypeTag>
    std::optional<typename StandardWell<TypeTag>::Scalar>
    StandardWell<TypeTag>::
    computeBhpAtThpLimitProdWithAlq(const Simulator& simulator,
                                    const SummaryState& summary_state,
                                    const Scalar alq_value,
                                    DeferredLogger& deferred_logger,
                                    bool iterate_if_no_solution) const
    {
        OPM_TIMEFUNCTION();
        // Make the frates() function.
        auto frates = [this, &simulator, &deferred_logger](const Scalar bhp) {
            // Not solving the well equations here, which means we are
            // calculating at the current Fg/Fw values of the
            // well. This does not matter unless the well is
            // crossflowing, and then it is likely still a good
            // approximation.
            std::vector<Scalar> rates(3);
            computeWellRatesWithBhp(simulator, bhp, rates, deferred_logger);
            this->adaptRatesForVFP(rates);
            return rates;
        };

        Scalar max_pressure = 0.0;
        for (int perf = 0; perf < this->number_of_local_perforations_; ++perf) {
            const int cell_idx = this->well_cells_[perf];
            const auto& int_quants = simulator.model().intensiveQuantities(cell_idx, /*timeIdx=*/ 0);
            const auto& fs = int_quants.fluidState();
            Scalar pressure_cell = this->getPerfCellPressure(fs).value();
            max_pressure = std::max(max_pressure, pressure_cell);
        }
        auto bhpAtLimit = WellBhpThpCalculator(*this).computeBhpAtThpLimitProd(frates,
                                                                               summary_state,
                                                                               max_pressure,
                                                                               this->connections_.rho(),
                                                                               alq_value,
                                                                               this->getTHPConstraint(summary_state),
                                                                               deferred_logger);

        if (bhpAtLimit) {
            auto v = frates(*bhpAtLimit);
            if (std::all_of(v.cbegin(), v.cend(), [](Scalar i){ return i <= 0; }) ) {
                return bhpAtLimit;
            }
        }

        if (!iterate_if_no_solution)
            return std::nullopt;

        auto fratesIter = [this, &simulator, &deferred_logger](const Scalar bhp) {
            // Solver the well iterations to see if we are
            // able to get a solution with an update
            // solution
            std::vector<Scalar> rates(3);
            computeWellRatesWithBhpIterations(simulator, bhp, rates, deferred_logger);
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


        if (bhpAtLimit) {
            // should we use fratesIter here since fratesIter is used in computeBhpAtThpLimitProd above?
            auto v = frates(*bhpAtLimit);
            if (std::all_of(v.cbegin(), v.cend(), [](Scalar i){ return i <= 0; }) ) {
                return bhpAtLimit;
            }
        }

        // we still don't get a valied solution.
        return std::nullopt;
    }



    template<typename TypeTag>
    std::optional<typename StandardWell<TypeTag>::Scalar>
    StandardWell<TypeTag>::
    computeBhpAtThpLimitInj(const Simulator& simulator,
                            const SummaryState& summary_state,
                            DeferredLogger& deferred_logger) const
    {
        // Make the frates() function.
        auto frates = [this, &simulator, &deferred_logger](const Scalar bhp) {
            // Not solving the well equations here, which means we are
            // calculating at the current Fg/Fw values of the
            // well. This does not matter unless the well is
            // crossflowing, and then it is likely still a good
            // approximation.
            std::vector<Scalar> rates(3);
            computeWellRatesWithBhp(simulator, bhp, rates, deferred_logger);
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
    iterateWellEqWithControl(const Simulator& simulator,
                             const double dt,
                             const Well::InjectionControls& inj_controls,
                             const Well::ProductionControls& prod_controls,
                             WellState<Scalar>& well_state,
                             const GroupState<Scalar>& group_state,
                             DeferredLogger& deferred_logger)
    {
        updatePrimaryVariables(simulator, well_state, deferred_logger);

        const int max_iter = this->param_.max_inner_iter_wells_;
        int it = 0;
        bool converged;
        bool relax_convergence = false;
        this->regularize_ = false;
        do {
            assembleWellEqWithoutIteration(simulator, dt, inj_controls, prod_controls, well_state, group_state, deferred_logger);

            if (it > this->param_.strict_inner_iter_wells_) {
                relax_convergence = true;
                this->regularize_ = true;
            }

            auto report = getWellConvergence(simulator, well_state, Base::B_avg_, deferred_logger, relax_convergence);

            converged = report.converged();
            if (converged) {
                break;
            }

            ++it;
            solveEqAndUpdateWellState(simulator, well_state, deferred_logger);

            // TODO: when this function is used for well testing purposes, will need to check the controls, so that we will obtain convergence
            // under the most restrictive control. Based on this converged results, we can check whether to re-open the well. Either we refactor
            // this function or we use different functions for the well testing purposes.
            // We don't allow for switching well controls while computing well potentials and testing wells
            // updateWellControl(simulator, well_state, deferred_logger);
        } while (it < max_iter);

        return converged;
    }


    template<typename TypeTag>
    bool
    StandardWell<TypeTag>::
    iterateWellEqWithSwitching(const Simulator& simulator,
                               const double dt,
                               const Well::InjectionControls& inj_controls,
                               const Well::ProductionControls& prod_controls,
                               WellState<Scalar>& well_state,
                               const GroupState<Scalar>& group_state,
                               DeferredLogger& deferred_logger,
                               const bool fixed_control /*false*/,
                               const bool fixed_status /*false*/)
    {
        updatePrimaryVariables(simulator, well_state, deferred_logger);

        const int max_iter = this->param_.max_inner_iter_wells_;
        int it = 0;
        bool converged = false;
        bool relax_convergence = false;
        this->regularize_ = false;
        const auto& summary_state = simulator.vanguard().summaryState();

        // Always take a few (more than one) iterations after a switch before allowing a new switch
        // The optimal number here is subject to further investigation, but it has been observerved
        // that unless this number is >1, we may get stuck in a cycle
        constexpr int min_its_after_switch = 4;
        int its_since_last_switch = min_its_after_switch;
        int switch_count= 0;
        // if we fail to solve eqs, we reset status/operability before leaving
        const auto well_status_orig = this->wellStatus_;
        const auto operability_orig = this->operability_status_;
        auto well_status_cur = well_status_orig;
        int status_switch_count = 0;
        // don't allow opening wells that are stopped from schedule or has a stopped well state
        const bool allow_open =  this->well_ecl_.getStatus() == WellStatus::OPEN &&
                                 well_state.well(this->index_of_well_).status == WellStatus::OPEN;
        // don't allow switcing for wells under zero rate target or requested fixed status and control
        const bool allow_switching =
            !this->wellUnderZeroRateTarget(simulator, well_state, deferred_logger) &&
            (!fixed_control || !fixed_status) && allow_open;

        bool changed = false;
        bool final_check = false;
        // well needs to be set operable or else solving/updating of re-opened wells is skipped
        this->operability_status_.resetOperability();
        this->operability_status_.solvable = true;
        do {
            its_since_last_switch++;
            if (allow_switching && its_since_last_switch >= min_its_after_switch){
                const Scalar wqTotal = this->primary_variables_.eval(WQTotal).value();
                changed = this->updateWellControlAndStatusLocalIteration(simulator, well_state, group_state,
                                                                         inj_controls, prod_controls, wqTotal,
                                                                         deferred_logger, fixed_control, fixed_status);
                if (changed){
                    its_since_last_switch = 0;
                    switch_count++;
                    if (well_status_cur != this->wellStatus_) {
                        well_status_cur = this->wellStatus_;
                        status_switch_count++;
                    }
                }
                if (!changed && final_check) {
                    break;
                } else {
                    final_check = false;
                }
            }

            assembleWellEqWithoutIteration(simulator, dt, inj_controls, prod_controls, well_state, group_state, deferred_logger);

            if (it > this->param_.strict_inner_iter_wells_) {
                relax_convergence = true;
                this->regularize_ = true;
            }

            auto report = getWellConvergence(simulator, well_state, Base::B_avg_, deferred_logger, relax_convergence);

            converged = report.converged();
            if (converged) {
                // if equations are sufficiently linear they might converge in less than min_its_after_switch
                // in this case, make sure all constraints are satisfied before returning
                if (switch_count > 0 && its_since_last_switch < min_its_after_switch) {
                    final_check = true;
                    its_since_last_switch = min_its_after_switch;
                } else {
                    break;
                }
            }

            ++it;
            solveEqAndUpdateWellState(simulator, well_state, deferred_logger);

        } while (it < max_iter);

        if (converged) {
            if (allow_switching){
                // update operability if status change
                const bool is_stopped = this->wellIsStopped();
                if (this->wellHasTHPConstraints(summary_state)){
                    this->operability_status_.can_obtain_bhp_with_thp_limit = !is_stopped;
                    this->operability_status_.obey_thp_limit_under_bhp_limit = !is_stopped;
                } else {
                    this->operability_status_.operable_under_only_bhp_limit = !is_stopped;
                }
            }
        } else {
            this->wellStatus_ = well_status_orig;
            this->operability_status_ = operability_orig;
            const std::string message = fmt::format("   Well {} did not converge in {} inner iterations ("
                                                    "{} switches, {} status changes).", this->name(), it, switch_count, status_switch_count);
            deferred_logger.debug(message);
            // add operability here as well ?
        }
        return converged;
    }

    template<typename TypeTag>
    std::vector<typename StandardWell<TypeTag>::Scalar>
    StandardWell<TypeTag>::
    computeCurrentWellRates(const Simulator& simulator,
                            DeferredLogger& deferred_logger) const
    {
        // Calculate the rates that follow from the current primary variables.
        std::vector<Scalar> well_q_s(this->num_components_, 0.);
        const EvalWell& bhp = this->primary_variables_.eval(Bhp);
        const bool allow_cf = this->getAllowCrossFlow() || openCrossFlowAvoidSingularity(simulator);
        for (int perf = 0; perf < this->number_of_local_perforations_; ++perf) {
            const int cell_idx = this->well_cells_[perf];
            const auto& intQuants = simulator.model().intensiveQuantities(cell_idx, /*timeIdx=*/ 0);
            std::vector<Scalar> mob(this->num_components_, 0.);
            getMobility(simulator, perf, mob, deferred_logger);
            std::vector<Scalar> cq_s(this->num_components_, 0.);
            Scalar trans_mult = simulator.problem().template wellTransMultiplier<Scalar>(intQuants,  cell_idx);
            const auto& wellstate_nupcol = simulator.problem().wellModel().nupcolWellState().well(this->index_of_well_);
            const std::vector<Scalar> Tw = this->wellIndex(perf, intQuants, trans_mult, wellstate_nupcol);
            PerforationRates<Scalar> perf_rates;
            computePerfRate(intQuants, mob, bhp.value(), Tw, perf, allow_cf,
                            cq_s, perf_rates, deferred_logger);
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
    std::vector<typename StandardWell<TypeTag>::Scalar>
    StandardWell<TypeTag>::
    getPrimaryVars() const
    {
        const int num_pri_vars = this->primary_variables_.numWellEq();
        std::vector<Scalar> retval(num_pri_vars);
        for (int ii = 0; ii < num_pri_vars; ++ii) {
            retval[ii] = this->primary_variables_.value(ii);
        }
        return retval;
    }





    template <typename TypeTag>
    int
    StandardWell<TypeTag>::
    setPrimaryVars(typename std::vector<Scalar>::const_iterator it)
    {
        const int num_pri_vars = this->primary_variables_.numWellEq();
        for (int ii = 0; ii < num_pri_vars; ++ii) {
            this->primary_variables_.setValue(ii, it[ii]);
        }
        return num_pri_vars;
    }


    template <typename TypeTag>
    typename StandardWell<TypeTag>::Eval
    StandardWell<TypeTag>::
    connectionRateEnergy(const Scalar maxOilSaturation,
                         const std::vector<EvalWell>& cq_s,
                         const IntensiveQuantities& intQuants,
                         DeferredLogger& deferred_logger) const
    {
        auto fs = intQuants.fluidState();
        Eval result = 0;
        for (unsigned phaseIdx = 0; phaseIdx < FluidSystem::numPhases; ++phaseIdx) {
            if (!FluidSystem::phaseIsActive(phaseIdx)) {
                continue;
            }

            // convert to reservoir conditions
            EvalWell cq_r_thermal(this->primary_variables_.numWellEq() + Indices::numEq, 0.);
            const unsigned activeCompIdx = Indices::canonicalToActiveComponentIndex(FluidSystem::solventComponentIndex(phaseIdx));
            const bool both_oil_gas = FluidSystem::phaseIsActive(FluidSystem::oilPhaseIdx) && FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx);
            if (!both_oil_gas || FluidSystem::waterPhaseIdx == phaseIdx) {
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
                    deferred_logger.debug(
                        fmt::format("Problematic d value {} obtained for well {}"
                                    " during calculateSinglePerf with rs {}"
                                    ", rv {}. Continue as if no dissolution (rs = 0) and"
                                    " vaporization (rv = 0) for this connection.",
                                    d, this->name(), fs.Rs(), fs.Rv()));
                    cq_r_thermal = cq_s[activeCompIdx] / this->extendEval(fs.invB(phaseIdx));
                } else {
                    if (FluidSystem::gasPhaseIdx == phaseIdx) {
                        cq_r_thermal = (cq_s[gasCompIdx] -
                                        this->extendEval(fs.Rs()) * cq_s[oilCompIdx]) /
                                        (d * this->extendEval(fs.invB(phaseIdx)) );
                    } else if (FluidSystem::oilPhaseIdx == phaseIdx) {
                        // q_or = 1 / (b_o * d) * (q_os - rv * q_gs)
                        cq_r_thermal = (cq_s[oilCompIdx] - this->extendEval(fs.Rv()) *
                                        cq_s[gasCompIdx]) /
                                       (d * this->extendEval(fs.invB(phaseIdx)) );
                    }
                }
            }

            // change temperature for injecting fluids
            if (this->isInjector() && !this->wellIsStopped() && cq_r_thermal > 0.0){
                // only handles single phase injection now
                assert(this->well_ecl_.injectorType() != InjectorType::MULTI);
                fs.setTemperature(this->well_ecl_.inj_temperature());
                typedef typename std::decay<decltype(fs)>::type::Scalar FsScalar;
                typename FluidSystem::template ParameterCache<FsScalar> paramCache;
                const unsigned pvtRegionIdx = intQuants.pvtRegionIndex();
                paramCache.setRegionIndex(pvtRegionIdx);
                paramCache.setMaxOilSat(maxOilSaturation);
                paramCache.updatePhase(fs, phaseIdx);

                const auto& rho = FluidSystem::density(fs, paramCache, phaseIdx);
                fs.setDensity(phaseIdx, rho);
                const auto& h = FluidSystem::enthalpy(fs, paramCache, phaseIdx);
                fs.setEnthalpy(phaseIdx, h);
                cq_r_thermal *= this->extendEval(fs.enthalpy(phaseIdx)) * this->extendEval(fs.density(phaseIdx));
                result += getValue(cq_r_thermal);
            } else if (cq_r_thermal > 0.0) {
                cq_r_thermal *= getValue(fs.enthalpy(phaseIdx)) * getValue(fs.density(phaseIdx));
                result += Base::restrictEval(cq_r_thermal);
            } else {
                // compute the thermal flux
                cq_r_thermal *= this->extendEval(fs.enthalpy(phaseIdx)) * this->extendEval(fs.density(phaseIdx));
                result += Base::restrictEval(cq_r_thermal);
            }
        }

        return result * this->well_efficiency_factor_;
    }
} // namespace Opm

#endif
