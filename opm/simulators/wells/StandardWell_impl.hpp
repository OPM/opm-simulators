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

#include <opm/simulators/utils/DeferredLoggingErrorHelpers.hpp>

namespace Opm
{

    template<typename TypeTag>
    StandardWell<TypeTag>::
    StandardWell(const Well2& well, const int time_step, const Wells* wells,
                 const ModelParameters& param,
                 const RateConverterType& rate_converter,
                 const int pvtRegionIdx,
                 const int num_components)
    : Base(well, time_step, wells, param, rate_converter, pvtRegionIdx, num_components)
    , perf_densities_(number_of_perforations_)
    , perf_pressure_diffs_(number_of_perforations_)
    , primary_variables_(numWellEq, 0.0)
    , primary_variables_evaluation_(numWellEq) // the number of the primary variables
    , F0_(numWellConservationEq)
    , ipr_a_(number_of_phases_)
    , ipr_b_(number_of_phases_)
    {
        assert(num_components_ == numWellConservationEq);

        duneB_.setBuildMode( OffDiagMatWell::row_wise );
        duneC_.setBuildMode( OffDiagMatWell::row_wise );
        invDuneD_.setBuildMode( DiagMatWell::row_wise );
    }





    template<typename TypeTag>
    void
    StandardWell<TypeTag>::
    init(const PhaseUsage* phase_usage_arg,
         const std::vector<double>& depth_arg,
         const double gravity_arg,
         const int num_cells)
    {
        Base::init(phase_usage_arg, depth_arg, gravity_arg, num_cells);

        perf_depth_.resize(number_of_perforations_, 0.);
        for (int perf = 0; perf < number_of_perforations_; ++perf) {
            const int cell_idx = well_cells_[perf];
            perf_depth_[perf] = depth_arg[cell_idx];
        }

        // setup sparsity pattern for the matrices
        //[A C^T    [x    =  [ res
        // B D] x_well]      res_well]
        // set the size of the matrices
        invDuneD_.setSize(1, 1, 1);
        duneB_.setSize(1, num_cells, number_of_perforations_);
        duneC_.setSize(1, num_cells, number_of_perforations_);

        for (auto row=invDuneD_.createbegin(), end = invDuneD_.createend(); row!=end; ++row) {
            // Add nonzeros for diagonal
            row.insert(row.index());
        }

        for (auto row = duneB_.createbegin(), end = duneB_.createend(); row!=end; ++row) {
            for (int perf = 0 ; perf < number_of_perforations_; ++perf) {
                const int cell_idx = well_cells_[perf];
                row.insert(cell_idx);
            }
        }

        // make the C^T matrix
        for (auto row = duneC_.createbegin(), end = duneC_.createend(); row!=end; ++row) {
            for (int perf = 0; perf < number_of_perforations_; ++perf) {
                const int cell_idx = well_cells_[perf];
                row.insert(cell_idx);
            }
        }

        resWell_.resize(1);

        // resize temporary class variables
        Bx_.resize( duneB_.N() );
        invDrw_.resize( invDuneD_.N() );
    }





    template<typename TypeTag>
    void StandardWell<TypeTag>::
    initPrimaryVariablesEvaluation() const
    {
        for (int eqIdx = 0; eqIdx < numWellEq; ++eqIdx) {
            assert( (size_t)eqIdx < primary_variables_.size() );

            primary_variables_evaluation_[eqIdx] = 0.0;
            primary_variables_evaluation_[eqIdx].setValue(primary_variables_[eqIdx]);
            primary_variables_evaluation_[eqIdx].setDerivative(numEq + eqIdx, 1.0);
        }
    }





    template<typename TypeTag>
    const typename StandardWell<TypeTag>::EvalWell&
    StandardWell<TypeTag>::
    getBhp() const
    {
        return primary_variables_evaluation_[Bhp];
    }





    template<typename TypeTag>
    const typename StandardWell<TypeTag>::EvalWell&
    StandardWell<TypeTag>::
    getWQTotal() const
    {
        return primary_variables_evaluation_[WQTotal];
    }





    template<typename TypeTag>
    typename StandardWell<TypeTag>::EvalWell
    StandardWell<TypeTag>::
    getQs(const int comp_idx) const
    {
        // Note: currently, the WQTotal definition is still depends on Injector/Producer.
        assert(comp_idx < num_components_);

        if (well_type_ == INJECTOR) { // only single phase injection
            // TODO: using comp_frac here is dangerous, it should be changed later
            // Most likely, it should be changed to use distr, or at least, we need to update comp_frac_ based on distr
            // while solvent might complicate the situation
            const auto pu = phaseUsage();
            const int legacyCompIdx = ebosCompIdxToFlowCompIdx(comp_idx);
            double comp_frac = 0.0;
            if (has_solvent && comp_idx == contiSolventEqIdx) { // solvent
                comp_frac = comp_frac_[pu.phase_pos[ Gas ]] * wsolvent();
            } else if (legacyCompIdx == pu.phase_pos[ Gas ]) {
                comp_frac = comp_frac_[legacyCompIdx];
                if (has_solvent) {
                    comp_frac *= (1.0 - wsolvent());
                }
            } else {
                comp_frac = comp_frac_[legacyCompIdx];
            }

            return comp_frac * primary_variables_evaluation_[WQTotal];
        } else { // producers
            return primary_variables_evaluation_[WQTotal] * wellVolumeFractionScaled(comp_idx);
        }
    }






    template<typename TypeTag>
    typename StandardWell<TypeTag>::EvalWell
    StandardWell<TypeTag>::
    wellVolumeFractionScaled(const int compIdx) const
    {

        const int legacyCompIdx = ebosCompIdxToFlowCompIdx(compIdx);
        const double scal = scalingFactor(legacyCompIdx);
        if (scal > 0)
            return  wellVolumeFraction(compIdx) / scal;

        // the scaling factor may be zero for RESV controlled wells.
        return wellVolumeFraction(compIdx);
    }





    template<typename TypeTag>
    typename StandardWell<TypeTag>::EvalWell
    StandardWell<TypeTag>::
    wellVolumeFraction(const unsigned compIdx) const
    {
        if (FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx) && compIdx == Indices::canonicalToActiveComponentIndex(FluidSystem::waterCompIdx)) {
            return primary_variables_evaluation_[WFrac];
        }

        if (FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx) && compIdx == Indices::canonicalToActiveComponentIndex(FluidSystem::gasCompIdx)) {
            return primary_variables_evaluation_[GFrac];
        }

        if (has_solvent && compIdx == (unsigned)contiSolventEqIdx) {
            return primary_variables_evaluation_[SFrac];
        }

        // Oil fraction
        EvalWell well_fraction = 1.0;
        if (FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx)) {
            well_fraction -= primary_variables_evaluation_[WFrac];
        }

        if (FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx)) {
            well_fraction -= primary_variables_evaluation_[GFrac];
        }
        if (has_solvent) {
            well_fraction -= primary_variables_evaluation_[SFrac];
        }
        return well_fraction;
    }





    template<typename TypeTag>
    typename StandardWell<TypeTag>::EvalWell
    StandardWell<TypeTag>::
    wellSurfaceVolumeFraction(const int compIdx) const
    {

        EvalWell sum_volume_fraction_scaled = 0.;
        for (int idx = 0; idx < num_components_; ++idx) {
            sum_volume_fraction_scaled += wellVolumeFractionScaled(idx);
        }

        assert(sum_volume_fraction_scaled.value() != 0.);

        return wellVolumeFractionScaled(compIdx) / sum_volume_fraction_scaled;
     }





    template<typename TypeTag>
    typename StandardWell<TypeTag>::EvalWell
    StandardWell<TypeTag>::
    extendEval(const Eval& in) const
    {
        EvalWell out = 0.0;
        out.setValue(in.value());
        for(int eqIdx = 0; eqIdx < numEq;++eqIdx) {
            out.setDerivative(eqIdx, in.derivative(eqIdx));
        }
        return out;
    }





    template<typename TypeTag>
    void
    StandardWell<TypeTag>::
    computePerfRate(const IntensiveQuantities& intQuants,
                    const std::vector<EvalWell>& mob,
                    const EvalWell& bhp,
                    const double Tw,
                    const int perf,
                    const bool allow_cf,
                    std::vector<EvalWell>& cq_s,
                    double& perf_dis_gas_rate,
                    double& perf_vap_oil_rate,
                    Opm::DeferredLogger& deferred_logger) const
    {
        const auto& fs = intQuants.fluidState();
        const EvalWell pressure = extendEval(fs.pressure(FluidSystem::oilPhaseIdx));
        const EvalWell rs = extendEval(fs.Rs());
        const EvalWell rv = extendEval(fs.Rv());
        std::vector<EvalWell> b_perfcells_dense(num_components_, 0.0);
        for (unsigned phaseIdx = 0; phaseIdx < FluidSystem::numPhases; ++phaseIdx) {
            if (!FluidSystem::phaseIsActive(phaseIdx)) {
                continue;
            }

            const unsigned compIdx = Indices::canonicalToActiveComponentIndex(FluidSystem::solventComponentIndex(phaseIdx));
            b_perfcells_dense[compIdx] = extendEval(fs.invB(phaseIdx));
        }
        if (has_solvent) {
            b_perfcells_dense[contiSolventEqIdx] = extendEval(intQuants.solventInverseFormationVolumeFactor());
        }

        // Pressure drawdown (also used to determine direction of flow)
        const EvalWell well_pressure = bhp + perf_pressure_diffs_[perf];
        const EvalWell drawdown = pressure - well_pressure;

        // producing perforations
        if ( drawdown.value() > 0 )  {
            //Do nothing if crossflow is not allowed
            if (!allow_cf && well_type_ == INJECTOR) {
                return;
            }

            // compute component volumetric rates at standard conditions
            for (int componentIdx = 0; componentIdx < num_components_; ++componentIdx) {
                const EvalWell cq_p = - Tw * (mob[componentIdx] * drawdown);
                cq_s[componentIdx] = b_perfcells_dense[componentIdx] * cq_p;
            }

            if (FluidSystem::phaseIsActive(FluidSystem::oilPhaseIdx) && FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx)) {
                const unsigned oilCompIdx = Indices::canonicalToActiveComponentIndex(FluidSystem::oilCompIdx);
                const unsigned gasCompIdx = Indices::canonicalToActiveComponentIndex(FluidSystem::gasCompIdx);
                const EvalWell cq_sOil = cq_s[oilCompIdx];
                const EvalWell cq_sGas = cq_s[gasCompIdx];
                const EvalWell dis_gas = rs * cq_sOil;
                const EvalWell vap_oil = rv * cq_sGas;

                cq_s[gasCompIdx] += dis_gas;
                cq_s[oilCompIdx] += vap_oil;

                // recording the perforation solution gas rate and solution oil rates
                if (well_type_ == PRODUCER) {
                    perf_dis_gas_rate = dis_gas.value();
                    perf_vap_oil_rate = vap_oil.value();
                }
            }

        } else {
            //Do nothing if crossflow is not allowed
            if (!allow_cf && well_type_ == PRODUCER) {
                return;
            }

            // Using total mobilities
            EvalWell total_mob_dense = mob[0];
            for (int componentIdx = 1; componentIdx < num_components_; ++componentIdx) {
                total_mob_dense += mob[componentIdx];
            }

            // injection perforations total volume rates
            const EvalWell cqt_i = - Tw * (total_mob_dense * drawdown);

            // surface volume fraction of fluids within wellbore
            std::vector<EvalWell> cmix_s(num_components_, 0.);
            for (int componentIdx = 0; componentIdx < num_components_; ++componentIdx) {
                cmix_s[componentIdx] = wellSurfaceVolumeFraction(componentIdx);
            }

            // compute volume ratio between connection at standard conditions
            EvalWell volumeRatio = 0.0;
            if (FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx)) {
                const unsigned waterCompIdx = Indices::canonicalToActiveComponentIndex(FluidSystem::waterCompIdx);
                volumeRatio += cmix_s[waterCompIdx] / b_perfcells_dense[waterCompIdx];
            }

            if (has_solvent) {
                volumeRatio += cmix_s[contiSolventEqIdx] / b_perfcells_dense[contiSolventEqIdx];
            }

            if (FluidSystem::phaseIsActive(FluidSystem::oilPhaseIdx) && FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx)) {
                const unsigned oilCompIdx = Indices::canonicalToActiveComponentIndex(FluidSystem::oilCompIdx);
                const unsigned gasCompIdx = Indices::canonicalToActiveComponentIndex(FluidSystem::gasCompIdx);
                // Incorporate RS/RV factors if both oil and gas active
                const EvalWell d = 1.0 - rv * rs;

                if (d.value() == 0.0) {
                    OPM_DEFLOG_THROW(Opm::NumericalIssue, "Zero d value obtained for well " << name() << " during flux calcuation"
                                                  << " with rs " << rs << " and rv " << rv, deferred_logger);
                }

                const EvalWell tmp_oil = (cmix_s[oilCompIdx] - rv * cmix_s[gasCompIdx]) / d;
                //std::cout << "tmp_oil " <<tmp_oil << std::endl;
                volumeRatio += tmp_oil / b_perfcells_dense[oilCompIdx];

                const EvalWell tmp_gas = (cmix_s[gasCompIdx] - rs * cmix_s[oilCompIdx]) / d;
                //std::cout << "tmp_gas " <<tmp_gas << std::endl;
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
            EvalWell cqt_is = cqt_i/volumeRatio;
            //std::cout << "volrat " << volumeRatio << " " << volrat_perf_[perf] << std::endl;
            for (int componentIdx = 0; componentIdx < num_components_; ++componentIdx) {
                cq_s[componentIdx] = cmix_s[componentIdx] * cqt_is; // * b_perfcells_dense[phase];
            }

            // calculating the perforation solution gas rate and solution oil rates
            if (well_type_ == PRODUCER) {
                if (FluidSystem::phaseIsActive(FluidSystem::oilPhaseIdx) && FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx)) {
                    const unsigned oilCompIdx = Indices::canonicalToActiveComponentIndex(FluidSystem::oilCompIdx);
                    const unsigned gasCompIdx = Indices::canonicalToActiveComponentIndex(FluidSystem::gasCompIdx);
                    // TODO: the formulations here remain to be tested with cases with strong crossflow through production wells
                    // s means standard condition, r means reservoir condition
                    // q_os = q_or * b_o + rv * q_gr * b_g
                    // q_gs = q_gr * g_g + rs * q_or * b_o
                    // d = 1.0 - rs * rv
                    // q_or = 1 / (b_o * d) * (q_os - rv * q_gs)
                    // q_gr = 1 / (b_g * d) * (q_gs - rs * q_os)

                    const double d = 1.0 - rv.value() * rs.value();
                    // vaporized oil into gas
                    // rv * q_gr * b_g = rv * (q_gs - rs * q_os) / d
                    perf_vap_oil_rate = rv.value() * (cq_s[gasCompIdx].value() - rs.value() * cq_s[oilCompIdx].value()) / d;
                    // dissolved of gas in oil
                    // rs * q_or * b_o = rs * (q_os - rv * q_gs) / d
                    perf_dis_gas_rate = rs.value() * (cq_s[oilCompIdx].value() - rv.value() * cq_s[gasCompIdx].value()) / d;
                }
            }
        }
    }





    template<typename TypeTag>
    void
    StandardWell<TypeTag>::
    assembleWellEq(const Simulator& ebosSimulator,
                   const std::vector<Scalar>& /* B_avg */,
                   const double dt,
                   WellState& well_state,
                   Opm::DeferredLogger& deferred_logger)
    {
        SummaryState summaryState;

        checkWellOperability(ebosSimulator, well_state, deferred_logger);

        if (!this->isOperable()) return;

        // clear all entries
        duneB_ = 0.0;
        duneC_ = 0.0;
        invDuneD_ = 0.0;
        resWell_ = 0.0;

        // TODO: it probably can be static member for StandardWell
        const double volume = 0.002831684659200; // 0.1 cu ft;

        const bool allow_cf = getAllowCrossFlow() || openCrossFlowAvoidSingularity(ebosSimulator);

        const EvalWell& bhp = getBhp();

        // the solution gas rate and solution oil rate needs to be reset to be zero for well_state.
        well_state.wellVaporizedOilRates()[index_of_well_] = 0.;
        well_state.wellDissolvedGasRates()[index_of_well_] = 0.;

        const int np = number_of_phases_;
        for (int p = 0; p < np; ++p) {
            well_state.productivityIndex()[np*index_of_well_ + p] = 0.;
        }

        for (int perf = 0; perf < number_of_perforations_; ++perf) {

            const int cell_idx = well_cells_[perf];
            const auto& intQuants = *(ebosSimulator.model().cachedIntensiveQuantities(cell_idx, /*timeIdx=*/ 0));
            std::vector<EvalWell> mob(num_components_, 0.0);
            getMobility(ebosSimulator, perf, mob, deferred_logger);

            std::vector<EvalWell> cq_s(num_components_, 0.0);
            double perf_dis_gas_rate = 0.;
            double perf_vap_oil_rate = 0.;
            double trans_mult = ebosSimulator.problem().template rockCompTransMultiplier<double>(intQuants,  cell_idx);
            const double Tw = well_index_[perf] * trans_mult;
            computePerfRate(intQuants, mob, bhp, Tw, perf, allow_cf,
                            cq_s, perf_dis_gas_rate, perf_vap_oil_rate, deferred_logger);

            // updating the solution gas rate and solution oil rate
            if (well_type_ == PRODUCER) {
                well_state.wellDissolvedGasRates()[index_of_well_] += perf_dis_gas_rate;
                well_state.wellVaporizedOilRates()[index_of_well_] += perf_vap_oil_rate;
            }

            if (has_energy) {
                connectionRates_[perf][contiEnergyEqIdx] = 0.0;
            }

            for (int componentIdx = 0; componentIdx < num_components_; ++componentIdx) {
                // the cq_s entering mass balance equations need to consider the efficiency factors.
                const EvalWell cq_s_effective = cq_s[componentIdx] * well_efficiency_factor_;

                connectionRates_[perf][componentIdx] = Base::restrictEval(cq_s_effective);

                // subtract sum of phase fluxes in the well equations.
                resWell_[0][componentIdx] += cq_s_effective.value();

                // assemble the jacobians
                for (int pvIdx = 0; pvIdx < numWellEq; ++pvIdx) {
                    // also need to consider the efficiency factor when manipulating the jacobians.
                    duneC_[0][cell_idx][pvIdx][componentIdx] -= cq_s_effective.derivative(pvIdx+numEq); // intput in transformed matrix
                    invDuneD_[0][0][componentIdx][pvIdx] += cq_s_effective.derivative(pvIdx+numEq);
                }

                for (int pvIdx = 0; pvIdx < numEq; ++pvIdx) {
                    duneB_[0][cell_idx][componentIdx][pvIdx] += cq_s_effective.derivative(pvIdx);
                }

                // Store the perforation phase flux for later usage.
                if (has_solvent && componentIdx == contiSolventEqIdx) {
                    well_state.perfRateSolvent()[first_perf_ + perf] = cq_s[componentIdx].value();
                } else {
                    well_state.perfPhaseRates()[(first_perf_ + perf) * np + ebosCompIdxToFlowCompIdx(componentIdx)] = cq_s[componentIdx].value();
                }
            }
            if (has_energy) {

                auto fs = intQuants.fluidState();
                for (unsigned phaseIdx = 0; phaseIdx < FluidSystem::numPhases; ++phaseIdx) {
                    if (!FluidSystem::phaseIsActive(phaseIdx)) {
                        continue;
                    }

                    const unsigned activeCompIdx = Indices::canonicalToActiveComponentIndex(FluidSystem::solventComponentIndex(phaseIdx));
                    // convert to reservoar conditions
                    EvalWell cq_r_thermal = 0.0;
                    if (FluidSystem::phaseIsActive(FluidSystem::oilPhaseIdx) && FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx)) {

                        if(FluidSystem::waterPhaseIdx == phaseIdx)
                             cq_r_thermal = cq_s[activeCompIdx] / extendEval(fs.invB(phaseIdx));

                        // remove dissolved gas and vapporized oil
                        const unsigned oilCompIdx = Indices::canonicalToActiveComponentIndex(FluidSystem::oilCompIdx);
                        const unsigned gasCompIdx = Indices::canonicalToActiveComponentIndex(FluidSystem::gasCompIdx);
                        // q_os = q_or * b_o + rv * q_gr * b_g
                        // q_gs = q_gr * g_g + rs * q_or * b_o
                        // d = 1.0 - rs * rv
                        const EvalWell d = extendEval(1.0 - fs.Rv() * fs.Rs());
                        // q_gr = 1 / (b_g * d) * (q_gs - rs * q_os)
                        if(FluidSystem::gasPhaseIdx == phaseIdx)
                            cq_r_thermal = (cq_s[gasCompIdx] - extendEval(fs.Rs()) * cq_s[oilCompIdx]) / (d * extendEval(fs.invB(phaseIdx)) );
                        // q_or = 1 / (b_o * d) * (q_os - rv * q_gs)
                        if(FluidSystem::oilPhaseIdx == phaseIdx)
                            cq_r_thermal = (cq_s[oilCompIdx] - extendEval(fs.Rv()) * cq_s[gasCompIdx]) / (d * extendEval(fs.invB(phaseIdx)) );

                    } else {
                        cq_r_thermal = cq_s[activeCompIdx] / extendEval(fs.invB(phaseIdx));
                    }

                    // change temperature for injecting fluids
                    if (well_type_ == INJECTOR && cq_s[activeCompIdx] > 0.0){
                        const auto& injProps = this->well_ecl_.injectionControls(summaryState);

                        // only handles single phase injection now
                        assert(injProps.injector_type != WellInjector::MULTI);

                        fs.setTemperature(injProps.temperature);
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
                    }
                    // compute the thermal flux
                    cq_r_thermal *= extendEval(fs.enthalpy(phaseIdx)) * extendEval(fs.density(phaseIdx));
                    connectionRates_[perf][contiEnergyEqIdx] += Base::restrictEval(cq_r_thermal);
                }
            }

            if (has_polymer) {
                // TODO: the application of well efficiency factor has not been tested with an example yet
                const unsigned waterCompIdx = Indices::canonicalToActiveComponentIndex(FluidSystem::waterCompIdx);
                EvalWell cq_s_poly = cq_s[waterCompIdx] * well_efficiency_factor_;
                if (well_type_ == INJECTOR) {
                    cq_s_poly *= wpolymer();
                } else {
                    cq_s_poly *= extendEval(intQuants.polymerConcentration() * intQuants.polymerViscosityCorrection());
                }
                connectionRates_[perf][contiPolymerEqIdx] = Base::restrictEval(cq_s_poly);

                if (this->has_polymermw) {
                    if (well_type_ == PRODUCER) {
                        EvalWell cq_s_polymw = cq_s_poly;
                        if (cq_s_polymw < 0.) {
                            cq_s_polymw *= extendEval(intQuants.polymerMoleWeight() );
                        } else {
                            // we do not consider the molecular weight from the polymer
                            // re-injecting back through producer
                            cq_s_polymw *= 0.;
                        }
                        connectionRates_[perf][this->contiPolymerMWEqIdx] = Base::restrictEval(cq_s_polymw);
                    }
                }
            }

            // Store the perforation pressure for later usage.
            well_state.perfPress()[first_perf_ + perf] = well_state.bhp()[index_of_well_] + perf_pressure_diffs_[perf];

            // Compute Productivity index if asked for
            const auto& pu = phaseUsage();
            const Opm::SummaryConfig& summaryConfig = ebosSimulator.vanguard().summaryConfig();
            const Opm::Schedule& schedule = ebosSimulator.vanguard().schedule();
            for (int p = 0; p < np; ++p) {
                if ( (pu.phase_pos[Water] == p && (summaryConfig.hasSummaryKey("WPIW:" + name()) || summaryConfig.hasSummaryKey("WPIL:" + name())))
                        || (pu.phase_pos[Oil] == p && (summaryConfig.hasSummaryKey("WPIO:" + name()) || summaryConfig.hasSummaryKey("WPIL:" + name())))
                        || (pu.phase_pos[Gas] == p && summaryConfig.hasSummaryKey("WPIG:" + name()))) {

                    const unsigned int compIdx = flowPhaseToEbosCompIdx(p);
                    const double drawdown = well_state.perfPress()[first_perf_ + perf] - intQuants.fluidState().pressure(FluidSystem::oilPhaseIdx).value();
                    const bool new_well = schedule.hasWellEvent(name(), ScheduleEvents::NEW_WELL, current_step_);
                    double productivity_index = cq_s[compIdx].value() / drawdown;
                    scaleProductivityIndex(perf, productivity_index, new_well, deferred_logger);
                    well_state.productivityIndex()[np*index_of_well_ + p] += productivity_index;
                }
            }

        }


        // add vol * dF/dt + Q to the well equations;
        for (int componentIdx = 0; componentIdx < numWellConservationEq; ++componentIdx) {
            // TODO: following the development in MSW, we need to convert the volume of the wellbore to be surface volume
            // since all the rates are under surface condition
            EvalWell resWell_loc = (wellSurfaceVolumeFraction(componentIdx) - F0_[componentIdx]) * volume / dt;
            resWell_loc -= getQs(componentIdx) * well_efficiency_factor_;
            for (int pvIdx = 0; pvIdx < numWellEq; ++pvIdx) {
                invDuneD_[0][0][componentIdx][pvIdx] += resWell_loc.derivative(pvIdx+numEq);
            }
            resWell_[0][componentIdx] += resWell_loc.value();
        }

        assembleControlEq(deferred_logger);

        // do the local inversion of D.
        try {
            Dune::ISTLUtility::invertMatrix(invDuneD_[0][0]);
        } catch( ... ) {
            OPM_DEFLOG_THROW(Opm::NumericalIssue,"Error when inverting local well equations for well " + name(), deferred_logger);
        }


    }





    template <typename TypeTag>
    void
    StandardWell<TypeTag>::
    assembleControlEq(Opm::DeferredLogger& deferred_logger)
    {
        EvalWell control_eq(0.0);
        switch (well_controls_get_current_type(well_controls_)) {
            case THP:
            {
                std::vector<EvalWell> rates(3, 0.);
                if (FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx)) {
                    rates[ Water ] = getQs(flowPhaseToEbosCompIdx(Water));
                }
                if (FluidSystem::phaseIsActive(FluidSystem::oilPhaseIdx)) {
                    rates[ Oil ] = getQs(flowPhaseToEbosCompIdx(Oil));
                }
                if (FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx)) {
                    rates[ Gas ] = getQs(flowPhaseToEbosCompIdx(Gas));
                }
                const int current = well_controls_get_current(well_controls_);
                control_eq = getBhp() - calculateBhpFromThp(rates, current, deferred_logger);
                break;
            }
            case BHP:
            {
                const double target_bhp = well_controls_get_current_target(well_controls_);
                control_eq = getBhp() - target_bhp;
                break;
            }
            case SURFACE_RATE:
            {
                const double target_rate = well_controls_get_current_target(well_controls_); // surface rate target
                if (well_type_ == INJECTOR) {
                    control_eq = getWQTotal() - target_rate;
                } else if (well_type_ == PRODUCER) {
                    if (target_rate != 0.) {
                        EvalWell rate_for_control(0.);
                        const EvalWell& g_total = getWQTotal();
                        // a variable to check if we are producing any targeting fluids
                        double sum_fraction = 0.;
                        const double* distr = well_controls_get_current_distr(well_controls_);
                        for (int phase = 0; phase < number_of_phases_; ++phase) {
                            if (distr[phase] > 0.) {
                                const EvalWell fraction_scaled = wellVolumeFractionScaled(flowPhaseToEbosCompIdx(phase));
                                rate_for_control += g_total * fraction_scaled;
                                sum_fraction += fraction_scaled.value();
                            }
                        }
                        if (sum_fraction > 0.) {
                            control_eq = rate_for_control - target_rate;
                        } else {
                            // we are not producing any fluids that specfied for a non-zero target
                            // which makes it a mission impossible, we will set all the rates to be zero for this case
                            const std::string msg = " Setting all rates to be zero for well " + name()
                                                  + " due to un-solvable situation. There is non-zero target for the phase "
                                                  + " that does not exist in the wellbore for the situation";
                            deferred_logger.warning("NON_SOLVABLE_WELL_SOLUTION", msg);

                            control_eq = getWQTotal() - target_rate;
                        }
                    } else {
                        // there is some special treatment for the zero rate control well
                        // 1. if the well can produce the specified phase, it means the well should not produce any fluid, this
                        //    is a fine situation.
                        // 2. if the well can not produce the specified phase, it cause a under-determined problem, we
                        //    basically assume the well not producing any fluid as a solution
                        // With both the situation, we can use the following well equation
                        control_eq = getWQTotal() - target_rate;
                    }
                }
                break;
            }
            case RESERVOIR_RATE:
            {
                const double target_rate = well_controls_get_current_target(well_controls_); // reservoir rate target
                if (well_type_ == INJECTOR) {
                    const double* distr = well_controls_get_current_distr(well_controls_);
                    for (int phase = 0; phase < number_of_phases_; ++phase) {
                        if (distr[phase] > 0.0) {
                            control_eq = getWQTotal() * scalingFactor(phase) - target_rate;
                            break;
                        }
                    }
                } else {
                    const EvalWell& g_total = getWQTotal();
                    EvalWell rate_for_control(0.0); // reservoir rate
                    for (int phase = 0; phase < number_of_phases_; ++phase) {
                        rate_for_control += g_total * wellVolumeFraction( flowPhaseToEbosCompIdx(phase) );
                    }
                    control_eq = rate_for_control - target_rate;
                }
                break;
            }
            default:
                OPM_DEFLOG_THROW(std::runtime_error, "Unknown well control control types for well " << name(), deferred_logger);
        }

        // using control_eq to update the matrix and residuals
        // TODO: we should use a different index system for the well equations
        resWell_[0][Bhp] = control_eq.value();
        for (int pv_idx = 0; pv_idx < numWellEq; ++pv_idx) {
            invDuneD_[0][0][Bhp][pv_idx] = control_eq.derivative(pv_idx + numEq);
        }
    }





    template<typename TypeTag>
    void
    StandardWell<TypeTag>::
    getMobility(const Simulator& ebosSimulator,
                const int perf,
                std::vector<EvalWell>& mob,
                Opm::DeferredLogger& deferred_logger) const
    {
        const int cell_idx = well_cells_[perf];
        assert (int(mob.size()) == num_components_);
        const auto& intQuants = *(ebosSimulator.model().cachedIntensiveQuantities(cell_idx, /*timeIdx=*/0));
        const auto& materialLawManager = ebosSimulator.problem().materialLawManager();

        // either use mobility of the perforation cell or calcualte its own
        // based on passing the saturation table index
        const int satid = saturation_table_number_[perf] - 1;
        const int satid_elem = materialLawManager->satnumRegionIdx(cell_idx);
        if( satid == satid_elem ) { // the same saturation number is used. i.e. just use the mobilty from the cell

            for (unsigned phaseIdx = 0; phaseIdx < FluidSystem::numPhases; ++phaseIdx) {
                if (!FluidSystem::phaseIsActive(phaseIdx)) {
                    continue;
                }

                const unsigned activeCompIdx = Indices::canonicalToActiveComponentIndex(FluidSystem::solventComponentIndex(phaseIdx));
                mob[activeCompIdx] = extendEval(intQuants.mobility(phaseIdx));
            }
            if (has_solvent) {
                mob[contiSolventEqIdx] = extendEval(intQuants.solventMobility());
            }
        } else {

            const auto& paramsCell = materialLawManager->connectionMaterialLawParams(satid, cell_idx);
            Eval relativePerms[3] = { 0.0, 0.0, 0.0 };
            MaterialLaw::relativePermeabilities(relativePerms, paramsCell, intQuants.fluidState());

            // reset the satnumvalue back to original
            materialLawManager->connectionMaterialLawParams(satid_elem, cell_idx);

            // compute the mobility
            for (unsigned phaseIdx = 0; phaseIdx < FluidSystem::numPhases; ++phaseIdx) {
                if (!FluidSystem::phaseIsActive(phaseIdx)) {
                    continue;
                }

                const unsigned activeCompIdx = Indices::canonicalToActiveComponentIndex(FluidSystem::solventComponentIndex(phaseIdx));
                mob[activeCompIdx] = extendEval(relativePerms[phaseIdx] / intQuants.fluidState().viscosity(phaseIdx));
            }

            // this may not work if viscosity and relperms has been modified?
            if (has_solvent) {
                OPM_DEFLOG_THROW(std::runtime_error, "individual mobility for wells does not work in combination with solvent", deferred_logger);
            }
        }

        // modify the water mobility if polymer is present
        if (has_polymer) {
            if (!FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx)) {
                OPM_DEFLOG_THROW(std::runtime_error, "Water is required when polymer is active", deferred_logger);
            }

            updateWaterMobilityWithPolymer(ebosSimulator, perf, mob, deferred_logger);
        }
    }





    template<typename TypeTag>
    void
    StandardWell<TypeTag>::
    updateWellState(const BVectorWell& dwells,
                    WellState& well_state,
                    Opm::DeferredLogger& deferred_logger) const
    {
        if (!this->isOperable()) return;

        updatePrimaryVariablesNewton(dwells, well_state);

        updateWellStateFromPrimaryVariables(well_state, deferred_logger);
    }





    template<typename TypeTag>
    void
    StandardWell<TypeTag>::
    updatePrimaryVariablesNewton(const BVectorWell& dwells,
                                 const WellState& /* well_state */) const
    {
        const double dFLimit = param_.dwell_fraction_max_;

        const std::vector<double> old_primary_variables = primary_variables_;

        // for injectors, very typical one of the fractions will be one, and it is easy to get zero value
        // fractions. not sure what is the best way to handle it yet, so we just use 1.0 here
        const double relaxation_factor_fractions = (well_type_ == PRODUCER) ?
                                         relaxationFactorFractionsProducer(old_primary_variables, dwells)
                                       : 1.0;

        // update the second and third well variable (The flux fractions)
        if (FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx)) {
            const int sign2 = dwells[0][WFrac] > 0 ? 1: -1;
            const double dx2_limited = sign2 * std::min(std::abs(dwells[0][WFrac] * relaxation_factor_fractions), dFLimit);
            // primary_variables_[WFrac] = old_primary_variables[WFrac] - dx2_limited;
            primary_variables_[WFrac] = old_primary_variables[WFrac] - dx2_limited;
        }

        if (FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx)) {
            const int sign3 = dwells[0][GFrac] > 0 ? 1: -1;
            const double dx3_limited = sign3 * std::min(std::abs(dwells[0][GFrac] * relaxation_factor_fractions), dFLimit);
            primary_variables_[GFrac] = old_primary_variables[GFrac] - dx3_limited;
        }

        if (has_solvent) {
            const int sign4 = dwells[0][SFrac] > 0 ? 1: -1;
            const double dx4_limited = sign4 * std::min(std::abs(dwells[0][SFrac]) * relaxation_factor_fractions, dFLimit);
            primary_variables_[SFrac] = old_primary_variables[SFrac] - dx4_limited;
        }

        processFractions();

        // updating the total rates Q_t
        const double relaxation_factor_rate = relaxationFactorRate(old_primary_variables, dwells);
        primary_variables_[WQTotal] = old_primary_variables[WQTotal] - dwells[0][WQTotal] * relaxation_factor_rate;

        // updating the bottom hole pressure
        {
            const double dBHPLimit = param_.dbhp_max_rel_;
            const int sign1 = dwells[0][Bhp] > 0 ? 1: -1;
            const double dx1_limited = sign1 * std::min(std::abs(dwells[0][Bhp]), std::abs(old_primary_variables[Bhp]) * dBHPLimit);
            // 1e5 to make sure bhp will not be below 1bar
            primary_variables_[Bhp] = std::max(old_primary_variables[Bhp] - dx1_limited, 1e5);
        }
    }





    template<typename TypeTag>
    void
    StandardWell<TypeTag>::
    processFractions() const
    {
        assert(FluidSystem::phaseIsActive(FluidSystem::oilPhaseIdx));
        const auto pu = phaseUsage();
        std::vector<double> F(number_of_phases_, 0.0);
        F[pu.phase_pos[Oil]] = 1.0;

        if (FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx)) {
            F[pu.phase_pos[Water]] = primary_variables_[WFrac];
            F[pu.phase_pos[Oil]] -= F[pu.phase_pos[Water]];
        }

        if (FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx)) {
            F[pu.phase_pos[Gas]] = primary_variables_[GFrac];
            F[pu.phase_pos[Oil]] -= F[pu.phase_pos[Gas]];
        }

        double F_solvent = 0.0;
        if (has_solvent) {
            F_solvent = primary_variables_[SFrac];
            F[pu.phase_pos[Oil]] -= F_solvent;
        }

        if (FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx)) {
            if (F[Water] < 0.0) {
                if (FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx)) {
                        F[pu.phase_pos[Gas]] /= (1.0 - F[pu.phase_pos[Water]]);
                }
                if (has_solvent) {
                    F_solvent /= (1.0 - F[pu.phase_pos[Water]]);
                }
                F[pu.phase_pos[Oil]] /= (1.0 - F[pu.phase_pos[Water]]);
                F[pu.phase_pos[Water]] = 0.0;
            }
        }

        if (FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx)) {
            if (F[pu.phase_pos[Gas]] < 0.0) {
                if (FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx)) {
                    F[pu.phase_pos[Water]] /= (1.0 - F[pu.phase_pos[Gas]]);
                }
                if (has_solvent) {
                    F_solvent /= (1.0 - F[pu.phase_pos[Gas]]);
                }
                F[pu.phase_pos[Oil]] /= (1.0 - F[pu.phase_pos[Gas]]);
                F[pu.phase_pos[Gas]] = 0.0;
            }
        }

        if (F[pu.phase_pos[Oil]] < 0.0) {
            if (FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx)) {
                F[pu.phase_pos[Water]] /= (1.0 - F[pu.phase_pos[Oil]]);
            }
            if (FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx)) {
                F[pu.phase_pos[Gas]] /= (1.0 - F[pu.phase_pos[Oil]]);
            }
            if (has_solvent) {
                F_solvent /= (1.0 - F[pu.phase_pos[Oil]]);
            }
            F[pu.phase_pos[Oil]] = 0.0;
        }

        if (FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx)) {
            primary_variables_[WFrac] = F[pu.phase_pos[Water]];
        }
        if (FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx)) {
            primary_variables_[GFrac] = F[pu.phase_pos[Gas]];
        }
        if(has_solvent) {
            primary_variables_[SFrac] = F_solvent;
        }
    }





    template<typename TypeTag>
    void
    StandardWell<TypeTag>::
    updateWellStateFromPrimaryVariables(WellState& well_state, Opm::DeferredLogger& deferred_logger) const
    {
        const PhaseUsage& pu = phaseUsage();
        assert( FluidSystem::phaseIsActive(FluidSystem::oilPhaseIdx) );
        const int oil_pos = pu.phase_pos[Oil];

        std::vector<double> F(number_of_phases_, 0.0);
        F[oil_pos] = 1.0;

        if ( FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx) ) {
            const int water_pos = pu.phase_pos[Water];
            F[water_pos] = primary_variables_[WFrac];
            F[oil_pos] -= F[water_pos];
        }

        if ( FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx) ) {
            const int gas_pos = pu.phase_pos[Gas];
            F[gas_pos] = primary_variables_[GFrac];
            F[oil_pos] -= F[gas_pos];
        }

        double F_solvent = 0.0;
        if (has_solvent) {
            F_solvent = primary_variables_[SFrac];
            F[oil_pos] -= F_solvent;
        }

        // convert the fractions to be Q_p / G_total to calculate the phase rates
        for (int p = 0; p < number_of_phases_; ++p) {
            const double scal = scalingFactor(p);
            // for injection wells, there should only one non-zero scaling factor
            if (scal > 0) {
                F[p] /= scal ;
            } else {
                // this should only happens to injection wells
                F[p] = 0.;
            }
        }

        // F_solvent is added to F_gas. This means that well_rate[Gas] also contains solvent.
        // More testing is needed to make sure this is correct for well groups and THP.
        if (has_solvent){
            F_solvent /= scalingFactor(contiSolventEqIdx);
            F[pu.phase_pos[Gas]] += F_solvent;
        }

        well_state.bhp()[index_of_well_] = primary_variables_[Bhp];

        // calculate the phase rates based on the primary variables
        // for producers, this is not a problem, while not sure for injectors here
        if (well_type_ == PRODUCER) {
            const double g_total = primary_variables_[WQTotal];
            for (int p = 0; p < number_of_phases_; ++p) {
                well_state.wellRates()[index_of_well_ * number_of_phases_ + p] = g_total * F[p];
            }
        } else { // injectors
            // TODO: using comp_frac_ here is very dangerous, since we do not update it based on the injection phase
            // Either we use distr (might conflict with RESV related) or we update comp_frac_ based on the injection phase
            for (int p = 0; p < number_of_phases_; ++p) {
                const double comp_frac = comp_frac_[p];
                well_state.wellRates()[index_of_well_ * number_of_phases_ + p] = comp_frac * primary_variables_[WQTotal];
            }
        }

        updateThp(well_state, deferred_logger);
    }





    template<typename TypeTag>
    void
    StandardWell<TypeTag>::
    updateThp(WellState& well_state, Opm::DeferredLogger& deferred_logger) const
    {
        // When there is no vaild VFP table provided, we set the thp to be zero.
        if (!this->isVFPActive(deferred_logger)) {
            well_state.thp()[index_of_well_] = 0.;
            return;
        }

        // avaiable VFP table is provided, we should update the thp value

        // if the well is under THP control, we should use its target value
        if (well_controls_get_current_type(well_controls_) == THP) {
            well_state.thp()[index_of_well_] = well_controls_get_current_target(well_controls_);
        } else {
            // the well is under other control types, we calculate the thp based on bhp and rates
            std::vector<double> rates(3, 0.0);

            const Opm::PhaseUsage& pu = phaseUsage();
            if (FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx)) {
                rates[ Water ] = well_state.wellRates()[index_of_well_ * number_of_phases_ + pu.phase_pos[ Water ] ];
            }
            if (FluidSystem::phaseIsActive(FluidSystem::oilPhaseIdx)) {
                rates[ Oil ] = well_state.wellRates()[index_of_well_ * number_of_phases_ + pu.phase_pos[ Oil ] ];
            }
            if (FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx)) {
                rates[ Gas ] = well_state.wellRates()[index_of_well_ * number_of_phases_ + pu.phase_pos[ Gas ] ];
            }

            const double bhp = well_state.bhp()[index_of_well_];

            well_state.thp()[index_of_well_] = calculateThpFromBhp(rates, bhp, deferred_logger);
        }
    }





    template<typename TypeTag>
    void
    StandardWell<TypeTag>::
    updateWellStateWithTarget(const Simulator& ebos_simulator,
                              WellState& well_state,
                              Opm::DeferredLogger& deferred_logger) const
    {
        
        // number of phases
        const int np = number_of_phases_;
        const int well_index = index_of_well_;
        const WellControls* wc = well_controls_;
        const int current = well_state.currentControls()[well_index];
        // Updating well state and primary variables.
        // Target values are used as initial conditions for BHP, THP, and SURFACE_RATE
        const double target = well_controls_iget_target(wc, current);
        const double* distr = well_controls_iget_distr(wc, current);
        switch (well_controls_iget_type(wc, current)) {
        case BHP:
            well_state.bhp()[well_index] = target;
            // TODO: similar to the way below to handle THP
            // we should not something related to thp here when there is thp constraint
            // or when can calculate the THP (table avaiable or requested for output?)
            // TODO: we should address this in a function updateWellStateWithBHPTarget.
            // TODO: however, the reason that this one minght not be that critical with
            // TODO: the effects remaining to be investigated.
            break;

        case THP: {
            // when a well can not work under THP target, it switches to BHP control
            if (this->operability_status_.isOperableUnderTHPLimit() ) {
                updateWellStateWithTHPTargetIPR(ebos_simulator, well_state, deferred_logger);
            } else { // go to BHP limit
                assert(this->operability_status_.isOperableUnderBHPLimit() );

                deferred_logger.info("well " + name() + " can not work with THP target, switching to BHP control");

                well_state.bhp()[well_index] = mostStrictBhpFromBhpLimits(deferred_logger);
            }
            break;
        }

        case RESERVOIR_RATE: // intentional fall-through
        case SURFACE_RATE:
            // TODO: something needs to be done with BHP and THP here
            // TODO: they should go to a separate function
            // checking the number of the phases under control
            int numPhasesWithTargetsUnderThisControl = 0;
            for (int phase = 0; phase < np; ++phase) {
                if (distr[phase] > 0.0) {
                    numPhasesWithTargetsUnderThisControl += 1;
                }
            }

            assert(numPhasesWithTargetsUnderThisControl > 0);

            if (well_type_ == INJECTOR) {
                // assign target value as initial guess for injectors
                // only handles single phase control at the moment
                assert(numPhasesWithTargetsUnderThisControl == 1);

                for (int phase = 0; phase < np; ++phase) {
                    if (distr[phase] > 0.) {
                        well_state.wellRates()[np*well_index + phase] = target / distr[phase];
                    } else {
                        well_state.wellRates()[np * well_index + phase] = 0.;
                    }
                }
            } else if (well_type_ == PRODUCER) {
                // update the rates of phases under control based on the target,
                // and also update rates of phases not under control to keep the rate ratio,
                // assuming the mobility ratio does not change for the production wells
                double original_rates_under_phase_control = 0.0;
                for (int phase = 0; phase < np; ++phase) {
                    if (distr[phase] > 0.0) {
                        original_rates_under_phase_control += well_state.wellRates()[np * well_index + phase] * distr[phase];
                    }
                }

                if (original_rates_under_phase_control != 0.0 ) {
                    const double scaling_factor = target / original_rates_under_phase_control;

                    for (int phase = 0; phase < np; ++phase) {
                        well_state.wellRates()[np * well_index + phase] *= scaling_factor;
                    }
                } else { // scaling factor is not well defined when original_rates_under_phase_control is zero
                    // separating targets equally between phases under control
                    const double target_rate_divided = target / numPhasesWithTargetsUnderThisControl;
                    for (int phase = 0; phase < np; ++phase) {
                        if (distr[phase] > 0.0) {
                            well_state.wellRates()[np * well_index + phase] = target_rate_divided / distr[phase];
                        } else {
                            // this only happens for SURFACE_RATE control
                            well_state.wellRates()[np * well_index + phase] = target_rate_divided;
                        }
                    }
                }
            } else {
                OPM_DEFLOG_THROW(std::logic_error, "Expected PRODUCER or INJECTOR type of well", deferred_logger);
            }

            break;
        } // end of switch
    }





    template<typename TypeTag>
    void
    StandardWell<TypeTag>::
    updateIPR(const Simulator& ebos_simulator, Opm::DeferredLogger& deferred_logger) const
    {
        // TODO: not handling solvent related here for now

        // TODO: it only handles the producers for now
        // the formular for the injectors are not formulated yet
        if (well_type_ == INJECTOR) {
            return;
        }

        // initialize all the values to be zero to begin with
        std::fill(ipr_a_.begin(), ipr_a_.end(), 0.);
        std::fill(ipr_b_.begin(), ipr_b_.end(), 0.);

        for (int perf = 0; perf < number_of_perforations_; ++perf) {
            std::vector<EvalWell> mob(num_components_, 0.0);
            // TODO: mabye we should store the mobility somewhere, so that we only need to calculate it one per iteration
            getMobility(ebos_simulator, perf, mob, deferred_logger);

            const int cell_idx = well_cells_[perf];
            const auto& int_quantities = *(ebos_simulator.model().cachedIntensiveQuantities(cell_idx, /*timeIdx=*/ 0));
            const auto& fs = int_quantities.fluidState();
            // the pressure of the reservoir grid block the well connection is in
            const double p_r = fs.pressure(FluidSystem::oilPhaseIdx).value();

            // calculating the b for the connection
            std::vector<double> b_perf(num_components_);
            for (size_t phase = 0; phase < FluidSystem::numPhases; ++phase) {
                if (!FluidSystem::phaseIsActive(phase)) {
                    continue;
                }
                const unsigned comp_idx = Indices::canonicalToActiveComponentIndex(FluidSystem::solventComponentIndex(phase));
                b_perf[comp_idx] = fs.invB(phase).value();
            }

            // the pressure difference between the connection and BHP
            const double h_perf = perf_pressure_diffs_[perf];
            const double pressure_diff = p_r - h_perf;

            // Let us add a check, since the pressure is calculated based on zero value BHP
            // it should not be negative anyway. If it is negative, we might need to re-formulate
            // to taking into consideration the crossflow here.
            if (pressure_diff <= 0.) {
                deferred_logger.warning("NON_POSITIVE_DRAWDOWN_IPR",
                                "non-positive drawdown found when updateIPR for well " + name());
            }

            // the well index associated with the connection
            const double tw_perf = well_index_[perf]*ebos_simulator.problem().template rockCompTransMultiplier<double>(int_quantities, cell_idx);

            // TODO: there might be some indices related problems here
            // phases vs components
            // ipr values for the perforation
            std::vector<double> ipr_a_perf(ipr_a_.size());
            std::vector<double> ipr_b_perf(ipr_b_.size());
            for (int p = 0; p < number_of_phases_; ++p) {
                const double tw_mob = tw_perf * mob[p].value() * b_perf[p];
                ipr_a_perf[p] += tw_mob * pressure_diff;
                ipr_b_perf[p] += tw_mob;
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

            for (int p = 0; p < number_of_phases_; ++p) {
                // TODO: double check the indices here
                ipr_a_[ebosCompIdxToFlowCompIdx(p)] += ipr_a_perf[p];
                ipr_b_[ebosCompIdxToFlowCompIdx(p)] += ipr_b_perf[p];
            }
        }
    }





    template<typename TypeTag>
    void
    StandardWell<TypeTag>::
    checkWellOperability(const Simulator& ebos_simulator,
                         const WellState& well_state,
                         Opm::DeferredLogger& deferred_logger
                         )
    {
        // focusing on PRODUCER for now
        if (well_type_ == INJECTOR) {
            return;
        }

        if (!this->underPredictionMode(deferred_logger) ) {
            return;
        }

        const bool old_well_operable = this->operability_status_.isOperable();

        updateWellOperability(ebos_simulator, well_state, deferred_logger);

        const bool well_operable = this->operability_status_.isOperable();

        if (!well_operable && old_well_operable) {
            deferred_logger.info(" well " + name() + " gets SHUT during iteration ");
        } else if (well_operable && !old_well_operable) {
            deferred_logger.info(" well " + name() + " gets REVIVED during iteration ");
        }
    }





    template<typename TypeTag>
    void
    StandardWell<TypeTag>::
    updateWellOperability(const Simulator& ebos_simulator,
                          const WellState& /* well_state */,
                          Opm::DeferredLogger& deferred_logger
                          )
    {
        this->operability_status_.reset();

        updateIPR(ebos_simulator, deferred_logger);

        // checking the BHP limit related
        checkOperabilityUnderBHPLimitProducer(ebos_simulator, deferred_logger);

        // checking whether the well can operate under the THP constraints.
        if (this->wellHasTHPConstraints()) {
            checkOperabilityUnderTHPLimitProducer(ebos_simulator, deferred_logger);
        }
    }





    template<typename TypeTag>
    void
    StandardWell<TypeTag>::
    checkOperabilityUnderBHPLimitProducer(const Simulator& ebos_simulator, Opm::DeferredLogger& deferred_logger)
    {
        const double bhp_limit = mostStrictBhpFromBhpLimits(deferred_logger);
        // Crude but works: default is one atmosphere.
        // TODO: a better way to detect whether the BHP is defaulted or not
        const bool bhp_limit_not_defaulted = bhp_limit > 1.5 * unit::barsa;
        if ( bhp_limit_not_defaulted || !this->wellHasTHPConstraints() ) {
            // if the BHP limit is not defaulted or the well does not have a THP limit
            // we need to check the BHP limit

            for (int p = 0; p < number_of_phases_; ++p) {
                const double temp = ipr_a_[p] - ipr_b_[p] * bhp_limit;
                if (temp < 0.) {
                    this->operability_status_.operable_under_only_bhp_limit = false;
                    break;
                }
            }

            // checking whether running under BHP limit will violate THP limit
            if (this->operability_status_.operable_under_only_bhp_limit && this->wellHasTHPConstraints()) {
                // option 1: calculate well rates based on the BHP limit.
                // option 2: stick with the above IPR curve
                // we use IPR here
                std::vector<double> well_rates_bhp_limit;
                computeWellRatesWithBhp(ebos_simulator, bhp_limit, well_rates_bhp_limit, deferred_logger);

                const double thp = calculateThpFromBhp(well_rates_bhp_limit, bhp_limit, deferred_logger);
                const double thp_limit = this->getTHPConstraint(deferred_logger);

                if (thp < thp_limit) {
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
    checkOperabilityUnderTHPLimitProducer(const Simulator& /* ebos_simulator */, Opm::DeferredLogger& deferred_logger)
    {
        const double obtain_bhp =  calculateBHPWithTHPTargetIPR(deferred_logger);

        if (obtain_bhp > 0.) {
            this->operability_status_.can_obtain_bhp_with_thp_limit = true;

            const double  bhp_limit = mostStrictBhpFromBhpLimits(deferred_logger);
            this->operability_status_.obey_bhp_limit_with_thp_limit = (obtain_bhp >= bhp_limit);

            const double thp_limit = this->getTHPConstraint(deferred_logger);
            if (obtain_bhp < thp_limit) {
                const std::string msg = " obtained bhp " + std::to_string(unit::convert::to(obtain_bhp, unit::barsa))
                                        + " bars is SMALLER than thp limit "
                                        + std::to_string(unit::convert::to(thp_limit, unit::barsa))
                                        + " bars as a producer for well " + name();
                deferred_logger.debug(msg);
            }
        } else {
            this->operability_status_.can_obtain_bhp_with_thp_limit = false;
            const double thp_limit = this->getTHPConstraint(deferred_logger);
            deferred_logger.debug(" COULD NOT find bhp value under thp_limit "
                          + std::to_string(unit::convert::to(thp_limit, unit::barsa))
                          + " bars for well " + name() + ", the well might need to be closed ");
            this->operability_status_.obey_bhp_limit_with_thp_limit = false;
        }
    }





    template<typename TypeTag>
    bool
    StandardWell<TypeTag>::
    allDrawDownWrongDirection(const Simulator& ebos_simulator) const
    {
        bool all_drawdown_wrong_direction = true;

        for (int perf = 0; perf < number_of_perforations_; ++perf) {
            const int cell_idx = well_cells_[perf];
            const auto& intQuants = *(ebos_simulator.model().cachedIntensiveQuantities(cell_idx, /*timeIdx=*/0));
            const auto& fs = intQuants.fluidState();

            const double pressure = (fs.pressure(FluidSystem::oilPhaseIdx)).value();
            const double bhp = getBhp().value();

            // Pressure drawdown (also used to determine direction of flow)
            const double well_pressure = bhp + perf_pressure_diffs_[perf];
            const double drawdown = pressure - well_pressure;

            // for now, if there is one perforation can produce/inject in the correct
            // direction, we consider this well can still produce/inject.
            // TODO: it can be more complicated than this to cause wrong-signed rates
            if ( (drawdown < 0. && well_type_ == INJECTOR) ||
                 (drawdown > 0. && well_type_ == PRODUCER) )  {
                all_drawdown_wrong_direction = false;
                break;
            }
        }

        return all_drawdown_wrong_direction;
    }




    template<typename TypeTag>
    bool
    StandardWell<TypeTag>::
    canProduceInjectWithCurrentBhp(const Simulator& ebos_simulator,
                                   const WellState& well_state,
                                   Opm::DeferredLogger& deferred_logger
                                   )
    {
        const double bhp = well_state.bhp()[index_of_well_];
        std::vector<double> well_rates;
        computeWellRatesWithBhp(ebos_simulator, bhp, well_rates, deferred_logger);

        const double sign = (well_type_ == PRODUCER) ? -1. : 1.;
        const double threshold = sign * std::numeric_limits<double>::min();

        bool can_produce_inject = false;
        for (const auto value : well_rates) {
            if (well_type_ == PRODUCER && value < threshold) {
                can_produce_inject = true;
                break;
            } else if (well_type_ == INJECTOR && value > threshold) {
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
        return !getAllowCrossFlow() && allDrawDownWrongDirection(ebos_simulator);
    }





    template<typename TypeTag>
    void
    StandardWell<TypeTag>::
    updateWellStateWithTHPTargetIPR(const Simulator& ebos_simulator,
                                    WellState& well_state,
                                    Opm::DeferredLogger& deferred_logger) const
    {
        if (well_type_ == PRODUCER) {
            updateWellStateWithTHPTargetIPRProducer(ebos_simulator,
                                                    well_state,
                                                    deferred_logger);
        }

        if (well_type_ == INJECTOR) {
            well_state.thp()[index_of_well_] = this->getTHPConstraint(deferred_logger);
            // TODO: more work needs to be done for the injectors here, while injectors
            // have been okay with the current strategy relying on well control equation directly.
        }
    }





    template<typename TypeTag>
    void
    StandardWell<TypeTag>::
    updateWellStateWithTHPTargetIPRProducer(const Simulator& ebos_simulator,
                                            WellState& well_state,
                                            Opm::DeferredLogger& deferred_logger) const
    {

        well_state.thp()[index_of_well_] = this->getTHPConstraint(deferred_logger);

        const double bhp = calculateBHPWithTHPTargetIPR(deferred_logger);

        assert(bhp > 0.0);

        well_state.bhp()[index_of_well_] = bhp;

        // TODO: explicit quantities are always tricky for this type of situation
        updatePrimaryVariables(well_state, deferred_logger);
        initPrimaryVariablesEvaluation();

        std::vector<double> rates;
        computeWellRatesWithBhp(ebos_simulator, bhp, rates, deferred_logger);

        // TODO: double checke the obtained rates
        // this is another places we might obtain negative rates

        for (int p = 0; p < number_of_phases_; ++p) {
            well_state.wellRates()[number_of_phases_ * index_of_well_ + p] = rates[p];
        }

        // TODO: there will be something need to be done for the cases not the defaulted 3 phases,
        // like 2 phases or solvent, polymer, etc. But we are not addressing them with THP control yet.
    }





    template<typename TypeTag>
    double
    StandardWell<TypeTag>::
    calculateBHPWithTHPTargetIPR(Opm::DeferredLogger& deferred_logger) const
    {
        const double thp_target = this->getTHPConstraint(deferred_logger);
        const double thp_control_index = this->getControlIndex(THP);
        const  int thp_table_id = well_controls_iget_vfp(well_controls_, thp_control_index);
        const double alq = well_controls_iget_alq(well_controls_, thp_control_index);

        // not considering injectors for now
        const double vfp_ref_depth = vfp_properties_->getProd()->getTable(thp_table_id)->getDatumDepth();

        // the density of the top perforation
        // TODO: make sure this is properly initialized
        // TODO: with IPR, it should be possible, even this well is a new well, we do not have
        // TODO: any information of previous rates. However, we are lacking the pressure though.
        // TODO: but let us not do more work related now to see what will happen
        const double rho = perf_densities_[0];

        // TODO: call a/the function for dp
        const double dp = (vfp_ref_depth - ref_depth_) * rho * gravity_;

        const double  bhp_limit = mostStrictBhpFromBhpLimits(deferred_logger);

        const double obtain_bhp = vfp_properties_->getProd()->calculateBhpWithTHPTarget(ipr_a_, ipr_b_,
                                             bhp_limit, thp_table_id, thp_target, alq, dp);

        return obtain_bhp;
    }





    template<typename TypeTag>
    void
    StandardWell<TypeTag>::
    computePropertiesForWellConnectionPressures(const Simulator& ebosSimulator,
                                                const WellState& well_state,
                                                std::vector<double>& b_perf,
                                                std::vector<double>& rsmax_perf,
                                                std::vector<double>& rvmax_perf,
                                                std::vector<double>& surf_dens_perf) const
    {
        const int nperf = number_of_perforations_;
        const PhaseUsage& pu = phaseUsage();
        b_perf.resize(nperf * num_components_);
        surf_dens_perf.resize(nperf * num_components_);
        const int w = index_of_well_;

        const bool waterPresent = FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx);
        const bool oilPresent = FluidSystem::phaseIsActive(FluidSystem::oilPhaseIdx);
        const bool gasPresent = FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx);

        //rs and rv are only used if both oil and gas is present
        if (oilPresent && gasPresent) {
            rsmax_perf.resize(nperf);
            rvmax_perf.resize(nperf);
        }

        // Compute the average pressure in each well block
        for (int perf = 0; perf < nperf; ++perf) {
            const int cell_idx = well_cells_[perf];
            const auto& intQuants = *(ebosSimulator.model().cachedIntensiveQuantities(cell_idx, /*timeIdx=*/0));
            const auto& fs = intQuants.fluidState();

            // TODO: this is another place to show why WellState need to be a vector of WellState.
            // TODO: to check why should be perf - 1
            const double p_above = perf == 0 ? well_state.bhp()[w] : well_state.perfPress()[first_perf_ + perf - 1];
            const double p_avg = (well_state.perfPress()[first_perf_ + perf] + p_above)/2;
            const double temperature = fs.temperature(FluidSystem::oilPhaseIdx).value();

            if (waterPresent) {
                const unsigned waterCompIdx = Indices::canonicalToActiveComponentIndex(FluidSystem::waterCompIdx);
                b_perf[ waterCompIdx + perf * num_components_] =
                FluidSystem::waterPvt().inverseFormationVolumeFactor(fs.pvtRegionIndex(), temperature, p_avg);
            }

            if (gasPresent) {
                const unsigned gasCompIdx = Indices::canonicalToActiveComponentIndex(FluidSystem::gasCompIdx);
                const int gaspos = gasCompIdx + perf * num_components_;
                const int gaspos_well = pu.phase_pos[Gas] + w * pu.num_phases;

                if (oilPresent) {
                    const int oilpos_well = pu.phase_pos[Oil] + w * pu.num_phases;
                    const double oilrate = std::abs(well_state.wellRates()[oilpos_well]); //in order to handle negative rates in producers
                    rvmax_perf[perf] = FluidSystem::gasPvt().saturatedOilVaporizationFactor(fs.pvtRegionIndex(), temperature, p_avg);
                    if (oilrate > 0) {
                        const double gasrate = std::abs(well_state.wellRates()[gaspos_well]) - well_state.solventWellRate(w);
                        double rv = 0.0;
                        if (gasrate > 0) {
                            rv = oilrate / gasrate;
                        }
                        rv = std::min(rv, rvmax_perf[perf]);

                        b_perf[gaspos] = FluidSystem::gasPvt().inverseFormationVolumeFactor(fs.pvtRegionIndex(), temperature, p_avg, rv);
                    }
                    else {
                        b_perf[gaspos] = FluidSystem::gasPvt().saturatedInverseFormationVolumeFactor(fs.pvtRegionIndex(), temperature, p_avg);
                    }

                } else {
                    b_perf[gaspos] = FluidSystem::gasPvt().saturatedInverseFormationVolumeFactor(fs.pvtRegionIndex(), temperature, p_avg);
                }
            }

            if (oilPresent) {
                const unsigned oilCompIdx = Indices::canonicalToActiveComponentIndex(FluidSystem::oilCompIdx);
                const int oilpos = oilCompIdx + perf * num_components_;
                const int oilpos_well = pu.phase_pos[Oil] + w * pu.num_phases;
                if (gasPresent) {
                    rsmax_perf[perf] = FluidSystem::oilPvt().saturatedGasDissolutionFactor(fs.pvtRegionIndex(), temperature, p_avg);
                    const int gaspos_well = pu.phase_pos[Gas] + w * pu.num_phases;
                    const double gasrate = std::abs(well_state.wellRates()[gaspos_well]) - well_state.solventWellRate(w);
                    if (gasrate > 0) {
                        const double oilrate = std::abs(well_state.wellRates()[oilpos_well]);
                        double rs = 0.0;
                        if (oilrate > 0) {
                            rs = gasrate / oilrate;
                        }
                        rs = std::min(rs, rsmax_perf[perf]);
                        b_perf[oilpos] = FluidSystem::oilPvt().inverseFormationVolumeFactor(fs.pvtRegionIndex(), temperature, p_avg, rs);
                    } else {
                        b_perf[oilpos] = FluidSystem::oilPvt().saturatedInverseFormationVolumeFactor(fs.pvtRegionIndex(), temperature, p_avg);
                    }
                } else {
                    b_perf[oilpos] = FluidSystem::oilPvt().saturatedInverseFormationVolumeFactor(fs.pvtRegionIndex(), temperature, p_avg);
                }
            }

            // Surface density.
            for (unsigned phaseIdx = 0; phaseIdx < FluidSystem::numPhases; ++phaseIdx) {
                if (!FluidSystem::phaseIsActive(phaseIdx)) {
                    continue;
                }

                const unsigned compIdx = Indices::canonicalToActiveComponentIndex(FluidSystem::solventComponentIndex(phaseIdx));
                surf_dens_perf[num_components_ * perf  + compIdx] = FluidSystem::referenceDensity( phaseIdx, fs.pvtRegionIndex() );
            }

            // We use cell values for solvent injector
            if (has_solvent) {
                b_perf[num_components_ * perf + contiSolventEqIdx] = intQuants.solventInverseFormationVolumeFactor().value();
                surf_dens_perf[num_components_ * perf + contiSolventEqIdx] = intQuants.solventRefDensity();
            }
        }
    }





    template<typename TypeTag>
    void
    StandardWell<TypeTag>::
    computeConnectionDensities(const std::vector<double>& perfComponentRates,
                               const std::vector<double>& b_perf,
                               const std::vector<double>& rsmax_perf,
                               const std::vector<double>& rvmax_perf,
                               const std::vector<double>& surf_dens_perf)
    {
        // Verify that we have consistent input.
        const int np = number_of_phases_;
        const int nperf = number_of_perforations_;
        const int num_comp = num_components_;

        // 1. Compute the flow (in surface volume units for each
        //    component) exiting up the wellbore from each perforation,
        //    taking into account flow from lower in the well, and
        //    in/out-flow at each perforation.
        std::vector<double> q_out_perf(nperf*num_comp);

        // TODO: investigate whether we should use the following techniques to calcuate the composition of flows in the wellbore
        // Iterate over well perforations from bottom to top.
        for (int perf = nperf - 1; perf >= 0; --perf) {
            for (int component = 0; component < num_comp; ++component) {
                if (perf == nperf - 1) {
                    // This is the bottom perforation. No flow from below.
                    q_out_perf[perf*num_comp+ component] = 0.0;
                } else {
                    // Set equal to flow from below.
                    q_out_perf[perf*num_comp + component] = q_out_perf[(perf+1)*num_comp + component];
                }
                // Subtract outflow through perforation.
                q_out_perf[perf*num_comp + component] -= perfComponentRates[perf*num_comp + component];
            }
        }

        // 2. Compute the component mix at each perforation as the
        //    absolute values of the surface rates divided by their sum.
        //    Then compute volume ratios (formation factors) for each perforation.
        //    Finally compute densities for the segments associated with each perforation.
        std::vector<double> mix(num_comp,0.0);
        std::vector<double> x(num_comp);
        std::vector<double> surf_dens(num_comp);

        for (int perf = 0; perf < nperf; ++perf) {
            // Find component mix.
            const double tot_surf_rate = std::accumulate(q_out_perf.begin() + num_comp*perf,
                                                         q_out_perf.begin() + num_comp*(perf+1), 0.0);
            if (tot_surf_rate != 0.0) {
                for (int component = 0; component < num_comp; ++component) {
                    mix[component] = std::fabs(q_out_perf[perf*num_comp + component]/tot_surf_rate);
                }
            } else {
                // No flow => use well specified fractions for mix.
                for (int component = 0; component < num_comp; ++component) {
                    if (component < np) {
                        mix[component] = comp_frac_[ ebosCompIdxToFlowCompIdx(component)];
                    }
                }
                // intialize 0.0 for comIdx >= np;
            }
            // Compute volume ratio.
            x = mix;

            // Subtract dissolved gas from oil phase and vapporized oil from gas phase
            if (FluidSystem::phaseIsActive(FluidSystem::gasCompIdx) && FluidSystem::phaseIsActive(FluidSystem::oilCompIdx)) {
                const unsigned gaspos = Indices::canonicalToActiveComponentIndex(FluidSystem::gasCompIdx);
                const unsigned oilpos = Indices::canonicalToActiveComponentIndex(FluidSystem::oilCompIdx);
                double rs = 0.0;
                double rv = 0.0;
                if (!rsmax_perf.empty() && mix[oilpos] > 0.0) {
                    rs = std::min(mix[gaspos]/mix[oilpos], rsmax_perf[perf]);
                }
                if (!rvmax_perf.empty() && mix[gaspos] > 0.0) {
                    rv = std::min(mix[oilpos]/mix[gaspos], rvmax_perf[perf]);
                }
                if (rs != 0.0) {
                    // Subtract gas in oil from gas mixture
                    x[gaspos] = (mix[gaspos] - mix[oilpos]*rs)/(1.0 - rs*rv);
                }
                if (rv != 0.0) {
                    // Subtract oil in gas from oil mixture
                    x[oilpos] = (mix[oilpos] - mix[gaspos]*rv)/(1.0 - rs*rv);;
                }
            }
            double volrat = 0.0;
            for (int component = 0; component < num_comp; ++component) {
                volrat += x[component] / b_perf[perf*num_comp+ component];
            }
            for (int component = 0; component < num_comp; ++component) {
                surf_dens[component] = surf_dens_perf[perf*num_comp+ component];
            }

            // Compute segment density.
            perf_densities_[perf] = std::inner_product(surf_dens.begin(), surf_dens.end(), mix.begin(), 0.0) / volrat;
        }
    }




    template<typename TypeTag>
    void
    StandardWell<TypeTag>::
    computeConnectionPressureDelta()
    {
        // Algorithm:

        // We'll assume the perforations are given in order from top to
        // bottom for each well.  By top and bottom we do not necessarily
        // mean in a geometric sense (depth), but in a topological sense:
        // the 'top' perforation is nearest to the surface topologically.
        // Our goal is to compute a pressure delta for each perforation.

        // 1. Compute pressure differences between perforations.
        //    dp_perf will contain the pressure difference between a
        //    perforation and the one above it, except for the first
        //    perforation for each well, for which it will be the
        //    difference to the reference (bhp) depth.

        const int nperf = number_of_perforations_;
        perf_pressure_diffs_.resize(nperf, 0.0);

        for (int perf = 0; perf < nperf; ++perf) {
            const double z_above = perf == 0 ? ref_depth_ : perf_depth_[perf - 1];
            const double dz = perf_depth_[perf] - z_above;
            perf_pressure_diffs_[perf] = dz * perf_densities_[perf] * gravity_;
        }

        // 2. Compute pressure differences to the reference point (bhp) by
        //    accumulating the already computed adjacent pressure
        //    differences, storing the result in dp_perf.
        //    This accumulation must be done per well.
        const auto beg = perf_pressure_diffs_.begin();
        const auto end = perf_pressure_diffs_.end();
        std::partial_sum(beg, end, beg);
    }





    template<typename TypeTag>
    ConvergenceReport
    StandardWell<TypeTag>::
    getWellConvergence(const std::vector<double>& B_avg,
                       Opm::DeferredLogger& deferred_logger) const
    {
        // the following implementation assume that the polymer is always after the w-o-g phases
        // For the polymer case and the energy case, there is one more mass balance equations of reservoir than wells
        assert((int(B_avg.size()) == num_components_) || has_polymer || has_energy);

        const double tol_wells = param_.tolerance_wells_;
        const double maxResidualAllowed = param_.max_residual_allowed_;

        std::vector<double> res(numWellEq);
        for (int eq_idx = 0; eq_idx < numWellEq; ++eq_idx) {
            // magnitude of the residual matters
            res[eq_idx] = std::abs(resWell_[0][eq_idx]);
        }

        std::vector<double> well_flux_residual(num_components_);

        // Finish computation
        for ( int compIdx = 0; compIdx < num_components_; ++compIdx )
        {
            well_flux_residual[compIdx] = B_avg[compIdx] * res[compIdx];
        }

        ConvergenceReport report;
        using CR = ConvergenceReport;
        CR::WellFailure::Type type = CR::WellFailure::Type::MassBalance;
        // checking if any NaN or too large residuals found
        for (unsigned phaseIdx = 0; phaseIdx < FluidSystem::numPhases; ++phaseIdx) {
            if (!FluidSystem::phaseIsActive(phaseIdx)) {
                continue;
            }

            const unsigned canonicalCompIdx = FluidSystem::solventComponentIndex(phaseIdx);
            const std::string& compName = FluidSystem::componentName(canonicalCompIdx);
            const int compIdx = Indices::canonicalToActiveComponentIndex(canonicalCompIdx);

            if (std::isnan(well_flux_residual[compIdx])) {
                report.setWellFailed({type, CR::Severity::NotANumber, compIdx, name()});
            } else if (well_flux_residual[compIdx] > maxResidualAllowed) {
                report.setWellFailed({type, CR::Severity::TooLarge, compIdx, name()});
            } else if (well_flux_residual[compIdx] > tol_wells) {
                report.setWellFailed({type, CR::Severity::Normal, compIdx, name()});
            }
        }

        checkConvergenceControlEq(report, deferred_logger);

        return report;
    }





    template<typename TypeTag>
    void
    StandardWell<TypeTag>::
    computeWellConnectionDensitesPressures(const WellState& well_state,
                                           const std::vector<double>& b_perf,
                                           const std::vector<double>& rsmax_perf,
                                           const std::vector<double>& rvmax_perf,
                                           const std::vector<double>& surf_dens_perf)
    {
        // Compute densities
        const int nperf = number_of_perforations_;
        const int np = number_of_phases_;
        std::vector<double> perfRates(b_perf.size(),0.0);

        for (int perf = 0; perf < nperf; ++perf) {
            for (int comp = 0; comp < np; ++comp) {
                perfRates[perf * num_components_ + comp] =  well_state.perfPhaseRates()[(first_perf_ + perf) * np + ebosCompIdxToFlowCompIdx(comp)];
            }
            if(has_solvent) {
                perfRates[perf * num_components_ + contiSolventEqIdx] =  well_state.perfRateSolvent()[first_perf_ + perf];
            }
        }

        computeConnectionDensities(perfRates, b_perf, rsmax_perf, rvmax_perf, surf_dens_perf);

        computeConnectionPressureDelta();

    }





    template<typename TypeTag>
    void
    StandardWell<TypeTag>::
    computeWellConnectionPressures(const Simulator& ebosSimulator,
                                   const WellState& well_state)
    {
         // 1. Compute properties required by computeConnectionPressureDelta().
         //    Note that some of the complexity of this part is due to the function
         //    taking std::vector<double> arguments, and not Eigen objects.
         std::vector<double> b_perf;
         std::vector<double> rsmax_perf;
         std::vector<double> rvmax_perf;
         std::vector<double> surf_dens_perf;
         computePropertiesForWellConnectionPressures(ebosSimulator, well_state, b_perf, rsmax_perf, rvmax_perf, surf_dens_perf);
         computeWellConnectionDensitesPressures(well_state, b_perf, rsmax_perf, rvmax_perf, surf_dens_perf);
    }





    template<typename TypeTag>
    void
    StandardWell<TypeTag>::
    solveEqAndUpdateWellState(WellState& well_state, Opm::DeferredLogger& deferred_logger)
    {
        if (!this->isOperable()) return;

        // We assemble the well equations, then we check the convergence,
        // which is why we do not put the assembleWellEq here.
        BVectorWell dx_well(1);
        invDuneD_.mv(resWell_, dx_well);

        updateWellState(dx_well, well_state, deferred_logger);
    }





    template<typename TypeTag>
    void
    StandardWell<TypeTag>::
    calculateExplicitQuantities(const Simulator& ebosSimulator,
                                const WellState& well_state,
                                Opm::DeferredLogger& deferred_logger)
    {
        updatePrimaryVariables(well_state, deferred_logger);
        initPrimaryVariablesEvaluation();
        computeWellConnectionPressures(ebosSimulator, well_state);
        computeAccumWell();
    }





    template<typename TypeTag>
    void
    StandardWell<TypeTag>::
    computeAccumWell()
    {
        for (int eq_idx = 0; eq_idx < numWellConservationEq; ++eq_idx) {
            F0_[eq_idx] = wellSurfaceVolumeFraction(eq_idx).value();
        }
    }





    template<typename TypeTag>
    void
    StandardWell<TypeTag>::
    apply(const BVector& x, BVector& Ax) const
    {
        if (!this->isOperable()) return;

        if ( param_.matrix_add_well_contributions_ )
        {
            // Contributions are already in the matrix itself
            return;
        }
        assert( Bx_.size() == duneB_.N() );
        assert( invDrw_.size() == invDuneD_.N() );

        // Bx_ = duneB_ * x
        duneB_.mv(x, Bx_);
        // invDBx = invDuneD_ * Bx_
        // TODO: with this, we modified the content of the invDrw_.
        // Is it necessary to do this to save some memory?
        BVectorWell& invDBx = invDrw_;
        invDuneD_.mv(Bx_, invDBx);

        // Ax = Ax - duneC_^T * invDBx
        duneC_.mmtv(invDBx,Ax);
    }




    template<typename TypeTag>
    void
    StandardWell<TypeTag>::
    apply(BVector& r) const
    {
        if (!this->isOperable()) return;

        assert( invDrw_.size() == invDuneD_.N() );

        // invDrw_ = invDuneD_ * resWell_
        invDuneD_.mv(resWell_, invDrw_);
        // r = r - duneC_^T * invDrw_
        duneC_.mmtv(invDrw_, r);
    }





    template<typename TypeTag>
    void
    StandardWell<TypeTag>::
    recoverSolutionWell(const BVector& x, BVectorWell& xw) const
    {
        if (!this->isOperable()) return;

        BVectorWell resWell = resWell_;
        // resWell = resWell - B * x
        duneB_.mmv(x, resWell);
        // xw = D^-1 * resWell
        invDuneD_.mv(resWell, xw);
    }





    template<typename TypeTag>
    void
    StandardWell<TypeTag>::
    recoverWellSolutionAndUpdateWellState(const BVector& x,
                                          WellState& well_state,
                                          Opm::DeferredLogger& deferred_logger) const
    {
        if (!this->isOperable()) return;

        BVectorWell xw(1);
        recoverSolutionWell(x, xw);
        updateWellState(xw, well_state, deferred_logger);
    }




    template<typename TypeTag>
    void
    StandardWell<TypeTag>::
    computeWellRatesWithBhp(const Simulator& ebosSimulator,
                            const double& bhp,
                            std::vector<double>& well_flux,
                            Opm::DeferredLogger& deferred_logger) const
    {

        const int np = number_of_phases_;
        well_flux.resize(np, 0.0);

        const bool allow_cf = getAllowCrossFlow();

        for (int perf = 0; perf < number_of_perforations_; ++perf) {
            const int cell_idx = well_cells_[perf];
            const auto& intQuants = *(ebosSimulator.model().cachedIntensiveQuantities(cell_idx, /*timeIdx=*/ 0));
            // flux for each perforation
            std::vector<EvalWell> mob(num_components_, 0.0);
            getMobility(ebosSimulator, perf, mob, deferred_logger);
            double trans_mult = ebosSimulator.problem().template rockCompTransMultiplier<double>(intQuants, cell_idx);
            const double Tw = well_index_[perf] * trans_mult;

            std::vector<EvalWell> cq_s(num_components_, 0.0);
            double perf_dis_gas_rate = 0.;
            double perf_vap_oil_rate = 0.;
            computePerfRate(intQuants, mob, EvalWell(bhp), Tw, perf, allow_cf,
                            cq_s, perf_dis_gas_rate, perf_vap_oil_rate, deferred_logger);

            for(int p = 0; p < np; ++p) {
                well_flux[ebosCompIdxToFlowCompIdx(p)] += cq_s[p].value();
            }
        }
    }



    template<typename TypeTag>
    void
    StandardWell<TypeTag>::
    computeWellRatesWithBhpPotential(const Simulator& ebosSimulator,
                            const std::vector<Scalar>& B_avg,
                            const double& bhp,
                            std::vector<double>& well_flux,
                            Opm::DeferredLogger& deferred_logger)
    {

        WellControls* wc = well_controls_;
        const int bhp_index = Base::getControlIndex(BHP);
        const double orig_bhp = well_controls_iget_target(wc, bhp_index);
        const auto orig_current = well_controls_get_current(wc);

        well_controls_iset_target(wc, bhp_index, bhp);
        well_controls_set_current(wc, bhp_index);

        // iterate to get a more accurate well density
        // create a copy of the well_state to use. If the operability checking is sucessful, we use this one
        // to replace the original one
        WellState well_state_copy = ebosSimulator.problem().wellModel().wellState();
        well_state_copy.currentControls()[index_of_well_] = bhp_index;

        bool converged = this->solveWellEqUntilConverged(ebosSimulator, B_avg, well_state_copy, deferred_logger);

        if (!converged) {
            const std::string msg = " well " + name() + " did not get converged during well potential calculations "
                                                        "returning zero values for the potential";
            deferred_logger.debug(msg);
            return;
        }
        updatePrimaryVariables(well_state_copy, deferred_logger);
        computeWellConnectionPressures(ebosSimulator, well_state_copy);
        initPrimaryVariablesEvaluation();


        computeWellRatesWithBhp(ebosSimulator, bhp, well_flux, deferred_logger);

        // reset bhp limit
        well_controls_iset_target(wc, bhp_index, orig_bhp);
        well_controls_set_current(wc, orig_current);
    }




    template<typename TypeTag>
    std::vector<double>
    StandardWell<TypeTag>::
    computeWellPotentialWithTHP(const Simulator& ebosSimulator,
                                const std::vector<Scalar>& B_avg,
                                const double initial_bhp, // bhp from BHP constraints
                                const std::vector<double>& initial_potential,
                                Opm::DeferredLogger& deferred_logger)
    {
        // TODO: pay attention to the situation that finally the potential is calculated based on the bhp control
        // TODO: should we consider the bhp constraints during the iterative process?
        const int np = number_of_phases_;

        assert( np == int(initial_potential.size()) );

        std::vector<double> potentials = initial_potential;
        std::vector<double> old_potentials = potentials; // keeping track of the old potentials

        double bhp = initial_bhp;
        double old_bhp = bhp;

        bool converged = false;
        const int max_iteration = 1000;
        const double bhp_tolerance = 1000.; // 1000 pascal

        int iteration = 0;

        while ( !converged && iteration < max_iteration ) {
            // for each iteration, we calculate the bhp based on the rates/potentials with thp constraints
            // with considering the bhp value from the bhp limits. At the beginning of each iteration,
            // we initialize the bhp to be the bhp value from the bhp limits. Then based on the bhp values calculated
            // from the thp constraints, we decide the effective bhp value for well potential calculation.
            bhp = initial_bhp;

            // The number of the well controls/constraints
            const int nwc = well_controls_get_num(well_controls_);

            for (int ctrl_index = 0; ctrl_index < nwc; ++ctrl_index) {
                if (well_controls_iget_type(well_controls_, ctrl_index) == THP) {
                    const Opm::PhaseUsage& pu = phaseUsage();

                    std::vector<double> rates(3, 0.0);
                    if (FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx)) {
                        rates[ Water ] = potentials[pu.phase_pos[ Water ] ];
                    }
                    if (FluidSystem::phaseIsActive(FluidSystem::oilPhaseIdx)) {
                        rates[ Oil ] = potentials[pu.phase_pos[ Oil ] ];
                    }
                    if (FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx)) {
                        rates[ Gas ] = potentials[pu.phase_pos[ Gas ] ];
                    }

                    const double bhp_calculated = calculateBhpFromThp(rates, ctrl_index, deferred_logger);

                    if (well_type_ == INJECTOR && bhp_calculated < bhp ) {
                        bhp = bhp_calculated;
                    }

                    if (well_type_ == PRODUCER && bhp_calculated > bhp) {
                        bhp = bhp_calculated;
                    }
                }
            }

            // there should be always some available bhp/thp constraints there
            if (std::isinf(bhp) || std::isnan(bhp)) {
                OPM_DEFLOG_THROW(std::runtime_error, "Unvalid bhp value obtained during the potential calculation for well " << name(), deferred_logger);
            }

            converged = std::abs(old_bhp - bhp) < bhp_tolerance;

            computeWellRatesWithBhpPotential(ebosSimulator, B_avg, bhp, potentials, deferred_logger);

            // checking whether the potentials have valid values
            for (const double value : potentials) {
                if (std::isinf(value) || std::isnan(value)) {
                    OPM_DEFLOG_THROW(std::runtime_error, "Unvalid potential value obtained during the potential calculation for well " << name(), deferred_logger);
                }
            }

            if (!converged) {
                old_bhp = bhp;
                for (int p = 0; p < np; ++p) {
                    // TODO: improve the interpolation, will it always be valid with the way below?
                    // TODO: finding better paramters, better iteration strategy for better convergence rate.
                    const double potential_update_damping_factor = 0.001;
                    potentials[p] = potential_update_damping_factor * potentials[p] + (1.0 - potential_update_damping_factor) * old_potentials[p];
                    old_potentials[p] = potentials[p];
                }
            }

            ++iteration;
        }

        if (!converged) {
            // TODO: temporary approach to recover the running of prediction case
            // more sophisicated treatment should be developed to handle the situation in a general way
            // OPM_DEFLOG_THROW(std::runtime_error, "Failed in getting converged for the potential calculation for well " << name(), deferred_logger);
            deferred_logger.warning("FAILURE_GETTING_CONVERGED_POTENTIAL",
                    "Failed in getting converged for the potential calculation for well " + name());
        }

        return potentials;
    }





    template<typename TypeTag>
    void
    StandardWell<TypeTag>::
    computeWellPotentials(const Simulator& ebosSimulator,
                          const std::vector<Scalar>& B_avg,
                          const WellState& well_state,
                          std::vector<double>& well_potentials,
                          Opm::DeferredLogger& deferred_logger) // const
    {

        updatePrimaryVariables(well_state, deferred_logger);
        computeWellConnectionPressures(ebosSimulator, well_state);

        // initialize the primary variables in Evaluation, which is used in computePerfRate for computeWellPotentials
        // TODO: for computeWellPotentials, no derivative is required actually
        initPrimaryVariablesEvaluation();

        const int np = number_of_phases_;
        well_potentials.resize(np, 0.0);

        // get the bhp value based on the bhp constraints
        const double bhp = mostStrictBhpFromBhpLimits(deferred_logger);

        // does the well have a THP related constraint?
        if ( !wellHasTHPConstraints() ) {
            assert(std::abs(bhp) != std::numeric_limits<double>::max());
            computeWellRatesWithBhpPotential(ebosSimulator, B_avg, bhp, well_potentials, deferred_logger);
        } else {
            // the well has a THP related constraint
            // checking whether a well is newly added, it only happens at the beginning of the report step
            if ( !well_state.effectiveEventsOccurred(index_of_well_) ) {
                for (int p = 0; p < np; ++p) {
                    // This is dangerous for new added well
                    // since we are not handling the initialization correctly for now
                    well_potentials[p] = well_state.wellRates()[index_of_well_ * np + p];
                }
            } else {
                // We need to generate a reasonable rates to start the iteration process
                computeWellRatesWithBhpPotential(ebosSimulator, B_avg, bhp, well_potentials, deferred_logger);
                for (double& value : well_potentials) {
                    // make the value a little safer in case the BHP limits are default ones
                    // TODO: a better way should be a better rescaling based on the investigation of the VFP table.
                    const double rate_safety_scaling_factor = 0.00001;
                    value *= rate_safety_scaling_factor;
                }
            }

            well_potentials = computeWellPotentialWithTHP(ebosSimulator, B_avg, bhp, well_potentials, deferred_logger);
        }
    }





    template<typename TypeTag>
    void
    StandardWell<TypeTag>::
    updatePrimaryVariables(const WellState& well_state, Opm::DeferredLogger& deferred_logger) const
    {

        if (!this->isOperable()) return;

        const int well_index = index_of_well_;
        const int np = number_of_phases_;

        // the weighted total well rate
        double total_well_rate = 0.0;
        for (int p = 0; p < np; ++p) {
            total_well_rate += scalingFactor(p) * well_state.wellRates()[np * well_index + p];
        }

        // Not: for the moment, the first primary variable for the injectors is not G_total. The injection rate
        // under surface condition is used here
        if (well_type_ == INJECTOR) {
            primary_variables_[WQTotal] = 0.;
            for (int p = 0; p < np; ++p) {
                // TODO: the use of comp_frac_ here is dangerous, since the injection phase can be different from
                // prefered phasse in WELSPECS, while comp_frac_ only reflect the one specified in WELSPECS
                primary_variables_[WQTotal] += well_state.wellRates()[np * well_index + p] * comp_frac_[p];
            }
        } else {
            for (int p = 0; p < np; ++p) {
                primary_variables_[WQTotal] = total_well_rate;
            }
        }


        const WellControls* wc = well_controls_;
        const double* distr = well_controls_get_current_distr(wc);
        const auto pu = phaseUsage();

        if(std::abs(total_well_rate) > 0.) {
            if (FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx)) {
                primary_variables_[WFrac] = scalingFactor(pu.phase_pos[Water]) * well_state.wellRates()[np*well_index + pu.phase_pos[Water]] / total_well_rate;
            }
            if (FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx)) {
                primary_variables_[GFrac] = scalingFactor(pu.phase_pos[Gas]) * (well_state.wellRates()[np*well_index + pu.phase_pos[Gas]] - well_state.solventWellRate(well_index)) / total_well_rate ;
            }
            if (has_solvent) {
                primary_variables_[SFrac] = scalingFactor(pu.phase_pos[Gas]) * well_state.solventWellRate(well_index) / total_well_rate ;
            }
        } else { // total_well_rate == 0
            if (well_type_ == INJECTOR) {
                // only single phase injection handled
                if (FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx)) {
                    if (distr[Water] > 0.0) {
                        primary_variables_[WFrac] = 1.0;
                    } else {
                        primary_variables_[WFrac] = 0.0;
                    }
                }

                if (FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx)) {
                    if (distr[pu.phase_pos[Gas]] > 0.0) {
                        primary_variables_[GFrac] = 1.0 - wsolvent();
                        if (has_solvent) {
                            primary_variables_[SFrac] = wsolvent();
                        }
                    } else {
                        primary_variables_[GFrac] = 0.0;
                    }
                }

                // TODO: it is possible to leave injector as a oil well,
                // when F_w and F_g both equals to zero, not sure under what kind of circumstance
                // this will happen.
            } else if (well_type_ == PRODUCER) { // producers
                // TODO: the following are not addressed for the solvent case yet
                if (FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx)) {
                    primary_variables_[WFrac] = 1.0 / np;
                }
                if (FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx)) {
                    primary_variables_[GFrac] = 1.0 / np;
                }
            } else {
                OPM_DEFLOG_THROW(std::logic_error, "Expected PRODUCER or INJECTOR type of well", deferred_logger);
            }
        }


        // BHP
        primary_variables_[Bhp] = well_state.bhp()[index_of_well_];
    }





    template<typename TypeTag>
    template<class ValueType>
    ValueType
    StandardWell<TypeTag>::
    calculateBhpFromThp(const std::vector<ValueType>& rates,
                        const int control_index,
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

        const int vfp        = well_controls_iget_vfp(well_controls_, control_index);
        const double& thp    = well_controls_iget_target(well_controls_, control_index);
        const double& alq    = well_controls_iget_alq(well_controls_, control_index);

        // pick the density in the top layer
        // TODO: it is possible it should be a Evaluation
        const double rho = perf_densities_[0];

        ValueType bhp = 0.;
        if (well_type_ == INJECTOR) {
            const double vfp_ref_depth = vfp_properties_->getInj()->getTable(vfp)->getDatumDepth();

            const double dp = wellhelpers::computeHydrostaticCorrection(ref_depth_, vfp_ref_depth, rho, gravity_);

            bhp = vfp_properties_->getInj()->bhp(vfp, aqua, liquid, vapour, thp) - dp;
         }
         else if (well_type_ == PRODUCER) {
             const double vfp_ref_depth = vfp_properties_->getProd()->getTable(vfp)->getDatumDepth();

             const double dp = wellhelpers::computeHydrostaticCorrection(ref_depth_, vfp_ref_depth, rho, gravity_);

             bhp = vfp_properties_->getProd()->bhp(vfp, aqua, liquid, vapour, thp, alq) - dp;
         }
         else {
             OPM_DEFLOG_THROW(std::logic_error, "Expected INJECTOR or PRODUCER well", deferred_logger);
         }

         return bhp;
    }





    template<typename TypeTag>
    double
    StandardWell<TypeTag>::
    calculateThpFromBhp(const std::vector<double>& rates,
                        const double bhp,
                        Opm::DeferredLogger& deferred_logger) const
    {
        assert(int(rates.size()) == 3); // the vfp related only supports three phases now.

        SummaryState summaryState;
        const double aqua = rates[Water];
        const double liquid = rates[Oil];
        const double vapour = rates[Gas];

        // pick the density in the top layer
        const double rho = perf_densities_[0];

        double thp = 0.0;
        if (well_type_ == INJECTOR) {
            const int table_id = well_ecl_.injectionControls(summaryState).vfp_table_number;
            const double vfp_ref_depth = vfp_properties_->getInj()->getTable(table_id)->getDatumDepth();
            const double dp = wellhelpers::computeHydrostaticCorrection(ref_depth_, vfp_ref_depth, rho, gravity_);

            thp = vfp_properties_->getInj()->thp(table_id, aqua, liquid, vapour, bhp + dp);
        }
        else if (well_type_ == PRODUCER) {
            const auto controls = well_ecl_.productionControls(summaryState);
            const int table_id = controls.vfp_table_number;
            const double alq = controls.alq_value;
            const double vfp_ref_depth = vfp_properties_->getProd()->getTable(table_id)->getDatumDepth();
            const double dp = wellhelpers::computeHydrostaticCorrection(ref_depth_, vfp_ref_depth, rho, gravity_);

            thp = vfp_properties_->getProd()->thp(table_id, aqua, liquid, vapour, bhp + dp, alq);
        }
        else {
            OPM_DEFLOG_THROW(std::logic_error, "Expected INJECTOR or PRODUCER well", deferred_logger);
        }

        return thp;
    }







    template<typename TypeTag>
    void
    StandardWell<TypeTag>::
    updateWaterMobilityWithPolymer(const Simulator& ebos_simulator,
                                   const int perf,
                                   std::vector<EvalWell>& mob,
                                   Opm::DeferredLogger& deferred_logger) const
    {
        const int cell_idx = well_cells_[perf];
        const auto& int_quant = *(ebos_simulator.model().cachedIntensiveQuantities(cell_idx, /*timeIdx=*/ 0));
        const EvalWell polymer_concentration = extendEval(int_quant.polymerConcentration());

        // TODO: not sure should based on the well type or injecting/producing peforations
        // it can be different for crossflow
        if (well_type_ == INJECTOR) {
            // assume fully mixing within injecting wellbore
            const auto& visc_mult_table = PolymerModule::plyviscViscosityMultiplierTable(int_quant.pvtRegionIndex());
            const unsigned waterCompIdx = Indices::canonicalToActiveComponentIndex(FluidSystem::waterCompIdx);
            mob[waterCompIdx] /= (extendEval(int_quant.waterViscosityCorrection()) * visc_mult_table.eval(polymer_concentration, /*extrapolate=*/true) );
        }

        if (PolymerModule::hasPlyshlog()) {
            // we do not calculate the shear effects for injection wells when they do not
            // inject polymer.
            if (well_type_ == INJECTOR && wpolymer() == 0.) {
                return;
            }
            // compute the well water velocity with out shear effects.
            // TODO: do we need to turn on crossflow here?
            const bool allow_cf = getAllowCrossFlow() || openCrossFlowAvoidSingularity(ebos_simulator);
            const EvalWell& bhp = getBhp();

            std::vector<EvalWell> cq_s(num_components_,0.0);
            double perf_dis_gas_rate = 0.;
            double perf_vap_oil_rate = 0.;
            double trans_mult = ebos_simulator.problem().template rockCompTransMultiplier<double>(int_quant, cell_idx);
            const double Tw = well_index_[perf] * trans_mult;
            computePerfRate(int_quant, mob, bhp, Tw, perf, allow_cf,
                            cq_s, perf_dis_gas_rate, perf_vap_oil_rate, deferred_logger);
            // TODO: make area a member
            const double area = 2 * M_PI * perf_rep_radius_[perf] * perf_length_[perf];
            const auto& material_law_manager = ebos_simulator.problem().materialLawManager();
            const auto& scaled_drainage_info =
                        material_law_manager->oilWaterScaledEpsInfoDrainage(cell_idx);
            const double swcr = scaled_drainage_info.Swcr;
            const EvalWell poro = extendEval(int_quant.porosity());
            const EvalWell sw = extendEval(int_quant.fluidState().saturation(FluidSystem::waterPhaseIdx));
            // guard against zero porosity and no water
            const EvalWell denom = Opm::max( (area * poro * (sw - swcr)), 1e-12);
            const unsigned waterCompIdx = Indices::canonicalToActiveComponentIndex(FluidSystem::waterCompIdx);
            EvalWell water_velocity = cq_s[waterCompIdx] / denom * extendEval(int_quant.fluidState().invB(FluidSystem::waterPhaseIdx));

            if (PolymerModule::hasShrate()) {
                // the equation for the water velocity conversion for the wells and reservoir are from different version
                // of implementation. It can be changed to be more consistent when possible.
                water_velocity *= PolymerModule::shrate( int_quant.pvtRegionIndex() ) / bore_diameters_[perf];
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
    StandardWell<TypeTag>::addWellContributions(Mat& mat) const
    {
        // We need to change matrx A as follows
        // A -= C^T D^-1 B
        // D is diagonal
        // B and C have 1 row, nc colums and nonzero
        // at (0,j) only if this well has a perforation at cell j.

        for ( auto colC = duneC_[0].begin(), endC = duneC_[0].end(); colC != endC; ++colC )
        {
            const auto row_index = colC.index();
            auto& row = mat[row_index];
            auto col = row.begin();

            for ( auto colB = duneB_[0].begin(), endB = duneB_[0].end(); colB != endB; ++colB )
            {
                const auto col_index = colB.index();
                // Move col to index col_index
                while ( col != row.end() && col.index() < col_index ) ++col;
                assert(col != row.end() && col.index() == col_index);

                Dune::FieldMatrix<Scalar, numWellEq, numEq> tmp;
                typename Mat::block_type tmp1;
                Dune::FMatrixHelp::multMatrix(invDuneD_[0][0],  (*colB), tmp);
                Detail::multMatrixTransposed((*colC), tmp, tmp1);
                (*col) -= tmp1;
            }
        }
    }





    template<typename TypeTag>
    double
    StandardWell<TypeTag>::
    relaxationFactorFraction(const double old_value,
                             const double dx)
    {
        assert(old_value >= 0. && old_value <= 1.0);

        double relaxation_factor = 1.;

        // updated values without relaxation factor
        const double possible_updated_value = old_value - dx;

        // 0.95 is an experimental value remains to be optimized
        if (possible_updated_value < 0.0) {
            relaxation_factor = std::abs(old_value / dx) * 0.95;
        } else if (possible_updated_value > 1.0) {
            relaxation_factor = std::abs((1. - old_value) / dx) * 0.95;
        }
        // if possible_updated_value is between 0. and 1.0, then relaxation_factor
        // remains to be one

        assert(relaxation_factor >= 0. && relaxation_factor <= 1.);

        return relaxation_factor;
    }





    template<typename TypeTag>
    double
    StandardWell<TypeTag>::
    relaxationFactorFractionsProducer(const std::vector<double>& primary_variables,
                                      const BVectorWell& dwells)
    {
        // TODO: not considering solvent yet
        // 0.95 is a experimental value, which remains to be optimized
        double relaxation_factor = 1.0;

        if (FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx)) {
            const double relaxation_factor_w = relaxationFactorFraction(primary_variables[WFrac], dwells[0][WFrac]);
            relaxation_factor = std::min(relaxation_factor, relaxation_factor_w);
        }

        if (FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx)) {
            const double relaxation_factor_g = relaxationFactorFraction(primary_variables[GFrac], dwells[0][GFrac]);
            relaxation_factor = std::min(relaxation_factor, relaxation_factor_g);
        }

        if (FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx) && FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx)) {
            // We need to make sure the even with the relaxation_factor, the sum of F_w and F_g is below one, so there will
            // not be negative oil fraction later
            const double original_sum = primary_variables[WFrac] + primary_variables[GFrac];
            const double relaxed_update = (dwells[0][WFrac] + dwells[0][GFrac]) * relaxation_factor;
            const double possible_updated_sum = original_sum - relaxed_update;

            if (possible_updated_sum > 1.0) {
                assert(relaxed_update != 0.);

                const double further_relaxation_factor = std::abs((1. - original_sum) / relaxed_update) * 0.95;
                relaxation_factor *= further_relaxation_factor;
            }
        }

        assert(relaxation_factor >= 0.0 && relaxation_factor <= 1.0);

        return relaxation_factor;
    }





    template<typename TypeTag>
    double
    StandardWell<TypeTag>::
    relaxationFactorRate(const std::vector<double>& primary_variables,
                         const BVectorWell& dwells)
    {
        double relaxation_factor = 1.0;

        // For injector, we only check the total rates to avoid sign change of rates
        const double original_total_rate = primary_variables[WQTotal];
        const double newton_update = dwells[0][WQTotal];
        const double possible_update_total_rate = primary_variables[WQTotal] - newton_update;

        // 0.8 here is a experimental value, which remains to be optimized
        // if the original rate is zero or possible_update_total_rate is zero, relaxation_factor will
        // always be 1.0, more thoughts might be needed.
        if (original_total_rate * possible_update_total_rate < 0.) { // sign changed
            relaxation_factor = std::abs(original_total_rate / newton_update) * 0.8;
        }

        assert(relaxation_factor >= 0.0 && relaxation_factor <= 1.0);

        return relaxation_factor;
    }





    template<typename TypeTag>
    void
    StandardWell<TypeTag>::
    wellTestingPhysical(const Simulator& ebos_simulator, const std::vector<double>& B_avg,
                        const double /* simulation_time */, const int /* report_step */,
                        WellState& well_state, WellTestState& welltest_state,
                        Opm::DeferredLogger& deferred_logger)
    {
        deferred_logger.debug(" well " + name() + " is being tested for physical limits");

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
        updatePrimaryVariables(well_state_copy, deferred_logger);
        initPrimaryVariablesEvaluation();

        const bool converged = this->solveWellEqUntilConverged(ebos_simulator, B_avg, well_state_copy, deferred_logger);

        if (!converged) {
            const std::string msg = " well " + name() + " did not get converged during well testing for physical reason";
            deferred_logger.debug(msg);
            return;
        }

        if (this->isOperable() ) {
            welltest_state.openWell(name() );
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
    StandardWell<TypeTag>::
    updateWaterThroughput(const double dt OPM_UNUSED, WellState& well_state OPM_UNUSED) const
    {
    }





    template<typename TypeTag>
    void
    StandardWell<TypeTag>::
    checkConvergenceControlEq(ConvergenceReport& report,
                              DeferredLogger& deferred_logger) const
    {
        double control_tolerance = 0.;
        using CR = ConvergenceReport;
        CR::WellFailure::Type ctrltype = CR::WellFailure::Type::Invalid;
        switch(well_controls_get_current_type(well_controls_)) {
            case THP:
                ctrltype = CR::WellFailure::Type::ControlTHP;
                control_tolerance = 1.e4; // 0.1 bar
                break;
            case BHP:  // pressure type of control
                ctrltype = CR::WellFailure::Type::ControlBHP;
                control_tolerance = 1.e3; // 0.01 bar
                break;
            case RESERVOIR_RATE:
            case SURFACE_RATE:
                ctrltype = CR::WellFailure::Type::ControlRate;
                control_tolerance = 1.e-4; // smaller tolerance for rate control
                break;
            default:
                OPM_DEFLOG_THROW(std::runtime_error, "Unknown well control control types for well " << name(), deferred_logger);
        }
        assert(ctrltype != CR::WellFailure::Type::Invalid);

        const double well_control_residual = std::abs(resWell_[0][numWellEq -1]);
        const int dummy_component = -1;
        const double max_residual_allowed = param_.max_residual_allowed_;
        if (std::isnan(well_control_residual)) {
            report.setWellFailed({ctrltype, CR::Severity::NotANumber, dummy_component, name()});
        } else if (well_control_residual > max_residual_allowed * 10.) {
            report.setWellFailed({ctrltype, CR::Severity::TooLarge, dummy_component, name()});
        } else if ( well_control_residual > control_tolerance) {
            report.setWellFailed({ctrltype, CR::Severity::Normal, dummy_component, name()});
        }
    }
}
