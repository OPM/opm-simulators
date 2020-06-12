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

#include <opm/common/utility/numeric/RootFinders.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/Well/WellInjectionProperties.hpp>
#include <opm/simulators/utils/DeferredLoggingErrorHelpers.hpp>
#include <opm/simulators/linalg/MatrixBlock.hpp>

namespace Opm
{

    template<typename TypeTag>
    StandardWell<TypeTag>::
    StandardWell(const Well& well, const int time_step,
                 const ModelParameters& param,
                 const RateConverterType& rate_converter,
                 const int pvtRegionIdx,
                 const int num_components,
                 const int num_phases,
                 const int index_of_well,
                 const int first_perf_index,
                 const std::vector<PerforationData>& perf_data)
        : Base(well, time_step, param, rate_converter, pvtRegionIdx, num_components, num_phases, index_of_well, first_perf_index, perf_data)
    , perf_densities_(number_of_perforations_)
    , perf_pressure_diffs_(number_of_perforations_)
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

        // counting/updating primary variable numbers
        if (this->has_polymermw && this->isInjector()) {
            // adding a primary variable for water perforation rate per connection
            numWellEq_ += number_of_perforations_;
            // adding a primary variable for skin pressure per connection
            numWellEq_ += number_of_perforations_;
        }

        // with the updated numWellEq_, we can initialize the primary variables and matrices now
        primary_variables_.resize(numWellEq_, 0.0);
        primary_variables_evaluation_.resize(numWellEq_, EvalWell{numWellEq_ + numEq, 0.0});

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
        // the block size is run-time determined now
        invDuneD_[0][0].resize(numWellEq_, numWellEq_);

        for (auto row = duneB_.createbegin(), end = duneB_.createend(); row!=end; ++row) {
            for (int perf = 0 ; perf < number_of_perforations_; ++perf) {
                const int cell_idx = well_cells_[perf];
                row.insert(cell_idx);
            }
        }

        for (int perf = 0 ; perf < number_of_perforations_; ++perf) {
            const int cell_idx = well_cells_[perf];
             // the block size is run-time determined now
             duneB_[0][cell_idx].resize(numWellEq_, numEq);
        }

        // make the C^T matrix
        for (auto row = duneC_.createbegin(), end = duneC_.createend(); row!=end; ++row) {
            for (int perf = 0; perf < number_of_perforations_; ++perf) {
                const int cell_idx = well_cells_[perf];
                row.insert(cell_idx);
            }
        }

        for (int perf = 0; perf < number_of_perforations_; ++perf) {
            const int cell_idx = well_cells_[perf];
            duneC_[0][cell_idx].resize(numWellEq_, numEq);
        }

        resWell_.resize(1);
        // the block size of resWell_ is also run-time determined now
        resWell_[0].resize(numWellEq_);

        // resize temporary class variables
        Bx_.resize( duneB_.N() );
        for (unsigned i = 0; i < duneB_.N(); ++i) {
            Bx_[i].resize(numWellEq_);
        }

        invDrw_.resize( invDuneD_.N() );
        for (unsigned i = 0; i < invDuneD_.N(); ++i) {
            invDrw_[i].resize(numWellEq_);
        }
    }





    template<typename TypeTag>
    void StandardWell<TypeTag>::
    initPrimaryVariablesEvaluation() const
    {
        for (int eqIdx = 0; eqIdx < numWellEq_; ++eqIdx) {
            primary_variables_evaluation_[eqIdx] =
                EvalWell::createVariable(numWellEq_ + numEq, primary_variables_[eqIdx], numEq + eqIdx);
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

        if (this->isInjector()) { // only single phase injection
            double inj_frac = 0.0;
            switch (this->wellEcl().injectorType()) {
            case InjectorType::WATER:
                if (comp_idx == int(Indices::canonicalToActiveComponentIndex(FluidSystem::waterCompIdx))) {
                    inj_frac = 1.0;
                }
                break;
            case InjectorType::GAS:
                if (has_solvent && comp_idx == contiSolventEqIdx) { // solvent
                    inj_frac = wsolvent();
                } else if (comp_idx == int(Indices::canonicalToActiveComponentIndex(FluidSystem::gasCompIdx))) {
                    inj_frac = has_solvent ? 1.0 - wsolvent() : 1.0;
                }
                break;
            case InjectorType::OIL:
                if (comp_idx == int(Indices::canonicalToActiveComponentIndex(FluidSystem::oilCompIdx))) {
                    inj_frac = 1.0;
                }
                break;
            case InjectorType::MULTI:
                // Not supported.
                // deferred_logger.warning("MULTI_PHASE_INJECTOR_NOT_SUPPORTED",
                //                         "Multi phase injectors are not supported, requested for well " + name());
                break;
            }
            return inj_frac * primary_variables_evaluation_[WQTotal];
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
        if (FluidSystem::numActivePhases() == 1) {
            return EvalWell(numWellEq_ + numEq, 1.0);
        }

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
        EvalWell well_fraction(numWellEq_ + numEq, 1.0);
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

        EvalWell sum_volume_fraction_scaled(numWellEq_ + numEq, 0.);
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
        EvalWell out(numWellEq_ + numEq, in.value());
        for(int eqIdx = 0; eqIdx < numEq;++eqIdx) {
            out.setDerivative(eqIdx, in.derivative(eqIdx));
        }
        return out;
    }





    template<typename TypeTag>
    typename StandardWell<TypeTag>::Eval
    StandardWell<TypeTag>::getPerfCellPressure(const typename StandardWell<TypeTag>::FluidState& fs) const
    {
        Eval pressure;
        if (Indices::oilEnabled) {
            pressure = fs.pressure(FluidSystem::oilPhaseIdx);
        } else {
            if (Indices::waterEnabled) {
                pressure = fs.pressure(FluidSystem::waterPhaseIdx);
            } else {
                pressure = fs.pressure(FluidSystem::gasPhaseIdx);
            }
        }
        return pressure;
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
        const EvalWell pressure = extendEval(getPerfCellPressure(fs));
        const EvalWell rs = extendEval(fs.Rs());
        const EvalWell rv = extendEval(fs.Rv());
        std::vector<EvalWell> b_perfcells_dense(num_components_, EvalWell{numWellEq_ + numEq, 0.0});
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
        EvalWell drawdown = pressure - well_pressure;

        if (this->has_polymermw && this->isInjector()) {
            const int pskin_index = Bhp + 1 + number_of_perforations_ + perf;
            const EvalWell& skin_pressure = primary_variables_evaluation_[pskin_index];
            drawdown += skin_pressure;
        }

        // producing perforations
        if ( drawdown.value() > 0 )  {
            //Do nothing if crossflow is not allowed
            if (!allow_cf && this->isInjector()) {
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
                if (this->isProducer()) {
                    perf_dis_gas_rate = dis_gas.value();
                    perf_vap_oil_rate = vap_oil.value();
                }
            }

        } else {
            //Do nothing if crossflow is not allowed
            if (!allow_cf && this->isProducer()) {
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
            std::vector<EvalWell> cmix_s(num_components_, EvalWell{numWellEq_ + numEq});
            for (int componentIdx = 0; componentIdx < num_components_; ++componentIdx) {
                cmix_s[componentIdx] = wellSurfaceVolumeFraction(componentIdx);
            }

            // compute volume ratio between connection at standard conditions
            EvalWell volumeRatio(numWellEq_ + numEq, 0.);
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
                const EvalWell d = EvalWell(numWellEq_ + numEq, 1.0) - rv * rs;

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
            if (this->isProducer()) {
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
                   const std::vector<Scalar>& B_avg,
                   const double dt,
                   WellState& well_state,
                   Opm::DeferredLogger& deferred_logger)
    {
        const bool use_inner_iterations = param_.use_inner_iterations_wells_;
        if (use_inner_iterations) {
            this->iterateWellEquations(ebosSimulator, B_avg, dt, well_state, deferred_logger);
        }

        // TODO: inj_controls and prod_controls are not used in the following function for now
        const auto& summary_state = ebosSimulator.vanguard().summaryState();
        const auto inj_controls = well_ecl_.isInjector() ? well_ecl_.injectionControls(summary_state) : Well::InjectionControls(0);
        const auto prod_controls = well_ecl_.isProducer() ? well_ecl_.productionControls(summary_state) : Well::ProductionControls(0);
        assembleWellEqWithoutIteration(ebosSimulator, dt, inj_controls, prod_controls, well_state, deferred_logger);
    }





    template<typename TypeTag>
    void
    StandardWell<TypeTag>::
    assembleWellEqWithoutIteration(const Simulator& ebosSimulator,
                                   const double dt,
                                   const Well::InjectionControls& /*inj_controls*/,
                                   const Well::ProductionControls& /*prod_controls*/,
                                   WellState& well_state,
                                   Opm::DeferredLogger& deferred_logger)
    {
        // TODO: only_wells should be put back to save some computation
        // for example, the matrices B C does not need to update if only_wells

        checkWellOperability(ebosSimulator, well_state, deferred_logger);

        if (!this->isOperable() && !this->wellIsStopped()) return;

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
            std::vector<EvalWell> mob(num_components_, {numWellEq_ + numEq, 0.});
            getMobility(ebosSimulator, perf, mob, deferred_logger);

            std::vector<EvalWell> cq_s(num_components_, {numWellEq_ + numEq, 0.});
            double perf_dis_gas_rate = 0.;
            double perf_vap_oil_rate = 0.;
            double trans_mult = ebosSimulator.problem().template rockCompTransMultiplier<double>(intQuants,  cell_idx);
            const double Tw = well_index_[perf] * trans_mult;
            computePerfRate(intQuants, mob, bhp, Tw, perf, allow_cf,
                            cq_s, perf_dis_gas_rate, perf_vap_oil_rate, deferred_logger);

            // better way to do here is that use the cq_s and then replace the cq_s_water here?
            if (has_polymer && this->has_polymermw && this->isInjector()) {
                handleInjectivityRateAndEquations(intQuants, well_state, perf, cq_s, deferred_logger);
            }

            // updating the solution gas rate and solution oil rate
            if (this->isProducer()) {
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
                for (int pvIdx = 0; pvIdx < numWellEq_; ++pvIdx) {
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
                    EvalWell cq_r_thermal(numWellEq_ + numEq, 0.);
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
                    }
                    // compute the thermal flux
                    cq_r_thermal *= extendEval(fs.enthalpy(phaseIdx)) * extendEval(fs.density(phaseIdx));
                    connectionRates_[perf][contiEnergyEqIdx] += Base::restrictEval(cq_r_thermal);
                }
            }

            if (has_polymer) {
                // TODO: the application of well efficiency factor has not been tested with an example yet
                const unsigned waterCompIdx = Indices::canonicalToActiveComponentIndex(FluidSystem::waterCompIdx);
                EvalWell cq_s_poly = cq_s[waterCompIdx];
                if (this->isInjector()) {
                    cq_s_poly *= wpolymer();
                } else {
                    cq_s_poly *= extendEval(intQuants.polymerConcentration() * intQuants.polymerViscosityCorrection());
                }
                // Note. Efficiency factor is handled in the output layer
                well_state.perfRatePolymer()[first_perf_ + perf] = cq_s_poly.value();

                cq_s_poly *= well_efficiency_factor_;
                connectionRates_[perf][contiPolymerEqIdx] = Base::restrictEval(cq_s_poly);

                if (this->has_polymermw) {
                    updateConnectionRatePolyMW(cq_s_poly, intQuants, well_state, perf, deferred_logger);
                }
            }

            if (has_foam) {
                // TODO: the application of well efficiency factor has not been tested with an example yet
                const unsigned gasCompIdx = Indices::canonicalToActiveComponentIndex(FluidSystem::gasCompIdx);
                EvalWell cq_s_foam = cq_s[gasCompIdx] * well_efficiency_factor_;
                if (this->isInjector()) {
                    cq_s_foam *= wfoam();
                } else {
                    cq_s_foam *= extendEval(intQuants.foamConcentration());
                }
                connectionRates_[perf][contiFoamEqIdx] = Base::restrictEval(cq_s_foam);
            }

            if (has_brine) {
                // TODO: the application of well efficiency factor has not been tested with an example yet
                const unsigned waterCompIdx = Indices::canonicalToActiveComponentIndex(FluidSystem::waterCompIdx);
                EvalWell cq_s_sm = cq_s[waterCompIdx];
                if (this->isInjector()) {
                    cq_s_sm *= wsalt();
                } else {
                    cq_s_sm *= extendEval(intQuants.fluidState().saltConcentration());
                }
                // Note. Efficiency factor is handled in the output layer
                well_state.perfRateBrine()[first_perf_ + perf] = cq_s_sm.value();

                cq_s_sm *= well_efficiency_factor_;
                connectionRates_[perf][contiBrineEqIdx] = Base::restrictEval(cq_s_sm);
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
                    const auto& fs = intQuants.fluidState();
                    Eval perf_pressure = getPerfCellPressure(fs);
                    const double drawdown  = well_state.perfPress()[first_perf_ + perf] - perf_pressure.value();
                    const bool new_well = schedule.hasWellGroupEvent(name(), ScheduleEvents::NEW_WELL, current_step_);
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
            EvalWell resWell_loc(numWellEq_ + numEq, 0.0);
            if (FluidSystem::numActivePhases() > 1) {
                assert(dt > 0);
                resWell_loc += (wellSurfaceVolumeFraction(componentIdx) - F0_[componentIdx]) * volume / dt;
            }
            resWell_loc -= getQs(componentIdx) * well_efficiency_factor_;
            for (int pvIdx = 0; pvIdx < numWellEq_; ++pvIdx) {
                invDuneD_[0][0][componentIdx][pvIdx] += resWell_loc.derivative(pvIdx+numEq);
            }
            resWell_[0][componentIdx] += resWell_loc.value();
        }

        const auto& summaryState = ebosSimulator.vanguard().summaryState();
        const Opm::Schedule& schedule = ebosSimulator.vanguard().schedule();

        if(ebosSimulator.model().linearizer().getLinearizationType().type == Opm::LinearizationType::seqtransport) {
            assembleControlEqSeqTrans(well_state);
        } else {
            assembleControlEq(well_state, schedule, summaryState, deferred_logger);
        }

        // do the local inversion of D.
        try {
            Dune::ISTLUtility::invertMatrix(invDuneD_[0][0]);
        } catch( ... ) {
            OPM_DEFLOG_THROW(Opm::NumericalIssue,"Error when inverting local well equations for well " + name(), deferred_logger);
        }


    }





    template <typename TypeTag>
    void
    StandardWell<TypeTag>::assembleControlEqSeqTrans(const WellState& /* well_state */)
    {
        // TODO: we can get the current bhp value from well_state
        // If the primary variables are well synchronized, they should have the same value.
        EvalWell control_eq(numWellEq_ + numEq, 0.0);
        const EvalWell& bhp = getBhp();
        control_eq = bhp - bhp.value();
        resWell_[0][Bhp] = control_eq.value();
        for (int pv_idx = 0; pv_idx < numWellEq_; ++pv_idx) {
            invDuneD_[0][0][Bhp][pv_idx] = control_eq.derivative(pv_idx + numEq);
        }
    }





    template <typename TypeTag>
    void
    StandardWell<TypeTag>::assembleControlEq(const WellState& well_state,
                                             const Opm::Schedule& schedule,
                                             const SummaryState& summaryState,
                                             Opm::DeferredLogger& deferred_logger)
    {
        EvalWell control_eq(numWellEq_ + numEq, 0.0);

        const auto& well = well_ecl_;

        auto getRates = [&]() {
            std::vector<EvalWell> rates(3, EvalWell(numWellEq_ + numEq, 0.0));
            if (FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx)) {
                rates[Water] = getQs(Indices::canonicalToActiveComponentIndex(FluidSystem::waterCompIdx));
            }
            if (FluidSystem::phaseIsActive(FluidSystem::oilPhaseIdx)) {
                rates[Oil] = getQs(Indices::canonicalToActiveComponentIndex(FluidSystem::oilCompIdx));
            }
            if (FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx)) {
                rates[Gas] = getQs(Indices::canonicalToActiveComponentIndex(FluidSystem::gasCompIdx));
            }
            return rates;
        };

        if (wellIsStopped_) {
            control_eq = getWQTotal();
        } else if (this->isInjector()) {
            // Find injection rate.
            const EvalWell injection_rate = getWQTotal();
            // Setup function for evaluation of BHP from THP (used only if needed).
            auto bhp_from_thp = [&]() {
                const auto rates = getRates();
                return calculateBhpFromThp(rates, well, summaryState, deferred_logger);
            };
            // Call generic implementation.
            const auto& inj_controls = well.injectionControls(summaryState);
            Base::assembleControlEqInj(well_state, schedule, summaryState, inj_controls, getBhp(), injection_rate, bhp_from_thp, control_eq, deferred_logger);
        } else {
            // Find rates.
            const auto rates = getRates();
            // Setup function for evaluation of BHP from THP (used only if needed).
            auto bhp_from_thp = [&]() {
                return calculateBhpFromThp(rates, well, summaryState, deferred_logger);
            };
            // Call generic implementation.
            const auto& prod_controls = well.productionControls(summaryState);
            Base::assembleControlEqProd(well_state, schedule, summaryState, prod_controls, getBhp(), rates, bhp_from_thp, control_eq, deferred_logger);
        }

        // using control_eq to update the matrix and residuals
        // TODO: we should use a different index system for the well equations
        resWell_[0][Bhp] = control_eq.value();
        for (int pv_idx = 0; pv_idx < numWellEq_; ++pv_idx) {
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
        if (!this->isOperable() && !this->wellIsStopped()) return;

        updatePrimaryVariablesNewton(dwells, well_state);

        updateWellStateFromPrimaryVariables(well_state, deferred_logger);
        Base::calculateReservoirRates(well_state);
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
        const double relaxation_factor_fractions = (this->isProducer()) ?
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

        updateExtraPrimaryVariables(dwells);

#ifndef NDEBUG
        for (double v : primary_variables_) {
            assert(Opm::isfinite(v));
        }
#endif

    }





    template<typename TypeTag>
    void
    StandardWell<TypeTag>::
    updateExtraPrimaryVariables(const BVectorWell& dwells) const
    {
        // for the water velocity and skin pressure
        if (this->has_polymermw && this->isInjector()) {
            for (int perf = 0; perf < number_of_perforations_; ++perf) {
                const int wat_vel_index = Bhp + 1 + perf;
                const int pskin_index = Bhp + 1 + number_of_perforations_ + perf;

                const double relaxation_factor = 0.9;
                const double dx_wat_vel = dwells[0][wat_vel_index];
                primary_variables_[wat_vel_index] -= relaxation_factor * dx_wat_vel;

                const double dx_pskin = dwells[0][pskin_index];
                primary_variables_[pskin_index] -= relaxation_factor * dx_pskin;
            }
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
        if (this->isProducer()) {
            const double g_total = primary_variables_[WQTotal];
            for (int p = 0; p < number_of_phases_; ++p) {
                well_state.wellRates()[index_of_well_ * number_of_phases_ + p] = g_total * F[p];
            }
        } else { // injectors
            for (int p = 0; p < number_of_phases_; ++p) {
                well_state.wellRates()[index_of_well_ * number_of_phases_ + p] = 0.0;
            }
            switch (this->wellEcl().injectorType()) {
            case InjectorType::WATER:
                well_state.wellRates()[index_of_well_ * number_of_phases_ + pu.phase_pos[Water]] = primary_variables_[WQTotal];
                break;
            case InjectorType::GAS:
                well_state.wellRates()[index_of_well_ * number_of_phases_ + pu.phase_pos[Gas]] = primary_variables_[WQTotal];
                break;
            case InjectorType::OIL:
                well_state.wellRates()[index_of_well_ * number_of_phases_ + pu.phase_pos[Oil]] = primary_variables_[WQTotal];
                break;
            case InjectorType::MULTI:
                // Not supported.
                deferred_logger.warning("MULTI_PHASE_INJECTOR_NOT_SUPPORTED",
                                        "Multi phase injectors are not supported, requested for well " + name());
                break;
            }
        }

        updateThp(well_state, deferred_logger);

        // other primary variables related to polymer injectivity study
        if (this->has_polymermw && this->isInjector()) {
            for (int perf = 0; perf < number_of_perforations_; ++perf) {
                well_state.perfWaterVelocity()[first_perf_ + perf] = primary_variables_[Bhp + 1 + perf];
                well_state.perfSkinPressure()[first_perf_ + perf] = primary_variables_[Bhp + 1 + number_of_perforations_ + perf];
            }
        }
    }





    template<typename TypeTag>
    void
    StandardWell<TypeTag>::
    updateThp(WellState& well_state, Opm::DeferredLogger& deferred_logger) const
    {
        // When there is no vaild VFP table provided, we set the thp to be zero.
        if (!this->isVFPActive(deferred_logger) || this->wellIsStopped()) {
            well_state.thp()[index_of_well_] = 0.;
            return;
        }

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





    template<typename TypeTag>
    void
    StandardWell<TypeTag>::
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

        if (wellIsStopped_) {
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
                throw("Expected WATER, OIL or GAS as type for injectors " + well.name());
            }

            const Opm::Well::InjectorCMode& current = well_state.currentInjectionControls()[well_index];

            switch(current) {
            case Well::InjectorCMode::RATE:
            {
                well_state.wellRates()[well_index*np + phasePos] = controls.surface_rate;
                break;
            }

            case Well::InjectorCMode::RESV:
            {
                std::vector<double> convert_coeff(number_of_phases_, 1.0);
                Base::rateConverter_.calcCoeff(/*fipreg*/ 0, Base::pvtRegionIdx_, convert_coeff);
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
                double bhp = calculateBhpFromThp(rates, well, summaryState, deferred_logger);
                well_state.bhp()[well_index] = bhp;
                break;
            }
            case Well::InjectorCMode::BHP:
            {
                well_state.bhp()[well_index] = controls.bhp_limit;
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
            const Well::ProducerCMode& current = well_state.currentProductionControls()[well_index];
            const auto& controls = well.productionControls(summaryState);

            switch (current) {
            case Well::ProducerCMode::ORAT:
            {
                double current_rate = -well_state.wellRates()[ well_index*np + pu.phase_pos[Oil] ];

                if (current_rate == 0.0)
                    break;

                for (int p = 0; p<np; ++p) {
                    well_state.wellRates()[well_index*np + p] *= controls.oil_rate/current_rate;
                }
                break;
            }
            case Well::ProducerCMode::WRAT:
            {
                double current_rate = -well_state.wellRates()[ well_index*np + pu.phase_pos[Water] ];

                if (current_rate == 0.0)
                    break;

                for (int p = 0; p<np; ++p) {
                    well_state.wellRates()[well_index*np + p] *= controls.water_rate/current_rate;
                }
                break;
            }
            case Well::ProducerCMode::GRAT:
            {
                double current_rate = -well_state.wellRates()[ well_index*np + pu.phase_pos[Gas] ];

                if (current_rate == 0.0)
                    break;

                for (int p = 0; p<np; ++p) {
                    well_state.wellRates()[well_index*np + p] *= controls.gas_rate/current_rate;
                }
                break;

            }
            case Well::ProducerCMode::LRAT:
            {
                double current_rate = -well_state.wellRates()[ well_index*np + pu.phase_pos[Water] ]
                        - well_state.wellRates()[ well_index*np + pu.phase_pos[Oil] ];

                if (current_rate == 0.0)
                    break;

                for (int p = 0; p<np; ++p) {
                    well_state.wellRates()[well_index*np + p] *= controls.liquid_rate/current_rate;
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
                Base::rateConverter_.calcCoeff(/*fipreg*/ 0, Base::pvtRegionIdx_, convert_coeff);
                double total_res_rate = 0.0;
                for (int p = 0; p<np; ++p) {
                    total_res_rate -= well_state.wellRates()[well_index*np + p] * convert_coeff[p];
                }
                if (total_res_rate == 0.0)
                    break;

                if (controls.prediction_mode) {
                    for (int p = 0; p<np; ++p) {
                        well_state.wellRates()[well_index*np + p] *= controls.resv_rate/total_res_rate;
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
                    Base::rateConverter_.calcReservoirVoidageRates(/*fipreg*/ 0, Base::pvtRegionIdx_, hrates, hrates_resv);
                    double target = std::accumulate(hrates_resv.begin(), hrates_resv.end(), 0.0);
                    for (int p = 0; p<np; ++p) {
                        well_state.wellRates()[well_index*np + p] *= target/total_res_rate;
                    }

                }
                break;
            }
            case Well::ProducerCMode::BHP:
            {
                well_state.bhp()[well_index] = controls.bhp_limit;
                break;
            }
            case Well::ProducerCMode::THP:
            {
                well_state.thp()[well_index] = controls.thp_limit;
                auto bhp = computeBhpAtThpLimitProd(ebos_simulator, summaryState, deferred_logger);
                if (bhp) {
                    well_state.bhp()[well_index] = *bhp;
                } else {
                    deferred_logger.warning("FAILURE_GETTING_CONVERGED_BHP",
                                            "Failed to find BHP when switching to THP control for well " + name());
                }
                break;
            }
            case Well::ProducerCMode::GRUP:
            {
                //do nothing at the moment
                break;
            }
            case Well::ProducerCMode::CMODE_UNDEFINED:
            {
                OPM_DEFLOG_THROW(std::runtime_error, "Well control must be specified for well "  + name(), deferred_logger );
            }
            case Well::ProducerCMode::NONE:
            {
                OPM_DEFLOG_THROW(std::runtime_error, "Well control must be specified for well "  + name() , deferred_logger);
            }

                break;
            } // end of switch
        }
    }





    template<typename TypeTag>
    void
    StandardWell<TypeTag>::
    updateIPR(const Simulator& ebos_simulator, Opm::DeferredLogger& deferred_logger) const
    {
        // TODO: not handling solvent related here for now

        // TODO: it only handles the producers for now
        // the formular for the injectors are not formulated yet
        if (this->isInjector()) {
            return;
        }

        // initialize all the values to be zero to begin with
        std::fill(ipr_a_.begin(), ipr_a_.end(), 0.);
        std::fill(ipr_b_.begin(), ipr_b_.end(), 0.);

        for (int perf = 0; perf < number_of_perforations_; ++perf) {
            std::vector<EvalWell> mob(num_components_, {numWellEq_ + numEq, 0.0});
            // TODO: mabye we should store the mobility somewhere, so that we only need to calculate it one per iteration
            getMobility(ebos_simulator, perf, mob, deferred_logger);

            const int cell_idx = well_cells_[perf];
            const auto& int_quantities = *(ebos_simulator.model().cachedIntensiveQuantities(cell_idx, /*timeIdx=*/ 0));
            const auto& fs = int_quantities.fluidState();
            // the pressure of the reservoir grid block the well connection is in
            Eval perf_pressure = getPerfCellPressure(fs);
            double p_r = perf_pressure.value();

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
    StandardWell<TypeTag>::
    updateWellOperability(const Simulator& ebos_simulator,
                          const WellState& /* well_state */,
                          Opm::DeferredLogger& deferred_logger)
    {
        this->operability_status_.reset();

        updateIPR(ebos_simulator, deferred_logger);

        // checking the BHP limit related
        checkOperabilityUnderBHPLimitProducer(ebos_simulator, deferred_logger);

        const auto& summaryState = ebos_simulator.vanguard().summaryState();

        // checking whether the well can operate under the THP constraints.
        if (this->wellHasTHPConstraints(summaryState)) {
            checkOperabilityUnderTHPLimitProducer(ebos_simulator, deferred_logger);
        }
    }





    template<typename TypeTag>
    void
    StandardWell<TypeTag>::
    checkOperabilityUnderBHPLimitProducer(const Simulator& ebos_simulator, Opm::DeferredLogger& deferred_logger)
    {
        const auto& summaryState = ebos_simulator.vanguard().summaryState();
        const double bhp_limit = mostStrictBhpFromBhpLimits(summaryState);
        // Crude but works: default is one atmosphere.
        // TODO: a better way to detect whether the BHP is defaulted or not
        const bool bhp_limit_not_defaulted = bhp_limit > 1.5 * unit::barsa;
        if ( bhp_limit_not_defaulted || !this->wellHasTHPConstraints(summaryState) ) {
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
            if (this->operability_status_.operable_under_only_bhp_limit && this->wellHasTHPConstraints(summaryState)) {
                // option 1: calculate well rates based on the BHP limit.
                // option 2: stick with the above IPR curve
                // we use IPR here
                std::vector<double> well_rates_bhp_limit;
                computeWellRatesWithBhp(ebos_simulator, bhp_limit, well_rates_bhp_limit, deferred_logger);

                const double thp = calculateThpFromBhp(well_rates_bhp_limit, bhp_limit, deferred_logger);
                const double thp_limit = this->getTHPConstraint(summaryState);

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
    checkOperabilityUnderTHPLimitProducer(const Simulator& ebos_simulator, Opm::DeferredLogger& deferred_logger)
    {
        const auto& summaryState = ebos_simulator.vanguard().summaryState();
        const auto obtain_bhp = computeBhpAtThpLimitProd(ebos_simulator, summaryState, deferred_logger);

        if (obtain_bhp) {
            this->operability_status_.can_obtain_bhp_with_thp_limit = true;

            const double  bhp_limit = mostStrictBhpFromBhpLimits(summaryState);
            this->operability_status_.obey_bhp_limit_with_thp_limit = (*obtain_bhp >= bhp_limit);

            const double thp_limit = this->getTHPConstraint(summaryState);
            if (*obtain_bhp < thp_limit) {
                const std::string msg = " obtained bhp " + std::to_string(unit::convert::to(*obtain_bhp, unit::barsa))
                                        + " bars is SMALLER than thp limit "
                                        + std::to_string(unit::convert::to(thp_limit, unit::barsa))
                                        + " bars as a producer for well " + name();
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
            if ( (drawdown < 0. && this->isInjector()) ||
                 (drawdown > 0. && this->isProducer()) )  {
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
                                   Opm::DeferredLogger& deferred_logger)
    {
        const double bhp = well_state.bhp()[index_of_well_];
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
        return !getAllowCrossFlow() && allDrawDownWrongDirection(ebos_simulator);
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
            const double saltConcentration = fs.saltConcentration().value();

            if (waterPresent) {
                const unsigned waterCompIdx = Indices::canonicalToActiveComponentIndex(FluidSystem::waterCompIdx);
                b_perf[ waterCompIdx + perf * num_components_] =
                FluidSystem::waterPvt().inverseFormationVolumeFactor(fs.pvtRegionIndex(), temperature, p_avg, saltConcentration);
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
                std::fill(mix.begin(), mix.end(), 0.0);
                // No flow => use well specified fractions for mix.
                if (this->isInjector()) {
                    switch (this->wellEcl().injectorType()) {
                    case InjectorType::WATER:
                        mix[FluidSystem::waterCompIdx] = 1.0;
                        break;
                    case InjectorType::GAS:
                        mix[FluidSystem::gasCompIdx] = 1.0;
                        break;
                    case InjectorType::OIL:
                        mix[FluidSystem::oilCompIdx] = 1.0;
                        break;
                    case InjectorType::MULTI:
                        // Not supported.
                        // deferred_logger.warning("MULTI_PHASE_INJECTOR_NOT_SUPPORTED",
                        //                         "Multi phase injectors are not supported, requested for well " + name());
                        break;
                    }
                } else {
                    assert(this->isProducer());
                    // Using the preferred phase to decide the mix initialization.
                    switch (this->wellEcl().getPreferredPhase()) {
                    case Phase::OIL:
                        mix[FluidSystem::oilCompIdx] = 1.0;
                        break;
                    case Phase::GAS:
                        mix[FluidSystem::gasCompIdx] = 1.0;
                        break;
                    case Phase::WATER:
                        mix[FluidSystem::waterCompIdx] = 1.0;
                        break;
                    default:
                        // No others supported.
                        break;
                    }

                }
            }
            // Compute volume ratio.
            x = mix;

            // Subtract dissolved gas from oil phase and vapporized oil from gas phase
            if (FluidSystem::phaseIsActive(FluidSystem::gasCompIdx) && FluidSystem::phaseIsActive(FluidSystem::oilCompIdx)) {
                const unsigned gaspos = Indices::canonicalToActiveComponentIndex(FluidSystem::gasCompIdx);
                const unsigned oilpos = Indices::canonicalToActiveComponentIndex(FluidSystem::oilCompIdx);
                double rs = 0.0;
                double rv = 0.0;
                if (!rsmax_perf.empty() && mix[oilpos] > 1e-12) {
                    rs = std::min(mix[gaspos]/mix[oilpos], rsmax_perf[perf]);
                }
                if (!rvmax_perf.empty() && mix[gaspos] > 1e-12) {
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
    getWellConvergence(const WellState& well_state,
                       const std::vector<double>& B_avg,
                       Opm::DeferredLogger& deferred_logger,
                       const bool /*relax_tolerance*/) const
    {
        // the following implementation assume that the polymer is always after the w-o-g phases
        // For the polymer, energy and foam cases, there is one more mass balance equations of reservoir than wells
        assert((int(B_avg.size()) == num_components_) || has_polymer || has_energy || has_foam || has_brine);

        const double tol_wells = param_.tolerance_wells_;
        const double maxResidualAllowed = param_.max_residual_allowed_;

        std::vector<double> res(numWellEq_);
        for (int eq_idx = 0; eq_idx < numWellEq_; ++eq_idx) {
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
            const int compIdx = Indices::canonicalToActiveComponentIndex(canonicalCompIdx);

            if (std::isnan(well_flux_residual[compIdx])) {
                report.setWellFailed({type, CR::Severity::NotANumber, compIdx, name()});
            } else if (well_flux_residual[compIdx] > maxResidualAllowed) {
                report.setWellFailed({type, CR::Severity::TooLarge, compIdx, name()});
            } else if (well_flux_residual[compIdx] > tol_wells) {
                report.setWellFailed({type, CR::Severity::Normal, compIdx, name()});
            }
        }

        checkConvergenceControlEq(well_state, report, deferred_logger);

        checkConvergenceExtraEqs(res, report);

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
        if (!this->isOperable() && !this->wellIsStopped()) return;

        // We assemble the well equations, then we check the convergence,
        // which is why we do not put the assembleWellEq here.
        BVectorWell dx_well(1);
        dx_well[0].resize(numWellEq_);
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
        if (!this->isOperable() && !this->wellIsStopped()) return;

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
        if (!this->isOperable() && !this->wellIsStopped()) return;

        assert( invDrw_.size() == invDuneD_.N() );

        // invDrw_ = invDuneD_ * resWell_
        invDuneD_.mv(resWell_, invDrw_);
        // r = r - duneC_^T * invDrw_
        duneC_.mmtv(invDrw_, r);
    }

#if HAVE_CUDA || HAVE_OPENCL
    template<typename TypeTag>
    void
    StandardWell<TypeTag>::
    addWellContribution(WellContributions& wellContribs) const
    {
        std::vector<int> colIndices;
        std::vector<double> nnzValues;
        colIndices.reserve(duneB_.nonzeroes());
        nnzValues.reserve(duneB_.nonzeroes()*numStaticWellEq * numEq);

        // duneC
        for ( auto colC = duneC_[0].begin(), endC = duneC_[0].end(); colC != endC; ++colC )
        {
            colIndices.emplace_back(colC.index());
            for (int i = 0; i < numStaticWellEq; ++i) {
                for (int j = 0; j < numEq; ++j) {
                    nnzValues.emplace_back((*colC)[i][j]);
                }
            }
        }
        wellContribs.addMatrix(WellContributions::MatrixType::C, colIndices.data(), nnzValues.data(), duneC_.nonzeroes());

        // invDuneD
        colIndices.clear();
        nnzValues.clear();
        colIndices.emplace_back(0);
        for (int i = 0; i < numStaticWellEq; ++i)
        {
            for (int j = 0; j < numStaticWellEq; ++j) {
                nnzValues.emplace_back(invDuneD_[0][0][i][j]);
            }
        }
        wellContribs.addMatrix(WellContributions::MatrixType::D, colIndices.data(), nnzValues.data(), 1);

        // duneB
        colIndices.clear();
        nnzValues.clear();
        for ( auto colB = duneB_[0].begin(), endB = duneB_[0].end(); colB != endB; ++colB )
        {
            colIndices.emplace_back(colB.index());
            for (int i = 0; i < numStaticWellEq; ++i) {
                for (int j = 0; j < numEq; ++j) {
                    nnzValues.emplace_back((*colB)[i][j]);
                }
            }
        }
        wellContribs.addMatrix(WellContributions::MatrixType::B, colIndices.data(), nnzValues.data(), duneB_.nonzeroes());
    }


    template<typename TypeTag>
    void
    StandardWell<TypeTag>::
    getNumBlocks(unsigned int& numBlocks) const
    {
        numBlocks = duneB_.nonzeroes();
    }
#endif


    template<typename TypeTag>
    void
    StandardWell<TypeTag>::
    recoverSolutionWell(const BVector& x, BVectorWell& xw) const
    {
        if (!this->isOperable() && !this->wellIsStopped()) return;

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
        if (!this->isOperable() && !this->wellIsStopped()) return;

        BVectorWell xw(1);
        xw[0].resize(numWellEq_);

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
            std::vector<EvalWell> mob(num_components_, {numWellEq_ + numEq, 0.});
            getMobility(ebosSimulator, perf, mob, deferred_logger);
            double trans_mult = ebosSimulator.problem().template rockCompTransMultiplier<double>(intQuants, cell_idx);
            const double Tw = well_index_[perf] * trans_mult;

            std::vector<EvalWell> cq_s(num_components_, {numWellEq_ + numEq, 0.});
            double perf_dis_gas_rate = 0.;
            double perf_vap_oil_rate = 0.;
            computePerfRate(intQuants, mob, EvalWell(numWellEq_ + numEq, bhp), Tw, perf, allow_cf,
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

        // iterate to get a more accurate well density
        // create a copy of the well_state to use. If the operability checking is sucessful, we use this one
        // to replace the original one
        WellState well_state_copy = ebosSimulator.problem().wellModel().wellState();

        //  Set current control to bhp, and bhp value in state, modify bhp limit in control object.
        if (well_ecl_.isInjector()) {
            well_state_copy.currentInjectionControls()[index_of_well_] = Well::InjectorCMode::BHP;
        } else {
            well_state_copy.currentProductionControls()[index_of_well_] = Well::ProducerCMode::BHP;
        }
        well_state_copy.bhp()[index_of_well_] = bhp;

        const double dt = ebosSimulator.timeStepSize();
        bool converged = this->iterateWellEquations(ebosSimulator, B_avg, dt, well_state_copy, deferred_logger);
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
    }




    template<typename TypeTag>
    std::vector<double>
    StandardWell<TypeTag>::
    computeWellPotentialWithTHP(const Simulator& ebos_simulator,
                                Opm::DeferredLogger& deferred_logger) const
    {
        std::vector<double> potentials(number_of_phases_, 0.0);
        const auto& summary_state = ebos_simulator.vanguard().summaryState();

        const auto& well = well_ecl_;
        if (well.isInjector()){
            auto bhp_at_thp_limit = computeBhpAtThpLimitInj(ebos_simulator, summary_state, deferred_logger);
            if (bhp_at_thp_limit) {
                const auto& controls = well_ecl_.injectionControls(summary_state);
                const double bhp = std::min(*bhp_at_thp_limit, controls.bhp_limit);
                computeWellRatesWithBhp(ebos_simulator, bhp, potentials, deferred_logger);
            } else {
                deferred_logger.warning("FAILURE_GETTING_CONVERGED_POTENTIAL",
                                        "Failed in getting converged thp based potential calculation for well "
                                        + name() + ". Instead the bhp based value is used");
                const auto& controls = well_ecl_.injectionControls(summary_state);
                const double bhp = controls.bhp_limit;
                computeWellRatesWithBhp(ebos_simulator, bhp, potentials, deferred_logger);
            }
        } else {
            auto bhp_at_thp_limit = computeBhpAtThpLimitProd(ebos_simulator, summary_state, deferred_logger);
            if (bhp_at_thp_limit) {
                const auto& controls = well_ecl_.productionControls(summary_state);
                const double bhp = std::max(*bhp_at_thp_limit, controls.bhp_limit);
                computeWellRatesWithBhp(ebos_simulator, bhp, potentials, deferred_logger);
            } else {
                deferred_logger.warning("FAILURE_GETTING_CONVERGED_POTENTIAL",
                                        "Failed in getting converged thp based potential calculation for well "
                                        + name() + ". Instead the bhp based value is used");
                const auto& controls = well_ecl_.productionControls(summary_state);
                const double bhp = controls.bhp_limit;
                computeWellRatesWithBhp(ebos_simulator, bhp, potentials, deferred_logger);
            }
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
        const int np = number_of_phases_;
        well_potentials.resize(np, 0.0);

        if (this->wellIsStopped()) {
            return;
        }

        // creating a copy of the well itself, to avoid messing up the explicit informations
        // during this copy, the only information not copied properly is the well controls
        StandardWell<TypeTag> well(*this);
        well.calculateExplicitQuantities(ebosSimulator, well_state, deferred_logger);

        // does the well have a THP related constraint?
        const auto& summaryState = ebosSimulator.vanguard().summaryState();
        const Well::ProducerCMode& current_control = well_state.currentProductionControls()[this->index_of_well_];
        if ( !well.Base::wellHasTHPConstraints(summaryState) || current_control == Well::ProducerCMode::BHP ) {
            // get the bhp value based on the bhp constraints
            const double bhp = well.mostStrictBhpFromBhpLimits(summaryState);
            assert(std::abs(bhp) != std::numeric_limits<double>::max());
            well.computeWellRatesWithBhpPotential(ebosSimulator, B_avg, bhp, well_potentials, deferred_logger);
        } else {
            // the well has a THP related constraint
            well_potentials = well.computeWellPotentialWithTHP(ebosSimulator, deferred_logger);
        }
    }





    template<typename TypeTag>
    void
    StandardWell<TypeTag>::
    updatePrimaryVariables(const WellState& well_state, Opm::DeferredLogger& deferred_logger) const
    {
        if (!this->isOperable() && !this->wellIsStopped()) return;

        const int well_index = index_of_well_;
        const int np = number_of_phases_;
        const auto& pu = phaseUsage();

        // the weighted total well rate
        double total_well_rate = 0.0;
        for (int p = 0; p < np; ++p) {
            total_well_rate += scalingFactor(p) * well_state.wellRates()[np * well_index + p];
        }

        // Not: for the moment, the first primary variable for the injectors is not G_total. The injection rate
        // under surface condition is used here
        if (this->isInjector()) {
            switch (this->wellEcl().injectorType()) {
            case InjectorType::WATER:
                primary_variables_[WQTotal] = well_state.wellRates()[np * well_index + pu.phase_pos[Water]];
                break;
            case InjectorType::GAS:
                primary_variables_[WQTotal] = well_state.wellRates()[np * well_index + pu.phase_pos[Gas]];
                break;
            case InjectorType::OIL:
                primary_variables_[WQTotal] = well_state.wellRates()[np * well_index + pu.phase_pos[Oil]];
                break;
            case InjectorType::MULTI:
                // Not supported.
                deferred_logger.warning("MULTI_PHASE_INJECTOR_NOT_SUPPORTED",
                                        "Multi phase injectors are not supported, requested for well " + name());
                break;
            }
        } else {
            for (int p = 0; p < np; ++p) {
                primary_variables_[WQTotal] = total_well_rate;
            }
        }

        if (std::abs(total_well_rate) > 0.) {
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
            if (this->isInjector()) {
                auto phase = well_ecl_.getInjectionProperties().injectorType;
                // only single phase injection handled
                if (FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx)) {
                    if (phase == InjectorType::WATER) {
                        primary_variables_[WFrac] = 1.0;
                    } else {
                        primary_variables_[WFrac] = 0.0;
                    }
                }

                if (FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx)) {
                    if (phase == InjectorType::GAS) {
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
            } else if (this->isProducer()) { // producers
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

        // other primary variables related to polymer injection
        if (this->has_polymermw && this->isInjector()) {
            for (int perf = 0; perf < number_of_perforations_; ++perf) {
                primary_variables_[Bhp + 1 + perf] = well_state.perfWaterVelocity()[first_perf_ + perf];
                primary_variables_[Bhp + 1 + number_of_perforations_ + perf] = well_state.perfSkinPressure()[first_perf_ + perf];
            }
        }
#ifndef NDEBUG
        for (double v : primary_variables_) {
            assert(Opm::isfinite(v));
        }
#endif
    }



    template<typename TypeTag>
    template<class ValueType>
    ValueType
    StandardWell<TypeTag>::
    calculateBhpFromThp(const std::vector<ValueType>& rates,
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

        // pick the density in the top layer
        // TODO: it is possible it should be a Evaluation
        const double rho = perf_densities_[0];

        if (this->isInjector() )
        {
            const auto& controls = well.injectionControls(summaryState);
            const double vfp_ref_depth = vfp_properties_->getInj()->getTable(controls.vfp_table_number)->getDatumDepth();
            const double dp = wellhelpers::computeHydrostaticCorrection(ref_depth_, vfp_ref_depth, rho, gravity_);
            return vfp_properties_->getInj()->bhp(controls.vfp_table_number, aqua, liquid, vapour, controls.thp_limit) - dp;
         }
         else if (this->isProducer()) {
             const auto& controls = well.productionControls(summaryState);
             const double vfp_ref_depth = vfp_properties_->getProd()->getTable(controls.vfp_table_number)->getDatumDepth();
             const double dp = wellhelpers::computeHydrostaticCorrection(ref_depth_, vfp_ref_depth, rho, gravity_);
             return vfp_properties_->getProd()->bhp(controls.vfp_table_number, aqua, liquid, vapour, controls.thp_limit, controls.alq_value) - dp;
         }
         else {
             OPM_DEFLOG_THROW(std::logic_error, "Expected INJECTOR or PRODUCER well", deferred_logger);
         }



    }



    template<typename TypeTag>
    double
    StandardWell<TypeTag>::
    calculateThpFromBhp(const std::vector<double>& rates,
                        const double bhp,
                        Opm::DeferredLogger& deferred_logger) const
    {
        assert(int(rates.size()) == 3); // the vfp related only supports three phases now.

        const double aqua = rates[Water];
        const double liquid = rates[Oil];
        const double vapour = rates[Gas];

        // pick the density in the top layer
        const double rho = perf_densities_[0];

        double thp = 0.0;
        if (this->isInjector()) {
            const int table_id = well_ecl_.vfp_table_number();
            const double vfp_ref_depth = vfp_properties_->getInj()->getTable(table_id)->getDatumDepth();
            const double dp = wellhelpers::computeHydrostaticCorrection(ref_depth_, vfp_ref_depth, rho, gravity_);

            thp = vfp_properties_->getInj()->thp(table_id, aqua, liquid, vapour, bhp + dp);
        }
        else if (this->isProducer()) {
            const int table_id = well_ecl_.vfp_table_number();
            const double alq = well_ecl_.alq_value();
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
        // for the cases related to polymer molecular weight, we assume fully mixing
        // as a result, the polymer and water share the same viscosity
        if (this->has_polymermw) {
            return;
        }
        const int cell_idx = well_cells_[perf];
        const auto& int_quant = *(ebos_simulator.model().cachedIntensiveQuantities(cell_idx, /*timeIdx=*/ 0));
        const EvalWell polymer_concentration = extendEval(int_quant.polymerConcentration());

        // TODO: not sure should based on the well type or injecting/producing peforations
        // it can be different for crossflow
        if (this->isInjector()) {
            // assume fully mixing within injecting wellbore
            const auto& visc_mult_table = PolymerModule::plyviscViscosityMultiplierTable(int_quant.pvtRegionIndex());
            const unsigned waterCompIdx = Indices::canonicalToActiveComponentIndex(FluidSystem::waterCompIdx);
            mob[waterCompIdx] /= (extendEval(int_quant.waterViscosityCorrection()) * visc_mult_table.eval(polymer_concentration, /*extrapolate=*/true) );
        }

        if (PolymerModule::hasPlyshlog()) {
            // we do not calculate the shear effects for injection wells when they do not
            // inject polymer.
            if (this->isInjector() && wpolymer() == 0.) {
                return;
            }
            // compute the well water velocity with out shear effects.
            // TODO: do we need to turn on crossflow here?
            const bool allow_cf = getAllowCrossFlow() || openCrossFlowAvoidSingularity(ebos_simulator);
            const EvalWell& bhp = getBhp();

            std::vector<EvalWell> cq_s(num_components_, {numWellEq_ + numEq, 0.});
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
    StandardWell<TypeTag>::addWellContributions(SparseMatrixAdapter& jacobian) const
    {
        // We need to change matrx A as follows
        // A -= C^T D^-1 B
        // D is diagonal
        // B and C have 1 row, nc colums and nonzero
        // at (0,j) only if this well has a perforation at cell j.
        typename SparseMatrixAdapter::MatrixBlock tmpMat;
        Dune::DynamicMatrix<Scalar> tmp;
        for ( auto colC = duneC_[0].begin(), endC = duneC_[0].end(); colC != endC; ++colC )
        {
            const auto row_index = colC.index();

            for ( auto colB = duneB_[0].begin(), endB = duneB_[0].end(); colB != endB; ++colB )
            {
                Detail::multMatrix(invDuneD_[0][0],  (*colB), tmp);
                Detail::negativeMultMatrixTransposed((*colC), tmp, tmpMat);
                jacobian.addToBlock( row_index, colB.index(), tmpMat );
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

        if (FluidSystem::numActivePhases() > 1) {
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
        }
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
        const bool converged = this->iterateWellEquations(ebos_simulator, B_avg, dt, well_state_copy, deferred_logger);

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
    typename StandardWell<TypeTag>::EvalWell
    StandardWell<TypeTag>::
    pskinwater(const double throughput,
               const EvalWell& water_velocity,
              Opm::DeferredLogger& deferred_logger) const
    {
        if (!this->has_polymermw) {
            OPM_DEFLOG_THROW(std::runtime_error, "Polymermw is not activated, "
                                          "while injecting skin pressure is requested for well " << name(), deferred_logger);
        }
        const int water_table_id = well_ecl_.getPolymerProperties().m_skprwattable;
        if (water_table_id <= 0) {
            OPM_DEFLOG_THROW(std::runtime_error, "Unused SKPRWAT table id used for well " << name(), deferred_logger);
        }
        const auto& water_table_func = PolymerModule::getSkprwatTable(water_table_id);
        const EvalWell throughput_eval(numWellEq_ + numEq, throughput);
        // the skin pressure when injecting water, which also means the polymer concentration is zero
        EvalWell pskin_water(numWellEq_ + numEq, 0.0);
        pskin_water = water_table_func.eval(throughput_eval, water_velocity);
        return pskin_water;
    }





    template<typename TypeTag>
    typename StandardWell<TypeTag>::EvalWell
    StandardWell<TypeTag>::
    pskin(const double throughput,
              const EvalWell& water_velocity,
              const EvalWell& poly_inj_conc,
              Opm::DeferredLogger& deferred_logger) const
    {
        if (!this->has_polymermw) {
            OPM_DEFLOG_THROW(std::runtime_error, "Polymermw is not activated, "
                                          "while injecting skin pressure is requested for well " << name(), deferred_logger);
        }
        const double sign = water_velocity >= 0. ? 1.0 : -1.0;
        const EvalWell water_velocity_abs = Opm::abs(water_velocity);
        if (poly_inj_conc == 0.) {
            return sign * pskinwater(throughput, water_velocity_abs, deferred_logger);
        }
        const int polymer_table_id = well_ecl_.getPolymerProperties().m_skprpolytable;
        if (polymer_table_id <= 0) {
            OPM_DEFLOG_THROW(std::runtime_error, "Unavailable SKPRPOLY table id used for well " << name(), deferred_logger);
        }
        const auto& skprpolytable = PolymerModule::getSkprpolyTable(polymer_table_id);
        const double reference_concentration = skprpolytable.refConcentration;
        const EvalWell throughput_eval(numWellEq_ + numEq, throughput);
        // the skin pressure when injecting water, which also means the polymer concentration is zero
        EvalWell pskin_poly(numWellEq_ + numEq, 0.0);
        pskin_poly = skprpolytable.table_func.eval(throughput_eval, water_velocity_abs);
        if (poly_inj_conc == reference_concentration) {
            return sign * pskin_poly;
        }
        // poly_inj_conc != reference concentration of the table, then some interpolation will be required
        const EvalWell pskin_water = pskinwater(throughput, water_velocity_abs, deferred_logger);
        const EvalWell pskin = pskin_water + (pskin_poly - pskin_water) / reference_concentration * poly_inj_conc;
        return sign * pskin;
    }





    template<typename TypeTag>
    typename StandardWell<TypeTag>::EvalWell
    StandardWell<TypeTag>::
    wpolymermw(const double throughput,
               const EvalWell& water_velocity,
              Opm::DeferredLogger& deferred_logger) const
    {
        if (!this->has_polymermw) {
            OPM_DEFLOG_THROW(std::runtime_error, "Polymermw is not activated, "
                                          "while injecting polymer molecular weight is requested for well " << name(), deferred_logger);
        }
        const int table_id = well_ecl_.getPolymerProperties().m_plymwinjtable;
        const auto& table_func = PolymerModule::getPlymwinjTable(table_id);
        const EvalWell throughput_eval(numWellEq_ + numEq, throughput);
        EvalWell molecular_weight(numWellEq_ + numEq, 0.);
        if (wpolymer() == 0.) { // not injecting polymer
            return molecular_weight;
        }
        molecular_weight = table_func.eval(throughput_eval, Opm::abs(water_velocity));
        return molecular_weight;
    }





    template<typename TypeTag>
    void
    StandardWell<TypeTag>::
    updateWaterThroughput(const double dt, WellState &well_state) const
    {
        if (this->has_polymermw && this->isInjector()) {
            for (int perf = 0; perf < number_of_perforations_; ++perf) {
                const double perf_water_vel = primary_variables_[Bhp + 1 + perf];
                // we do not consider the formation damage due to water flowing from reservoir into wellbore
                if (perf_water_vel > 0.) {
                    well_state.perfThroughput()[first_perf_ + perf] += perf_water_vel * dt;
                }
            }
        }
    }





    template<typename TypeTag>
    void
    StandardWell<TypeTag>::
    handleInjectivityRateAndEquations(const IntensiveQuantities& int_quants,
                                      const WellState& well_state,
                                      const int perf,
                                      std::vector<EvalWell>& cq_s,
                                      Opm::DeferredLogger& deferred_logger)
    {
        const unsigned water_comp_idx = Indices::canonicalToActiveComponentIndex(FluidSystem::waterCompIdx);
        const EvalWell& water_flux_s = cq_s[water_comp_idx];
        const auto& fs = int_quants.fluidState();
        const EvalWell b_w = extendEval(fs.invB(FluidSystem::waterPhaseIdx));
        const EvalWell water_flux_r = water_flux_s / b_w;
        const double area = M_PI * bore_diameters_[perf] * perf_length_[perf];
        const EvalWell water_velocity = water_flux_r / area;
        const int wat_vel_index = Bhp + 1 + perf;

        // equation for the water velocity
        const EvalWell eq_wat_vel = primary_variables_evaluation_[wat_vel_index] - water_velocity;
        resWell_[0][wat_vel_index] = eq_wat_vel.value();

        const double throughput = well_state.perfThroughput()[first_perf_ + perf];
        const int pskin_index = Bhp + 1 + number_of_perforations_ + perf;

        EvalWell poly_conc(numWellEq_ + numEq, 0.0);
        poly_conc.setValue(wpolymer());

        // equation for the skin pressure
        const EvalWell eq_pskin = primary_variables_evaluation_[pskin_index]
                                  - pskin(throughput, primary_variables_evaluation_[wat_vel_index], poly_conc, deferred_logger);

        resWell_[0][pskin_index] = eq_pskin.value();
        for (int pvIdx = 0; pvIdx < numWellEq_; ++pvIdx) {
            invDuneD_[0][0][wat_vel_index][pvIdx] = eq_wat_vel.derivative(pvIdx+numEq);
            invDuneD_[0][0][pskin_index][pvIdx] = eq_pskin.derivative(pvIdx+numEq);
        }

        // water rate is update to use the form from water velocity, since water velocity is
        // a primary variable now
        cq_s[water_comp_idx] = area * primary_variables_evaluation_[wat_vel_index] * b_w;

        // the water velocity is impacted by the reservoir primary varaibles. It needs to enter matrix B
        const int cell_idx = well_cells_[perf];
        for (int pvIdx = 0; pvIdx < numEq; ++pvIdx) {
            duneB_[0][cell_idx][wat_vel_index][pvIdx] = eq_wat_vel.derivative(pvIdx);
        }
    }





    template<typename TypeTag>
    void
    StandardWell<TypeTag>::
    checkConvergenceControlEq(const WellState& well_state,
                              ConvergenceReport& report,
                              DeferredLogger& deferred_logger) const
    {
        double control_tolerance = 0.;
        using CR = ConvergenceReport;
        CR::WellFailure::Type ctrltype = CR::WellFailure::Type::Invalid;

        const int well_index = index_of_well_;
        if (wellIsStopped_) {
            ctrltype = CR::WellFailure::Type::ControlRate;
            control_tolerance = 1.e-6; // use smaller tolerance for zero control?
        }
        else if (this->isInjector() )
        {
            const Opm::Well::InjectorCMode& current = well_state.currentInjectionControls()[well_index];
            switch(current) {
            case Well::InjectorCMode::THP:
                ctrltype = CR::WellFailure::Type::ControlTHP;
                control_tolerance = 1.e4; // 0.1 bar
                break;
            case Well::InjectorCMode::BHP:
                ctrltype = CR::WellFailure::Type::ControlBHP;
                control_tolerance = 1.e3; // 0.01 bar
                break;
            case Well::InjectorCMode::RATE:
            case Well::InjectorCMode::RESV:
                ctrltype = CR::WellFailure::Type::ControlRate;
                control_tolerance = 1.e-4; //
                break;
            case Well::InjectorCMode::GRUP:
                ctrltype = CR::WellFailure::Type::ControlRate;
                control_tolerance = 1.e-6; //
                break;
            default:
                OPM_DEFLOG_THROW(std::runtime_error, "Unknown well control control types for well " << name(), deferred_logger);
            }
        }
        else if (this->isProducer() )
        {
            const Well::ProducerCMode& current = well_state.currentProductionControls()[well_index];
            switch(current) {
            case Well::ProducerCMode::THP:
                ctrltype = CR::WellFailure::Type::ControlTHP;
                control_tolerance = 1.e4; // 0.1 bar
                break;
            case Well::ProducerCMode::BHP:
                ctrltype = CR::WellFailure::Type::ControlBHP;
                control_tolerance = 1.e3; // 0.01 bar
                break;
            case Well::ProducerCMode::ORAT:
            case Well::ProducerCMode::WRAT:
            case Well::ProducerCMode::GRAT:
            case Well::ProducerCMode::LRAT:
            case Well::ProducerCMode::RESV:
            case Well::ProducerCMode::CRAT:
                ctrltype = CR::WellFailure::Type::ControlRate;
                control_tolerance = 1.e-4; // smaller tolerance for rate control
                break;
            case Well::ProducerCMode::GRUP:
                ctrltype = CR::WellFailure::Type::ControlRate;
                control_tolerance = 1.e-6; // smaller tolerance for rate control
                break;
            default:
                OPM_DEFLOG_THROW(std::runtime_error, "Unknown well control control types for well " << name(), deferred_logger);
            }
        }

        const double well_control_residual = std::abs(resWell_[0][Bhp]);
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





    template<typename TypeTag>
    void
    StandardWell<TypeTag>::
    checkConvergenceExtraEqs(const std::vector<double>& res,
                             ConvergenceReport& report) const
    {
        // if different types of extra equations are involved, this function needs to be refactored further

        // checking the convergence of the extra equations related to polymer injectivity
        if (this->has_polymermw && this->isInjector()) {
            //  checking the convergence of the perforation rates
            const double wat_vel_tol = 1.e-8;
            const int dummy_component = -1;
            const double maxResidualAllowed = param_.max_residual_allowed_;
            using CR = ConvergenceReport;
            const auto wat_vel_failure_type = CR::WellFailure::Type::MassBalance;
            for (int perf = 0; perf < number_of_perforations_; ++perf) {
                const double wat_vel_residual = res[Bhp + 1 + perf];
                if (std::isnan(wat_vel_residual)) {
                    report.setWellFailed({wat_vel_failure_type, CR::Severity::NotANumber, dummy_component, name()});
                } else if (wat_vel_residual > maxResidualAllowed * 10.) {
                    report.setWellFailed({wat_vel_failure_type, CR::Severity::TooLarge, dummy_component, name()});
                } else if (wat_vel_residual > wat_vel_tol) {
                    report.setWellFailed({wat_vel_failure_type, CR::Severity::Normal, dummy_component, name()});
                }
            }

            // checking the convergence of the skin pressure
            const double pskin_tol = 1000.; // 1000 pascal
            const auto pskin_failure_type = CR::WellFailure::Type::Pressure;
            for (int perf = 0; perf < number_of_perforations_; ++perf) {
                const double pskin_residual = res[Bhp + 1 + perf + number_of_perforations_];
                if (std::isnan(pskin_residual)) {
                    report.setWellFailed({pskin_failure_type, CR::Severity::NotANumber, dummy_component, name()});
                } else if (pskin_residual > maxResidualAllowed * 10.) {
                    report.setWellFailed({pskin_failure_type, CR::Severity::TooLarge, dummy_component, name()});
                } else if (pskin_residual > pskin_tol) {
                    report.setWellFailed({pskin_failure_type, CR::Severity::Normal, dummy_component, name()});
                }
            }
        }
    }





    template<typename TypeTag>
    void
    StandardWell<TypeTag>::
    updateConnectionRatePolyMW(const EvalWell& cq_s_poly,
                               const IntensiveQuantities& int_quants,
                               const WellState& well_state,
                               const int perf,
                               DeferredLogger& deferred_logger)
    {
        // the source term related to transport of molecular weight
        EvalWell cq_s_polymw = cq_s_poly;
        if (this->isInjector()) {
            const int wat_vel_index = Bhp + 1 + perf;
            const EvalWell water_velocity = primary_variables_evaluation_[wat_vel_index];
            if (water_velocity > 0.) { // injecting
                const double throughput = well_state.perfThroughput()[first_perf_ + perf];
                const EvalWell molecular_weight = wpolymermw(throughput, water_velocity, deferred_logger);
                cq_s_polymw *= molecular_weight;
            } else {
                // we do not consider the molecular weight from the polymer
                // going-back to the wellbore through injector
                cq_s_polymw *= 0.;
            }
        } else if (this->isProducer()) {
            if (cq_s_polymw < 0.) {
                cq_s_polymw *= extendEval(int_quants.polymerMoleWeight() );
            } else {
                // we do not consider the molecular weight from the polymer
                // re-injecting back through producer
                cq_s_polymw *= 0.;
            }
        }
        connectionRates_[perf][this->contiPolymerMWEqIdx] = Base::restrictEval(cq_s_polymw);
    }






    template<typename TypeTag>
    std::optional<double>
    StandardWell<TypeTag>::
    computeBhpAtThpLimitProd(const Simulator& ebos_simulator,
                             const SummaryState& summary_state,
                             DeferredLogger& deferred_logger) const
    {
        // Given a VFP function returning bhp as a function of phase
        // rates and thp:
        //     fbhp(rates, thp),
        // a function extracting the particular flow rate used for VFP
        // lookups:
        //     flo(rates)
        // and the inflow function (assuming the reservoir is fixed):
        //     frates(bhp)
        // we want to solve the equation:
        //     fbhp(frates(bhp, thplimit)) - bhp = 0
        // for bhp.
        //
        // This may result in 0, 1 or 2 solutions. If two solutions,
        // the one corresponding to the lowest bhp (and therefore
        // highest rate) is returned.
        //
        // In order to detect these situations, we will find piecewise
        // linear approximations both to the inverse of the frates
        // function and to the fbhp function.
        //
        // We first take the FLO sample points of the VFP curve, and
        // find the corresponding bhp values by solving the equation:
        //     flo(frates(bhp)) - flo_sample = 0
        // for bhp, for each flo_sample. The resulting (flo_sample,
        // bhp_sample) values give a piecewise linear approximation to
        // the true inverse inflow function, at the same flo values as
        // the VFP data.
        //
        // Then we extract a piecewise linear approximation from the
        // multilinear fbhp() by evaluating it at the flo_sample
        // points, with fractions given by the frates(bhp_sample)
        // values.
        //
        // When we have both piecewise linear curves defined on the
        // same flo_sample points, it is easy to distinguish between
        // the 0, 1 or 2 solution cases, and obtain the right interval
        // in which to solve for the solution we want (with highest
        // flow in case of 2 solutions).

        // Make the fbhp() function.
        const auto& controls = well_ecl_.productionControls(summary_state);
        const auto& table = *(vfp_properties_->getProd()->getTable(controls.vfp_table_number));
        const double vfp_ref_depth = table.getDatumDepth();
        const double rho = perf_densities_[0]; // Use the density at the top perforation.
        const double dp = wellhelpers::computeHydrostaticCorrection(ref_depth_, vfp_ref_depth, rho, gravity_);
        auto fbhp = [this, &controls, dp](const std::vector<double>& rates) {
            assert(rates.size() == 3);
            return this->vfp_properties_->getProd()
            ->bhp(controls.vfp_table_number, rates[Water], rates[Oil], rates[Gas], controls.thp_limit, controls.alq_value) - dp;
        };

        // Make the flo() function.
        auto flo_type = table.getFloType();
        auto flo = [flo_type](const std::vector<double>& rates) {
            return detail::getFlo(rates[Water], rates[Oil], rates[Gas], flo_type);
        };

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

        // Get the flo samples, add extra samples at low rates and bhp
        // limit point if necessary. Then the sign must be flipped
        // since the VFP code expects that production flo values are
        // negative.
        std::vector<double> flo_samples = table.getFloAxis();
        if (flo_samples[0] > 0.0) {
            const double f0 = flo_samples[0];
            flo_samples.insert(flo_samples.begin(), { f0/20.0, f0/10.0, f0/5.0, f0/2.0 });
        }
        const double flo_bhp_limit = -flo(frates(controls.bhp_limit));
        if (flo_samples.back() < flo_bhp_limit) {
            flo_samples.push_back(flo_bhp_limit);
        }
        for (double& x : flo_samples) {
            x = -x;
        }

        // Find bhp values for inflow relation corresponding to flo samples.
        std::vector<double> bhp_samples;
        for (double flo_sample : flo_samples) {
            if (flo_sample < -flo_bhp_limit) {
                // We would have to go under the bhp limit to obtain a
                // flow of this magnitude. We associate all such flows
                // with simply the bhp limit. The first one
                // encountered is considered valid, the rest not. They
                // are therefore skipped.
                bhp_samples.push_back(controls.bhp_limit);
                break;
            }
            auto eq = [&flo, &frates, flo_sample](double bhp) {
                return flo(frates(bhp)) - flo_sample;
            };
            // TODO: replace hardcoded low/high limits.
            const double low = 10.0 * unit::barsa;
            const double high = 600.0 * unit::barsa;
            const int max_iteration = 50;
            const double flo_tolerance = 1e-6 * std::fabs(flo_samples.back());
            int iteration = 0;
            try {
                const double solved_bhp = RegulaFalsiBisection<>::
                    solve(eq, low, high, max_iteration, flo_tolerance, iteration);
                bhp_samples.push_back(solved_bhp);
            }
            catch (...) {
                // Use previous value (or max value if at start) if we failed.
                bhp_samples.push_back(bhp_samples.empty() ? high : bhp_samples.back());
                deferred_logger.warning("FAILED_ROBUST_BHP_THP_SOLVE_EXTRACT_SAMPLES",
                                        "Robust bhp(thp) solve failed extracting bhp values at flo samples for well " + name());
            }
        }

        // Find bhp values for VFP relation corresponding to flo samples.
        const int num_samples = bhp_samples.size(); // Note that this can be smaller than flo_samples.size()
        std::vector<double> fbhp_samples(num_samples);
        for (int ii = 0; ii < num_samples; ++ii) {
            fbhp_samples[ii] = fbhp(frates(bhp_samples[ii]));
        }
// #define EXTRA_THP_DEBUGGING
#ifdef EXTRA_THP_DEBUGGING
        std::string dbgmsg;
        dbgmsg += "flo: ";
        for (int ii = 0; ii < num_samples; ++ii) {
            dbgmsg += "  " + std::to_string(flo_samples[ii]);
        }
        dbgmsg += "\nbhp: ";
        for (int ii = 0; ii < num_samples; ++ii) {
            dbgmsg += "  " + std::to_string(bhp_samples[ii]);
        }
        dbgmsg += "\nfbhp: ";
        for (int ii = 0; ii < num_samples; ++ii) {
            dbgmsg += "  " + std::to_string(fbhp_samples[ii]);
        }
        OpmLog::debug(dbgmsg);
#endif // EXTRA_THP_DEBUGGING

        // Look for sign changes for the (fbhp_samples - bhp_samples) piecewise linear curve.
        // We only look at the valid
        int sign_change_index = -1;
        for (int ii = 0; ii < num_samples - 1; ++ii) {
            const double curr = fbhp_samples[ii] - bhp_samples[ii];
            const double next = fbhp_samples[ii + 1] - bhp_samples[ii + 1];
            if (curr * next < 0.0) {
                // Sign change in the [ii, ii + 1] interval.
                sign_change_index = ii; // May overwrite, thereby choosing the highest-flo solution.
            }
        }

        // Handle the no solution case.
        if (sign_change_index == -1) {
            return std::optional<double>();
        }

        // Solve for the proper solution in the given interval.
        auto eq = [&fbhp, &frates](double bhp) {
            return fbhp(frates(bhp)) - bhp;
        };
        // TODO: replace hardcoded low/high limits.
        const double low = bhp_samples[sign_change_index + 1];
        const double high = bhp_samples[sign_change_index];
        const int max_iteration = 50;
        const double bhp_tolerance = 0.01 * unit::barsa;
        int iteration = 0;
        if (low == high) {
            // We are in the high flow regime where the bhp_samples
            // are all equal to the bhp_limit.
            assert(low == controls.bhp_limit);
            deferred_logger.warning("FAILED_ROBUST_BHP_THP_SOLVE",
                                    "Robust bhp(thp) solve failed for well " + name());
            return std::optional<double>();
        }
        try {
            const double solved_bhp = RegulaFalsiBisection<>::
                solve(eq, low, high, max_iteration, bhp_tolerance, iteration);
#ifdef EXTRA_THP_DEBUGGING
            OpmLog::debug("*****    " + name() + "    solved_bhp = " + std::to_string(solved_bhp)
                          + "    flo_bhp_limit = " + std::to_string(flo_bhp_limit));
#endif // EXTRA_THP_DEBUGGING
            return solved_bhp;
        }
        catch (...) {
            deferred_logger.warning("FAILED_ROBUST_BHP_THP_SOLVE",
                                    "Robust bhp(thp) solve failed for well " + name());
            return std::optional<double>();
        }

    }



    template<typename TypeTag>
    std::optional<double>
    StandardWell<TypeTag>::
    computeBhpAtThpLimitInj(const Simulator& ebos_simulator,
                            const SummaryState& summary_state,
                            DeferredLogger& deferred_logger) const
    {
        // Given a VFP function returning bhp as a function of phase
        // rates and thp:
        //     fbhp(rates, thp),
        // a function extracting the particular flow rate used for VFP
        // lookups:
        //     flo(rates)
        // and the inflow function (assuming the reservoir is fixed):
        //     frates(bhp)
        // we want to solve the equation:
        //     fbhp(frates(bhp, thplimit)) - bhp = 0
        // for bhp.
        //
        // This may result in 0, 1 or 2 solutions. If two solutions,
        // the one corresponding to the lowest bhp (and therefore
        // highest rate) is returned.
        //
        // In order to detect these situations, we will find piecewise
        // linear approximations both to the inverse of the frates
        // function and to the fbhp function.
        //
        // We first take the FLO sample points of the VFP curve, and
        // find the corresponding bhp values by solving the equation:
        //     flo(frates(bhp)) - flo_sample = 0
        // for bhp, for each flo_sample. The resulting (flo_sample,
        // bhp_sample) values give a piecewise linear approximation to
        // the true inverse inflow function, at the same flo values as
        // the VFP data.
        //
        // Then we extract a piecewise linear approximation from the
        // multilinear fbhp() by evaluating it at the flo_sample
        // points, with fractions given by the frates(bhp_sample)
        // values.
        //
        // When we have both piecewise linear curves defined on the
        // same flo_sample points, it is easy to distinguish between
        // the 0, 1 or 2 solution cases, and obtain the right interval
        // in which to solve for the solution we want (with highest
        // flow in case of 2 solutions).

        // Make the fbhp() function.
        const auto& controls = well_ecl_.injectionControls(summary_state);
        const auto& table = *(vfp_properties_->getInj()->getTable(controls.vfp_table_number));
        const double vfp_ref_depth = table.getDatumDepth();
        const double rho = perf_densities_[0]; // Use the density at the top perforation.
        const double dp = wellhelpers::computeHydrostaticCorrection(ref_depth_, vfp_ref_depth, rho, gravity_);
        auto fbhp = [this, &controls, dp](const std::vector<double>& rates) {
            assert(rates.size() == 3);
            return this->vfp_properties_->getInj()
                    ->bhp(controls.vfp_table_number, rates[Water], rates[Oil], rates[Gas], controls.thp_limit) - dp;
        };

        // Make the flo() function.
        auto flo_type = table.getFloType();
        auto flo = [flo_type](const std::vector<double>& rates) {
            return detail::getFlo(rates[Water], rates[Oil], rates[Gas], flo_type);
        };

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

        // Get the flo samples, add extra samples at low rates and bhp
        // limit point if necessary.
        std::vector<double> flo_samples = table.getFloAxis();
        if (flo_samples[0] > 0.0) {
            const double f0 = flo_samples[0];
            flo_samples.insert(flo_samples.begin(), { f0/20.0, f0/10.0, f0/5.0, f0/2.0 });
        }
        const double flo_bhp_limit = flo(frates(controls.bhp_limit));
        if (flo_samples.back() < flo_bhp_limit) {
            flo_samples.push_back(flo_bhp_limit);
        }

        // Find bhp values for inflow relation corresponding to flo samples.
        std::vector<double> bhp_samples;
        for (double flo_sample : flo_samples) {
            if (flo_sample > flo_bhp_limit) {
                // We would have to go over the bhp limit to obtain a
                // flow of this magnitude. We associate all such flows
                // with simply the bhp limit. The first one
                // encountered is considered valid, the rest not. They
                // are therefore skipped.
                bhp_samples.push_back(controls.bhp_limit);
                break;
            }
            auto eq = [&flo, &frates, flo_sample](double bhp) {
                return flo(frates(bhp)) - flo_sample;
            };
            // TODO: replace hardcoded low/high limits.
            const double low = 10.0 * unit::barsa;
            const double high = 800.0 * unit::barsa;
            const int max_iteration = 50;
            const double flo_tolerance = 1e-6 * std::fabs(flo_samples.back());
            int iteration = 0;
            try {
                const double solved_bhp = RegulaFalsiBisection<>::
                        solve(eq, low, high, max_iteration, flo_tolerance, iteration);
                bhp_samples.push_back(solved_bhp);
            }
            catch (...) {
                // Use previous value (or max value if at start) if we failed.
                bhp_samples.push_back(bhp_samples.empty() ? low : bhp_samples.back());
                deferred_logger.warning("FAILED_ROBUST_BHP_THP_SOLVE_EXTRACT_SAMPLES",
                                        "Robust bhp(thp) solve failed extracting bhp values at flo samples for well " + name());
            }
        }

        // Find bhp values for VFP relation corresponding to flo samples.
        const int num_samples = bhp_samples.size(); // Note that this can be smaller than flo_samples.size()
        std::vector<double> fbhp_samples(num_samples);
        for (int ii = 0; ii < num_samples; ++ii) {
            fbhp_samples[ii] = fbhp(frates(bhp_samples[ii]));
        }
// #define EXTRA_THP_DEBUGGING
#ifdef EXTRA_THP_DEBUGGING
        std::string dbgmsg;
        dbgmsg += "flo: ";
        for (int ii = 0; ii < num_samples; ++ii) {
            dbgmsg += "  " + std::to_string(flo_samples[ii]);
        }
        dbgmsg += "\nbhp: ";
        for (int ii = 0; ii < num_samples; ++ii) {
            dbgmsg += "  " + std::to_string(bhp_samples[ii]);
        }
        dbgmsg += "\nfbhp: ";
        for (int ii = 0; ii < num_samples; ++ii) {
            dbgmsg += "  " + std::to_string(fbhp_samples[ii]);
        }
        OpmLog::debug(dbgmsg);
#endif // EXTRA_THP_DEBUGGING

        // Look for sign changes for the (fbhp_samples - bhp_samples) piecewise linear curve.
        // We only look at the valid
        int sign_change_index = -1;
        for (int ii = 0; ii < num_samples - 1; ++ii) {
            const double curr = fbhp_samples[ii] - bhp_samples[ii];
            const double next = fbhp_samples[ii + 1] - bhp_samples[ii + 1];
            if (curr * next < 0.0) {
                // Sign change in the [ii, ii + 1] interval.
                sign_change_index = ii; // May overwrite, thereby choosing the highest-flo solution.
            }
        }

        // Handle the no solution case.
        if (sign_change_index == -1) {
            return std::optional<double>();
        }

        // Solve for the proper solution in the given interval.
        auto eq = [&fbhp, &frates](double bhp) {
            return fbhp(frates(bhp)) - bhp;
        };
        // TODO: replace hardcoded low/high limits.
        const double low = bhp_samples[sign_change_index + 1];
        const double high = bhp_samples[sign_change_index];
        const int max_iteration = 50;
        const double bhp_tolerance = 0.01 * unit::barsa;
        int iteration = 0;
        if (low == high) {
            // We are in the high flow regime where the bhp_samples
            // are all equal to the bhp_limit.
            assert(low == controls.bhp_limit);
            deferred_logger.warning("FAILED_ROBUST_BHP_THP_SOLVE",
                                    "Robust bhp(thp) solve failed for well " + name());
            return std::optional<double>();
        }
        try {
            const double solved_bhp = RegulaFalsiBisection<>::
                    solve(eq, low, high, max_iteration, bhp_tolerance, iteration);
#ifdef EXTRA_THP_DEBUGGING
            OpmLog::debug("*****    " + name() + "    solved_bhp = " + std::to_string(solved_bhp)
                          + "    flo_bhp_limit = " + std::to_string(flo_bhp_limit));
#endif // EXTRA_THP_DEBUGGING
            return solved_bhp;
        }
        catch (...) {
            deferred_logger.warning("FAILED_ROBUST_BHP_THP_SOLVE",
                                    "Robust bhp(thp) solve failed for well " + name());
            return std::optional<double>();
        }

    }





    template<typename TypeTag>
    bool
    StandardWell<TypeTag>::
    iterateWellEqWithControl(const Simulator& ebosSimulator,
                             const std::vector<double>& B_avg,
                             const double dt,
                             const Well::InjectionControls& inj_controls,
                             const Well::ProductionControls& prod_controls,
                             WellState& well_state,
                             Opm::DeferredLogger& deferred_logger)
    {
        const int max_iter = param_.max_inner_iter_wells_;
        int it = 0;
        bool converged;
        do {
            assembleWellEqWithoutIteration(ebosSimulator, dt, inj_controls, prod_controls, well_state, deferred_logger);

            auto report = getWellConvergence(well_state, B_avg, deferred_logger);

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


} // namespace Opm
