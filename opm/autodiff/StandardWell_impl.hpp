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


namespace Opm
{
    template<typename TypeTag>
    StandardWell<TypeTag>::
    StandardWell(const Well* well, const int time_step, const Wells* wells, const ModelParameters& param)
    : Base(well, time_step, wells, param)
    , perf_densities_(number_of_perforations_)
    , perf_pressure_diffs_(number_of_perforations_)
    , primary_variables_(numWellEq, 0.0)
    , primary_variables_evaluation_(numWellEq) // the number of the primary variables
    , F0_(numWellEq)
    {
        duneB_.setBuildMode( OffDiagMatWell::row_wise );
        duneC_.setBuildMode( OffDiagMatWell::row_wise );
        invDuneD_.setBuildMode( DiagMatWell::row_wise );
    }





    template<typename TypeTag>
    void
    StandardWell<TypeTag>::
    init(const PhaseUsage* phase_usage_arg,
         const std::vector<bool>* active_arg,
         const std::vector<double>& depth_arg,
         const double gravity_arg,
         const int num_cells)
    {
        Base::init(phase_usage_arg, active_arg,
                   depth_arg, gravity_arg, num_cells);

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
        // TODO: using numComp here is only to make the 2p + dummy phase work
        // TODO: in theory, we should use numWellEq here.
        // for (int eqIdx = 0; eqIdx < numWellEq; ++eqIdx) {
        for (int eqIdx = 0; eqIdx < numComponents(); ++eqIdx) {
            assert( (size_t)eqIdx < primary_variables_.size() );

            primary_variables_evaluation_[eqIdx] = 0.0;
            primary_variables_evaluation_[eqIdx].setValue(primary_variables_[eqIdx]);
            primary_variables_evaluation_[eqIdx].setDerivative(numEq + eqIdx, 1.0);
        }
    }





    template<typename TypeTag>
    typename StandardWell<TypeTag>::EvalWell
    StandardWell<TypeTag>::
    getBhp() const
    {
        const WellControls* wc = well_controls_;
        if (well_controls_get_current_type(wc) == BHP) {
            EvalWell bhp = 0.0;
            const double target_rate = well_controls_get_current_target(wc);
            bhp.setValue(target_rate);
            return bhp;
        } else if (well_controls_get_current_type(wc) == THP) {
            const int control = well_controls_get_current(wc);

            const Opm::PhaseUsage& pu = phaseUsage();
            std::vector<EvalWell> rates(3, 0.0);
            if (active()[ Water ]) {
                rates[ Water ]= getQs(pu.phase_pos[ Water]);
            }
            if (active()[ Oil ]) {
                rates[ Oil ] = getQs(pu.phase_pos[ Oil ]);
            }
            if (active()[ Gas ]) {
                rates[ Gas ] = getQs(pu.phase_pos[ Gas ]);
            }
            return calculateBhpFromThp(rates, control);
        }

        return primary_variables_evaluation_[XvarWell];
    }





    template<typename TypeTag>
    typename StandardWell<TypeTag>::EvalWell
    StandardWell<TypeTag>::
    getQs(const int comp_idx) const // TODO: phase or component?
    {
        EvalWell qs = 0.0;

        const WellControls* wc = well_controls_;
        const int np = number_of_phases_;
        const double target_rate = well_controls_get_current_target(wc);

        assert(comp_idx < numComponents());
        const auto pu = phaseUsage();

        // TODO: the formulation for the injectors decides it only work with single phase
        // surface rate injection control. Improvement will be required.
        if (well_type_ == INJECTOR) {
            if (has_solvent) {
                // TODO: investigate whether the use of the comp_frac is justified.
                // The usage of the comp_frac is not correct, which should be changed later.
                double comp_frac = 0.0;
                if (has_solvent && comp_idx == contiSolventEqIdx) { // solvent
                    comp_frac = comp_frac_[pu.phase_pos[ Gas ]] * wsolvent();
                } else if (comp_idx == pu.phase_pos[ Gas ]) {
                    comp_frac = comp_frac_[comp_idx] * (1.0 - wsolvent());
                } else {
                    comp_frac = comp_frac_[comp_idx];
                }
                if (comp_frac == 0.0) {
                    return qs; //zero
                }

                if (well_controls_get_current_type(wc) == BHP || well_controls_get_current_type(wc) == THP) {
                    return comp_frac * primary_variables_evaluation_[XvarWell];
                }

                qs.setValue(comp_frac * target_rate);
                return qs;
            }

            const double comp_frac = comp_frac_[comp_idx];
            if (comp_frac == 0.0) {
                return qs;
            }

            if (well_controls_get_current_type(wc) == BHP || well_controls_get_current_type(wc) == THP) {
                return primary_variables_evaluation_[XvarWell];
            }
            qs.setValue(target_rate);
            return qs;
        }

        // Producers
        if (well_controls_get_current_type(wc) == BHP || well_controls_get_current_type(wc) == THP ) {
            return primary_variables_evaluation_[XvarWell] * wellVolumeFractionScaled(comp_idx);
        }

        if (well_controls_get_current_type(wc) == SURFACE_RATE) {
            // checking how many phases are included in the rate control
            // to decide wheter it is a single phase rate control or not
            const double* distr = well_controls_get_current_distr(wc);
            int num_phases_under_rate_control = 0;
            for (int phase = 0; phase < np; ++phase) {
                if (distr[phase] > 0.0) {
                    num_phases_under_rate_control += 1;
                }
            }

            // there should be at least one phase involved
            assert(num_phases_under_rate_control > 0);

            // when it is a single phase rate limit
            if (num_phases_under_rate_control == 1) {

                // looking for the phase under control
                int phase_under_control = -1;
                for (int phase = 0; phase < np; ++phase) {
                    if (distr[phase] > 0.0) {
                        phase_under_control = phase;
                        break;
                    }
                }

                assert(phase_under_control >= 0);

                EvalWell wellVolumeFractionScaledPhaseUnderControl = wellVolumeFractionScaled(phase_under_control);
                if (has_solvent && phase_under_control == Gas) {
                    // for GRAT controlled wells solvent is included in the target
                    wellVolumeFractionScaledPhaseUnderControl += wellVolumeFractionScaled(contiSolventEqIdx);
                }

                if (comp_idx == phase_under_control) {
                    if (has_solvent && phase_under_control == Gas) {
                        qs.setValue(target_rate * wellVolumeFractionScaled(Gas).value() / wellVolumeFractionScaledPhaseUnderControl.value() );
                        return qs;
                    }
                    qs.setValue(target_rate);
                    return qs;
                }

                // TODO: not sure why the single phase under control will have near zero fraction
                const double eps = 1e-6;
                if (wellVolumeFractionScaledPhaseUnderControl < eps) {
                    return qs;
                }
                return (target_rate * wellVolumeFractionScaled(comp_idx) / wellVolumeFractionScaledPhaseUnderControl);
            }

            // when it is a combined two phase rate limit, such like LRAT
            // we neec to calculate the rate for the certain phase
            if (num_phases_under_rate_control == 2) {
                EvalWell combined_volume_fraction = 0.;
                for (int p = 0; p < np; ++p) {
                    if (distr[p] == 1.0) {
                        combined_volume_fraction += wellVolumeFractionScaled(p);
                    }
                }
                return (target_rate * wellVolumeFractionScaled(comp_idx) / combined_volume_fraction);
            }

            // TODO: three phase surface rate control is not tested yet
            if (num_phases_under_rate_control == 3) {
                return target_rate * wellSurfaceVolumeFraction(comp_idx);
            }
        } else if (well_controls_get_current_type(wc) == RESERVOIR_RATE) {
            // ReservoirRate
            return target_rate * wellVolumeFractionScaled(comp_idx);
        } else {
            OPM_THROW(std::logic_error, "Unknown control type for well " << name());
        }

        // avoid warning of condition reaches end of non-void function
        return qs;
    }






    template<typename TypeTag>
    typename StandardWell<TypeTag>::EvalWell
    StandardWell<TypeTag>::
    wellVolumeFractionScaled(const int compIdx) const
    {

        const double scal = scalingFactor(compIdx);
        if (scal > 0)
            return  wellVolumeFraction(compIdx) / scal;

        // the scaling factor may be zero for RESV controlled wells.
        return wellVolumeFraction(compIdx);
    }





    template<typename TypeTag>
    typename StandardWell<TypeTag>::EvalWell
    StandardWell<TypeTag>::
    wellVolumeFraction(const int compIdx) const
    {
        const auto& pu = phaseUsage();
        if (active()[Water] && compIdx == pu.phase_pos[Water]) {
            return primary_variables_evaluation_[WFrac];
        }

        if (active()[Gas] && compIdx == pu.phase_pos[Gas]) {
            return primary_variables_evaluation_[GFrac];
        }

        if (has_solvent && compIdx == contiSolventEqIdx) {
            return primary_variables_evaluation_[SFrac];
        }

        // Oil fraction
        EvalWell well_fraction = 1.0;
        if (active()[Water]) {
            well_fraction -= primary_variables_evaluation_[WFrac];
        }

        if (active()[Gas]) {
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
        const int numComp = numComponents();
        for (int idx = 0; idx < numComp; ++idx) {
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
                    const std::vector<EvalWell>& mob_perfcells_dense,
                    const double Tw, const EvalWell& bhp, const double& cdp,
                    const bool& allow_cf, std::vector<EvalWell>& cq_s) const
    {
        const Opm::PhaseUsage& pu = phaseUsage();
        const int np = number_of_phases_;
        const int numComp = numComponents();
        std::vector<EvalWell> cmix_s(numComp,0.0);
        for (int componentIdx = 0; componentIdx < numComp; ++componentIdx) {
            cmix_s[componentIdx] = wellSurfaceVolumeFraction(componentIdx);
        }
        const auto& fs = intQuants.fluidState();

        const EvalWell pressure = extendEval(fs.pressure(FluidSystem::oilPhaseIdx));
        const EvalWell rs = extendEval(fs.Rs());
        const EvalWell rv = extendEval(fs.Rv());
        std::vector<EvalWell> b_perfcells_dense(numComp, 0.0);
        for (int phase = 0; phase < np; ++phase) {
            const int ebosPhaseIdx = flowPhaseToEbosPhaseIdx(phase);
            b_perfcells_dense[phase] = extendEval(fs.invB(ebosPhaseIdx));
        }
        if (has_solvent) {
            b_perfcells_dense[contiSolventEqIdx] = extendEval(intQuants.solventInverseFormationVolumeFactor());
        }

        // Pressure drawdown (also used to determine direction of flow)
        const EvalWell well_pressure = bhp + cdp;
        const EvalWell drawdown = pressure - well_pressure;

        // producing perforations
        if ( drawdown.value() > 0 )  {
            //Do nothing if crossflow is not allowed
            if (!allow_cf && well_type_ == INJECTOR) {
                return;
            }

            // compute component volumetric rates at standard conditions
            for (int componentIdx = 0; componentIdx < numComp; ++componentIdx) {
                const EvalWell cq_p = - Tw * (mob_perfcells_dense[componentIdx] * drawdown);
                cq_s[componentIdx] = b_perfcells_dense[componentIdx] * cq_p;
            }

            if (active()[Oil] && active()[Gas]) {
                const int oilpos = pu.phase_pos[Oil];
                const int gaspos = pu.phase_pos[Gas];
                const EvalWell cq_sOil = cq_s[oilpos];
                const EvalWell cq_sGas = cq_s[gaspos];
                cq_s[gaspos] += rs * cq_sOil;
                cq_s[oilpos] += rv * cq_sGas;
            }

        } else {
            //Do nothing if crossflow is not allowed
            if (!allow_cf && well_type_ == PRODUCER) {
                return;
            }

            // Using total mobilities
            EvalWell total_mob_dense = mob_perfcells_dense[0];
            for (int componentIdx = 1; componentIdx < numComp; ++componentIdx) {
                total_mob_dense += mob_perfcells_dense[componentIdx];
            }

            // injection perforations total volume rates
            const EvalWell cqt_i = - Tw * (total_mob_dense * drawdown);

            // compute volume ratio between connection at standard conditions
            EvalWell volumeRatio = 0.0;
            if (active()[Water]) {
                const int watpos = pu.phase_pos[Water];
                volumeRatio += cmix_s[watpos] / b_perfcells_dense[watpos];
            }

            if (has_solvent) {
                volumeRatio += cmix_s[contiSolventEqIdx] / b_perfcells_dense[contiSolventEqIdx];
            }

            if (active()[Oil] && active()[Gas]) {
                const int oilpos = pu.phase_pos[Oil];
                const int gaspos = pu.phase_pos[Gas];

                // Incorporate RS/RV factors if both oil and gas active
                const EvalWell d = 1.0 - rv * rs;

                if (d.value() == 0.0) {
                    OPM_THROW(Opm::NumericalProblem, "Zero d value obtained for well " << name() << " during flux calcuation"
                                                  << " with rs " << rs << " and rv " << rv);
                }

                const EvalWell tmp_oil = (cmix_s[oilpos] - rv * cmix_s[gaspos]) / d;
                //std::cout << "tmp_oil " <<tmp_oil << std::endl;
                volumeRatio += tmp_oil / b_perfcells_dense[oilpos];

                const EvalWell tmp_gas = (cmix_s[gaspos] - rs * cmix_s[oilpos]) / d;
                //std::cout << "tmp_gas " <<tmp_gas << std::endl;
                volumeRatio += tmp_gas / b_perfcells_dense[gaspos];
            }
            else {
                if (active()[Oil]) {
                    const int oilpos = pu.phase_pos[Oil];
                    volumeRatio += cmix_s[oilpos] / b_perfcells_dense[oilpos];
                }
                if (active()[Gas]) {
                    const int gaspos = pu.phase_pos[Gas];
                    volumeRatio += cmix_s[gaspos] / b_perfcells_dense[gaspos];
                }
            }

            // injecting connections total volumerates at standard conditions
            EvalWell cqt_is = cqt_i/volumeRatio;
            //std::cout << "volrat " << volumeRatio << " " << volrat_perf_[perf] << std::endl;
            for (int componentIdx = 0; componentIdx < numComp; ++componentIdx) {
                cq_s[componentIdx] = cmix_s[componentIdx] * cqt_is; // * b_perfcells_dense[phase];
            }
        }
    }





    template<typename TypeTag>
    void
    StandardWell<TypeTag>::
    assembleWellEq(Simulator& ebosSimulator,
                   const double dt,
                   WellState& well_state,
                   bool only_wells)
    {
        const int numComp = numComponents();
        const int np = number_of_phases_;

        // clear all entries
        if (!only_wells) {
            duneB_ = 0.0;
            duneC_ = 0.0;
        }
        invDuneD_ = 0.0;
        resWell_ = 0.0;

        auto& ebosJac = ebosSimulator.model().linearizer().matrix();
        auto& ebosResid = ebosSimulator.model().linearizer().residual();

        // TODO: it probably can be static member for StandardWell
        const double volume = 0.002831684659200; // 0.1 cu ft;

        const bool allow_cf = crossFlowAllowed(ebosSimulator);

        const EvalWell& bhp = getBhp();

        for (int perf = 0; perf < number_of_perforations_; ++perf) {

            const int cell_idx = well_cells_[perf];
            const auto& intQuants = *(ebosSimulator.model().cachedIntensiveQuantities(cell_idx, /*timeIdx=*/ 0));
            std::vector<EvalWell> cq_s(numComp,0.0);
            std::vector<EvalWell> mob(numComp, 0.0);
            getMobility(ebosSimulator, perf, mob);
            computePerfRate(intQuants, mob, well_index_[perf], bhp, perf_pressure_diffs_[perf], allow_cf, cq_s);

            for (int componentIdx = 0; componentIdx < numComp; ++componentIdx) {
                // the cq_s entering mass balance equations need to consider the efficiency factors.
                const EvalWell cq_s_effective = cq_s[componentIdx] * well_efficiency_factor_;

                if (!only_wells) {
                    // subtract sum of component fluxes in the reservoir equation.
                    // need to consider the efficiency factor
                    ebosResid[cell_idx][flowPhaseToEbosCompIdx(componentIdx)] -= cq_s_effective.value();
                }

                // subtract sum of phase fluxes in the well equations.
                resWell_[0][componentIdx] -= cq_s_effective.value();

                // assemble the jacobians
                for (int pvIdx = 0; pvIdx < numWellEq; ++pvIdx) {
                    if (!only_wells) {
                        // also need to consider the efficiency factor when manipulating the jacobians.
                        duneC_[0][cell_idx][pvIdx][flowPhaseToEbosCompIdx(componentIdx)] -= cq_s_effective.derivative(pvIdx+numEq); // intput in transformed matrix
                    }
                    invDuneD_[0][0][componentIdx][pvIdx] -= cq_s_effective.derivative(pvIdx+numEq);
                }

                for (int pvIdx = 0; pvIdx < numEq; ++pvIdx) {
                    if (!only_wells) {
                        // also need to consider the efficiency factor when manipulating the jacobians.
                        ebosJac[cell_idx][cell_idx][flowPhaseToEbosCompIdx(componentIdx)][pvIdx] -= cq_s_effective.derivative(pvIdx);
                        duneB_[0][cell_idx][componentIdx][pvIdx] -= cq_s_effective.derivative(pvIdx);
                    }
                }

                // Store the perforation phase flux for later usage.
                if (has_solvent && componentIdx == contiSolventEqIdx) {// if (flowPhaseToEbosCompIdx(componentIdx) == Solvent)
                    well_state.perfRateSolvent()[first_perf_ + perf] = cq_s[componentIdx].value();
                } else {
                    well_state.perfPhaseRates()[(first_perf_ + perf) * np + componentIdx] = cq_s[componentIdx].value();
                }
            }

            if (has_polymer) {
                EvalWell cq_s_poly = cq_s[Water];
                if (well_type_ == INJECTOR) {
                    cq_s_poly *= wpolymer();
                } else {
                    cq_s_poly *= extendEval(intQuants.polymerConcentration() * intQuants.polymerViscosityCorrection());
                }
                if (!only_wells) {
                    // TODO: we need to consider the efficiency here.
                    for (int pvIdx = 0; pvIdx < numEq; ++pvIdx) {
                        ebosJac[cell_idx][cell_idx][contiPolymerEqIdx][pvIdx] -= cq_s_poly.derivative(pvIdx);
                    }
                    ebosResid[cell_idx][contiPolymerEqIdx] -= cq_s_poly.value();
                }
            }

            // Store the perforation pressure for later usage.
            well_state.perfPress()[first_perf_ + perf] = well_state.bhp()[index_of_well_] + perf_pressure_diffs_[perf];
        }

        // add vol * dF/dt + Q to the well equations;
        for (int componentIdx = 0; componentIdx < numComp; ++componentIdx) {
            EvalWell resWell_loc = (wellSurfaceVolumeFraction(componentIdx) - F0_[componentIdx]) * volume / dt;
            resWell_loc += getQs(componentIdx) * well_efficiency_factor_;
            for (int pvIdx = 0; pvIdx < numWellEq; ++pvIdx) {
                invDuneD_[0][0][componentIdx][pvIdx] += resWell_loc.derivative(pvIdx+numEq);
            }
            resWell_[0][componentIdx] += resWell_loc.value();
        }

        // do the local inversion of D.
        invDuneD_[0][0].invert();
    }





    template<typename TypeTag>
    bool
    StandardWell<TypeTag>::
    crossFlowAllowed(const Simulator& ebosSimulator) const
    {
        if (getAllowCrossFlow()) {
            return true;
        }

        // TODO: investigate the justification of the following situation

        // check for special case where all perforations have cross flow
        // then the wells must allow for cross flow
        for (int perf = 0; perf < number_of_perforations_; ++perf) {
            const int cell_idx = well_cells_[perf];
            const auto& intQuants = *(ebosSimulator.model().cachedIntensiveQuantities(cell_idx, /*timeIdx=*/0));
            const auto& fs = intQuants.fluidState();
            EvalWell pressure = extendEval(fs.pressure(FluidSystem::oilPhaseIdx));
            EvalWell bhp = getBhp();

            // Pressure drawdown (also used to determine direction of flow)
            EvalWell well_pressure = bhp + perf_pressure_diffs_[perf];
            EvalWell drawdown = pressure - well_pressure;

            if (drawdown.value() < 0 && well_type_ == INJECTOR)  {
                return false;
            }

            if (drawdown.value() > 0 && well_type_ == PRODUCER)  {
                return false;
            }
        }
        return true;
    }





    template<typename TypeTag>
    void
    StandardWell<TypeTag>::
    getMobility(const Simulator& ebosSimulator,
                const int perf,
                std::vector<EvalWell>& mob) const
    {
        const int np = number_of_phases_;
        const int cell_idx = well_cells_[perf];
        assert (int(mob.size()) == numComponents());
        const auto& intQuants = *(ebosSimulator.model().cachedIntensiveQuantities(cell_idx, /*timeIdx=*/0));
        const auto& materialLawManager = ebosSimulator.problem().materialLawManager();

        // either use mobility of the perforation cell or calcualte its own
        // based on passing the saturation table index
        const int satid = saturation_table_number_[perf] - 1;
        const int satid_elem = materialLawManager->satnumRegionIdx(cell_idx);
        if( satid == satid_elem ) { // the same saturation number is used. i.e. just use the mobilty from the cell

            for (int phase = 0; phase < np; ++phase) {
                int ebosPhaseIdx = flowPhaseToEbosPhaseIdx(phase);
                mob[phase] = extendEval(intQuants.mobility(ebosPhaseIdx));
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
            for (int phase = 0; phase < np; ++phase) {
                int ebosPhaseIdx = flowPhaseToEbosPhaseIdx(phase);
                mob[phase] = extendEval(relativePerms[ebosPhaseIdx] / intQuants.fluidState().viscosity(ebosPhaseIdx));
            }

            // this may not work if viscosity and relperms has been modified?
            if (has_solvent) {
                OPM_THROW(std::runtime_error, "individual mobility for wells does not work in combination with solvent");
            }
        }

        // modify the water mobility if polymer is present
        if (has_polymer) {
            // assume fully mixture for wells.
            EvalWell polymerConcentration = extendEval(intQuants.polymerConcentration());

            if (well_type_ == INJECTOR) {
                const auto& viscosityMultiplier = PolymerModule::plyviscViscosityMultiplierTable(intQuants.pvtRegionIndex());
                mob[ Water ] /= (extendEval(intQuants.waterViscosityCorrection()) * viscosityMultiplier.eval(polymerConcentration, /*extrapolate=*/true) );
            }

            if (PolymerModule::hasPlyshlog()) {
                // compute the well water velocity with out shear effects.
                const int numComp = numComponents();
                const bool allow_cf = crossFlowAllowed(ebosSimulator);
                const EvalWell& bhp = getBhp();
                std::vector<EvalWell> cq_s(numComp,0.0);
                computePerfRate(intQuants, mob, well_index_[perf], bhp, perf_pressure_diffs_[perf], allow_cf, cq_s);
                // TODO: make area a member
                double area = 2 * M_PI * perf_rep_radius_[perf] * perf_length_[perf];
                const auto& materialLawManager = ebosSimulator.problem().materialLawManager();
                const auto& scaledDrainageInfo =
                        materialLawManager->oilWaterScaledEpsInfoDrainage(cell_idx);
                const Scalar& Swcr = scaledDrainageInfo.Swcr;
                const EvalWell poro = extendEval(intQuants.porosity());
                const EvalWell Sw = extendEval(intQuants.fluidState().saturation(flowPhaseToEbosPhaseIdx(Water)));
                // guard against zero porosity and no water
                const EvalWell denom = Opm::max( (area * poro * (Sw - Swcr)), 1e-12);
                EvalWell waterVelocity = cq_s[ Water ] / denom * extendEval(intQuants.fluidState().invB(flowPhaseToEbosPhaseIdx(Water)));

                if (PolymerModule::hasShrate()) {
                    // TODO Use the same conversion as for the reservoar equations.
                    // Need the "permeability" of the well?
                    // For now use the same formula as in legacy.
                    waterVelocity *= PolymerModule::shrate( intQuants.pvtRegionIndex() ) / bore_diameters_[perf];
                }
                EvalWell polymerConcentration = extendEval(intQuants.polymerConcentration());
                EvalWell shearFactor = PolymerModule::computeShearFactor(polymerConcentration,
                                                                         intQuants.pvtRegionIndex(),
                                                                         waterVelocity);

                // modify the mobility with the shear factor and recompute the well fluxes.
                mob[ Water ] /= shearFactor;
            }
        }
    }





    template<typename TypeTag>
    void
    StandardWell<TypeTag>::
    updateWellState(const BVectorWell& dwells,
                    WellState& well_state) const
    {
        const int np = number_of_phases_;
        const double dBHPLimit = param_.dbhp_max_rel_;
        const double dFLimit = param_.dwell_fraction_max_;
        const auto pu = phaseUsage();

        const std::vector<double> xvar_well_old = primary_variables_;

        // update the second and third well variable (The flux fractions)
        std::vector<double> F(np,0.0);
        if (active()[ Water ]) {
            const int sign2 = dwells[0][WFrac] > 0 ? 1: -1;
            const double dx2_limited = sign2 * std::min(std::abs(dwells[0][WFrac]),dFLimit);
            primary_variables_[WFrac] = xvar_well_old[WFrac] - dx2_limited;
        }

        if (active()[ Gas ]) {
            const int sign3 = dwells[0][GFrac] > 0 ? 1: -1;
            const double dx3_limited = sign3 * std::min(std::abs(dwells[0][GFrac]),dFLimit);
            primary_variables_[GFrac] = xvar_well_old[GFrac] - dx3_limited;
        }

        if (has_solvent) {
            const int sign4 = dwells[0][SFrac] > 0 ? 1: -1;
            const double dx4_limited = sign4 * std::min(std::abs(dwells[0][SFrac]),dFLimit);
            primary_variables_[SFrac] = xvar_well_old[SFrac] - dx4_limited;
        }

        assert(active()[ Oil ]);
        F[pu.phase_pos[Oil]] = 1.0;

        if (active()[ Water ]) {
            F[pu.phase_pos[Water]] = primary_variables_[WFrac];
            F[pu.phase_pos[Oil]] -= F[pu.phase_pos[Water]];
        }

        if (active()[ Gas ]) {
            F[pu.phase_pos[Gas]] = primary_variables_[GFrac];
            F[pu.phase_pos[Oil]] -= F[pu.phase_pos[Gas]];
        }

        double F_solvent = 0.0;
        if (has_solvent) {
            F_solvent = primary_variables_[SFrac];
            F[pu.phase_pos[Oil]] -= F_solvent;
        }

        if (active()[ Water ]) {
            if (F[Water] < 0.0) {
                if (active()[ Gas ]) {
                        F[pu.phase_pos[Gas]] /= (1.0 - F[pu.phase_pos[Water]]);
                }
                if (has_solvent) {
                    F_solvent /= (1.0 - F[pu.phase_pos[Water]]);
                }
                F[pu.phase_pos[Oil]] /= (1.0 - F[pu.phase_pos[Water]]);
                F[pu.phase_pos[Water]] = 0.0;
            }
        }

        if (active()[ Gas ]) {
            if (F[pu.phase_pos[Gas]] < 0.0) {
                if (active()[ Water ]) {
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
            if (active()[ Water ]) {
                F[pu.phase_pos[Water]] /= (1.0 - F[pu.phase_pos[Oil]]);
            }
            if (active()[ Gas ]) {
                F[pu.phase_pos[Gas]] /= (1.0 - F[pu.phase_pos[Oil]]);
            }
            if (has_solvent) {
                F_solvent /= (1.0 - F[pu.phase_pos[Oil]]);
            }
            F[pu.phase_pos[Oil]] = 0.0;
        }

        if (active()[ Water ]) {
            primary_variables_[WFrac] = F[pu.phase_pos[Water]];
        }
        if (active()[ Gas ]) {
            primary_variables_[GFrac] = F[pu.phase_pos[Gas]];
        }
        if(has_solvent) {
            primary_variables_[SFrac] = F_solvent;
        }

        // F_solvent is added to F_gas. This means that well_rate[Gas] also contains solvent.
        // More testing is needed to make sure this is correct for well groups and THP.
        if (has_solvent){
            F[pu.phase_pos[Gas]] += F_solvent;
        }

            // The interpretation of the first well variable depends on the well control
        const WellControls* wc = well_controls_;

        // TODO: we should only maintain one current control either from the well_state or from well_controls struct.
        // Either one can be more favored depending on the final strategy for the initilzation of the well control
        const int current = well_state.currentControls()[index_of_well_];
        const double target_rate = well_controls_iget_target(wc, current);

        for (int p = 0; p < np; ++p) {
            const double scal = scalingFactor(p);
            if (scal > 0) {
                F[p] /= scal ;
            } else {
                F[p] = 0.;
            }
        }

        switch (well_controls_iget_type(wc, current)) {
            case THP: // The BHP and THP both uses the total rate as first well variable.
            case BHP:
            {
                primary_variables_[XvarWell] = xvar_well_old[XvarWell] - dwells[0][XvarWell];

                switch (well_type_) {
                case INJECTOR:
                    for (int p = 0; p < np; ++p) {
                        const double comp_frac = comp_frac_[p];
                        well_state.wellRates()[index_of_well_ * np + p] = comp_frac * primary_variables_[XvarWell];
                    }
                    break;
                case PRODUCER:
                    for (int p = 0; p < np; ++p) {
                        well_state.wellRates()[index_of_well_ * np + p] = primary_variables_[XvarWell] * F[p];
                    }
                    break;
                }

                if (well_controls_iget_type(wc, current) == THP) {

                    // Calculate bhp from thp control and well rates
                    std::vector<double> rates(3, 0.0); // the vfp related only supports three phases for the moment

                    const Opm::PhaseUsage& pu = phaseUsage();
                    if (active()[ Water ]) {
                        rates[ Water ] = well_state.wellRates()[index_of_well_ * np + pu.phase_pos[ Water ] ];
                    }
                    if (active()[ Oil ]) {
                        rates[ Oil ]= well_state.wellRates()[index_of_well_ * np + pu.phase_pos[ Oil ] ];
                    }
                    if (active()[ Gas ]) {
                        rates[ Gas ]= well_state.wellRates()[index_of_well_ * np + pu.phase_pos[ Gas ] ];
                    }

                    well_state.bhp()[index_of_well_] = calculateBhpFromThp(rates, current);
                }
            }
                break;
            case SURFACE_RATE: // Both rate controls use bhp as first well variable
            case RESERVOIR_RATE:
            {
                const int sign1 = dwells[0][XvarWell] > 0 ? 1: -1;
                const double dx1_limited = sign1 * std::min(std::abs(dwells[0][XvarWell]),std::abs(xvar_well_old[XvarWell])*dBHPLimit);
                primary_variables_[XvarWell] = std::max(xvar_well_old[XvarWell] - dx1_limited,1e5);
                well_state.bhp()[index_of_well_] = primary_variables_[XvarWell];

                if (well_controls_iget_type(wc, current) == SURFACE_RATE) {
                    if (well_type_ == PRODUCER) {

                        const double* distr = well_controls_iget_distr(wc, current);

                        double F_target = 0.0;
                        for (int p = 0; p < np; ++p) {
                            F_target += distr[p] * F[p];
                        }
                        for (int p = 0; p < np; ++p) {
                            well_state.wellRates()[np * index_of_well_ + p] = F[p] * target_rate / F_target;
                        }
                    } else {

                        for (int p = 0; p < np; ++p) {
                            well_state.wellRates()[index_of_well_ * np + p] = comp_frac_[p] * target_rate;
                        }
                    }
                } else { // RESERVOIR_RATE
                    for (int p = 0; p < np; ++p) {
                        well_state.wellRates()[np * index_of_well_ + p] = F[p] * target_rate;
                    }
                }
            }
            break;
        } // end of switch (well_controls_iget_type(wc, current))

        // for the wells having a THP constaint, we should update their thp value
        // If it is under THP control, it will be set to be the target value. Otherwise,
        // the thp value will be calculated based on the bhp value, assuming the bhp value is correctly calculated.
        const int nwc = well_controls_get_num(wc);
        // Looping over all controls until we find a THP constraint
        int ctrl_index = 0;
        for ( ; ctrl_index < nwc; ++ctrl_index) {
            if (well_controls_iget_type(wc, ctrl_index) == THP) {
                // the current control
                const int current = well_state.currentControls()[index_of_well_];
                // If under THP control at the moment
                if (current == ctrl_index) {
                    const double thp_target = well_controls_iget_target(wc, current);
                    well_state.thp()[index_of_well_] = thp_target;
                } else { // otherwise we calculate the thp from the bhp value
                    const Opm::PhaseUsage& pu = phaseUsage();
                    std::vector<double> rates(3, 0.0);

                    if (active()[ Water ]) {
                        rates[ Water ] = well_state.wellRates()[index_of_well_*np + pu.phase_pos[ Water ] ];
                    }
                    if (active()[ Oil ]) {
                        rates[ Oil ] = well_state.wellRates()[index_of_well_*np + pu.phase_pos[ Oil ] ];
                    }
                    if (active()[ Gas ]) {
                        rates[ Gas ] = well_state.wellRates()[index_of_well_*np + pu.phase_pos[ Gas ] ];
                    }

                    const double bhp = well_state.bhp()[index_of_well_];

                    well_state.thp()[index_of_well_] = calculateThpFromBhp(rates, ctrl_index, bhp);
                }

                // the THP control is found, we leave the loop now
                break;
            }
        }  // end of for loop for seaching THP constraints

        // no THP constraint found
        if (ctrl_index == nwc) { // not finding a THP contstraints
            well_state.thp()[index_of_well_] = 0.0;
        }
    }





    template<typename TypeTag>
    void
    StandardWell<TypeTag>::
    updateWellStateWithTarget(const int current,
                              WellState& xw) const
    {
        // number of phases
        const int np = number_of_phases_;
        const int well_index = index_of_well_;
        const WellControls* wc = well_controls_;
        // Updating well state and primary variables.
        // Target values are used as initial conditions for BHP, THP, and SURFACE_RATE
        const double target = well_controls_iget_target(wc, current);
        const double* distr = well_controls_iget_distr(wc, current);
        switch (well_controls_iget_type(wc, current)) {
        case BHP:
            xw.bhp()[well_index] = target;
            // TODO: similar to the way below to handle THP
            // we should not something related to thp here when there is thp constraint
            break;

        case THP: {
            xw.thp()[well_index] = target;

            const Opm::PhaseUsage& pu = phaseUsage();
            std::vector<double> rates(3, 0.0);
            if (active()[ Water ]) {
                rates[ Water ] = xw.wellRates()[well_index*np + pu.phase_pos[ Water ] ];
            }
            if (active()[ Oil ]) {
                 rates[ Oil ] = xw.wellRates()[well_index*np + pu.phase_pos[ Oil ] ];
            }
            if (active()[ Gas ]) {
                rates[ Gas ] = xw.wellRates()[well_index*np + pu.phase_pos[ Gas ] ];
            }

            xw.bhp()[well_index] = calculateBhpFromThp(rates, current);
            break;
        }

        case RESERVOIR_RATE: // intentional fall-through
        case SURFACE_RATE:
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
                        xw.wellRates()[np*well_index + phase] = target / distr[phase];
                    } else {
                        xw.wellRates()[np * well_index + phase] = 0.;
                    }
                }
            } else if (well_type_ == PRODUCER) {
                // update the rates of phases under control based on the target,
                // and also update rates of phases not under control to keep the rate ratio,
                // assuming the mobility ratio does not change for the production wells
                double original_rates_under_phase_control = 0.0;
                for (int phase = 0; phase < np; ++phase) {
                    if (distr[phase] > 0.0) {
                        original_rates_under_phase_control += xw.wellRates()[np * well_index + phase] * distr[phase];
                    }
                }

                if (original_rates_under_phase_control != 0.0 ) {
                    double scaling_factor = target / original_rates_under_phase_control;

                    for (int phase = 0; phase < np; ++phase) {
                        xw.wellRates()[np * well_index + phase] *= scaling_factor;
                    }
                } else { // scaling factor is not well defied when original_rates_under_phase_control is zero
                    // separating targets equally between phases under control
                    const double target_rate_divided = target / numPhasesWithTargetsUnderThisControl;
                    for (int phase = 0; phase < np; ++phase) {
                        if (distr[phase] > 0.0) {
                            xw.wellRates()[np * well_index + phase] = target_rate_divided / distr[phase];
                        } else {
                            // this only happens for SURFACE_RATE control
                            xw.wellRates()[np * well_index + phase] = target_rate_divided;
                        }
                    }
                }
            } else {
                OPM_THROW(std::logic_error, "Expected PRODUCER or INJECTOR type of well");
            }

            break;
        } // end of switch

        updatePrimaryVariables(xw);
    }





    template<typename TypeTag>
    void
    StandardWell<TypeTag>::
    computePropertiesForWellConnectionPressures(const Simulator& ebosSimulator,
                                                const WellState& xw,
                                                std::vector<double>& b_perf,
                                                std::vector<double>& rsmax_perf,
                                                std::vector<double>& rvmax_perf,
                                                std::vector<double>& surf_dens_perf) const
    {
        const int nperf = number_of_perforations_;
        // TODO: can make this a member?
        const int numComp = numComponents();
        const PhaseUsage& pu = phaseUsage();
        b_perf.resize(nperf*numComp);
        surf_dens_perf.resize(nperf*numComp);
        const int w = index_of_well_;

        //rs and rv are only used if both oil and gas is present
        if (pu.phase_used[BlackoilPhases::Vapour] && pu.phase_used[BlackoilPhases::Liquid]) {
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
            const double p_above = perf == 0 ? xw.bhp()[w] : xw.perfPress()[first_perf_ + perf - 1];
            const double p_avg = (xw.perfPress()[first_perf_ + perf] + p_above)/2;
            const double temperature = fs.temperature(FluidSystem::oilPhaseIdx).value();

            if (pu.phase_used[BlackoilPhases::Aqua]) {
                b_perf[ pu.phase_pos[BlackoilPhases::Aqua] + perf * numComp] =
                FluidSystem::waterPvt().inverseFormationVolumeFactor(fs.pvtRegionIndex(), temperature, p_avg);
            }

            if (pu.phase_used[BlackoilPhases::Vapour]) {
                const int gaspos = pu.phase_pos[BlackoilPhases::Vapour] + perf * numComp;
                const int gaspos_well = pu.phase_pos[BlackoilPhases::Vapour] + w * pu.num_phases;

                if (pu.phase_used[BlackoilPhases::Liquid]) {
                    const int oilpos_well = pu.phase_pos[BlackoilPhases::Liquid] + w * pu.num_phases;
                    const double oilrate = std::abs(xw.wellRates()[oilpos_well]); //in order to handle negative rates in producers
                    rvmax_perf[perf] = FluidSystem::gasPvt().saturatedOilVaporizationFactor(fs.pvtRegionIndex(), temperature, p_avg);
                    if (oilrate > 0) {
                        const double gasrate = std::abs(xw.wellRates()[gaspos_well]) - xw.solventWellRate(w);
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

            if (pu.phase_used[BlackoilPhases::Liquid]) {
                const int oilpos = pu.phase_pos[BlackoilPhases::Liquid] + perf * numComp;
                const int oilpos_well = pu.phase_pos[BlackoilPhases::Liquid] + w * pu.num_phases;
                if (pu.phase_used[BlackoilPhases::Vapour]) {
                    rsmax_perf[perf] = FluidSystem::oilPvt().saturatedGasDissolutionFactor(fs.pvtRegionIndex(), temperature, p_avg);
                    const int gaspos_well = pu.phase_pos[BlackoilPhases::Vapour] + w * pu.num_phases;
                    const double gasrate = std::abs(xw.wellRates()[gaspos_well]) - xw.solventWellRate(w);
                    if (gasrate > 0) {
                        const double oilrate = std::abs(xw.wellRates()[oilpos_well]);
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
            for (int p = 0; p < pu.num_phases; ++p) {
                 surf_dens_perf[numComp*perf + p] = FluidSystem::referenceDensity( flowPhaseToEbosPhaseIdx( p ), fs.pvtRegionIndex());
            }

            // We use cell values for solvent injector
            if (has_solvent) {
                b_perf[numComp*perf + contiSolventEqIdx] = intQuants.solventInverseFormationVolumeFactor().value();
                surf_dens_perf[numComp*perf + contiSolventEqIdx] = intQuants.solventRefDensity();
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
        const int num_comp = numComponents();
        const PhaseUsage& phase_usage = phaseUsage();

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
        const int gaspos = phase_usage.phase_pos[BlackoilPhases::Vapour];
        const int oilpos = phase_usage.phase_pos[BlackoilPhases::Liquid];
        std::vector<double> mix(num_comp,0.0);
        std::vector<double> x(num_comp);
        std::vector<double> surf_dens(num_comp);
        std::vector<double> dens(nperf);

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
                for (int phase = 0; phase < np; ++phase) {
                    mix[phase] = comp_frac_[phase];
                }
                // intialize 0.0 for comIdx >= np;
            }
            // Compute volume ratio.
            x = mix;
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
    typename StandardWell<TypeTag>::ConvergenceReport
    StandardWell<TypeTag>::
    getWellConvergence(const std::vector<double>& B_avg) const
    {
        typedef double Scalar;
        typedef std::vector< Scalar > Vector;

        const int np = number_of_phases_;
        const int numComp = numComponents();

        // the following implementation assume that the polymer is always after the w-o-g phases
        // For the polymer case, there is one more mass balance equations of reservoir than wells
        assert((int(B_avg.size()) == numComp) || has_polymer);

        const double tol_wells = param_.tolerance_wells_;
        const double maxResidualAllowed = param_.max_residual_allowed_;

        // TODO: it should be the number of numWellEq
        // using numComp here for flow_ebos running 2p case.
        std::vector<Scalar> res(numComp);
        for (int comp = 0; comp < numComp; ++comp) {
            // magnitude of the residual matters
            res[comp] = std::abs(resWell_[0][comp]);
        }

        Vector well_flux_residual(numComp);

        // Finish computation
        for ( int compIdx = 0; compIdx < numComp; ++compIdx )
        {
            well_flux_residual[compIdx] = B_avg[compIdx] * res[compIdx];
        }

        ConvergenceReport report;
        // checking if any NaN or too large residuals found
        // TODO: not understand why phase here while component in other places.
        for (int phaseIdx = 0; phaseIdx < np; ++phaseIdx) {
            const auto& phaseName = FluidSystem::phaseName(flowPhaseToEbosPhaseIdx(phaseIdx));

            if (std::isnan(well_flux_residual[phaseIdx])) {
                report.nan_residual_found = true;
                const typename ConvergenceReport::ProblemWell problem_well = {name(), phaseName};
                report.nan_residual_wells.push_back(problem_well);
            } else {
                if (well_flux_residual[phaseIdx] > maxResidualAllowed) {
                    report.too_large_residual_found = true;
                    const typename ConvergenceReport::ProblemWell problem_well = {name(), phaseName};
                    report.too_large_residual_wells.push_back(problem_well);
                }
            }
        }

        if ( !(report.nan_residual_found || report.too_large_residual_found) ) { // no abnormal residual value found
            // check convergence
            for ( int compIdx = 0; compIdx < numComp; ++compIdx )
            {
                report.converged = report.converged && (well_flux_residual[compIdx] < tol_wells);
            }
        } else { // abnormal values found and no need to check the convergence
            report.converged = false;
        }

        return report;
    }





    template<typename TypeTag>
    void
    StandardWell<TypeTag>::
    computeWellConnectionDensitesPressures(const WellState& xw,
                                           const std::vector<double>& b_perf,
                                           const std::vector<double>& rsmax_perf,
                                           const std::vector<double>& rvmax_perf,
                                           const std::vector<double>& surf_dens_perf)
    {
        // Compute densities
        const int nperf = number_of_perforations_;
        const int numComponent = numComponents();
        const int np = number_of_phases_;
        std::vector<double> perfRates(b_perf.size(),0.0);

        for (int perf = 0; perf < nperf; ++perf) {
            for (int phase = 0; phase < np; ++phase) {
                perfRates[perf*numComponent + phase] =  xw.perfPhaseRates()[(first_perf_ + perf) * np + phase];
            }
            if(has_solvent) {
                perfRates[perf*numComponent + contiSolventEqIdx] =  xw.perfRateSolvent()[first_perf_ + perf];
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
    solveEqAndUpdateWellState(WellState& well_state)
    {
        // We assemble the well equations, then we check the convergence,
        // which is why we do not put the assembleWellEq here.
        BVectorWell dx_well(1);
        invDuneD_.mv(resWell_, dx_well);

        updateWellState(dx_well, well_state);
    }





    template<typename TypeTag>
    void
    StandardWell<TypeTag>::
    calculateExplicitQuantities(const Simulator& ebosSimulator,
                                const WellState& well_state)
    {
        computeWellConnectionPressures(ebosSimulator, well_state);
        computeAccumWell();
    }





    template<typename TypeTag>
    void
    StandardWell<TypeTag>::
    computeAccumWell()
    {
        // TODO: it should be num_comp, while it also bring problem for
        // the polymer case.
        for (int eq_idx = 0; eq_idx < numWellEq; ++eq_idx) {
            F0_[eq_idx] = wellSurfaceVolumeFraction(eq_idx).value();
        }
    }





    template<typename TypeTag>
    void
    StandardWell<TypeTag>::
    apply(const BVector& x, BVector& Ax) const
    {
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
                                          WellState& well_state) const
    {
        BVectorWell xw(1);
        recoverSolutionWell(x, xw);
        updateWellState(xw, well_state);
    }





    template<typename TypeTag>
    void
    StandardWell<TypeTag>::
    computeWellRatesWithBhp(const Simulator& ebosSimulator,
                            const EvalWell& bhp,
                            std::vector<double>& well_flux) const
    {
        const int np = number_of_phases_;
        const int numComp = numComponents();
        well_flux.resize(np, 0.0);

        const bool allow_cf = crossFlowAllowed(ebosSimulator);

        for (int perf = 0; perf < number_of_perforations_; ++perf) {
            const int cell_idx = well_cells_[perf];
            const auto& intQuants = *(ebosSimulator.model().cachedIntensiveQuantities(cell_idx, /*timeIdx=*/ 0));
            // flux for each perforation
            std::vector<EvalWell> cq_s(numComp, 0.0);
            std::vector<EvalWell> mob(numComp, 0.0);
            getMobility(ebosSimulator, perf, mob);
            computePerfRate(intQuants, mob, well_index_[perf], bhp, perf_pressure_diffs_[perf], allow_cf, cq_s);

            for(int p = 0; p < np; ++p) {
                well_flux[p] += cq_s[p].value();
            }
        }
    }





    template<typename TypeTag>
    std::vector<double>
    StandardWell<TypeTag>::
    computeWellPotentialWithTHP(const Simulator& ebosSimulator,
                                const double initial_bhp, // bhp from BHP constraints
                                const std::vector<double>& initial_potential) const
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
                    if (active()[ Water ]) {
                        rates[ Water ] = potentials[pu.phase_pos[ Water ] ];
                    }
                    if (active()[ Oil ]) {
                        rates[ Oil ] = potentials[pu.phase_pos[ Oil ] ];
                    }
                    if (active()[ Gas ]) {
                        rates[ Gas ] = potentials[pu.phase_pos[ Gas ] ];
                    }

                    const double bhp_calculated = calculateBhpFromThp(rates, ctrl_index);

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
                OPM_THROW(std::runtime_error, "Unvalid bhp value obtained during the potential calculation for well " << name());
            }

            converged = std::abs(old_bhp - bhp) < bhp_tolerance;

            computeWellRatesWithBhp(ebosSimulator, bhp, potentials);

            // checking whether the potentials have valid values
            for (const double value : potentials) {
                if (std::isinf(value) || std::isnan(value)) {
                    OPM_THROW(std::runtime_error, "Unvalid potential value obtained during the potential calculation for well " << name());
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
            OPM_THROW(std::runtime_error, "Failed in getting converged for the potential calculation for well " << name());
        }

        return potentials;
    }





    template<typename TypeTag>
    void
    StandardWell<TypeTag>::
    computeWellPotentials(const Simulator& ebosSimulator,
                          const WellState& well_state,
                          std::vector<double>& well_potentials) // const
    {
        updatePrimaryVariables(well_state);
        computeWellConnectionPressures(ebosSimulator, well_state);

        // initialize the primary variables in Evaluation, which is used in computePerfRate for computeWellPotentials
        // TODO: for computeWellPotentials, no derivative is required actually
        initPrimaryVariablesEvaluation();

        const int np = number_of_phases_;
        well_potentials.resize(np, 0.0);

        // get the bhp value based on the bhp constraints
        const double bhp = mostStrictBhpFromBhpLimits();

        // does the well have a THP related constraint?
        if ( !wellHasTHPConstraints() ) {
            assert(std::abs(bhp) != std::numeric_limits<double>::max());

            computeWellRatesWithBhp(ebosSimulator, bhp, well_potentials);
        } else {
            // the well has a THP related constraint
            // checking whether a well is newly added, it only happens at the beginning of the report step
            if ( !well_state.isNewWell(index_of_well_) ) {
                for (int p = 0; p < np; ++p) {
                    // This is dangerous for new added well
                    // since we are not handling the initialization correctly for now
                    well_potentials[p] = well_state.wellRates()[index_of_well_ * np + p];
                }
            } else {
                // We need to generate a reasonable rates to start the iteration process
                computeWellRatesWithBhp(ebosSimulator, bhp, well_potentials);
                for (double& value : well_potentials) {
                    // make the value a little safer in case the BHP limits are default ones
                    // TODO: a better way should be a better rescaling based on the investigation of the VFP table.
                    const double rate_safety_scaling_factor = 0.00001;
                    value *= rate_safety_scaling_factor;
                }
            }

            well_potentials = computeWellPotentialWithTHP(ebosSimulator, bhp, well_potentials);
        }
    }





    template<typename TypeTag>
    void
    StandardWell<TypeTag>::
    updatePrimaryVariables(const WellState& well_state) const
    {
        const int np = number_of_phases_;
        const int well_index = index_of_well_;
        const WellControls* wc = well_controls_;
        const double* distr = well_controls_get_current_distr(wc);
        const auto pu = phaseUsage();

        switch (well_controls_get_current_type(wc)) {
        case THP:
        case BHP: {
            primary_variables_[XvarWell] = 0.0;
            if (well_type_ == INJECTOR) {
                for (int p = 0; p < np; ++p) {
                    primary_variables_[XvarWell] += well_state.wellRates()[np*well_index + p] * comp_frac_[p];
                }
            } else {
                for (int p = 0; p < np; ++p) {
                    primary_variables_[XvarWell] += scalingFactor(p) * well_state.wellRates()[np*well_index + p];
                }
            }
            break;
        }
        case RESERVOIR_RATE: // Intentional fall-through
        case SURFACE_RATE:
            primary_variables_[XvarWell] = well_state.bhp()[well_index];
            break;
        } // end of switch

        double tot_well_rate = 0.0;
        for (int p = 0; p < np; ++p)  {
            tot_well_rate += scalingFactor(p) * well_state.wellRates()[np*well_index + p];
        }
        if(std::abs(tot_well_rate) > 0) {
            if (active()[ Water ]) {
                primary_variables_[WFrac] = scalingFactor(pu.phase_pos[Water]) * well_state.wellRates()[np*well_index + pu.phase_pos[Water]] / tot_well_rate;
            }
            if (active()[ Gas ]) {
                primary_variables_[GFrac] = scalingFactor(pu.phase_pos[Gas]) * (well_state.wellRates()[np*well_index + pu.phase_pos[Gas]] - well_state.solventWellRate(well_index)) / tot_well_rate ;
            }
            if (has_solvent) {
                primary_variables_[SFrac] = scalingFactor(pu.phase_pos[Gas]) * well_state.solventWellRate(well_index) / tot_well_rate ;
            }
        } else { // tot_well_rate == 0
            if (well_type_ == INJECTOR) {
                // only single phase injection handled
                if (active()[Water]) {
                    if (distr[Water] > 0.0) {
                        primary_variables_[WFrac] = 1.0;
                    } else {
                        primary_variables_[WFrac] = 0.0;
                    }
                }

                if (active()[Gas]) {
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
                if (active()[Water]) {
                    primary_variables_[WFrac] = 1.0 / np;
                }
                if (active()[Gas]) {
                    primary_variables_[GFrac] = 1.0 / np;
                }
            } else {
                OPM_THROW(std::logic_error, "Expected PRODUCER or INJECTOR type of well");
            }
        }
    }





    template<typename TypeTag>
    template<class ValueType>
    ValueType
    StandardWell<TypeTag>::
    calculateBhpFromThp(const std::vector<ValueType>& rates,
                        const int control_index) const
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
             OPM_THROW(std::logic_error, "Expected INJECTOR or PRODUCER well");
         }

         return bhp;
    }





    template<typename TypeTag>
    double
    StandardWell<TypeTag>::
    calculateThpFromBhp(const std::vector<double>& rates,
                        const int control_index,
                        const double bhp) const
    {
        assert(int(rates.size()) == 3); // the vfp related only supports three phases now.

        const double aqua = rates[Water];
        const double liquid = rates[Oil];
        const double vapour = rates[Gas];

        const int vfp        = well_controls_iget_vfp(well_controls_, control_index);
        const double& alq    = well_controls_iget_alq(well_controls_, control_index);

        // pick the density in the top layer
        const double rho = perf_densities_[0];

        double thp = 0.0;
        if (well_type_ == INJECTOR) {
            const double vfp_ref_depth = vfp_properties_->getInj()->getTable(vfp)->getDatumDepth();

            const double dp = wellhelpers::computeHydrostaticCorrection(ref_depth_, vfp_ref_depth, rho, gravity_);

            thp = vfp_properties_->getInj()->thp(vfp, aqua, liquid, vapour, bhp + dp);
         }
         else if (well_type_ == PRODUCER) {
             const double vfp_ref_depth = vfp_properties_->getProd()->getTable(vfp)->getDatumDepth();

             const double dp = wellhelpers::computeHydrostaticCorrection(ref_depth_, vfp_ref_depth, rho, gravity_);

             thp = vfp_properties_->getProd()->thp(vfp, aqua, liquid, vapour, bhp + dp, alq);
         }
         else {
             OPM_THROW(std::logic_error, "Expected INJECTOR or PRODUCER well");
         }

         return thp;
    }

    template<typename TypeTag>
    double
    StandardWell<TypeTag>::scalingFactor(const int phaseIdx) const
    {
        const WellControls* wc = well_controls_;
        const double* distr = well_controls_get_current_distr(wc);

        if (well_controls_get_current_type(wc) == RESERVOIR_RATE) {
            if (has_solvent && phaseIdx == contiSolventEqIdx )
                   OPM_THROW(std::runtime_error, "RESERVOIR_RATE control in combination with solvent is not implemented");

            return distr[phaseIdx];
        }
        const auto& pu = phaseUsage();
        if (active()[Water] && pu.phase_pos[Water] == phaseIdx)
            return 1.0;
        if (active()[Oil] && pu.phase_pos[Oil] == phaseIdx)
            return 1.0;
        if (active()[Gas] && pu.phase_pos[Gas] == phaseIdx)
            return 0.01;
        if (has_solvent && phaseIdx == contiSolventEqIdx )
            return 0.01;

        // we should not come this far
        assert(false);
        return 1.0;
    }


}
