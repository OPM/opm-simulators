/*
  Copyright 2017 SINTEF ICT, Applied Mathematics.
  Copyright 2017 Statoil ASA.

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
    StandardWell(const Well* well, const int time_step, const Wells* wells)
    : WellInterface<TypeTag>(well, time_step, wells)
    , perf_densities_(numberOfPerforations())
    , perf_pressure_diffs_(numberOfPerforations())
    , well_variables_(numWellEq) // the number of the primary variables
    {
        duneB_.setBuildMode( Mat::row_wise );
        duneC_.setBuildMode( Mat::row_wise );
        invDuneD_.setBuildMode( Mat::row_wise );
    }





    template<typename TypeTag>
    void
    StandardWell<TypeTag>::
    init(const PhaseUsage* phase_usage_arg,
         const std::vector<bool>* active_arg,
         const VFPProperties* vfp_properties_arg,
         const double gravity_arg,
         const int num_cells)
    {
        WellInterface<TypeTag>(phase_usage_arg, active_arg,
                               vfp_properties_arg, gravity_arg, num_cells);

        // setup sparsity pattern for the matrices
        // TODO: C and B are opposite compared with the notations used in the paper.
        //[A B^T    [x    =  [ res
        // C D] x_well]      res_well]
        // set the size of the matrices
        invDuneD_.setSize(1, 1, 1);
        duneC_.setSize(1, num_cells, numberOfPerforations());
        duneB_.setSize(1, num_cells, numberOfPerforations());

        for (auto row=invDuneD_.createbegin(), end = invDuneD_.createend(); row!=end; ++row) {
            // Add nonzeros for diagonal
            row.insert(row.index());
        }

        for (auto row = duneC_.createbegin(), end = duneC_.createend(); row!=end; ++row) {
            // Add nonzeros for diagonal
            for (int perf = 0 ; perf < numberOfPerforations(); ++perf) {
                const int cell_idx = wellCells()[perf];
                row.insert(cell_idx);
            }
        }

        // make the B^T matrix
        for (auto row = duneB_.createbegin(), end = duneB_.createend(); row!=end; ++row) {
            for (int perf = 0; perf < numberOfPerforations(); ++perf) {
                const int cell_idx = wellCells()[perf];
                row.insert(cell_idx);
            }
        }

        resWell_.resize(1);

        // resize temporary class variables
        Cx_.resize( duneC_.N() );
        invDrw_.resize( invDuneD_.N() );
    }





    template<typename TypeTag>
    const std::vector<double>&
    StandardWell<TypeTag>::
    perfDensities() const
    {
        return perf_densities_;
    }





    template<typename TypeTag>
    std::vector<double>&
    StandardWell<TypeTag>::
    perfDensities()
    {
        return perf_densities_;
    }





    template<typename TypeTag>
    const std::vector<double>&
    StandardWell<TypeTag>::
    perfPressureDiffs() const
    {
        return perf_pressure_diffs_;
    }





    template<typename TypeTag>
    std::vector<double>&
    StandardWell<TypeTag>::
    perfPressureDiffs()
    {
        return perf_pressure_diffs_;
    }





    template<typename TypeTag>
    void StandardWell<TypeTag>::
    setWellVariables(const WellState& well_state)
    {
        const int nw = well_state.bhp().size();
        const int numComp = numComponents();
        for (int eqIdx = 0; eqIdx < numComp; ++eqIdx) {
            const unsigned int idx = nw * eqIdx + indexOfWell();
            assert( eqIdx < well_variables_.size() );
            assert( idx < well_state.wellSolutions().size() );

            well_variables_[eqIdx] = 0.0;
            well_variables_[eqIdx].setValue(well_state.wellSolutions()[idx]);
            well_variables_[eqIdx].setDerivative(numEq + eqIdx, 1.0);
        }
    }





    template<typename TypeTag>
    typename StandardWell<TypeTag>::EvalWell
    StandardWell<TypeTag>::
    getBhp() const
    {
        const WellControls* wc = wellControls();
        if (well_controls_get_current_type(wc) == BHP) {
            EvalWell bhp = 0.0;
            const double target_rate = well_controls_get_current_target(wc);
            bhp.setValue(target_rate);
            return bhp;
        } else if (well_controls_get_current_type(wc) == THP) {
            const int control = well_controls_get_current(wc);
            const double thp = well_controls_get_current_target(wc);
            const double alq = well_controls_iget_alq(wc, control);
            const int table_id = well_controls_iget_vfp(wc, control);
            EvalWell aqua = 0.0;
            EvalWell liquid = 0.0;
            EvalWell vapour = 0.0;
            EvalWell bhp = 0.0;
            double vfp_ref_depth = 0.0;

            const Opm::PhaseUsage& pu = phaseUsage();

            if (active()[ Water ]) {
                aqua = getQs(pu.phase_pos[ Water]);
            }
            if (active()[ Oil ]) {
                liquid = getQs(pu.phase_pos[ Oil ]);
            }
            if (active()[ Gas ]) {
                vapour = getQs(pu.phase_pos[ Gas ]);
            }
            if (wellType() == INJECTOR) {
                bhp = vfp_properties_->getInj()->bhp(table_id, aqua, liquid, vapour, thp);
                vfp_ref_depth = vfp_properties_->getInj()->getTable(table_id)->getDatumDepth();
            } else {
                bhp = vfp_properties_->getProd()->bhp(table_id, aqua, liquid, vapour, thp, alq);
                vfp_ref_depth = vfp_properties_->getProd()->getTable(table_id)->getDatumDepth();
            }

            // pick the density in the top layer
            const double rho = perf_densities_[0];
            // TODO: not sure whether it is always correct
            const double well_ref_depth = perfDepth()[0];
            const double dp = wellhelpers::computeHydrostaticCorrection(well_ref_depth, vfp_ref_depth, rho, gravity_);
            bhp -= dp;
            return bhp;
        }

        return well_variables_[XvarWell];
    }





    template<typename TypeTag>
    typename StandardWell<TypeTag>::EvalWell
    StandardWell<TypeTag>::
    getQs(const int phase) const
    {
        EvalWell qs = 0.0;

        const WellControls* wc = wellControls();
        const int np = numberOfPhases();
        const double target_rate = well_controls_get_current_target(wc);

        // TODO: we need to introduce numComponents() for StandardWell
        // assert(phase < numComponents());
        const auto pu = phaseUsage();

        // TODO: the formulation for the injectors decides it only work with single phase
        // surface rate injection control. Improvement will be required.
        if (wellType() == INJECTOR) {
            // TODO: adding the handling related to solvent
            /* if (has_solvent_ ) {
                // TODO: investigate whether the use of the comp_frac is justified.
                double comp_frac = 0.0;
                if (compIdx == solventCompIdx) { // solvent
                    comp_frac = wells().comp_frac[np*wellIdx + pu.phase_pos[ Gas ]] * wsolvent(wellIdx);
                } else if (compIdx == pu.phase_pos[ Gas ]) {
                    comp_frac = wells().comp_frac[np*wellIdx + compIdx] * (1.0 - wsolvent(wellIdx));
                } else {
                    comp_frac = wells().comp_frac[np*wellIdx + compIdx];
                }
                if (comp_frac == 0.0) {
                    return qs; //zero
                }

                if (well_controls_get_current_type(wc) == BHP || well_controls_get_current_type(wc) == THP) {
                    return comp_frac * well_variables_[nw*XvarWell + wellIdx];
                }

                qs.setValue(comp_frac * target_rate);
                return qs;
            } */
            const double comp_frac = compFrac()[phase];
            if (comp_frac == 0.0) {
                return qs;
            }

            if (well_controls_get_current_type(wc) == BHP || well_controls_get_current_type(wc) == THP) {
                return well_variables_[XvarWell];
            }
            qs.setValue(target_rate);
            return qs;
        }

        // Producers
        if (well_controls_get_current_type(wc) == BHP || well_controls_get_current_type(wc) == THP ) {
            return well_variables_[XvarWell] * wellVolumeFractionScaled(phase);
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
                // TODO: handling solvent related later
                /* if (has_solvent_ && phase_under_control == Gas) {
                    // for GRAT controlled wells solvent is included in the target
                    wellVolumeFractionScaledPhaseUnderControl += wellVolumeFractionScaled(solventCompIdx);
                } */

                if (phase == phase_under_control) {
                    /* if (has_solvent_ && phase_under_control == Gas) {
                        qs.setValue(target_rate * wellVolumeFractionScaled(Gas).value() / wellVolumeFractionScaledPhaseUnderControl.value() );
                        return qs;
                    } */
                    qs.setValue(target_rate);
                    return qs;
                }

                // TODO: not sure why the single phase under control will have near zero fraction
                const double eps = 1e-6;
                if (wellVolumeFractionScaledPhaseUnderControl < eps) {
                    return qs;
                }
                return (target_rate * wellVolumeFractionScaled(phase) / wellVolumeFractionScaledPhaseUnderControl);
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
                return (target_rate * wellVolumeFractionScaled(phase) / combined_volume_fraction);
            }

            // TODO: three phase surface rate control is not tested yet
            if (num_phases_under_rate_control == 3) {
                return target_rate * wellSurfaceVolumeFraction(phase);
            }
        } else if (well_controls_get_current_type(wc) == RESERVOIR_RATE) {
            // ReservoirRate
            return target_rate * wellVolumeFractionScaled(phase);
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
        // TODO: we should be able to set the g for the well based on the control type
        // instead of using explicit code for g all the times
        const WellControls* wc = wellControls();
        if (well_controls_get_current_type(wc) == RESERVOIR_RATE) {

            if (has_solvent && compIdx == solventCompIdx) {
                return wellVolumeFraction(compIdx);
            }
            const double* distr = well_controls_get_current_distr(wc);
            assert(compIdx < 3);
            if (distr[compIdx] > 0.) {
                return wellVolumeFraction(compIdx) / distr[compIdx];
            } else {
                // TODO: not sure why return EvalWell(0.) causing problem here
                // Probably due to the wrong Jacobians.
                return wellVolumeFraction(compIdx);
            }
        }

        std::vector<double> g = {1, 1, 0.01, 0.01};
        return (wellVolumeFraction(compIdx) / g[compIdx]);
    }





    template<typename TypeTag>
    typename StandardWell<TypeTag>::EvalWell
    StandardWell<TypeTag>::
    wellVolumeFraction(const int compIdx) const
    {
        if (compIdx == Water) {
            return well_variables_[WFrac];
        }

        if (compIdx == Gas) {
            return well_variables_[GFrac];
        }

        if (compIdx == solventCompIdx) {
            return well_variables_[SFrac];
        }

        // Oil fraction
        EvalWell well_fraction = 1.0;
        if (active()[Water]) {
            well_fraction -= well_variables_[WFrac];
        }

        if (active()[Gas]) {
            well_fraction -= well_variables_[GFrac];
        }
        if (has_solvent) {
            well_fraction -= well_variables_[SFrac];
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
            out.setDerivative(eqIdx, in.derivative(flowToEbosPvIdx(eqIdx)));
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
        const int np = numPhases();
        const int numComp = numComponents();
        std::vector<EvalWell> cmix_s(numComp,0.0);
        for (int componentIdx = 0; componentIdx < numComp; ++componentIdx) {
            cmix_s[componentIdx] = wellSurfaceVolumeFraction(componentIdx);
        }
        auto& fs = intQuants.fluidState();

        EvalWell pressure = extendEval(fs.pressure(FluidSystem::oilPhaseIdx));
        EvalWell rs = extendEval(fs.Rs());
        EvalWell rv = extendEval(fs.Rv());
        std::vector<EvalWell> b_perfcells_dense(numComp, 0.0);
        for (int phase = 0; phase < np; ++phase) {
            int ebosPhaseIdx = flowPhaseToEbosPhaseIdx(phase);
            b_perfcells_dense[phase] = extendEval(fs.invB(ebosPhaseIdx));
        }
        if (has_solvent) {
            b_perfcells_dense[solventCompIdx] = extendEval(intQuants.solventInverseFormationVolumeFactor());
        }

        // Pressure drawdown (also used to determine direction of flow)
        EvalWell well_pressure = bhp + cdp;
        EvalWell drawdown = pressure - well_pressure;

        // producing perforations
        if ( drawdown.value() > 0 )  {
            //Do nothing if crossflow is not allowed
            if (!allow_cf && wellType() == INJECTOR) {
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
            if (!allow_cf && wellType() == PRODUCER) {
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
                volumeRatio += cmix_s[solventCompIdx] / b_perfcells_dense[solventCompIdx];
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
        // TODO: accessing well_state information is the only place to use nw at the moment
        const int nw = well_state.bhp().size();
        const int numComp = numComponents();
        const int np = numPhases();

        // clear all entries
        duneB_ = 0.0;
        duneC_ = 0.0;
        invDuneD_ = 0.0;
        resWell_ = 0.0;

        auto& ebosJac = ebosSimulator.model().linearizer().matrix();
        auto& ebosResid = ebosSimulator.model().linearizer().residual();

        // TODO: it probably can be static member for StandardWell
        const double volume = 0.002831684659200; // 0.1 cu ft;

        const bool allow_cf = allow_cross_flow(ebosSimulator);

        const EvalWell& bhp = getBhp();

        for (int perf = 0; perf < numberOfPerforations(); ++perf) {

            const int cell_idx = wellCells()[perf];
            const auto& intQuants = *(ebosSimulator.model().cachedIntensiveQuantities(cell_idx, /*timeIdx=*/ 0));
            std::vector<EvalWell> cq_s(numComp,0.0);
            std::vector<EvalWell> mob(numComp, 0.0);
            getMobility(ebosSimulator, perf, mob);
            computePerfRate(intQuants, mob, wellIndex()[perf], bhp, perfPressureDiffs()[perf], allow_cf, cq_s);

                for (int componentIdx = 0; componentIdx < numComp; ++componentIdx) {
                    // the cq_s entering mass balance equations need to consider the efficiency factors.
                    const EvalWell cq_s_effective = cq_s[componentIdx] * well_efficiency_factor_;

                    if (!only_wells) {
                        // subtract sum of component fluxes in the reservoir equation.
                        // need to consider the efficiency factor
                        ebosResid[cell_idx][flowPhaseToEbosCompIdx(componentIdx)] -= cq_s_effective.value();
                    }

                    // subtract sum of phase fluxes in the well equations.
                    resWell_[0][componentIdx] -= cq_s[componentIdx].value();

                    // assemble the jacobians
                    for (int pvIdx = 0; pvIdx < numWellEq; ++pvIdx) {
                        if (!only_wells) {
                            // also need to consider the efficiency factor when manipulating the jacobians.
                            ebosJac[cell_idx][cell_idx][flowPhaseToEbosCompIdx(componentIdx)][flowToEbosPvIdx(pvIdx)] -= cq_s_effective.derivative(pvIdx);
                            duneB_[0][cell_idx][pvIdx][flowPhaseToEbosCompIdx(componentIdx)] -= cq_s_effective.derivative(pvIdx+numEq); // intput in transformed matrix
                            duneC_[0][cell_idx][componentIdx][flowToEbosPvIdx(pvIdx)] -= cq_s_effective.derivative(pvIdx);
                        }
                        invDuneD_[0][0][componentIdx][pvIdx] -= cq_s[componentIdx].derivative(pvIdx+numEq);
                    }

                    // add trivial equation for 2p cases (Only support water + oil)
                    if (numComp == 2) {
                        assert(!active()[ Gas ]);
                        invDuneD_[0][0][Gas][Gas] = 1.0;
                    }

                    // Store the perforation phase flux for later usage.
                    if (componentIdx == solventCompIdx) {// if (flowPhaseToEbosCompIdx(componentIdx) == Solvent)
                        well_state.perfRateSolvent()[perf] = cq_s[componentIdx].value();
                    } else {
                        well_state.perfPhaseRates()[perf*np + componentIdx] = cq_s[componentIdx].value();
                    }
                }

                // Store the perforation pressure for later usage.
                well_state.perfPress()[perf] = well_state.bhp()[indexOfWell()] + perfPressureDiffs()[perf];
            }

            // add vol * dF/dt + Q to the well equations;
            for (int componentIdx = 0; componentIdx < numComp; ++componentIdx) {
                // TODO: the F0_ here is not initialized yet here, which should happen in the first iteration, so it should happen in the assemble function
                EvalWell resWell_loc = (wellSurfaceVolumeFraction(componentIdx) - F0_[componentIdx]) * volume / dt;
                resWell_loc += getQs(componentIdx);
                for (int pvIdx = 0; pvIdx < numWellEq; ++pvIdx) {
                    invDuneD_[0][0][componentIdx][pvIdx] += resWell_loc.derivative(pvIdx+numEq);
                }
                resWell_[0][componentIdx] += resWell_loc.value();
            }

        // do the local inversion of D.
        localInvert( invDuneD_ );
    }





    template<typename TypeTag>
    bool
    StandardWell<TypeTag>::
    allow_cross_flow(const Simulator& ebosSimulator) const
    {
        if (allowCrossFlow()) {
            return true;
        }

        // TODO: investigate the justification of the following situation

        // check for special case where all perforations have cross flow
        // then the wells must allow for cross flow
        for (int perf = 0; perf < numberOfPerforations(); ++perf) {
            const int cell_idx = wellCells()[perf];
            const auto& intQuants = *(ebosSimulator.model().cachedIntensiveQuantities(cell_idx, /*timeIdx=*/0));
            const auto& fs = intQuants.fluidState();
            EvalWell pressure = extendEval(fs.pressure(FluidSystem::oilPhaseIdx));
            EvalWell bhp = getBhp();

            // Pressure drawdown (also used to determine direction of flow)
            EvalWell well_pressure = bhp + perfPressureDiffs()[perf];
            EvalWell drawdown = pressure - well_pressure;

            if (drawdown.value() < 0 && wellType() == INJECTOR)  {
                return false;
            }

            if (drawdown.value() > 0 && wellType() == PRODUCER)  {
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
        const int np = numberOfPhases();
        const int cell_idx = wellCells()[perf];
        assert (int(mob.size()) == numComponents());
        const auto& intQuants = *(ebosSimulator.model().cachedIntensiveQuantities(cell_idx, /*timeIdx=*/0));
        const auto& materialLawManager = ebosSimulator.problem().materialLawManager();

        // either use mobility of the perforation cell or calcualte its own
        // based on passing the saturation table index
        const int satid = saturationTableNumber()[perf] - 1;
        const int satid_elem = materialLawManager->satnumRegionIdx(cell_idx);
        if( satid == satid_elem ) { // the same saturation number is used. i.e. just use the mobilty from the cell

            for (int phase = 0; phase < np; ++phase) {
                int ebosPhaseIdx = flowPhaseToEbosPhaseIdx(phase);
                mob[phase] = extendEval(intQuants.mobility(ebosPhaseIdx));
            }
            if (has_solvent) {
                mob[solventCompIdx] = extendEval(intQuants.solventMobility());
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
    }





    template<typename TypeTag>
    void
    StandardWell<TypeTag>::
    updateWellState(const BVector& dwells,
                    const BlackoilModelParameters& param,
                    WellState& well_state) const
    {
        const int np = numberOfPhases();
        const int nw = well_state.bhp().size();
        const double dFLimit = param.dbhp_max_rel_;
        const double dBHPLimit = param.dwell_fraction_max_;

        std::vector<double> xvar_well_old(numWellEq);
        // TODO: better way to handle this?
        for (int i = 0; i < numWellEq; ++i) {
            xvar_well_old[i] = well_state.wellSolutions()[i * nw + indexOfWell()];
        }

        // update the second and third well variable (The flux fractions)
        std::vector<double> F(np,0.0);
        if (active()[ Water ]) {
            const int sign2 = dwells[0][WFrac] > 0 ? 1: -1;
            const double dx2_limited = sign2 * std::min(std::abs(dwells[0][WFrac]),dFLimit);
            well_state.wellSolutions()[WFrac * nw + indexOfWell()] = xvar_well_old[WFrac] - dx2_limited;
        }

        if (active()[ Gas ]) {
            const int sign3 = dwells[0][GFrac] > 0 ? 1: -1;
            const double dx3_limited = sign3 * std::min(std::abs(dwells[0][GFrac]),dFLimit);
            well_state.wellSolutions()[GFrac*nw + indexOfWell()] = xvar_well_old[GFrac] - dx3_limited;
        }

        if (has_solvent) {
            const int sign4 = dwells[0][SFrac] > 0 ? 1: -1;
            const double dx4_limited = sign4 * std::min(std::abs(dwells[0][SFrac]),dFLimit);
            well_state.wellSolutions()[SFrac*nw + indexOfWell()] = xvar_well_old[SFrac] - dx4_limited;
        }

        assert(active()[ Oil ]);
        F[Oil] = 1.0;

        if (active()[ Water ]) {
            F[Water] = well_state.wellSolutions()[WFrac*nw + indexOfWell()];
            F[Oil] -= F[Water];
        }

        if (active()[ Gas ]) {
            F[Gas] = well_state.wellSolutions()[GFrac*nw + indexOfWell()];
            F[Oil] -= F[Gas];
        }

        double F_solvent = 0.0;
        if (has_solvent) {
            F_solvent = well_state.wellSolutions()[SFrac*nw + indexOfWell()];
            F[Oil] -= F_solvent;
        }

        if (active()[ Water ]) {
            if (F[Water] < 0.0) {
                if (active()[ Gas ]) {
                        F[Gas] /= (1.0 - F[Water]);
                }
                if (has_solvent) {
                    F_solvent /= (1.0 - F[Water]);
                }
                F[Oil] /= (1.0 - F[Water]);
                F[Water] = 0.0;
            }
        }

        if (active()[ Gas ]) {
            if (F[Gas] < 0.0) {
                if (active()[ Water ]) {
                    F[Water] /= (1.0 - F[Gas]);
                }
                if (has_solvent) {
                    F_solvent /= (1.0 - F[Gas]);
                }
                F[Oil] /= (1.0 - F[Gas]);
                F[Gas] = 0.0;
            }
        }

        if (F[Oil] < 0.0) {
            if (active()[ Water ]) {
                F[Water] /= (1.0 - F[Oil]);
            }
            if (active()[ Gas ]) {
                F[Gas] /= (1.0 - F[Oil]);
            }
            if (has_solvent) {
                F_solvent /= (1.0 - F[Oil]);
            }
            F[Oil] = 0.0;
        }

        if (active()[ Water ]) {
            well_state.wellSolutions()[WFrac*nw + indexOfWell()] = F[Water];
        }
        if (active()[ Gas ]) {
            well_state.wellSolutions()[GFrac*nw + indexOfWell()] = F[Gas];
        }
        if(has_solvent) {
            well_state.wellSolutions()[SFrac*nw + indexOfWell()] = F_solvent;
        }

        // F_solvent is added to F_gas. This means that well_rate[Gas] also contains solvent.
        // More testing is needed to make sure this is correct for well groups and THP.
        if (has_solvent){
            F[Gas] += F_solvent;
        }

            // The interpretation of the first well variable depends on the well control
        const WellControls* wc = wellControls();

        // TODO: we should only maintain one current control either from the well_state or from well_controls struct.
        // Either one can be more favored depending on the final strategy for the initilzation of the well control
        const int current = well_state.currentControls()[indexOfWell()];
        const double target_rate = well_controls_iget_target(wc, current);

        std::vector<double> g = {1,1,0.01};
        if (well_controls_iget_type(wc, current) == RESERVOIR_RATE) {
            const double* distr = well_controls_iget_distr(wc, current);
            for (int p = 0; p < np; ++p) {
                if (distr[p] > 0.) { // For injection wells, there only one non-zero distr value
                    F[p] /= distr[p];
                } else {
                    F[p] = 0.;
                }
            }
        } else {
            for (int p = 0; p < np; ++p) {
                F[p] /= g[p];
            }
        }

        switch (well_controls_iget_type(wc, current)) {
            case THP: // The BHP and THP both uses the total rate as first well variable.
            case BHP:
            {
                well_state.wellSolutions()[nw*XvarWell + indexOfWell()] = xvar_well_old[XvarWell] - dwells[0][XvarWell];

                switch (wellType()) {
                case INJECTOR:
                    for (int p = 0; p < np; ++p) {
                        const double comp_frac = compFrac()[p];
                        well_state.wellRates()[indexOfWell() * np + p] = comp_frac * well_state.wellSolutions()[nw*XvarWell + indexOfWell()];
                    }
                    break;
                case PRODUCER:
                    for (int p = 0; p < np; ++p) {
                        well_state.wellRates()[indexOfWell() * np + p] = well_state.wellSolutions()[nw*XvarWell + indexOfWell()] * F[p];
                    }
                    break;
                }

                if (well_controls_iget_type(wc, current) == THP) {

                    // Calculate bhp from thp control and well rates
                    double aqua = 0.0;
                    double liquid = 0.0;
                    double vapour = 0.0;

                    const Opm::PhaseUsage& pu = phaseUsage();

                    if (active()[ Water ]) {
                        aqua = well_state.wellRates()[indexOfWell() * np + pu.phase_pos[ Water ] ];
                    }
                    if (active()[ Oil ]) {
                        liquid = well_state.wellRates()[indexOfWell() * np + pu.phase_pos[ Oil ] ];
                    }
                    if (active()[ Gas ]) {
                        vapour = well_state.wellRates()[indexOfWell() * np + pu.phase_pos[ Gas ] ];
                    }

                    const int vfp        = well_controls_iget_vfp(wc, current);
                    const double& thp    = well_controls_iget_target(wc, current);
                    const double& alq    = well_controls_iget_alq(wc, current);

                    // Set *BHP* target by calculating bhp from THP
                    const WellType& well_type = wellType();
                    // pick the density in the top layer
                    const double rho = perf_densities_[0];
                    const double well_ref_depth = perfDepth()[0];

                    if (well_type == INJECTOR) {
                        const double vfp_ref_depth = vfp_properties_->getInj()->getTable(vfp)->getDatumDepth();

                        const double dp = wellhelpers::computeHydrostaticCorrection(well_ref_depth, vfp_ref_depth, rho, gravity_);

                        well_state.bhp()[indexOfWell()] = vfp_properties_->getInj()->bhp(vfp, aqua, liquid, vapour, thp) - dp;
                    }
                    else if (well_type == PRODUCER) {
                        const double vfp_ref_depth = vfp_properties_->getProd()->getTable(vfp)->getDatumDepth();

                        const double dp = wellhelpers::computeHydrostaticCorrection(well_ref_depth, vfp_ref_depth, rho, gravity_);

                        well_state.bhp()[indexOfWell()] = vfp_properties_->getProd()->bhp(vfp, aqua, liquid, vapour, thp, alq) - dp;
                    }
                    else {
                        OPM_THROW(std::logic_error, "Expected INJECTOR or PRODUCER well");
                    }
                }
            }
                break;
            case SURFACE_RATE: // Both rate controls use bhp as first well variable
            case RESERVOIR_RATE:
            {
                const int sign1 = dwells[0][XvarWell] > 0 ? 1: -1;
                const double dx1_limited = sign1 * std::min(std::abs(dwells[0][XvarWell]),std::abs(xvar_well_old[nw*XvarWell + indexOfWell()])*dBHPLimit);
                well_state.wellSolutions()[nw*XvarWell + indexOfWell()] = std::max(xvar_well_old[nw*XvarWell + indexOfWell()] - dx1_limited,1e5);
                well_state.bhp()[indexOfWell()] = well_state.wellSolutions()[nw*XvarWell + indexOfWell()];

                if (well_controls_iget_type(wc, current) == SURFACE_RATE) {
                    if (wellType() == PRODUCER) {

                        const double* distr = well_controls_iget_distr(wc, current);

                        double F_target = 0.0;
                        for (int p = 0; p < np; ++p) {
                            F_target += distr[p] * F[p];
                        }
                        for (int p = 0; p < np; ++p) {
                            well_state.wellRates()[np * indexOfWell() + p] = F[p] * target_rate / F_target;
                        }
                    } else {

                        for (int p = 0; p < np; ++p) {
                            well_state.wellRates()[indexOfWell() * np + p] = compFrac()[p] * target_rate;
                        }
                    }
                } else { // RESERVOIR_RATE
                    for (int p = 0; p < np; ++p) {
                        well_state.wellRates()[np * indexOfWell() + p] = F[p] * target_rate;
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
                const int current = well_state.currentControls()[indexOfWell()];
                // If under THP control at the moment
                if (current == ctrl_index) {
                    const double thp_target = well_controls_iget_target(wc, current);
                    well_state.thp()[indexOfWell()] = thp_target;
                } else { // otherwise we calculate the thp from the bhp value
                    double aqua = 0.0;
                    double liquid = 0.0;
                    double vapour = 0.0;

                    const Opm::PhaseUsage& pu = phaseUsage();

                    if (active()[ Water ]) {
                        aqua = well_state.wellRates()[indexOfWell()*np + pu.phase_pos[ Water ] ];
                    }
                    if (active()[ Oil ]) {
                        liquid = well_state.wellRates()[indexOfWell()*np + pu.phase_pos[ Oil ] ];
                    }
                    if (active()[ Gas ]) {
                        vapour = well_state.wellRates()[indexOfWell()*np + pu.phase_pos[ Gas ] ];
                    }

                    const double alq = well_controls_iget_alq(wc, ctrl_index);
                    const int table_id = well_controls_iget_vfp(wc, ctrl_index);

                    const WellType& well_type = wellType();
                    const double rho = perf_densities_[0];
                    const double well_ref_depth = perfDepth()[0];
                    if (well_type == INJECTOR) {
                        const double vfp_ref_depth = vfp_properties_->getInj()->getTable(table_id)->getDatumDepth();

                        const double dp = wellhelpers::computeHydrostaticCorrection(well_ref_depth, vfp_ref_depth, rho, gravity_);

                        const double bhp = well_state.bhp()[indexOfWell()];

                        well_state.thp()[indexOfWell()] = vfp_properties_->getInj()->thp(table_id, aqua, liquid, vapour, bhp + dp);
                    } else if (well_type == PRODUCER) {
                        const double vfp_ref_depth = vfp_properties_->getProd()->getTable(table_id)->getDatumDepth();

                        const double dp = wellhelpers::computeHydrostaticCorrection(well_ref_depth, vfp_ref_depth, rho, gravity_);

                        const double bhp = well_state.bhp()[indexOfWell()];

                        well_state.thp()[indexOfWell()] = vfp_properties_->getProd()->thp(table_id, aqua, liquid, vapour, bhp + dp, alq);
                    } else {
                        OPM_THROW(std::logic_error, "Expected INJECTOR or PRODUCER well");
                    }
                }

                // the THP control is found, we leave the loop now
                break;
            }
        }  // end of for loop for seaching THP constraints

        // no THP constraint found
        if (ctrl_index == nwc) { // not finding a THP contstraints
            well_state.thp()[indexOfWell()] = 0.0;
        }
    }





    template<typename TypeTag>
    void
    StandardWell<TypeTag>::
    localInvert(Mat& istlA) const
    {
        for (auto row = istlA.begin(), rowend = istlA.end(); row != rowend; ++row ) {
            for (auto col = row->begin(), colend = row->end(); col != colend; ++col ) {
                //std::cout << (*col) << std::endl;
                (*col).invert();
            }
        }
    }

}
