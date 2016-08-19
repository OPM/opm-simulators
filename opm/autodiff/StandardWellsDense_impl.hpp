/*
  Copyright 2016 SINTEF ICT, Applied Mathematics.
  Copyright 2016 Statoil ASA.
  Copyright 2016 IRIS AS.


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


#include <opm/autodiff/StandardWellsDense.hpp>
#include <opm/autodiff/WellDensitySegmented.hpp>

#include <opm/autodiff/VFPInjProperties.hpp>
#include <opm/autodiff/VFPProdProperties.hpp>
#include <opm/autodiff/WellHelpers.hpp>
#include <opm/autodiff/BlackoilModelEnums.hpp>





namespace Opm
{


    StandardWellsDense::
    WellOps::WellOps(const Wells* wells)
      : w2p(),
        p2w(),
        well_cells()
    {
        if( wells )
        {
            w2p = Eigen::SparseMatrix<double>(wells->well_connpos[ wells->number_of_wells ], wells->number_of_wells);
            p2w = Eigen::SparseMatrix<double>(wells->number_of_wells, wells->well_connpos[ wells->number_of_wells ]);

            const int        nw   = wells->number_of_wells;
            const int* const wpos = wells->well_connpos;

            typedef Eigen::Triplet<double> Tri;

            std::vector<Tri> scatter, gather;
            scatter.reserve(wpos[nw]);
            gather .reserve(wpos[nw]);

            for (int w = 0, i = 0; w < nw; ++w) {
                for (; i < wpos[ w + 1 ]; ++i) {
                    scatter.push_back(Tri(i, w, 1.0));
                    gather .push_back(Tri(w, i, 1.0));
                }
            }

            w2p.setFromTriplets(scatter.begin(), scatter.end());
            p2w.setFromTriplets(gather .begin(), gather .end());

            well_cells.assign(wells->well_cells, wells->well_cells + wells->well_connpos[wells->number_of_wells]);
        }
    }




    StandardWellsDense::StandardWellsDense(const Wells* wells_arg)
      : wells_active_(wells_arg!=nullptr)
      , wells_(wells_arg)
      , wops_(wells_arg)
      , fluid_(nullptr)
      , active_(nullptr)
      , phase_condition_(nullptr)
      , vfp_properties_(nullptr)
      , well_perforation_densities_(Vector())
      , well_perforation_pressure_diffs_(Vector())
      , store_well_perforation_fluxes_(false)
    {
    }





    void
    StandardWellsDense::init(const BlackoilPropsAdInterface* fluid_arg,
                        const std::vector<bool>* active_arg,
                        const std::vector<PhasePresence>* pc_arg,
                        const VFPProperties*  vfp_properties_arg,
                        const double gravity_arg,
                        const Vector& depth_arg)
    {
        fluid_ = fluid_arg;
        active_ = active_arg;
        phase_condition_ = pc_arg;
        vfp_properties_ = vfp_properties_arg;
        gravity_ = gravity_arg;
        perf_cell_depth_ = subset(depth_arg, wellOps().well_cells);;
    }





    const Wells& StandardWellsDense::wells() const
    {
        assert(wells_ != 0);
        return *(wells_);
    }


    const Wells* StandardWellsDense::wellsPointer() const
    {
        return wells_;
    }



    bool StandardWellsDense::wellsActive() const
    {
        return wells_active_;
    }





    void StandardWellsDense::setWellsActive(const bool wells_active)
    {
        wells_active_ = wells_active;
    }





    bool StandardWellsDense::localWellsActive() const
    {
        return wells_ ? (wells_->number_of_wells > 0 ) : false;
    }





    int
    StandardWellsDense::numWellVars() const
    {
        if ( !localWellsActive() )
        {
            return 0;
        }

        // For each well, we have a bhp variable, and one flux per phase.
        const int nw = wells().number_of_wells;
        return (numPhases() + 1) * nw;
    }





    const StandardWellsDense::WellOps&
    StandardWellsDense::wellOps() const
    {
        return wops_;
    }





    StandardWellsDense::Vector& StandardWellsDense::wellPerforationDensities()
    {
        return well_perforation_densities_;
    }





    const StandardWellsDense::Vector&
    StandardWellsDense::wellPerforationDensities() const
    {
        return well_perforation_densities_;
    }





    StandardWellsDense::Vector&
    StandardWellsDense::wellPerforationPressureDiffs()
    {
        return well_perforation_pressure_diffs_;
    }





    const StandardWellsDense::Vector&
    StandardWellsDense::wellPerforationPressureDiffs() const
    {
        return well_perforation_pressure_diffs_;
    }




    template<class SolutionState, class WellState>
    void
    StandardWellsDense::
    computePropertiesForWellConnectionPressures(const SolutionState& state,
                                                const WellState& xw,
                                                std::vector<double>& b_perf,
                                                std::vector<double>& rsmax_perf,
                                                std::vector<double>& rvmax_perf,
                                                std::vector<double>& surf_dens_perf)
    {
        const int nperf = wells().well_connpos[wells().number_of_wells];
        const int nw = wells().number_of_wells;

        // Compute the average pressure in each well block
        const Vector perf_press = Eigen::Map<const Vector>(xw.perfPress().data(), nperf);
        Vector avg_press = perf_press*0;
        for (int w = 0; w < nw; ++w) {
            for (int perf = wells().well_connpos[w]; perf < wells().well_connpos[w+1]; ++perf) {
                const double p_above = perf == wells().well_connpos[w] ? state.bhp.value()[w] : perf_press[perf - 1];
                const double p_avg = (perf_press[perf] + p_above)/2;
                avg_press[perf] = p_avg;
            }
        }

        const std::vector<int>& well_cells = wellOps().well_cells;

        // Use cell values for the temperature as the wells don't knows its temperature yet.
        const ADB perf_temp = subset(state.temperature, well_cells);

        // Compute b, rsmax, rvmax values for perforations.
        // Evaluate the properties using average well block pressures
        // and cell values for rs, rv, phase condition and temperature.
        const ADB avg_press_ad = ADB::constant(avg_press);
        std::vector<PhasePresence> perf_cond(nperf);
        // const std::vector<PhasePresence>& pc = phaseCondition();
        for (int perf = 0; perf < nperf; ++perf) {
            perf_cond[perf] = (*phase_condition_)[well_cells[perf]];
        }
        const PhaseUsage& pu = fluid_->phaseUsage();
        DataBlock b(nperf, pu.num_phases);
        if (pu.phase_used[BlackoilPhases::Aqua]) {
            const Vector bw = fluid_->bWat(avg_press_ad, perf_temp, well_cells).value();
            b.col(pu.phase_pos[BlackoilPhases::Aqua]) = bw;
        }
        assert((*active_)[Oil]);
        const Vector perf_so =  subset(state.saturation[pu.phase_pos[Oil]].value(), well_cells);
        if (pu.phase_used[BlackoilPhases::Liquid]) {
            const ADB perf_rs = (state.rs.size() > 0) ? subset(state.rs, well_cells) : ADB::null();
            const Vector bo = fluid_->bOil(avg_press_ad, perf_temp, perf_rs, perf_cond, well_cells).value();
            b.col(pu.phase_pos[BlackoilPhases::Liquid]) = bo;
        }
        if (pu.phase_used[BlackoilPhases::Vapour]) {
            const ADB perf_rv = (state.rv.size() > 0) ? subset(state.rv, well_cells) : ADB::null();
            const Vector bg = fluid_->bGas(avg_press_ad, perf_temp, perf_rv, perf_cond, well_cells).value();
            b.col(pu.phase_pos[BlackoilPhases::Vapour]) = bg;
        }
        if (pu.phase_used[BlackoilPhases::Liquid] && pu.phase_used[BlackoilPhases::Vapour]) {
            const Vector rssat = fluid_->rsSat(ADB::constant(avg_press), ADB::constant(perf_so), well_cells).value();
            rsmax_perf.assign(rssat.data(), rssat.data() + nperf);

            const Vector rvsat = fluid_->rvSat(ADB::constant(avg_press), ADB::constant(perf_so), well_cells).value();
            rvmax_perf.assign(rvsat.data(), rvsat.data() + nperf);
        }

        // b is row major, so can just copy data.
        b_perf.assign(b.data(), b.data() + nperf * pu.num_phases);

        // Surface density.
        // The compute density segment wants the surface densities as
        // an np * number of wells cells array
        Vector rho = superset(fluid_->surfaceDensity(0 , well_cells), Span(nperf, pu.num_phases, 0), nperf*pu.num_phases);
        for (int phase = 1; phase < pu.num_phases; ++phase) {
            rho += superset(fluid_->surfaceDensity(phase , well_cells), Span(nperf, pu.num_phases, phase), nperf*pu.num_phases);
        }
        surf_dens_perf.assign(rho.data(), rho.data() + nperf * pu.num_phases);

    }





    template <class WellState>
    void
    StandardWellsDense::
    computeWellConnectionDensitesPressures(const WellState& xw,
                                           const std::vector<double>& b_perf,
                                           const std::vector<double>& rsmax_perf,
                                           const std::vector<double>& rvmax_perf,
                                           const std::vector<double>& surf_dens_perf,
                                           const std::vector<double>& depth_perf,
                                           const double grav)
    {
        // Compute densities
        std::vector<double> cd =
                WellDensitySegmented::computeConnectionDensities(
                        wells(), xw, fluid_->phaseUsage(),
                        b_perf, rsmax_perf, rvmax_perf, surf_dens_perf);

        const int nperf = wells().well_connpos[wells().number_of_wells];

        // Compute pressure deltas
        std::vector<double> cdp =
                WellDensitySegmented::computeConnectionPressureDelta(
                        wells(), depth_perf, cd, grav);

        // Store the results
        well_perforation_densities_ = Eigen::Map<const Vector>(cd.data(), nperf);
        well_perforation_pressure_diffs_ = Eigen::Map<const Vector>(cdp.data(), nperf);
    }





    template <class SolutionState, class WellState>
    void
    StandardWellsDense::
    computeWellConnectionPressures(const SolutionState& state,
                                   const WellState& xw)
    {
        if( ! localWellsActive() ) return ;
        // 1. Compute properties required by computeConnectionPressureDelta().
        //    Note that some of the complexity of this part is due to the function
        //    taking std::vector<double> arguments, and not Eigen objects.
        std::vector<double> b_perf;
        std::vector<double> rsmax_perf;
        std::vector<double> rvmax_perf;
        std::vector<double> surf_dens_perf;
        computePropertiesForWellConnectionPressures(state, xw, b_perf, rsmax_perf, rvmax_perf, surf_dens_perf);

        const Vector& pdepth = perf_cell_depth_;
        const int nperf = wells().well_connpos[wells().number_of_wells];
        const std::vector<double> depth_perf(pdepth.data(), pdepth.data() + nperf);

        computeWellConnectionDensitesPressures(xw, b_perf, rsmax_perf, rvmax_perf, surf_dens_perf, depth_perf, gravity_);

    }





    template <class ReservoirResidualQuant, class SolutionState>
    void
    StandardWellsDense::
    extractWellPerfProperties(const SolutionState& /* state */,
                              const std::vector<ReservoirResidualQuant>& rq,
                              std::vector<ADB>& mob_perfcells,
                              std::vector<ADB>& b_perfcells) const
    {
        // If we have wells, extract the mobilities and b-factors for
        // the well-perforated cells.
        if ( !localWellsActive() ) {
            mob_perfcells.clear();
            b_perfcells.clear();
            return;
        } else {
            const std::vector<int>& well_cells = wellOps().well_cells;
            const int num_phases = wells().number_of_phases;
            mob_perfcells.resize(num_phases, ADB::null());
            b_perfcells.resize(num_phases, ADB::null());
            for (int phase = 0; phase < num_phases; ++phase) {
                mob_perfcells[phase] = subset(rq[phase].mob, well_cells);
                b_perfcells[phase] = subset(rq[phase].b, well_cells);
            }
        }
    }

    template <class SolutionState>
    void
    StandardWellsDense::
    computeWellFluxDense(const SolutionState& state,
                    const std::vector<ADB>& mob_perfcells,
                    const std::vector<ADB>& b_perfcells,
                    std::vector<ADB>& cq_s) const
    {
        if( ! localWellsActive() ) return ;

        const int np = wells().number_of_phases;
        const int nw = wells().number_of_wells;
        const int nperf = wells().well_connpos[nw];
        const Opm::PhaseUsage& pu = fluid_->phaseUsage();
        Vector Tw = Eigen::Map<const Vector>(wells().WI, nperf);
        const std::vector<int>& well_cells = wellOps().well_cells;
        std::vector<int> well_id(nperf);

        // pressure diffs computed already (once per step, not changing per iteration)
        const Vector& cdp = wellPerforationPressureDiffs();

        std::vector<std::vector<EvalWell>> cq_s_dense(np, std::vector<EvalWell>(nperf,0.0));
        std::vector<ADB> cmix_s_ADB = wellVolumeFractions(state);

        const int oilpos = pu.phase_pos[Oil];
        ADB perfpressure = (wellOps().w2p * state.bhp) + cdp;
        const ADB rsSat = fluid_->rsSat(perfpressure, cmix_s_ADB[oilpos], well_cells);
        const ADB rvSat = fluid_->rvSat(perfpressure, cmix_s_ADB[oilpos], well_cells);

        for (int w = 0; w < nw; ++w) {

            EvalWell bhp = extractDenseADWell(state.bhp,w);


            // TODO: fix for 2-phase case
            std::vector<EvalWell> cmix_s(np,0.0);
            for (int phase = 0; phase < np; ++phase) {
                cmix_s[phase] = extractDenseADWell(cmix_s_ADB[phase],w);
            }

            //std::cout <<"cmix gas "<< w<< " "<<cmix_s[Gas] << std::endl;

            for (int perf = wells().well_connpos[w] ; perf < wells().well_connpos[w+1]; ++perf) {
                const int cell_idx = well_cells[perf];
                well_id[perf] = w;
                EvalWell pressure = extractDenseAD(state.pressure, cell_idx, cell_idx);
                EvalWell rs = extractDenseAD(state.rs, cell_idx, cell_idx);
                EvalWell rv = extractDenseAD(state.rv, cell_idx, cell_idx);
                std::vector<EvalWell> b_perfcells_dense(np, 0.0);
                std::vector<EvalWell> mob_perfcells_dense(np, 0.0);
                for (int phase = 0; phase < np; ++phase) {
                    b_perfcells_dense[phase] = extractDenseAD(b_perfcells[phase], perf, cell_idx);
                    mob_perfcells_dense[phase] = extractDenseAD(mob_perfcells[phase], perf, cell_idx);

                }

                // Pressure drawdown (also used to determine direction of flow)
                EvalWell drawdown = pressure - bhp - cdp[perf];

                // injection perforations
                if ( drawdown.value > 0 )  {

                    //Do nothing if crossflow is not allowed
                    if (!wells().allow_cf[w] && wells().type[w] == INJECTOR)
                        continue;
                    // compute phase volumetric rates at standard conditions
                    std::vector<EvalWell> cq_ps(np, 0.0);
                    for (int phase = 0; phase < np; ++phase) {
                        const EvalWell cq_p = - Tw[perf] * (mob_perfcells_dense[phase] * drawdown);
                        cq_ps[phase] = b_perfcells_dense[phase] * cq_p;
                    }


                    if ((*active_)[Oil] && (*active_)[Gas]) {
                        const int oilpos = pu.phase_pos[Oil];
                        const int gaspos = pu.phase_pos[Gas];
                        const EvalWell cq_psOil = cq_ps[oilpos];
                        const EvalWell cq_psGas = cq_ps[gaspos];
                        cq_ps[gaspos] += rs * cq_psOil;
                        cq_ps[oilpos] += rv * cq_psGas;
                    }


                    // map to ADB
                    for (int phase = 0; phase < np; ++phase) {
                        cq_s_dense[phase][perf] = cq_ps[phase];
                    }

                } else {
                    //Do nothing if crossflow is not allowed
                    if (!wells().allow_cf[w] && wells().type[w] == PRODUCER)
                        continue;

                    // Using total mobilities
                    EvalWell total_mob_dense = mob_perfcells_dense[0];
                    for (int phase = 1; phase < np; ++phase) {
                        total_mob_dense += mob_perfcells_dense[phase];
                    }
                    // injection perforations total volume rates
                    const EvalWell cqt_i = - Tw[perf] * (total_mob_dense * drawdown);

                    // compute volume ratio between connection at standard conditions
                    EvalWell volumeRatio = 0.0;
                    if ((*active_)[Water]) {
                        const int watpos = pu.phase_pos[Water];
                        volumeRatio += cmix_s[watpos] / b_perfcells_dense[watpos];
                    }

                    if ((*active_)[Oil] && (*active_)[Gas]) {

                        const int oilpos = pu.phase_pos[Oil];
                        const int gaspos = pu.phase_pos[Gas];
                        EvalWell rvPerf = 0.0;
                        if (cmix_s[gaspos] > 0)
                            rvPerf = cmix_s[oilpos] / cmix_s[gaspos];

                        if (rvPerf.value > rvSat.value()[w]) {
                            rvPerf = 0.0;
                            rvPerf.value = rvSat.value()[w];
                        }

                        EvalWell rsPerf = 0.0;
                        if (cmix_s[oilpos] > 0)
                            rsPerf = cmix_s[gaspos] / cmix_s[oilpos];

                        if (rsPerf.value > rsSat.value()[w]) {
                            rsPerf = 0.0;
                            rsPerf.value = rsSat.value()[w];
                        }

                        // Incorporate RS/RV factors if both oil and gas active
                        const EvalWell d = 1.0 - rvPerf * rsPerf;

                        const EvalWell tmp_oil = (cmix_s[oilpos] - rvPerf * cmix_s[gaspos]) / d;
                        //std::cout << "tmp_oil " <<tmp_oil << std::endl;
                        volumeRatio += tmp_oil / b_perfcells_dense[oilpos];

                        const EvalWell tmp_gas = (cmix_s[gaspos] - rsPerf * cmix_s[oilpos]) / d;
                        //std::cout << "tmp_gas " <<tmp_gas << std::endl;
                        volumeRatio += tmp_gas / b_perfcells_dense[gaspos];
                    }
                    else {
                        if ((*active_)[Oil]) {
                            const int oilpos = pu.phase_pos[Oil];
                            volumeRatio += cmix_s[oilpos] / b_perfcells_dense[oilpos];
                        }
                        if ((*active_)[Gas]) {
                            const int gaspos = pu.phase_pos[Gas];
                            volumeRatio += cmix_s[gaspos] / b_perfcells_dense[gaspos];
                        }
                    }
                    // injecting connections total volumerates at standard conditions
                    EvalWell cqt_is = cqt_i/volumeRatio;
                    //std::cout << "volrat " << volumeRatio << " " << volrat_perf_[perf] << std::endl;
                    for (int phase = 0; phase < np; ++phase) {
                        cq_s_dense[phase][perf] = cmix_s[phase] * cqt_is; // * b_perfcells_dense[phase];
                    }
                }
            }
        }
        cq_s.resize(np, ADB::null());
        for (int phase = 0; phase < np; ++phase) {
            cq_s[phase] = convertToADB(cq_s_dense[phase], well_cells, state.pressure.size(), well_id, nw, state.bhp.numBlocks());
            //std::cout << "cq_s " <<cq_s[phase] << std::endl;
        }
        //std::cout << aliveWells << std::endl;
        //std::vector<ADB> cq_s2;
        //Vector aliveWells;
        //computeWellFlux(state,mob_perfcells,b_perfcells, aliveWells,cq_s2);

        //for (int phase = 0; phase < np; ++phase) {
            //if( !(((cq_s[phase].value() - cq_s2[phase].value()).abs()<1e-10).all()) ) {
                //std::cout << "phase " << phase << std::endl;
                //std::cout << cq_s2[phase].value() << std::endl;
                //std::cout << cq_s[phase].value() << std::endl;
            //}
        //}
    }






    template <class SolutionState>
    void
    StandardWellsDense::
    computeWellFlux(const SolutionState& state,
                    const std::vector<ADB>& mob_perfcells,
                    const std::vector<ADB>& b_perfcells,
                    Vector& aliveWells,
                    std::vector<ADB>& cq_s) const
    {
        if( ! localWellsActive() ) return ;

        const int np = wells().number_of_phases;
        const int nw = wells().number_of_wells;
        const int nperf = wells().well_connpos[nw];
        Vector Tw = Eigen::Map<const Vector>(wells().WI, nperf);
        const std::vector<int>& well_cells = wellOps().well_cells;

        // pressure diffs computed already (once per step, not changing per iteration)
        const Vector& cdp = wellPerforationPressureDiffs();
        // Extract needed quantities for the perforation cells
        const ADB& p_perfcells = subset(state.pressure, well_cells);

        // Perforation pressure
        const ADB perfpressure = (wellOps().w2p * state.bhp) + cdp;

        // Pressure drawdown (also used to determine direction of flow)
        const ADB drawdown =  p_perfcells - perfpressure;

        // Compute vectors with zero and ones that
        // selects the wanted quantities.

        // selects injection perforations
        Vector selectInjectingPerforations = Vector::Zero(nperf);
        // selects producing perforations
        Vector selectProducingPerforations = Vector::Zero(nperf);
        for (int c = 0; c < nperf; ++c){
            if (drawdown.value()[c] < 0)
                selectInjectingPerforations[c] = 1;
            else
                selectProducingPerforations[c] = 1;
        }
        std::vector<V> distri(np);
        for (int p = 0; p < np; ++p) {
            distri[p] = V::Zero(nw);
        }


        Vector isResv = Vector::Zero(nw);
        for (int w = 0; w < nw; ++w) {
            const WellControls* wc = wells().ctrls[w];
            // The current control in the well state overrides
            // the current control set in the Wells struct, which
            // is instead treated as a default.
            const double* distr = well_controls_get_current_distr(wc);
            if (well_controls_get_current_type(wc) == RESERVOIR_RATE) {
                isResv[w] = 1;
            }
            for (int p = 0; p < np; ++p) {
                distri[p][w] = distr[p];
                //std::cout << "distr " << distr[p] << " comp_frac " << comp_frac << std::endl;
            }

        }


        // Handle cross flow
        const Vector numInjectingPerforations = (wellOps().p2w * ADB::constant(selectInjectingPerforations)).value();
        const Vector numProducingPerforations = (wellOps().p2w * ADB::constant(selectProducingPerforations)).value();
        for (int w = 0; w < nw; ++w) {
            if (!wells().allow_cf[w]) {
                for (int perf = wells().well_connpos[w] ; perf < wells().well_connpos[w+1]; ++perf) {
                    // Crossflow is not allowed; reverse flow is prevented.
                    // At least one of the perforation must be open in order to have a meeningful
                    // equation to solve. For the special case where all perforations have reverse flow,
                    // and the target rate is non-zero all of the perforations are keept open.
                    if (wells().type[w] == INJECTOR && numInjectingPerforations[w] > 0) {
                        selectProducingPerforations[perf] = 0.0;
                    } else if (wells().type[w] == PRODUCER && numProducingPerforations[w] > 0 ){
                        selectInjectingPerforations[perf] = 0.0;
                    }
                }
            }
        }

        // HANDLE FLOW INTO WELLBORE
        // compute phase volumetric rates at standard conditions
        std::vector<ADB> cq_p(np, ADB::null());
        std::vector<ADB> cq_ps(np, ADB::null());
        for (int phase = 0; phase < np; ++phase) {
            cq_p[phase] = -(selectProducingPerforations * Tw) * (mob_perfcells[phase] * drawdown);
            cq_ps[phase] = b_perfcells[phase] * cq_p[phase];
        }
        const Opm::PhaseUsage& pu = fluid_->phaseUsage();
        if ((*active_)[Oil] && (*active_)[Gas]) {
            const int oilpos = pu.phase_pos[Oil];
            const int gaspos = pu.phase_pos[Gas];
            const ADB cq_psOil = cq_ps[oilpos];
            const ADB cq_psGas = cq_ps[gaspos];
            const ADB& rv_perfcells = subset(state.rv, well_cells);
            const ADB& rs_perfcells = subset(state.rs, well_cells);
            cq_ps[gaspos] += rs_perfcells * cq_psOil;
            cq_ps[oilpos] += rv_perfcells * cq_psGas;
        }

        // HANDLE FLOW OUT FROM WELLBORE
        // Using total mobilities
        ADB total_mob = mob_perfcells[0];
        for (int phase = 1; phase < np; ++phase) {
            total_mob += mob_perfcells[phase];
        }
        // injection perforations total volume rates
        const ADB cqt_i = -(selectInjectingPerforations * Tw) * (total_mob * drawdown);

        // Store well perforation total fluxes (reservor volumes) if requested.
        if (store_well_perforation_fluxes_) {
            // Ugly const-cast, but unappealing alternatives.
            Vector& wf = const_cast<Vector&>(well_perforation_fluxes_);
            wf = cqt_i.value();
            for (int phase = 0; phase < np; ++phase) {
                wf += cq_p[phase].value();
            }
        }

        // compute wellbore mixture for injecting perforations
        // The wellbore mixture depends on the inflow from the reservoar
        // and the well injection rates.

        // compute avg. and total wellbore phase volumetric rates at standard conds
//        const DataBlock compi = Eigen::Map<const DataBlock>(wells().comp_frac, nw, np);
//        std::vector<ADB> wbq(np, ADB::null());
//        ADB wbqt = ADB::constant(Vector::Zero(nw));
//        for (int phase = 0; phase < np; ++phase) {
//            const ADB& q_ps = wellOps().p2w * cq_ps[phase];
//            const ADB& q_s = subset(state.qs, Span(nw, 1, phase*nw));
//            Selector<double> injectingPhase_selector(q_s.value(), Selector<double>::GreaterZero);
//            const int pos = pu.phase_pos[phase];
//            wbq[phase] = (compi.col(pos) * injectingPhase_selector.select(q_s,ADB::constant(Vector::Zero(nw))))  - q_ps;
//            wbqt += wbq[phase];
//        }
        // compute wellbore mixture at standard conditions.


        std::vector<ADB> cmix_s(np, ADB::null());
        Vector ones = Vector::Constant(nperf,1.0);
        cmix_s[Water] = wellOps().w2p * subset(state.wellVariables, Span(nw, 1, 1*nw));
        cmix_s[Gas] = wellOps().w2p * subset(state.wellVariables, Span(nw, 1, 2*nw));
        cmix_s[Oil] = (ones  - cmix_s[Water] - cmix_s[Gas]);
//        std::cout << cmix_s[Gas].value() << std::endl;

//        std::vector<Vector> g(np, Vector::Constant(nw,1.0));
//        g[Gas] = Vector::Constant(nw,0.01);
//        const DataBlock compi = Eigen::Map<const DataBlock>(wells().comp_frac, nw, np);
//        std::vector<ADB> wbq(np, ADB::null());
//        ADB wbqt = ADB::constant(Vector::Zero(nw));
//        for (int phase = 0; phase < np; ++phase) {
//            const ADB& q_s = subset(state.qs, Span(nw, 1, phase*nw));
//            wbq[phase] = q_s * ( (Vector::Constant(nw,1.0) - isResv) * g[phase] + isResv * distri[phase]);
//            wbqt += wbq[phase];
//        }
//        Selector<double> notDeadWells_selector(wbqt.value(), Selector<double>::Zero);
//        std::vector<ADB> cmix_s(np, ADB::null());
//        for (int phase = 0; phase < np; ++phase) {
//            const int pos = pu.phase_pos[phase];
//            cmix_s[phase] = wellOps().w2p * notDeadWells_selector.select(ADB::constant(compi.col(pos)), wbq[phase]/wbqt);
//        }


        // compute volume ratio between connection at standard conditions
        ADB volumeRatio = ADB::constant(Vector::Zero(nperf));

        if ((*active_)[Water]) {
            const int watpos = pu.phase_pos[Water];
            volumeRatio += cmix_s[watpos] / b_perfcells[watpos];
        }


        if ((*active_)[Oil] && (*active_)[Gas]) {
            // Incorporate RS/RV factors if both oil and gas active
//            const ADB& rv = subset(state.rv, well_cells);
//            const ADB& rs = subset(state.rs, well_cells);
//            const ADB d = Vector::Constant(nperf,1.0) - rv* rs;




            const int oilpos = pu.phase_pos[Oil];
            const int gaspos = pu.phase_pos[Gas];

            Selector<double> noGas_selector(cmix_s[gaspos].value(), Selector<double>::GreaterZero);
            Selector<double> noOil_selector(cmix_s[oilpos].value(), Selector<double>::GreaterZero);
            ADB rv = noGas_selector.select(cmix_s[oilpos]/cmix_s[gaspos], ADB::constant(Vector::Constant(nperf,0.0)));
            ADB rs = noOil_selector.select(cmix_s[gaspos]/cmix_s[oilpos], ADB::constant(Vector::Constant(nperf,0.0)));

            const ADB rsSat = fluid_->rsSat(perfpressure, cmix_s[oilpos], well_cells);
            const ADB rvSat = fluid_->rvSat(perfpressure, cmix_s[oilpos], well_cells);

            Selector<double> maxRs_selector(rs.value() - rsSat.value(), Selector<double>::GreaterZero);
            Selector<double> maxRv_selector(rv.value() - rvSat.value(), Selector<double>::GreaterZero);
            rv = maxRv_selector.select(rvSat, rv);
            rs = maxRs_selector.select(rsSat, rs);

            const ADB d = Vector::Constant(nperf,1.0) - rv * rs;


            const ADB tmp_oil = (cmix_s[oilpos] - rv * cmix_s[gaspos]) / d;
            volumeRatio += tmp_oil / b_perfcells[oilpos];

            const ADB tmp_gas = (cmix_s[gaspos] - rs * cmix_s[oilpos]) / d;
            volumeRatio += tmp_gas / b_perfcells[gaspos];
        }
        else {
            if ((*active_)[Oil]) {
                const int oilpos = pu.phase_pos[Oil];
                volumeRatio += cmix_s[oilpos] / b_perfcells[oilpos];
            }
            if ((*active_)[Gas]) {
                const int gaspos = pu.phase_pos[Gas];
                volumeRatio += cmix_s[gaspos] / b_perfcells[gaspos];
            }
        }


        // injecting connections total volumerates at standard conditions
        //std::cout << volumeRatio.value() << std::endl;
        ADB cqt_is = cqt_i/volumeRatio;

        // connection phase volumerates at standard conditions
        cq_s.resize(np, ADB::null());
        for (int phase = 0; phase < np; ++phase) {
            cq_s[phase] = cq_ps[phase] + cmix_s[phase]*cqt_is; //*b_perfcells[phase];
        }

        // check for dead wells (used in the well controll equations)
        aliveWells = Vector::Constant(nw, 1.0);
//        for (int w = 0; w < nw; ++w) {
//            if (wbqt.value()[w] == 0) {
//                aliveWells[w] = 0.0;
//            }
//        }
    }

    typedef DenseAd::Evaluation<double, /*size=*/6> EvalWell;
    EvalWell
    StandardWellsDense::
    extractDenseAD(const ADB& data, int i, int j) const
    {
        EvalWell output = 0.0;
        output.value = data.value()[i];
        const int np = wells().number_of_phases;
        const std::vector<Opm::AutoDiffMatrix>& jac = data.derivative();
        //std::cout << jac.size() << std::endl;
        int numblocs = jac.size();
        for (int b = 0; b < numblocs; ++b) {
            if (b < np) { // don't copy well blocks)
                //std::cout << jac[b].coeff(i,j) << std::endl;
                output.derivatives[b] = jac[b].coeff(i,j);
            }
        }
        return output;
    }

    typedef DenseAd::Evaluation<double, /*size=*/6> EvalWell;
    EvalWell
    StandardWellsDense::
    extractDenseADWell(const ADB& data, int i) const
    {
        EvalWell output = 0.0;
        output.value = data.value()[i];
        const int nw = wells().number_of_wells;
        const int np = wells().number_of_phases;
        const std::vector<Opm::AutoDiffMatrix>& jac = data.derivative();
        //std::cout << jac.size() << std::endl;
        int numblocs = jac.size();
        for (int b = 0; b < np; ++b) {
            output.derivatives[b+np] = jac[numblocs-1].coeff(i, b*nw + i);
        }
        return output;
    }

    const AutoDiffBlock<double> StandardWellsDense::convertToADB(const std::vector<EvalWell>& local, const std::vector<int>& well_cells, const int nc, const std::vector<int>& well_id, const int nw, const int numVars) const
    {
        typedef typename ADB::M  M;
        const int nLocal = local.size();
        typename ADB::V value( nLocal );
        //const int numVars = 5;
        const int np = wells().number_of_phases;

        std::vector<Eigen::SparseMatrix<double>> mat(np, Eigen::SparseMatrix<double>(nLocal,nc));
        Eigen::SparseMatrix<double> matFlux(nLocal,np*nw);
        Eigen::SparseMatrix<double> matBHP(nLocal,nw);

        for( int i=0; i<nLocal; ++i )
        {
            value[ i ] = local[ i ].value;
            for( int d=0; d<np; ++d ) {
                //std::cout << i << " " <<d << " "<<local[i].derivatives[d] << std::endl;
                mat[d].insert(i, well_cells[i]) = local[i].derivatives[d];
            }

            for (int phase = 0; phase < np; ++phase) {
                //std::cout << "well: "<< i << " " << phase << " " << local[i].derivatives[np + phase] << std::endl;

                matFlux.insert(i, nw*phase + well_id[i]) = local[i].derivatives[np + phase];
            }
            //matBHP.insert(i,well_id[i]) = local[i].derivatives[2*np];
        }

        std::vector< M > jacs( numVars );
        if (numVars == 4) {
            for( int d=0; d<np; ++d ) {
                //Eigen::DiagonalMatrix<double>(deri[d]);
                jacs[ d ] = M(mat[d]);
            }

            jacs[3] = M(matFlux);
            //jacs[4] = M(matBHP);
        }
        else if (numVars == 1) {
            jacs[0] = M(matFlux);
            //jacs[1] = M(matBHP);
        }
        //std::cout << numVars << std::endl;

        return ADB::function( std::move( value ), std::move( jacs ));
    }




    template <class SolutionState, class WellState>
    void
    StandardWellsDense::
    updatePerfPhaseRatesAndPressures(const std::vector<ADB>& cq_s,
                                     const SolutionState& state,
                                     WellState& xw) const
    {
        if ( !localWellsActive() )
        {
            // If there are no wells in the subdomain of the proces then
            // cq_s has zero size and will cause a segmentation fault below.
            return;
        }

        // Update the perforation phase rates (used to calculate the pressure drop in the wellbore).
        const int np = wells().number_of_phases;
        const int nw = wells().number_of_wells;
        const int nperf = wells().well_connpos[nw];

        Vector cq = superset(cq_s[0].value(), Span(nperf, np, 0), nperf*np);
        for (int phase = 1; phase < np; ++phase) {
            cq += superset(cq_s[phase].value(), Span(nperf, np, phase), nperf*np);
        }
        xw.perfPhaseRates().assign(cq.data(), cq.data() + nperf*np);

        // Update the perforation pressures.
        const Vector& cdp = wellPerforationPressureDiffs();
        const Vector perfpressure = (wellOps().w2p * state.bhp.value().matrix()).array() + cdp;
        xw.perfPress().assign(perfpressure.data(), perfpressure.data() + nperf);
    }







    template <class WellState>
    void
    StandardWellsDense::
    updateWellState(const Vector& dwells,
                    const double dpmaxrel,
                    WellState& well_state)
    {
        if( localWellsActive() )
        {
            const int np = wells().number_of_phases;
            const int nw = wells().number_of_wells;

            // Extract parts of dwells corresponding to each part.
            int varstart = 0;
            const Vector dxvar_well = subset(dwells, Span(np*nw, 1, varstart));
            //const Vector dqs = subset(dwells, Span(np*nw, 1, varstart));
            varstart += dxvar_well.size();
            //const Vector dbhp = subset(dwells, Span(nw, 1, varstart));
            //varstart += dbhp.size();
            assert(varstart == dwells.size());

            // Qs update.
            // Since we need to update the wellrates, that are ordered by wells,
            // from dqs which are ordered by phase, the simplest is to compute
            // dwr, which is the data from dqs but ordered by wells.
            const Vector xvar_well_old = Eigen::Map<const Vector>(&well_state.wellSolutions()[0], nw*np);

            double dFLimit = 0.2;
            double dBHPLimit = 2;
            double dTotalRateLimit = 0.5;
            //std::cout << "dxvar_well "<<dxvar_well << std::endl;


            for (int w = 0; w < nw; ++w) {
                const WellControls* wc = wells().ctrls[w];
                // The current control in the well state overrides
                // the current control set in the Wells struct, which
                // is instead treated as a default.
                const int current = well_state.currentControls()[w];
                const double target_rate = well_controls_iget_target(wc, current);
                const double* distr = well_controls_iget_distr(wc, current);
                std::vector<double> F(np,0.0);
                const int sign2 = dxvar_well[nw + w] > 0 ? 1: -1;
                const double dx2_limited = sign2 * std::min(std::abs(dxvar_well[nw + w]),dFLimit);
                well_state.wellSolutions()[nw + w] = xvar_well_old[nw + w] - dx2_limited;
                const int sign3 = dxvar_well[2*nw + w] > 0 ? 1: -1;
                const double dx3_limited = sign3 * std::min(std::abs(dxvar_well[2*nw + w]),dFLimit);
                well_state.wellSolutions()[2*nw + w] = xvar_well_old[2*nw + w] - dx3_limited;
                F[Water] = well_state.wellSolutions()[nw + w];
                F[Gas] = well_state.wellSolutions()[2*nw + w];
                F[Oil] = 1.0 - F[Water] - F[Gas];

//                const double dFw = dxvar_well[nw + w];
//                const double dFg = dxvar_well[nw*2 + w];
//                const double dFo = - dFw - dFg;
//                //std::cout << w << " "<< F[Water] << " "  << F[Oil] << " " << F[Gas] << std::endl;
//                double step = dFLimit / std::max(std::abs(dFw),std::max(std::abs(dFg),std::abs(dFo))); //)) / dFLimit;
//                step = std::min(step, 1.0);
//                //std::cout << step << std::endl;
//                F[Water] = xvar_well_old[nw + w] - step*dFw;
//                F[Gas] = xvar_well_old[2*nw + w] - step*dFg;
//                F[Oil] = (1.0 - xvar_well_old[2*nw + w] - xvar_well_old[nw + w]) - step * dFo;

                if (F[Water] < 0.0) {
                    F[Gas] /= (1.0 - F[Water]);
                    F[Oil] /= (1.0 - F[Water]);
                    F[Water] = 0.0;
                }
                if (F[Gas] < 0.0) {
                    F[Water] /= (1.0 - F[Gas]);
                    F[Oil] /= (1.0 - F[Gas]);
                    F[Gas] = 0.0;
                }
                if (F[Oil] < 0.0) {
                    F[Water] /= (1.0 - F[Oil]);
                    F[Gas] /= (1.0 - F[Oil]);
                    F[Oil] = 0.0;
                }
                well_state.wellSolutions()[nw + w] = F[Water];
                well_state.wellSolutions()[2*nw + w] = F[Gas];

                //std::cout << wells().name[w] << " "<< F[Water] << " "  << F[Oil] << " " << F[Gas] << std::endl;


                std::vector<double> g = {1,1,0.01};

                if (well_controls_iget_type(wc, current) == RESERVOIR_RATE) {
                    for (int p = 0; p < np; ++p) {
                        F[p] /= distr[p];
                    }
                } else {
                    for (int p = 0; p < np; ++p) {
                        F[p] /= g[p];
                    }
                }
                //std::cout << w << " "<< F[Water] << " "  << F[Oil] << " " << F[Gas] << std::endl;






//                const double dFw = dxvar_well[nw + w];
//                const double dFg = dxvar_well[nw*2 + w];
//                const double dFo = - dFw - dFg;
                    //std::cout << w << " "<< F[Water] << " "  << F[Oil] << " " << F[Gas] << std::endl;
//                double step = dFLimit / std::max(std::abs(dFw),std::max(std::abs(dFg),std::abs(dFo))); //)) / dFLimit;
//                step = std::min(step, 1.0);
//                std::cout << step << std::endl;
//                F[Water] = xvar_well_old[nw + w] - step*dFw;
//                F[Gas] = xvar_well_old[2*nw + w] - step*dFg;
//                F[Oil] = (1.0 - xvar_well_old[2*nw + w] - xvar_well_old[nw + w]) - step * dFo;
//                double sumF = F[Water]+F[Gas]+F[Oil];
//                F[Water] /= sumF;
//                F[Gas] /= sumF;
//                F[Oil] /= sumF;
//                well_state.wellSolutions()[nw + w] = F[Water];
//                well_state.wellSolutions()[2 * nw + w] = F[Gas];

                switch (well_controls_iget_type(wc, current)) {
                case BHP:
                {
                    //const int sign1 = dxvar_well[w] > 0 ? 1: -1;
                    //const double dx1_limited = sign1 * std::min(std::abs(dxvar_well[w]),std::abs(xvar_well_old[w])*dTotalRateLimit);
                    well_state.wellSolutions()[w] = xvar_well_old[w] - dxvar_well[w];

                    switch (wells().type[w]) {
                    case INJECTOR:
                        for (int p = 0; p < np; ++p) {
                            const double comp_frac = wells().comp_frac[np*w + p];
                            //if (comp_frac > 0) {
                            well_state.wellRates()[w*np + p] = comp_frac * well_state.wellSolutions()[w];
                            //}

                        }
                        break;
                    case PRODUCER:
                        for (int p = 0; p < np; ++p) {
                            well_state.wellRates()[w*np + p] = well_state.wellSolutions()[w] * F[p];
                        }
                        break;
                    }
                }
                    break;
                case SURFACE_RATE:
                {
                    const int sign1 = dxvar_well[w] > 0 ? 1: -1;
                    const double dx1_limited = sign1 * std::min(std::abs(dxvar_well[w]),std::abs(xvar_well_old[w])*dBHPLimit);
                    well_state.wellSolutions()[w] = xvar_well_old[w] - dx1_limited;
                    //const int sign = (dxvar_well1[w] < 0) ? -1 : 1;
                    //well_state.bhp()[w] -= sign * std::min( std::abs(dxvar_well1[w]), std::abs(well_state.bhp()[w])*dpmaxrel) ;
                    well_state.bhp()[w] = well_state.wellSolutions()[w];

                    if (wells().type[w]==PRODUCER) {

                        double F_target = 0.0;
                        for (int p = 0; p < np; ++p) {
                            F_target += wells().comp_frac[np*w + p] * F[p];
                        }
                        for (int p = 0; p < np; ++p) {
                            //std::cout << F[p] << std::endl;
                            well_state.wellRates()[np*w + p] = F[p] * target_rate /F_target;
                        }
                    } else {

                        for (int p = 0; p < np; ++p) {
                            //std::cout << wells().comp_frac[np*w + p] << " " <<distr[p] << std::endl;
                            well_state.wellRates()[w*np + p] = wells().comp_frac[np*w + p] * target_rate;
                        }
                    }


                }
                    break;
                case RESERVOIR_RATE: {
                    const int sign1 = dxvar_well[w] > 0 ? 1: -1;
                    const double dx1_limited = sign1 * std::min(std::abs(dxvar_well[w]),std::abs(xvar_well_old[w])*dBHPLimit);
                    well_state.wellSolutions()[w] = xvar_well_old[w] - dx1_limited;
                    //const int sign = (dxvar_well1[w] < 0) ? -1 : 1;
                    //well_state.bhp()[w] -= sign * std::min( std::abs(dxvar_well1[w]), std::abs(well_state.bhp()[w])*dpmaxrel) ;
                    well_state.bhp()[w] = well_state.wellSolutions()[w];
                    for (int p = 0; p < np; ++p) {
                        well_state.wellRates()[np*w + p] = F[p] * target_rate;
                    }
                }
                    break;

                }
            }





            const Opm::PhaseUsage& pu = fluid_->phaseUsage();
            //Loop over all wells
#pragma omp parallel for schedule(static)
            for (int w = 0; w < nw; ++w) {
                const WellControls* wc = wells().ctrls[w];
                const int nwc = well_controls_get_num(wc);
                //Loop over all controls until we find a THP control
                //that specifies what we need...
                //Will only update THP for wells with THP control
                for (int ctrl_index=0; ctrl_index < nwc; ++ctrl_index) {
                    if (well_controls_iget_type(wc, ctrl_index) == THP) {
                        double aqua = 0.0;
                        double liquid = 0.0;
                        double vapour = 0.0;

                        if ((*active_)[ Water ]) {
                            aqua = well_state.wellRates()[w*np + pu.phase_pos[ Water ] ];
                        }
                        if ((*active_)[ Oil ]) {
                            liquid = well_state.wellRates()[w*np + pu.phase_pos[ Oil ] ];
                        }
                        if ((*active_)[ Gas ]) {
                            vapour = well_state.wellRates()[w*np + pu.phase_pos[ Gas ] ];
                        }

                        double alq = well_controls_iget_alq(wc, ctrl_index);
                        int table_id = well_controls_iget_vfp(wc, ctrl_index);

                        const WellType& well_type = wells().type[w];
                        if (well_type == INJECTOR) {
                            double dp = wellhelpers::computeHydrostaticCorrection(
                                    wells(), w, vfp_properties_->getInj()->getTable(table_id)->getDatumDepth(),
                                    wellPerforationDensities(), gravity_);

                            well_state.thp()[w] = vfp_properties_->getInj()->thp(table_id, aqua, liquid, vapour, well_state.bhp()[w] + dp);
                        }
                        else if (well_type == PRODUCER) {
                            double dp = wellhelpers::computeHydrostaticCorrection(
                                    wells(), w, vfp_properties_->getProd()->getTable(table_id)->getDatumDepth(),
                                    wellPerforationDensities(), gravity_);

                            well_state.thp()[w] = vfp_properties_->getProd()->thp(table_id, aqua, liquid, vapour, well_state.bhp()[w] + dp, alq);
                        }
                        else {
                            OPM_THROW(std::logic_error, "Expected INJECTOR or PRODUCER well");
                        }

                        //Assume only one THP control specified for each well
                        break;
                    }
                }
            }
        }
    }





    template <class WellState>
    void
    StandardWellsDense::
    updateWellControls(WellState& xw)
    {
        if( !localWellsActive() ) return ;

        std::string modestring[4] = { "BHP", "THP", "RESERVOIR_RATE", "SURFACE_RATE" };
        // Find, for each well, if any constraints are broken. If so,
        // switch control to first broken constraint.
        const int np = wells().number_of_phases;
        const int nw = wells().number_of_wells;
#pragma omp parallel for schedule(dynamic)
        for (int w = 0; w < nw; ++w) {
            WellControls* wc = wells().ctrls[w];
            // The current control in the well state overrides
            // the current control set in the Wells struct, which
            // is instead treated as a default.
            int current = xw.currentControls()[w];
            // Loop over all controls except the current one, and also
            // skip any RESERVOIR_RATE controls, since we cannot
            // handle those.
            const int nwc = well_controls_get_num(wc);
            int ctrl_index = 0;
            for (; ctrl_index < nwc; ++ctrl_index) {
                if (ctrl_index == current) {
                    // This is the currently used control, so it is
                    // used as an equation. So this is not used as an
                    // inequality constraint, and therefore skipped.
                    continue;
                }
                if (wellhelpers::constraintBroken(
                        xw.bhp(), xw.thp(), xw.wellRates(),
                        w, np, wells().type[w], wc, ctrl_index)) {
                    // ctrl_index will be the index of the broken constraint after the loop.
                    break;
                }
            }
            if (ctrl_index != nwc) {
                // Constraint number ctrl_index was broken, switch to it.
                // We disregard terminal_ouput here as with it only messages
                // for wells on one process will be printed.
                std::ostringstream ss;
                ss << "Switching control mode for well " << wells().name[w]
                      << " from " << modestring[well_controls_iget_type(wc, current)]
                      << " to " << modestring[well_controls_iget_type(wc, ctrl_index)] << std::endl;
                OpmLog::info(ss.str());
                xw.currentControls()[w] = ctrl_index;
                current = xw.currentControls()[w];
                well_controls_set_current( wc, current);



                // Updating well state and primary variables if constraint is broken

                // Target values are used as initial conditions for BHP, THP, and SURFACE_RATE
                const double target = well_controls_iget_target(wc, current);
                const double* distr = well_controls_iget_distr(wc, current);
                switch (well_controls_iget_type(wc, current)) {
                case BHP:
                    xw.bhp()[w] = target;
                    break;

                case THP: {
                    double aqua = 0.0;
                    double liquid = 0.0;
                    double vapour = 0.0;

                    const Opm::PhaseUsage& pu = fluid_->phaseUsage();

                    if ((*active_)[ Water ]) {
                        aqua = xw.wellRates()[w*np + pu.phase_pos[ Water ] ];
                    }
                    if ((*active_)[ Oil ]) {
                        liquid = xw.wellRates()[w*np + pu.phase_pos[ Oil ] ];
                    }
                    if ((*active_)[ Gas ]) {
                        vapour = xw.wellRates()[w*np + pu.phase_pos[ Gas ] ];
                    }

                    const int vfp        = well_controls_iget_vfp(wc, current);
                    const double& thp    = well_controls_iget_target(wc, current);
                    const double& alq    = well_controls_iget_alq(wc, current);

                    //Set *BHP* target by calculating bhp from THP
                    const WellType& well_type = wells().type[w];

                    if (well_type == INJECTOR) {
                        double dp = wellhelpers::computeHydrostaticCorrection(
                                    wells(), w, vfp_properties_->getInj()->getTable(vfp)->getDatumDepth(),
                                    wellPerforationDensities(), gravity_);

                        xw.bhp()[w] = vfp_properties_->getInj()->bhp(vfp, aqua, liquid, vapour, thp) - dp;
                    }
                    else if (well_type == PRODUCER) {
                        double dp = wellhelpers::computeHydrostaticCorrection(
                                    wells(), w, vfp_properties_->getProd()->getTable(vfp)->getDatumDepth(),
                                    wellPerforationDensities(), gravity_);

                        xw.bhp()[w] = vfp_properties_->getProd()->bhp(vfp, aqua, liquid, vapour, thp, alq) - dp;
                    }
                    else {
                        OPM_THROW(std::logic_error, "Expected PRODUCER or INJECTOR type of well");
                    }
                    break;
                }

                case RESERVOIR_RATE:
                    // No direct change to any observable quantity at
                    // surface condition.  In this case, use existing
                    // flow rates as initial conditions as reservoir
                    // rate acts only in aggregate.
                    break;

                case SURFACE_RATE:
                    // assign target value as initial guess for injectors and
                    // single phase producers (orat, grat, wrat)
                    const WellType& well_type = wells().type[w];
                    if (well_type == INJECTOR) {
                        for (int phase = 0; phase < np; ++phase) {
                            const double& compi = wells().comp_frac[np * w + phase];
                            //if (compi > 0.0) {
                                xw.wellRates()[np*w + phase] = target * compi;
                            //}
                        }
                    } else if (well_type == PRODUCER) {

                        // only set target as initial rates for single phase
                        // producers. (orat, grat and wrat, and not lrat)
                        // lrat will result in numPhasesWithTargetsUnderThisControl == 2
                        int numPhasesWithTargetsUnderThisControl = 0;
                        for (int phase = 0; phase < np; ++phase) {
                            if (distr[phase] > 0.0) {
                                numPhasesWithTargetsUnderThisControl += 1;
                            }
                        }
                        for (int phase = 0; phase < np; ++phase) {
                            if (distr[phase] > 0.0 && numPhasesWithTargetsUnderThisControl < 2 ) {
                                xw.wellRates()[np*w + phase] = target * distr[phase];
                            }
                        }
                    } else {
                        OPM_THROW(std::logic_error, "Expected PRODUCER or INJECTOR type of well");
                    }


                    break;
                }
                std::vector<double> g = {1,1,0.01};
                if (well_controls_iget_type(wc, current) == RESERVOIR_RATE) {
                    const double* distr = well_controls_iget_distr(wc, current);
                    for (int phase = 0; phase < np; ++phase) {
                        g[phase] = distr[phase];
                    }
                }
                switch (well_controls_iget_type(wc, current)) {
                case BHP:
                {
                    const WellType& well_type = wells().type[w];
                    xw.wellSolutions()[w] = 0.0;
                    if (well_type == INJECTOR) {
                        for (int p = 0; p < np; ++p)  {
                            xw.wellSolutions()[w] += xw.wellRates()[np*w + p] * wells().comp_frac[np*w + p];
                        }
                    } else {
                        for (int p = 0; p < np; ++p)  {
                            xw.wellSolutions()[w] += g[p] * xw.wellRates()[np*w + p];
                        }
                    }
                }
                    break;


                case RESERVOIR_RATE: // Intentional fall-through
                case SURFACE_RATE:
                {
                    xw.wellSolutions()[w] = xw.bhp()[w];
                }
                    break;
                }

                double tot_well_rate = 0.0;
                for (int p = 0; p < np; ++p)  {
                    tot_well_rate += g[p] * xw.wellRates()[np*w + p];
                }
                if(std::abs(tot_well_rate) > 0) {
                    xw.wellSolutions()[nw + w] = g[Water] * xw.wellRates()[np*w + Water] / tot_well_rate; //wells->comp_frac[np*w + Water]; // Water;
                    xw.wellSolutions()[2*nw + w] = g[Gas] * xw.wellRates()[np*w + Gas] / tot_well_rate ; //wells->comp_frac[np*w + Gas]; //Gas
                } else {
                    //xw.wellSolutions()[nw + w] =  wells().comp_frac[np*w + Water];
                    //xw.wellSolutions()[2 * nw + w] =  wells().comp_frac[np*w + Gas];
                }
            }
        }
    }





    template <class SolutionState>
    void
    StandardWellsDense::
    addWellFluxEq(const std::vector<ADB>& cq_s,
                  const SolutionState& state,
                  const double dt,
                  LinearisedBlackoilResidual& residual)
    {
        if( !localWellsActive() )
        {
            // If there are no wells in the subdomain of the proces then
            // cq_s has zero size and will cause a segmentation fault below.
            return;
        }

        const int np = wells().number_of_phases;
        const int nw = wells().number_of_wells;


        double volume = 0.002831684659200; // 0.1 cu ft;
        const Vector vol_dt = Vector::Constant(nw,volume/dt);
        std::vector<ADB> F = wellVolumeFractions(state);
        //std::cout << F0_[0] << std::endl;
        //std::cout << F[0] << std::endl;

        ADB qs = state.qs;
        for (int phase = 0; phase < np; ++phase) {
            qs -= superset(wellOps().p2w * cq_s[phase], Span(nw, 1, phase*nw), nw*np);
            qs += superset((F[phase]-F0_[phase]) * vol_dt, Span(nw,1,phase*nw), nw*np);
        }


        residual.well_flux_eq = qs;
        //std::cout << "etter dense " <<qs << std::endl;
    }





    template <class SolutionState, class WellState>
    void
    StandardWellsDense::addWellControlEq(const SolutionState& state,
                                    const WellState& xw,
                                    const Vector& aliveWells,
                                    LinearisedBlackoilResidual& residual)
    {
        if( ! localWellsActive() ) return;

        const int np = wells().number_of_phases;
        const int nw = wells().number_of_wells;

        ADB aqua   = ADB::constant(Vector::Zero(nw));
        ADB liquid = ADB::constant(Vector::Zero(nw));
        ADB vapour = ADB::constant(Vector::Zero(nw));

        if ((*active_)[Water]) {
            aqua += subset(state.qs, Span(nw, 1, BlackoilPhases::Aqua*nw));
        }
        if ((*active_)[Oil]) {
            liquid += subset(state.qs, Span(nw, 1, BlackoilPhases::Liquid*nw));
        }
        if ((*active_)[Gas]) {
            vapour += subset(state.qs, Span(nw, 1, BlackoilPhases::Vapour*nw));
        }

        //THP calculation variables
        std::vector<int> inj_table_id(nw, -1);
        std::vector<int> prod_table_id(nw, -1);
        Vector thp_inj_target_v = Vector::Zero(nw);
        Vector thp_prod_target_v = Vector::Zero(nw);
        Vector alq_v = Vector::Zero(nw);

        //Hydrostatic correction variables
        Vector rho_v = Vector::Zero(nw);
        Vector vfp_ref_depth_v = Vector::Zero(nw);

        //Target vars
        Vector bhp_targets  = Vector::Zero(nw);
        Vector rate_targets = Vector::Zero(nw);
        Eigen::SparseMatrix<double> rate_distr(nw, np*nw);

        //Selection variables
        std::vector<int> bhp_elems;
        std::vector<int> thp_inj_elems;
        std::vector<int> thp_prod_elems;
        std::vector<int> rate_elems;

        //Run through all wells to calculate BHP/RATE targets
        //and gather info about current control
        for (int w = 0; w < nw; ++w) {
            auto wc = wells().ctrls[w];

            // The current control in the well state overrides
            // the current control set in the Wells struct, which
            // is instead treated as a default.
            const int current = xw.currentControls()[w];

            switch (well_controls_iget_type(wc, current)) {
            case BHP:
            {
                bhp_elems.push_back(w);
                bhp_targets(w)  = well_controls_iget_target(wc, current);
                rate_targets(w) = -1e100;
            }
            break;

            case THP:
            {
                const int perf = wells().well_connpos[w];
                rho_v[w] = wellPerforationDensities()[perf];

                const int table_id = well_controls_iget_vfp(wc, current);
                const double target = well_controls_iget_target(wc, current);

                const WellType& well_type = wells().type[w];
                if (well_type == INJECTOR) {
                    inj_table_id[w]  = table_id;
                    thp_inj_target_v[w] = target;
                    alq_v[w]     = -1e100;

                    vfp_ref_depth_v[w] = vfp_properties_->getInj()->getTable(table_id)->getDatumDepth();

                    thp_inj_elems.push_back(w);
                }
                else if (well_type == PRODUCER) {
                    prod_table_id[w]  = table_id;
                    thp_prod_target_v[w] = target;
                    alq_v[w]      = well_controls_iget_alq(wc, current);

                    vfp_ref_depth_v[w] =  vfp_properties_->getProd()->getTable(table_id)->getDatumDepth();

                    thp_prod_elems.push_back(w);
                }
                else {
                    OPM_THROW(std::logic_error, "Expected INJECTOR or PRODUCER type well");
                }
                bhp_targets(w)  = -1e100;
                rate_targets(w) = -1e100;
            }
            break;

            case RESERVOIR_RATE: // Intentional fall-through
            case SURFACE_RATE:
            {
                rate_elems.push_back(w);
                // RESERVOIR and SURFACE rates look the same, from a
                // high-level point of view, in the system of
                // simultaneous linear equations.

                const double* const distr =
                    well_controls_iget_distr(wc, current);

                for (int p = 0; p < np; ++p) {
                    rate_distr.insert(w, p*nw + w) = distr[p];
                }

                bhp_targets(w)  = -1.0e100;
                rate_targets(w) = well_controls_iget_target(wc, current);
            }
            break;
            }
        }

        //Calculate BHP target from THP
        const ADB thp_inj_target = ADB::constant(thp_inj_target_v);
        const ADB thp_prod_target = ADB::constant(thp_prod_target_v);
        const ADB alq = ADB::constant(alq_v);
        const ADB bhp_from_thp_inj = vfp_properties_->getInj()->bhp(inj_table_id, aqua, liquid, vapour, thp_inj_target);
        const ADB bhp_from_thp_prod = vfp_properties_->getProd()->bhp(prod_table_id, aqua, liquid, vapour, thp_prod_target, alq);

        //Perform hydrostatic correction to computed targets
        const Vector dp_v = wellhelpers::computeHydrostaticCorrection(wells(), vfp_ref_depth_v, wellPerforationDensities(), gravity_);
        const ADB dp = ADB::constant(dp_v);
        const ADB dp_inj = superset(subset(dp, thp_inj_elems), thp_inj_elems, nw);
        const ADB dp_prod = superset(subset(dp, thp_prod_elems), thp_prod_elems, nw);

        //Calculate residuals
        const ADB thp_inj_residual = state.bhp - bhp_from_thp_inj + dp_inj;
        const ADB thp_prod_residual = state.bhp - bhp_from_thp_prod + dp_prod;
        const ADB bhp_residual = state.bhp - bhp_targets;
        const ADB rate_residual = rate_distr * state.qs - rate_targets;

        //Select the right residual for each well
        residual.well_eq = superset(subset(bhp_residual, bhp_elems), bhp_elems, nw) +
                superset(subset(thp_inj_residual, thp_inj_elems), thp_inj_elems, nw) +
                superset(subset(thp_prod_residual, thp_prod_elems), thp_prod_elems, nw) +
                superset(subset(rate_residual, rate_elems), rate_elems, nw);

        // For wells that are dead (not flowing), and therefore not communicating
        // with the reservoir, we set the equation to be equal to the well's total
        // flow. This will be a solution only if the target rate is also zero.
        Eigen::SparseMatrix<double> rate_summer(nw, np*nw);
        for (int w = 0; w < nw; ++w) {
            for (int phase = 0; phase < np; ++phase) {
                rate_summer.insert(w, phase*nw + w) = 1.0;
            }
        }
        Selector<double> alive_selector(aliveWells, Selector<double>::NotEqualZero);
        residual.well_eq = alive_selector.select(residual.well_eq, rate_summer * state.qs);
        // OPM_AD_DUMP(residual_.well_eq);
    }





    template <class SolutionState, class WellState>
    void
    StandardWellsDense::computeWellPotentials(const std::vector<ADB>& mob_perfcells,
                                         const std::vector<ADB>& b_perfcells,
                                         SolutionState& state0,
                                         WellState& well_state)
    {
        const int nw = wells().number_of_wells;
        const int np = wells().number_of_phases;
        const Opm::PhaseUsage& pu = fluid_->phaseUsage();

        Vector bhps = Vector::Zero(nw);
        for (int w = 0; w < nw; ++w) {
            const WellControls* ctrl = wells().ctrls[w];
            const int nwc = well_controls_get_num(ctrl);
            //Loop over all controls until we find a BHP control
            //or a THP control that specifies what we need.
            //Pick the value that gives the most restrictive flow
            for (int ctrl_index=0; ctrl_index < nwc; ++ctrl_index) {

                if (well_controls_iget_type(ctrl, ctrl_index) == BHP) {
                    bhps[w] = well_controls_iget_target(ctrl, ctrl_index);
                }

                if (well_controls_iget_type(ctrl, ctrl_index) == THP) {
                    double aqua = 0.0;
                    double liquid = 0.0;
                    double vapour = 0.0;

                    if ((*active_)[ Water ]) {
                        aqua = well_state.wellRates()[w*np + pu.phase_pos[ Water ] ];
                    }
                    if ((*active_)[ Oil ]) {
                        liquid = well_state.wellRates()[w*np + pu.phase_pos[ Oil ] ];
                    }
                    if ((*active_)[ Gas ]) {
                        vapour = well_state.wellRates()[w*np + pu.phase_pos[ Gas ] ];
                    }

                    const int vfp        = well_controls_iget_vfp(ctrl, ctrl_index);
                    const double& thp    = well_controls_iget_target(ctrl, ctrl_index);
                    const double& alq    = well_controls_iget_alq(ctrl, ctrl_index);

                    //Set *BHP* target by calculating bhp from THP
                    const WellType& well_type = wells().type[w];

                    if (well_type == INJECTOR) {
                        double dp = wellhelpers::computeHydrostaticCorrection(
                                    wells(), w, vfp_properties_->getInj()->getTable(vfp)->getDatumDepth(),
                                    wellPerforationDensities(), gravity_);
                        const double bhp = vfp_properties_->getInj()->bhp(vfp, aqua, liquid, vapour, thp) - dp;
                        // apply the strictest of the bhp controlls i.e. smallest bhp for injectors
                        if ( bhp < bhps[w]) {
                            bhps[w] = bhp;
                        }
                    }
                    else if (well_type == PRODUCER) {
                        double dp = wellhelpers::computeHydrostaticCorrection(
                                    wells(), w, vfp_properties_->getProd()->getTable(vfp)->getDatumDepth(),
                                    wellPerforationDensities(), gravity_);

                        const double bhp = vfp_properties_->getProd()->bhp(vfp, aqua, liquid, vapour, thp, alq) - dp;
                        // apply the strictest of the bhp controlls i.e. largest bhp for producers
                        if ( bhp > bhps[w]) {
                            bhps[w] = bhp;
                        }
                    }
                    else {
                        OPM_THROW(std::logic_error, "Expected PRODUCER or INJECTOR type of well");
                    }
                }
            }

        }

        // use bhp limit from control
        state0.bhp = ADB::constant(bhps);

        // compute well potentials
        Vector aliveWells;
        std::vector<ADB> well_potentials;
        computeWellFlux(state0, mob_perfcells,  b_perfcells, aliveWells, well_potentials);

        // store well potentials in the well state
        // transform to a single vector instead of separate vectors pr phase
        const int nperf = wells().well_connpos[nw];
        Vector cq = superset(well_potentials[0].value(), Span(nperf, np, 0), nperf*np);
        for (int phase = 1; phase < np; ++phase) {
            cq += superset(well_potentials[phase].value(), Span(nperf, np, phase), nperf*np);
        }
        well_state.wellPotentials().assign(cq.data(), cq.data() + nperf*np);
    }





    void
    StandardWellsDense::variableStateWellIndices(std::vector<int>& indices,
                                            int& next) const
    {
        indices[Qs] = next++;
        //indices[Bhp] = next++;
    }




    template <class SolutionState, class WellState>
    void
    StandardWellsDense::
    variableStateExtractWellsVars(const std::vector<int>& indices,
                                  std::vector<ADB>& vars,
                                  SolutionState& state,
                                  WellState &xw) const
    {

        state.wellVariables = std::move(vars[indices[Qs]]);
        //std::cout << "state.wellVariables " << state.wellVariables.value() << std::endl;

        const int nw = wells().number_of_wells;
        const int np = wells().number_of_phases;
        const Opm::PhaseUsage& pu = fluid_->phaseUsage();

        Vector isBHPControlled = Vector::Zero(nw);
        Vector isReservoirRateControlled = Vector::Zero(nw);
        Vector isSurfaceRateControlled = Vector::Zero(nw);

        Vector isInjector = Vector::Zero(nw);
        for (int w = 0; w < nw; ++w) {
            const WellControls* wc = wells().ctrls[w];
            // The current control in the well state overrides
            // the current control set in the Wells struct, which
            // is instead treated as a default.
            const int current = xw.currentControls()[w];

            switch (well_controls_iget_type(wc, current)) {

            case BHP:
                isBHPControlled[w] = 1;
                break;

            case RESERVOIR_RATE:
                isReservoirRateControlled[w] = 1;

                break;
            case SURFACE_RATE:
                isSurfaceRateControlled[w] = 1;
                break;
            }

            if (wells().type[w] == INJECTOR) {
                isInjector[w] = 1;
            }
        }

        Vector ones = Vector::Constant(nw,1);


        const ADB& xvar_well1 = subset(state.wellVariables,Span(nw,1,0));
        //std::cout << xvar_well1.value() << std::endl;
        const V bhp = Eigen::Map<const V>(& xw.bhp()[0], xw.bhp().size());

        state.bhp = isBHPControlled * bhp + (1-isBHPControlled) * xvar_well1;

        const ADB Fw = subset(state.wellVariables,Span(nw,1,nw));
        const ADB Fg = subset(state.wellVariables,Span(nw,1,nw*2));
        const ADB Fo = ones - Fw - Fg;

        const DataBlock compi = Eigen::Map<const DataBlock>(wells().comp_frac, nw, np);

        V target_rates = V::Zero(nw);
        V isNoFlow = V::Zero(nw);
        std::vector<V> distri(np);
        for (int p = 0; p < np; ++p) {
            distri[p] = V::Zero(nw);
        }


        for (int w = 0; w < nw; ++w) {
            const WellControls* wc = wells().ctrls[w];
            // The current control in the well state overrides
            // the current control set in the Wells struct, which
            // is instead treated as a default.
            const int current = xw.currentControls()[w];
            const double* distr = well_controls_iget_distr(wc, current);


            target_rates[w] = well_controls_iget_target(wc, current);
            if (target_rates[w]==0)
                isNoFlow[w] = 1;
            for (int p = 0; p < np; ++p) {
                const double comp_frac = wells().comp_frac[np*w + p];
                distri[p][w] = distr[p];
                //std::cout << "distr " << distr[p] << " comp_frac " << comp_frac << std::endl;
            }

        }

        std::vector<ADB> F(np, ADB::null());
        F[Water] = Fw;
        F[Oil] = Fo;
        F[Gas] = Fg; // / Vector::Constant(nw,0.01);

        std::vector<Vector> g(np,Vector::Constant(nw,1.0));
        g[Gas] = Vector::Constant(nw,0.01);
//        g[Water] = Vector::Constant(nw,1.0);
//        g[Oil] = Vector::Constant(nw,1.0);
        for (int p = 0; p < np; ++p) {
            const V tmp = (ones - isReservoirRateControlled) * g[p] + isReservoirRateControlled * distri[p];
            F[p] = F[p] / tmp;
        }

        ADB targetF = compi.col(pu.phase_pos[0])*F[0];
        for (int p = 1; p < np; ++p) {
            targetF += compi.col(pu.phase_pos[p])*F[p];
        }


        ADB sumF = distri[0] * F[0];
        for (int p = 1; p < np; ++p) {
            sumF += distri[p] * (F[p]);
        }
        Selector<double> noTargetFSelector(targetF.value(), Selector<double>::Zero);
        Selector<double> noSumFSelector(sumF.value(), Selector<double>::Zero);


        //std::cout << target_rates << std::endl;
        ADB qs = ADB::constant(V::Zero(nw*np));
        for (int p = 0; p < np; ++p) {

            const int pos = pu.phase_pos[p];
            ADB Q = isSurfaceRateControlled * ADB::constant(compi.col(pos) * target_rates);
            //DUMPVAL(Q);
            Q += isBHPControlled * isInjector * xvar_well1 * ADB::constant(compi.col(pos));
            //DUMPVAL(Q);
            Q += isBHPControlled * (ones - isInjector) * xvar_well1 * F[p];
            //DUMPVAL(Q);
            Q += (ones - isInjector) * isSurfaceRateControlled  * noTargetFSelector.select( ADB::constant(V::Zero(nw)), (ones - compi.col(pos)) * target_rates * F[p] / targetF);

            //Q += isNoFlow * -1e-24 * F[p];

            //Q += isSurfaceRateControlled * isInjector * -1e-12 * F[p] * (ones - compi.col(pos));

            Q += isReservoirRateControlled * F[p] * target_rates;
            //DUMPVAL(Q);
            qs += superset(Q, Span(nw, 1, p*nw), nw*np);
        }
        state.qs = qs;
        //if (isNoFlow.sum() > 0)
            //std::cout << state.qs << std::endl;

        // Bhp.
        //state.bhp =
        //state.qs =
    }





    std::vector<int>
    StandardWellsDense::variableWellStateIndices() const
    {
        // Black oil model standard is 5 equation.
        // For the pure well solve, only the well equations are picked.
        std::vector<int> indices(5, -1);
        int next = 0;

        variableStateWellIndices(indices, next);

        assert(next == 1);
        return indices;
    }





    template <class WellState>
    void
    StandardWellsDense::variableWellStateInitials(const WellState& xw,
                                             std::vector<Vector>& vars0) const
    {
        // Initial well rates.
        if ( localWellsActive() )
        {
            // Need to reshuffle well rates, from phase running fastest
            // to wells running fastest.
            const int nw = wells().number_of_wells;
            const int np = wells().number_of_phases;

            // The transpose() below switches the ordering.
            //const DataBlock wrates = Eigen::Map<const DataBlock>(& xw.wellSolutions()[0], nw, np).transpose();
            //const Vector qs = Eigen::Map<const V>(wrates.data(), nw*np);
            const Vector qs = Eigen::Map<const V>(& xw.wellSolutions()[0], xw.wellSolutions().size());
            vars0.push_back(qs);

            // Initial well bottom-hole pressure.
//            assert (not xw.bhp().empty());
//            const Vector bhp = Eigen::Map<const V>(& xw.bhp()[0], xw.bhp().size());
//            vars0.push_back(bhp);

        }
        else
        {
            // push null states for qs and bhp
            vars0.push_back(V());
            //vars0.push_back(V());
        }
    }

    template <class SolutionState>
    std::vector<AutoDiffBlock<double>> StandardWellsDense::wellVolumeFractions(const SolutionState& state) const
    {
        const int np = wells().number_of_phases;
        const int nw = wells().number_of_wells;
        std::vector<ADB> F(np,ADB::null());
        F[Water] = subset(state.wellVariables,Span(nw,1,nw));
        F[Gas] = subset(state.wellVariables,Span(nw,1,2 * nw));
        F[Oil] = Vector::Constant(nw,1.0) - F[Water] - F[Gas];
        return F;
    }


    template <class SolutionState>
    void StandardWellsDense::computeAccumWells(const SolutionState& state) {
        F0_ = wellVolumeFractions(state);
    }



    void
    StandardWellsDense::setStoreWellPerforationFluxesFlag(const bool store_fluxes)
    {
        store_well_perforation_fluxes_ = store_fluxes;
    }





    const StandardWellsDense::Vector&
    StandardWellsDense::getStoredWellPerforationFluxes() const
    {
        assert(store_well_perforation_fluxes_);
        return well_perforation_fluxes_;
    }






   template<class WellState>
   void
   StandardWellsDense::
   updateListEconLimited(ScheduleConstPtr schedule,
                         const int current_step,
                         const Wells* wells_struct,
                         const WellState& well_state,
                         DynamicListEconLimited& list_econ_limited) const
   {
       const int nw = wells_struct->number_of_wells;

       for (int w = 0; w < nw; ++w) {
           // flag to check if the mim oil/gas rate limit is violated
           bool rate_limit_violated = false;
           const std::string& well_name = wells_struct->name[w];
           const Well* well_ecl = schedule->getWell(well_name);
           const WellEconProductionLimits& econ_production_limits = well_ecl->getEconProductionLimits(current_step);

           // economic limits only apply for production wells.
           if (wells_struct->type[w] != PRODUCER) {
               continue;
           }

           // if no limit is effective here, then continue to the next well
           if ( !econ_production_limits.onAnyEffectiveLimit() ) {
               continue;
           }
           // for the moment, we only handle rate limits, not handling potential limits
           // the potential limits should not be difficult to add
           const WellEcon::QuantityLimitEnum& quantity_limit = econ_production_limits.quantityLimit();
           if (quantity_limit == WellEcon::POTN) {
               const std::string msg = std::string("POTN limit for well ") + well_name + std::string(" is not supported for the moment. \n")
                                     + std::string("All the limits will be evaluated based on RATE. ");
               OpmLog::warning("NOT_SUPPORTING_POTN", msg);
           }

           const WellMapType& well_map = well_state.wellMap();
           const typename WellMapType::const_iterator i_well = well_map.find(well_name);
           assert(i_well != well_map.end()); // should always be found?
           const WellMapEntryType& map_entry = i_well->second;
           const int well_number = map_entry[0];

           if (econ_production_limits.onAnyRateLimit()) {
               rate_limit_violated = checkRateEconLimits(econ_production_limits, well_state, well_number);
           }

           if (rate_limit_violated) {
               if (econ_production_limits.endRun()) {
                   const std::string warning_message = std::string("ending run after well closed due to economic limits is not supported yet \n")
                                                     + std::string("the program will keep running after ") + well_name + std::string(" is closed");
                   OpmLog::warning("NOT_SUPPORTING_ENDRUN", warning_message);
               }

               if (econ_production_limits.validFollowonWell()) {
                   OpmLog::warning("NOT_SUPPORTING_FOLLOWONWELL", "opening following on well after well closed is not supported yet");
               }

               if (well_ecl->getAutomaticShutIn()) {
                   list_econ_limited.addShutWell(well_name);
                   const std::string msg = std::string("well ") + well_name + std::string(" will be shut in due to economic limit");
                   OpmLog::info(msg);
               } else {
                   list_econ_limited.addStoppedWell(well_name);
                   const std::string msg = std::string("well ") + well_name + std::string(" will be stopped due to economic limit");
                   OpmLog::info(msg);
               }
               // the well is closed, not need to check other limits
               continue;
           }

           // checking for ratio related limits, mostly all kinds of ratio.
           bool ratio_limits_violated = false;
           RatioCheckTuple ratio_check_return;

           if (econ_production_limits.onAnyRatioLimit()) {
               ratio_check_return = checkRatioEconLimits(econ_production_limits, well_state, map_entry);
               ratio_limits_violated = std::get<0>(ratio_check_return);
           }

           if (ratio_limits_violated) {
               const bool last_connection = std::get<1>(ratio_check_return);
               const int worst_offending_connection = std::get<2>(ratio_check_return);

               const int perf_start = map_entry[1];

               assert((worst_offending_connection >= 0) && (worst_offending_connection <  map_entry[2]));

               const int cell_worst_offending_connection = wells_struct->well_cells[perf_start + worst_offending_connection];
               list_econ_limited.addClosedConnectionsForWell(well_name, cell_worst_offending_connection);
               const std::string msg = std::string("Connection ") + std::to_string(worst_offending_connection) + std::string(" for well ")
                                     + well_name + std::string(" will be closed due to economic limit");
               OpmLog::info(msg);

               if (last_connection) {
                   list_econ_limited.addShutWell(well_name);
                   const std::string msg2 = well_name + std::string(" will be shut due to the last connection closed");
                   OpmLog::info(msg2);
               }
           }

       }
   }





    template <class WellState>
    bool
    StandardWellsDense::
    checkRateEconLimits(const WellEconProductionLimits& econ_production_limits,
                        const WellState& well_state,
                        const int well_number) const
    {
        const Opm::PhaseUsage& pu = fluid_->phaseUsage();
        const int np = well_state.numPhases();

        if (econ_production_limits.onMinOilRate()) {
            assert((*active_)[Oil]);
            const double oil_rate = well_state.wellRates()[well_number * np + pu.phase_pos[ Oil ] ];
            const double min_oil_rate = econ_production_limits.minOilRate();
            if (std::abs(oil_rate) < min_oil_rate) {
                return true;
            }
        }

        if (econ_production_limits.onMinGasRate() ) {
            assert((*active_)[Gas]);
            const double gas_rate = well_state.wellRates()[well_number * np + pu.phase_pos[ Gas ] ];
            const double min_gas_rate = econ_production_limits.minGasRate();
            if (std::abs(gas_rate) < min_gas_rate) {
                return true;
            }
        }

        if (econ_production_limits.onMinLiquidRate() ) {
            assert((*active_)[Oil]);
            assert((*active_)[Water]);
            const double oil_rate = well_state.wellRates()[well_number * np + pu.phase_pos[ Oil ] ];
            const double water_rate = well_state.wellRates()[well_number * np + pu.phase_pos[ Water ] ];
            const double liquid_rate = oil_rate + water_rate;
            const double min_liquid_rate = econ_production_limits.minLiquidRate();
            if (std::abs(liquid_rate) < min_liquid_rate) {
                return true;
            }
        }

        if (econ_production_limits.onMinReservoirFluidRate()) {
            OpmLog::warning("NOT_SUPPORTING_MIN_RESERVOIR_FLUID_RATE", "Minimum reservoir fluid production rate limit is not supported yet");
        }

        return false;
    }





    template <class WellState>
    StandardWellsDense::RatioCheckTuple
    StandardWellsDense::
    checkRatioEconLimits(const WellEconProductionLimits& econ_production_limits,
                         const WellState& well_state,
                         const WellMapEntryType& map_entry) const
    {
        // TODO: not sure how to define the worst-offending connection when more than one
        //       ratio related limit is violated.
        //       The defintion used here is that we define the violation extent based on the
        //       ratio between the value and the corresponding limit.
        //       For each violated limit, we decide the worst-offending connection separately.
        //       Among the worst-offending connections, we use the one has the biggest violation
        //       extent.

        bool any_limit_violated = false;
        bool last_connection = false;
        int worst_offending_connection = INVALIDCONNECTION;
        double violation_extent = -1.0;

        if (econ_production_limits.onMaxWaterCut()) {
            const RatioCheckTuple water_cut_return = checkMaxWaterCutLimit(econ_production_limits, well_state, map_entry);
            bool water_cut_violated = std::get<0>(water_cut_return);
            if (water_cut_violated) {
                any_limit_violated = true;
                const double violation_extent_water_cut = std::get<3>(water_cut_return);
                if (violation_extent_water_cut > violation_extent) {
                    violation_extent = violation_extent_water_cut;
                    worst_offending_connection = std::get<2>(water_cut_return);
                    last_connection = std::get<1>(water_cut_return);
                }
            }
        }

        if (econ_production_limits.onMaxGasOilRatio()) {
            OpmLog::warning("NOT_SUPPORTING_MAX_GOR", "the support for max Gas-Oil ratio is not implemented yet!");
        }

        if (econ_production_limits.onMaxWaterGasRatio()) {
            OpmLog::warning("NOT_SUPPORTING_MAX_WGR", "the support for max Water-Gas ratio is not implemented yet!");
        }

        if (econ_production_limits.onMaxGasLiquidRatio()) {
            OpmLog::warning("NOT_SUPPORTING_MAX_GLR", "the support for max Gas-Liquid ratio is not implemented yet!");
        }

        if (any_limit_violated) {
            assert(worst_offending_connection >=0);
            assert(violation_extent > 1.);
        }

        return std::make_tuple(any_limit_violated, last_connection, worst_offending_connection, violation_extent);
    }





    template <class WellState>
    StandardWellsDense::RatioCheckTuple
    StandardWellsDense::
    checkMaxWaterCutLimit(const WellEconProductionLimits& econ_production_limits,
                          const WellState& well_state,
                          const WellMapEntryType& map_entry) const
    {
        bool water_cut_limit_violated = false;
        int worst_offending_connection = INVALIDCONNECTION;
        bool last_connection = false;
        double violation_extent = -1.0;

        const int np = well_state.numPhases();
        const Opm::PhaseUsage& pu = fluid_->phaseUsage();
        const int well_number = map_entry[0];

        assert((*active_)[Oil]);
        assert((*active_)[Water]);

        const double oil_rate = well_state.wellRates()[well_number * np + pu.phase_pos[ Oil ] ];
        const double water_rate = well_state.wellRates()[well_number * np + pu.phase_pos[ Water ] ];
        const double liquid_rate = oil_rate + water_rate;
        double water_cut;
        if (std::abs(liquid_rate) != 0.) {
            water_cut = water_rate / liquid_rate;
        } else {
            water_cut = 0.0;
        }

        const double max_water_cut_limit = econ_production_limits.maxWaterCut();
        if (water_cut > max_water_cut_limit) {
            water_cut_limit_violated = true;
        }

        if (water_cut_limit_violated) {
            // need to handle the worst_offending_connection
            const int perf_start = map_entry[1];
            const int perf_number = map_entry[2];

            std::vector<double> water_cut_perf(perf_number);
            for (int perf = 0; perf < perf_number; ++perf) {
                const int i_perf = perf_start + perf;
                const double oil_perf_rate = well_state.perfPhaseRates()[i_perf * np + pu.phase_pos[ Oil ] ];
                const double water_perf_rate = well_state.perfPhaseRates()[i_perf * np + pu.phase_pos[ Water ] ];
                const double liquid_perf_rate = oil_perf_rate + water_perf_rate;
                if (std::abs(liquid_perf_rate) != 0.) {
                    water_cut_perf[perf] = water_perf_rate / liquid_perf_rate;
                } else {
                    water_cut_perf[perf] = 0.;
                }
            }

            last_connection = (perf_number == 1);
            if (last_connection) {
                worst_offending_connection = 0;
                violation_extent = water_cut_perf[0] / max_water_cut_limit;
                return std::make_tuple(water_cut_limit_violated, last_connection, worst_offending_connection, violation_extent);
            }

            double max_water_cut_perf = 0.;
            for (int perf = 0; perf < perf_number; ++perf) {
                if (water_cut_perf[perf] > max_water_cut_perf) {
                    worst_offending_connection = perf;
                    max_water_cut_perf = water_cut_perf[perf];
                }
            }

            assert(max_water_cut_perf != 0.);
            assert((worst_offending_connection >= 0) && (worst_offending_connection < perf_number));

            violation_extent = max_water_cut_perf / max_water_cut_limit;
        }

        return std::make_tuple(water_cut_limit_violated, last_connection, worst_offending_connection, violation_extent);
    }


} // namespace Opm
