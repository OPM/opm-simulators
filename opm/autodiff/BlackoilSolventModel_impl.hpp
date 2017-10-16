/*
  Copyright 2015 IRIS AS

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

#ifndef OPM_BLACKOILSOLVENTMODEL_IMPL_HEADER_INCLUDED
#define OPM_BLACKOILSOLVENTMODEL_IMPL_HEADER_INCLUDED

#include <opm/autodiff/BlackoilSolventModel.hpp>

#include <opm/autodiff/AutoDiffBlock.hpp>
#include <opm/autodiff/AutoDiffHelpers.hpp>
#include <opm/autodiff/GridHelpers.hpp>
#include <opm/autodiff/BlackoilPropsAdFromDeck.hpp>
#include <opm/autodiff/GeoProps.hpp>
#include <opm/autodiff/WellDensitySegmented.hpp>

#include <opm/core/grid.h>
#include <opm/core/linalg/LinearSolverInterface.hpp>
#include <opm/core/linalg/ParallelIstlInformation.hpp>
#include <opm/core/props/rock/RockCompressibility.hpp>
#include <opm/common/ErrorMacros.hpp>
#include <opm/common/Exceptions.hpp>
#include <opm/parser/eclipse/Units/Units.hpp>
#include <opm/core/well_controls.h>
#include <opm/core/utility/parameters/ParameterGroup.hpp>

#include <cassert>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <limits>

namespace Opm {



    namespace detail {

        template <class PU>
        int solventPos(const PU& pu)
        {
            const int maxnp = Opm::BlackoilPhases::MaxNumPhases;
            int pos = 0;
            for (int phase = 0; phase < maxnp; ++phase) {
                if (pu.phase_used[phase]) {
                    pos++;
                }
            }

            return pos;
        }

    } // namespace detail



    template <class Grid>
    BlackoilSolventModel<Grid>::BlackoilSolventModel(const typename Base::ModelParameters&   param,
                                                     const Grid&                             grid,
                                                     const BlackoilPropsAdFromDeck&         fluid,
                                                     const DerivedGeology&                   geo,
                                                     const RockCompressibility*              rock_comp_props,
                                                     const SolventPropsAdFromDeck&           solvent_props,
                                                     const StandardWellsSolvent&             well_model,
                                                     const NewtonIterationBlackoilInterface& linsolver,
                                                     std::shared_ptr< const EclipseState >   eclState,
                                                     const bool                              has_disgas,
                                                     const bool                              has_vapoil,
                                                     const bool                              terminal_output,
                                                     const bool                              has_solvent,
                                                     const bool                              is_miscible)
        : Base(param, grid, fluid, geo, rock_comp_props, well_model, linsolver,
               eclState, has_disgas, has_vapoil, terminal_output),
          has_solvent_(has_solvent),
          solvent_pos_(detail::solventPos(fluid.phaseUsage())),
          solvent_props_(solvent_props),
          is_miscible_(is_miscible)

    {
        if (has_solvent_) {

            // If deck has solvent, residual_ should contain solvent equation.
            sd_.rq.resize(fluid_.numPhases() + 1);
            residual_.material_balance_eq.resize(fluid_.numPhases() + 1, ADB::null());
            Base::material_name_.push_back("Solvent");
            assert(solvent_pos_ == fluid_.numPhases());
            residual_.matbalscale.resize(fluid_.numPhases() + 1, 0.0031); // use the same as gas

            wellModel().initSolvent(&solvent_props_, solvent_pos_, has_solvent_);
        }
        if (is_miscible_) {
            mu_eff_.resize(fluid_.numPhases() + 1, ADB::null());
            b_eff_.resize(fluid_.numPhases() + 1, ADB::null());
        }
    }





    template <class Grid>
    void
    BlackoilSolventModel<Grid>::makeConstantState(SolutionState& state) const
    {
        Base::makeConstantState(state);
        state.solvent_saturation = ADB::constant(state.solvent_saturation.value());
    }





    template <class Grid>
    std::vector<V>
    BlackoilSolventModel<Grid>::variableStateInitials(const ReservoirState& x,
                                                      const WellState& xw) const
    {
        std::vector<V> vars0 = Base::variableStateInitials(x, xw);
        assert(int(vars0.size()) == fluid_.numPhases() + 2);

        // Initial solvent concentration.
        if (has_solvent_) {
            const auto& solvent_saturation = x.getCellData( BlackoilState::SSOL );
            const int nc = solvent_saturation.size();
            const V ss = Eigen::Map<const V>(solvent_saturation.data() , nc);

            // This is some insanely detailed flickety flackety code;
            // Solvent belongs after other reservoir vars but before well vars.
            auto solvent_pos = vars0.begin() + fluid_.numPhases();
            assert (not solvent_saturation.empty());
            assert(solvent_pos == vars0.end() - 2);
            vars0.insert(solvent_pos, ss);
        }
        return vars0;
    }





    template <class Grid>
    std::vector<int>
    BlackoilSolventModel<Grid>::variableStateIndices() const
    {
        std::vector<int> ind = Base::variableStateIndices();
        assert(ind.size() == 5);
        if (has_solvent_) {
            ind.resize(6);
            // Solvent belongs after other reservoir vars but before well vars.
            ind[Solvent] = fluid_.numPhases();
            // Solvent is pushing back the well vars.
            ++ind[Qs];
            ++ind[Bhp];
        }
        return ind;
    }




    template <class Grid>
    typename BlackoilSolventModel<Grid>::SolutionState
    BlackoilSolventModel<Grid>::variableStateExtractVars(const ReservoirState& x,
                                                         const std::vector<int>& indices,
                                                         std::vector<ADB>& vars) const
    {
        // This is more or less a copy of the base class. Refactoring is needed in the base class
        // to avoid unnecessary copying.

        //using namespace Opm::AutoDiffGrid;
        const int nc = Opm::AutoDiffGrid::numCells(grid_);
        const Opm::PhaseUsage pu = fluid_.phaseUsage();

        SolutionState state(fluid_.numPhases());

        // Pressure.
        state.pressure = std::move(vars[indices[Pressure]]);

        // Temperature cannot be a variable at this time (only constant).
        const V temp = Eigen::Map<const V>(& x.temperature()[0], x.temperature().size());
        state.temperature = ADB::constant(temp);

        // Saturations
        {
            ADB so = ADB::constant(V::Ones(nc, 1));

            if (active_[ Water ]) {
                state.saturation[pu.phase_pos[ Water ]] = std::move(vars[indices[Sw]]);
                const ADB& sw = state.saturation[pu.phase_pos[ Water ]];
                so -= sw;
            }
            if (has_solvent_) {
                state.solvent_saturation = std::move(vars[indices[Solvent]]);
                so -= state.solvent_saturation;
            }

            if (active_[ Gas ]) {
                // Define Sg Rs and Rv in terms of xvar.
                // Xvar is only defined if gas phase is active
                const ADB& xvar = vars[indices[Xvar]];
                ADB& sg = state.saturation[ pu.phase_pos[ Gas ] ];
                sg = Base::isSg_*xvar + Base::isRv_*so;
                so -= sg;

                if (active_[ Oil ]) {
                    // RS and RV is only defined if both oil and gas phase are active.
                    state.canonical_phase_pressures = computePressures(state.pressure, state.saturation[pu.phase_pos[ Water ]], so, sg, state.solvent_saturation);
                    sd_.rsSat = fluidRsSat(state.canonical_phase_pressures[ Oil ], so , cells_);
                    if (has_disgas_) {
                        state.rs = (1-Base::isRs_)*sd_.rsSat + Base::isRs_*xvar;
                    } else {
                        state.rs = sd_.rsSat;
                    }
                    sd_.rvSat = fluidRvSat(state.canonical_phase_pressures[ Gas ], so , cells_);
                    if (has_vapoil_) {
                        state.rv = (1-Base::isRv_)*sd_.rvSat + Base::isRv_*xvar;
                    } else {
                        state.rv = sd_.rvSat;
                    }
                }
            }

            if (active_[ Oil ]) {
                // Note that so is never a primary variable.
                state.saturation[pu.phase_pos[ Oil ]] = std::move(so);
            }
        }
        // wells
        wellModel().variableStateExtractWellsVars(indices, vars, state);
        return state;

    }





    template <class Grid>
    void
    BlackoilSolventModel<Grid>::computeAccum(const SolutionState& state,
                                             const int            aix  )
    {

        if (is_miscible_) {
            computeEffectiveProperties(state);
        }
        Base::computeAccum(state, aix);

        // Compute accumulation of the solvent
        if (has_solvent_) {
            const ADB& press = state.pressure;
            const ADB& ss = state.solvent_saturation;
            const ADB pv_mult = poroMult(press); // also computed in Base::computeAccum, could be optimized.
            const Opm::PhaseUsage& pu = fluid_.phaseUsage();

            const ADB& pg = state.canonical_phase_pressures[pu.phase_pos[Gas]];
            const std::vector<PhasePresence>& cond = phaseCondition();
            sd_.rq[solvent_pos_].b = fluidReciprocFVF(Solvent, pg, state.temperature, state.rs, state.rv,cond);
            sd_.rq[solvent_pos_].accum[aix] = pv_mult * sd_.rq[solvent_pos_].b * ss;
        }
    }





    template <class Grid>
    void
    BlackoilSolventModel<Grid>::
    assembleMassBalanceEq(const SolutionState& state)
    {

        Base::assembleMassBalanceEq(state);

        if (has_solvent_) {
            residual_.material_balance_eq[ solvent_pos_ ] =
                pvdt_ * (sd_.rq[solvent_pos_].accum[1] - sd_.rq[solvent_pos_].accum[0])
                + ops_.div*sd_.rq[solvent_pos_].mflux;
        }

    }

    template <class Grid>
    void
    BlackoilSolventModel<Grid>::updateEquationsScaling()
    {
        Base::updateEquationsScaling();
        assert(MaxNumPhases + 1 == residual_.matbalscale.size());
        if (has_solvent_) {
            const ADB& temp_b = sd_.rq[solvent_pos_].b;
            ADB::V B = 1. / temp_b.value();
#if HAVE_MPI
            if ( linsolver_.parallelInformation().type() == typeid(ParallelISTLInformation) )
            {
                const ParallelISTLInformation& real_info =
                    boost::any_cast<const ParallelISTLInformation&>(linsolver_.parallelInformation());
                double B_global_sum = 0;
                real_info.computeReduction(B, Reduction::makeGlobalSumFunctor<double>(), B_global_sum);
                residual_.matbalscale[solvent_pos_] = B_global_sum / Base::global_nc_;
            }
            else
#endif
            {
                residual_.matbalscale[solvent_pos_] = B.mean();
            }
        }
    }

    template <class Grid>
    void BlackoilSolventModel<Grid>::addWellContributionToMassBalanceEq(const std::vector<ADB>& cq_s,
                                                                        const SolutionState& state,
                                                                        WellState& xw)

    {

        // Add well contributions to solvent mass balance equation

        Base::addWellContributionToMassBalanceEq(cq_s, state, xw);

        if (has_solvent_) {
            const int nperf = wells().well_connpos[wells().number_of_wells];
            const int nc = Opm::AutoDiffGrid::numCells(grid_);

            const Opm::PhaseUsage& pu = fluid_.phaseUsage();
            const ADB zero = ADB::constant(V::Zero(nc));
            const ADB& ss = state.solvent_saturation;
            const ADB& sg = (active_[ Gas ]
                             ? state.saturation[ pu.phase_pos[ Gas ] ]
                             : zero);

            const std::vector<int> well_cells(wells().well_cells, wells().well_cells + nperf);
            Selector<double> zero_selector(ss.value() + sg.value(), Selector<double>::Zero);
            ADB F_solvent = subset(zero_selector.select(ss, ss / (ss + sg)),well_cells);

            const int nw = wells().number_of_wells;
            V injectedSolventFraction = Eigen::Map<const V>(&xw.solventFraction()[0], nperf);

            V isProducer = V::Zero(nperf);
            V ones = V::Constant(nperf,1.0);
            for (int w = 0; w < nw; ++w) {
                if(wells().type[w] == PRODUCER) {
                    for (int perf = wells().well_connpos[w]; perf < wells().well_connpos[w+1]; ++perf) {
                        isProducer[perf] = 1;
                    }
                }
            }

            const ADB& rs_perfcells = subset(state.rs, well_cells);
            const ADB& rv_perfcells = subset(state.rv, well_cells);
            int gas_pos = fluid_.phaseUsage().phase_pos[Gas];
            int oil_pos = fluid_.phaseUsage().phase_pos[Oil];

            // remove contribution from the dissolved gas.
            const ADB cq_s_solvent = (isProducer * F_solvent + (ones - isProducer) * injectedSolventFraction) * (cq_s[gas_pos] - rs_perfcells * cq_s[oil_pos]) / (ones - rs_perfcells * rv_perfcells);

            // Solvent contribution to the mass balance equation is given as a fraction
            // of the gas contribution.
            residual_.material_balance_eq[solvent_pos_] -= superset(cq_s_solvent, well_cells, nc);

            // The gas contribution must be reduced accordingly for the total contribution to be
            // the same.
            residual_.material_balance_eq[gas_pos] += superset(cq_s_solvent, well_cells, nc);

        }
    }





    template <class Grid>
    void
    BlackoilSolventModel<Grid>::
    updateState(const V& dx,
                ReservoirState& reservoir_state,
                WellState& well_state)
    {
        //TODO:
        // This is basicly a copy of updateState in the Base class
        // The convergence is very sensitive to details in the appelyard process
        // and the hydrocarbonstate detection. Further changes may occur, refactoring
        // to reuse more of the base class is planned when the code mature a bit more.
        using namespace Opm::AutoDiffGrid;
        const int np = fluid_.numPhases();
        const int nc = numCells(grid_);
        const V null;
        assert(null.size() == 0);
        const V zero = V::Zero(nc);
        const V ones = V::Constant(nc,1.0);

        // Extract parts of dx corresponding to each part.
        const V dp = subset(dx, Span(nc));
        int varstart = nc;
        const V dsw = active_[Water] ? subset(dx, Span(nc, 1, varstart)) : null;
        varstart += dsw.size();

        const V dxvar = active_[Gas] ? subset(dx, Span(nc, 1, varstart)): null;
        varstart += dxvar.size();

        const V dss = has_solvent_ ? subset(dx, Span(nc, 1, varstart)) : null;
        varstart += dss.size();

        // Extract well parts np phase rates + bhp
        const V dwells = subset(dx, Span(wellModel().numWellVars(), 1, varstart));
        varstart += dwells.size();

        assert(varstart == dx.size());

        // Pressure update.
        const double dpmaxrel = dpMaxRel();
        const V p_old = Eigen::Map<const V>(&reservoir_state.pressure()[0], nc, 1);
        const V absdpmax = dpmaxrel*p_old.abs();
        const V dp_limited = sign(dp) * dp.abs().min(absdpmax);
        const V p = (p_old - dp_limited).max(zero);
        std::copy(&p[0], &p[0] + nc, reservoir_state.pressure().begin());


        // Saturation updates.
        const Opm::PhaseUsage& pu = fluid_.phaseUsage();
        const DataBlock s_old = Eigen::Map<const DataBlock>(& reservoir_state.saturation()[0], nc, np);
        const double dsmax = dsMax();

        // initialize with zeros
        // if the phase is active the saturation are overwritten.
        V so;
        V sw = zero;
        V sg = zero;
        V ss = zero;

        // Appleyard chop process.
        // We chop too large updates in the saturations
        {
            V maxVal = zero;
            V dso = zero;
            if (active_[Water]){
                maxVal = dsw.abs().max(maxVal);
                dso = dso - dsw;
            }

            V dsg;
            if (active_[Gas]){
                dsg = Base::isSg_ * dxvar - Base::isRv_ * dsw;
                maxVal = dsg.abs().max(maxVal);
                dso = dso - dsg;
            }
            if (has_solvent_){
                maxVal = dss.abs().max(maxVal);
                // solvent is not added note that the so calculated
                // here is overwritten later
                //dso = dso - dss;
            }

            maxVal = dso.abs().max(maxVal);

            V step = dsmax/maxVal;
            step = step.min(1.);

            if (active_[Water]) {
                const int pos = pu.phase_pos[ Water ];
                const V sw_old = s_old.col(pos);
                sw = sw_old - step * dsw;
            }

            if (active_[Gas]) {
                const int pos = pu.phase_pos[ Gas ];
                const V sg_old = s_old.col(pos);
                sg = sg_old - step * dsg;
            }
            if (has_solvent_) {
                auto& solvent_saturation = reservoir_state.getCellData( reservoir_state.SSOL );
                const V ss_old = Eigen::Map<const V>(&solvent_saturation[0], nc, 1);
                ss  = ss_old - step * dss;
            }

            const int pos = pu.phase_pos[ Oil ];
            const V so_old = s_old.col(pos);
            so = so_old - step * dso;
        }

        // solvent is not included in the adjustment for negative saturation
        auto ixg = sg < 0;
        for (int c = 0; c < nc; ++c) {
            if (ixg[c]) {
                sw[c] = sw[c] / (1-sg[c]);
                so[c] = so[c] / (1-sg[c]);
                sg[c] = 0;
            }
        }

        auto ixo = so < 0;
        for (int c = 0; c < nc; ++c) {
            if (ixo[c]) {
                sw[c] = sw[c] / (1-so[c]);
                sg[c] = sg[c] / (1-so[c]);
                so[c] = 0;
            }
        }

        auto ixw = sw < 0;
        for (int c = 0; c < nc; ++c) {
            if (ixw[c]) {
                so[c] = so[c] / (1-sw[c]);
                sg[c] = sg[c] / (1-sw[c]);
                sw[c] = 0;
            }
        }

        auto ixs = ss < 0;
        for (int c = 0; c < nc; ++c) {
            if (ixs[c]) {
                ss[c] = 0;
            }
        }


        // The oil saturation is defined to
        // fill the rest of the pore space.
        // For convergence reasons oil saturations
        // is included in the appelyard chopping
        so = V::Constant(nc,1.0) - sw - sg - ss;

        // Update rs and rv
        const double drmaxrel = drMaxRel();
        V rs;
        if (has_disgas_) {
            const V rs_old = Eigen::Map<const V>(&reservoir_state.gasoilratio()[0], nc);
            const V drs = Base::isRs_ * dxvar;
            const V drs_limited = sign(drs) * drs.abs().min( (rs_old.abs()*drmaxrel).max( ones*1.0));
            rs = rs_old - drs_limited;
            rs = rs.max(zero);
        }
        V rv;
        if (has_vapoil_) {
            const V rv_old = Eigen::Map<const V>(&reservoir_state.rv()[0], nc);
            const V drv = Base::isRv_ * dxvar;
            const V drv_limited = sign(drv) * drv.abs().min( (rv_old.abs()*drmaxrel).max( ones*1e-3));
            rv = rv_old - drv_limited;
            rv = rv.max(zero);
        }

        sd_.soMax = fluid_.satOilMax();

        // Sg is used as primal variable for water only cells.
        const double epsilon = std::sqrt(std::numeric_limits<double>::epsilon());
        auto watOnly = sw >  (1 - epsilon);

        // phase translation sg <-> rs
        std::vector<HydroCarbonState>& hydroCarbonState = reservoir_state.hydroCarbonState();
        std::fill(hydroCarbonState.begin(), hydroCarbonState.end(), HydroCarbonState::GasAndOil);

        if (has_disgas_) {
            const V rsSat0 = fluidRsSat(p_old, s_old.col(pu.phase_pos[Oil]), cells_);
            const V rsSat = fluidRsSat(p, so, cells_);
            sd_.rsSat = ADB::constant(rsSat);
            // The obvious case
            auto hasGas = (sg > 0 && Base::isRs_ == 0);

            // Set oil saturated if previous rs is sufficiently large
            const V rs_old = Eigen::Map<const V>(&reservoir_state.gasoilratio()[0], nc);
            auto gasVaporized =  ( (rs > rsSat * (1+epsilon) && Base::isRs_ == 1 ) && (rs_old > rsSat0 * (1-epsilon)) );
            auto useSg = watOnly || hasGas || gasVaporized;
            for (int c = 0; c < nc; ++c) {
                if (useSg[c]) {
                    rs[c] = rsSat[c];
                    if (watOnly[c]) {
                        so[c] = 0;
                        sg[c] = 0;
                        ss[c] = 0;
                        rs[c] = 0;
                    }

                } else {
                    hydroCarbonState[c] = HydroCarbonState::OilOnly;
                }
            }
            rs = rs.min(rsSat);
        }

        // phase transitions so <-> rv
        if (has_vapoil_) {

            // The gas pressure is needed for the rvSat calculations
            const V gaspress_old = computeGasPressure(p_old, s_old.col(Water), s_old.col(Oil), s_old.col(Gas));
            const V gaspress = computeGasPressure(p, sw, so, sg);
            const V rvSat0 = fluidRvSat(gaspress_old, s_old.col(pu.phase_pos[Oil]), cells_);
            const V rvSat = fluidRvSat(gaspress, so, cells_);
            sd_.rvSat = ADB::constant(rvSat);

            // The obvious case
            auto hasOil = (so > 0 && Base::isRv_ == 0);

            // Set oil saturated if previous rv is sufficiently large
            const V rv_old = Eigen::Map<const V>(&reservoir_state.rv()[0], nc);
            auto oilCondensed = ( (rv > rvSat * (1+epsilon) && Base::isRv_ == 1) && (rv_old > rvSat0 * (1-epsilon)) );
            auto useSg = watOnly || hasOil || oilCondensed;
            for (int c = 0; c < nc; ++c) {
                if (useSg[c]) {
                    rv[c] = rvSat[c];
                    if (watOnly[c]) {
                        so[c] = 0;
                        sg[c] = 0;
                        ss[c] = 0;
                        rv[c] = 0;
                    }

                } else {
                    hydroCarbonState[c] = HydroCarbonState::GasOnly;
                }
            }
            rv = rv.min(rvSat);

        }

        // Update the reservoir_state
        if (has_solvent_) {
            auto& solvent_saturation = reservoir_state.getCellData( reservoir_state.SSOL );
            std::copy(&ss[0], &ss[0] + nc, solvent_saturation.begin());
        }

        for (int c = 0; c < nc; ++c) {
            reservoir_state.saturation()[c*np + pu.phase_pos[ Water ]] = sw[c];
        }

        for (int c = 0; c < nc; ++c) {
            reservoir_state.saturation()[c*np + pu.phase_pos[ Gas ]] = sg[c];
        }

        if (active_[ Oil ]) {
            const int pos = pu.phase_pos[ Oil ];
            for (int c = 0; c < nc; ++c) {
                reservoir_state.saturation()[c*np + pos] = so[c];
            }
        }

        // Update the reservoir_state
        if (has_disgas_) {
            std::copy(&rs[0], &rs[0] + nc, reservoir_state.gasoilratio().begin());
        }

        if (has_vapoil_) {
            std::copy(&rv[0], &rv[0] + nc, reservoir_state.rv().begin());
        }

        wellModel().updateWellState(dwells, Base::dbhpMaxRel(), well_state);

        for( auto w = 0; w < wells().number_of_wells; ++w ) {
            if (wells().type[w] == INJECTOR) {
                continue;
            }

            for (int perf = wells().well_connpos[w]; perf < wells().well_connpos[w+1]; ++perf ) {
                int wc = wells().well_cells[perf];
                if ( (ss[wc] + sg[wc]) > 0) {
                    well_state.solventFraction()[perf] = ss[wc] / (ss[wc] + sg[wc]);
                } else {
                    well_state.solventFraction()[perf] = 0.0;
                }
            }
        }



        // Update phase conditions used for property calculations.
        updatePhaseCondFromPrimalVariable(reservoir_state);
    }



    template <class Grid>
    void
    BlackoilSolventModel<Grid>::computeMassFlux(const int               actph ,
                                                const V&                transi,
                                                const ADB&              kr    ,
                                                const ADB&              mu    ,
                                                const ADB&              rho    ,
                                                const ADB&              phasePressure,
                                                const SolutionState&    state)
    {

        const int canonicalPhaseIdx = canph_[ actph ];
        // make a copy to make it possible to modify it
        ADB kr_mod = kr;
        if (canonicalPhaseIdx == Gas) {
            if (has_solvent_) {
                const int  nc   = Opm::UgGridHelpers::numCells(grid_);
                const Opm::PhaseUsage& pu = fluid_.phaseUsage();
                const ADB zero = ADB::constant(V::Zero(nc));
                const V ones = V::Constant(nc, 1.0);

                const ADB& ss = state.solvent_saturation;
                const ADB& sg = (active_[ Gas ]
                                 ? state.saturation[ pu.phase_pos[ Gas ] ]
                                 : zero);
                Selector<double> zero_selector(ss.value() + sg.value(), Selector<double>::Zero);
                const ADB F_solvent = zero_selector.select(zero, ss / (ss + sg));

                // Compute solvent properties
                const std::vector<PhasePresence>& cond = phaseCondition();
                ADB mu_s = fluidViscosity(Solvent, phasePressure,state.temperature, state.rs, state.rv, cond);
                ADB rho_s = fluidDensity(Solvent,sd_.rq[solvent_pos_].b, state.rs, state.rv);

                // Compute solvent relperm and mass flux
                ADB krs = solvent_props_.solventRelPermMultiplier(F_solvent, cells_) * kr_mod;
                Base::computeMassFlux(solvent_pos_, transi, krs, mu_s, rho_s, phasePressure, state);

                // Modify gas relperm
                kr_mod = solvent_props_.gasRelPermMultiplier( (ones - F_solvent) , cells_) * kr_mod;

            }
        }
        // Compute mobility and flux
        Base::computeMassFlux(actph, transi, kr_mod, mu, rho, phasePressure, state);

    }

    template <class Grid>
    ADB
    BlackoilSolventModel<Grid>::fluidViscosity(const int               phase,
                                                            const ADB&              p    ,
                                                            const ADB&              temp ,
                                                            const ADB&              rs   ,
                                                            const ADB&              rv   ,
                                                            const std::vector<PhasePresence>& cond) const
    {
        const Opm::PhaseUsage& pu = fluid_.phaseUsage();
        if (phase == Solvent) {
            if (!is_miscible_) {
                return solvent_props_.muSolvent(p, cells_);
            } else {
                return mu_eff_[solvent_pos_];
            }

        } else {
            if (!is_miscible_) {
                return Base::fluidViscosity(phase, p, temp, rs, rv, cond);
            } else {
                return mu_eff_[pu.phase_pos[ phase ]];
            }
        }
     }

    template <class Grid>
    ADB
    BlackoilSolventModel<Grid>::fluidReciprocFVF(const int               phase,
                                                 const ADB&              p    ,
                                                 const ADB&              temp ,
                                                 const ADB&              rs   ,
                                                 const ADB&              rv   ,
                                                 const std::vector<PhasePresence>& cond) const
    {
        const Opm::PhaseUsage& pu = fluid_.phaseUsage();
        if (phase == Solvent) {
            if (!is_miscible_) {
                return solvent_props_.bSolvent(p, cells_);
            } else {
                return b_eff_[solvent_pos_];
            }

        } else {
            if (!is_miscible_) {
                return Base::fluidReciprocFVF(phase, p, temp, rs, rv, cond);
            } else {
                return b_eff_[pu.phase_pos[ phase ]];
            }
        }
    }

    template <class Grid>
    ADB
    BlackoilSolventModel<Grid>::fluidDensity(const int  phase,
                                                          const ADB& b,
                                                          const ADB& rs,
                                                          const ADB& rv) const
    {
        if (phase == Solvent && has_solvent_) {
            return solvent_props_.solventSurfaceDensity(cells_) * b;
        } else {
            return Base::fluidDensity(phase, b, rs, rv);
        }
    }

    template <class Grid>
    std::vector<ADB>
    BlackoilSolventModel<Grid>::computeRelPerm(const SolutionState& state) const
    {
        using namespace Opm::AutoDiffGrid;
        const int               nc   = numCells(grid_);

        const ADB zero = ADB::constant(V::Zero(nc));

        const Opm::PhaseUsage& pu = fluid_.phaseUsage();
        const ADB& sw = (active_[ Water ]
                         ? state.saturation[ pu.phase_pos[ Water ] ]
                         : zero);

        const ADB& so = (active_[ Oil ]
                         ? state.saturation[ pu.phase_pos[ Oil ] ]
                         : zero);

        const ADB& sg = (active_[ Gas ]
                         ? state.saturation[ pu.phase_pos[ Gas ] ]
                         : zero);

        if (has_solvent_) {
            const ADB& ss = state.solvent_saturation;
            if (is_miscible_) {

                assert(active_[ Oil ]);
                assert(active_[ Gas ]);

                std::vector<ADB> relperm = fluid_.relperm(sw, so, sg+ss, cells_);

                Selector<double> zero_selector(ss.value() + sg.value(), Selector<double>::Zero);
                ADB F_solvent = zero_selector.select(ss, ss / (ss + sg));
                const ADB& po = state.canonical_phase_pressures[ Oil ];
                const ADB misc = solvent_props_.miscibilityFunction(F_solvent, cells_)
                        * solvent_props_.pressureMiscibilityFunction(po, cells_);

                const ADB sn = ss + so + sg;

                // adjust endpoints
                const V sgcr = fluid_.scaledCriticalGasSaturations(cells_);
                const V sogcr = fluid_.scaledCriticalOilinGasSaturations(cells_);
                const ADB sorwmis = solvent_props_.miscibleResidualOilSaturationFunction(sw, cells_);
                const ADB sgcwmis = solvent_props_.miscibleCriticalGasSaturationFunction(sw, cells_);

                const V ones = V::Constant(nc, 1.0);
                ADB sor = misc * sorwmis + (ones - misc) * sogcr;
                ADB sgc = misc * sgcwmis + (ones - misc) * sgcr;

                ADB ssg = ss + sg - sgc;
                const ADB sn_eff = sn - sor - sgc;

                // avoid negative values
                Selector<double> negSsg_selector(ssg.value(), Selector<double>::LessZero);
                ssg = negSsg_selector.select(zero, ssg);

                // avoid negative value and division on zero
                Selector<double> zeroSn_selector(sn_eff.value(), Selector<double>::LessEqualZero);
                const ADB F_totalGas = zeroSn_selector.select( zero, ssg / sn_eff);

                const ADB mkrgt = solvent_props_.miscibleSolventGasRelPermMultiplier(F_totalGas, cells_) * solvent_props_.misicibleHydrocarbonWaterRelPerm(sn, cells_);
                const ADB mkro = solvent_props_.miscibleOilRelPermMultiplier(ones - F_totalGas, cells_) * solvent_props_.misicibleHydrocarbonWaterRelPerm(sn, cells_);

                const V eps = V::Constant(nc, 1e-7);
                Selector<double> noOil_selector(so.value()-eps, Selector<double>::LessEqualZero);
                relperm[Gas] = (ones - misc) * relperm[Gas] + misc * mkrgt;
                relperm[Oil] = noOil_selector.select(relperm[Oil], (ones - misc) * relperm[Oil] + misc * mkro);
                return relperm;
            } else {
                return fluid_.relperm(sw, so, sg+ss, cells_);
            }
        } else {
            return fluid_.relperm(sw, so, sg, cells_);
        }

    }


    template <class Grid>
    void
    BlackoilSolventModel<Grid>::computeEffectiveProperties(const SolutionState&    state)
    {
        // Viscosity
        const Opm::PhaseUsage& pu = fluid_.phaseUsage();
        const int np = fluid_.numPhases();
        const int nc   = Opm::UgGridHelpers::numCells(grid_);
        const ADB zero = ADB::constant(V::Zero(nc));

        const ADB& pw = state.canonical_phase_pressures[pu.phase_pos[Water]];
        const ADB& po = state.canonical_phase_pressures[pu.phase_pos[Oil]];
        const ADB& pg = state.canonical_phase_pressures[pu.phase_pos[Gas]];
        const std::vector<PhasePresence>& cond = phaseCondition();

        const ADB mu_w = fluid_.muWat(pw, state.temperature, cells_);
        const ADB mu_o = fluid_.muOil(po, state.temperature, state.rs, cond, cells_);
        const ADB mu_g = fluid_.muGas(pg, state.temperature, state.rv, cond, cells_);
        const ADB mu_s = solvent_props_.muSolvent(pg,cells_);
        std::vector<ADB> viscosity(np + 1, ADB::null());
        viscosity[pu.phase_pos[Oil]] = mu_o;
        viscosity[pu.phase_pos[Gas]] = mu_g;
        viscosity[pu.phase_pos[Water]] = mu_w;
        viscosity[solvent_pos_] = mu_s;

        // Density
        const ADB bw = fluid_.bWat(pw, state.temperature, cells_);
        const ADB bo = fluid_.bOil(po, state.temperature, state.rs, cond, cells_);
        const ADB bg = fluid_.bGas(pg, state.temperature, state.rv, cond, cells_);
        const ADB bs = solvent_props_.bSolvent(pg, cells_);

        std::vector<ADB> density(np + 1, ADB::null());
        density[pu.phase_pos[Oil]] = fluidDensity(Oil, bo, state.rs, state.rv);
        density[pu.phase_pos[Gas]] = fluidDensity(Gas, bg, state.rs, state.rv);
        density[pu.phase_pos[Water]] = fluidDensity(Water, bw, state.rs, state.rv);
        density[solvent_pos_] = fluidDensity(Solvent, bs, state.rs, state.rv);

        // Effective saturations
        const ADB& ss = state.solvent_saturation;
        const ADB& so = state.saturation[ pu.phase_pos[ Oil ] ];
        const ADB& sg = (active_[ Gas ]
                         ? state.saturation[ pu.phase_pos[ Gas ] ]
                         : zero);
        const ADB& sw = (active_[ Water ]
                         ? state.saturation[ pu.phase_pos[ Water ] ]
                         : zero);

        const ADB sorwmis = solvent_props_.miscibleResidualOilSaturationFunction(sw, cells_);
        const ADB sgcwmis = solvent_props_.miscibleCriticalGasSaturationFunction(sw, cells_);

        std::vector<ADB> effective_saturations (np + 1, ADB::null());
        effective_saturations[pu.phase_pos[Oil]] = so - sorwmis;
        effective_saturations[pu.phase_pos[Gas]] = sg - sgcwmis;
        effective_saturations[pu.phase_pos[Water]] = sw;
        effective_saturations[solvent_pos_] = ss - sgcwmis;

        // Compute effective viscosities and densities
        computeToddLongstaffMixing(viscosity, density, effective_saturations, po, pu);

        // compute the volume factors from the densities
        const ADB b_eff_o = density[pu.phase_pos[ Oil ]] / (fluid_.surfaceDensity(pu.phase_pos[ Oil ],  cells_) + fluid_.surfaceDensity(pu.phase_pos[ Gas ], cells_) * state.rs);
        const ADB b_eff_g = density[pu.phase_pos[ Gas ]] / (fluid_.surfaceDensity(pu.phase_pos[ Gas ],  cells_) + fluid_.surfaceDensity(pu.phase_pos[ Oil ], cells_) * state.rv);
        const ADB b_eff_s = density[solvent_pos_] / solvent_props_.solventSurfaceDensity(cells_);

        // account for pressure effects and store the computed volume factors and viscosities
        const V ones = V::Constant(nc, 1.0);
        const ADB pmisc = solvent_props_.pressureMiscibilityFunction(po, cells_);

        b_eff_[pu.phase_pos[ Oil ]] = pmisc * b_eff_o + (ones - pmisc) * bo;
        b_eff_[pu.phase_pos[ Gas ]] = pmisc * b_eff_g + (ones - pmisc) * bg;
        b_eff_[solvent_pos_] = pmisc * b_eff_s + (ones - pmisc) * bs;

        // keep the mu*b interpolation
        mu_eff_[pu.phase_pos[ Oil ]] = b_eff_[pu.phase_pos[ Oil ]] / (pmisc * b_eff_o / viscosity[pu.phase_pos[ Oil ]] + (ones - pmisc) * bo / mu_o);
        mu_eff_[pu.phase_pos[ Gas ]] = b_eff_[pu.phase_pos[ Gas ]] / (pmisc * b_eff_g / viscosity[pu.phase_pos[ Gas ]] + (ones - pmisc) * bg / mu_g);
        mu_eff_[solvent_pos_] = b_eff_[solvent_pos_] / (pmisc * b_eff_s / viscosity[solvent_pos_] + (ones - pmisc) * bs / mu_s);

        // for water the pure values are used
        mu_eff_[pu.phase_pos[ Water ]] = mu_w;
        b_eff_[pu.phase_pos[ Water ]] = bw;
    }

    template <class Grid>
    void
    BlackoilSolventModel<Grid>::computeToddLongstaffMixing(std::vector<ADB>& viscosity, std::vector<ADB>& density, const std::vector<ADB>& saturations, const ADB po, const Opm::PhaseUsage pu)
    {
        const int  nc = cells_.size();
        const V ones = V::Constant(nc, 1.0);
        const ADB zero = ADB::constant(V::Zero(nc));

        // Calculation of effective saturations
        ADB so_eff = saturations[pu.phase_pos[ Oil ]];
        ADB sg_eff = saturations[pu.phase_pos[ Gas ]];
        ADB ss_eff = saturations[solvent_pos_];

        // Avoid negative values
        Selector<double> negative_selectorSo(so_eff.value(), Selector<double>::LessZero);
        Selector<double> negative_selectorSg(sg_eff.value(), Selector<double>::LessZero);
        Selector<double> negative_selectorSs(ss_eff.value(), Selector<double>::LessZero);
        so_eff = negative_selectorSo.select(zero, so_eff);
        sg_eff = negative_selectorSg.select(zero, sg_eff);
        ss_eff = negative_selectorSs.select(zero, ss_eff);

        // Make copy of the pure viscosities
        const ADB mu_o = viscosity[pu.phase_pos[ Oil ]];
        const ADB mu_g = viscosity[pu.phase_pos[ Gas ]];
        const ADB mu_s = viscosity[solvent_pos_];

        const ADB sn_eff =  so_eff + sg_eff + ss_eff;
        const ADB sos_eff = so_eff + ss_eff;
        const ADB ssg_eff = ss_eff + sg_eff;

        // Avoid division by zero
        Selector<double> zero_selectorSos(sos_eff.value(), Selector<double>::Zero);
        Selector<double> zero_selectorSsg(ssg_eff.value(), Selector<double>::Zero);
        Selector<double> zero_selectorSn(sn_eff.value(), Selector<double>::Zero);

        const ADB mu_s_pow = pow(mu_s, 0.25);
        const ADB mu_o_pow = pow(mu_o, 0.25);
        const ADB mu_g_pow = pow(mu_g, 0.25);

        const ADB mu_mos = zero_selectorSos.select(mu_o, mu_o * mu_s / pow( ( (so_eff / sos_eff) * mu_s_pow) + ( (ss_eff / sos_eff) * mu_o_pow) , 4.0));
        const ADB mu_msg = zero_selectorSsg.select(mu_g, mu_g * mu_s / pow( ( (sg_eff / ssg_eff) * mu_s_pow) + ( (ss_eff / ssg_eff) * mu_g_pow) , 4.0));
        const ADB mu_m = zero_selectorSn.select(mu_s, mu_o * mu_s * mu_g / pow( ( (so_eff / sn_eff) * mu_s_pow *  mu_g_pow)
                                                       + ( (ss_eff / sn_eff) * mu_o_pow *  mu_g_pow) + ( (sg_eff / sn_eff) * mu_s_pow * mu_o_pow), 4.0));
        // Mixing parameter for viscosity
        // The pressureMixingParameter represent the miscibility of the solvent while the mixingParameterViscosity the effect of the porous media.
        // The pressureMixingParameter is not implemented in ecl100.
        const ADB mix_param_mu = solvent_props_.mixingParameterViscosity(cells_) * solvent_props_.pressureMixingParameter(po, cells_);

        Selector<double> zero_selectorSs(ss_eff.value(), Selector<double>::Zero);
        // Update viscosities, use pure values if solvent saturation is zero
        viscosity[pu.phase_pos[ Oil ]] = zero_selectorSs.select(mu_o, pow(mu_o,ones - mix_param_mu) * pow(mu_mos, mix_param_mu));
        viscosity[pu.phase_pos[ Gas ]] = zero_selectorSs.select(mu_g, pow(mu_g,ones - mix_param_mu) * pow(mu_msg, mix_param_mu));
        viscosity[solvent_pos_] =  zero_selectorSs.select(mu_s, pow(mu_s,ones - mix_param_mu) * pow(mu_m, mix_param_mu));

        // Density
        ADB& rho_o = density[pu.phase_pos[ Oil ]];
        ADB& rho_g = density[pu.phase_pos[ Gas ]];
        ADB& rho_s = density[solvent_pos_];

        // mixing parameter for density
        const ADB mix_param_rho = solvent_props_.mixingParameterDensity(cells_) * solvent_props_.pressureMixingParameter(po, cells_);

        // compute effective viscosities for density calculations. These have to
        // be recomputed as a different mixing parameter may be used.
        const ADB mu_o_eff = pow(mu_o,ones - mix_param_rho) * pow(mu_mos, mix_param_rho);
        const ADB mu_g_eff = pow(mu_g,ones - mix_param_rho) * pow(mu_msg, mix_param_rho);
        const ADB mu_s_eff = pow(mu_s,ones - mix_param_rho) * pow(mu_m, mix_param_rho);

        const ADB sog_eff = so_eff + sg_eff;
        // Avoid division by zero
        Selector<double> zero_selectorSog_eff(sog_eff.value(), Selector<double>::Zero);
        const ADB sof = zero_selectorSog_eff.select(zero , so_eff / sog_eff);
        const ADB sgf = zero_selectorSog_eff.select(zero , sg_eff / sog_eff);

        // Effective densities
        const ADB mu_sog_pow = mu_s_pow * ( (sgf * mu_o_pow) + (sof * mu_g_pow) );
        const ADB mu_o_eff_pow = pow(mu_o_eff, 0.25);
        const ADB mu_g_eff_pow = pow(mu_g_eff, 0.25);
        const ADB mu_s_eff_pow = pow(mu_s_eff, 0.25);       
        const ADB sfraction_oe = (mu_o_pow * (mu_o_eff_pow - mu_s_pow)) / (mu_o_eff_pow * (mu_o_pow - mu_s_pow));
        const ADB sfraction_ge = (mu_g_pow * (mu_s_pow - mu_g_eff_pow)) / (mu_g_eff_pow * (mu_s_pow - mu_g_pow));
        const ADB sfraction_se = (mu_sog_pow - ( mu_o_pow * mu_g_pow * mu_s_pow / mu_s_eff_pow) ) / ( mu_sog_pow - (mu_o_pow * mu_g_pow));
        const ADB rho_o_eff = (rho_o * sfraction_oe) + (rho_s * (ones - sfraction_oe));
        const ADB rho_g_eff = (rho_g * sfraction_ge) + (rho_s * (ones - sfraction_ge));
        const ADB rho_s_eff = (rho_s * sfraction_se) + (rho_g * sgf * (ones - sfraction_se)) + (rho_o * sof * (ones - sfraction_se));

        // Avoid division by zero for equal mobilities. For equal mobilities the effecitive density is calculated
        // based on the saturation fraction directly.
        Selector<double> unitGasSolventMobilityRatio_selector(mu_s.value() - mu_g.value(), Selector<double>::Zero);
        Selector<double> unitOilSolventMobilityRatio_selector(mu_s.value() - mu_o.value(), Selector<double>::Zero);

        // Effective densities when the mobilities are equal
        const ADB rho_m = zero_selectorSn.select(zero, (rho_o * so_eff / sn_eff) + (rho_g * sg_eff / sn_eff) + (rho_s * ss_eff / sn_eff));
        const ADB rho_o_eff_simple = ((ones - mix_param_rho) * rho_o) + (mix_param_rho * rho_m);
        const ADB rho_g_eff_simple = ((ones - mix_param_rho) * rho_g) + (mix_param_rho * rho_m);
        //const ADB rho_s_eff_simple = ((ones - mix_param_rho) * rho_s) + (mix_param_rho * rho_m);

        // Update densities, use pure values if solvent saturation is zero
        rho_o = zero_selectorSs.select(rho_o, unitOilSolventMobilityRatio_selector.select(rho_o_eff_simple, rho_o_eff) );
        rho_g = zero_selectorSs.select(rho_g, unitGasSolventMobilityRatio_selector.select(rho_g_eff_simple, rho_g_eff) );
        rho_s = zero_selectorSs.select(rho_s, rho_s_eff);


    }


    template <class Grid>
    std::vector<ADB>
    BlackoilSolventModel<Grid>::
    computePressures(const ADB& po,
                     const ADB& sw,
                     const ADB& so,
                     const ADB& sg,
                     const ADB& ss) const
    {
        std::vector<ADB> pressures = Base::computePressures(po, sw, so, sg);

        if (has_solvent_) {
            // The imiscible capillary pressure is evaluated using the total gas saturation (sg + ss)
            std::vector<ADB> pressures_imisc = Base::computePressures(po, sw, so, sg + ss);

            // Pressure effects on capillary pressure miscibility
            const ADB pmisc = solvent_props_.pressureMiscibilityFunction(po, cells_);
            // Only the pcog is effected by the miscibility. Since pg = po + pcog, changing pg is eqvivalent
            // to changing the gas pressure directly.
            const int  nc = cells_.size();
            const V ones = V::Constant(nc, 1.0);
            pressures[Gas] = ( pmisc * pressures[Gas] + ((ones - pmisc) * pressures_imisc[Gas]));
        }
        return pressures;
    }




    template <class Grid>
    std::vector<std::vector<double> >
    BlackoilSolventModel<Grid>::
    computeFluidInPlace(const ReservoirState& x,
                        const std::vector<int>& fipnum)
    {
        if (has_solvent_ && is_miscible_ && b_eff_[0].size() == 0) {
            // A hack to avoid trouble for initial fluid in place, due to usage of b_eff_.
            WellState xw, xwdummy;
            const Opm::PhaseUsage& pu = fluid_.phaseUsage();
            xw.init(&wells(), x, xwdummy, pu);
            SolutionState solstate = variableState(x, xw);
            computeEffectiveProperties(solstate);
        }
        return Base::computeFluidInPlace(x, fipnum);
    }



}


#endif // OPM_BLACKOILSOLVENT_IMPL_HEADER_INCLUDED
