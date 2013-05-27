/*
  Copyright 2013 SINTEF ICT, Applied Mathematics.

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

#include <opm/autodiff/FullyImplicitBlackoilSolver.hpp>


#include <opm/autodiff/AutoDiffBlock.hpp>
#include <opm/autodiff/AutoDiffHelpers.hpp>
#include <opm/autodiff/BlackoilPropsAdInterface.hpp>
#include <opm/autodiff/GeoProps.hpp>

#include <opm/core/grid.h>
#include <opm/core/linalg/LinearSolverInterface.hpp>
#include <opm/core/simulator/BlackoilState.hpp>
#include <opm/core/simulator/WellState.hpp>
#include <opm/core/utility/ErrorMacros.hpp>

#include <cassert>
#include <iomanip>


typedef AutoDiff::ForwardBlock<double> ADB;
typedef ADB::V V;
typedef ADB::M M;
typedef Eigen::Array<double,
                     Eigen::Dynamic,
                     Eigen::Dynamic,
                     Eigen::RowMajor> DataBlock;


namespace {


    std::vector<int>
    buildAllCells(const int nc)
    {
        std::vector<int> all_cells(nc);

        for (int c = 0; c < nc; ++c) { all_cells[c] = c; }

        return all_cells;
    }



    template <class GeoProps>
    AutoDiff::ForwardBlock<double>::M
    gravityOperator(const UnstructuredGrid& grid,
                    const HelperOps&        ops ,
                    const GeoProps&         geo )
    {
        const int nc = grid.number_of_cells;

        std::vector<int> f2hf(2 * grid.number_of_faces, -1);
        for (int c = 0, i = 0; c < nc; ++c) {
            for (; i < grid.cell_facepos[c + 1]; ++i) {
                const int f = grid.cell_faces[ i ];
                const int p = 0 + (grid.face_cells[2*f + 0] != c);

                f2hf[2*f + p] = i;
            }
        }

        typedef AutoDiff::ForwardBlock<double>::V V;
        typedef AutoDiff::ForwardBlock<double>::M M;

        const V& gpot  = geo.gravityPotential();
        const V& trans = geo.transmissibility();

        const HelperOps::IFaces::Index ni = ops.internal_faces.size();

        typedef Eigen::Triplet<double> Tri;
        std::vector<Tri> grav;  grav.reserve(2 * ni);
        for (HelperOps::IFaces::Index i = 0; i < ni; ++i) {
            const int f  = ops.internal_faces[ i ];
            const int c1 = grid.face_cells[2*f + 0];
            const int c2 = grid.face_cells[2*f + 1];

            assert ((c1 >= 0) && (c2 >= 0));

            const double dG1 = gpot[ f2hf[2*f + 0] ];
            const double dG2 = gpot[ f2hf[2*f + 1] ];
            const double t   = trans[ f ];

            grav.push_back(Tri(i, c1,   t * dG1));
            grav.push_back(Tri(i, c2, - t * dG2));
        }

        M G(ni, nc);  G.setFromTriplets(grav.begin(), grav.end());

        return G;
    }



    template <class PU>
    std::vector<bool>
    activePhases(const PU& pu)
    {
        const int maxnp = Opm::BlackoilPhases::MaxNumPhases;
        std::vector<bool> active(maxnp, false);

        for (int p = 0; p < pu.MaxNumPhases; ++p) {
            active[ p ] = pu.phase_used[ p ] != 0;
        }

        return active;
    }



    template <class PU>
    std::vector<int>
    active2Canonical(const PU& pu)
    {
        const int maxnp = Opm::BlackoilPhases::MaxNumPhases;
        std::vector<int> act2can(maxnp, -1);

        for (int phase = 0; phase < maxnp; ++phase) {
            if (pu.phase_used[ phase ]) {
                act2can[ pu.phase_pos[ phase ] ] = phase;
            }
        }

        return act2can;
    }


} // Anonymous namespace



namespace Opm {


    FullyImplicitBlackoilSolver::
    FullyImplicitBlackoilSolver(const UnstructuredGrid&         grid ,
                                const BlackoilPropsAdInterface& fluid,
                                const DerivedGeology&           geo  ,
                                const Wells&                    wells,
                                const LinearSolverInterface&    linsolver)
        : grid_  (grid)
        , fluid_ (fluid)
        , geo_   (geo)
        , wells_ (wells)
        , linsolver_ (linsolver)
        , active_(activePhases(fluid.phaseUsage()))
        , canph_ (active2Canonical(fluid.phaseUsage()))
        , cells_ (buildAllCells(grid.number_of_cells))
        , ops_   (grid)
        , wops_  (wells)
        , grav_  (gravityOperator(grid_, ops_, geo_))
        , rq_    (fluid.numPhases())
        , residual_ ( { std::vector<ADB>(fluid.numPhases(), ADB::null()),
                        ADB::null() } )
    {
    }





    void
    FullyImplicitBlackoilSolver::
    step(const double   dt,
         BlackoilState& x ,
         WellState&     xw)
    {
        const V dtpv = geo_.poreVolume() / dt;

        {
            const SolutionState state = constantState(x, xw);
            computeAccum(state, 0);
        }

        const double atol  = 1.0e-12;
        const double rtol  = 5.0e-8;
        const int    maxit = 15;

        assemble(dtpv, x, xw);

        const double r0  = residualNorm();
        int          it  = 0;
        std::cout << "\nIteration         Residual\n"
                  << std::setw(9) << it << std::setprecision(9)
                  << std::setw(18) << r0 << std::endl;
        bool resTooLarge = r0 > atol;
        while (resTooLarge && (it < maxit)) {
            solveJacobianSystem(x, xw);

            assemble(dtpv, x, xw);

            const double r = residualNorm();

            resTooLarge = (r > atol) && (r > rtol*r0);

            it += 1;
            std::cout << std::setw(9) << it << std::setprecision(9)
                      << std::setw(18) << r << std::endl;
        }

        if (resTooLarge) {
            std::cerr << "Failed to compute converged solution in " << it << " iterations. Ignoring!\n";
            // THROW("Failed to compute converged solution in " << it << " iterations.");
        }
    }





    FullyImplicitBlackoilSolver::ReservoirResidualQuant::ReservoirResidualQuant()
        : accum(2, ADB::null())
        , mflux(   ADB::null())
        , b    (   ADB::null())
        , head (   ADB::null())
        , mob  (   ADB::null())
    {
    }





    FullyImplicitBlackoilSolver::SolutionState::SolutionState(const int np)
        : pressure  (    ADB::null())
        , saturation(np, ADB::null())
        , Rs        (    ADB::null())
        , bhp       (    ADB::null())
    {
    }





    FullyImplicitBlackoilSolver::
    WellOps::WellOps(const Wells& wells)
        : w2p(wells.well_connpos[ wells.number_of_wells ],
              wells.number_of_wells)
        , p2w(wells.number_of_wells,
              wells.well_connpos[ wells.number_of_wells ])
    {
        const int        nw   = wells.number_of_wells;
        const int* const wpos = wells.well_connpos;

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
    }





    FullyImplicitBlackoilSolver::SolutionState
    FullyImplicitBlackoilSolver::constantState(const BlackoilState& x,
                                               const WellState&     xw)
    {
        const int nc = grid_.number_of_cells;
        const int np = x.numPhases();

        // The block pattern assumes the following primary variables:
        //    pressure
        //    water saturation (if water present)
        //    gas saturation (if gas present)
        //    gas solution factor (if both gas and oil present)
        //    well bottom-hole pressure
        // Note that oil is assumed to always be present, but is never
        // a primary variable.
        ASSERT(active_[ Oil ]);
        std::vector<int> bpat(np, nc);
        const bool gasandoil = (active_[ Oil ] && active_[ Gas ]);
        if (gasandoil) {
            bpat.push_back(nc);
        }
        bpat.push_back(xw.bhp().size());

        SolutionState state(np);

        // Pressure.
        assert (not x.pressure().empty());
        const V p = Eigen::Map<const V>(& x.pressure()[0], nc, 1);
        state.pressure = ADB::constant(p, bpat);

        // Saturation.
        assert (not x.saturation().empty());
        const DataBlock s = Eigen::Map<const DataBlock>(& x.saturation()[0], nc, np);
        const Opm::PhaseUsage pu = fluid_.phaseUsage();
        {
            V so = V::Ones(nc, 1);
            if (active_[ Water ]) {
                const int pos = pu.phase_pos[ Water ];
                const V   sw  = s.col(pos);
                so -= sw;

                state.saturation[pos] = ADB::constant(sw, bpat);
            }
            if (active_[ Gas ]) {
                const int pos = pu.phase_pos[ Gas ];
                const V   sg  = s.col(pos);
                so -= sg;

                state.saturation[pos] = ADB::constant(sg, bpat);
            }
            if (active_[ Oil ]) {
                const int pos = pu.phase_pos[ Oil ];
                state.saturation[pos] = ADB::constant(so, bpat);
            }
        }

        // Gas solution factor (Rs).
        if (active_[ Oil ] && active_[ Gas ]) {
            const Eigen::Map<const DataBlock> z(&x.surfacevol()[0], nc, np);
            const Opm::PhaseUsage pu = fluid_.phaseUsage();
            const int pos_oil = pu.phase_pos[ Oil ];
            const int pos_gas = pu.phase_pos[ Gas ];
            // This may fail unless residual oil > 0.
            const V rs_from_z = z.col(pos_gas) / z.col(pos_oil);
            const V rs_max = fluidRsMax(state.pressure, cells_).value();
            const V rs = rs_from_z.min(rs_max);
            state.Rs = ADB::constant(rs, bpat);
        } else {
            const V Rs = V::Zero(nc, 1);
            state.Rs = ADB::constant(Rs, bpat);
        }

        // Well bottom-hole pressure.
        assert (not xw.bhp().empty());
        const V bhp = Eigen::Map<const V>(& xw.bhp()[0], xw.bhp().size());
        state.bhp = ADB::constant(bhp, bpat);

        return state;
    }





    FullyImplicitBlackoilSolver::SolutionState
    FullyImplicitBlackoilSolver::variableState(const BlackoilState& x,
                                               const WellState&     xw)
    {
        const int nc = grid_.number_of_cells;
        const int np = x.numPhases();

        std::vector<V> vars0;
        vars0.reserve(active_[Oil] && active_[Gas] ? np + 2 : np + 1); // Rs is primary if oil and gas present.

        // Initial pressure.
        assert (not x.pressure().empty());
        const V p = Eigen::Map<const V>(& x.pressure()[0], nc, 1);
        vars0.push_back(p);

        // Initial saturation.
        assert (not x.saturation().empty());
        const DataBlock s = Eigen::Map<const DataBlock>(& x.saturation()[0], nc, np);
        const Opm::PhaseUsage pu = fluid_.phaseUsage();
        // We do not handle a Water/Gas situation correctly, guard against it.
        ASSERT (active_[ Oil]);
        if (active_[ Water ]) {
            const V sw = s.col(pu.phase_pos[ Water ]);
            vars0.push_back(sw);
        }
        if (active_[ Gas ]) {
            const V sg = s.col(pu.phase_pos[ Gas ]);
            vars0.push_back(sg);
        }

        // Initial gas solution factor (Rs).
        if (active_[ Oil ] && active_[ Gas ]) {
            const Eigen::Map<const DataBlock> z(&x.surfacevol()[0], nc, np);
            const Opm::PhaseUsage pu = fluid_.phaseUsage();
            const int pos_oil = pu.phase_pos[ Oil ];
            const int pos_gas = pu.phase_pos[ Gas ];
            // This may fail unless residual oil > 0.
            const V rs_from_z = z.col(pos_gas) / z.col(pos_oil);
            std::vector<int> dummy_bpat;
            const V rs_max = fluidRsMax(ADB::constant(p, dummy_bpat), cells_).value();
            const V rs = rs_from_z.min(rs_max);
            vars0.push_back(rs);
        }

        // Initial well bottom-hole pressure.
        assert (not xw.bhp().empty());
        const V bhp = Eigen::Map<const V>(& xw.bhp()[0], xw.bhp().size());
        vars0.push_back(bhp);

        std::vector<ADB> vars = ADB::variables(vars0);

        SolutionState state(np);

        // Pressure.
        int nextvar = 0;
        state.pressure = vars[ nextvar++ ];

        // Saturation.
        const std::vector<int>& bpat = vars[0].blockPattern();
        {
            ADB so  = ADB::constant(V::Ones(nc, 1), bpat);

            if (active_[ Water ]) {
                ADB& sw = vars[ nextvar++ ];
                state.saturation[ pu.phase_pos[ Water ] ] = sw;

                so = so - sw;
            }
            if (active_[ Gas ]) {
                ADB& sg = vars[ nextvar++ ];
                state.saturation[ pu.phase_pos[ Gas ] ] = sg;

                so = so - sg;
            }
            if (active_[ Oil ]) {
                // Note that so is never a primary variable.
                state.saturation[ pu.phase_pos[ Oil ] ] = so;
            }
        }

        // Rs
        if (active_[ Oil ] && active_[ Gas ]) {
            state.Rs = vars[ nextvar++ ];
        } else {
            state.Rs = ADB::constant(V::Zero(nc), bpat);
        }

        // Bhp.
        state.bhp = vars[ nextvar ++];

        ASSERT(nextvar == int(vars.size()));

        return state;
    }





    void
    FullyImplicitBlackoilSolver::computeAccum(const SolutionState& state,
                                              const int            aix  )
    {
        const Opm::PhaseUsage& pu = fluid_.phaseUsage();

        const ADB&              press = state.pressure;
        const std::vector<ADB>& sat   = state.saturation;

        const int maxnp = Opm::BlackoilPhases::MaxNumPhases;
        for (int phase = 0; phase < maxnp; ++phase) {
            if (active_[ phase ]) {
                const int pos = pu.phase_pos[ phase ];
                rq_[pos].b = fluidReciprocFVF(phase, press, cells_);
                rq_[pos].accum[aix] = rq_[pos].b * sat[pos];
            }
        }

        if (active_[ Oil ] && active_[ Gas ]) {
            // Account for gas dissolved in oil.
            const int po = pu.phase_pos[ Oil ];
            const int pg = pu.phase_pos[ Gas ];

            rq_[pg].accum[aix] += state.Rs * rq_[po].accum[aix];
        }
    }





    void
    FullyImplicitBlackoilSolver::
    assemble(const V&             dtpv,
             const BlackoilState& x   ,
             const WellState&     xw  )
    {
        // Create the primary variables.
        const SolutionState state = variableState(x, xw);

        // -------- Mass balance equations --------

        // Compute b_p and the accumulation term b_p*s_p for each phase,
        // except gas. For gas, we compute b_g*s_g + Rs*b_o*s_o.
        // These quantities are stored in rq_[phase].accum[1].
        // The corresponding accumulation terms from the start of
        // the timestep (b^0_p*s^0_p etc.) were already computed
        // in step() and stored in rq_[phase].accum[0].
        computeAccum(state, 1);

        // Set up the common parts of the mass balance equations
        // for each active phase.
        const V transi = subset(geo_.transmissibility(), ops_.internal_faces);
        const std::vector<ADB> kr = computeRelPerm(state);
        for (int phase = 0; phase < fluid_.numPhases(); ++phase) {
            computeMassFlux(phase, transi, kr, state);

            residual_.mass_balance[ phase ] =
                dtpv*(rq_[phase].accum[1] - rq_[phase].accum[0])
                + ops_.div*rq_[phase].mflux;
        }

        // Add the extra (flux) terms to the gas mass balance equations
        // from gas dissolved in the oil phase.
        // The extra terms in the accumulation part of the equation
        if (active_[ Oil ] && active_[ Gas ]) {
            const int po = fluid_.phaseUsage().phase_pos[ Oil ];
            const UpwindSelector<double> upwind(grid_, ops_,
                                                rq_[po].head.value());
            const ADB Rs = upwind.select(state.Rs);

            residual_.mass_balance[ Gas ] += ops_.div * (Rs * rq_[po].mflux);
        }

        // -------- Well equation, and well contributions to the mass balance equations --------

        // Contribution to mass balance will have to wait.

        const int nc = grid_.number_of_cells;
        const int np = wells_.number_of_phases;
        const int nw = wells_.number_of_wells;
        const int nperf = wells_.well_connpos[nw];

        const std::vector<int> cells = buildAllCells(nc);
        const std::vector<int> well_cells(wells_.well_cells, wells_.well_cells + nperf);
        const V transw = Eigen::Map<const V>(wells_.WI, nperf);

        const ADB& bhp = state.bhp;

        const DataBlock well_s = wops_.w2p * Eigen::Map<const DataBlock>(wells_.comp_frac, nw, np).matrix();

        // Extract variables for perforation cell pressures
        // and corresponding perforation well pressures.
        const ADB p_perfcell = subset(state.pressure, well_cells);
        // Finally construct well perforation pressures and well flows.
        const V well_perf_dp_ = V::Zero(nperf);
        const ADB p_perfwell = wops_.w2p * bhp + well_perf_dp_;
        const ADB nkgradp_well = transw * (p_perfcell - p_perfwell);
        const Selector<double> cell_to_well_selector(nkgradp_well.value());
        ADB qs = ADB::constant(V::Zero(nw*np), state.bhp.blockPattern());
        const std::vector<ADB> well_kr = computeRelPermWells(state, well_s, well_cells);
        for (int phase = 0; phase < np; ++phase) {
            const ADB& cell_b = rq_[phase].b;
            const ADB well_b = fluidReciprocFVF(canph_[phase], p_perfwell, well_cells);
            const ADB perf_b = cell_to_well_selector.select(subset(cell_b, well_cells), well_b);

            const ADB& cell_mob = rq_[phase].mob;

            const ADB well_mu = fluidViscosity(canph_[phase], p_perfwell, well_cells);
            const ADB well_mob = well_kr[phase] / well_mu;
            const ADB perf_mob = cell_to_well_selector.select(subset(cell_mob, well_cells), well_mob);

            const ADB perf_flux = perf_mob * (nkgradp_well); // No gravity term for perforations.
            const ADB well_rates = wops_.p2w * (perf_flux*perf_b);
            qs = qs + superset(well_rates, Span(nw, 1, phase*nw), nw*np);

            const ADB well_contrib = superset(perf_flux*perf_b, well_cells, nc);
            residual_.mass_balance[phase] += well_contrib;
        }
        // Handling BHP and SURFACE_RATE wells.
        V bhp_targets(nw);
        V rate_targets(nw);
        M rate_distr(nw, np*nw);
        for (int w = 0; w < nw; ++w) {
            const WellControls* wc = wells_.ctrls[w];
            if (wc->type[wc->current] == BHP) {
                bhp_targets[w] = wc->target[wc->current];
                rate_targets[w] = -1e100;
            } else if (wc->type[wc->current] == SURFACE_RATE) {
                bhp_targets[w] = -1e100;
                rate_targets[w] = wc->target[wc->current];
                for (int phase = 0; phase < np; ++phase) {
                    rate_distr.insert(w, phase*nw + w) = wc->distr[phase];
                }
            } else {
                THROW("Can only handle BHP and SURFACE_RATE type controls.");
            }
        }
        const ADB bhp_residual = bhp - bhp_targets;
        const ADB rate_residual = rate_distr * qs - rate_targets;
        // Choose bhp residual for positive bhp targets.
        Selector<double> bhp_selector(bhp_targets);
        residual_.well_eq = bhp_selector.select(bhp_residual, rate_residual);
    }





    void FullyImplicitBlackoilSolver::solveJacobianSystem(BlackoilState& state,
                                                          WellState& well_state) const
    {
        const int np = fluid_.numPhases();
        const ADB mass_res = (np == 2) ?
            vertcat(residual_.mass_balance[0], residual_.mass_balance[1])
            : vertcat(vertcat(residual_.mass_balance[0], residual_.mass_balance[1]), 
                      residual_.mass_balance[2]);
        const ADB total_residual = collapseJacs(vertcat(mass_res, residual_.well_eq));

        const Eigen::SparseMatrix<double, Eigen::RowMajor> matr = total_residual.derivative()[0];

        V dx(V::Zero(total_residual.size()));
        Opm::LinearSolverInterface::LinearSolverReport rep
            = linsolver_.solve(matr.rows(), matr.nonZeros(),
                               matr.outerIndexPtr(), matr.innerIndexPtr(), matr.valuePtr(),
                               total_residual.value().data(), dx.data());
        if (!rep.converged) {
            THROW("ImpesTPFAAD::solve(): Linear solver convergence failure.");
        }

        const int nc = grid_.number_of_cells;
        const int nw = wells_.number_of_wells;

        // Pressure update.
        const V p_old = Eigen::Map<const V>(&state.pressure()[0], nc, 1);
        const V dp = subset(dx, Span(nc));
        const V p = p_old - dp;
        std::copy(&p[0], &p[0] + nc, state.pressure().begin());

        // Saturation updates.
        const DataBlock s_old = Eigen::Map<const DataBlock>(& state.saturation()[0], nc, np);
        int varstart = nc;
        V so = V::Constant(nc, 1.0);
        const Opm::PhaseUsage& pu = fluid_.phaseUsage();
        if (active_[ Water ]) {
            const int pos = pu.phase_pos[ Water ];
            const V sw_old = s_old.col(pos);
            const V dsw = subset(dx, Span(nc, 1, varstart));
            const V sw = sw_old - dsw;
            so -= sw;
            varstart += nc;
            for (int c = 0; c < nc; ++c) {
                state.saturation()[c*np + pos] = sw[c];
            }
        }
        if (active_[ Gas ]) {
            const int pos = pu.phase_pos[ Gas ];
            const V sg_old = s_old.col(pos);
            const V dsg = subset(dx, Span(nc, 1, varstart));
            const V sg = sg_old - dsg;
            so -= sg;
            varstart += nc;
            for (int c = 0; c < nc; ++c) {
                state.saturation()[c*np + pos] = sg[c];
            }
        }
        if (active_[ Oil ]) {
            const int pos = pu.phase_pos[ Oil ];
            for (int c = 0; c < nc; ++c) {
                state.saturation()[c*np + pos] = so[c];
            }
        }
        ASSERT(varstart + nw == total_residual.size());

        // Bhp update.
        const V bhp_old = Eigen::Map<const V>(&well_state.bhp()[0], nw, 1);
        const V dbhp = subset(dx, Span(nw, 1, varstart));
        const V bhp = bhp_old - dbhp;
        std::copy(&bhp[0], &bhp[0] + nw, well_state.bhp().begin());
    }





    std::vector<ADB>
    FullyImplicitBlackoilSolver::computeRelPerm(const SolutionState& state) const
    {
        const int               nc   = grid_.number_of_cells;
        const std::vector<int>& bpat = state.pressure.blockPattern();

        const ADB null = ADB::constant(V::Zero(nc, 1), bpat);

        const Opm::PhaseUsage& pu = fluid_.phaseUsage();
        const ADB sw = (active_[ Water ]
                        ? state.saturation[ pu.phase_pos[ Water ] ]
                        : null);

        const ADB so = (active_[ Oil ]
                        ? state.saturation[ pu.phase_pos[ Oil ] ]
                        : null);

        const ADB sg = (active_[ Gas ]
                        ? state.saturation[ pu.phase_pos[ Gas ] ]
                        : null);

        return fluid_.relperm(sw, so, sg, cells_);
    }





    std::vector<ADB>
    FullyImplicitBlackoilSolver::computeRelPermWells(const SolutionState& state,
                                                     const DataBlock& well_s,
                                                     const std::vector<int>& well_cells) const
    {
        const int nw = wells_.number_of_wells;
        const int nperf = wells_.well_connpos[nw];
        const std::vector<int>& bpat = state.pressure.blockPattern();

        const ADB null = ADB::constant(V::Zero(nperf), bpat);

        const Opm::PhaseUsage& pu = fluid_.phaseUsage();
        const ADB sw = (active_[ Water ]
                        ? ADB::constant(well_s.col(pu.phase_pos[ Water ]), bpat)
                        : null);

        const ADB so = (active_[ Oil ]
                        ? ADB::constant(well_s.col(pu.phase_pos[ Oil ]), bpat)
                        : null);

        const ADB sg = (active_[ Gas ]
                        ? ADB::constant(well_s.col(pu.phase_pos[ Gas ]), bpat)
                        : null);

        return fluid_.relperm(sw, so, sg, well_cells);
    }





    void
    FullyImplicitBlackoilSolver::computeMassFlux(const int               actph ,
                                                 const V&                transi,
                                                 const std::vector<ADB>& kr    ,
                                                 const SolutionState&    state )
    {
        const int phase = canph_[ actph ];
        const ADB mu    = fluidViscosity(phase, state.pressure, cells_);

        rq_[ actph ].mob = kr[ phase ] / mu;

        const ADB rho   = fluidDensity(phase, state.pressure, cells_);
        const ADB gflux = grav_ * rho;

        ADB& head = rq_[ actph ].head;
        head      = transi*(ops_.ngrad * state.pressure) + gflux;

        UpwindSelector<double> upwind(grid_, ops_, head.value());

        const ADB& b       = rq_[ actph ].b;
        const ADB& mob     = rq_[ actph ].mob;
        rq_[ actph ].mflux = upwind.select(b * mob) * head;
    }





    double
    FullyImplicitBlackoilSolver::residualNorm() const
    {
        double r = 0;
        for (std::vector<ADB>::const_iterator
                 b = residual_.mass_balance.begin(),
                 e = residual_.mass_balance.end();
             b != e; ++b)
        {
            r = std::max(r, (*b).value().matrix().norm());
        }
        r = std::max(r, residual_.well_eq.value().matrix().norm());

        return r;
    }





    ADB
    FullyImplicitBlackoilSolver::fluidViscosity(const int               phase,
                                                const ADB&              p    ,
                                                const std::vector<int>& cells) const
    {
        switch (phase) {
        case Water:
            return fluid_.muWat(p, cells);
        case Oil: {
            ADB dummy_rs = V::Zero(p.size(), 1) * p;
            return fluid_.muOil(p, dummy_rs, cells);
        }
        case Gas:
            return fluid_.muGas(p, cells);
        default:
            THROW("Unknown phase index " << phase);
        }
    }





    ADB
    FullyImplicitBlackoilSolver::fluidReciprocFVF(const int               phase,
                                                  const ADB&              p    ,
                                                  const std::vector<int>& cells) const
    {
        switch (phase) {
        case Water:
            return fluid_.bWat(p, cells);
        case Oil: {
            ADB dummy_rs = V::Zero(p.size(), 1) * p;
            return fluid_.bOil(p, dummy_rs, cells);
        }
        case Gas:
            return fluid_.bGas(p, cells);
        default:
            THROW("Unknown phase index " << phase);
        }
    }





    ADB
    FullyImplicitBlackoilSolver::fluidDensity(const int               phase,
                                              const ADB&              p    ,
                                              const std::vector<int>& cells) const
    {
        const double* rhos = fluid_.surfaceDensity();
        ADB b   = fluidReciprocFVF(phase, p, cells);
        ADB rho = V::Constant(p.size(), 1, rhos[phase]) * b;
        return rho;
    }





    ADB
    FullyImplicitBlackoilSolver::fluidRsMax(const ADB&              p,
                                            const std::vector<int>& cells) const
    {
        return fluid_.rsMax(p, cells);
    }


} // namespace Opm

