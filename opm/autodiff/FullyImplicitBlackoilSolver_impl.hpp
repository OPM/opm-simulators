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
#include <opm/autodiff/GridHelpers.hpp>
#include <opm/autodiff/BlackoilPropsAdInterface.hpp>
#include <opm/autodiff/GeoProps.hpp>
#include <opm/autodiff/WellDensitySegmented.hpp>
#include <opm/autodiff/WellStateFullyImplicitBlackoil.hpp>

#include <opm/core/grid.h>
#include <opm/core/linalg/LinearSolverInterface.hpp>
#include <opm/core/props/rock/RockCompressibility.hpp>
#include <opm/core/simulator/BlackoilState.hpp>
#include <opm/core/utility/ErrorMacros.hpp>
#include <opm/core/utility/Exceptions.hpp>
#include <opm/core/utility/Units.hpp>
#include <opm/core/well_controls.h>
#include <opm/core/utility/parameters/ParameterGroup.hpp>

#include <cassert>
#include <cmath>
#include <iostream>
#include <iomanip>
//#include <fstream>

// A debugging utility.
#define DUMP(foo)                                                       \
    do {                                                                \
        std::cout << "==========================================\n"     \
                  << #foo ":\n"                                         \
                  << collapseJacs(foo) << std::endl;                    \
    } while (0)

#define DUMPVAL(foo)                                                    \
    do {                                                                \
        std::cout << "==========================================\n"     \
                  << #foo ":\n"                                         \
                  << foo.value() << std::endl;                          \
    } while (0)

#define DISKVAL(foo)                                                    \
    do {                                                                \
        std::ofstream os(#foo);                                         \
        os.precision(16);                                               \
        os << foo.value() << std::endl;                                 \
    } while (0)


namespace Opm {

typedef AutoDiffBlock<double> ADB;
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

    template <class GeoProps, class Grid>
    AutoDiffBlock<double>::M
    gravityOperator(const Grid&             grid,
                    const HelperOps&        ops ,
                    const GeoProps&         geo )
    {
        using namespace Opm::AutoDiffGrid;
        const int nc = numCells(grid);
        typedef typename Opm::UgGridHelpers::Cell2FacesTraits<Grid>::Type Cell2Faces;
        Cell2Faces c2f = cell2Faces(grid);

        std::vector<int> f2hf(2 * numFaces(grid), -1);
        typename ADFaceCellTraits<Grid>::Type
            face_cells = faceCells(grid);
        for (int c = 0, i = 0; c < nc; ++c) {
            typename Cell2Faces::row_type faces=c2f[c];
            typedef typename Cell2Faces::row_type::iterator Iter;
            for (Iter f=faces.begin(), end=faces.end(); f!=end; ++f) {
                const int p = 0 + (face_cells(*f, 0) != c);

                f2hf[2*(*f) + p] = i;
            }
        }

        typedef AutoDiffBlock<double>::V V;
        typedef AutoDiffBlock<double>::M M;

        const V& gpot  = geo.gravityPotential();
        const V& trans = geo.transmissibility();

        const HelperOps::IFaces::Index ni = ops.internal_faces.size();

        typedef Eigen::Triplet<double> Tri;
        std::vector<Tri> grav;  grav.reserve(2 * ni);
        for (HelperOps::IFaces::Index i = 0; i < ni; ++i) {
            const int f  = ops.internal_faces[ i ];
            const int c1 = faceCells(grid)(f, 0);
            const int c2 = faceCells(grid)(f, 1);

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


    template<class Grid>
    V computePerfPress(const Grid& grid, const Wells& wells, const V& rho, const double grav)
    {
        using namespace Opm::AutoDiffGrid;
        const int nw = wells.number_of_wells;
        const int nperf = wells.well_connpos[nw];
        const int dim = dimensions(grid);
        V wdp = V::Zero(nperf,1);
        assert(wdp.size() == rho.size());

        // Main loop, iterate over all perforations,
        // using the following formula:
        //    wdp(perf) = g*(perf_z - well_ref_z)*rho(perf)
        // where the total density rho(perf) is taken to be
        //    sum_p (rho_p*saturation_p) in the perforation cell.
        // [although this is computed on the outside of this function].
        for (int w = 0; w < nw; ++w) {
            const double ref_depth = wells.depth_ref[w];
            for (int j = wells.well_connpos[w]; j < wells.well_connpos[w + 1]; ++j) {
                const int cell = wells.well_cells[j];
                const double cell_depth = cellCentroid(grid, cell)[dim - 1];
                wdp[j] = rho[j]*grav*(cell_depth - ref_depth);
            }
        }
        return wdp;
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



    template<class T>
    FullyImplicitBlackoilSolver<T>::
    FullyImplicitBlackoilSolver(const parameter::ParameterGroup& param,
                                const Grid&                     grid ,
                                const BlackoilPropsAdInterface& fluid,
                                const DerivedGeology&           geo  ,
                                const RockCompressibility*      rock_comp_props,
                                const Wells&                    wells,
                                const NewtonIterationBlackoilInterface&    linsolver,
                                const bool has_disgas,
                                const bool has_vapoil)
        : grid_  (grid)
        , fluid_ (fluid)
        , geo_   (geo)
        , rock_comp_props_(rock_comp_props)
        , wells_ (wells)
        , linsolver_ (linsolver)
        , active_(activePhases(fluid.phaseUsage()))
        , canph_ (active2Canonical(fluid.phaseUsage()))
        , cells_ (buildAllCells(Opm::AutoDiffGrid::numCells(grid)))
        , ops_   (grid)
        , wops_  (wells)
        , grav_  (gravityOperator(grid_, ops_, geo_))
        , has_disgas_(has_disgas)
        , has_vapoil_(has_vapoil)
        , dp_max_rel_ (1.0e9)
        , ds_max_ (0.2)
        , drs_max_rel_ (1.0e9)
        , relax_type_ (DAMPEN)
        , relax_max_ (0.5)
        , relax_increment_ (0.1)
        , relax_rel_tol_ (0.2)
        , max_iter_ (15)
        , rq_    (fluid.numPhases())
        , phaseCondition_(AutoDiffGrid::numCells(grid))
        , residual_ ( { std::vector<ADB>(fluid.numPhases(), ADB::null()),
                        ADB::null(),
                        ADB::null() } )
    {
        dp_max_rel_ = param.getDefault("dp_max_rel", dp_max_rel_);
        ds_max_ = param.getDefault("ds_max", ds_max_);
        drs_max_rel_ = param.getDefault("drs_max_rel", drs_max_rel_);
        relax_max_ = param.getDefault("relax_max", relax_max_);
        max_iter_ = param.getDefault("max_iter", max_iter_);

        std::string relaxation_type = param.getDefault("relax_type", std::string("dampen"));
        if (relaxation_type == "dampen") {
            relax_type_ = DAMPEN;
        } else if (relaxation_type == "sor") {
            relax_type_ = SOR;
        } else {
            OPM_THROW(std::runtime_error, "Unknown Relaxtion Type " << relaxation_type);
        }
    }


    template<class T>
    void
    FullyImplicitBlackoilSolver<T>::
    step(const double   dt,
         BlackoilState& x ,
         WellStateFullyImplicitBlackoil& xw)
    {
        const V pvdt = geo_.poreVolume() / dt;

        if (active_[Gas]) { updatePrimalVariableFromState(x); }

        {
            const SolutionState state = constantState(x, xw);
            computeAccum(state, 0);
            computeWellConnectionPressures(state, xw);
        }


        std::vector<std::vector<double>> residual_history;

        assemble(pvdt, x, xw);


        bool converged = false;
        double omega = 1.;
        const double r0  = residualNorm();

        residual_history.push_back(residuals());

        converged = getConvergence(dt);

        int          it  = 0;
        std::cout << "\nIteration         Residual\n"
                  << std::setw(9) << it << std::setprecision(9)
                  << std::setw(18) << r0 << std::endl;

        const int sizeNonLinear = residual_.sizeNonLinear();

        V dxOld = V::Zero(sizeNonLinear);

        bool isOscillate = false;
        bool isStagnate = false;
        const enum RelaxType relaxtype = relaxType();

        while ((!converged) && (it < maxIter())) {
            V dx = solveJacobianSystem();

            detectNewtonOscillations(residual_history, it, relaxRelTol(), isOscillate, isStagnate);

            if (isOscillate) {
                omega -= relaxIncrement();
                omega = std::max(omega, relaxMax());
                std::cout << " Oscillating behavior detected: Relaxation set to " << omega << std::endl;
            }

            stablizeNewton(dx, dxOld, omega, relaxtype);

            updateState(dx, x, xw);

            assemble(pvdt, x, xw);

            const double r = residualNorm();

            residual_history.push_back(residuals());

            converged = getConvergence(dt);

            it += 1;
            std::cout << std::setw(9) << it << std::setprecision(9)
                      << std::setw(18) << r << std::endl;
        }

        if (!converged) {
            std::cerr << "Failed to compute converged solution in " << it << " iterations. Ignoring!\n";
            // OPM_THROW(std::runtime_error, "Failed to compute converged solution in " << it << " iterations.");
        }
    }





    template<class T>
    FullyImplicitBlackoilSolver<T>::ReservoirResidualQuant::ReservoirResidualQuant()
        : accum(2, ADB::null())
        , mflux(   ADB::null())
        , b    (   ADB::null())
        , head (   ADB::null())
        , mob  (   ADB::null())
    {
    }





    template<class T>
    FullyImplicitBlackoilSolver<T>::SolutionState::SolutionState(const int np)
        : pressure  (    ADB::null())
        , saturation(np, ADB::null())
        , rs        (    ADB::null())
        , rv        (    ADB::null())
        , qs        (    ADB::null())
        , bhp       (    ADB::null())
    {
    }





    template<class T>
    FullyImplicitBlackoilSolver<T>::
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





    template<class T>
    typename FullyImplicitBlackoilSolver<T>::SolutionState
    FullyImplicitBlackoilSolver<T>::constantState(const BlackoilState& x,
                                                  const WellStateFullyImplicitBlackoil&     xw)
    {
        auto state = variableState(x, xw);

        // HACK: throw away the derivatives. this may not be the most
        // performant way to do things, but it will make the state
        // automatically consistent with variableState() (and doing
        // things automatically is all the rage in this module ;)
        state.pressure = ADB::constant(state.pressure.value());
        state.rs = ADB::constant(state.rs.value());
        state.rv = ADB::constant(state.rv.value());
        for (int phaseIdx= 0; phaseIdx < x.numPhases(); ++ phaseIdx)
            state.saturation[phaseIdx] = ADB::constant(state.saturation[phaseIdx].value());
        state.qs = ADB::constant(state.qs.value());
        state.bhp = ADB::constant(state.bhp.value());

        return state;
    }





    template<class T>
    typename FullyImplicitBlackoilSolver<T>::SolutionState
    FullyImplicitBlackoilSolver<T>::variableState(const BlackoilState& x,
                                                  const WellStateFullyImplicitBlackoil&     xw)
    {
        using namespace Opm::AutoDiffGrid;
        const int nc = numCells(grid_);
        const int np = x.numPhases();

        std::vector<V> vars0;
        // p, Sw and Rs, Rv or Sg is used as primary depending on solution conditions
        vars0.reserve(np + 1);
        // Initial pressure.
        assert (not x.pressure().empty());
        const V p = Eigen::Map<const V>(& x.pressure()[0], nc, 1);
        vars0.push_back(p);

        // Initial saturation.
        assert (not x.saturation().empty());
        const DataBlock s = Eigen::Map<const DataBlock>(& x.saturation()[0], nc, np);
        const Opm::PhaseUsage pu = fluid_.phaseUsage();
        // We do not handle a Water/Gas situation correctly, guard against it.
        assert (active_[ Oil]);
        if (active_[ Water ]) {
            const V sw = s.col(pu.phase_pos[ Water ]);
            vars0.push_back(sw);
        }

        // store cell status in vectors
        V isRs = V::Zero(nc,1);
        V isRv = V::Zero(nc,1);
        V isSg = V::Zero(nc,1);

        if (active_[ Gas ]){
            for (int c = 0; c < nc ; c++ ) {
                if ( primalVariable_[c] == PrimalVariables::RS ) {
                    isRs[c] = 1;
                }
                else if ( primalVariable_[c] == PrimalVariables::RV ) {
                    isRv[c] = 1;
                }
                else {
                    isSg[c] = 1;
                }
            }


            // define new primary variable xvar depending on solution condition
            V xvar(nc);
            const V sg = s.col(pu.phase_pos[ Gas ]);
            const V rs = Eigen::Map<const V>(& x.gasoilratio()[0], x.gasoilratio().size());
            const V rv = Eigen::Map<const V>(& x.rv()[0], x.rv().size());
            xvar = isRs*rs + isRv*rv + isSg*sg;
            vars0.push_back(xvar);
        }



        // Initial well rates.
        assert (not xw.wellRates().empty());
        // Need to reshuffle well rates, from phase running fastest
        // to wells running fastest.
        const int nw = wells_.number_of_wells;
        // The transpose() below switches the ordering.
        const DataBlock wrates = Eigen::Map<const DataBlock>(& xw.wellRates()[0], nw, np).transpose();
        const V qs = Eigen::Map<const V>(wrates.data(), nw*np);
        vars0.push_back(qs);

        // Initial well bottom-hole pressure.
        assert (not xw.bhp().empty());
        const V bhp = Eigen::Map<const V>(& xw.bhp()[0], xw.bhp().size());
        vars0.push_back(bhp);

        std::vector<ADB> vars = ADB::variables(vars0);

        SolutionState state(np);

        // Pressure.
        int nextvar = 0;
        state.pressure = vars[ nextvar++ ];

        // Saturations
        const std::vector<int>& bpat = vars[0].blockPattern();
        {
            ADB so  = ADB::constant(V::Ones(nc, 1), bpat);

            if (active_[ Water ]) {
                ADB& sw = vars[ nextvar++ ];
                state.saturation[pu.phase_pos[ Water ]] = sw;
                so = so - sw;
            }

            if (active_[ Gas]) {
                // Define Sg Rs and Rv in terms of xvar.
                ADB rsSat = fluidRsSat(state.pressure, cells_);
                ADB rvSat = fluidRvSat(state.pressure, cells_);
                ADB xvar = vars[ nextvar++ ];
                ADB sg = isSg*xvar + isRv* so;
                state.saturation[ pu.phase_pos[ Gas ] ] = sg;
                so = so - sg;

                if (has_disgas_) {
                    state.rs = (1-isRs) * rsSat + isRs*xvar;
                } else {
                    state.rs = rsSat;
                }
                if (has_vapoil_) {
                    state.rv = (1-isRv) * rvSat + isRv*xvar;
                } else {
                    state.rv = rvSat;
                }
            }

            if (active_[ Oil ]) {
                // Note that so is never a primary variable.
                state.saturation[ pu.phase_pos[ Oil ] ] = so;
            }
        }

        // Qs.
        state.qs = vars[ nextvar++ ];

        // Bhp.
        state.bhp = vars[ nextvar++ ];

        assert(nextvar == int(vars.size()));

        return state;
    }





    template<class T>
    void
    FullyImplicitBlackoilSolver<T>::computeAccum(const SolutionState& state,
                                              const int            aix  )
    {
        const Opm::PhaseUsage& pu = fluid_.phaseUsage();

        const ADB&              press = state.pressure;
        const std::vector<ADB>& sat   = state.saturation;
        const ADB&              rs    = state.rs;
        const ADB&              rv    = state.rv;

        const std::vector<PhasePresence> cond = phaseCondition();

        const ADB pv_mult = poroMult(press);

        const int maxnp = Opm::BlackoilPhases::MaxNumPhases;
        for (int phase = 0; phase < maxnp; ++phase) {
            if (active_[ phase ]) {
                const int pos = pu.phase_pos[ phase ];
                rq_[pos].b = fluidReciprocFVF(phase, press, rs, rv, cond, cells_);
                rq_[pos].accum[aix] = pv_mult * rq_[pos].b * sat[pos];
                // DUMP(rq_[pos].b);
                // DUMP(rq_[pos].accum[aix]);
            }
        }

        if (active_[ Oil ] && active_[ Gas ]) {
            // Account for gas dissolved in oil and vaporized oil
            const int po = pu.phase_pos[ Oil ];
            const int pg = pu.phase_pos[ Gas ];

            rq_[pg].accum[aix] += state.rs * rq_[po].accum[aix];
            rq_[po].accum[aix] += state.rv * rq_[pg].accum[aix];
            //DUMP(rq_[pg].accum[aix]);
        }
    }





    template<class T>
    void FullyImplicitBlackoilSolver<T>::computeWellConnectionPressures(const SolutionState& state,
                                                                        const WellStateFullyImplicitBlackoil& xw)
    {
        using namespace Opm::AutoDiffGrid;
        // 1. Compute properties required by computeConnectionPressureDelta().
        //    Note that some of the complexity of this part is due to the function
        //    taking std::vector<double> arguments, and not Eigen objects.
        const int nperf = wells_.well_connpos[wells_.number_of_wells];
        const std::vector<int> well_cells(wells_.well_cells, wells_.well_cells + nperf);
        // Compute b, rsmax, rvmax values for perforations.
        const ADB perf_press = subset(state.pressure, well_cells);
        std::vector<PhasePresence> perf_cond(nperf);
        const std::vector<PhasePresence>& pc = phaseCondition();
        for (int perf = 0; perf < nperf; ++perf) {
            perf_cond[perf] = pc[well_cells[perf]];
        }
        const PhaseUsage& pu = fluid_.phaseUsage();
        DataBlock b(nperf, pu.num_phases);
        std::vector<double> rssat_perf(nperf, 0.0);
        std::vector<double> rvsat_perf(nperf, 0.0);
        if (pu.phase_used[BlackoilPhases::Aqua]) {
            const ADB bw = fluid_.bWat(perf_press, well_cells);
            b.col(pu.phase_pos[BlackoilPhases::Aqua]) = bw.value();
        }
        if (pu.phase_used[BlackoilPhases::Liquid]) {
            const ADB perf_rs = subset(state.rs, well_cells);
            const ADB bo = fluid_.bOil(perf_press, perf_rs, perf_cond, well_cells);
            b.col(pu.phase_pos[BlackoilPhases::Liquid]) = bo.value();
            const V rssat = fluidRsSat(perf_press.value(), well_cells);
            rssat_perf.assign(rssat.data(), rssat.data() + nperf);
        }
        if (pu.phase_used[BlackoilPhases::Vapour]) {
            const ADB perf_rv = subset(state.rv, well_cells);
            const ADB bg = fluid_.bGas(perf_press, perf_rv, perf_cond, well_cells);
            b.col(pu.phase_pos[BlackoilPhases::Vapour]) = bg.value();
            const V rvsat = fluidRvSat(perf_press.value(), well_cells);
            rvsat_perf.assign(rvsat.data(), rvsat.data() + nperf);
        }
        // b is row major, so can just copy data.
        std::vector<double> b_perf(b.data(), b.data() + nperf * pu.num_phases);
        // Extract well connection depths.
        const V depth = cellCentroidsZ(grid_);
        const V pdepth = subset(depth, well_cells);
        std::vector<double> perf_depth(pdepth.data(), pdepth.data() + nperf);
        // Surface density.
        std::vector<double> surf_dens(fluid_.surfaceDensity(), fluid_.surfaceDensity() + pu.num_phases);
        // Gravity
        double grav = 0.0;
        const double* g = geo_.gravity();
        const int dim = dimensions(grid_);
        if (g) {
            // Guard against gravity in anything but last dimension.
            for (int dd = 0; dd < dim - 1; ++dd) {
                assert(g[dd] == 0.0);
            }
            grav = g[dim - 1];
        }

        // 2. Compute pressure deltas, and store the results.
        std::vector<double> cdp = WellDensitySegmented
            ::computeConnectionPressureDelta(wells_, xw, fluid_.phaseUsage(),
                                             b_perf, rssat_perf, rvsat_perf, perf_depth,
                                             surf_dens, grav);
        well_perforation_pressure_diffs_ = Eigen::Map<const V>(cdp.data(), nperf);
    }





    template<class T>
    void
    FullyImplicitBlackoilSolver<T>::
    assemble(const V&             pvdt,
             const BlackoilState& x   ,
             WellStateFullyImplicitBlackoil& xw  )
    {
        using namespace Opm::AutoDiffGrid;
        // Create the primary variables.
        SolutionState state = variableState(x, xw);

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
        const std::vector<ADB> pressures = computePressures(state);
        for (int phaseIdx = 0; phaseIdx < fluid_.numPhases(); ++phaseIdx) {
            computeMassFlux(phaseIdx, transi, kr[canph_[phaseIdx]], pressures[canph_[phaseIdx]], state);
            // std::cout << "===== kr[" << phase << "] = \n" << std::endl;
            // std::cout << kr[phase];
            // std::cout << "===== rq_[" << phase << "].mflux = \n" << std::endl;
            // std::cout << rq_[phase].mflux;

            residual_.material_balance_eq[ phaseIdx ] =
                pvdt*(rq_[phaseIdx].accum[1] - rq_[phaseIdx].accum[0])
                + ops_.div*rq_[phaseIdx].mflux;


            // DUMP(ops_.div*rq_[phase].mflux);
            // DUMP(residual_.material_balance_eq[phase]);
        }

        // -------- Extra (optional) rs and rv contributions to the mass balance equations --------

        // Add the extra (flux) terms to the mass balance equations
        // From gas dissolved in the oil phase (rs) and oil vaporized in the gas phase (rv)
        // The extra terms in the accumulation part of the equation are already handled.
        if (active_[ Oil ] && active_[ Gas ]) {
            const int po = fluid_.phaseUsage().phase_pos[ Oil ];
            const UpwindSelector<double> upwindOil(grid_, ops_,
                                                rq_[po].head.value());
            const ADB rs_face = upwindOil.select(state.rs);

            residual_.material_balance_eq[ Gas ] += ops_.div * (rs_face * rq_[po].mflux);

            const int pg = fluid_.phaseUsage().phase_pos[ Gas ];
            const UpwindSelector<double> upwindGas(grid_, ops_,
                                                rq_[pg].head.value());
            const ADB rv_face = upwindGas.select(state.rv);

            residual_.material_balance_eq[ Oil ] += ops_.div * (rv_face * rq_[pg].mflux);

            // DUMP(residual_.material_balance_eq[ Gas ]);

        }

        // Note: updateWellControls() can change all its arguments if
        // a well control is switched.
        updateWellControls(state.bhp, state.qs, xw);
        V aliveWells;
        addWellEq(state, xw, aliveWells);
        addWellControlEq(state, xw, aliveWells);
    }





    template <class T>
    void FullyImplicitBlackoilSolver<T>::addWellEq(const SolutionState& state,
                                                   WellStateFullyImplicitBlackoil& xw,
                                                   V& aliveWells)
    {
        const int nc = Opm::AutoDiffGrid::numCells(grid_);
        const int np = wells_.number_of_phases;
        const int nw = wells_.number_of_wells;
        const int nperf = wells_.well_connpos[nw];
        const Opm::PhaseUsage& pu = fluid_.phaseUsage();
        V Tw = Eigen::Map<const V>(wells_.WI, nperf);
        const std::vector<int> well_cells(wells_.well_cells, wells_.well_cells + nperf);

        // pressure diffs computed already (once per step, not changing per iteration)
        const V& cdp = well_perforation_pressure_diffs_;

        // Extract variables for perforation cell pressures
        // and corresponding perforation well pressures.
        const ADB p_perfcell = subset(state.pressure, well_cells);

        // DUMPVAL(p_perfcell);
        // DUMPVAL(state.bhp);
        // DUMPVAL(ADB::constant(cdp));

        // Pressure drawdown (also used to determine direction of flow)
        const ADB drawdown =  p_perfcell - (wops_.w2p * state.bhp + cdp);

        // current injecting connections
        auto connInjInx = drawdown.value() < 0;

        // injector == 1, producer == 0
        V isInj = V::Zero(nw);
        for (int w = 0; w < nw; ++w) {
            if (wells_.type[w] == INJECTOR) {
                isInj[w] = 1;
            }
        }

//        // A cross-flow connection is defined as a connection which has opposite
//        // flow-direction to the well total flow
//         V isInjPerf = (wops_.w2p * isInj);
//        auto crossFlowConns = (connInjInx != isInjPerf);

//        bool allowCrossFlow = true;

//        if (not allowCrossFlow) {
//            auto closedConns = crossFlowConns;
//            for (int c = 0; c < nperf; ++c) {
//                if (closedConns[c]) {
//                    Tw[c] = 0;
//                }
//            }
//            connInjInx = !closedConns;
//        }
//        TODO: not allow for crossflow


        V isInjInx = V::Zero(nperf);
        V isNotInjInx = V::Zero(nperf);
        for (int c = 0; c < nperf; ++c){
            if (connInjInx[c])
                isInjInx[c] = 1;
            else
                isNotInjInx[c] = 1;
        }


        // HANDLE FLOW INTO WELLBORE

        // compute phase volumerates standard conditions
        std::vector<ADB> cq_ps(np, ADB::null());
        for (int phase = 0; phase < np; ++phase) {
            const ADB& wellcell_mob = subset ( rq_[phase].mob, well_cells);
            const ADB cq_p = -(isNotInjInx * Tw) * (wellcell_mob * drawdown);
            cq_ps[phase] = subset(rq_[phase].b,well_cells) * cq_p;
        }
        if (active_[Oil] && active_[Gas]) {
            const int oilpos = pu.phase_pos[Oil];
            const int gaspos = pu.phase_pos[Gas];
            ADB cq_psOil = cq_ps[oilpos];
            ADB cq_psGas = cq_ps[gaspos];
            cq_ps[gaspos] += subset(state.rs,well_cells) * cq_psOil;
            cq_ps[oilpos] += subset(state.rv,well_cells) * cq_psGas;
        }

        // phase rates at std. condtions
        std::vector<ADB> q_ps(np, ADB::null());
        for (int phase = 0; phase < np; ++phase) {
            q_ps[phase] = wops_.p2w * cq_ps[phase];
        }

        // total rates at std
        ADB qt_s = ADB::constant(V::Zero(nw), state.bhp.blockPattern());
        for (int phase = 0; phase < np; ++phase) {
            qt_s += subset(state.qs, Span(nw, 1, phase*nw));
        }

        // compute avg. and total wellbore phase volumetric rates at std. conds
        const DataBlock compi = Eigen::Map<const DataBlock>(wells_.comp_frac, nw, np);
        std::vector<ADB> wbq(np, ADB::null());
        ADB wbqt = ADB::constant(V::Zero(nw), state.pressure.blockPattern());
        for (int phase = 0; phase < np; ++phase) {
            const int pos = pu.phase_pos[phase];
            wbq[phase] = (isInj * compi.col(pos)) * qt_s - q_ps[phase];
            wbqt += wbq[phase];
        }
        // DUMPVAL(wbqt);

        // check for dead wells
        aliveWells = V::Constant(nw, 1.0);
        for (int w = 0; w < nw; ++w) {
            if (wbqt.value()[w] == 0) {
                aliveWells[w] = 0.0;
            }
        }
        // compute wellbore mixture at std conds
        Selector<double> notDeadWells_selector(wbqt.value(), Selector<double>::Zero);
        std::vector<ADB> mix_s(np, ADB::null());
        for (int phase = 0; phase < np; ++phase) {
            const int pos = pu.phase_pos[phase];
            mix_s[phase] = notDeadWells_selector.select(ADB::constant(compi.col(pos)), wbq[phase]/wbqt);
        }


        // HANDLE FLOW OUT FROM WELLBORE

        // Total mobilities
        ADB mt = subset(rq_[0].mob,well_cells);
        for (int phase = 1; phase < np; ++phase) {
            mt += subset(rq_[phase].mob,well_cells);
        }

        // DUMPVAL(ADB::constant(isInjInx));
        // DUMPVAL(ADB::constant(Tw));
        // DUMPVAL(mt);
        // DUMPVAL(drawdown);

        // injection connections total volumerates
        ADB cqt_i = -(isInjInx * Tw) * (mt * drawdown);

        // compute volume ratio between connection at standard conditions
        ADB volRat = ADB::constant(V::Zero(nperf), state.pressure.blockPattern());
        std::vector<ADB> cmix_s(np, ADB::null());
        for (int phase = 0; phase < np; ++phase) {
            cmix_s[phase] = wops_.w2p * mix_s[phase];
        }

        ADB well_rv = subset(state.rv,well_cells);
        ADB well_rs = subset(state.rs,well_cells);
        ADB d = V::Constant(nperf,1.0) -  well_rv * well_rs;

        for (int phase = 0; phase < np; ++phase) {
            ADB tmp = cmix_s[phase];

            if (phase == Oil && active_[Gas]) {
                const int gaspos = pu.phase_pos[Gas];
                tmp = tmp - subset(state.rv,well_cells) * cmix_s[gaspos] / d;
            }
            if (phase == Gas && active_[Oil]) {
                const int oilpos = pu.phase_pos[Oil];
                tmp = tmp - subset(state.rs,well_cells) * cmix_s[oilpos] / d;
            }
            volRat += tmp / subset(rq_[phase].b,well_cells);
        }

        // DUMPVAL(cqt_i);
        // DUMPVAL(volRat);

        // injecting connections total volumerates at std cond
        ADB cqt_is = cqt_i/volRat;

        // connection phase volumerates at std cond
        std::vector<ADB> cq_s(np, ADB::null());
        for (int phase = 0; phase < np; ++phase) {
            cq_s[phase] = cq_ps[phase] + (wops_.w2p * mix_s[phase])*cqt_is;
        }

        // DUMPVAL(mix_s[2]);
        // DUMPVAL(cq_ps[2]);

        // Add well contributions to mass balance equations
        for (int phase = 0; phase < np; ++phase) {
            residual_.material_balance_eq[phase] -= superset(cq_s[phase],well_cells,nc);
        }


        // Add WELL EQUATIONS
        ADB qs = state.qs;
        for (int phase = 0; phase < np; ++phase) {
            qs -= superset(wops_.p2w * cq_s[phase], Span(nw, 1, phase*nw), nw*np);

        }


        V cq = superset(cq_s[0].value(), Span(nperf, np, 0), nperf*np);
        for (int phase = 1; phase < np; ++phase) {
            cq += superset(cq_s[phase].value(), Span(nperf, np, phase), nperf*np);
        }

        std::vector<double> cq_d(cq.data(), cq.data() + nperf*np);
        xw.perfPhaseRates() = cq_d;

        residual_.well_flux_eq = qs;
    }





    namespace
    {
        double rateToCompare(const ADB& well_phase_flow_rate,
                             const int well,
                             const int num_phases,
                             const double* distr)
        {
            const int num_wells = well_phase_flow_rate.size() / num_phases;
            double rate = 0.0;
            for (int phase = 0; phase < num_phases; ++phase) {
                // Important: well_phase_flow_rate is ordered with all rates for first
                // phase coming first, then all for second phase etc.
                rate += well_phase_flow_rate.value()[well + phase*num_wells] * distr[phase];
            }
            return rate;
        }

        bool constraintBroken(const ADB& bhp,
                              const ADB& well_phase_flow_rate,
                              const int well,
                              const int num_phases,
                              const WellType& well_type,
                              const WellControls* wc,
                              const int ctrl_index)
        {
            const WellControlType ctrl_type = well_controls_iget_type(wc, ctrl_index);
            const double target = well_controls_iget_target(wc, ctrl_index);
            const double* distr = well_controls_iget_distr(wc, ctrl_index);
            switch (well_type) {
            case INJECTOR:
                switch (ctrl_type) {
                case BHP:
                    return bhp.value()[well] > target;
                case SURFACE_RATE:
                    return rateToCompare(well_phase_flow_rate, well, num_phases, distr) > target;
                case RESERVOIR_RATE:
                    // Intentional fall-through.
                default:
                    OPM_THROW(std::logic_error, "Can only handle BHP and SURFACE_RATE controls.");
                }
                break;
            case PRODUCER:
                switch (ctrl_type) {
                case BHP:
                    return bhp.value()[well] < target;
                case SURFACE_RATE:
                    // Note that the rates compared below are negative,
                    // so breaking the constraints means: too high flow rate
                    // (as for injection).
                    return rateToCompare(well_phase_flow_rate, well, num_phases, distr) < target;
                case RESERVOIR_RATE:
                    // Intentional fall-through.
                default:
                    OPM_THROW(std::logic_error, "Can only handle BHP and SURFACE_RATE controls.");
                }
                break;
            default:
                OPM_THROW(std::logic_error, "Can only handle INJECTOR and PRODUCER wells.");
            }
        }
    } // anonymous namespace




    template<class T>
    void FullyImplicitBlackoilSolver<T>::updateWellControls(ADB& bhp,
                                                            ADB& well_phase_flow_rate,
                                                            WellStateFullyImplicitBlackoil& xw) const
    {
        std::string modestring[3] = { "BHP", "RESERVOIR_RATE", "SURFACE_RATE" };
        // Find, for each well, if any constraints are broken. If so,
        // switch control to first broken constraint.
        const int np = wells_.number_of_phases;
        const int nw = wells_.number_of_wells;
        bool bhp_changed = false;
        bool rates_changed = false;
        for (int w = 0; w < nw; ++w) {
            const WellControls* wc = wells_.ctrls[w];
            // The current control in the well state overrides
            // the current control set in the Wells struct, which
            // is instead treated as a default.
            const int current = xw.currentControls()[w];
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
                if (well_controls_iget_type(wc, ctrl_index) == RESERVOIR_RATE) {
                    // We cannot handle this yet.
#ifdef OPM_VERBOSE
                    std::cout << "Warning: a RESERVOIR_RATE well control exists for well "
                              << wells_.name[w] << ", but will never be checked." << std::endl;
#endif
                    continue;
                }
                if (constraintBroken(bhp, well_phase_flow_rate, w, np, wells_.type[w], wc, ctrl_index)) {
                    // ctrl_index will be the index of the broken constraint after the loop.
                    break;
                }
            }
            if (ctrl_index != nwc) {
                // Constraint number ctrl_index was broken, switch to it.
                std::cout << "Switching control mode for well " << wells_.name[w]
                          << " from " << modestring[well_controls_iget_type(wc, current)]
                          << " to " << modestring[well_controls_iget_type(wc, ctrl_index)] << std::endl;
                xw.currentControls()[w] = ctrl_index;
                // Also updating well state and primary variables.
                // We can only be switching to BHP and SURFACE_RATE
                // controls since we do not support RESERVOIR_RATE.
                const double target = well_controls_iget_target(wc, ctrl_index);
                const double* distr = well_controls_iget_distr(wc, ctrl_index);
                switch (well_controls_iget_type(wc, ctrl_index)) {
                case BHP:
                    xw.bhp()[w] = target;
                    bhp_changed = true;
                    break;
                case SURFACE_RATE:
                    for (int phase = 0; phase < np; ++phase) {
                        if (distr[phase] > 0.0) {
                            xw.wellRates()[np*w + phase] = target * distr[phase];
                        }
                    }
                    rates_changed = true;
                    break;
                default:
                    OPM_THROW(std::logic_error, "Programmer error: should not have switched "
                                                "to anything but BHP or SURFACE_RATE.");
                }
            }
        }

        // Update primary variables, if necessary.
        if (bhp_changed) {
            ADB::V new_bhp = Eigen::Map<ADB::V>(xw.bhp().data(), nw);
            bhp = ADB::function(new_bhp, bhp.derivative());
        }
        if (rates_changed) {
            // Need to reshuffle well rates, from phase running fastest
            // to wells running fastest.
            // The transpose() below switches the ordering.
            const DataBlock wrates = Eigen::Map<const DataBlock>(xw.wellRates().data(), nw, np).transpose();
            const ADB::V new_qs = Eigen::Map<const V>(wrates.data(), nw*np);
            well_phase_flow_rate = ADB::function(new_qs, well_phase_flow_rate.derivative());
        }
    }




    template<class T>
    void FullyImplicitBlackoilSolver<T>::addWellControlEq(const SolutionState& state,
                                                          const WellStateFullyImplicitBlackoil& xw,
                                                          const V& aliveWells)
    {
        // Handling BHP and SURFACE_RATE wells.

        const int np = wells_.number_of_phases;
        const int nw = wells_.number_of_wells;

        V bhp_targets(nw);
        V rate_targets(nw);
        M rate_distr(nw, np*nw);
        for (int w = 0; w < nw; ++w) {
            const WellControls* wc = wells_.ctrls[w];
            // The current control in the well state overrides
            // the current control set in the Wells struct, which
            // is instead treated as a default.
            const int current = xw.currentControls()[w];
            if (well_controls_iget_type(wc, current) == BHP) {
                bhp_targets[w] = well_controls_iget_target(wc, current);
                rate_targets[w] = -1e100;
            } else if (well_controls_iget_type(wc, current) == SURFACE_RATE) {
                bhp_targets[w] = -1e100;
                rate_targets[w] = well_controls_iget_target(wc, current);
                {
                    const double * distr = well_controls_iget_distr(wc, current);
                    for (int phase = 0; phase < np; ++phase) {
                        rate_distr.insert(w, phase*nw + w) = distr[phase];
                    }
                }
            } else {
                OPM_THROW(std::runtime_error, "Can only handle BHP and SURFACE_RATE type controls.");
            }
        }
        const ADB bhp_residual = state.bhp - bhp_targets;
        const ADB rate_residual = rate_distr * state.qs - rate_targets;
        // Choose bhp residual for positive bhp targets.
        Selector<double> bhp_selector(bhp_targets);
        residual_.well_eq = bhp_selector.select(bhp_residual, rate_residual);
        // For wells that are dead (not flowing), and therefore not communicating
        // with the reservoir, we set the equation to be equal to the well's total
        // flow. This will be a solution only if the target rate is also zero.
        M rate_summer(nw, np*nw);
        for (int w = 0; w < nw; ++w) {
            for (int phase = 0; phase < np; ++phase) {
                rate_summer.insert(w, phase*nw + w) = 1.0;
            }
        }
        Selector<double> alive_selector(aliveWells, Selector<double>::NotEqualZero);
        residual_.well_eq = alive_selector.select(residual_.well_eq, rate_summer * state.qs);
        // DUMP(residual_.well_eq);
    }





    template<class T>
    V FullyImplicitBlackoilSolver<T>::solveJacobianSystem() const
    {
        return linsolver_.computeNewtonIncrement(residual_);
    }





    namespace {
        struct Chop01 {
            double operator()(double x) const { return std::max(std::min(x, 1.0), 0.0); }
        };
    }





    template<class T>
    void FullyImplicitBlackoilSolver<T>::updateState(const V& dx,
                                                  BlackoilState& state,
                                                  WellStateFullyImplicitBlackoil& well_state)
    {
        using namespace Opm::AutoDiffGrid;
        const int np = fluid_.numPhases();
        const int nc = numCells(grid_);
        const int nw = wells_.number_of_wells;
        const V null;
        assert(null.size() == 0);
        const V zero = V::Zero(nc);
        const V one = V::Constant(nc, 1.0);

        // store cell status in vectors
        V isRs = V::Zero(nc,1);
        V isRv = V::Zero(nc,1);
        V isSg = V::Zero(nc,1);
        if (active_[Gas]) {
            for (int c = 0; c < nc ; c++ ) {
                if ( primalVariable_[c] == PrimalVariables::RS ) {
                    isRs[c] = 1;
                }
                else if ( primalVariable_[c] == PrimalVariables::RV ) {
                    isRv[c] = 1;
                }
                else {
                    isSg[c] = 1;
                }
            }
        }

        // Extract parts of dx corresponding to each part.
        const V dp = subset(dx, Span(nc));
        int varstart = nc;
        const V dsw = active_[Water] ? subset(dx, Span(nc, 1, varstart)) : null;
        varstart += dsw.size();

        const V dxvar = active_[Gas] ? subset(dx, Span(nc, 1, varstart)): null;
        varstart += dxvar.size();

        const V dqs = subset(dx, Span(np*nw, 1, varstart));
        varstart += dqs.size();
        const V dbhp = subset(dx, Span(nw, 1, varstart));
        varstart += dbhp.size();
        assert(varstart == dx.size());

        // Pressure update.
        const double dpmaxrel = dpMaxRel();
        const V p_old = Eigen::Map<const V>(&state.pressure()[0], nc, 1);
        const V absdpmax = dpmaxrel*p_old.abs();
        const V dp_limited = sign(dp) * dp.abs().min(absdpmax);
        const V p = (p_old - dp_limited).max(zero);
        std::copy(&p[0], &p[0] + nc, state.pressure().begin());


        // Saturation updates.
        const Opm::PhaseUsage& pu = fluid_.phaseUsage();
        const DataBlock s_old = Eigen::Map<const DataBlock>(& state.saturation()[0], nc, np);
        const double dsmax = dsMax();
        V so = one;
        V sw;
        V sg;

        {
            V maxVal = zero;
            V dso = zero;
            if (active_[Water]){
                maxVal = dsw.abs().max(maxVal);
                dso = dso - dsw;
            }

            V dsg;
            if (active_[Gas]){
                dsg = isSg * dxvar - isRv * dsw;
                maxVal = dsg.abs().max(maxVal);
                dso = dso - dsg;
            }

            maxVal = dso.abs().max(maxVal);

            V step = dsmax/maxVal;
            step = step.min(1.);

            if (active_[Water]) {
                const int pos = pu.phase_pos[ Water ];
                const V sw_old = s_old.col(pos);
                sw = sw_old - step * dsw;
                so -= sw;
            }

            if (active_[Gas]) {
                const int pos = pu.phase_pos[ Gas ];
                const V sg_old = s_old.col(pos);
                sg = sg_old - step * dsg;
                so -= sg;
            }
        }

        const double drsmaxrel = drsMaxRel();
        const double drvmax = 1e9;//% same as in Mrst
        V rs;
        if (has_disgas_) {
            const V rs_old = Eigen::Map<const V>(&state.gasoilratio()[0], nc);
            const V drs = isRs * dxvar;
            const V drs_limited = sign(drs) * drs.abs().min(rs_old.abs()*drsmaxrel);
            rs = rs_old - drs_limited;
        }
        V rv;
        if (has_vapoil_) {
            const V rv_old = Eigen::Map<const V>(&state.rv()[0], nc);
            const V drv = isRv * dxvar;
            const V drv_limited = sign(drv) * drv.abs().min(drvmax);
            rv = rv_old - drv_limited;
        }

        // Appleyard chop process.
        const double epsilon = std::sqrt(std::numeric_limits<double>::epsilon());
        auto watOnly = sw >  (1 - epsilon);


        // phase translation sg <-> rs
        const V rsSat0 = fluidRsSat(p_old, cells_);
        const V rsSat = fluidRsSat(p, cells_);

        // reset the phase conditions
        std::vector<PhasePresence> cond(nc);

        std::fill(primalVariable_.begin(), primalVariable_.end(), PrimalVariables::Sg);

        if (has_disgas_) {
            // The obvious case
            auto hasGas = (sg > 0 && isRs == 0);

            // keep oil saturated if previous sg is sufficient large:
            const int pos = pu.phase_pos[ Gas ];
            auto hadGas = (sg < 0 && s_old.col(pos) > epsilon);
            // Set oil saturated if previous rs is sufficiently large
            const V rs_old = Eigen::Map<const V>(&state.gasoilratio()[0], nc);
            auto gasVaporized =  ( (rs > rsSat * (1+epsilon) && isRs == 1 ) && (rs_old > rsSat0 * (1-epsilon)) );

            auto useSg = watOnly || hasGas || hadGas  || gasVaporized;
            for (int c = 0; c < nc; ++c) {
                if (useSg[c]) { rs[c] = rsSat[c];}
                else { primalVariable_[c] = PrimalVariables::RS; }

            }
        }

        // phase transitions so <-> rv
        const V rvSat0 = fluidRvSat(p_old, cells_);
        const V rvSat = fluidRvSat(p, cells_);

        if (has_vapoil_) {
            // The obvious case
            auto hasOil = (so > 0 && isRv == 0);

            // keep oil saturated if previous so is sufficient large:
            const int pos = pu.phase_pos[ Oil ];
            auto hadOil = (so < 0 && s_old.col(pos) > epsilon );
            // Set oil saturated if previous rv is sufficiently large
            const V rv_old = Eigen::Map<const V>(&state.rv()[0], nc);
            auto oilCondensed = ( (rv > rvSat * (1+epsilon) && isRv == 1) && (rv_old > rvSat0 * (1-epsilon)) );
            auto useSg = watOnly || hasOil || hadOil || oilCondensed;
            for (int c = 0; c < nc; ++c) {
                if (useSg[c]) { rv[c] = rvSat[c]; }
                else {primalVariable_[c] = PrimalVariables::RV; }

            }

        }

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
                sg[c] = sg[c] / (1-so[c]);
                sw[c] = 0;
            }
        }


        // Update saturations

        for (int c = 0; c < nc; ++c) {
            state.saturation()[c*np + pu.phase_pos[ Water ]] = sw[c];
        }

        for (int c = 0; c < nc; ++c) {
            state.saturation()[c*np + pu.phase_pos[ Gas ]] = sg[c];
        }

        if (active_[ Oil ]) {
            const int pos = pu.phase_pos[ Oil ];
            for (int c = 0; c < nc; ++c) {
                state.saturation()[c*np + pos] = so[c];
            }
        }

        // Rs and Rv updates
        if (has_disgas_)
            std::copy(&rs[0], &rs[0] + nc, state.gasoilratio().begin());

        if (has_vapoil_)
            std::copy(&rv[0], &rv[0] + nc, state.rv().begin());



        // Qs update.
        // Since we need to update the wellrates, that are ordered by wells,
        // from dqs which are ordered by phase, the simplest is to compute
        // dwr, which is the data from dqs but ordered by wells.
        const DataBlock wwr = Eigen::Map<const DataBlock>(dqs.data(), np, nw).transpose();
        const V dwr = Eigen::Map<const V>(wwr.data(), nw*np);
        const V wr_old = Eigen::Map<const V>(&well_state.wellRates()[0], nw*np);
        const V wr = wr_old - dwr;
        std::copy(&wr[0], &wr[0] + wr.size(), well_state.wellRates().begin());

        // Bhp update.
        const V bhp_old = Eigen::Map<const V>(&well_state.bhp()[0], nw, 1);
        const V dbhp_limited = sign(dbhp) * dbhp.abs().min(bhp_old.abs()*dpmaxrel);
        const V bhp = bhp_old - dbhp_limited;
        std::copy(&bhp[0], &bhp[0] + bhp.size(), well_state.bhp().begin());

    }





    template<class T>
    std::vector<ADB>
    FullyImplicitBlackoilSolver<T>::computeRelPerm(const SolutionState& state) const
    {
        using namespace Opm::AutoDiffGrid;
        const int               nc   = numCells(grid_);
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


    template<class T>
    std::vector<ADB>
    FullyImplicitBlackoilSolver<T>::computePressures(const SolutionState& state) const
    {
        using namespace Opm::AutoDiffGrid;
        const int               nc   = numCells(grid_);
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

        // convert the pressure offsets to the capillary pressures
        std::vector<ADB> pressure = fluid_.capPress(sw, so, sg, cells_);
        for (int phaseIdx = 0; phaseIdx < BlackoilPhases::MaxNumPhases; ++phaseIdx) {
            // The reference pressure is always the liquid phase (oil) pressure.
            if (phaseIdx == BlackoilPhases::Liquid)
                continue;
            pressure[phaseIdx] = pressure[phaseIdx] - pressure[BlackoilPhases::Liquid];
        }

        // Since pcow = po - pw, but pcog = pg - po,
        // we have
        //   pw = po - pcow
        //   pg = po + pcgo
        // This is an unfortunate inconsistency, but a convention we must handle.
        for (int phaseIdx = 0; phaseIdx < BlackoilPhases::MaxNumPhases; ++phaseIdx) {
            if (phaseIdx == BlackoilPhases::Aqua) {
                pressure[phaseIdx] = state.pressure - pressure[phaseIdx];
            } else {
                pressure[phaseIdx] += state.pressure;
            }
        }

        return pressure;
    }




    template<class T>
    std::vector<ADB>
    FullyImplicitBlackoilSolver<T>::computeRelPermWells(const SolutionState& state,
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





    template<class T>
    void
    FullyImplicitBlackoilSolver<T>::computeMassFlux(const int               actph ,
                                                 const V&                transi,
                                                 const ADB&              kr    ,
                                                 const ADB&              phasePressure,
                                                 const SolutionState&    state)
    {
        const int canonicalPhaseIdx = canph_[ actph ];

        const std::vector<PhasePresence> cond = phaseCondition();

        const ADB tr_mult = transMult(state.pressure);
        const ADB mu    = fluidViscosity(canonicalPhaseIdx, phasePressure, state.rs, state.rv,cond, cells_);

        rq_[ actph ].mob = tr_mult * kr / mu;

        const ADB rho   = fluidDensity(canonicalPhaseIdx, phasePressure, state.rs, state.rv,cond, cells_);

        ADB& head = rq_[ actph ].head;

        // compute gravity potensial using the face average as in eclipse and MRST
        const ADB rhoavg = ops_.caver * rho;

        const ADB dp = ops_.ngrad * phasePressure - geo_.gravity()[2] * (rhoavg * (ops_.ngrad * geo_.z().matrix()));

        head = transi*dp;
        //head      = transi*(ops_.ngrad * phasePressure) + gflux;

        UpwindSelector<double> upwind(grid_, ops_, head.value());

        const ADB& b       = rq_[ actph ].b;
        const ADB& mob     = rq_[ actph ].mob;
        rq_[ actph ].mflux = upwind.select(b * mob) * head;
        // DUMP(rq_[ actph ].mob);
        // DUMP(rq_[ actph ].mflux);
    }





    template<class T>
    double
    FullyImplicitBlackoilSolver<T>::residualNorm() const
    {
        double globalNorm = 0;
        std::vector<ADB>::const_iterator quantityIt = residual_.material_balance_eq.begin();
        const std::vector<ADB>::const_iterator endQuantityIt = residual_.material_balance_eq.end();
        for (; quantityIt != endQuantityIt; ++quantityIt) {
            const double quantityResid = (*quantityIt).value().matrix().norm();
            if (!std::isfinite(quantityResid)) {
                OPM_THROW(Opm::NumericalProblem,
                          "Encountered a non-finite residual");
            }
            globalNorm = std::max(globalNorm, quantityResid);
        }
        globalNorm = std::max(globalNorm, residual_.well_flux_eq.value().matrix().norm());
        globalNorm = std::max(globalNorm, residual_.well_eq.value().matrix().norm());

        return globalNorm;
    }


    template<class T>
    std::vector<double>
    FullyImplicitBlackoilSolver<T>::residuals() const
    {
        std::vector<double> residual;

        std::vector<ADB>::const_iterator massBalanceIt = residual_.material_balance_eq.begin();
        const std::vector<ADB>::const_iterator endMassBalanceIt = residual_.material_balance_eq.end();

        for (; massBalanceIt != endMassBalanceIt; ++massBalanceIt) {
            const double massBalanceResid = (*massBalanceIt).value().matrix().template lpNorm<Eigen::Infinity>();
            if (!std::isfinite(massBalanceResid)) {
                OPM_THROW(Opm::NumericalProblem,
                          "Encountered a non-finite residual");
            }
            residual.push_back(massBalanceResid);
        }

        // the following residuals are not used in the oscillation detection now
        const double wellFluxResid = residual_.well_flux_eq.value().matrix().template lpNorm<Eigen::Infinity>();
        if (!std::isfinite(wellFluxResid)) {
           OPM_THROW(Opm::NumericalProblem,
               "Encountered a non-finite residual");
        }
        residual.push_back(wellFluxResid);

        const double wellResid = residual_.well_eq.value().matrix().template lpNorm<Eigen::Infinity>();
        if (!std::isfinite(wellResid)) {
           OPM_THROW(Opm::NumericalProblem,
               "Encountered a non-finite residual");
        }
        residual.push_back(wellResid);

        return residual;
    }

    template<class T>
    void
    FullyImplicitBlackoilSolver<T>::detectNewtonOscillations(const std::vector<std::vector<double>>& residual_history,
                                                             const int it, const double relaxRelTol,
                                                             bool& oscillate, bool& stagnate) const
    {
        // The detection of oscillation in two primary variable results in the report of the detection
        // of oscillation for the solver.
        // Only the saturations are used for oscillation detection for the black oil model.
        // Stagnate is not used for any treatment here.

        if ( it < 2 ) {
            oscillate = false;
            stagnate = false;
            return;
        }

        stagnate = true;
        int oscillatePhase = 0;

        for (int phaseIdx= 0; phaseIdx < fluid_.numPhases(); ++ phaseIdx){
            if (active_[phaseIdx]) {
                double relChange1 = std::fabs((residual_history[it][phaseIdx] - residual_history[it - 2][phaseIdx]) /
                                               residual_history[it][phaseIdx]);
                double relChange2 = std::fabs((residual_history[it][phaseIdx] - residual_history[it - 1][phaseIdx]) /
                                               residual_history[it][phaseIdx]);
                oscillatePhase += (relChange1 < relaxRelTol) && (relChange2 > relaxRelTol);

                double relChange3 = std::fabs((residual_history[it - 1][phaseIdx] - residual_history[it - 2][phaseIdx]) /
                                               residual_history[it - 2][phaseIdx]);
                if (relChange3 > 1.e-3) {
                    stagnate = false;
                }
            }
        }

        oscillate = (oscillatePhase > 1);
    }


    template<class T>
    void
    FullyImplicitBlackoilSolver<T>::stablizeNewton(V& dx, V& dxOld, const double omega,
                                                    const RelaxType relax_type) const
    {
        // The dxOld is updated with dx.
        // If omega is equal to 1., no relaxtion will be appiled.

        const V tempDxOld = dxOld;
        dxOld = dx;

        switch (relax_type) {
            case DAMPEN:
                if (omega == 1.) {
                    return;
                }
                dx = dx*omega;
                return;
            case SOR:
                if (omega == 1.) {
                    return;
                }
                dx = dx*omega + (1.-omega)*tempDxOld;
                return;
            default:
                OPM_THROW(std::runtime_error, "Can only handle DAMPEN and SOR relaxation type.");
        }

        return;
    }

    template<class T>
    bool
    FullyImplicitBlackoilSolver<T>::getConvergence(const double dt)
    {
        const double tol_mb = 1.0e-7;
        const double tol_cnv = 1.0e-3;

        const int nc = Opm::AutoDiffGrid::numCells(grid_);
        const Opm::PhaseUsage& pu = fluid_.phaseUsage();

        const V pv = geo_.poreVolume();
        const double pvSum = pv.sum();

        const std::vector<PhasePresence> cond = phaseCondition();

        double CNVW = 0.;
        double CNVO = 0.;
        double CNVG = 0.;

        double RW_sum = 0.;
        double RO_sum = 0.;
        double RG_sum = 0.;

        double BW_avg = 0.;
        double BO_avg = 0.;
        double BG_avg = 0.;

        if (active_[Water]) {
            const int pos = pu.phase_pos[Water];
            const ADB& tempBW = rq_[pos].b;
            V BW = 1./tempBW.value();
            V RW = residual_.material_balance_eq[Water].value();
            BW_avg = BW.sum()/nc;
            const V tempV = RW.abs()/pv;

            CNVW = BW_avg * dt * tempV.maxCoeff();
            RW_sum = RW.sum();
        }

        if (active_[Oil]) {
            // Omit the disgas here. We should add it.
            const int pos = pu.phase_pos[Oil];
            const ADB& tempBO = rq_[pos].b;
            V BO = 1./tempBO.value();
            V RO = residual_.material_balance_eq[Oil].value();
            BO_avg = BO.sum()/nc;
            const V tempV = RO.abs()/pv;

            CNVO = BO_avg * dt * tempV.maxCoeff();
            RO_sum = RO.sum();
        }

        if (active_[Gas]) {
            // Omit the vapoil here. We should add it.
            const int pos = pu.phase_pos[Gas];
            const ADB& tempBG = rq_[pos].b;
            V BG = 1./tempBG.value();
            V RG = residual_.material_balance_eq[Gas].value();
            BG_avg = BG.sum()/nc;
            const V tempV = RG.abs()/pv;

            CNVG = BG_avg * dt * tempV.maxCoeff();
            RG_sum = RG.sum();
        }

        double tempValue = tol_mb * pvSum /dt;

        bool converged_MB = (fabs(BW_avg*RW_sum) < tempValue)
                         && (fabs(BO_avg*RO_sum) < tempValue)
                         && (fabs(BG_avg*RG_sum) < tempValue);

        bool converged_CNV = (CNVW < tol_cnv) && (CNVO < tol_cnv) && (CNVG < tol_cnv);

        double residualWellFlux = residual_.well_flux_eq.value().matrix().template lpNorm<Eigen::Infinity>();
        double residualWell = residual_.well_eq.value().matrix().template lpNorm<Eigen::Infinity>();

        bool converged_Well = (residualWellFlux < 1./Opm::unit::day) && (residualWell < Opm::unit::barsa);

        bool converged = converged_MB && converged_CNV && converged_Well;

#ifdef OPM_VERBOSE
        std::cout << " CNVW " << CNVW << " CNVO " << CNVO << " CNVG " << CNVG << std::endl;
        std::cout << " converged_MB " << converged_MB << " converged_CNV " << converged_CNV
                  << " converged_Well " << converged_Well << " converged " << converged << std::endl;
#endif
        return converged;
    }


    template<class T>
    ADB
    FullyImplicitBlackoilSolver<T>::fluidViscosity(const int               phase,
                                                const ADB&              p    ,
                                                const ADB&              rs   ,
                                                const ADB&              rv   ,
                                                const std::vector<PhasePresence>& cond,
                                                const std::vector<int>& cells) const
    {
        switch (phase) {
        case Water:
            return fluid_.muWat(p, cells);
        case Oil: {
            return fluid_.muOil(p, rs, cond, cells);
        }
        case Gas:
            return fluid_.muGas(p, rv, cond, cells);
        default:
            OPM_THROW(std::runtime_error, "Unknown phase index " << phase);
        }
    }





    template<class T>
    ADB
    FullyImplicitBlackoilSolver<T>::fluidReciprocFVF(const int               phase,
                                                  const ADB&              p    ,
                                                  const ADB&              rs   ,
                                                  const ADB&              rv   ,
                                                  const std::vector<PhasePresence>& cond,
                                                  const std::vector<int>& cells) const
    {
        switch (phase) {
        case Water:
            return fluid_.bWat(p, cells);
        case Oil: {
            return fluid_.bOil(p, rs, cond, cells);
        }
        case Gas:
            return fluid_.bGas(p, rv, cond, cells);
        default:
            OPM_THROW(std::runtime_error, "Unknown phase index " << phase);
        }
    }





    template<class T>
    ADB
    FullyImplicitBlackoilSolver<T>::fluidDensity(const int               phase,
                                              const ADB&              p    ,
                                              const ADB&              rs   ,
                                              const ADB&              rv   ,
                                              const std::vector<PhasePresence>& cond,
                                              const std::vector<int>& cells) const
    {
        const double* rhos = fluid_.surfaceDensity();
        ADB b = fluidReciprocFVF(phase, p, rs, rv, cond, cells);
        ADB rho = V::Constant(p.size(), 1, rhos[phase]) * b;
        if (phase == Oil && active_[Gas]) {
            // It is correct to index into rhos with canonical phase indices.
            rho += V::Constant(p.size(), 1, rhos[Gas]) * rs * b;
        }
        if (phase == Gas && active_[Oil]) {
            // It is correct to index into rhos with canonical phase indices.
            rho += V::Constant(p.size(), 1, rhos[Oil]) * rv * b;
        }
        return rho;
    }





    template<class T>
    V
    FullyImplicitBlackoilSolver<T>::fluidRsSat(const V&                p,
                                            const std::vector<int>& cells) const
    {
        return fluid_.rsSat(p, cells);
    }





    template<class T>
    ADB
    FullyImplicitBlackoilSolver<T>::fluidRsSat(const ADB&              p,
                                            const std::vector<int>& cells) const
    {
        return fluid_.rsSat(p, cells);
    }

    template<class T>
    V
    FullyImplicitBlackoilSolver<T>::fluidRvSat(const V&                p,
                                            const std::vector<int>& cells) const
    {
        return fluid_.rvSat(p, cells);
    }





    template<class T>
    ADB
    FullyImplicitBlackoilSolver<T>::fluidRvSat(const ADB&              p,
                                            const std::vector<int>& cells) const
    {
        return fluid_.rvSat(p, cells);
    }



    template<class T>
    ADB
    FullyImplicitBlackoilSolver<T>::poroMult(const ADB& p) const
    {
        const int n = p.size();
        if (rock_comp_props_ && rock_comp_props_->isActive()) {
            V pm(n);
            V dpm(n);
            for (int i = 0; i < n; ++i) {
                pm[i] = rock_comp_props_->poroMult(p.value()[i]);
                dpm[i] = rock_comp_props_->poroMultDeriv(p.value()[i]);
            }
            ADB::M dpm_diag = spdiag(dpm);
            const int num_blocks = p.numBlocks();
            std::vector<ADB::M> jacs(num_blocks);
            for (int block = 0; block < num_blocks; ++block) {
                jacs[block] = dpm_diag * p.derivative()[block];
            }
            return ADB::function(pm, jacs);
        } else {
            return ADB::constant(V::Constant(n, 1.0), p.blockPattern());
        }
    }





    template<class T>
    ADB
    FullyImplicitBlackoilSolver<T>::transMult(const ADB& p) const
    {
        const int n = p.size();
        if (rock_comp_props_ && rock_comp_props_->isActive()) {
            V tm(n);
            V dtm(n);
            for (int i = 0; i < n; ++i) {
                tm[i] = rock_comp_props_->transMult(p.value()[i]);
                dtm[i] = rock_comp_props_->transMultDeriv(p.value()[i]);
            }
            ADB::M dtm_diag = spdiag(dtm);
            const int num_blocks = p.numBlocks();
            std::vector<ADB::M> jacs(num_blocks);
            for (int block = 0; block < num_blocks; ++block) {
                jacs[block] = dtm_diag * p.derivative()[block];
            }
            return ADB::function(tm, jacs);
        } else {
            return ADB::constant(V::Constant(n, 1.0), p.blockPattern());
        }
    }


    /*
    template<class T>
    void
    FullyImplicitBlackoilSolver<T>::
    classifyCondition(const SolutionState&        state,
                      std::vector<PhasePresence>& cond ) const
    {
        const PhaseUsage& pu = fluid_.phaseUsage();

        if (active_[ Gas ]) {
            // Oil/Gas or Water/Oil/Gas system
            const int po = pu.phase_pos[ Oil ];
            const int pg = pu.phase_pos[ Gas ];

            const V&  so = state.saturation[ po ].value();
            const V&  sg = state.saturation[ pg ].value();

            cond.resize(sg.size());

            for (V::Index c = 0, e = sg.size(); c != e; ++c) {
                if (so[c] > 0)        { cond[c].setFreeOil  (); }
                if (sg[c] > 0)        { cond[c].setFreeGas  (); }
                if (active_[ Water ]) { cond[c].setFreeWater(); }
            }
        }
        else {
            // Water/Oil system
            assert (active_[ Water ]);

            const int po = pu.phase_pos[ Oil ];
            const V&  so = state.saturation[ po ].value();

            cond.resize(so.size());

            for (V::Index c = 0, e = so.size(); c != e; ++c) {
                cond[c].setFreeWater();

                if (so[c] > 0) { cond[c].setFreeOil(); }
            }
        }
    } */


    template<class T>
    void
    FullyImplicitBlackoilSolver<T>::classifyCondition(const BlackoilState& state)
    {
        using namespace Opm::AutoDiffGrid;
        const int nc = numCells(grid_);
        const int np = state.numPhases();

        const PhaseUsage& pu = fluid_.phaseUsage();
        const DataBlock s = Eigen::Map<const DataBlock>(& state.saturation()[0], nc, np);
        if (active_[ Gas ]) {
            // Oil/Gas or Water/Oil/Gas system
            const V so = s.col(pu.phase_pos[ Oil ]);
            const V sg = s.col(pu.phase_pos[ Gas ]);

            for (V::Index c = 0, e = sg.size(); c != e; ++c) {
                if (so[c] > 0)        { phaseCondition_[c].setFreeOil  (); }
                if (sg[c] > 0)        { phaseCondition_[c].setFreeGas  (); }
                if (active_[ Water ]) { phaseCondition_[c].setFreeWater(); }
            }
        }
        else {
            // Water/Oil system
            assert (active_[ Water ]);

            const V so = s.col(pu.phase_pos[ Oil ]);


            for (V::Index c = 0, e = so.size(); c != e; ++c) {
                phaseCondition_[c].setFreeWater();

                if (so[c] > 0) { phaseCondition_[c].setFreeOil(); }
            }
        }


    }

    template<class T>
    void
    FullyImplicitBlackoilSolver<T>::updatePrimalVariableFromState(const BlackoilState& state)
    {
        using namespace Opm::AutoDiffGrid;
        const int nc = numCells(grid_);
        const int np = state.numPhases();

        const PhaseUsage& pu = fluid_.phaseUsage();
        const DataBlock s = Eigen::Map<const DataBlock>(& state.saturation()[0], nc, np);

        assert (active_[ Gas ]);

        // reset the primary variables if RV and RS is not set Sg is used as primary variable.
        primalVariable_.resize(nc);
        std::fill(primalVariable_.begin(), primalVariable_.end(), PrimalVariables::Sg);

        if (has_disgas_) {
            // Oil/Gas or Water/Oil/Gas system
            const V sg = s.col(pu.phase_pos[ Gas ]);
            const V so = s.col(pu.phase_pos[ Oil ]);

            for (V::Index c = 0, e = sg.size(); c != e; ++c) {
                if ( sg[c] <= 0 && so[c] > 0 ) {primalVariable_[c] = PrimalVariables::RS; }
            }
        }

        if (has_vapoil_) {
            // Oil/Gas or Water/Oil/Gas system
            const V sg = s.col(pu.phase_pos[ Gas ]);
            const V so = s.col(pu.phase_pos[ Oil ]);

            for (V::Index c = 0, e = so.size(); c != e; ++c) {
                if (so[c] <= 0 && sg[c] > 0) {primalVariable_[c] = PrimalVariables::RV; }
            }
        }

    }



} // namespace Opm
