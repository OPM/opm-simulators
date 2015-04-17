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
#include <config.h>

#include <opm/autodiff/ImpesTPFAAD.hpp>
#include <opm/autodiff/GeoProps.hpp>
#include <opm/autodiff/GridHelpers.hpp>

#include <opm/core/simulator/BlackoilState.hpp>
#include <opm/core/simulator/WellState.hpp>
#include <opm/core/utility/ErrorMacros.hpp>
#include <opm/core/utility/Exceptions.hpp>
#include <opm/core/linalg/LinearSolverInterface.hpp>
#include <opm/core/wells.h>

#include <iostream>
#include <iomanip>

namespace Opm {

// Repeated from inside ImpesTPFAAD for convenience.
typedef AutoDiffBlock<double> ADB;
typedef ADB::V V;
typedef ADB::M M;


namespace {
    std::vector<int>
    buildAllCells(const int nc)
    {
        std::vector<int> all_cells(nc);

        for (int c = 0; c < nc; ++c) { all_cells[c] = c; }

        return all_cells;
    }

    template <class GeoProps>
    AutoDiffBlock<double>::M
    gravityOperator(const UnstructuredGrid& grid,
                    const HelperOps&        ops ,
                    const GeoProps&         geo )
    {
        using namespace Opm::AutoDiffGrid;
        const int nc = numCells(grid);

        std::vector<int> f2hf(2 * numFaces(grid), -1);
        Eigen::Array<int, Eigen::Dynamic, 2, Eigen::RowMajor>
            face_cells = faceCellsToEigen(grid);

        typedef typename Opm::UgGridHelpers::Cell2FacesTraits<UnstructuredGrid>::Type
            Cell2Faces;
        Cell2Faces c2f=cell2Faces(grid);
        for (int c = 0; c < nc; ++c) {
            typename Cell2Faces::row_type
                cell_faces = c2f[c];
            typedef typename Cell2Faces::row_type::iterator Iter;
            for (Iter f=cell_faces.begin(), end=cell_faces.end();
                 f!=end; ++f) {
                const int p = 0 + (face_cells(*f,0) != c);
                f2hf[2*(*f) + p] = f-c2f[0].begin();
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
            const int c1 = face_cells(f,0);
            const int c2 = face_cells(f,1);

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

    V computePerfPress(const UnstructuredGrid& grid, const Wells& wells, const V& rho, const double grav)
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
                const double cell_depth = grid.cell_centroids[dim * cell + dim - 1];
                wdp[j] = rho[j]*grav*(cell_depth - ref_depth);
            }
        }
        return wdp;
    }

} // anonymous namespace





    ImpesTPFAAD::ImpesTPFAAD(const UnstructuredGrid&         grid,
                             const BlackoilPropsAdInterface& fluid,
                             const DerivedGeology&           geo,
                             const Wells&                    wells,
                             const LinearSolverInterface&    linsolver)
        : grid_     (grid)
        , fluid_    (fluid)
        , geo_      (geo)
        , wells_    (wells)
        , linsolver_(linsolver)
          // , pdepfdata_(grid.number_of_cells, fluid)
        , ops_      (grid)
        , grav_     (gravityOperator(grid_, ops_, geo_))
        , cell_residual_ (ADB::null())
        , well_flow_residual_ ()
        , well_residual_ (ADB::null())
        , total_residual_ (ADB::null())
        , qs_ (ADB::null())
    {
    }





    void
    ImpesTPFAAD::solve(const double   dt,
                       BlackoilState& state,
                       WellState& well_state)
    {
        using namespace Opm::AutoDiffGrid;
        const int nc = numCells(grid_);
        const int np = state.numPhases();

        well_flow_residual_.resize(np, ADB::null());

        // Compute dynamic data that are treated explicitly.
        computeExplicitData(dt, state, well_state);
        // Compute relperms once and for all (since saturations are explicit).
        DataBlock s = Eigen::Map<const DataBlock>(state.saturation().data(), nc, np);
        assert(np == 2);
        kr_ = fluidRelperm(s.col(0), s.col(1), V::Zero(nc,1), buildAllCells(nc));
        // Compute relperms for wells. This must be revisited for crossflow.
        const int nw = wells_.number_of_wells;
        const int nperf = wells_.well_connpos[nw];
        DataBlock well_s(nperf, np);
        for (int w = 0; w < nw; ++w) {
            const double* comp_frac = &wells_.comp_frac[np*w];
            for (int j = wells_.well_connpos[w]; j < wells_.well_connpos[w+1]; ++j) {
                well_s.row(j) = Eigen::Map<const DataBlock>(comp_frac, 1, np);
            }
        }
        const std::vector<int> well_cells(wells_.well_cells,
                                          wells_.well_cells + nperf);
        well_kr_ = fluidRelperm(well_s.col(0), well_s.col(1), V::Zero(nperf,1), well_cells);

        const double atol  = 1.0e-10;
        const double rtol  = 5.0e-6;
        const int    maxit = 15;

        assemble(dt, state, well_state);

        const double r0  = residualNorm();
        int          it  = 0;
        std::cout << "\nIteration         Residual\n"
                  << std::setw(9) << it << std::setprecision(9)
                  << std::setw(18) << r0 << std::endl;
        bool resTooLarge = r0 > atol;
        while (resTooLarge && (it < maxit)) {
            solveJacobianSystem(state, well_state);

            assemble(dt, state, well_state);

            const double r = residualNorm();

            resTooLarge = (r > atol) && (r > rtol*r0);

            it += 1;

            std::cout << std::setw(9) << it << std::setprecision(9)
                      << std::setw(18) << r << std::endl;
        }

        if (resTooLarge) {
            OPM_THROW(std::runtime_error, "Failed to compute converged pressure solution");
        }
        else {
            computeFluxes(state, well_state);
        }
    }





    void
    ImpesTPFAAD::computeExplicitData(const double         dt,
                                     const BlackoilState& state,
                                     const WellState& well_state)
    {
        using namespace Opm::AutoDiffGrid;
        // Suppress warnings about "unused parameters".
        static_cast<void>(dt);
        static_cast<void>(well_state);

        const int nc = numCells(grid_);
        const int np = state.numPhases();
        const int nw = wells_.number_of_wells;
        const int nperf = wells_.well_connpos[nw];
        const int dim = dimensions(grid_);

        const std::vector<int> cells = buildAllCells(nc);

        // Compute relperms.
        DataBlock s = Eigen::Map<const DataBlock>(state.saturation().data(), nc, np);
        assert(np == 2);
        kr_ = fluidRelperm(s.col(0), s.col(1), V::Zero(nc,1), buildAllCells(nc));

        // Compute relperms for wells. This must be revisited for crossflow.
        DataBlock well_s(nperf, np);
        for (int w = 0; w < nw; ++w) {
            const double* comp_frac = &wells_.comp_frac[np*w];
            for (int j = wells_.well_connpos[w]; j < wells_.well_connpos[w+1]; ++j) {
                well_s.row(j) = Eigen::Map<const DataBlock>(comp_frac, 1, np);
            }
        }
        const std::vector<int> well_cells(wells_.well_cells,
                                          wells_.well_cells + nperf);
        well_kr_ = fluidRelperm(well_s.col(0), well_s.col(1), V::Zero(nperf,1), well_cells);

        // Compute well pressure differentials.
        // Construct pressure difference vector for wells.
        const double* g = geo_.gravity();
        if (g) {
            // Guard against gravity in anything but last dimension.
            for (int dd = 0; dd < dim - 1; ++dd) {
                assert(g[dd] == 0.0);
            }
        }
        V cell_rho_total = V::Zero(nc,1);
        const Eigen::Map<const V> p(state.pressure().data(), nc, 1);
        const Eigen::Map<const V> T(state.temperature().data(), nc, 1);
        for (int phase = 0; phase < np; ++phase) {
            const V cell_rho = fluidRho(phase, p, T, cells);
            const V cell_s = s.col(phase);
            cell_rho_total += cell_s * cell_rho;
        }
        V rho_perf = subset(cell_rho_total, well_cells);
        well_perf_dp_ = computePerfPress(grid_, wells_, rho_perf, g ? g[dim-1] : 0.0);
    }





    void
    ImpesTPFAAD::assemble(const double         dt,
                          const BlackoilState& state,
                          const WellState& well_state)
    {
        using namespace Opm::AutoDiffGrid;
        const V& pv = geo_.poreVolume();
        const int nc = numCells(grid_);              ;
        const int np = state.numPhases();
        const int nw = wells_.number_of_wells;
        const int nperf = wells_.well_connpos[nw];

        const std::vector<int> cells = buildAllCells(nc);

        const Eigen::Map<const DataBlock> z0all(&state.surfacevol()[0], nc, np);
        const DataBlock qall = DataBlock::Zero(nc, np);
        const V delta_t = dt * V::Ones(nc, 1);
        const V transi = subset(geo_.transmissibility(),
                                ops_.internal_faces);
        const std::vector<int> well_cells(wells_.well_cells,
                                          wells_.well_cells + nperf);
        const V transw = Eigen::Map<const V>(wells_.WI, nperf, 1);

        // Initialize AD variables: p (cell pressures) and bhp (well bhp).
        const V p0 = Eigen::Map<const V>(&state.pressure()[0], nc, 1);
        const V T0 = Eigen::Map<const V>(&state.temperature()[0], nc, 1);
        const V bhp0 = Eigen::Map<const V>(&well_state.bhp()[0], nw, 1);
        std::vector<V> vars0 = { p0, bhp0 };
        std::vector<ADB> vars = ADB::variables(vars0);
        const ADB& p = vars[0];
        const ADB T = ADB::constant(T0);
        const ADB& bhp = vars[1];

        // Compute T_ij * (p_i - p_j).
        const ADB nkgradp = transi * (ops_.ngrad * p);

        // Extract variables for perforation cell pressures
        // and corresponding perforation well pressures.
        const ADB p_perfcell = subset(p, well_cells);
        const ADB T_perfcell = subset(T, well_cells);
        // Construct matrix to map wells->perforations.
        M well_to_perf(well_cells.size(), nw);
        typedef Eigen::Triplet<double> Tri;
        std::vector<Tri> w2p;
        for (int w = 0; w < nw; ++w) {
            for (int perf = wells_.well_connpos[w]; perf < wells_.well_connpos[w+1]; ++perf) {
                w2p.emplace_back(perf, w, 1.0);
            }
        }
        well_to_perf.setFromTriplets(w2p.begin(), w2p.end());
        const M perf_to_well = well_to_perf.transpose();
        // Finally construct well perforation pressures and well flows.
        const ADB p_perfwell = well_to_perf*bhp + well_perf_dp_;
        const ADB nkgradp_well = transw * (p_perfcell - p_perfwell);
        const Selector<double> cell_to_well_selector(nkgradp_well.value());

        cell_residual_ = ADB::constant(pv);
        well_residual_ = ADB::constant(V::Zero(nw,1));
        ADB divcontrib_sum = ADB::constant(V::Zero(nc,1));
        qs_ = ADB::constant(V::Zero(nw*np, 1));
        for (int phase = 0; phase < np; ++phase) {
            const ADB cell_b = fluidFvf(phase, p, T, cells);
            const ADB cell_rho = fluidRho(phase, p, T, cells);
            const ADB well_b = fluidFvf(phase, p_perfwell, T_perfcell, well_cells);
            const V kr = fluidKr(phase);
            // Explicitly not asking for derivatives of viscosity,
            // since they are not available yet.
            const V mu = fluidMu(phase, p.value(), T.value(), cells);
            const V cell_mob = kr / mu;
            const ADB head_diff_grav = (grav_ * cell_rho);
            const ADB head = nkgradp + (grav_ * cell_rho);
            const UpwindSelector<double> upwind(grid_, ops_, head.value());
            const V face_mob = upwind.select(cell_mob);
            const V well_kr = fluidKrWell(phase);
            const V well_mu = fluidMu(phase, p_perfwell.value(), T_perfcell.value(), well_cells);
            const V well_mob = well_kr / well_mu;
            const V perf_mob = cell_to_well_selector.select(subset(cell_mob, well_cells), well_mob);
            const ADB flux = face_mob * head;
            const ADB perf_flux = perf_mob * (nkgradp_well); // No gravity term for perforations.
            const ADB face_b = upwind.select(cell_b);
            const ADB perf_b = cell_to_well_selector.select(subset(cell_b, well_cells), well_b);
            const V z0 = z0all.block(0, phase, nc, 1);
            const V q  = qall .block(0, phase, nc, 1);
            const ADB well_contrib = superset(perf_flux*perf_b, well_cells, nc);
            const ADB divcontrib = delta_t * (ops_.div * (flux * face_b) + well_contrib);
            const V qcontrib = delta_t * q;
            const ADB pvcontrib = ADB::constant(pv*z0);
            const ADB component_contrib = pvcontrib + qcontrib;
            divcontrib_sum = divcontrib_sum - divcontrib/cell_b;
            cell_residual_ = cell_residual_ - (component_contrib/cell_b);
            const ADB well_rates = perf_to_well * (perf_flux*perf_b);
            qs_ = qs_ +  superset(well_rates, Span(nw, 1, phase*nw), nw*np);
        }
        cell_residual_ = cell_residual_ + divcontrib_sum;
        // Handling BHP and SURFACE_RATE wells.
        V bhp_targets(nw,1);
        V rate_targets(nw,1);
        M rate_distr(nw, np*nw);
        for (int w = 0; w < nw; ++w) {
            const WellControls* wc = wells_.ctrls[w];
            if (well_controls_get_current_type(wc) == BHP) {
                bhp_targets[w] = well_controls_get_current_target( wc );
                rate_targets[w] = -1e100;
            } else if (well_controls_get_current_type(wc) == SURFACE_RATE) {
                bhp_targets[w] = -1e100;
                rate_targets[w] = well_controls_get_current_target( wc );
                {
                    const double * distr = well_controls_get_current_distr( wc );
                    for (int phase = 0; phase < np; ++phase) {
                        rate_distr.insert(w, phase*nw + w) = distr[phase];
                    }
                }
            } else {
                OPM_THROW(std::runtime_error, "Can only handle BHP and SURFACE_RATE type controls.");
            }
        }
        const ADB bhp_residual = bhp - bhp_targets;
        const ADB rate_residual = rate_distr * qs_ - rate_targets;
        // Choose bhp residual for positive bhp targets.
        Selector<double> bhp_selector(bhp_targets);
        well_residual_ = bhp_selector.select(bhp_residual, rate_residual);

        // Build full residual by concatenation of residual arrays and
        // jacobian matrices.
        total_residual_ = collapseJacs(vertcat(cell_residual_, well_residual_));

        // std::cout.precision(16);
        // std::cout << total_residual_;
    }





    void
    ImpesTPFAAD::solveJacobianSystem(BlackoilState& state,
                                     WellState& well_state) const
    {
        using namespace Opm::AutoDiffGrid;
        const int nc = numCells(grid_);
        const int nw = wells_.number_of_wells;
        // const int np = state.numPhases();

        Eigen::SparseMatrix<double, Eigen::RowMajor> matr = total_residual_.derivative()[0];

        V dx(V::Zero(total_residual_.size()));
        Opm::LinearSolverInterface::LinearSolverReport rep
            = linsolver_.solve(matr.rows(), matr.nonZeros(),
                               matr.outerIndexPtr(), matr.innerIndexPtr(), matr.valuePtr(),
                               total_residual_.value().data(), dx.data());
        if (!rep.converged) {
            OPM_THROW(LinearSolverProblem, "ImpesTPFAAD::solve(): Linear solver convergence failure.");
        }
        const V p0 = Eigen::Map<const V>(&state.pressure()[0], nc, 1);
        const V dp = subset(dx, Span(nc));
        const V p = p0 - dp;
        std::copy(&p[0], &p[0] + nc, state.pressure().begin());
        const V bhp0 = Eigen::Map<const V>(&well_state.bhp()[0], nw, 1);
        Span bhp_dofs(nw, 1, nc);
        const V dbhp = subset(dx, bhp_dofs);
        const V bhp = bhp0 - dbhp;
        std::copy(&bhp[0], &bhp[0] + nw, well_state.bhp().begin());
    }





    double
    ImpesTPFAAD::residualNorm() const
    {
        return total_residual_.value().matrix().norm();
    }





    void
    ImpesTPFAAD::computeFluxes(BlackoilState& state,
                               WellState& well_state) const
    {
        using namespace Opm::AutoDiffGrid;
        // This method computes state.faceflux(),
        // well_state.perfRates() and well_state.perfPress().
        const int nc = numCells(grid_);
        const int np = state.numPhases();
        const int nw = wells_.number_of_wells;
        const int nperf = wells_.well_connpos[nw];

        // Build cell sets.
        const std::vector<int> cells = buildAllCells(nc);
        const std::vector<int> well_cells(wells_.well_cells,
                                          wells_.well_cells + nperf);
        // Construct matrix to map wells->perforations.
        M well_to_perf(well_cells.size(), nw);
        typedef Eigen::Triplet<double> Tri;
        std::vector<Tri> w2p;
        for (int w = 0; w < nw; ++w) {
            for (int perf = wells_.well_connpos[w]; perf < wells_.well_connpos[w+1]; ++perf) {
                w2p.emplace_back(perf, w, 1.0);
            }
        }
        well_to_perf.setFromTriplets(w2p.begin(), w2p.end());
        const M perf_to_well = well_to_perf.transpose();
        const V transw = Eigen::Map<const V>(wells_.WI, nperf, 1);

        const V p = Eigen::Map<const V>(&state.pressure()[0], nc, 1);
        const V T = Eigen::Map<const V>(&state.temperature()[0], nc, 1);
        const V bhp = Eigen::Map<const V>(&well_state.bhp()[0], nw, 1);

        const V p_perfcell = subset(p, well_cells);
        const V T_perfcell = subset(T, well_cells);

        const V transi  = subset(geo_.transmissibility(),
                                 ops_.internal_faces);
        const V nkgradp = transi * (ops_.ngrad * p.matrix()).array();

        const V p_perfwell = (well_to_perf*bhp.matrix()).array() + well_perf_dp_;
        const V nkgradp_well = transw * (p_perfcell - p_perfwell);
        const Selector<double> cell_to_well_selector(nkgradp_well);

        V flux = V::Zero(ops_.internal_faces.size(), 1);
        V perf_flux = V::Zero(nperf, 1);

        for (int phase = 0; phase < np; ++phase) {
            const V cell_rho = fluidRho(phase, p, T, cells);
            const V head = nkgradp + (grav_ * cell_rho.matrix()).array();
            const UpwindSelector<double> upwind(grid_, ops_, head);
            const V kr = fluidKr(phase);
            const V mu = fluidMu(phase, p, T, cells);
            const V cell_mob = kr / mu;
            const V face_mob = upwind.select(cell_mob);
            const V well_kr = fluidKrWell(phase);
            const V well_mu = fluidMu(phase, p_perfwell, T_perfcell, well_cells);
            const V well_mob = well_kr / well_mu;
            const V perf_mob = cell_to_well_selector.select(subset(cell_mob, well_cells), well_mob);

            perf_flux += perf_mob * (nkgradp_well); // No gravity term for perforations.
            flux += face_mob * head;
        }

        V all_flux = superset(flux, ops_.internal_faces, numFaces(grid_));
        std::copy(all_flux.data(), all_flux.data() + numFaces(grid_), state.faceflux().begin());

        perf_flux = -perf_flux; // well_state.perfRates() assumed to be inflows.
        std::copy(perf_flux.data(), perf_flux.data() + nperf, well_state.perfRates().begin());

        std::copy(p_perfwell.data(), p_perfwell.data() + nperf, well_state.perfPress().begin());
        std::copy(qs_.value().data(), qs_.value().data() + np*nw, &well_state.wellRates()[0]);
    }





    V ImpesTPFAAD::fluidMu(const int phase, const V& p, const V& T, const std::vector<int>& cells) const
    {
        return fluidMu(phase, ADB::constant(p), ADB::constant(T), cells).value();
    }





    ADB ImpesTPFAAD::fluidMu(const int phase, const ADB& p, const ADB& T, const std::vector<int>& cells) const
    {
        switch (phase) {
        case Water:
            return fluid_.muWat(p, T, cells);
        case Oil: {
            ADB dummy_rs = V::Zero(p.size(), 1) * p;
            std::vector<PhasePresence> cond(dummy_rs.size());
            return fluid_.muOil(p, T, dummy_rs, cond, cells);
        }
        case Gas: {
            ADB dummy_rv = V::Zero(p.size(), 1) * p;
            std::vector<PhasePresence> cond(dummy_rv.size());
            return fluid_.muGas(p, T, dummy_rv, cond, cells);
        }
        default:
            OPM_THROW(std::runtime_error, "Unknown phase index " << phase);
        }
    }





    V ImpesTPFAAD::fluidFvf(const int phase, const V& p, const V& T, const std::vector<int>& cells) const
    {
        return fluidFvf(phase, ADB::constant(p), ADB::constant(T), cells).value();
    }





    ADB ImpesTPFAAD::fluidFvf(const int phase, const ADB& p, const ADB& T, const std::vector<int>& cells) const
    {
        switch (phase) {
        case Water:
            return fluid_.bWat(p, T, cells);
        case Oil: {
            ADB dummy_rs = V::Zero(p.size(), 1) * p;
            std::vector<PhasePresence> cond(dummy_rs.size());
            return fluid_.bOil(p, T, dummy_rs, cond, cells);
        }
        case Gas: {
            ADB dummy_rv = V::Zero(p.size(), 1) * p;
            std::vector<PhasePresence> cond(dummy_rv.size());
            return fluid_.bGas(p, T, dummy_rv, cond, cells);
        }
        default:
            OPM_THROW(std::runtime_error, "Unknown phase index " << phase);
        }
    }





    V ImpesTPFAAD::fluidRho(const int phase, const V& p, const V& T, const std::vector<int>& cells) const
    {
        const double* rhos = fluid_.surfaceDensity();
        V b = fluidFvf(phase, p, T, cells);
        V rho = V::Constant(p.size(), 1, rhos[phase]) * b;
        return rho;
    }





    ADB ImpesTPFAAD::fluidRho(const int phase, const ADB& p, const ADB& T, const std::vector<int>& cells) const
    {
        const double* rhos = fluid_.surfaceDensity();
        ADB b = fluidFvf(phase, p, T, cells);
        ADB rho = V::Constant(p.size(), 1, rhos[phase]) * b;
        return rho;
    }





    std::vector<V> ImpesTPFAAD::fluidRelperm(const V& sw,
                                             const V& so,
                                             const V& sg,
                                             const std::vector<int>& cells) const
    {
        std::vector<ADB> kr_ad = fluid_.relperm(ADB::constant(sw), ADB::constant(so), ADB::constant(sg), cells);
        std::vector<V> kr = { kr_ad[0].value(), kr_ad[1].value(), kr_ad[2].value() };
        return kr;
    }





    V ImpesTPFAAD::fluidKr(const int phase) const
    {
        return kr_[phase];
    }





    V ImpesTPFAAD::fluidKrWell(const int phase) const
    {
        return well_kr_[phase];
    }

} // namespace Opm

