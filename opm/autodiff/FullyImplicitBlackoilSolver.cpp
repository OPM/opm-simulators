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

#include <opm/core/simulator/BlackoilState.hpp>
#include <opm/core/simulator/WellState.hpp>
#include <opm/core/grid.h>
#include <opm/core/utility/ErrorMacros.hpp>


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
                                const Wells&                    wells)
        : grid_  (grid)
        , fluid_ (fluid)
        , geo_   (geo)
        , wells_ (wells)
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
            const SolutionState state = constantState(x);
            computeAccum(state, 0);
        }

#if 0
        const double atol  = 1.0e-15;
        const double rtol  = 5.0e-10;
        const int    maxit = 15;
#endif

        assemble(dtpv, x, xw);

#if 0
        const double r0  = residualNorm();
        int          it  = 0;
        bool resTooLarge = r0 > atol;
        while (resTooLarge && (it < maxit)) {
            solveJacobianSystem(x);

            assemble(dtpv, x, xw);

            const double r = residualNorm();

            resTooLarge = (r > atol) && (r > rtol*r0);

            it += 1;
        }

        if (resTooLarge) {
            THROW("Failed to compute converge solution");
        }
#endif
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
    FullyImplicitBlackoilSolver::constantState(const BlackoilState& x)
    {
        const int nc = grid_.number_of_cells;
        const int np = x.numPhases();

        const std::vector<int> bpat(np, nc);

        SolutionState state(np);

        assert (not x.pressure().empty());
        const V p = Eigen::Map<const V>(& x.pressure()[0], nc, 1);
        state.pressure = ADB::constant(p, bpat);

        assert (not x.saturation().empty());
        const DataBlock s =
            Eigen::Map<const DataBlock>(& x.saturation()[0], nc, np);

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

        // Ignore miscibility effects (no dissolved gas) for now!
        const V Rs = V::Zero(nc, 1);
        state.Rs = ADB::constant(Rs, bpat);

        return state;
    }

    FullyImplicitBlackoilSolver::SolutionState
    FullyImplicitBlackoilSolver::variableState(const BlackoilState& x)
    {
        const int nc = grid_.number_of_cells;
        const int np = x.numPhases();

        std::vector<V> vars0;
        vars0.reserve(np);

        assert (not x.pressure().empty());
        const V p = Eigen::Map<const V>(& x.pressure()[0], nc, 1);
        vars0.push_back(p);

        assert (not x.saturation().empty());
        const DataBlock s =
            Eigen::Map<const DataBlock>(& x.saturation()[0], nc, np);

        const Opm::PhaseUsage pu = fluid_.phaseUsage();

        if (active_[ Water ]) {
            const V sw = s.col(pu.phase_pos[ Water ]);
            vars0.push_back(sw);
        }
        if (active_[ Gas ]) {
            const V sg = s.col(pu.phase_pos[ Gas ]);
            vars0.push_back(sg);
        }

        std::vector<ADB> vars = ADB::variables(vars0);
            
        SolutionState state(np);
        state.pressure = vars[0];

        const std::vector<int>& bpat = vars[0].blockPattern();
        {
            ADB so  = ADB::constant(V::Ones(nc, 1), bpat);
            int off = 1;    // First saturation variable at offset 1.

            if (active_[ Water ]) {
                ADB& sw = vars[ off++ ];
                state.saturation[ pu.phase_pos[ Water ] ] = sw;

                so = so - sw;
            }
            if (active_[ Gas ]) {
                ADB& sg = vars[ off++ ];
                state.saturation[ pu.phase_pos[ Gas ] ] = sg;

                so = so - sg;
            }
            if (active_[ Oil ]) {
                state.saturation[ pu.phase_pos[ Oil ] ] = so;
            }
        }

        // Ignore miscibility effects (no dissolved gas) for now!
        V Rs = V::Zero(nc, 1);
        state.Rs = ADB::constant(Rs, bpat);

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
        const V transi = subset(geo_.transmissibility(),
                                ops_.internal_faces);

        const SolutionState    state = variableState(x);
        const std::vector<ADB> kr    = computeRelPerm(state);

        computeAccum(state, 1);

        for (int phase = 0; phase < fluid_.numPhases(); ++phase) {
            computeMassFlux(phase, transi, kr, state);

            residual_.mass_balance[ phase ] =
                dtpv*(rq_[phase].accum[1] - rq_[phase].accum[0])
                + ops_.div*rq_[phase].mflux;
        }

        if (active_[ Oil ] && active_[ Gas ]) {
            const int po = fluid_.phaseUsage().phase_pos[ Oil ];
            const UpwindSelector<double> upwind(grid_, ops_,
                                                rq_[po].head.value());
            const ADB Rs = upwind.select(state.Rs);

            residual_.mass_balance[ Gas ] += ops_.div * (Rs * rq_[po].mflux);
        }
    }

    std::vector<ADB>
    FullyImplicitBlackoilSolver::computeRelPerm(const SolutionState& state)
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

} // namespace Opm

