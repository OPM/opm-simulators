/*
  Copyright 2013 SINTEF ICT, Applied Mathematics.
  Copyright 2013 Statoil ASA.

  This file is part of the Open Porous Media Project (OPM).

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

#ifndef OPM_IMPESTPFAAD_HEADER_INCLUDED
#define OPM_IMPESTPFAAD_HEADER_INCLUDED

#include <opm/autodiff/AutoDiffBlock.hpp>
#include <opm/autodiff/AutoDiffHelpers.hpp>

#include <opm/core/simulator/BlackoilState.hpp>
#include <opm/core/simulator/WellState.hpp>
#include <opm/core/utility/ErrorMacros.hpp>
#include <opm/core/linalg/LinearSolverInterface.hpp>
#include <opm/core/wells.h>

#include <algorithm>
#include <cassert>
#include <vector>

#include <boost/shared_ptr.hpp>

struct UnstructuredGrid;

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
}

namespace Opm {

#if 0
    template <typename Scalar, class BOFluid>
    class PressureDependentFluidData {
    public:
        typedef AutoDiff::ForwardBlock<Scalar> ADB;
        typedef typename AutoDiff::ForwardBlock<Scalar>::V V;
        typedef typename AutoDiff::ForwardBlock<Scalar>::M M;

        PressureDependentFluidData(const int      nc,
                                   const BOFluid& fluid)
            : nc_   (nc)
            , np_   (fluid.numPhases())
            , cells_(buildAllCells(nc))
            , fluid_(fluid)
            , A_    (nc_, np_ * np_)
            , dA_   (nc_, np_ * np_)
            , mu_   (nc_, np_      )
            , dmu_  (nc_, np_      )
            , kr_   (nc_, np_      )
            , zero_ (ADB::V::Zero(nc, 1))
            , one_  (ADB::V::Ones(nc, 1))
        {
        }

        void
        computeSatQuant(const BlackoilState& state)
        {
            const std::vector<double>& s = state.saturation();

            assert (s.size() == std::vector<double>::size_type(nc_ * np_));

            double* dkrds = 0;  // Ignore rel-perm derivatives
            fluid_.relperm(nc_, & s[0], & cells_[0],
                           kr_.data(), dkrds);
        }

        void
        computePressQuant(const BlackoilState& state)
        {
            const std::vector<double>& p = state.pressure();
            const std::vector<double>& z = state.surfacevol();

            assert (p.size() == std::vector<double>::size_type(nc_ * 1  ));
            assert (z.size() == std::vector<double>::size_type(nc_ * np_));

            fluid_.matrix   (nc_, & p[0], & z[0], & cells_[0],
                             A_ .data(), dA_ .data());

            fluid_.viscosity(nc_, & p[0], & z[0], & cells_[0],
                             mu_.data(), /*dmu_.data()*/ 0);
            dmu_.fill(0.0);
        }

        ADB
        fvf(const int phase, const ADB& p) const
        {
            assert (0     <= phase);
            assert (phase <  np_  );

            typedef typename ADB::V V;

            const V A  = A_ .block(0, phase * (np_ + 1), nc_, 1);
            const V dA = dA_.block(0, phase * (np_ + 1), nc_, 1);

            std::vector<typename ADB::M> jac(p.numBlocks());
            assert(p.numBlocks() == 2);
            jac[0] = spdiag(dA);
            assert(jac[0].cols() == p.blockPattern()[0]);
            jac[1] = M(A.rows(), p.blockPattern()[1]);
            return one_ / ADB::function(A, jac);
        }

        typename ADB::V
        phaseRelPerm(const int phase) const
        {
            const typename ADB::V kr = kr_.block(0, phase, nc_, 1);

            return kr;
        }

        ADB
        phaseViscosity(const int phase, const ADB& p) const
        {
            assert (0     <= phase);
            assert (phase <  np_  );

            typedef typename ADB::V V;

            const V mu  = mu_ .block(0, phase, nc_, 1);
            const V dmu = dmu_.block(0, phase, nc_, 1);

            std::vector<typename ADB::M> jac(p.numBlocks());
            assert(p.numBlocks() == 2);
            jac[0] = spdiag(dmu);
            assert(jac[0].cols() == p.blockPattern()[0]);
            jac[1] = M(mu.rows(), p.blockPattern()[1]);

            return ADB::function(mu, jac);
        }

        ADB
        phaseDensity(const int phase, const ADB& p) const
        {
            typedef typename ADB::V V;

            const double* rho0 = fluid_.surfaceDensity();

            V rho  = V::Zero(nc_, 1);
            V drho = V::Zero(nc_, 1);
            for (int i = 0; i < np_; ++i) {
                rho  += rho0[i] * A_ .block(0, phase*np_ + i, nc_, 1);
                drho += rho0[i] * dA_.block(0, phase*np_ + i, nc_, 1);
            }

            assert (p.numBlocks() == 2);
            std::vector<typename ADB::M> jac(p.numBlocks());
            jac[0] = spdiag(drho);
            jac[1] = M(rho.rows(), p.blockPattern()[1]);

            assert (jac[0].cols() == p.blockPattern()[0]);

            return ADB::function(rho, jac);
        }

    private:
        const int nc_;
        const int np_;

        const std::vector<int> cells_;
        const BOFluid&         fluid_;

        typedef Eigen::Array<Scalar,
                             Eigen::Dynamic,
                             Eigen::Dynamic,
                             Eigen::RowMajor> DerivedQuant;

        // Pressure dependent quantities (essentially B and \mu)
        DerivedQuant A_  ;
        DerivedQuant dA_ ;
        DerivedQuant mu_ ;
        DerivedQuant dmu_;

        // Saturation dependent quantities (rel-perm only)
        DerivedQuant kr_;

        const typename ADB::V zero_;
        const typename ADB::V one_ ;
    };
#endif




    template <class BOFluid, class GeoProps>
    class ImpesTPFAAD {
    public:
        ImpesTPFAAD(const UnstructuredGrid& grid ,
                    const BOFluid&          fluid,
                    const GeoProps&         geo  ,
                    const Wells&            wells,
                    const LinearSolverInterface& linsolver)
            : grid_     (grid)
            , fluid_    (fluid)
            , geo_      (geo)
            , wells_    (wells)
            , linsolver_(linsolver)
            // , pdepfdata_(grid.number_of_cells, fluid)
            , ops_      (grid)
            , grav_     (gravityOperator(grid_, ops_, geo_))
            , cell_residual_ (ADB::null())
            , well_residual_ (ADB::null())
        {
        }

        void
        solve(const double   dt,
              BlackoilState& state,
              WellState& well_state)
        {
            // pdepfdata_.computeSatQuant(state);
            const int nc = grid_.number_of_cells;
            const int np = state.numPhases();
            DataBlock s = Eigen::Map<const DataBlock>(state.saturation().data(), nc, np);
            ASSERT(np == 2);
            kr_ = fluid_.relperm(s.col(0), s.col(1), V::Zero(nc,1), buildAllCells(nc));

            const double atol  = 1.0e-15;
            const double rtol  = 5.0e-10;
            const int    maxit = 15;

            assemble(dt, state, well_state);

            const double r0  = residualNorm();
            int          it  = 0;
            bool resTooLarge = r0 > atol;
            while (resTooLarge && (it < maxit)) {
                solveJacobianSystem(state);

                assemble(dt, state, well_state);

                const double r = residualNorm();

                resTooLarge = (r > atol) && (r > rtol*r0);

                it += 1;
            }

            if (resTooLarge) {
                THROW("Failed to compute converged pressure solution");
            }
            else {
                computeFluxes(state);
            }
        }

    private:
        // Disallow copying and assignment
        ImpesTPFAAD(const ImpesTPFAAD& rhs);
        ImpesTPFAAD& operator=(const ImpesTPFAAD& rhs);

        // typedef PressureDependentFluidData<double, BOFluid> PDepFData;
        // typedef typename PDepFData::ADB ADB;
        typedef AutoDiff::ForwardBlock<double> ADB;
        typedef ADB::V V;
        typedef ADB::M M;
        typedef Eigen::Array<double,
                             Eigen::Dynamic,
                             Eigen::Dynamic,
                             Eigen::RowMajor> DataBlock;

        const UnstructuredGrid& grid_;
        const BOFluid&          fluid_;
        const GeoProps&         geo_ ;
        const Wells&            wells_;
        const LinearSolverInterface& linsolver_;
        // PDepFData               pdepfdata_;
        HelperOps               ops_;
        const M                 grav_;
        ADB                     cell_residual_;
        ADB                     well_residual_;
        std::vector<V> kr_;

        enum { Water = BOFluid::Water,
               Oil = BOFluid::Oil,
               Gas = BOFluid::Gas };

        void
        assemble(const double         dt,
                 const BlackoilState& state,
                 const WellState& well_state)
        {

            const V& pv = geo_.poreVolume();
            const int nc = grid_.number_of_cells;
            const int np = state.numPhases();
            const int nw = wells_.number_of_wells;

            // pdepfdata_.computePressQuant(state);

            const std::vector<int> cells = buildAllCells(nc);

            const Eigen::Map<const DataBlock> z0all(&state.surfacevol()[0], nc, np);
            const DataBlock qall = DataBlock::Zero(nc, np);
            const V delta_t = dt * V::Ones(nc, 1);
            const V transi = subset(geo_.transmissibility(),
                                    ops_.internal_faces);
            const int num_perf = wells_.well_connpos[nw];
            const std::vector<int> well_cells(wells_.well_cells,
                                              wells_.well_cells + num_perf);
            const V transw = Eigen::Map<const V>(wells_.WI, num_perf, 1);

            // Initialize AD variables: p (cell pressures) and bhp (well bhp).
            const V p0 = Eigen::Map<const V>(&state.pressure()[0], nc, 1);
            const V bhp0 = Eigen::Map<const V>(&well_state.bhp()[0], nw, 1);
            std::vector<V> vars0 = { p0, bhp0 };
            std::vector<ADB> vars= ADB::variables(vars0);
            const ADB& p = vars[0];
            const ADB& bhp = vars[1];
            std::vector<int> bpat = p.blockPattern();

            // Compute T_ij * (p_i - p_j) and use for upwinding.
            const ADB nkgradp = transi * (ops_.ngrad * p);
            const UpwindSelector<double> upwind(grid_, ops_, nkgradp.value());

            // Extract variables for perforation cell pressures
            // and corresponding perforation well pressures.
            const ADB p_perfcell = subset(p, well_cells);
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
            // Construct pressure difference vector for wells.
            const V well_perf_dp = V::Zero(well_cells.size()); // No gravity yet!
            // Finally construct well perforation pressures.
            const ADB p_perfwell = well_to_perf*bhp + well_perf_dp;

            const ADB nkgradp_well = transw * (p_perfcell - p_perfwell);

            cell_residual_ = ADB::constant(pv, bpat);
#define COMPENSATE_FLOAT_PRECISION 0
#if COMPENSATE_FLOAT_PRECISION
            ADB divcontrib_sum = ADB::constant(V::Zero(nc,1), bpat);
#endif
            for (int phase = 0; phase < np; ++phase) {
                // const ADB cell_B = pdepfdata_.fvf(phase, p);
                // const ADB cell_rho = pdepfdata_.phaseDensity(phase, p);
                // const V   kr = pdepfdata_.phaseRelPerm(phase);
                // const ADB mu = pdepfdata_.phaseViscosity(phase, p);
                const ADB cell_b = fluidFvf(phase, p, cells);
                const ADB cell_rho = fluidRho(phase, p, cells);
                const V   kr = fluidKr(phase);
                const ADB mu = fluidMu(phase, p, cells);

                const ADB mf = upwind.select(kr / mu);
                const ADB flux = mf * (nkgradp + (grav_ * cell_rho));

                const ADB face_b = upwind.select(cell_b);

                const V z0 = z0all.block(0, phase, nc, 1);
                const V q  = qall .block(0, phase, nc, 1);

#if COMPENSATE_FLOAT_PRECISION
                const ADB divcontrib = delta_t * (ops_.div * (flux * face_b));
                const V qcontrib = delta_t * q;
                const ADB pvcontrib = ADB::constant(pv*z0, bpat);
                const ADB component_contrib = pvcontrib + qcontrib;
                divcontrib_sum = divcontrib_sum - divcontrib/cell_b;
                cell_residual_ = cell_residual_ - (component_contrib/cell_b);
#else
                const ADB component_contrib = pv*z0 + delta_t*(q - (ops_.div * (flux * face_b)));
                cell_residual_ = cell_residual_ - (component_contrib / cell_b);
#endif
            }
#if COMPENSATE_FLOAT_PRECISION
            cell_residual_ = cell_residual_ + divcontrib_sum;
#endif
        }

        void
        solveJacobianSystem(BlackoilState& state) const
        {
            const int nc = grid_.number_of_cells;
            Eigen::SparseMatrix<double, Eigen::RowMajor> matr = cell_residual_.derivative()[0];

#if HACK_INCOMPRESSIBLE_GRAVITY
            matr.coeffRef(0, 0) *= 2;
#endif

            V dp(nc);
            const V p0 = Eigen::Map<const V>(&state.pressure()[0], nc, 1);
            Opm::LinearSolverInterface::LinearSolverReport rep
                = linsolver_.solve(nc, matr.nonZeros(),
                                   matr.outerIndexPtr(), matr.innerIndexPtr(), matr.valuePtr(),
                                   cell_residual_.value().data(), dp.data());
            if (!rep.converged) {
                THROW("ImpesTPFAAD::solve(): Linear solver convergence failure.");
            }
            const V p = p0 - dp;
            std::copy(&p[0], &p[0] + nc, state.pressure().begin());
        }

        double
        residualNorm() const
        {
            return cell_residual_.value().matrix().norm();
        }

        void
        computeFluxes(BlackoilState& state) const
        {
            const int nc = grid_.number_of_cells;
            const int np = state.numPhases();

            const std::vector<int> cells = buildAllCells(nc);

            const V   p0 = Eigen::Map<const V>(&state.pressure()[0], nc, 1);
            const ADB p  = ADB::constant(p0, cell_residual_.blockPattern());

            const V transi  = subset(geo_.transmissibility(),
                                     ops_.internal_faces);
            const V nkgradp = transi * (ops_.ngrad * p0.matrix()).array();

            V flux = V::Zero(ops_.internal_faces.size(), 1);

            for (int phase = 0; phase < np; ++phase) {
                // const ADB cell_rho = pdepfdata_.phaseDensity(phase, p);
                const ADB cell_rho = fluidRho(phase, p, cells);
                const V   head     = nkgradp +
                    (grav_ * cell_rho.value().matrix()).array();

                const UpwindSelector<double> upwind(grid_, ops_, head);

                // const V   kr = pdepfdata_.phaseRelPerm(phase);
                // const ADB mu = pdepfdata_.phaseViscosity(phase, p);
                const V   kr = fluidKr(phase);
                const ADB mu = fluidMu(phase, p, cells);
                const V   mf = upwind.select(kr / mu.value());

                flux += mf * head;
            }
            V all_flux = superset(flux, ops_.internal_faces, grid_.number_of_faces);
            std::copy(all_flux.data(), all_flux.data() + grid_.number_of_faces, state.faceflux().data());
        }


        ADB fluidMu(const int phase, const ADB& p, const std::vector<int>& cells) const
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

        ADB fluidFvf(const int phase, const ADB& p, const std::vector<int>& cells) const
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

        ADB fluidRho(const int phase, const ADB& p, const std::vector<int>& cells) const
        {
            const double* rhos = fluid_.surfaceDensity();
            ADB b = fluidFvf(phase, p, cells);
            ADB rho = V::Constant(p.size(), 1, rhos[phase]) * b;
            return rho;
        }


        V fluidKr(const int phase) const
        {
            return kr_[phase];
        }
    };
} // namespace Opm

#endif  /* OPM_IMPESTPFAAD_HEADER_INCLUDED */
