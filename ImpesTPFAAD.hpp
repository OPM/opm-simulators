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

#include "AutoDiffBlock.hpp"
#include "AutoDiffHelpers.hpp"

#include <opm/core/simulator/BlackoilState.hpp>
#include <opm/core/simulator/WellState.hpp>
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
}

namespace Opm {
    class LinearSolverInterface;

    template <typename Scalar, class BOFluid>
    class PressureDependentFluidData {
    public:
        typedef AutoDiff::ForwardBlock<Scalar> ADB;

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
        }

        ADB
        fvf(const int phase) const
        {
            assert (0     <= phase);
            assert (phase <  np_  );

            typename ADB::V A  = A_ .block(0, phase * (np_ + 1), nc_, 1);
            typename ADB::V dA = dA_.block(0, phase * (np_ + 1), nc_, 1);

            std::vector<typename ADB::M> jac(1, spdiag(dA));

            return one_ / ADB::function(A, jac);
        }

        typename ADB::V
        phaseRelPerm(const int phase) const
        {
            typename ADB::V kr = kr_.block(0, phase, nc_, 1);

            return kr;
        }

        ADB
        phaseViscosity(const int phase) const
        {
            assert (0     <= phase);
            assert (phase <  np_  );

            typename ADB::V mu  = mu_ .block(0, phase, nc_, 1);
            typename ADB::V dmu = dmu_.block(0, phase, nc_, 1);

            std::vector<typename ADB::M> jac(1, spdiag(dmu));

            return ADB::function(mu, jac);
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

    template <class BOFluid, class GeoProps>
    class ImpesTPFAAD {
    public:
        ImpesTPFAAD(const UnstructuredGrid& grid ,
                    const BOFluid&          fluid,
                    const GeoProps&         geo  ,
                    const Wells&            wells)
            : grid_     (grid)
            , geo_      (geo)
            , wells_    (wells)
            , pdepfdata_(grid.number_of_cells, fluid)
            , ops_      (grid)
        {
        }

        void
        solve(const double   dt,
              BlackoilState& state,
              WellState& well_state)
        {
            pdepfdata_.computeSatQuant(state);

            assemble(dt, state, well_state);
        }

    private:
        // Disallow copying and assignment
        ImpesTPFAAD(const ImpesTPFAAD& rhs);
        ImpesTPFAAD& operator=(const ImpesTPFAAD& rhs);

        typedef PressureDependentFluidData<double, BOFluid> PDepFData;
        typedef typename PDepFData::ADB ADB;

        const UnstructuredGrid& grid_;
        const GeoProps&         geo_ ;
        const Wells&            wells_;
        PDepFData               pdepfdata_;
        HelperOps               ops_;

        void
        assemble(const double         dt,
                 const BlackoilState& state,
                 const WellState& well_state)
        {
            typedef typename ADB::V V;

            typedef Eigen::Array<double,
                                 Eigen::Dynamic,
                                 Eigen::Dynamic,
                                 Eigen::RowMajor> DataBlock;

            const V& pv = geo_.poreVolume();
            const int nc = grid_.number_of_cells;
            const int np = state.numPhases();
            const int nw = wells_.number_of_wells;

            pdepfdata_.computePressQuant(state);

            const Eigen::Map<const DataBlock> z0all(&state.surfacevol()[0], nc, np);
            const DataBlock qall = DataBlock::Zero(nc, np);

            const V transi = subset(geo_.transmissibility(),
                                    ops_.internal_faces);
            const V p0 = Eigen::Map<const V>(&state.pressure()[0], nc, 1);
            const V bhp = Eigen::Map<const V>(&well_state.bhp()[0], nw, 1);

            const std::vector<int> bpat(1, nc);
            ADB p = ADB::variable(0, p0, bpat);
            const ADB nkgradp = transi * (ops_.ngrad * p);

            const UpwindSelector<double> upwind(grid_, ops_, nkgradp.value());

            const V delta_t = dt * V::Ones(nc, 1);

            ADB residual = ADB::constant(pv, bpat);
            for (int phase = 0; phase < np; ++phase) {
                const ADB cell_B = pdepfdata_.fvf(phase);

                const V   kr = pdepfdata_.phaseRelPerm(phase);
                const ADB mu = pdepfdata_.phaseViscosity(phase);
                const ADB mf = upwind.select(kr / mu);
                const ADB flux = mf * nkgradp;

                const ADB face_B = upwind.select(cell_B);

                const V z0 = z0all.block(0, phase, nc, 1);
                const V q  = qall .block(0, phase, nc, 1);

                ADB component_contrib = pv*z0 + delta_t*(q - (ops_.div * (flux / face_B)));
                residual = residual - (cell_B * component_contrib);
            }
        }
    };
} // namespace Opm

#endif  /* OPM_IMPESTPFAAD_HEADER_INCLUDED */
